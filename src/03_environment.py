"""Environmental layers per tract: Landsat LST/NDVI, MODIS day/night LST, ESA WorldCover,
cooling-centre access, HOLC grades, and population-weighted exposure from 2020 census blocks."""
import base64, gzip, io, json
import numpy as np, pandas as pd, geopandas as gpd, rasterio
from rasterio.features import rasterize
from rasterio.warp import reproject, Resampling
from rasterio.transform import from_origin
from common import RAW, PROC, CRS, GRID

tr = gpd.read_file(PROC / "tracts_social_health.gpkg")[["GEOID", "geometry"]].to_crs(CRS).reset_index(drop=True)
tr["tid"] = np.arange(1, len(tr) + 1)
T30 = from_origin(GRID["x0"], GRID["y0"], 30, 30)
shape30 = (GRID["h"], GRID["w"])
ids30 = rasterize(zip(tr.geometry, tr.tid), out_shape=shape30, transform=T30, fill=0, dtype="int32")
out = tr[["GEOID", "tid"]].copy()


def zonal(arr, ids, name, stats=("mean",)):
    s = pd.DataFrame({"tid": ids.ravel(), "v": arr.ravel()})
    s = s[(s.tid > 0) & np.isfinite(s.v)]
    g = s.groupby("tid")["v"]
    res = {}
    if "mean" in stats: res[name] = g.mean()
    if "p90" in stats: res[name + "_p90"] = g.quantile(0.9)
    if "n" in stats: res[name + "_npix"] = g.size()
    return pd.DataFrame(res)


def read(p):
    with rasterio.open(p) as s:
        return s.read(1).astype("float64")

# ---- Landsat composites ----
lst = read(PROC / "lst_median.tif")
out = out.merge(zonal(lst, ids30, "lst_c", ("mean", "p90", "n")), left_on="tid", right_index=True, how="left")
out = out.merge(zonal(read(PROC / "lst_anom_mean.tif"), ids30, "lst_anom_c"), left_on="tid", right_index=True, how="left")
out = out.merge(zonal(read(PROC / "ndvi_median.tif"), ids30, "ndvi"), left_on="tid", right_index=True, how="left")
for y in range(2019, 2026):
    p = PROC / f"lst_median_{y}.tif"
    if p.exists():
        out = out.merge(zonal(read(p), ids30, f"lst_{y}"), left_on="tid", right_index=True, how="left")

# ---- MODIS 8-day LST (Terra 10:30/22:30, Aqua 13:30/01:30 local solar time) ----
M = json.load(gzip.open(RAW / "modis.json.gz"))
dates = dict(M["items"])
acc = {k: [] for k in ["MOD_day", "MOD_night", "MYD_day", "MYD_night"]}
mprof = None
for iid, b64 in M["tifs"].items():
    start = dates[iid][:10]
    if not ("06-01" <= start[5:] <= "08-31"):
        continue
    with rasterio.MemoryFile(base64.b64decode(b64)) as mf, mf.open() as s:
        a = s.read().astype("float64"); mprof = s.profile
    sat = iid[:3]
    for band, qcb, tag in [(0, 2, "day"), (1, 3, "night")]:
        v, qc = a[band], a[qcb].astype(np.uint16)
        good = (v > 0) & (((qc & 3) == 0) | (((qc & 3) == 1) & (((qc >> 6) & 3) <= 1)))
        acc[f"{sat}_{tag}"].append(np.where(good, v * 0.02 - 273.15, np.nan))
modis = {}
for k, lst_ in acc.items():
    st = np.stack(lst_)
    modis[k] = np.nanmedian(st, axis=0)
    print(k, "composites:", len(lst_), "median valid per cell:", int(np.median(np.isfinite(st).sum(0))))
for k, v in modis.items():
    dst = np.full(shape30, np.nan)
    reproject(v, dst, src_transform=mprof["transform"], src_crs=mprof["crs"], dst_transform=T30, dst_crs=CRS,
              resampling=Resampling.nearest, src_nodata=np.nan, dst_nodata=np.nan)
    out = out.merge(zonal(dst, ids30, f"modis_{k.lower()}_c"), left_on="tid", right_index=True, how="left")
    prof = dict(driver="GTiff", height=shape30[0], width=shape30[1], count=1, dtype="float32", crs=CRS, transform=T30, nodata=np.nan)
    with rasterio.open(PROC / f"modis_{k.lower()}.tif", "w", **prof) as d:
        d.write(dst.astype("float32"), 1)

# ---- ESA WorldCover 2021 (10 m) ----
with rasterio.open(RAW / "landsat" / "bhs_worldcover2021.tif") as s:
    wc = s.read(1); wt = s.transform
ids10 = rasterize(zip(tr.geometry, tr.tid), out_shape=wc.shape, transform=wt, fill=0, dtype="int32")
m = ids10 > 0
d = pd.DataFrame({"tid": ids10[m], "c": wc[m]})
d = d[d.c > 0]
land = d[d.c != 80]
frac = land.groupby("tid")["c"].agg(pct_tree=lambda x: 100 * (x == 10).mean(),
                                    pct_built=lambda x: 100 * (x == 50).mean(),
                                    pct_grass=lambda x: 100 * x.isin([20, 30]).mean())
out = out.merge(frac, left_on="tid", right_index=True, how="left")

# ---- misc bundle: blocks, cooling centres, ZCTAs, ACIS ----
B = json.load(open(RAW / "bhs_misc_bundle.json"))
blocks = pd.read_csv(io.StringIO(B["blocks_csv"]), dtype={"GEOID": str})
blocks = gpd.GeoDataFrame(blocks, geometry=gpd.points_from_xy(blocks.INTPTLON, blocks.INTPTLAT), crs=4326).to_crs(CRS)
blocks["tract"] = blocks.GEOID.str[:11]
col = ((GRID["x0"] - 0) * 0 + (blocks.geometry.x - GRID["x0"]) // 30).astype(int)
row = ((GRID["y0"] - blocks.geometry.y) // 30).astype(int)
inside = (col >= 0) & (col < GRID["w"]) & (row >= 0) & (row < GRID["h"])
blocks["lst"] = np.nan
blocks.loc[inside, "lst"] = lst[row[inside], col[inside]]
# block value = mean LST within 90 m (3x3) to avoid single-pixel noise
pad = np.pad(lst, 1, constant_values=np.nan)
win = np.nanmean(np.stack([pad[1 + dy:1 + dy + lst.shape[0], 1 + dx:1 + dx + lst.shape[1]]
                           for dy in (-1, 0, 1) for dx in (-1, 0, 1)]), axis=0)
blocks.loc[inside, "lst3"] = win[row[inside], col[inside]]
bp = blocks[(blocks.POP100 > 0) & blocks.lst3.notna()]
pw = bp.groupby("tract").apply(lambda g: np.average(g.lst3, weights=g.POP100), include_groups=False).rename("lst_popw_c")
out = out.merge(pw, left_on="GEOID", right_index=True, how="left")
blocks.drop(columns="geometry").to_csv(PROC / "blocks_lst.csv", index=False)

cool = gpd.GeoDataFrame.from_features(B["cooling"]["features"], crs=4326).to_crs(CRS)
cool.to_file(PROC / "cooling_centers.gpkg", driver="GPKG")
# population-weighted mean straight-line distance (km) from residents to nearest cooling centre
from scipy.spatial import cKDTree
kd = cKDTree(np.c_[cool.geometry.x, cool.geometry.y])
dist, _ = kd.query(np.c_[blocks.geometry.x, blocks.geometry.y])
blocks["dist_cool_km"] = dist / 1000
bpop = blocks[blocks.POP100 > 0]
dc = bpop.groupby("tract").apply(lambda g: np.average(g.dist_cool_km, weights=g.POP100), include_groups=False).rename("dist_cool_km")
out = out.merge(dc, left_on="GEOID", right_index=True, how="left")

# block -> ZCTA (point in polygon) for the outcome crosswalk
zc = gpd.GeoDataFrame.from_features(B["zcta"]["features"], crs=4326).to_crs(CRS)
zc.to_file(PROC / "zcta2020.gpkg", driver="GPKG")
bz = gpd.sjoin(blocks[["GEOID", "tract", "POP100", "geometry"]], zc[["ZCTA5", "geometry"]], how="left", predicate="within")
bz.drop(columns="geometry").to_csv(PROC / "block_zcta.csv", index=False)

# ACIS station daily data
rows = []
for sid in ["BWI", "DMH"]:
    for dte, tmax, tmin in B[f"acis_{sid}"]["data"]:
        rows.append((sid, dte, tmax, tmin))
ac = pd.DataFrame(rows, columns=["station", "date", "tmax_f", "tmin_f"])
for c in ["tmax_f", "tmin_f"]:
    ac[c] = pd.to_numeric(ac[c], errors="coerce")
ac.to_csv(PROC / "acis_daily.csv", index=False)

# ---- HOLC grades (Mapping Inequality v3, 2020-tract crosswalk) ----
h = gpd.read_file(RAW / "holc_baltimore_tract2020.gpkg")
h = h[h.GEOID.str.startswith("24510") & h.grade.isin(list("ABCD"))]
share = h.pivot_table(index="GEOID", columns="grade", values="pct_tract", aggfunc="sum").fillna(0)
share.columns = [f"holc_{c}" for c in share.columns]
share["holc_cover"] = share.sum(axis=1).clip(upper=1)
wsum = sum(share.get(f"holc_{g}", 0) * w for g, w in zip("ABCD", [1, 2, 3, 4]))
share["holc_score"] = wsum / share[[f"holc_{g}" for g in "ABCD"]].sum(axis=1)
share["holc_grade"] = share[[f"holc_{g}" for g in "ABCD"]].idxmax(axis=1).str[-1]
share.loc[share.holc_cover < 0.5, "holc_grade"] = None
out = out.merge(share, left_on="GEOID", right_index=True, how="left")

out.drop(columns="tid").to_csv(PROC / "tracts_environment.csv", index=False)
print(out.drop(columns=["tid"]).describe().T[["count", "mean", "min", "max"]].round(2).to_string())
