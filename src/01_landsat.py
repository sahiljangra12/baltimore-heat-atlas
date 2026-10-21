"""Landsat 8/9 C2 L2 summer (Jun-Aug) 2019-2025: cloud-masked LST and NDVI composites."""
import json, re, glob
import numpy as np, rasterio
from common import RAW, PROC

QA_REJECT = (1 << 1) | (1 << 2) | (1 << 3) | (1 << 4) | (1 << 7)  # dilated cloud, cirrus, cloud, shadow, water
files = sorted(glob.glob(str(RAW / "landsat" / "bhs_ls_*.tif")))
stack_lst, stack_ndvi, meta, dates = [], [], None, []
for f in files:
    with rasterio.open(f) as s:
        a = s.read().astype("float64"); meta = s.profile
    st, qa, red, nir = a[0], a[1].astype(np.uint16), a[2], a[3]
    ok = (st > 0) & ((qa & QA_REJECT) == 0) & (qa != 1)
    lst = np.where(ok, st * 0.00341802 + 149.0 - 273.15, np.nan)
    # residual cloud / cloud-edge screen: summer daytime land surfaces in Baltimore are never
    # >15 C below the scene median or below 15 C; such pixels are undetected cloud or haze
    if np.isfinite(lst).any():
        med = np.nanmedian(lst)
        bad = (lst < max(15.0, med - 15.0))
        ok = ok & ~bad
        lst = np.where(ok, lst, np.nan)
    r = red * 0.0000275 - 0.2; n = nir * 0.0000275 - 0.2
    with np.errstate(invalid="ignore", divide="ignore"):
        ndvi = np.where(ok & (red > 0) & (nir > 0), (n - r) / (n + r), np.nan)
    frac = np.isfinite(lst).mean()
    d = re.search(r"_(\d{8})_", f).group(1)
    dates.append({"scene": f.split("bhs_ls_")[1][:-4], "date": d, "clear_frac": round(float(frac), 3),
                  "median_lst": round(float(np.nanmedian(lst)), 2) if frac > 0 else None})
    stack_lst.append(lst); stack_ndvi.append(ndvi)
L = np.stack(stack_lst); N = np.stack(stack_ndvi)
clear = np.array([d["clear_frac"] for d in dates])
use = clear >= 0.30          # scenes with at least 30% of the box clear
L, N = L[use], N[use]
# scene anomaly: pixel minus that scene's city-box median (removes day-to-day weather)
anom = L - np.nanmedian(L, axis=(1, 2), keepdims=True)
out = {
    "lst_median": np.nanmedian(L, axis=0),
    "lst_anom_mean": np.nanmean(anom, axis=0),
    "ndvi_median": np.nanmedian(N, axis=0),
    "n_clear": np.isfinite(L).sum(0).astype("float64"),
}
# per-year medians for stability check
years = np.array([d["date"][:4] for d, u in zip(dates, use) if u])
for y in sorted(set(years)):
    out[f"lst_median_{y}"] = np.nanmedian(L[years == y], axis=0)
prof = meta.copy(); prof.update(count=1, dtype="float32", nodata=np.nan)
for k, v in out.items():
    with rasterio.open(PROC / f"{k}.tif", "w", **prof) as dst:
        dst.write(v.astype("float32"), 1)
for d, u in zip(dates, use): d["used"] = bool(u)
json.dump(dates, open(PROC / "landsat_scenes.json", "w"), indent=1)
print(f"{len(files)} scenes, {use.sum()} used; years {sorted(set(years))}")
print("median clear obs per pixel:", np.nanmedian(out["n_clear"]))
print("city-box median LST:", round(float(np.nanmedian(out['lst_median'])), 2))
