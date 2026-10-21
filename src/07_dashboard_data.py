"""Compact JSON for the interactive dashboard."""
import json, geopandas as gpd, pandas as pd, numpy as np
from common import PROC, ROOT
g = gpd.read_file(PROC / "baltimore_heat_tracts.gpkg").to_crs(4326)
names = gpd.read_file(PROC / "tracts_social_health.gpkg")[["GEOID", "NAME"]]
g = g.drop(columns=["NAME"], errors="ignore").merge(names, on="GEOID", how="left")
g["geometry"] = g.geometry.simplify(0.00008, preserve_topology=True)
# d3-geo expects clockwise exterior rings (the opposite of RFC 7946); orient them so each tract draws as itself
from shapely.geometry.polygon import orient
from shapely.geometry import MultiPolygon, Polygon
def cw(geom):
    if isinstance(geom, Polygon): return orient(geom, sign=-1.0)
    return MultiPolygon([orient(p, sign=-1.0) for p in geom.geoms])
g["geometry"] = g.geometry.apply(cw)
cols = {'GEOID': 'id', 'NAME': 'name', 'csa': 'csa', 'pop_acs': 'pop', 'hvi': 'hvi', 'hvi_rank': 'rank', 'rank_lo90': 'rlo', 'rank_hi90': 'rhi',
        'p_top20': 'p20', 'hvi_class': 'cls', 'score_exposure': 'sE', 'score_sensitivity': 'sS', 'score_low_capacity': 'sC',
        'lst_popw_c': 'lst', 'modis_myd_night_c': 'night', 'pct_tree': 'tree', 'pct_poverty': 'pov', 'pct_65plus': 'age65',
        'pl_CHD': 'chd', 'pl_COPD': 'copd', 'pl_CASTHMA': 'asthma', 'pl_DIABETES': 'diab', 'pct_no_vehicle': 'noveh',
        'pct_65_alone': 'alone', 'pl_SHUTUTILITY': 'shut', 'dist_cool_km': 'cool', 'holc_grade': 'holc', 'lisa': 'lisa',
        'pct_black_nh': 'black', 'median_income': 'inc', 'svi_themes': 'svi'}
d = g[list(cols) + ['geometry']].rename(columns=cols)
def rnd(v):
    if v is None or (isinstance(v, float) and np.isnan(v)): return None
    if isinstance(v, (float, np.floating)): return round(float(v), 3)
    if isinstance(v, np.integer): return int(v)
    return v
def rc(c): return [rc(x) for x in c] if isinstance(c, list) else round(c, 5)
feats = []
for _, r in d.iterrows():
    geom = json.loads(gpd.GeoSeries([r.geometry]).to_json())['features'][0]['geometry']
    geom['coordinates'] = rc(geom['coordinates'])
    p = {k: rnd(r[k]) for k in cols.values()}
    p['name'] = (p['name'] or '').replace('Census Tract ', '')
    feats.append({'type': 'Feature', 'properties': p, 'geometry': geom})
cool = gpd.read_file(PROC / "cooling_centers.gpkg").to_crs(4326)
cc = [{'n': r.USER_Locat, 'a': r.USER_Stree, 't': r.USER_Type, 'h': r.USER_Hours, 'x': round(r.geometry.x, 5), 'y': round(r.geometry.y, 5)} for _, r in cool.iterrows()]
R = json.load(open(PROC / "results.json"))
Z = pd.read_csv(PROC / "zcta_validation.csv", dtype={'ZCTA5': str})
val = [{'z': r.ZCTA5, 'hvi': round(r.hvi, 3), 'ems': int(r.PR_HRI), 'exp': round(r.score_exposure, 3), 'cap': round(r.score_low_capacity, 3)} for _, r in Z.iterrows()]
keys = ['n_scored', 'lst_range', 'night_range', 'moran_hvi', 'lisa_counts', 'lisa_hh_pop', 'veryhigh_n', 'veryhigh_pop', 'holc_lst_by_grade',
        'holc_tree_by_grade', 'holc_n_by_grade', 'val_spearman', 'val_hvi_boot_ci', 'zcta_n', 'spearman_lst_tree', 'ols_lst_tree_built',
        'rank_stability_spearman', 'n_robust_top20', 'n_top20', 'mc_draws', 'station_summers', 'dmh_minus_bwi_tmin_f', 'landsat_scenes_used',
        'excess_nested_F', 'all_mean_dist_cool', 'veryhigh_mean_dist_cool', 'n_cooling_centers', 'exp_day_black_nh', 'exp_day_white_nh',
        'exp_day_poverty', 'exp_day_all', 'holc_D_minus_A_or_B', 'paired_days_dmh_bwi']
out = {'tracts': {'type': 'FeatureCollection', 'features': feats}, 'cool': cc, 'val': val, 'R': {k: R[k] for k in keys}}
s = json.dumps(out, separators=(',', ':'), default=lambda o: o.item() if hasattr(o, 'item') else str(o))
(ROOT / "dashboard").mkdir(exist_ok=True)
open(ROOT / "dashboard" / "data.json", "w").write(s)
print(round(len(s) / 1024), "KB;", cc[0])
