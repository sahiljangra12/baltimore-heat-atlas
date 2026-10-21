"""Main analysis: heat vulnerability index, uncertainty, spatial structure, drivers, equity,
and validation against heat-related EMS activations (CDC HHI, ZCTA level)."""
import json, warnings
import numpy as np, pandas as pd, geopandas as gpd
from scipy import stats
import statsmodels.formula.api as smf
from sklearn.decomposition import PCA
from libpysal.weights import Queen
from esda.moran import Moran, Moran_Local
from common import RAW, PROC
warnings.filterwarnings("ignore")
rng = np.random.default_rng(20261021)
R = {}  # every number reported in the paper/dashboard goes in here

g = gpd.read_file(PROC / "tracts_social_health.gpkg")
env = pd.read_csv(PROC / "tracts_environment.csv", dtype={"GEOID": str})
g = g.merge(env, on="GEOID", how="left")
g["scored"] = g["pop2020"] >= 100
# Community Statistical Area (BNIA 2020) with the largest overlap, for readable neighbourhood names
import gzip
csa = gpd.read_file(gzip.open(RAW / "csa2020.geojson.gz")).to_crs(32618)[["CSA2020", "geometry"]]
ov = gpd.overlay(g[["GEOID", "geometry"]].to_crs(32618), csa, how="intersection")
ov["a"] = ov.area
g = g.merge(ov.sort_values("a").drop_duplicates("GEOID", keep="last")[["GEOID", "CSA2020"]], on="GEOID", how="left")
g = g.rename(columns={"CSA2020": "csa"})
R["n_tracts"] = int(len(g)); R["n_scored"] = int(g.scored.sum())
R["pop2020_total"] = int(g.pop2020.sum()); R["pop_acs_total"] = int(g.pop_acs.sum())

# ---------------- indicator framework ----------------
IND = {
    "exposure": {
        "lst_popw_c": "Daytime land surface temperature where people live (Landsat, pop-weighted, °C)",
        "modis_myd_night_c": "Night-time land surface temperature, ~01:30 (MODIS Aqua, °C)",
    },
    "sensitivity": {
        "pct_65plus": "Residents aged 65+ (%)",
        "pct_under5": "Children under 5 (%)",
        "pl_CHD": "Coronary heart disease, adults (%)",
        "pl_COPD": "COPD, adults (%)",
        "pl_CASTHMA": "Current asthma, adults (%)",
        "pl_DIABETES": "Diabetes, adults (%)",
        "pl_STROKE": "Stroke, adults (%)",
        "pl_MHLTH": "Frequent mental distress, adults (%)",
        "pct_disability": "Any disability (%)",
    },
    "low_capacity": {
        "pct_poverty": "Below poverty line (%)",
        "pct_no_vehicle": "Households with no vehicle (%)",
        "pct_65_alone": "Households: one person aged 65+ living alone (%)",
        "pct_limited_english": "Limited-English households (%)",
        "pct_uninsured": "No health insurance (%)",
        "pct_no_internet": "Households with no internet (%)",
        "pl_SHUTUTILITY": "Threat of utility shut-off in past year, adults (%)",
        "pct_pre1940": "Housing built before 1940 (%)",
        "dist_cool_km": "Distance to nearest cooling centre (km, pop-weighted)",
        "low_tree": "Absence of tree cover (100 - % tree, ESA WorldCover)",
    },
}
g["low_tree"] = 100 - g["pct_tree"]
ALL = [c for comp in IND.values() for c in comp]
S = g[g.scored].copy()


def pct_rank(df):
    return df.rank(pct=True)


def build_index(df, weights=None):
    P = pct_rank(df[ALL])
    comps = {k: P[list(v)].mean(axis=1) for k, v in IND.items()}
    C = pd.DataFrame(comps)
    w = weights if weights is not None else np.ones(3) / 3
    return C, (C.values * w).sum(axis=1)

C, hvi = build_index(S)
for k in C: S[f"score_{k}"] = C[k]
S["hvi"] = hvi
S["hvi_rank"] = S.hvi.rank(ascending=False, method="min").astype(int)
S["hvi_pct"] = S.hvi.rank(pct=True)
S["hvi_class"] = pd.qcut(S.hvi, 5, labels=["Very low", "Low", "Moderate", "High", "Very high"])

# indicator correlations with each other (redundancy check)
R["indicator_missing"] = {c: int(S[c].isna().sum()) for c in ALL}
R["component_spearman"] = {f"{a}~{b}": round(stats.spearmanr(S[f"score_{a}"], S[f"score_{b}"]).statistic, 3)
                           for a, b in [("exposure", "sensitivity"), ("exposure", "low_capacity"), ("sensitivity", "low_capacity")]}

# ---------------- alternative constructions ----------------
Z = S[ALL].apply(lambda s: (s - s.mean()) / s.std()).fillna(0)
pca = PCA().fit(Z)
k = int((pca.explained_variance_ > 1).sum())
R["pca_components_kaiser"] = k
R["pca_var_explained"] = round(float(pca.explained_variance_ratio_[:k].sum()), 3)
scores = pca.transform(Z)[:, :k]
# orient each component so that it correlates positively with the equal-weight HVI
for j in range(k):
    if np.corrcoef(scores[:, j], S.hvi)[0, 1] < 0: scores[:, j] *= -1
S["hvi_pca"] = (scores * pca.explained_variance_ratio_[:k]).sum(axis=1)
P = pct_rank(S[ALL])
S["hvi_geo"] = np.exp(np.log(S[[f"score_{c}" for c in IND]].clip(lower=0.01)).mean(axis=1))
S["hvi_flat"] = P.mean(axis=1)  # every indicator equal (no components)
R["alt_index_spearman"] = {n: round(stats.spearmanr(S.hvi, S[n]).statistic, 3) for n in ["hvi_pca", "hvi_geo", "hvi_flat"]}
R["svi_vs_hvi_spearman"] = round(stats.spearmanr(S.hvi, S.svi_themes, nan_policy="omit").statistic, 3)

# ---------------- uncertainty: Monte Carlo over data error + weights ----------------
N = 1000
ranks_data, ranks_w, ranks_both = [], [], []
se = {}
for c in ALL:
    if c + "_moe" in S: se[c] = S[c + "_moe"] / 1.645
    elif c.startswith("pl_") and c + "_se" in S: se[c] = S[c + "_se"]
yrs = [f"lst_{y}" for y in range(2019, 2026) if f"lst_{y}" in S]
# LST year-to-year spread in tract anomaly -> SE of a 7-summer mean
anom_years = S[yrs].sub(S[yrs].mean(axis=0), axis=1)
se["lst_popw_c"] = anom_years.std(axis=1) / np.sqrt(len(yrs))
for it in range(N):
    D = S[ALL].copy()
    for c, s in se.items():
        D[c] = D[c] + rng.normal(0, 1, len(D)) * s.fillna(0).values
    w_rand = rng.dirichlet([4, 4, 4])
    for store, use_data, use_w in [(ranks_data, True, False), (ranks_w, False, True), (ranks_both, True, True)]:
        src = D if use_data else S
        _, h = build_index(src, w_rand if use_w else None)
        store.append(pd.Series(h, index=S.index).rank(ascending=False).values)
RB = np.array(ranks_both)
S["rank_lo90"] = np.percentile(RB, 5, axis=0); S["rank_hi90"] = np.percentile(RB, 95, axis=0)
S["p_top20"] = (RB <= round(0.2 * len(S))).mean(axis=0)
R["mc_draws"] = N
R["median_rank_interval_width"] = float(np.median(S.rank_hi90 - S.rank_lo90))
R["n_robust_top20"] = int((S.p_top20 >= 0.9).sum())
R["n_top20"] = int(round(0.2 * len(S)))
R["rank_stability_spearman"] = {
    "data_error": round(float(np.median([stats.spearmanr(r, S.hvi_rank).statistic for r in ranks_data])), 3),
    "weights": round(float(np.median([stats.spearmanr(r, S.hvi_rank).statistic for r in ranks_w])), 3),
    "both": round(float(np.median([stats.spearmanr(r, S.hvi_rank).statistic for r in ranks_both])), 3)}

# ---------------- exposure facts ----------------
R["lst_range"] = [round(float(S.lst_popw_c.min()), 1), round(float(S.lst_popw_c.max()), 1)]
R["lst_iqr"] = [round(float(S.lst_popw_c.quantile(.25)), 1), round(float(S.lst_popw_c.quantile(.75)), 1)]
R["night_range"] = [round(float(S.modis_myd_night_c.min()), 1), round(float(S.modis_myd_night_c.max()), 1)]
R["lst_year_corr_min"] = round(float(S[yrs].corr(method="spearman").values[np.triu_indices(len(yrs), 1)].min()), 3)
R["lst_year_corr_median"] = round(float(np.median(S[yrs].corr(method="spearman").values[np.triu_indices(len(yrs), 1)])), 3)
R["day_night_spearman"] = round(stats.spearmanr(S.lst_popw_c, S.modis_myd_night_c).statistic, 3)
R["lst_vs_popw_spearman"] = round(stats.spearmanr(S.lst_c, S.lst_popw_c).statistic, 3)
scn = json.load(open(PROC / "landsat_scenes.json"))
R["landsat_scenes_total"] = len(scn); R["landsat_scenes_used"] = sum(s["used"] for s in scn)
R["landsat_dates_used"] = [s["date"] for s in scn if s["used"]]

# ---------------- drivers of surface heat ----------------
m1 = smf.ols("lst_popw_c ~ pct_tree + pct_built", data=S).fit(cov_type="HC3")
R["ols_lst_tree_built"] = {"r2": round(m1.rsquared, 3), "b_tree": round(m1.params.pct_tree, 3),
                           "ci_tree": [round(x, 3) for x in m1.conf_int().loc["pct_tree"]],
                           "b_built": round(m1.params.pct_built, 3),
                           "ci_built": [round(x, 3) for x in m1.conf_int().loc["pct_built"]]}
R["spearman_lst_tree"] = round(stats.spearmanr(S.lst_popw_c, S.pct_tree).statistic, 3)
R["spearman_lst_ndvi"] = round(stats.spearmanr(S.lst_popw_c, S.ndvi).statistic, 3)

# ---------------- equity ----------------
def corr(a, b):
    r = stats.spearmanr(S[a], S[b], nan_policy="omit"); return [round(r.statistic, 3), float(r.pvalue)]
R["eq_spearman"] = {f"{x}~{y}": corr(x, y) for x in ["lst_popw_c", "modis_myd_night_c", "hvi"]
                    for y in ["pct_poverty", "pct_black_nh", "median_income", "pct_tree"]}
m2 = smf.ols("lst_popw_c ~ pct_poverty", data=S).fit(cov_type="HC3")
R["lst_poverty_slope_per10"] = [round(10 * m2.params.pct_poverty, 2), [round(10 * x, 2) for x in m2.conf_int().loc["pct_poverty"]]]
# group-weighted exposure (tract shares x tract population)
for grp, col in [("black_nh", "pct_black_nh"), ("white_nh", "pct_white_nh"), ("hispanic", "pct_hispanic"),
                 ("poverty", "pct_poverty"), ("age65", "pct_65plus")]:
    ok_ = S[col].notna()
    w = (S.pop_acs * S[col] / 100)[ok_]
    R[f"exp_day_{grp}"] = round(float(np.average(S.lst_popw_c[ok_], weights=w)), 2)
    R[f"exp_night_{grp}"] = round(float(np.average(S.modis_myd_night_c[ok_], weights=w)), 2)
R["exp_day_all"] = round(float(np.average(S.lst_popw_c, weights=S.pop_acs)), 2)
R["exp_night_all"] = round(float(np.average(S.modis_myd_night_c, weights=S.pop_acs)), 2)
# HOLC
H = S[S.holc_grade.notna()]
R["holc_n_by_grade"] = H.holc_grade.value_counts().sort_index().to_dict()
R["holc_lst_by_grade"] = H.groupby("holc_grade").lst_popw_c.mean().round(2).to_dict()
R["holc_night_by_grade"] = H.groupby("holc_grade").modis_myd_night_c.mean().round(3).to_dict()
R["holc_tree_by_grade"] = H.groupby("holc_grade").pct_tree.mean().round(1).to_dict()
R["holc_hvi_by_grade"] = H.groupby("holc_grade").hvi.mean().round(3).to_dict()
kw = stats.kruskal(*[H[H.holc_grade == x].lst_popw_c for x in "ABCD" if (H.holc_grade == x).any()])
R["holc_kruskal_lst"] = [round(kw.statistic, 2), float(kw.pvalue)]
m3 = smf.ols("lst_popw_c ~ holc_score", data=S[S.holc_cover >= 0.5]).fit(cov_type="HC3")
R["holc_score_slope"] = [round(m3.params.holc_score, 2), [round(x, 2) for x in m3.conf_int().loc["holc_score"]], float(m3.pvalues.holc_score)]
R["holc_D_minus_A_or_B"] = round(float(H[H.holc_grade == "D"].lst_popw_c.mean() - H[H.holc_grade.isin(["A", "B"])].lst_popw_c.mean()), 2)

# ---------------- spatial structure ----------------
W = Queen.from_dataframe(S, use_index=False, silence_warnings=True); W.transform = "r"
for v in ["hvi", "lst_popw_c", "modis_myd_night_c", "score_sensitivity", "score_low_capacity"]:
    mi = Moran(S[v].values, W, permutations=999)
    R[f"moran_{v}"] = [round(mi.I, 3), float(mi.p_sim)]
lisa = Moran_Local(S.hvi.values, W, permutations=999, seed=1)
lab = np.array(["Not significant"] * len(S), dtype=object)
sig = lisa.p_sim < 0.05
lab[sig & (lisa.q == 1)] = "High-High"; lab[sig & (lisa.q == 3)] = "Low-Low"
lab[sig & (lisa.q == 2)] = "Low-High"; lab[sig & (lisa.q == 4)] = "High-Low"
S["lisa"] = lab
R["lisa_counts"] = pd.Series(lab).value_counts().to_dict()
R["lisa_hh_pop"] = int(S.loc[S.lisa == "High-High", "pop_acs"].sum())

# ---------------- cooling-centre coverage ----------------
hi = S[S.hvi_class == "Very high"]
R["veryhigh_n"] = int(len(hi)); R["veryhigh_pop"] = int(hi.pop_acs.sum())
R["veryhigh_mean_dist_cool"] = round(float(np.average(hi.dist_cool_km, weights=hi.pop_acs)), 2)
R["all_mean_dist_cool"] = round(float(np.average(S.dist_cool_km, weights=S.pop_acs)), 2)
blk = pd.read_csv(PROC / "blocks_lst.csv", dtype={"GEOID": str, "tract": str})
blk = blk[blk.POP100 > 0]
vh = blk[blk.tract.isin(hi.GEOID)]
from scipy.spatial import cKDTree
cool = gpd.read_file(PROC / "cooling_centers.gpkg")
R["n_cooling_centers"] = int(len(cool))

# ---------------- validation vs heat-related EMS (CDC HHI, ZCTA) ----------------
hh = pd.read_excel(RAW / "hhi" / "HHI Data 2024 United States.xlsx", dtype={"ZCTA": str})
hh = hh.replace(-999, np.nan)
bz = pd.read_csv(PROC / "block_zcta.csv", dtype={"GEOID": str, "tract": str, "ZCTA5": str})
bz = bz[bz.POP100 > 0].merge(S[["GEOID", "hvi", "lst_popw_c", "modis_myd_night_c", "score_exposure", "score_sensitivity",
                                 "score_low_capacity", "svi_themes"]].rename(columns={"GEOID": "tract"}), on="tract", how="inner")
def wavg(d):
    w = d.POP100
    return pd.Series({c: np.average(d[c].fillna(d[c].mean()), weights=w) for c in
                      ["hvi", "lst_popw_c", "modis_myd_night_c", "score_exposure", "score_sensitivity", "score_low_capacity", "svi_themes"]}
                     | {"pop_city": w.sum()})
zz = bz.groupby("ZCTA5").apply(wavg, include_groups=False).reset_index()
zz = zz.merge(hh[["ZCTA", "POP", "PR_HRI", "LOW_EMS", "P_NEHD", "OVERALL_RANK"]], left_on="ZCTA5", right_on="ZCTA", how="inner")
zz["city_share"] = zz.pop_city / zz.POP
Zv = zz[(zz.pop_city >= 5000) & (zz.city_share >= 0.5) & (zz.LOW_EMS == 0)].copy()
R["zcta_n"] = int(len(Zv)); R["zcta_list"] = Zv.ZCTA5.tolist()
def sp(a, b):
    r = stats.spearmanr(Zv[a], Zv[b]); return [round(r.statistic, 3), round(float(r.pvalue), 4)]
R["val_spearman"] = {v: sp(v, "PR_HRI") for v in ["hvi", "lst_popw_c", "modis_myd_night_c", "score_exposure",
                                                    "score_sensitivity", "score_low_capacity", "svi_themes", "OVERALL_RANK"]}
# bootstrap CI for HVI ~ EMS
bs = []
for _ in range(2000):
    s = Zv.sample(len(Zv), replace=True, random_state=int(rng.integers(1e9)))
    bs.append(stats.spearmanr(s.hvi, s.PR_HRI).statistic)
R["val_hvi_boot_ci"] = [round(float(np.nanpercentile(bs, 2.5)), 3), round(float(np.nanpercentile(bs, 97.5)), 3)]
# excess harm: EMS percentile not explained by exposure
mx = smf.ols("PR_HRI ~ score_exposure", data=Zv).fit()
Zv["excess_ems"] = mx.resid
R["excess_model_r2"] = round(mx.rsquared, 3)
R["excess_vs"] = {v: sp(v, "excess_ems") for v in ["score_sensitivity", "score_low_capacity", "svi_themes"]}
mfull = smf.ols("PR_HRI ~ score_exposure + score_sensitivity + score_low_capacity", data=Zv).fit()
R["excess_full_r2"] = round(mfull.rsquared, 3); R["excess_full_adj_r2"] = round(mfull.rsquared_adj, 3)
R["excess_full_coef"] = {k: [round(mfull.params[k], 2), [round(x, 2) for x in mfull.conf_int().loc[k]], round(float(mfull.pvalues[k]), 4)] for k in mfull.params.index}
R["excess_exposure_only_adj_r2"] = round(mx.rsquared_adj, 3)
# nested F-test: do sensitivity + capacity add to exposure?
ft = mfull.compare_f_test(mx); R["excess_nested_F"] = [round(ft[0], 2), round(float(ft[1]), 4)]
# leave-one-out prediction check
loo = []
for i in Zv.index:
    tr_ = Zv.drop(i); m_ = smf.ols("PR_HRI ~ score_exposure + score_sensitivity + score_low_capacity", data=tr_).fit()
    m0 = smf.ols("PR_HRI ~ score_exposure", data=tr_).fit()
    loo.append((Zv.loc[i, "PR_HRI"], m_.predict(Zv.loc[[i]]).iloc[0], m0.predict(Zv.loc[[i]]).iloc[0]))
loo = np.array(loo)
R["loo_mae_full"] = round(float(np.abs(loo[:, 0] - loo[:, 1]).mean()), 2)
R["loo_mae_exposure"] = round(float(np.abs(loo[:, 0] - loo[:, 2]).mean()), 2)
Zv.to_csv(PROC / "zcta_validation.csv", index=False)

# ---------------- station climate context ----------------
ac = pd.read_csv(PROC / "acis_daily.csv", parse_dates=["date"])
ac = ac[ac.date.dt.month.isin([6, 7, 8])]
cl = ac.groupby([ac.date.dt.year, "station"]).agg(n_valid=("tmax_f", lambda x: int(x.notna().sum())),
                                                    days90=("tmax_f", lambda x: int((x >= 90).sum())),
                                                    nights75=("tmin_f", lambda x: int((x >= 75).sum())),
                                                    mean_tmin=("tmin_f", "mean")).round(1)
cl = cl[cl.n_valid >= 85]  # keep summers with >= 85 of 92 days reported
R["station_summers"] = {f"{y}_{s}": v for (y, s), v in cl.to_dict("index").items()}
R["paired_days_dmh_bwi"] = None
both = ac.pivot_table(index="date", columns="station", values="tmin_f").dropna()
R["dmh_minus_bwi_tmin_f"] = round(float((both.DMH - both.BWI).mean()), 2)
R["paired_days_dmh_bwi"] = int(len(both))
bothx = ac.pivot_table(index="date", columns="station", values="tmax_f").dropna()
R["dmh_minus_bwi_tmax_f"] = round(float((bothx.DMH - bothx.BWI).mean()), 2)

# ---------------- top tracts ----------------
top = S.sort_values("hvi", ascending=False).head(15)
R["top15"] = top[["GEOID", "NAME", "csa", "hvi", "hvi_rank", "rank_lo90", "rank_hi90", "p_top20", "lst_popw_c", "modis_myd_night_c",
                  "pct_65plus", "pct_poverty", "pl_CHD", "pct_tree", "lisa"]].round(3).to_dict("records")

# ---------------- save ----------------
keep_cols = [c for c in S.columns if c != "geometry"]
out = g[["GEOID", "geometry"]].merge(S[keep_cols], on="GEOID", how="left")
out.to_file(PROC / "baltimore_heat_tracts.gpkg", driver="GPKG")
out.drop(columns="geometry").to_csv(PROC / "baltimore_heat_tracts.csv", index=False)
json.dump(R, open(PROC / "results.json", "w"), indent=1, default=lambda o: o.item() if hasattr(o, "item") else str(o))
json.dump(IND, open(PROC / "indicators.json", "w"), indent=1)
print(json.dumps({k: v for k, v in R.items() if k not in ("top15", "landsat_dates_used", "station_summers")}, indent=1, default=str))
