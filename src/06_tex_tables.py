"""Write LaTeX macros (every number quoted in the paper) and tables straight from the results,
so no figure in the text is typed by hand."""
import json
import numpy as np, pandas as pd
from common import PROC, ROOT

R = json.load(open(PROC / "results.json"))
IND = json.load(open(PROC / "indicators.json"))
d = pd.read_csv(PROC / "baltimore_heat_tracts.csv", dtype={"GEOID": str})
S = d[d.hvi.notna()].copy()
TEX = ROOT / "paper"; TEX.mkdir(exist_ok=True)
mac = {}


def f(x, n=1):
    return f"{x:,.{n}f}"


def p(x):
    return r"$<$0.001" if x < 0.001 else f"{x:.3f}"


mac.update({
    "NTracts": R["n_tracts"], "NScored": R["n_scored"], "PopTwenty": f"{R['pop2020_total']:,}", "PopACS": f"{R['pop_acs_total']:,}",
    "NScenes": R["landsat_scenes_total"], "NScenesUsed": R["landsat_scenes_used"],
    "LSTmin": f(R["lst_range"][0]), "LSTmax": f(R["lst_range"][1]), "LSTqOne": f(R["lst_iqr"][0]), "LSTqThree": f(R["lst_iqr"][1]),
    "LSTspread": f(R["lst_range"][1] - R["lst_range"][0]),
    "NightMin": f(R["night_range"][0]), "NightMax": f(R["night_range"][1]),
    "YearCorrMin": f(R["lst_year_corr_min"], 2), "YearCorrMed": f(R["lst_year_corr_median"], 2),
    "DayNightRho": f(R["day_night_spearman"], 2), "PopwRho": f(R["lst_vs_popw_spearman"], 2),
    "TreeRho": f(R["spearman_lst_tree"], 2), "NDVIRho": f(R["spearman_lst_ndvi"], 2),
    "OLSRtwo": f(R["ols_lst_tree_built"]["r2"], 2), "BTree": f(R["ols_lst_tree_built"]["b_tree"], 3),
    "BTreeLo": f(R["ols_lst_tree_built"]["ci_tree"][0], 3), "BTreeHi": f(R["ols_lst_tree_built"]["ci_tree"][1], 3),
    "BTreeTen": f(-10 * R["ols_lst_tree_built"]["b_tree"], 2),
    "BBuilt": f(R["ols_lst_tree_built"]["b_built"], 3), "BBuiltLo": f(R["ols_lst_tree_built"]["ci_built"][0], 3),
    "BBuiltHi": f(R["ols_lst_tree_built"]["ci_built"][1], 3),
    "CompES": f(R["component_spearman"]["exposure~sensitivity"], 2), "CompEC": f(R["component_spearman"]["exposure~low_capacity"], 2),
    "CompSC": f(R["component_spearman"]["sensitivity~low_capacity"], 2),
    "PCAk": R["pca_components_kaiser"], "PCAvar": f(100 * R["pca_var_explained"], 0),
    "AltPCA": f(R["alt_index_spearman"]["hvi_pca"], 2), "AltGeo": f(R["alt_index_spearman"]["hvi_geo"], 2),
    "AltFlat": f(R["alt_index_spearman"]["hvi_flat"], 2), "SVIrho": f(R["svi_vs_hvi_spearman"], 2),
    "MCdraws": f"{R['mc_draws']:,}", "RankWidth": f(R["median_rank_interval_width"], 0),
    "NRobust": R["n_robust_top20"], "NTopTwenty": R["n_top20"],
    "StabData": f(R["rank_stability_spearman"]["data_error"], 2), "StabW": f(R["rank_stability_spearman"]["weights"], 2),
    "StabBoth": f(R["rank_stability_spearman"]["both"], 2),
    "PovSlope": f(R["lst_poverty_slope_per10"][0], 2), "PovSlopeLo": f(R["lst_poverty_slope_per10"][1][0], 2),
    "PovSlopeHi": f(R["lst_poverty_slope_per10"][1][1], 2),
    "RhoLSTPov": f(R["eq_spearman"]["lst_popw_c~pct_poverty"][0], 2), "RhoLSTBlack": f(R["eq_spearman"]["lst_popw_c~pct_black_nh"][0], 2),
    "PLSTBlack": f(R["eq_spearman"]["lst_popw_c~pct_black_nh"][1], 2),
    "RhoNightPov": f(R["eq_spearman"]["modis_myd_night_c~pct_poverty"][0], 2), "RhoNightBlack": f(R["eq_spearman"]["modis_myd_night_c~pct_black_nh"][0], 2),
    "RhoHVIPov": f(R["eq_spearman"]["hvi~pct_poverty"][0], 2), "RhoHVIBlack": f(R["eq_spearman"]["hvi~pct_black_nh"][0], 2),
    "RhoHVIInc": f(R["eq_spearman"]["hvi~median_income"][0], 2),
    "ExpDayAll": f(R["exp_day_all"], 2), "ExpDayBlack": f(R["exp_day_black_nh"], 2), "ExpDayWhite": f(R["exp_day_white_nh"], 2),
    "ExpDayHisp": f(R["exp_day_hispanic"], 2), "ExpDayPov": f(R["exp_day_poverty"], 2), "ExpDayOld": f(R["exp_day_age65"], 2),
    "ExpNightAll": f(R["exp_night_all"], 2), "ExpNightBlack": f(R["exp_night_black_nh"], 2), "ExpNightWhite": f(R["exp_night_white_nh"], 2),
    "ExpNightHisp": f(R["exp_night_hispanic"], 2), "ExpNightPov": f(R["exp_night_poverty"], 2), "ExpNightOld": f(R["exp_night_age65"], 2),
    "HolcA": f(R["holc_lst_by_grade"]["A"], 1), "HolcB": f(R["holc_lst_by_grade"]["B"], 1), "HolcC": f(R["holc_lst_by_grade"]["C"], 1),
    "HolcD": f(R["holc_lst_by_grade"]["D"], 1), "HolcDminusAB": f(R["holc_D_minus_A_or_B"], 1),
    "HolcTreeA": f(R["holc_tree_by_grade"]["A"], 0), "HolcTreeB": f(R["holc_tree_by_grade"]["B"], 0),
    "HolcTreeC": f(R["holc_tree_by_grade"]["C"], 0), "HolcTreeD": f(R["holc_tree_by_grade"]["D"], 0),
    "HolcNA": R["holc_n_by_grade"]["A"], "HolcNB": R["holc_n_by_grade"]["B"], "HolcNC": R["holc_n_by_grade"]["C"], "HolcND": R["holc_n_by_grade"]["D"],
    "HolcKW": f(R["holc_kruskal_lst"][0], 1), "HolcSlope": f(R["holc_score_slope"][0], 2),
    "HolcSlopeLo": f(R["holc_score_slope"][1][0], 2), "HolcSlopeHi": f(R["holc_score_slope"][1][1], 2),
    "HolcNightD": f(R["holc_night_by_grade"]["D"], 1), "HolcNightB": f(R["holc_night_by_grade"]["B"], 1),
    "MoranHVI": f(R["moran_hvi"][0], 2), "MoranLST": f(R["moran_lst_popw_c"][0], 2), "MoranNight": f(R["moran_modis_myd_night_c"][0], 2),
    "MoranSens": f(R["moran_score_sensitivity"][0], 2), "MoranCap": f(R["moran_score_low_capacity"][0], 2),
    "LisaHH": R["lisa_counts"].get("High-High", 0), "LisaLL": R["lisa_counts"].get("Low-Low", 0),
    "LisaHL": R["lisa_counts"].get("High-Low", 0), "LisaHHpop": f"{R['lisa_hh_pop']:,}",
    "VHn": R["veryhigh_n"], "VHpop": f"{R['veryhigh_pop']:,}", "VHdist": f(R["veryhigh_mean_dist_cool"], 2),
    "AllDist": f(R["all_mean_dist_cool"], 2), "NCool": R["n_cooling_centers"],
    "NZcta": R["zcta_n"],
    "ValHVI": f(R["val_spearman"]["hvi"][0], 2), "ValHVIp": p(R["val_spearman"]["hvi"][1]),
    "ValHVIlo": f(R["val_hvi_boot_ci"][0], 2), "ValHVIhi": f(R["val_hvi_boot_ci"][1], 2),
    "ValLST": f(R["val_spearman"]["lst_popw_c"][0], 2), "ValLSTp": p(R["val_spearman"]["lst_popw_c"][1]),
    "ValNight": f(R["val_spearman"]["modis_myd_night_c"][0], 2), "ValNightp": p(R["val_spearman"]["modis_myd_night_c"][1]),
    "ValExp": f(R["val_spearman"]["score_exposure"][0], 2), "ValExpp": p(R["val_spearman"]["score_exposure"][1]),
    "ValSens": f(R["val_spearman"]["score_sensitivity"][0], 2), "ValSensp": p(R["val_spearman"]["score_sensitivity"][1]),
    "ValCap": f(R["val_spearman"]["score_low_capacity"][0], 2), "ValCapp": p(R["val_spearman"]["score_low_capacity"][1]),
    "ValSVI": f(R["val_spearman"]["svi_themes"][0], 2), "ValSVIp": p(R["val_spearman"]["svi_themes"][1]),
    "ValHHI": f(R["val_spearman"]["OVERALL_RANK"][0], 2), "ValHHIp": p(R["val_spearman"]["OVERALL_RANK"][1]),
    "ExcessRtwo": f(R["excess_model_r2"], 2), "ExcessAdjExp": f(R["excess_exposure_only_adj_r2"], 2),
    "FullRtwo": f(R["excess_full_r2"], 2), "FullAdj": f(R["excess_full_adj_r2"], 2),
    "NestedF": f(R["excess_nested_F"][0], 1), "NestedFp": p(R["excess_nested_F"][1]),
    "LooFull": f(R["loo_mae_full"], 1), "LooExp": f(R["loo_mae_exposure"], 1),
    "CapCoef": f(R["excess_full_coef"]["score_low_capacity"][0], 0), "CapCoefLo": f(R["excess_full_coef"]["score_low_capacity"][1][0], 0),
    "CapCoefHi": f(R["excess_full_coef"]["score_low_capacity"][1][1], 0),
    "SensCoef": f(R["excess_full_coef"]["score_sensitivity"][0], 0), "SensCoefLo": f(R["excess_full_coef"]["score_sensitivity"][1][0], 0),
    "SensCoefHi": f(R["excess_full_coef"]["score_sensitivity"][1][1], 0),
    "DMHtmin": f(R["dmh_minus_bwi_tmin_f"], 1), "DMHtmax": f(R["dmh_minus_bwi_tmax_f"], 1),
    "DMHtminC": f(R["dmh_minus_bwi_tmin_f"] * 5 / 9, 1), "DMHtmaxC": f(R["dmh_minus_bwi_tmax_f"] * 5 / 9, 1),
    "PairedDays": R["paired_days_dmh_bwi"],
})
ss = R["station_summers"]
mac["BWInightsTwentyFive"] = ss["2025_BWI"]["nights75"]; mac["DMHnightsTwentyFive"] = ss["2025_DMH"]["nights75"]
mac["BWIdaysTwentyFive"] = ss["2025_BWI"]["days90"]; mac["DMHdaysTwentyFive"] = ss["2025_DMH"]["days90"]
# very-high vs rest
vh = S[S.hvi_class == "Very high"]; rest = S[S.hvi_class != "Very high"]
def wa(df, c):
    m = df[c].notna(); return np.average(df.loc[m, c], weights=df.loc[m, "pop_acs"])
for c, k, n in [("lst_popw_c", "LST", 1), ("modis_myd_night_c", "Night", 1), ("pct_poverty", "Pov", 0), ("pl_CHD", "CHD", 1),
                ("pct_tree", "Tree", 0), ("pct_black_nh", "Black", 0), ("pct_no_vehicle", "NoVeh", 0), ("pl_SHUTUTILITY", "Shut", 0)]:
    mac[f"VH{k}"] = f(wa(vh, c), n); mac[f"Rest{k}"] = f(wa(rest, c), n)
top = S.sort_values("hvi", ascending=False)
csas = top.head(R["veryhigh_n"]).csa.value_counts()
mac["TopCSAone"] = csas.index[0]; mac["TopCSAtwo"] = csas.index[1]; mac["TopCSAthree"] = top.csa.iloc[0]
hq = S[S.lst_popw_c >= S.lst_popw_c.quantile(.75)]
mac["NHotQ"] = len(hq); mac["NHotWhite"] = int((hq.pct_white_nh > 50).sum())
with open(TEX / "macros.tex", "w") as fh:
    for k, v in mac.items():
        fh.write(f"\\newcommand{{\\{k}}}{{{str(v).replace('&', '\\&')}}}\n")

# ---- Table: indicators with descriptive statistics ----
lab = {}
for comp, items in IND.items():
    for c, l in items.items(): lab[c] = (comp, l)
rows = []
for c, (comp, l) in lab.items():
    v = S[c]
    rows.append((comp.replace("low_capacity", "Low capacity").capitalize() if comp != "low_capacity" else "Low capacity",
                 l.replace("%", r"\%").replace("&", r"\&"), f(v.median(), 1), f"{f(v.quantile(.25),1)}--{f(v.quantile(.75),1)}",
                 f(v.min(), 1), f(v.max(), 1)))
src = {"lst_popw_c": "Landsat 8/9 C2 L2, 2019--25", "modis_myd_night_c": "MODIS MYD11A2 v061, 2019--25"}
with open(TEX / "tab_indicators.tex", "w") as fh:
    last = None
    for comp, l, med, iqr, mn, mx in rows:
        if comp != last and last is not None: fh.write("\\midrule\n")
        fh.write(f"{comp if comp != last else ''} & {l} & {med} & {iqr} & {mn} & {mx} \\\\\n"); last = comp

# ---- Table: top 15 tracts ----
with open(TEX / "tab_top.tex", "w") as fh:
    for _, r in top.head(15).iterrows():
        fh.write(f"{int(r.hvi_rank)} & {r.NAME.replace('Census Tract ', '')} & {r.csa.replace('&', '\\&')} & {r.hvi:.3f} & "
                 f"{int(r.rank_lo90)}--{int(np.ceil(r.rank_hi90))} & {r.p_top20:.2f} & {r.lst_popw_c:.1f} & {r.modis_myd_night_c:.1f} & "
                 f"{r.pct_poverty:.0f} & {r.pl_CHD:.1f} & {r.pct_tree:.0f} \\\\\n")

# ---- Table: validation ----
names = {"hvi": "Heat Vulnerability Index (this study)", "score_exposure": "\\quad Exposure component",
         "score_sensitivity": "\\quad Sensitivity component", "score_low_capacity": "\\quad Low-capacity component",
         "lst_popw_c": "Daytime LST alone", "modis_myd_night_c": "Night-time LST alone",
         "svi_themes": "CDC/ATSDR SVI 2022", "OVERALL_RANK": "CDC Heat \\& Health Index 2024 (overall)"}
with open(TEX / "tab_validation.tex", "w") as fh:
    for k, n in names.items():
        rho, pv = R["val_spearman"][k]
        fh.write(f"{n} & {rho:.2f} & {p(pv)} \\\\\n")
print(len(mac), "macros written")
