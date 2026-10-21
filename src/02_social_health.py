"""Build tract-level social (ACS 2020-2024), health (CDC PLACES 2025) and SVI 2022 tables.

ACS shares carry margins of error (MOE, 90%) using the Census Bureau's formulas for
derived proportions (ACS General Handbook, ch. 8): MOE_sum = sqrt(sum MOE_i^2);
MOE_p = sqrt(MOE_num^2 - p^2 * MOE_den^2) / den  (ratio formula used when the radicand is negative).
"""
import json
import numpy as np
import pandas as pd
import geopandas as gpd
from common import RAW, PROC

J = json.load(open(RAW / "bhs_acs2024_tracts.json"))
est, moe = {}, {}
for grp in J.values():
    for geo, tabs in grp["data"].items():
        g = geo.replace("14000US", "")
        for t, v in tabs.items():
            est.setdefault(g, {}).update(v["estimate"])
            moe.setdefault(g, {}).update(v["error"])
E = pd.DataFrame(est).T.astype(float)
M = pd.DataFrame(moe).T.astype(float)
E.index.name = M.index.name = "GEOID"


def cols(prefix, nums):
    return [f"{prefix}{n:03d}" for n in nums]


def share(num_cols, den_col):
    num = E[num_cols].sum(axis=1)
    mnum = np.sqrt((M[num_cols] ** 2).sum(axis=1))
    den = E[den_col]; mden = M[den_col]
    p = num / den
    rad = mnum ** 2 - p ** 2 * mden ** 2
    rad = np.where(rad < 0, mnum ** 2 + p ** 2 * mden ** 2, rad)
    mp = np.sqrt(rad) / den
    p = p.where(den > 0); mp = pd.Series(mp, index=E.index).where(den > 0)
    return 100 * p, 100 * mp


acs = pd.DataFrame(index=E.index)
acs["pop_acs"] = E["B01003001"]
spec = {
    "pct_65plus": (cols("B01001", list(range(20, 26)) + list(range(44, 50))), "B01001001"),
    "pct_under5": (cols("B01001", [3, 27]), "B01001001"),
    "pct_poverty": (["B17001002"], "B17001001"),
    "pct_no_vehicle": (["B25044003", "B25044010"], "B25044001"),
    "pct_65_alone": (["B11007003"], "B11007001"),
    "pct_black_nh": (["B03002004"], "B03002001"),
    "pct_white_nh": (["B03002003"], "B03002001"),
    "pct_hispanic": (["B03002012"], "B03002001"),
    "pct_pre1940": (["B25034011"], "B25034001"),
    "pct_renter": (["B25003003"], "B25003001"),
    "pct_no_hs": (cols("B15003", range(2, 17)), "B15003001"),
    "pct_disability": (cols("B18101", [4, 7, 10, 13, 16, 19, 23, 26, 29, 32, 35, 38]), "B18101001"),
    "pct_unemployed": (["B23025005"], "B23025003"),
    "pct_limited_english": (cols("C16002", [4, 7, 10, 13]), "C16002001"),
    "pct_rent_burden50": (["B25070010"], "B25070001"),
    "pct_no_internet": (["B28002013"], "B28002001"),
    "pct_uninsured": (cols("B27010", [17, 33, 50, 66]), "B27010001"),
    "pct_crowded": (cols("B25014", [5, 6, 7, 11, 12, 13]), "B25014001"),
}
for name, (num, den) in spec.items():
    acs[name], acs[name + "_moe"] = share(num, den)
acs["median_income"] = E["B19013001"]
acs["median_income_moe"] = M["B19013001"]
acs = acs.replace([np.inf, -np.inf], np.nan)

# ---- CDC PLACES 2025 release (BRFSS 2023 for most measures) ----
pl = pd.read_csv(RAW / "bhs_places_2025.csv", dtype={"GEOID": str})
wide = pl.pivot(index="GEOID", columns="measure", values="value")
lo = pl.pivot(index="GEOID", columns="measure", values="lo")
hi = pl.pivot(index="GEOID", columns="measure", values="hi")
keep = ["CHD", "COPD", "CASTHMA", "DIABETES", "KIDNEY", "BPHIGH", "OBESITY", "STROKE",
        "MHLTH", "DISABILITY", "MOBILITY", "COGNITION", "INDEPLIVE", "SELFCARE", "LONELINESS",
        "EMOTIONSPT", "LACKTRPT", "HOUSINSECU", "FOODINSECU", "SHUTUTILITY", "GHLTH", "CSMOKING", "ACCESS2"]
keep = [k for k in keep if k in wide.columns]
places = wide[keep].add_prefix("pl_")
for k in keep:
    places[f"pl_{k}_se"] = (hi[k] - lo[k]) / (2 * 1.96)

# ---- SVI 2022 (Maryland, state-ranked) ----
svi = pd.read_csv(RAW / "bhs_svi2022_maryland.csv", dtype={"FIPS": str})
svi = svi[svi["FIPS"].str.startswith("24510")].set_index("FIPS")
svi = svi[["RPL_THEMES", "RPL_THEME1", "RPL_THEME2", "RPL_THEME3", "RPL_THEME4"]]
svi = svi.where(svi >= 0).rename(columns=lambda c: "svi_" + c.lower().replace("rpl_", ""))
svi.index.name = "GEOID"

# ---- Tracts (2020) ----
tr = gpd.read_file(RAW / "bhs_tracts2020.geojson")[["GEOID", "NAME", "AREALAND", "POP100", "geometry"]]
tr = tr.rename(columns={"POP100": "pop2020", "AREALAND": "aland_m2"})
tr = tr[tr["aland_m2"] > 0]
df = tr.merge(acs.reset_index(), on="GEOID", how="left") \
       .merge(places.reset_index(), on="GEOID", how="left") \
       .merge(svi.reset_index(), on="GEOID", how="left")
df.to_file(PROC / "tracts_social_health.gpkg", driver="GPKG")
print(df.shape, "tracts; ACS missing:", df["pop_acs"].isna().sum(), "PLACES missing:", df["pl_CHD"].isna().sum(),
      "SVI missing:", df["svi_themes"].isna().sum())
print(df[["pct_65plus", "pct_poverty", "pl_CHD", "svi_themes"]].describe().round(2))
