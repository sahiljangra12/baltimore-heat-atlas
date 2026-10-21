"""Publication figures (PDF + PNG) for the paper and website."""
import json, gzip
import numpy as np, pandas as pd, geopandas as gpd, rasterio
import matplotlib as mpl, matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, ListedColormap, BoundaryNorm
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from common import RAW, PROC, FIG, CRS, GRID

mpl.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8.5, "axes.titlesize": 9.5, "axes.titleweight": "bold",
                     "axes.titlelocation": "left", "axes.spines.top": False, "axes.spines.right": False,
                     "axes.edgecolor": "#8a8780", "axes.labelcolor": "#3d3c39", "xtick.color": "#3d3c39",
                     "ytick.color": "#3d3c39", "savefig.dpi": 300, "figure.dpi": 150, "pdf.fonttype": 42})
INK, INK2, GRID_C, NOD = "#1f1f1e", "#5b5a56", "#e9e7e1", "#dcdad3"
HEAT = LinearSegmentedColormap.from_list("heat", ["#fff1e8", "#fbc4a6", "#f28a5b", "#d65a24", "#9a3a12", "#5c1f07"])
BLUE = LinearSegmentedColormap.from_list("blue", ["#eaf2fd", "#b7d3f6", "#6da7ec", "#2a78d6", "#184f95", "#0d366b"])
NIGHT = LinearSegmentedColormap.from_list("night", ["#f3effc", "#cfc5f3", "#9c8ae6", "#6250d6", "#3d2f99", "#211a5c"])
GREEN = LinearSegmentedColormap.from_list("green", ["#f1f8ec", "#c8e6b5", "#8ccb6f", "#4a9d34", "#256d18", "#123d0b"])

R = json.load(open(PROC / "results.json"))
g = gpd.read_file(PROC / "baltimore_heat_tracts.gpkg").to_crs(CRS)
city = g.dissolve()
cool = gpd.read_file(PROC / "cooling_centers.gpkg").to_crs(CRS)
S = g[g.hvi.notna()]
ext = city.total_bounds


def base(ax):
    ax.set_axis_off(); ax.set_aspect("equal")
    ax.set_xlim(ext[0] - 300, ext[2] + 300); ax.set_ylim(ext[1] - 300, ext[3] + 300)


def scalebar(ax, km=2):
    x0, y0 = ext[0] + 200, ext[1] + 200
    ax.plot([x0, x0 + km * 1000], [y0, y0], color=INK, lw=2, solid_capstyle="butt")
    ax.text(x0 + km * 500, y0 + 250, f"{km} km", ha="center", fontsize=7, color=INK)


def cbar(fig, ax, mappable, label):
    cb = fig.colorbar(mappable, ax=ax, orientation="horizontal", fraction=0.04, pad=0.01, shrink=0.75)
    cb.set_label(label, fontsize=7.5, color=INK2); cb.ax.tick_params(labelsize=7); cb.outline.set_visible(False)
    return cb


def save(fig, name):
    fig.savefig(FIG / f"{name}.pdf", bbox_inches="tight"); fig.savefig(FIG / f"{name}.png", bbox_inches="tight")
    plt.close(fig)


# ---- Fig 1: 30 m daytime LST raster + night LST ----
with rasterio.open(PROC / "lst_median.tif") as s:
    lst = s.read(1); b = s.bounds
with rasterio.open(PROC / "modis_myd_night.tif") as s:
    night = s.read(1)
from rasterio.features import geometry_mask
with rasterio.open(PROC / "lst_median.tif") as s:
    outside = geometry_mask(city.geometry, out_shape=lst.shape, transform=s.transform)
lst_m = np.where(outside, np.nan, lst); night_m = np.where(outside, np.nan, night)
fig, axs = plt.subplots(1, 2, figsize=(7.2, 4.6))
for ax, arr, cmap, lab, ttl in [(axs[0], lst_m, HEAT, "Land surface temperature (°C)", "a  Daytime surface heat, late morning"),
                                (axs[1], night_m, NIGHT, "Land surface temperature (°C)", "b  Night-time surface heat, ~1:30 am solar time")]:
    lo, hi = np.nanpercentile(arr, [5, 99.5])
    im = ax.imshow(arr, extent=(b.left, b.right, b.bottom, b.top), cmap=cmap, vmin=lo, vmax=hi, interpolation="nearest")
    city.boundary.plot(ax=ax, color=INK, lw=0.6)
    base(ax); ax.set_title(ttl); cbar(fig, ax, im, lab)
cool.plot(ax=axs[0], marker="^", color="white", edgecolor=INK, markersize=14, lw=0.6, zorder=5)
axs[0].legend(handles=[Line2D([], [], marker="^", ls="", mfc="white", mec=INK, label="Cooling centre (2026)")],
              loc="upper left", fontsize=7, frameon=False)
scalebar(axs[0])
fig.text(0.01, -0.01, f"Landsat 8/9 Collection 2 surface temperature, per-pixel median of {R['landsat_scenes_used']} cloud-screened summer scenes (Jun–Aug 2019–2025), 30 m.\n"
         "MODIS Aqua MYD11A2 8-day night LST, per-cell median of summer composites 2019–2025, 1 km.", fontsize=6.5, color=INK2)
save(fig, "fig1_lst_day_night")

# ---- Fig 2: components ----
fig, axs = plt.subplots(1, 3, figsize=(7.4, 3.6))
for ax, (c, cmap, t) in zip(axs, [("score_exposure", HEAT, "a  Exposure"), ("score_sensitivity", BLUE, "b  Sensitivity"),
                                  ("score_low_capacity", BLUE, "c  Low adaptive capacity")]):
    g.plot(ax=ax, color=NOD, edgecolor="white", lw=0.2)
    S.plot(ax=ax, column=c, cmap=cmap, vmin=0, vmax=1, edgecolor="white", lw=0.2)
    base(ax); ax.set_title(t)
    sm = mpl.cm.ScalarMappable(cmap=cmap, norm=mpl.colors.Normalize(0, 1)); cbar(fig, ax, sm, "Mean percentile (0–1)")
save(fig, "fig2_components")

# ---- Fig 3: HVI + probability of being top 20% ----
fig, axs = plt.subplots(1, 2, figsize=(7.2, 4.6))
cls = ["Very low", "Low", "Moderate", "High", "Very high"]
cols5 = ["#fde4d3", "#f6b48f", "#ec7f4e", "#c8501c", "#7a2b0b"]
g.plot(ax=axs[0], color=NOD, edgecolor="white", lw=0.2)
for c, col in zip(cls, cols5):
    S[S.hvi_class == c].plot(ax=axs[0], color=col, edgecolor="white", lw=0.2)
axs[0].legend(handles=[Patch(color=c, label=l) for c, l in zip(cols5, cls)][::-1], title="HVI quintile", fontsize=7,
              title_fontsize=7.5, loc="lower right", frameon=True, facecolor="white", edgecolor="none", framealpha=0.95)
base(axs[0]); axs[0].set_title("a  Heat Vulnerability Index"); scalebar(axs[0])
g.plot(ax=axs[1], color=NOD, edgecolor="white", lw=0.2)
S.plot(ax=axs[1], column="p_top20", cmap=HEAT, vmin=0, vmax=1, edgecolor="white", lw=0.2)
robust = S[S.p_top20 >= 0.9]
robust.boundary.plot(ax=axs[1], color=INK, lw=0.9)
base(axs[1]); axs[1].set_title("b  Probability of ranking in the top 20%")
sm = mpl.cm.ScalarMappable(cmap=HEAT, norm=mpl.colors.Normalize(0, 1)); cbar(fig, axs[1], sm, "Share of 1,000 Monte Carlo draws")
axs[1].legend(handles=[Patch(fc="none", ec=INK, lw=0.9, label="≥ 90% of draws")], loc="lower right", fontsize=7, frameon=True, facecolor="white", edgecolor="none", framealpha=0.95)
save(fig, "fig3_hvi_uncertainty")

# ---- Fig 4: LISA ----
fig, ax = plt.subplots(figsize=(6.0, 4.6))
lc = {"High-High": "#b83c0f", "Low-Low": "#2a78d6", "High-Low": "#f2a37c", "Low-High": "#9ec5f4", "Not significant": "#ecebe6"}
g.plot(ax=ax, color=NOD, edgecolor="white", lw=0.2)
for k, c in lc.items():
    sub = S[S.lisa == k]
    if len(sub): sub.plot(ax=ax, color=c, edgecolor="white", lw=0.2)
ax.legend(handles=[Patch(color=c, label=f"{k} ({R['lisa_counts'].get(k, 0)})") for k, c in lc.items()], fontsize=7,
          loc="center left", bbox_to_anchor=(1.0, 0.5), frameon=False, title="Local Moran's I, p < 0.05", title_fontsize=7.5)
base(ax); ax.set_title("Clusters of heat vulnerability")
save(fig, "fig4_lisa")

# ---- Fig 5: HOLC ----
H = S[S.holc_grade.notna()]
fig, axs = plt.subplots(1, 2, figsize=(7.0, 3.0))
holc_c = {"A": "#76a865", "B": "#7cb5bd", "C": "#ffff00", "D": "#d9838d"}
for ax, v, lab in [(axs[0], "lst_popw_c", "Daytime surface temperature (°C)"), (axs[1], "pct_tree", "Tree cover (%)")]:
    data = [H[H.holc_grade == k][v].values for k in "ABCD"]
    bp = ax.boxplot(data, widths=0.55, patch_artist=True, medianprops=dict(color=INK, lw=1.2), whiskerprops=dict(color=INK2),
                    capprops=dict(color=INK2), flierprops=dict(marker="o", ms=2.5, mfc=INK2, mec="none"))
    for p, k in zip(bp["boxes"], "ABCD"):
        p.set_facecolor(holc_c[k]); p.set_edgecolor(INK2); p.set_alpha(0.85)
    ax.set_xticks(range(1, 5), [f"{k}\n(n={len(d)})" for k, d in zip("ABCD", data)])
    ax.set_ylabel(lab); ax.grid(axis="y", color=GRID_C, lw=0.6); ax.set_axisbelow(True)
axs[0].set_title("a  Surface heat by 1937 HOLC grade"); axs[1].set_title("b  Tree cover by 1937 HOLC grade")
fig.text(0.01, -0.06, "Tracts assigned the HOLC grade covering most of their area (≥ 50% of tract graded). HOLC colours follow the original maps.",
         fontsize=6.5, color=INK2)
fig.tight_layout(); save(fig, "fig5_holc")

# ---- Fig 6: drivers scatter ----
fig, axs = plt.subplots(1, 2, figsize=(7.0, 3.0))
ax = axs[0]
ax.scatter(S.pct_tree, S.lst_popw_c, s=12, color="#d65a24", edgecolor="white", lw=0.5)
z = np.polyfit(S.pct_tree, S.lst_popw_c, 1); xs = np.linspace(0, S.pct_tree.max(), 50)
ax.plot(xs, np.polyval(z, xs), color=INK, lw=1.2)
ax.set_xlabel("Tree cover (%, ESA WorldCover 2021)"); ax.set_ylabel("Daytime surface temperature (°C)")
ax.set_title(f"a  Trees and surface heat (ρ = {R['spearman_lst_tree']:.2f})"); ax.grid(color=GRID_C, lw=0.6); ax.set_axisbelow(True)
ax = axs[1]
ax.scatter(S.lst_popw_c, S.modis_myd_night_c, s=12, color="#6250d6", edgecolor="white", lw=0.5)
ax.set_xlabel("Daytime surface temperature (°C, Landsat)"); ax.set_ylabel("Night surface temperature (°C, MODIS)")
ax.set_title(f"b  Day vs night heat (ρ = {R['day_night_spearman']:.2f})"); ax.grid(color=GRID_C, lw=0.6); ax.set_axisbelow(True)
fig.tight_layout(); save(fig, "fig6_drivers")

# ---- Fig 7: validation against EMS ----
Z = pd.read_csv(PROC / "zcta_validation.csv", dtype={"ZCTA5": str})
fig, axs = plt.subplots(1, 3, figsize=(7.4, 2.7), sharey=True)
for ax, v, lab, c in [(axs[0], "hvi", "Heat Vulnerability Index", "#d65a24"), (axs[1], "score_exposure", "Exposure score", "#9a3a12"),
                      (axs[2], "score_low_capacity", "Low-capacity score", "#2a78d6")]:
    ax.scatter(Z[v], Z.PR_HRI, s=18, color=c, edgecolor="white", lw=0.5)
    zz = np.polyfit(Z[v], Z.PR_HRI, 1); xs = np.linspace(Z[v].min(), Z[v].max(), 20)
    ax.plot(xs, np.polyval(zz, xs), color=INK, lw=1)
    rho, p = R["val_spearman"][v]
    ptxt = "p < 0.001" if p < 0.001 else f"p = {p:.3f}"
    ax.set_title(f"ρ = {rho:.2f} ({ptxt})", fontsize=8.5)
    ax.set_xlabel(f"{lab}\n(ZCTA, population-weighted)"); ax.grid(color=GRID_C, lw=0.6); ax.set_axisbelow(True)
axs[0].set_ylabel("Heat-related EMS activations\n(national percentile, CDC HHI)")
fig.suptitle(f"Validation against heat-related EMS activations, {R['zcta_n']} Baltimore ZIP code areas", x=0.01, ha="left",
             fontsize=9.5, fontweight="bold", y=1.04)
fig.tight_layout(); save(fig, "fig7_validation")

# ---- Fig 8: station climate ----
ss = R["station_summers"]
yrs = list(range(2019, 2026))
fig, axs = plt.subplots(1, 2, figsize=(7.0, 2.7))
for ax, k, t in [(axs[0], "days90", "a  Days ≥ 90°F (32.2°C)"), (axs[1], "nights75", "b  Nights ≥ 75°F (23.9°C)")]:
    for j, (st, c, lab) in enumerate([("BWI", "#8a8780", "BWI Airport"), ("DMH", "#d65a24", "Inner Harbor (DMH)")]):
        vals = [ss.get(f"{y}_{st}", {}).get(k, np.nan) for y in yrs]
        ax.bar(np.arange(len(yrs)) + (j - 0.5) * 0.38, vals, width=0.36, color=c, label=lab)
    ax.set_xticks(range(len(yrs)), yrs); ax.set_title(t); ax.grid(axis="y", color=GRID_C, lw=0.6); ax.set_axisbelow(True)
axs[0].legend(fontsize=7, frameon=False)
fig.text(0.01, -0.05, "June–August. Missing bars = station summer with fewer than 85 of 92 days reported. Source: NOAA RCC-ACIS.",
         fontsize=6.5, color=INK2)
fig.tight_layout(); save(fig, "fig8_station_climate")

# ---- Fig 9: rank uncertainty for top 40 ----
T = S.sort_values("hvi_rank").head(40)
fig, ax = plt.subplots(figsize=(7.0, 3.2))
x = np.arange(len(T))
ax.vlines(x, T.rank_lo90, T.rank_hi90, color="#f28a5b", lw=2.2)
ax.scatter(x, T.hvi_rank, s=14, color="#7a2b0b", zorder=3)
ax.axhline(R["n_top20"], color=INK2, lw=0.8, ls="--"); ax.text(len(T) - 0.5, R["n_top20"] + 2, "top-20% cut-off", ha="right", fontsize=7, color=INK2)
ax.set_xticks(x, [n.replace("Census Tract ", "") for n in T.NAME], rotation=90, fontsize=6.5)
ax.invert_yaxis(); ax.set_ylabel("HVI rank (1 = most vulnerable)"); ax.set_xlabel("Census tract")
ax.set_title("Rank and 90% uncertainty interval, 40 highest-ranked tracts"); ax.grid(axis="y", color=GRID_C, lw=0.6); ax.set_axisbelow(True)
fig.tight_layout(); save(fig, "fig9_rank_uncertainty")
print("figures saved:", sorted(p.name for p in FIG.glob("*.png")))
