from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
PROC = ROOT / "data" / "processed"
FIG = ROOT / "outputs" / "figures"
for d in (PROC, FIG):
    d.mkdir(parents=True, exist_ok=True)
CRS = "EPSG:32618"
GRID = dict(x0=351435, y0=4360395, res=30, w=589, h=715)
