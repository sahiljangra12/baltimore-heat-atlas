#!/bin/bash
# inline the auto-generated tables (\input inside tabular breaks \bottomrule), then compile
set -e
cd "$(dirname "$0")"
python3 - <<'PY'
import re
t = open("main.tex").read()
t = re.sub(r"\\input\{(tab_[a-z]+\.tex)\}", lambda m: open(m.group(1)).read().rstrip("\n"), t)
open("paper_build.tex", "w").write(t)
PY
latexmk -pdf -interaction=nonstopmode -quiet paper_build.tex >/dev/null 2>&1 || true
cp paper_build.pdf Baltimore_Heat_Vulnerability_Sahil_2026.pdf
