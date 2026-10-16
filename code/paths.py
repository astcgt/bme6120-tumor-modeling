"""專案路徑：程式在 code/，輸出分成 output/data（CSV／JSON）與 output/plots（PNG／SVG）。

Project paths: code lives in code/; outputs split into output/data (CSV / JSON) and output/plots (PNG / SVG).
"""

from __future__ import annotations

from pathlib import Path

# 專案根目錄是 code/ 的上一層。
# The project root is the parent of code/.
ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "output" / "data"
PLOTS = ROOT / "output" / "plots"


# 把資料目錄對應到同名的圖目錄，例如 output/data/phase_plots → output/plots/phase_plots。
# Map a data directory to the matching plot directory, e.g. output/data/phase_plots → output/plots/phase_plots.
def plots_for(directory) -> Path:
    directory = Path(directory).resolve()
    try:
        return PLOTS / directory.relative_to(DATA)
    except ValueError:
        return directory
