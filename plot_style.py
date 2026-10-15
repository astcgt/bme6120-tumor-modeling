"""繪圖共用設定：matplotlib 後端、配色、字型樣式與存檔。

每個 set_*_style 先把 matplotlib 還原成初始設定，再套用該組圖的樣式；
因此 main.py 在同一個行程依序畫所有圖時，前一組圖的字型不會影響下一組。
"""

from __future__ import annotations

import os
from pathlib import Path

# 圖形快取留在專案中，不寫入 conda 或使用者的 matplotlib 設定目錄。
os.environ.setdefault(
    "MPLCONFIGDIR", str(Path(__file__).resolve().parent / ".mplconfig")
)
import matplotlib

# 使用不需要視窗的繪圖後端，適合直接輸出簡報圖片。
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager

# 相圖與基因型圖的雙語版本。
LANGS = ("zh", "en")
# 相圖、基因型圖共用的文字與座標軸顏色。
INK, GRAY = "#173042", "#687783"
# 集中定義展示圖的配色，讓各情境與面板使用一致的視覺樣式。
COLORS = ["#167D9A", "#ED7953", "#6C5B9C", "#57A773", "#D6A43B", "#BC5B86"]


def reset_style():
    """還原成 matplotlib 載入時的設定（保留 Agg 後端）。"""
    original = {k: v for k, v in matplotlib.rcParamsOrig.items() if k != "backend"}
    plt.rcParams.update(original)


def set_plot_style():
    """統一簡報用視覺樣式；英文圖標籤避免依賴中文系統字型。"""
    reset_style()
    # 統一字型、線寬、背景與淡色網格，並移除上方及右側邊框。
    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 11,
        "axes.titlesize": 13,
        "axes.titleweight": "bold",
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.labelcolor": "#334155",
        "text.color": "#172B42",
        "axes.grid": True,
        "grid.alpha": 0.16,
        "figure.facecolor": "#FAFBFD",
        "axes.facecolor": "white",
        "savefig.facecolor": "#FAFBFD",
        "lines.linewidth": 2.3,
    })


def save_figure(fig, output: Path, name: str):
    # 建立輸出目錄；同時支援呼叫端傳入字串或 Path 路徑。
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    # 同時輸出高解析度點陣 PNG 與可縮放向量 SVG，緊密裁切多餘留白。
    for extension in ("png", "svg"):
        fig.savefig(output / f"{name}.{extension}", dpi=240, bbox_inches="tight")
    # 儲存後釋放圖形資源，避免重複繪圖時累積記憶體。
    plt.close(fig)


def set_model_style():
    """數學圖：DejaVu Sans 與 STIX 數學字型；SVG 保留可編輯文字。"""
    reset_style()
    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "mathtext.fontset": "stix",
        "svg.fonttype": "none",
        "axes.unicode_minus": False,
    })


def set_phase_style():
    """中文標籤及 Computer Modern 數學字型；SVG 保留字形路徑。"""
    reset_style()
    available = {f.name for f in font_manager.fontManager.ttflist}
    candidates = (
        "Arial Unicode MS",
        "Microsoft JhengHei",
        "Noto Sans CJK TC",
        "PingFang TC",
    )
    family = next((f for f in candidates if f in available), None)
    if family is None:
        raise RuntimeError(
            "需要支援繁體中文的字型：Arial Unicode MS 或 Noto Sans CJK TC。"
        )
    plt.rcParams.update({
        "font.family": family,
        "font.size": 13,
        "mathtext.fontset": "cm",
        "axes.unicode_minus": False,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.labelcolor": INK,
        "text.color": INK,
        "axes.edgecolor": GRAY,
        "xtick.color": GRAY,
        "ytick.color": GRAY,
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "svg.fonttype": "path",
    })


def set_style(lang: str):
    """中文沿用既有相圖字型；英文用 DejaVu Sans，兩者數學式皆為 Computer Modern。"""
    reset_style()
    if lang == "zh":
        set_phase_style()
    else:
        plt.rcParams.update({
            "font.family": "DejaVu Sans",
            "font.size": 13,
            "mathtext.fontset": "cm",
            "axes.unicode_minus": True,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.labelcolor": INK,
            "text.color": INK,
            "axes.edgecolor": GRAY,
            "xtick.color": GRAY,
            "ytick.color": GRAY,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "svg.fonttype": "path",
        })


def save_lang(fig, directory: Path, name: str, lang: str):
    """雙語圖：檔名加上 _zh／_en，PNG 300 dpi 與 SVG 各一份。"""
    directory.mkdir(parents=True, exist_ok=True)
    for extension in ("png", "svg"):
        fig.savefig(
            directory / f"{name}_{lang}.{extension}", dpi=300, facecolor="white"
        )
    plt.close(fig)
