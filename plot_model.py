"""數學圖：單細胞規則、從 M 矩陣到流向、附錄（完整轉移核與程式流程）。

矩陣與向量場由 analysis.py 計算；這裡只負責畫圖與存檔。
"""

from __future__ import annotations

import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
import numpy as np

from analysis import CELL_UNIT, lineage_extinction
from experiments import OUTPUTS

OUT = OUTPUTS / "mathematical_model"
INK, BLUE, RED, GRAY = "#182634", "#21658B", "#A54C38", "#66717A"


# 將座標軸設為隱藏刻度的單位畫布，加入面板編號與標題。
def panel(ax, letter, title, gap=0.055):
    ax.set(xlim=(0, 1), ylim=(0, 1))
    ax.axis("off")
    ax.text(0, 1, letter, fontsize=17, weight="bold", va="top", color=BLUE)
    ax.text(gap, 0.995, title, fontsize=13, weight="bold", va="top", color=INK)


# 以統一的預設字級、顏色與頂端對齊方式，在指定位置加入文字。
def label(ax, x, y, text, size=12, **kw):
    return ax.text(
        x, y, text, fontsize=size, color=kw.pop("color", INK), va="top", **kw
    )


# 在指定中心繪製圓角流程節點並加入置中文字。
def node(ax, xy, text, width=0.25, height=0.115, color=BLUE):
    x, y = xy
    patch = FancyBboxPatch(
        (x - width / 2, y - height / 2),
        width,
        height,
        boxstyle="round,pad=0.008,rounding_size=0.015",
        facecolor="#F3F6F8",
        edgecolor=color,
        linewidth=1.1,
    )
    ax.add_patch(patch)
    ax.text(x, y, text, ha="center", va="center", fontsize=11, color=INK)


# 在兩點間加入可彎曲的箭頭，並依選項標示說明文字。
def arrow(ax, a, b, text=None, at=None, rad=0, color=GRAY):
    ax.add_patch(
        FancyArrowPatch(
            a,
            b,
            arrowstyle="-|>",
            mutation_scale=12,
            linewidth=1.2,
            color=color,
            connectionstyle=f"arc3,rad={rad}",
        )
    )
    if text:
        x, y = at if at else ((a[0] + b[0]) / 2, (a[1] + b[1]) / 2)
        ax.text(
            x,
            y,
            text,
            fontsize=11,
            ha="center",
            va="center",
            color=INK,
            bbox=dict(facecolor="white", edgecolor="none", pad=2),
        )


# 以列至欄的慣例繪製數值矩陣，並依選項顯示列欄標籤。
def matrix_display(
    ax,
    values,
    x,
    y,
    width,
    height,
    rows=None,
    cols=None,
    fmt=".6g",
    size=12,
    colors=(BLUE, GRAY),
    lw=1.3,
):
    """Draw an actual numeric matrix, with conventional row-to-column orientation."""
    nr, nc = values.shape
    cw, ch = width / nc, height / nr
    for i in range(nr):
        for j in range(nc):
            ax.text(
                x + (j + 0.5) * cw,
                y - (i + 0.5) * ch,
                format(values[i, j], fmt),
                ha="center",
                va="center",
                fontsize=size,
                color=colors[0] if values[i, j] else colors[1],
            )
    for bx, inward in ((x, 0.015), (x + width, -0.015)):
        ax.plot(
            [bx + inward, bx, bx, bx + inward],
            [y, y, y - height, y - height],
            color=INK,
            lw=lw,
        )
    if rows:
        for i, name in enumerate(rows):
            ax.text(
                x - 0.035,
                y - (i + 0.5) * ch,
                name,
                ha="right",
                va="center",
                fontsize=size,
            )
    if cols:
        for j, name in enumerate(cols):
            ax.text(
                x + (j + 0.5) * cw,
                y + 0.04,
                name,
                ha="center",
                va="bottom",
                fontsize=size,
            )


# 將模型圖儲存為 300 dpi PNG 與 SVG，然後關閉圖形。
def save(fig, name):
    for suffix in ("png", "svg"):
        fig.savefig(
            OUT / f"{name}.{suffix}",
            dpi=300,
            facecolor="white",
            bbox_inches="tight",
            pad_inches=0.3,
        )
    plt.close(fig)


# 一個細胞一代會怎樣：子代分布、符號說明，以及四種代表性細胞的數字。
def single_cell_figure(config, offspring, phenotypes):
    """一個細胞一代會怎樣：子代分布、符號說明，以及四種代表性細胞的數字。"""
    fig = plt.figure(figsize=(15, 11.2))
    fig.text(
        0.04,
        0.965,
        "Single-cell rules: what one cell does in one generation",
        fontsize=22,
        weight="bold",
        color=INK,
    )

    ax = fig.add_axes([0.04, 0.52, 0.44, 0.36])
    panel(ax, "A", "Offspring distribution of a mutation-free cell")
    label(ax, 0.015, 0.86, r"$d=0.30,\quad q=0.04,\quad\mathrm{no\ therapy}$", 14)
    node(ax, (0.16, 0.54), "1 parent", 0.23, 0.15)
    for y, n, p in zip((0.73, 0.50, 0.27), (0, 2, 3), offspring[[0, 2, 3]]):
        node(ax, (0.80, y), f"{n} offspring", 0.29, 0.14, RED if n == 0 else BLUE)
        arrow(ax, (0.285, 0.54), (0.645, y), f"{p:.3f}", (0.48, (y + 0.54) / 2 + 0.035))
    label(ax, 0.015, 0.12, r"$f(z)=0.300+0.672z^2+0.028z^3$", 16)
    label(
        ax,
        0.015,
        -0.005,
        r"$m=f'(1)=1.428,\qquad s=f(s)\Rightarrow s=%.3f$"
        % lineage_extinction(config.death_rate, config.triple_probability),
        15,
    )

    # B：符號說明，讓 A 與 C 的每個變數都有文字定義。
    ax = fig.add_axes([0.52, 0.52, 0.44, 0.36])
    panel(ax, "B", "What the symbols mean")
    terms = [
        (r"$d$", "death probability per generation (cell leaves no offspring)"),
        (r"$q$", "probability that a surviving cell makes 3 daughters (else 2)"),
        (r"$k$", "extra kill probability from therapy (sensitive 0.90, resistant 0.04)"),
        (r"$P(n)$", "probability that one parent leaves n offspring"),
        (r"$f(z)$", "offspring generating function; coefficient of $z^n$ = P(n)"),
        (r"$m$", "mean offspring per cell = (1−d)(1−k)(2+q); m > 1 grows"),
        (r"$s$", "probability that one cell's lineage eventually dies out"),
    ]
    for i, (sym, text) in enumerate(terms):
        y = 0.84 - i * 0.12
        label(ax, 0.02, y, sym, 14)
        label(ax, 0.16, y - 0.005, text, 11.5, color=GRAY)

    ax = fig.add_axes([0.04, 0.09, 0.92, 0.34])
    panel(ax, "C", "Numbers that decide the outcome (default parameters)", gap=0.025)
    columns = [
        "Parent phenotype",
        "d\ndeath prob.",
        "q\nP(3 daughters\n| survives)",
        "P(0)\nno offspring",
        "P(2)\n2 daughters",
        "P(3)\n3 daughters",
        "m, no therapy\nmean offspring",
        "m, on therapy\nmean offspring",
    ]
    # 滅絕機率不放進表格；只保留決定成長方向的 d、q、子代機率與 m。
    table = ax.table(
        cellText=[row[: len(columns)] for row in phenotypes],
        colLabels=columns,
        colWidths=[0.22, 0.09, 0.12, 0.10, 0.10, 0.10, 0.135, 0.135],
        bbox=[0, 0.08, 1, 0.80],
        cellLoc="center",
    )
    table.auto_set_font_size(False)
    table.set_fontsize(11)
    for (r, c), cell in table.get_celld().items():
        cell.set_edgecolor("#CED5DA")
        cell.set_linewidth(0.5)
        cell.set_text_props(color=INK)
        if r == 0:
            cell.set_facecolor("#E9EFF3")
            cell.set_text_props(weight="bold", fontsize=10, linespacing=1.25)
        elif c in (6, 7):
            cell.set_text_props(
                weight="bold",
                color=RED if float(cell.get_text().get_text()) < 1 else BLUE,
            )
    save(fig, "single_cell_rules")


# 期刊風格配色：Okabe–Ito 朱紅標示代表試驗，其他試驗用淺灰；向量場用低飽和連續色階。
ACCENT, TRIAL_GRAY = "#D55E00", "#B8C0C6"
JOURNAL_RC = {
    "font.family": ["Arial", "DejaVu Sans"],
    "mathtext.fontset": "custom",
    "mathtext.rm": "Arial",
    "mathtext.it": "Arial:italic",
    "mathtext.bf": "Arial:bold",
    "font.size": 7,
    "axes.linewidth": 0.6,
    "axes.labelsize": 7,
    "axes.titlesize": 7.5,
    "xtick.labelsize": 6.5,
    "ytick.labelsize": 6.5,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "xtick.major.size": 2.5,
    "ytick.major.size": 2.5,
    "xtick.direction": "out",
    "ytick.direction": "out",
    "legend.fontsize": 6.5,
    "svg.fonttype": "none",
    "axes.unicode_minus": False,
}


# 實際模擬路徑（完整基因型模型）：未治療面板畫治療前那段，治療面板畫治療開始之後那段。
def draw_paths(ax, trials, events, highlight, therapy_phase):
    """實際模擬路徑（完整基因型模型）：未治療面板畫治療前那段，治療面板畫治療開始之後那段。"""
    for rep, states in trials.items():
        start = events[rep]["therapy_generation"]
        part = states[: start + 1] if not therapy_phase else states[start:]
        xy = np.array([[r["sensitive"], r["resistant"]] for r in part]) / 1e4
        if rep == highlight:
            continue
        ax.plot(*xy.T, color=TRIAL_GRAY, lw=0.5, alpha=0.8, zorder=2)
    start = events[highlight]["therapy_generation"]
    states = trials[highlight]
    part = states[: start + 1] if not therapy_phase else states[start:]
    xy = np.array([[r["sensitive"], r["resistant"]] for r in part]) / 1e4
    ax.plot(*xy.T, color=ACCENT, lw=1.1, zorder=5)
    ax.scatter(*xy.T, s=5, color=ACCENT, edgecolor="white", lw=0.3, zorder=6)
    ax.scatter(*xy[0 if therapy_phase else -1], s=16, marker="D", color=INK, zorder=7)


# 線性座標的一步期望向量場；兩面板共用倍率與色階，不做單支箭頭正規化。
def draw_field(ax, field, scale, norm, title, subtitle, show_ylabel=True):
    """線性座標的一步期望向量場；兩面板共用倍率與色階，不做單支箭頭正規化。"""
    quiver = ax.quiver(
        field["x"],
        field["y"],
        field["display_u"],
        field["display_v"],
        field["magnitude_cells"],
        angles="xy",
        scale_units="xy",
        scale=scale,
        pivot="tail",
        cmap="cividis",
        norm=norm,
        width=0.006,
        headwidth=3.2,
        headlength=3.8,
        headaxislength=3.4,
        minlength=0,
        minshaft=1,
    )
    ax.scatter(0, 0, s=6, color=INK, zorder=4)
    ax.set_aspect("equal", adjustable="box")
    ax.set(xlim=(-0.08, 3.2), ylim=(-0.08, 3.2))
    ax.set_xticks(np.arange(0, 3.1, 1.0))
    ax.set_yticks(np.arange(0, 3.1, 1.0))
    ax.set_xticks(np.arange(0, 3.1, 0.5), minor=True)
    ax.set_yticks(np.arange(0, 3.1, 0.5), minor=True)
    ax.tick_params(which="minor", length=1.5, width=0.5)
    ax.set_xlabel(r"Sensitive cells, $S$ ($\times10^4$)")
    if show_ylabel:
        ax.set_ylabel(r"Resistant cells, $R$ ($\times10^4$)")
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.spines["left"].set_position(("outward", 3))
    ax.spines["bottom"].set_position(("outward", 3))
    ax.set_title(title, loc="left", fontsize=7.5, weight="bold", color=INK, pad=12)
    ax.text(0, 1.02, subtitle, transform=ax.transAxes, fontsize=6.5, color=GRAY)
    return quiver


# M 矩陣（數字）與向量場（形狀）放同一頁；附條件期望公式。期刊風格：雙欄寬、7 pt 字。
def matrices_to_flow_figure(
    config, qmatrix, before, during, fields, trials, events, highlight
):
    """M 矩陣（數字）與向量場（形狀）放同一頁；附條件期望公式。期刊風格：雙欄寬、7 pt 字。"""
    with plt.rc_context(JOURNAL_RC):
        _matrices_to_flow(config, qmatrix, before, during, fields, trials, events, highlight)


def _matrices_to_flow(config, qmatrix, before, during, fields, trials, events, highlight):
    scale = 9.0  # 11×11 網格比原本 19×19 稀疏，箭頭可以畫長一點；兩面板共用。
    fig = plt.figure(figsize=(7.2, 4.0))
    fig.text(
        0.012,
        0.975,
        "From the mean matrix M to the flow of the tumor",
        fontsize=9,
        weight="bold",
        color=INK,
        va="top",
    )
    fig.text(
        0.012,
        0.93,
        "Two-type reduction: sensitive (S) and resistant (R). The numbers on the left"
        " are the arrows on the right.",
        fontsize=7,
        color=GRAY,
        va="top",
    )
    letter = dict(fontsize=9, weight="bold", color=INK, va="bottom")

    # a：矩陣與條件期望。
    ax = fig.add_axes([0.012, 0.22, 0.30, 0.62])
    ax.set(xlim=(0, 1), ylim=(0, 1))
    ax.axis("off")
    ax.text(0, 1.0, "a", **letter)
    ax.text(0.07, 1.0, "Mutation and mean-offspring matrices", fontsize=7.5,
            weight="bold", color=INK, va="bottom")
    ax.text(0.0, 0.93, "Reduced slice: resistance bit only; other growth effects off.",
            fontsize=6.3, color=GRAY, va="top")
    mat = dict(size=7, colors=(INK, "#9AA3AA"), lw=0.7)
    ax.text(0.0, 0.73, r"$Q=$", fontsize=8.5, va="center")
    matrix_display(ax, qmatrix, 0.22, 0.79, 0.42, 0.12, ["S", "R"], ["S", "R"], **mat)
    ax.text(0.72, 0.73, r"$\mu=0.001$", fontsize=7.5, va="center")
    ax.text(0.0, 0.56, r"$M_{gh}=m_g\,Q_{gh}$", fontsize=8, va="center")
    ax.text(0.0, 0.445, r"$M_{\mathrm{off}}=$", fontsize=8, va="center")
    matrix_display(ax, before, 0.25, 0.49, 0.66, 0.09, fmt=".7g", **mat)
    ax.text(0.0, 0.315, r"$M_{\mathrm{on}}=$", fontsize=8, va="center")
    matrix_display(ax, during, 0.25, 0.36, 0.66, 0.09, fmt=".7g", **mat)
    ax.plot([0, 0.95], [0.235, 0.235], color="#D5DBDF", lw=0.5)
    for y, text in (
        (0.18, r"$\mathrm{E}[S_{t+1}\mid s]=(1-\mu)\,B_S(s)$"),
        (0.10, r"$\mathrm{E}[R_{t+1}\mid s]=B_R(s)+\mu\,B_S(s)$"),
        (0.02, r"$B_S=\sum_{g\in\mathcal{S}}n_gm_g,\quad B_R=\sum_{g\in\mathcal{R}}n_gm_g$"),
    ):
        ax.text(0.0, y, text, fontsize=7.5, va="center", color=INK)
    fig.text(
        0.012,
        0.16,
        "Full model: m_g uses each genotype's own d, q and therapy kill.\n"
        "M is NOT a probability matrix; row sums are m.",
        fontsize=6.3,
        color=GRAY,
        va="top",
    )

    # b、c：向量場與實際軌跡，共用倍率與色階。
    axes = [
        fig.add_axes([0.375, 0.30, 0.25, 0.45]),
        fig.add_axes([0.655, 0.30, 0.25, 0.45]),
    ]
    color_ax = fig.add_axes([0.925, 0.34, 0.007, 0.36])
    norm = Normalize(0, max(np.max(f["magnitude_cells"]) for f in fields))
    quivers = []
    for i, (axis, field) in enumerate(zip(axes, fields)):
        m_s, m_r = np.sum(field["matrix"][0]), np.sum(field["matrix"][1])
        title = "No therapy" if i == 0 else "Continuous therapy"
        subtitle = rf"$m_S={m_s:.4g},\ m_R={m_r:.4g}$"
        quivers.append(
            draw_field(axis, field, scale, norm, title, subtitle, show_ylabel=i == 0)
        )
        draw_paths(axis, trials, events, highlight, therapy_phase=i == 1)
        axis.text(-0.17, 1.1, "bc"[i], transform=axis.transAxes, **letter)
        a, b, c = (
            field["matrix"][0, 0] - 1,
            field["matrix"][0, 1],
            field["matrix"][1, 1] - 1,
        )
        axis.text(
            0.0,
            -0.30,
            rf"$\Delta S={a:.4g}\,S$" + "\n" + rf"$\Delta R={b:.4g}\,S+{c:.4g}\,R$",
            transform=axis.transAxes,
            fontsize=7,
            color=INK,
            va="top",
            linespacing=1.5,
        )
    bar = fig.colorbar(quivers[0], cax=color_ax)
    bar.set_label(r"$\|\Delta\|$ ($10^4$ cells / generation)", fontsize=6.5)
    bar.ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{v / 1e4:g}"))
    bar.ax.tick_params(labelsize=6, width=0.5, length=2)
    bar.outline.set_linewidth(0.5)
    axes[0].quiverkey(
        quivers[0],
        X=0.375,
        Y=0.075,
        U=1,
        label=f"{CELL_UNIT:,} cells / generation",
        labelpos="E",
        coordinates="figure",
        fontproperties={"size": 6.5},
        color=INK,
    )
    fig.legend(
        handles=[
            Line2D([], [], color=TRIAL_GRAY, lw=0.8, label=f"{len(trials)} simulated trials"),
            Line2D(
                [], [], color=ACCENT, lw=1.1, marker="o", ms=2.5,
                label=f"Trial {highlight} (one dot per generation)",
            ),
            Line2D([], [], color=INK, marker="D", ms=3, lw=0, label="Therapy starts"),
        ],
        loc="center left",
        bbox_to_anchor=(0.53, 0.075),
        ncol=3,
        columnspacing=1.4,
        frameon=False,
        handlelength=1.8,
        labelspacing=0.3,
    )
    fig.text(
        0.012,
        -0.005,
        "Arrows = one-generation expected change E[(S',R') | s] − (S,R) on the baseline"
        " slice; shared length scale. Lines = actual trials of the full model (with"
        " oncogenes):\nb shows each trial before therapy, c from therapy start on."
        " Trials run along the arrows: rightward before therapy, then left to S = 0 and"
        " up the R axis (beyond the plotted range).\nOn linear axes the paths hug the"
        " edges because R stays below a few hundred cells before therapy; phase"
        " portrait 01 shows the same paths in log coordinates.",
        fontsize=5.8,
        color=GRAY,
        va="top",
        linespacing=1.4,
    )
    assert quivers[0].scale == quivers[1].scale
    save(fig, "matrices_to_flow")


# 附錄：完整基因型計數過程的 PGF 轉移核，以及程式的執行狀態機。
def appendix_figure(config):
    """附錄：完整基因型計數過程的 PGF 轉移核，以及程式的執行狀態機。"""
    fig = plt.figure(figsize=(15, 14))
    fig.text(
        0.04,
        0.972,
        "Appendix A1 | Full Markov kernel and execution controller",
        fontsize=22,
        weight="bold",
        color=INK,
    )
    fig.text(
        0.04,
        0.948,
        "Exact population transition kernel; the controller decides when to record,"
        " stop and start therapy",
        fontsize=12,
        color=GRAY,
    )
    ax = fig.add_axes([0.04, 0.64, 0.92, 0.28])
    panel(ax, "A", "Full genotype-count process", gap=0.025)
    label(
        ax,
        0.015,
        0.84,
        r"$S_t=(\mathbf{N}_t,A_t),\quad \mathbf{N}_t\in\mathbb{N}_0^{2^{22}},\quad"
        r" A_t\in\{0,1\}$",
        17,
    )
    label(
        ax,
        0.015,
        0.65,
        r"$\alpha=A_t\vee\mathbf{1}\{\mathrm{therapy\ enabled},\ \sum_g"
        r" N_{g,t}\geq20{,}000\}$",
        15,
    )
    label(
        ax,
        0.015,
        0.47,
        r"$\delta_g=1-(1-d_g)(1-\alpha k_g),\qquad"
        r" u_g(\mathbf{z})=\sum_hQ_{gh}z_h$",
        15,
    )
    label(
        ax,
        0.015,
        0.29,
        "Fixed Config; kernel applies before stopping.",
        11,
        color=GRAY,
    )
    label(
        ax,
        0.57,
        0.82,
        r"$F_g(\mathbf{z}\mid"
        r" S)=\delta_g+(1-\delta_g)\left[(1-q_g)u_g^2+q_gu_g^3\right]$",
        15,
    )
    label(
        ax,
        0.57,
        0.61,
        r"$P_{S,S'}=\mathbf{1}\{A'=\alpha\}\,[\mathbf{z}^{\mathbf{N}'}]\prod_gF_g(\mathbf{z}\mid"
        r" S)^{N_g}$",
        16,
    )
    label(
        ax,
        0.57,
        0.39,
        r"$Q_{gh}=\mu^{|h|-|g|}(1-\mu)^{22-|h|}\quad\mathrm{if}\ g\subseteq h;\quad"
        r" 0\ \mathrm{otherwise}$",
        14,
    )
    label(
        ax,
        0.57,
        0.19,
        "Coefficient extraction gives the exact next-population probability.\nQ is"
        " inherited, irreversible, independent-bit mutation; mu = 0.001.",
        11,
    )

    ax = fig.add_axes([0.04, 0.07, 0.52, 0.52])
    panel(ax, "B", "Execution state machine")
    node(ax, (0.50, 0.87), "Initialize\nN = 100, A = 0", 0.45, 0.09)
    node(
        ax, (0.50, 0.73), "Record current state\nCheck stopping conditions", 0.62, 0.10
    )
    arrow(ax, (0.50, 0.822), (0.50, 0.785))
    node(ax, (0.50, 0.54), "Continue\nNo stop condition holds", 0.53, 0.09)
    arrow(ax, (0.50, 0.675), (0.50, 0.591), "otherwise", (0.50, 0.637))
    node(ax, (0.50, 0.32), "Update treatment\nalpha = A OR (N >= 20,000)", 0.68, 0.11)
    arrow(ax, (0.50, 0.491), (0.50, 0.380))
    node(
        ax,
        (0.50, 0.14),
        "Draw offspring + mutations\nN <- N', A <- alpha; t <- t+1",
        0.68,
        0.105,
    )
    arrow(ax, (0.50, 0.26), (0.50, 0.197))
    arrow(ax, (0.151, 0.14), (0.05, 0.14))
    ax.plot([0.05, 0.05], [0.14, 0.73], color=GRAY, lw=1.2)
    arrow(ax, (0.05, 0.73), (0.18, 0.73))
    node(ax, (0.88, 0.62), "STOP", 0.18, 0.055, RED)
    arrow(ax, (0.814, 0.73), (0.88, 0.655), rad=-0.12, color=RED)
    label(
        ax,
        0.09,
        0.445,
        "A = 0: untreated  |  A = 1: persistent therapy",
        11,
        bbox=dict(facecolor="white", edgecolor="none", pad=2),
    )
    ax = fig.add_axes([0.62, 0.07, 0.34, 0.52])
    ax.axis("off")
    notes = [
        "Stop priority:\nextinction > target_size > max_cells > final generation.",
        "Treatment update is shown for therapy=True;\notherwise alpha = A.",
        "Therapy never switches off once started (A is absorbing).",
        (
            "N in the controller is total burden; the transition\ndraws use the full"
            " genotype vector."
        ),
        (
            "Relapse is a trajectory classification made after the\nrun, not a"
            " controller state."
        ),
    ]
    for i, text in enumerate(notes):
        label(ax, 0, 0.86 - i * 0.16, text, 12)
    save(fig, "markov_kernel_and_controller")
