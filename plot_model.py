"""數學圖：06 單細胞規則、07 從 M 矩陣到流向、A1 附錄（完整轉移核與程式流程）。

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


def panel(ax, letter, title, gap=0.055):
    ax.set(xlim=(0, 1), ylim=(0, 1))
    ax.axis("off")
    ax.text(0, 1, letter, fontsize=17, weight="bold", va="top", color=BLUE)
    ax.text(gap, 0.995, title, fontsize=13, weight="bold", va="top", color=INK)


def label(ax, x, y, text, size=12, **kw):
    return ax.text(
        x, y, text, fontsize=size, color=kw.pop("color", INK), va="top", **kw
    )


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


def matrix_display(ax, values, x, y, width, height, rows=None, cols=None, fmt=".6g"):
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
                fontsize=12,
                color=BLUE if values[i, j] else GRAY,
            )
    for bx, inward in ((x, 0.015), (x + width, -0.015)):
        ax.plot(
            [bx + inward, bx, bx, bx + inward],
            [y, y, y - height, y - height],
            color=INK,
            lw=1.3,
        )
    if rows:
        for i, name in enumerate(rows):
            ax.text(
                x - 0.035,
                y - (i + 0.5) * ch,
                name,
                ha="right",
                va="center",
                fontsize=12,
            )
    if cols:
        for j, name in enumerate(cols):
            ax.text(
                x + (j + 0.5) * cw,
                y + 0.04,
                name,
                ha="center",
                va="bottom",
                fontsize=12,
            )


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


def single_cell_figure(config, matrix, offspring, phenotypes):
    """一個細胞一代會怎樣：子代分布、有限計數鏈，以及四種代表性細胞的數字。"""
    fig = plt.figure(figsize=(15, 11.2))
    fig.text(
        0.04,
        0.965,
        "Single-cell rules: what one cell does in one generation",
        fontsize=22,
        weight="bold",
        color=INK,
    )
    fig.text(
        0.04,
        0.93,
        "Every simulated outcome is built from these per-cell offspring laws",
        fontsize=12,
        color=GRAY,
    )

    ax = fig.add_axes([0.04, 0.52, 0.40, 0.36])
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

    ax = fig.add_axes([0.49, 0.52, 0.47, 0.36])
    panel(ax, "B", "Same rule as a finite Markov matrix")
    label(
        ax,
        0.015,
        0.86,
        r"$P_{ij}=[z^j]f(z)^i;\qquad \boldsymbol{\pi}_{t+1}=\boldsymbol{\pi}_tP$",
        15,
    )
    labels = ["0", "1", "2", "3", r"$\dagger$"]
    matrix_display(ax, matrix, 0.12, 0.65, 0.83, 0.47, labels, labels)
    label(
        ax,
        0.015,
        0.08,
        r"$\dagger$: first entry into $N\geq4$; bookkeeping stop, not a carrying"
        r" capacity.",
        11,
    )
    label(
        ax,
        0.015,
        -0.005,
        "Rows: current cells. Columns: next generation. Each row sums to 1; 0 is"
        " absorbing.",
        11,
        color=GRAY,
    )

    ax = fig.add_axes([0.04, 0.09, 0.92, 0.34])
    panel(ax, "C", "Numbers that decide the outcome (default parameters)", gap=0.025)
    columns = [
        "Parent phenotype",
        "d",
        "q",
        "P(0)",
        "P(2)",
        "P(3)",
        "m: OFF",
        "m: ON",
        "Extinct: OFF",
        "Extinct: ON",
    ]
    table = ax.table(
        cellText=phenotypes,
        colLabels=columns,
        colWidths=[0.20, 0.06, 0.06, 0.075, 0.075, 0.075, 0.09, 0.09, 0.135, 0.135],
        bbox=[0, 0.22, 1, 0.66],
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
            cell.set_text_props(weight="bold")
        elif c in (6, 7):
            cell.set_text_props(
                weight="bold",
                color=RED if float(cell.get_text().get_text()) < 1 else BLUE,
            )
    label(
        ax,
        0.015,
        0.15,
        "P(0/2/3) are offspring probabilities without therapy. m = expected offspring ="
        " (1−d)(1−k)(2+q); m > 1 grows, m < 1 shrinks.",
        11.5,
    )
    label(
        ax,
        0.015,
        0.06,
        "Therapy (ON) adds extra kill k = 0.90 for sensitive and 0.04 for resistant"
        " cells. Resistance is neutral without therapy (same m = 1.428).",
        11.5,
    )
    fig.text(
        0.04,
        0.035,
        "Panels A–B drop mutations, genotype advantages and therapy; the threshold 4 in"
        " B is illustrative. "
        "Extinct = probability that one cell's lineage eventually dies out (smallest"
        " root of s = f(s)).\n"
        "Starting from 100 baseline cells, total extinction without therapy has"
        " probability s^100 ≈ "
        f"{lineage_extinction(config.death_rate, config.triple_probability)**100:.0e}.",
        fontsize=10,
        color=GRAY,
    )
    save(fig, "06_single_cell_rules")


def draw_paths(ax, trials, events, highlight, therapy_phase):
    """實際模擬路徑（完整基因型模型）：未治療面板畫治療前那段，治療面板畫治療開始之後那段。"""
    for rep, states in trials.items():
        start = events[rep]["therapy_generation"]
        part = states[: start + 1] if not therapy_phase else states[start:]
        xy = np.array([[r["sensitive"], r["resistant"]] for r in part]) / 1e4
        if rep == highlight:
            continue
        ax.plot(*xy.T, color="#9AA7B0", lw=0.9, alpha=0.55, zorder=2)
    start = events[highlight]["therapy_generation"]
    states = trials[highlight]
    part = states[: start + 1] if not therapy_phase else states[start:]
    xy = np.array([[r["sensitive"], r["resistant"]] for r in part]) / 1e4
    ax.plot(*xy.T, color="#D9822B", lw=2.4, zorder=5)
    ax.scatter(*xy.T, s=14, color="#D9822B", edgecolor="white", lw=0.5, zorder=6)
    ax.scatter(*xy[0 if therapy_phase else -1], s=60, marker="D", color=INK, zorder=7)


def draw_field(ax, field, scale, norm, title, show_ylabel=True):
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
        cmap="viridis",
        norm=norm,
        width=0.0045,
        headwidth=3.4,
        headlength=4.4,
        headaxislength=3.8,
        minlength=0,
        minshaft=1,
    )
    ax.scatter(0, 0, s=18, color=INK, zorder=4)
    ax.set_aspect("equal", adjustable="box")
    ax.set(xlim=(-0.08, 3.2), ylim=(-0.08, 3.2))
    ax.set_xticks(np.arange(0, 3.1, 0.5))
    ax.set_yticks(np.arange(0, 3.1, 0.5))
    ax.set_xlabel(r"Sensitive cells $S$ ($\times10^4$)", fontsize=12)
    if show_ylabel:
        ax.set_ylabel(r"Resistant cells $R$ ($\times10^4$)", fontsize=12)
    ax.grid(color="#E1E7EB", lw=0.65, alpha=0.75)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.set_title(title, loc="left", fontsize=13, weight="bold", color=INK, pad=10)
    return quiver


def matrices_to_flow_figure(
    config, qmatrix, before, during, fields, trials, events, highlight
):
    """M 矩陣（數字）與向量場（形狀）放同一頁；附條件期望公式。"""
    scale = 9.0  # 11×11 網格比原本 19×19 稀疏，箭頭可以畫長一點；兩面板共用。
    fig = plt.figure(figsize=(17, 10.2))
    fig.text(
        0.035,
        0.965,
        "From the mean matrix M to the flow of the tumor",
        fontsize=22,
        weight="bold",
        color=INK,
    )
    fig.text(
        0.035,
        0.93,
        "Two-type reduction: sensitive (S) and resistant (R). The numbers on the left"
        " are the arrows on the right.",
        fontsize=12,
        color=GRAY,
    )
    ax = fig.add_axes([0.035, 0.10, 0.34, 0.78])
    panel(ax, "A", "Mutation and mean-offspring matrices")
    label(
        ax,
        0.015,
        0.90,
        "Reduced slice: resistance bit only; other growth effects and IL11 off.",
        10,
        color=GRAY,
    )
    label(ax, 0.015, 0.80, r"$Q=$", 17)
    matrix_display(ax, qmatrix, 0.20, 0.81, 0.52, 0.12, ["S", "R"], ["S", "R"])
    label(ax, 0.78, 0.79, r"$\mu=0.001$", 14)
    label(ax, 0.015, 0.62, r"$M_{gh}=m_gQ_{gh}$", 15)
    label(ax, 0.015, 0.52, r"$M_{\mathrm{off}}=$", 14)
    matrix_display(ax, before, 0.25, 0.52, 0.70, 0.085, fmt=".7g")
    label(ax, 0.015, 0.40, r"$M_{\mathrm{on}}=$", 14)
    matrix_display(ax, during, 0.25, 0.40, 0.70, 0.085, fmt=".7g")
    label(
        ax, 0.015, 0.26, r"$\mathrm{E}[S_{t+1}\mid s]=(1-\mu)\,B_S(s)$", 15, color=BLUE
    )
    label(
        ax,
        0.015,
        0.19,
        r"$\mathrm{E}[R_{t+1}\mid s]=B_R(s)+\mu\,B_S(s)$",
        15,
        color=BLUE,
    )
    label(
        ax,
        0.015,
        0.12,
        r"$B_S=\sum_{g\in\mathcal{S}}n_gm_g,\quad B_R=\sum_{g\in\mathcal{R}}n_gm_g$",
        14,
        color=BLUE,
    )
    label(
        ax,
        0.015,
        0.045,
        "Full model: m_g uses each genotype's own d, q and therapy kill.\n"
        "M is NOT a probability matrix; row sums are m.",
        10.5,
        color=GRAY,
    )

    axes = [
        fig.add_axes([0.42, 0.30, 0.24, 0.52]),
        fig.add_axes([0.69, 0.30, 0.24, 0.52]),
    ]
    color_ax = fig.add_axes([0.945, 0.35, 0.009, 0.42])
    norm = Normalize(0, max(np.max(f["magnitude_cells"]) for f in fields))
    quivers = []
    for i, (axis, field) in enumerate(zip(axes, fields)):
        m_s, m_r = np.sum(field["matrix"][0]), np.sum(field["matrix"][1])
        title = (
            "B  No therapy" if i == 0 else "C  Continuous therapy"
        ) + f"\nm_S = {m_s:.4g}, m_R = {m_r:.4g}"
        quivers.append(draw_field(axis, field, scale, norm, title, show_ylabel=i == 0))
        draw_paths(axis, trials, events, highlight, therapy_phase=i == 1)
        a, b, c = (
            field["matrix"][0, 0] - 1,
            field["matrix"][0, 1],
            field["matrix"][1, 1] - 1,
        )
        fig.text(
            0.42 + 0.27 * i,
            0.245,
            rf"$\Delta S={a:.4g}\,S$" + "\n" + rf"$\Delta R={b:.4g}\,S+{c:.4g}\,R$",
            fontsize=13.5,
            color=BLUE,
            va="top",
        )
    bar = fig.colorbar(quivers[0], cax=color_ax)
    bar.set_label(r"$\|\Delta\|$ (cells / generation)", fontsize=11)
    bar.ax.ticklabel_format(style="sci", scilimits=(4, 4))
    bar.ax.tick_params(labelsize=9)
    axes[0].quiverkey(
        quivers[0],
        X=0.43,
        Y=0.135,
        U=1,
        label=f"{CELL_UNIT:,} cells / generation",
        labelpos="E",
        coordinates="figure",
        fontproperties={"size": 11},
        color=BLUE,
    )
    axes[1].legend(
        handles=[
            Line2D(
                [], [], color="#9AA7B0", lw=1.2, label=f"{len(trials)} simulated trials"
            ),
            Line2D(
                [],
                [],
                color="#D9822B",
                lw=2.4,
                marker="o",
                ms=4,
                label=f"Trial {highlight} (one dot per generation)",
            ),
            Line2D([], [], color=INK, marker="D", lw=0, label="Therapy starts"),
        ],
        loc="upper left",
        bbox_to_anchor=(0.69, 0.19),
        bbox_transform=fig.transFigure,
        fontsize=10,
        frameon=False,
    )
    fig.text(
        0.42,
        0.07,
        "Arrows = one-generation expected change E[(S',R') | s] − (S,R) on the baseline"
        " slice; shared length scale. Lines = actual trials of the full model (with"
        " oncogenes):\nB shows each trial before therapy, C from therapy start on."
        " Trials run along the arrows: rightward before therapy, then left to S = 0 and"
        " up the R axis (beyond the plotted range).\nOn linear axes the paths hug the"
        " edges because R stays below a few hundred cells before therapy; phase"
        " portrait 01 shows the same paths in log coordinates.",
        fontsize=10.5,
        color=GRAY,
    )
    assert quivers[0].scale == quivers[1].scale
    save(fig, "07_matrices_to_flow")


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
        r"$b(\mathbf{N})=0.06\,f/(0.02+f),\qquad f=\sum_{g:\,g_6=1}N_g/\sum_gN_g$",
        15,
    )
    label(
        ax,
        0.015,
        0.29,
        r"$\delta_g=1-(1-d_g(\mathbf{N}))(1-\alpha k_g),\qquad"
        r" u_g(\mathbf{z})=\sum_hQ_{gh}z_h$",
        15,
    )
    label(
        ax,
        0.015,
        0.10,
        "IL11 OFF: b = 0. Empty population: b = 0. Fixed Config; kernel applies before"
        " stopping.",
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
    save(fig, "A1_markov_kernel_and_controller")
