"""相圖：實際軌跡與期望下一步（phase_portrait）、理論 vs 經驗相場（plot_pair）。

期望步伐與箭頭資料由 analysis.py、experiments.py 計算；這裡只負責畫圖與存檔。
"""

from __future__ import annotations

import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch
import numpy as np

from analysis import MIN_TRANSITIONS, event_generations, logit_rho

INK, BLUE, ORANGE, RED, GRAY = "#173042", "#236B8E", "#BC7844", "#A74F3D", "#687783"


# ---------------------------------------------------------------- 實際軌跡相圖 ------
# 將逐代計數轉為敏感與抗藥對數座標，或總負荷對數與抗藥百分比座標。
def phase_coordinates(states, composition=False):
    total = np.array([r["total"] for r in states], dtype=float)
    if composition:
        # 滅絕時比例未定義，用 NaN 讓圖斷開，而非畫成 0% 抗藥。
        fraction = np.divide(
            [r["resistant"] for r in states],
            total,
            out=np.full(len(states), np.nan),
            where=total > 0,
        )
        return np.column_stack((np.log10(total + 1), 100 * fraction))
    return np.log10(np.array([[r["sensitive"], r["resistant"]] for r in states]) + 1)


# 依相座標種類設定標題、網格、軸標籤及抗藥比例的顯示範圍。
def decorate_axis(ax, title, composition=False):
    ax.set_title(title, loc="left", fontsize=17, weight="bold", pad=16)
    ax.grid(color="#DDE4E8", alpha=0.6, linewidth=0.7)
    ax.set_axisbelow(True)
    if composition:
        ax.set_xlabel(r"總負荷 $\log_{10}(T+1)$", labelpad=12)
        ax.set_ylabel(r"抗藥比例 $100R/T$（%）", labelpad=12)
        ax.set_ylim(-4, 107)
    else:
        ax.set_xlabel(r"敏感細胞 $\log_{10}(S+1)$", labelpad=12)
        ax.set_ylabel(r"抗藥細胞 $\log_{10}(R+1)$", labelpad=12)
        ax.set_xlim(-0.16, 4.85)
        ax.set_ylim(-0.18, 6.55)


# 僅在起終點皆為有限數值時，於座標軸加入指定樣式的箭頭。
def arrow(ax, a, b, color, size=12, linewidth=2, zorder=6):
    if np.isfinite([*a, *b]).all():
        ax.add_patch(
            FancyArrowPatch(
                a,
                b,
                arrowstyle="-|>",
                mutation_scale=size,
                color=color,
                linewidth=linewidth,
                shrinkA=0,
                shrinkB=0,
                zorder=zorder,
            )
        )


# 繪製實際試驗軌跡，突出代表案例並疊加由完整基因型算出的下一代期望箭頭。
def plot_phase_portrait(trials, events, replicate, config, result, steps):
    """單一相平面：30 條實際軌跡、突出案例，以及由完整基因型狀態算出的期望下一代箭頭。

    取代舊的 01（軌跡兩格）與 02（條件期望）：同一試驗、同一座標不再畫兩次。
    負荷 × 抗藥比例的座標改看 phase_field_burden_logit（logit 軸）；推導公式在 matrices_to_flow。
    """
    fig = plt.figure(figsize=(15.5, 10.2))
    ax = fig.add_axes((0.07, 0.2, 0.5, 0.58))
    side = fig.add_axes((0.62, 0.2, 0.35, 0.58))
    side.axis("off")
    fig.text(
        0.055,
        0.965,
        "相平面：負荷下降不等於抗藥消失",
        fontsize=25,
        weight="bold",
        va="top",
    )
    event, states = events[replicate], trials[replicate]
    fig.text(
        0.055,
        0.912,
        f"{len(trials)} 次試驗的實際軌跡；突出案例 replicate"
        f" {replicate}（seed={event['seed']}）。"
        "黑色箭頭＝由當代完整基因型組成算出的期望下一代位置。",
        fontsize=13.5,
        color=GRAY,
    )
    start, relapse = event["therapy_generation"], event["relapse_generation"]
    # fmt: off
    colors = [
        BLUE if start is None or g < start
        else RED if relapse is not None and g >= relapse
        else ORANGE
        for g in range(len(states) - 1)
    ]
    # fmt: on
    for path in trials.values():
        ax.plot(
            *phase_coordinates(path).T, color="#BCC7CE", alpha=0.48, lw=0.85, zorder=1
        )
    xy = phase_coordinates(states)
    for g, color in enumerate(colors):
        ax.plot(*xy[g : g + 2].T, color=color, lw=2.6, zorder=3)
    ax.scatter(*xy.T, s=9, color=INK, alpha=0.55, zorder=4, linewidths=0)
    # 期望箭頭：先求原始計數期望再取 log；箭頭頭端與下一個實際點的距離就是隨機偏差。
    for g in sorted(
        {0, 5, 10, 13, 14, 15, 16, 17, 18, 19, 20, 22, 24, 28, 32, len(steps) - 1}
    ):
        if 0 <= g < len(steps):
            row = steps[g]
            arrow(
                ax,
                xy[g],
                np.log10(
                    [row["expected_sensitive"] + 1, row["expected_resistant"] + 1]
                ),
                INK,
                size=11,
                linewidth=1.5,
                zorder=6,
            )
    sensitive = np.r_[0, np.geomspace(1, config.therapy_start_size, 250)]
    resistant = np.maximum(0, config.therapy_start_size - sensitive)
    ax.plot(
        np.log10(sensitive + 1),
        np.log10(resistant + 1),
        color=GRAY,
        ls="--",
        lw=1,
        alpha=0.7,
    )
    ax.text(
        0.65,
        4.37,
        rf"$S+R={config.therapy_start_size:,}$（治療門檻）",
        color=GRAY,
        fontsize=11.5,
    )
    positions = {
        "初始": (2.75, 0.56),
        "治療啟動": (3.64, 3.5),
        "負荷谷底": (1.05, 2.08),
        "首次復發": (1.27, 5.15),
    }
    points = event_generations(states, event)
    for label, g in points.items():
        marker = {"初始": "o", "治療啟動": "D", "負荷谷底": "v", "首次復發": "s"}[label]
        ax.scatter(
            *xy[g],
            s=62,
            marker=marker,
            c=INK,
            edgecolors="white",
            linewidth=0.9,
            zorder=7,
        )
        ax.annotate(
            f"{label}\n$t={g}$",
            xy=xy[g],
            xytext=positions[label],
            fontsize=11.5,
            ha="center",
            arrowprops={"arrowstyle": "-", "color": GRAY, "lw": 0.9},
            bbox={"facecolor": "white", "edgecolor": "none", "pad": 2},
            zorder=8,
        )
    ax.scatter(*xy[-1], s=50, facecolor="white", edgecolor=RED, linewidth=1.5, zorder=7)
    decorate_axis(ax, "敏感 × 抗藥細胞（log 座標）")
    handles = [
        Line2D([], [], color="#BCC7CE", lw=1, label=f"全部 {len(trials)} 次試驗"),
        Line2D([], [], color=BLUE, lw=3, label="治療前"),
        Line2D([], [], color=ORANGE, lw=3, label="治療至首次復發"),
        Line2D([], [], color=RED, lw=3, label="復發後"),
        Line2D([], [], color=INK, lw=1.5, marker=">", markersize=7, label="期望下一代"),
    ]
    fig.legend(
        handles=handles,
        loc="upper left",
        bbox_to_anchor=(0.055, 0.875),
        frameon=False,
        ncol=5,
        fontsize=12,
    )
    # 右側：案例數字。期望與實際並列，說明大族群幾乎沒有偏差、小族群偏差變大。
    y = 0.97
    side.text(0, y, "讀圖", fontsize=16, weight="bold", va="top")
    y -= 0.08
    for text in (
        "治療前：S、R 一起長，路徑往右上。",
        "治療後：S 每代少約一個數量級，\n路徑幾乎水平往左。",
        "S 歸零後：沿縱軸直直往上，\n純抗藥族群在生長。",
    ):
        side.text(0, y, text, fontsize=12.5, va="top")
        y -= 0.035 + 0.05 * (1 + text.count("\n"))
    y -= 0.02
    if start is not None and "負荷谷底" in points:
        nadir = states[points["負荷谷底"]]
        side.text(
            0,
            y,
            "抗藥細胞在治療中從未減少",
            fontsize=14,
            weight="bold",
            va="top",
            color=RED,
        )
        y -= 0.065
        side.text(
            0,
            y,
            f"治療起點 R = {states[start]['resistant']:,}\n谷底 R ="
            f" {nadir['resistant']:,}"
            f"（占 {100*nadir['resistant']/nadir['total']:.1f}%）",
            fontsize=12.5,
            va="top",
        )
        y -= 0.14
    nadir_g = points.get("負荷谷底")
    for g, title in ((start, "治療啟動"), (nadir_g, "負荷谷底")):
        if g is not None and g < len(steps):
            row = steps[g]
            side.text(
                0,
                y,
                f"{title} t={g}：(S, R) = ({row['sensitive']:,}, {row['resistant']:,})",
                fontsize=12.5,
                weight="bold",
                va="top",
            )
            y -= 0.055
            side.text(
                0,
                y,
                f"期望下一代 ({row['expected_sensitive']:,.1f},"
                f" {row['expected_resistant']:,.1f})\n實際下一代"
                f" ({row['observed_next_sensitive']:,},"
                f" {row['observed_next_resistant']:,})",
                fontsize=12,
                va="top",
                color=GRAY,
            )
            y -= 0.12
    fig.text(
        0.07,
        0.115,
        r"箭頭端點 $=\phi(\mathrm{E}[S_{t+1}\mid s],\mathrm{E}[R_{t+1}\mid s])$，"
        r"$\phi(S,R)=(\log_{10}(S+1),\log_{10}(R+1))$；公式推導見矩陣流向圖（matrices_to_flow）。",
        fontsize=12.5,
    )
    fig.text(
        0.07,
        0.07,
        "族群越小，實際下一步偏離期望越多（人口學隨機性）。座標含 +1，"
        "因此零細胞可顯示；空心終點是計算上限，不是平衡點。",
        fontsize=11.5,
        color=GRAY,
    )
    fig.text(
        0.07,
        0.035,
        "相同 (S, R) 底下的 driver／TSG 組成可能不同，因此這是完整狀態的投影，"
        "不是封閉的二維向量場。",
        fontsize=11.5,
        color=GRAY,
    )
    return fig


# 將相圖儲存為 300 dpi PNG 與 SVG，然後關閉圖形。
def save_phase_figure(fig, destination, name):
    for extension in ("png", "svg"):
        fig.savefig(destination / f"{name}.{extension}", dpi=300, facecolor="white")
    plt.close(fig)


# ---------------------------------------------------------------- 理論 vs 經驗相場 --
THEORY, EMPIRICAL_CMAP, TRACE = "#2a78d6", "Oranges", "#C3CCD2"
RHO_TICKS = (0.001, 0.01, 0.1, 0.5, 0.9, 0.99, 0.999, 0.9999, 0.99999, 0.999999)
TEXT = {
    "zh": {
        "title": {
            "logN": "負荷 × 抗藥組成：理論截面與實際流向",
            "fitness": "適應度 × 抗藥組成：演化與選擇的分工",
        },
        "subtitle": (
            "上排：固定隱藏基因型的一步期望；下排：{t} 次實際試驗、"
            "起點同格的平均一步位移（≥{m} 次轉移、≥{k} 次試驗）。"
        ),
        "cols": ("未治療", "持續治療"),
        "rows": ("理論", "經驗"),
        "x": {
            "logN": r"總負荷 $\log_{10}(N+1)$",
            "fitness": r"平均無治療適應度 $\bar w_{\mathrm{off}}$（每細胞期望子代）",
        },
        "y": "抗藥比例 ρ（logit 刻度）",
        "S0": "S = 0 邊界",
        "R0": "R = 0 邊界",
        "legend": ("理論期望位移", "經驗平均位移（顏色＝轉移數）", "實際試驗路徑"),
        "count": "格內轉移數",
        "dx": {"logN": "Δ log₁₀ N", "fitness": "Δ w̄_off"},
        "slice": {
            "logN": "理論截面：敏感細胞全為基因型 0，抗藥細胞只帶抗藥位元。",
            "fitness": (
                "理論截面：N=治療門檻 {n:,}；oncogene 比例 f 由 w̄ 換算，與抗藥獨立，無 TSG 失活。"
            ),
        },
        "foot": (
            "治療中理論 logit ρ 每代上升 {d:.2f} 個 10 倍（log₁₀ w_R/w_S）；ρ→1"
            " 是組成上的穩定方向，但 log N 無不動點。"
        ),
        "foot2": (
            "經驗箭頭是 E[φ(X′)]−φ(X)，理論是 φ(E[X′])−φ(X)；大族群時兩者接近。"
            "計算上限不是平衡點。"
        ),
    },
    "en": {
        "title": {
            "logN": "Burden × resistant composition: theory slice vs observed flow",
            "fitness": "Fitness × resistant composition: evolution versus selection",
        },
        "subtitle": (
            "Top: one-step expectation on a fixed hidden-genotype slice. Bottom: mean"
            " one-step displacement of {t} trials from each cell (≥{m} transitions,"
            " ≥{k} trials)."
        ),
        "cols": ("No therapy", "Continuous therapy"),
        "rows": ("Theory", "Empirical"),
        "x": {
            "logN": r"Burden $\log_{10}(N+1)$",
            "fitness": (
                r"Mean no-therapy fitness $\bar w_{\mathrm{off}}$ (offspring per cell)"
            ),
        },
        "y": "Resistant fraction ρ (logit scale)",
        "S0": "S = 0 boundary",
        "R0": "R = 0 boundary",
        "legend": (
            "Theoretical expected step",
            "Empirical mean step (color = transitions)",
            "Observed trials",
        ),
        "count": "Transitions per cell",
        "dx": {"logN": "Δ log₁₀ N", "fitness": "Δ w̄_off"},
        "slice": {
            "logN": (
                "Theory slice: sensitive cells are genotype 0; resistant cells carry"
                " only the resistance bit."
            ),
            "fitness": (
                "Theory slice: N = therapy threshold {n:,}; oncogene fraction f from w̄,"
                " independent of resistance, no TSG loss."
            ),
        },
        "foot": (
            "On therapy, theoretical logit ρ rises {d:.2f} decades per generation"
            " (log₁₀ w_R/w_S): ρ→1 is a stable compositional direction, but log N has"
            " no fixed point."
        ),
        "foot2": (
            "Empirical arrows are E[φ(X′)]−φ(X); theory arrows are φ(E[X′])−φ(X); they"
            " agree for large populations. The computational cap is not an equilibrium."
        ),
    },
}


# 選取不大於輸入值的最大 1、2、5 乘十次方刻度，作為箭頭比例鍵。
def nice(value):
    """取 1、2、5 × 10^k 中不大於 value 的最大值，作為比例鍵。"""
    exponent = np.floor(np.log10(value))
    return max(m for m in (1, 2, 5) if m * 10**exponent <= value) * 10**exponent


# 將預設抗藥比例轉為 logit 刻度位置及百分比標籤。
def rho_ticks(lang):
    positions = [np.log10(p / (1 - p)) for p in RHO_TICKS]
    labels = [
        f"{100 * p:.{max(0, int(-np.floor(np.log10(min(p, 1 - p)))) - 2)}f}%"
        for p in RHO_TICKS
    ]
    return positions, labels


# 依同一治療狀態的理論及經驗位移，計算共用箭頭倍率與比例鍵。
def column_scale(fields, therapy, x_edges, y_edges):
    """同一治療狀態的理論與經驗共用比例：95 百分位箭頭 ≈ 0.9 格；長度仍與位移成正比。"""
    span = np.array([x_edges[-1] - x_edges[0], y_edges[-1] - y_edges[0]])
    cell = np.hypot(*(span / [len(x_edges) - 1, len(y_edges) - 1] / span))
    lengths = [
        np.hypot(f["dx"] / span[0], f["dy"] / span[1])
        for f in fields
        if f["therapy"] == therapy
    ]
    return 0.9 * cell / np.percentile(lengths, 95)


# 繪製指定相座標下未治療及治療中的理論與經驗向量場四面板比較圖。
def plot_pair(pair, rows, fields, steps, x_edges, y_edges, config, lang):
    t = TEXT[lang]
    fig, axes = plt.subplots(2, 2, figsize=(16.8, 13.2), sharex=True, sharey=True)
    fig.subplots_adjust(
        left=0.1, right=0.9, bottom=0.14, top=0.84, hspace=0.13, wspace=0.07
    )
    fig.text(0.05, 0.975, t["title"][pair], fontsize=24, weight="bold", va="top")
    norm = Normalize(
        MIN_TRANSITIONS,
        max(f["transitions"] for f in fields if f["kind"] == "empirical"),
    )
    cmap = plt.get_cmap(EMPIRICAL_CMAP)
    positions, labels = rho_ticks(lang)
    by_trial = {}
    for s in steps:
        by_trial.setdefault((s["replicate"], s["therapy"]), []).append(s)
    for c, therapy in enumerate((False, True)):
        k = column_scale(fields, therapy, x_edges, y_edges)
        for r, kind in enumerate(("theory", "empirical")):
            ax = axes[r, c]
            for segment in by_trial.values():
                if segment[0]["therapy"] == therapy:
                    xs = [s["x"] for s in segment] + [
                        segment[-1]["x"] + segment[-1]["dx"]
                    ]
                    ys = [s["y"] for s in segment] + [
                        segment[-1]["y"] + segment[-1]["dy"]
                    ]
                    ax.plot(xs, ys, color=TRACE, lw=0.7, alpha=0.55, zorder=1)
            chosen = [
                f for f in fields if f["kind"] == kind and f["therapy"] == therapy
            ]
            if chosen:
                X, Y, U, V = (
                    np.array([f[key] for f in chosen]) for key in ("x", "y", "dx", "dy")
                )
                color = (
                    THEORY
                    if kind == "theory"
                    else cmap(norm([f["transitions"] for f in chosen]) * 0.65 + 0.35)
                )
                ax.quiver(
                    X,
                    Y,
                    U * k,
                    V * k,
                    angles="xy",
                    scale_units="xy",
                    scale=1,
                    color=color,
                    width=0.0042,
                    headwidth=3.6,
                    headlength=4.2,
                    headaxislength=3.8,
                    zorder=4,
                )
            if pair == "logN":
                grid = np.linspace(x_edges[0], x_edges[-1], 100)
                total = 10**grid - 1
                ax.plot(grid, logit_rho(total, 0), color=GRAY, ls="--", lw=1)
                ax.plot(grid, logit_rho(0, total), color=GRAY, ls=":", lw=1)
                if r == 1 and c == 1:
                    ax.text(
                        grid[6],
                        logit_rho(10 ** grid[6] - 1, 0) + 0.3,
                        t["S0"],
                        color=GRAY,
                        fontsize=10.5,
                    )
                    ax.text(
                        grid[6],
                        logit_rho(0, 10 ** grid[6] - 1) + 0.3,
                        t["R0"],
                        color=GRAY,
                        fontsize=10.5,
                    )
            ax.set_xlim(x_edges[0], x_edges[-1])
            ax.set_ylim(y_edges[0], y_edges[-1])
            ax.set_yticks(positions, labels)
            ax.grid(color="#E6EBEE", lw=0.7)
            ax.set_axisbelow(True)
            if r == 0:
                ax.set_title(t["cols"][c], fontsize=16, weight="bold", pad=10)
            if c == 0:
                ax.set_ylabel(f"{t['rows'][r]}\n{t['y']}", fontsize=13)
            if r == 1:
                ax.set_xlabel(t["x"][pair], fontsize=13)
        # 箭頭比例鍵：x、y 各一支，取接近該欄位移 95 百分位的整齊數值；長度與資料位移成正比。
        chosen = [f for f in fields if f["therapy"] == therapy]
        ax = axes[0, c]
        span_x, span_y = x_edges[-1] - x_edges[0], y_edges[-1] - y_edges[0]
        x0, y0 = x_edges[0] + 0.05 * span_x, y_edges[-1] - 0.08 * span_y
        for axis, values in (
            ("x", [abs(f["dx"]) for f in chosen]),
            ("y", [abs(f["dy"]) for f in chosen]),
        ):
            value = nice(np.percentile(values, 95))
            u, v = (value * k, 0) if axis == "x" else (0, -value * k)
            ax.quiver(
                [x0],
                [y0],
                [u],
                [v],
                angles="xy",
                scale_units="xy",
                scale=1,
                color=INK,
                width=0.005,
                zorder=6,
            )
            label = (t["dx"][pair] if axis == "x" else "Δ logit ρ") + f" = {value:g}"
            where = (
                (x0 + u + 0.015 * span_x, y0)
                if axis == "x"
                else (x0 + 0.015 * span_x, y0 + v / 2)
            )
            ax.text(
                *where,
                label,
                fontsize=10.5,
                va="center",
                bbox={"facecolor": "white", "edgecolor": "none", "pad": 1.5},
                zorder=7,
            )
            if axis == "x":
                y0 -= 0.07 * span_y
    fig.legend(
        handles=[
            Line2D([], [], color=THEORY, lw=3, label=t["legend"][0]),
            Line2D([], [], color=cmap(0.8), lw=3, label=t["legend"][1]),
            Line2D([], [], color=TRACE, lw=2, label=t["legend"][2]),
        ],
        loc="upper left",
        bbox_to_anchor=(0.05, 0.905),
        ncol=3,
        frameon=False,
        fontsize=12,
    )
    bar = fig.colorbar(
        plt.cm.ScalarMappable(norm, cmap),
        ax=axes[1, :].tolist(),
        fraction=0.02,
        pad=0.015,
    )
    bar.set_label(t["count"])
    bar.outline.set_visible(False)
    return fig
