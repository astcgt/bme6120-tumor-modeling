"""主結果圖：腫瘤總覽（治療反應）、死亡率與異質性、one-hit vs two-hit。

資料由 experiments.py 計算；這裡只負責畫圖與存檔。
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import numpy as np

from experiments import HIT_SCENARIOS
from plot_style import COLORS, INK, save_figure, set_plot_style
from simulation import Config, Result, suppressor_alleles

HIT_SHORT = {
    "No TSG advantage": "Control\n(no TSG advantage)",
    "TSG one-hit": "One-hit",
    "TSG two-hit": "Two-hit",
}


# 腫瘤總覽：負荷；逐基因頻率改由基因頻率圖呈現。
def plot_therapy_ensemble(trials: list[Result], summary: list[dict], output: Path):
    """腫瘤總覽：負荷；逐基因頻率改由基因頻率圖呈現。"""
    length = len(summary)
    keys = ("total", "sensitive", "resistant")
    set_plot_style()
    fig, ax = plt.subplots(figsize=(10, 6.4), layout="constrained")
    axes = [ax]
    fig.suptitle(
        f"Tumor evolution under continuous therapy | {len(trials)} stochastic trials",
        fontsize=21,
        fontweight="bold",
    )
    x = np.arange(length)
    starts = [t.therapy_generation for t in trials if t.therapy_generation is not None]
    for key, color in zip(keys, ("#172B42", COLORS[0], COLORS[1])):
        axes[0].plot(
            x, [r[f"{key}_median"] for r in summary], label=key.title(), color=color
        )
        axes[0].fill_between(
            x,
            [r[f"{key}_p10"] for r in summary],
            [r[f"{key}_p90"] for r in summary],
            color=color,
            alpha=0.12,
        )
    axes[0].set(title="Population response", ylabel="Cells (symlog)")
    axes[0].set_yscale("symlog", linthresh=1)
    # 持續治療：從最早啟動的世代起一路到觀察結束都以淡色背景標示；起點因試驗而異。
    if starts:
        for ax in axes:
            ax.axvspan(
                min(starts),
                length - 1,
                color=COLORS[2],
                alpha=0.10,
                lw=0,
                label=f"Continuous therapy (starts gen {min(starts)}–{max(starts)})",
            )
    axes[0].legend(frameon=False, fontsize=9)
    for ax in axes:
        ax.set_xlabel("Generation")
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    save_figure(fig, output, "tumor_dashboard")


# 同樣 5,000 個細胞：死亡率越高越異質；B 顯示機制是需要更多代（更多分裂）才長到同樣大小。
def plot_death_rate(
    shannon: dict,
    tests: dict,
    config: Config,
    repeats: int,
    target: int,
    seed: int,
    output: Path,
):
    """同樣 5,000 個細胞：死亡率越高越異質；x 軸標出每代淨成長，說明需要更多代才長到同樣大小。"""
    set_plot_style()
    fig, ax = plt.subplots(figsize=(9, 6), layout="constrained")
    axes = [ax]
    fig.suptitle(
        "Higher cell turnover → more heterogeneity at the same tumor size",
        fontsize=17,
        fontweight="bold",
    )
    deaths = list(shannon)
    for ax, values, key in ((axes[0], shannon, "shannon"),):
        for index, death in enumerate(deaths):
            data = values[death]
            if not data:
                continue
            # 散點左右微幅偏移只改善重疊顯示；黑色橫線標示該組中位數。
            jitter = np.random.default_rng(seed + index).uniform(-0.10, 0.10, len(data))
            ax.scatter(index + jitter, data, color=COLORS[index], s=34, alpha=0.7)
            median = float(np.median(data))
            ax.plot([index - 0.18, index + 0.18], [median] * 2, color="#172B42", lw=3)
            ax.annotate(
                f"{median:.2f}" if key == "shannon" else f"{median:g} gen",
                (index + 0.22, median),
                va="center",
                fontsize=10.5,
                color="#172B42",
            )
            # 標示達到目標的試驗數，提醒統計是以達標為條件。
            ax.text(
                index,
                0.98,
                f"{len(data)}/{repeats} reached",
                ha="center",
                va="top",
                transform=ax.get_xaxis_transform(),
                fontsize=9,
            )
    # 無突變細胞的每代淨倍增 (1-d)(2+q)：死亡率越高，越接近 1，達標需要越多代。
    labels = [
        f"d = {d:.2f}\nm = ×{(1 - d) * (2 + config.triple_probability):.2f} / generation"
        for d in deaths
    ]
    axes[0].set(
        title=f"Diversity in {target:,}-cell samples",
        ylabel="Shannon diversity (natural log)",
    )
    for ax in axes:
        ax.set(
            xticks=range(len(deaths)),
            xticklabels=labels,
            xlabel=(
                "Death probability per cell per generation, d\n"
                "(m = expected offspring of a mutation-free cell = (1 − d)(2 + q))"
            ),
        )
        ax.margins(y=0.2, x=0.20)
        # 留出頂端空間給「達標數」文字與 p 值括號，避免壓到散點。
        low, high = ax.get_ylim()
        span = high - low
        ax.set_ylim(low, high + 0.42 * span)
    # 兩兩置換檢定：相鄰組在下排，跨組在上排。
    index = {d: i for i, d in enumerate(deaths)}
    for (a, b), p in tests.items():
        i, k = index[a], index[b]
        y = high + (0.04 if k - i == 1 else 0.17) * span
        top = y + 0.02 * span
        left, right = i + 0.04, k - 0.04
        axes[0].plot([left, left, right, right], [y, top, top, y], color=INK, lw=1)
        axes[0].text(
            (i + k) / 2,
            y + 0.03 * span,
            f"p = {p:.4f}" if p >= 1e-4 else "p < 0.0001",
            ha="center",
            va="bottom",
            fontsize=10,
        )
    save_figure(fig, output, "death_rate_heterogeneity")


# one-hit vs two-hit：對照組、one-hit、two-hit 的生長、TSG 達標比例與第 18 代負荷。
def plot_hit_comparison(hits: dict, config: Config, output: Path):
    """one-hit vs two-hit：對照組、one-hit、two-hit 的生長、TSG 達標比例與第 18 代負荷。"""
    # 圖 A 採線性尺度凸顯絕對數量差異；圖 B 顯示機制差異；圖 C 保留個別試驗。
    set_plot_style()
    output = Path(output)
    fig, axes = plt.subplots(1, 3, figsize=(17, 5.8), layout="constrained")
    last, repeats, order = hits["last"], hits["repeats"], hits["order"]
    fig.suptitle(
        "One-hit vs two-hit | tumor suppressor loss", fontsize=21, fontweight="bold"
    )
    styles = {
        HIT_SCENARIOS[0]: (COLORS[3], ":"),
        HIT_SCENARIOS[1]: (COLORS[1], "-"),
        HIT_SCENARIOS[2]: (COLORS[2], "--"),
    }
    medians = []
    for index, label in enumerate(order):
        color, linestyle = styles[label]
        total, fraction = hits["series"][label]
        for ax, values in zip(axes[:2], (total, fraction)):
            low, median, high = np.percentile(values, [10, 50, 90], axis=0)
            ax.plot(range(last + 1), median, color=color, label=label, ls=linestyle)
            ax.fill_between(
                range(last + 1),
                low,
                high,
                color=color,
                alpha=0.10 if index == 0 else 0.15,
            )
        # 固定終點比較避免拿不同世代的負荷相比；散點左右位移僅為減少遮蔽。
        endpoint = total[:, -1]
        jitter = np.random.default_rng(900 + index).uniform(-0.10, 0.10, repeats)
        axes[2].scatter(index + jitter, endpoint, color=color, alpha=0.7, s=40)
        median = float(np.median(endpoint))
        medians.append(median)
        axes[2].plot([index - 0.2, index + 0.2], [median] * 2, color=color, lw=3)
    # TSG 比例差距可能跨數量級；symlog 可同時保留零值和低頻 two-hit 細胞。
    axes[0].set(
        title="A  Tumor growth (linear scale)",
        xlabel="Generation",
        ylabel="Total cells",
    )
    axes[1].set(
        title="B  Cells meeting the TSG threshold",
        xlabel="Generation",
        ylabel="Cells with TSG loss (%)",
    )
    axes[1].set_yscale("symlog", linthresh=0.01)
    axes[1].text(
        0.03,
        0.97,
        "Symlog; linear below 0.01%",
        transform=axes[1].transAxes,
        va="top",
        fontsize=9,
    )
    axes[2].set(
        title=f"C  Tumor burden at generation {last}",
        ylabel="Total cells",
        xticks=range(len(order)),
        xticklabels=[
            f"{HIT_SHORT[s]}\nmedian {m:,.0f}" for s, m in zip(order, medians)
        ],
        xlim=(-0.5, len(order) - 0.5),
    )
    axes[2].margins(y=0.25)
    # 置換檢定結果：在 C 面板上方以括號標出兩兩比較的 p 值；相鄰組在下、跨組在上。
    low, high = axes[2].get_ylim()
    span = high - low
    axes[2].set_ylim(low, high + 0.32 * span)
    index = {label: i for i, label in enumerate(order)}
    pairs = sorted(hits["tests"], key=lambda ab: index[ab[1]] - index[ab[0]])
    for a, b in pairs:
        i, k = index[a], index[b]
        y = high + (0.04 if k - i == 1 else 0.17) * span
        top = y + 0.02 * span
        left, right = i + 0.04, k - 0.04
        axes[2].plot([left, left, right, right], [y, top, top, y], color=INK, lw=1)
        p = hits["tests"][(a, b)]
        axes[2].text(
            (i + k) / 2,
            y + 0.03 * span,
            f"p = {p:.4f}" if p >= 1e-4 else "p < 0.0001",
            ha="center",
            va="bottom",
            fontsize=10,
        )
    for ax in axes[:2]:
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))
        ax.legend(frameon=False, loc="lower right", fontsize=9)
    save_figure(fig, output, "one_hit_vs_two_hit")


# 輸出生長、基因型組成、多樣性與突變頻率四面板圖。
def plot_dashboard(result: Result, output: Path):
    """輸出生長、基因型組成、多樣性與突變頻率四面板圖。"""
    # 建立 2×2 面板，使用自動排版協調標題、座標標籤與色條的空間。
    set_plot_style()
    fig, axes = plt.subplots(2, 2, figsize=(15, 11), layout="constrained")
    fig.suptitle(
        "Tumor evolution | treatment & resistance", fontsize=21, fontweight="bold"
    )
    t = np.array([r["generation"] for r in result.history])
    # 面板 A：並列總細胞、敏感細胞與抗藥細胞的逐代數量。
    ax = axes[0, 0]
    for key, label, color in (
        ("total", "Total", "#172B42"),
        ("sensitive", "Sensitive", COLORS[0]),
        ("resistant", "Resistant", COLORS[1]),
    ):
        ax.plot(t, [r[key] for r in result.history], label=label, color=color)
    # symlog 同時呈現零細胞及跨數量級的增長，1 以下保留線性尺度。
    ax.set(
        yscale="symlog",
        ylabel="Cells (symlog; linear below 1)",
        title="A  Population response",
    )
    ax.set_yscale("symlog", linthresh=1)
    # 以垂直虛線標示治療起點，淡色背景表示持續接受治療的區間。
    if result.therapy_generation is not None:
        ax.axvspan(result.therapy_generation, t[-1], color=COLORS[1], alpha=0.09)
        ax.axvline(
            result.therapy_generation,
            color=COLORS[1],
            ls="--",
            lw=1.5,
            label="Therapy begins",
        )
    ax.legend(frameon=False, fontsize=9)
    # 面板 B：累加各世代的基因型比例，選出整段歷程中主要的五種基因型。
    ax = axes[0, 1]
    totals = np.array([r["total"] for r in result.history])
    scores = {}
    for pop, total in zip(result.populations, totals):
        for g, n in pop.items():
            scores[g] = scores.get(g, 0) + n / max(total, 1)
    top = sorted(scores, key=scores.get, reverse=True)[:5]
    # 將主要基因型轉為逐代比例，其餘合併為 Other；滅絕時總比例為零。
    bands = np.array([
        [p.get(g, 0) / max(n, 1) for p, n in zip(result.populations, totals)]
        for g in top
    ])
    others = np.maximum(0, (totals > 0).astype(float) - bands.sum(axis=0))
    # 以十六進位顯示基因型 bitmask；含抗藥基因者附上 R，方便辨識。
    labels = [
        f"G{g:06X}" + (" (R)" if g & (1 << result.config.resistance_gene) else "")
        for g in top
    ]
    ax.stackplot(t, *bands, others, labels=labels + ["Other"], colors=COLORS, alpha=0.9)
    ax.set(ylim=(0, 1), ylabel="Fraction of cells", title="B  Genotype composition")
    ax.legend(loc="upper left", fontsize=8, frameon=False, ncol=2)
    # 面板 C：exp(Shannon H) 將熵轉為有效基因型數，滅絕時另外定義為零。
    ax = axes[1, 0]
    ax.plot(
        t,
        [np.exp(r["shannon"]) if r["total"] else 0 for r in result.history],
        color=COLORS[2],
    )
    ax.set(
        ylabel="Effective genotypes: exp(Shannon H)", title="C  Genetic heterogeneity"
    )
    # 面板 D：TSG 顯示任一等位基因失活的細胞比例，不代表已達功能喪失門檻。
    ax = axes[1, 1]
    frequency = np.array([
        [
            sum(
                n
                for g, n in p.items()
                if (
                    suppressor_alleles(g, gene, result.config) > 0
                    if gene in result.config.suppressor_genes
                    else bool(g & (1 << gene))
                )
            )
            / max(total, 1)
            for p, total in zip(result.populations, totals)
        ]
        for gene in range(result.config.n_genes)
    ])
    # 用網格邊界繪製熱圖，固定色階 0–1，便於比較不同基因的突變頻率。
    mesh = ax.pcolormesh(
        np.arange(len(t) + 1) - 0.5,
        np.arange(result.config.n_genes + 1) - 0.5,
        frequency,
        cmap="YlGnBu",
        vmin=0,
        vmax=1,
        shading="flat",
    )
    # 顯示時將零起算索引改成 Gene 1 起算，並標出 driver 與 resistance。
    config = result.config

    # 依基因索引回傳 driver、抗藥或 TSG 的圖例文字，其他基因回傳空字串。
    def role(g: int) -> str:
        if g in config.driver_genes:
            return " · driver"
        if g == config.resistance_gene:
            return " · resistance"
        if g in config.suppressor_genes:
            return " · TSG (any hit)"
        return ""

    labels = [f"Gene {g + 1}{role(g)}" for g in range(config.n_genes)]
    ax.set(
        yticks=range(result.config.n_genes),
        yticklabels=labels,
        title="D  Mutation frequencies",
    )
    ax.tick_params(axis="y", labelsize=8)
    ax.grid(False)
    fig.colorbar(mesh, ax=ax, label="Mutant fraction", pad=0.02)
    # 所有面板共用世代標籤及整數刻度，最後一次輸出整張展示圖。
    for ax in axes.flat:
        ax.set_xlabel("Generation")
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    save_figure(fig, output, "tumor_dashboard")
