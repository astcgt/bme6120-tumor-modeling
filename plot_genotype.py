"""基因型演化圖：Muller 圖（Top-k 實際基因型）與每基因突變頻率軌跡。

Muller 圖以「基因型巢狀」排列：子基因型畫在已呈現的最大子集基因型內。
同一基因型可能來自不同祖先，因此這不是重建的親緣樹。
巢狀版面與頻率由 analysis.py 計算；這裡只負責畫圖。
"""

from __future__ import annotations

import matplotlib.pyplot as plt
from matplotlib.colors import to_rgb
from matplotlib.patches import Patch
import numpy as np

from analysis import (
    TOP_K,
    choose_exemplar,
    depth,
    genotype_label,
    muller_bands,
    muller_layout,
)
from plot_style import GRAY, INK
from simulation import suppressor_alleles, suppressor_loss

# 依功能配色（Tableau 10 系列，論文常用、色盲與灰階下仍可分辨）：
# oncogene 暖色、TSG 綠色、兩者皆有為紫色、無功能突變為淺灰；抗藥另加斑點。
# 同一色組內依序混白或混黑，讓同類但不同的基因型仍好分辨。
ONCOGENE_COLORS = ("#E15759", "#F28E2B", "#EDC948")  # ONC1 珊瑚紅、ONC2 橘、ONC3 金黃
TSG_COLORS = ("#59A14F", "#2E9C84")  # TSG1 綠、TSG2 青綠
BOTH_COLOR = "#A374B5"  # oncogene + TSG：紫
NEUTRAL_COLOR = "#C3C9CF"  # 無功能突變：淺灰
SINGLE_ALLELE_PALE = 0.5  # TSG 只有一個等位基因突變：與白色混合的比例
SHADES = (
    (None, 0),
    ("white", 0.4),
    ("black", 0.28),
    ("white", 0.65),
    ("black", 0.45),
    ("white", 0.2),
    ("black", 0.15),
)
# 灰色組先變深：起始基因型 WT 淺灰，第二個（通常是只有 RES）明顯較深，兩者不靠斑點也分得開。
NEUTRAL_SHADES = (
    (None, 0),
    ("black", 0.3),
    ("white", 0.55),
    ("black", 0.5),
    ("white", 0.3),
    ("black", 0.15),
)
SPOTS = "."
TEXT = {
    "zh": {
        "title": "基因型組成的演變：治療篩選出預先存在的抗藥基因型",
        "subtitle": (
            "replicate {r}（seed={s}）：每代前 {k}"
            " 個實際基因型（以峰值頻率選取）；其餘細胞併入最接近的已呈現祖先基因型。"
        ),
        "burden": r"負荷 $\log_{10}(N+1)$",
        "freq": "細胞比例",
        "gen": "世代",
        "therapy": "治療開始",
        "nadir": "谷底",
        "relapse": "復發",
        "genotypes": "基因型（主圖，縮排 = 巢狀層級）",
        "key_title": "顏色規則",
        "key": {
            "onc": "暖色：oncogene\n（ONC1 珊瑚紅、ONC2 橘、ONC3 金黃）",
            "tsg": "綠色：TSG 失去功能\n（TSG1 綠、TSG2 青綠）",
            "tsg_one": "淡綠：TSG 只有一個等位基因突變",
            "both": "紫色：oncogene + TSG",
            "neutral": "灰色：無功能突變（WT、passenger）",
            "res": "斑點：帶抗藥突變（RES）",
        },
        "small": "其他試驗（顏色規則同主圖）",
        "rep": "rep {r}：復發 t={g}",
        "norelapse": "rep {r}：未復發",
        "foot": (
            "巢狀規則：子基因型位於「已呈現且位元為其子集的最大基因型」內（同大小時取該基因型出現時"
            "細胞較多者）；這是基因型包含關係，不是重建的細胞譜系。"
        ),
        "foot2": (
            "標籤：ONC=oncogene，TSGa/b=同一抑癌基因的兩個等位基因，RES=抗藥，"
            "P=passenger（基因編號由 1 起算）。"
        ),
        "freq_title": "每個基因的突變細胞比例（{n} 次試驗，以治療開始對齊）",
        "freq_sub": (
            "線：中位數；帶：10–90 百分位（試驗間變異，"
            "非信賴區間）；只呈現所有試驗都有觀察的世代。"
        ),
        "panels": (
            "A  Oncogene",
            "B  抑癌基因（任一 / 雙等位基因）",
            "C  抗藥基因",
            "D  Passenger（中性對照）",
        ),
        "rel": "距治療開始的世代",
        "mutfrac": "突變細胞比例",
        "any": "任一等位基因",
        "bi": "雙等位基因失活",
        "pmean": "14 個 passenger 平均",
        "freq_foot": (
            "縱軸 symlog（10⁻⁴ 以下線性），保留 0；passenger 無適應度效果，"
            "用來對照純突變累積與搭便車（hitchhiking）。"
        ),
    },
    "en": {
        "title": (
            "Genotype composition over time: therapy selects pre-existing resistant"
            " genotypes"
        ),
        "subtitle": (
            "Replicate {r} (seed={s}): top {k} observed genotypes by peak frequency;"
            " all other cells merged into their nearest shown ancestor genotype."
        ),
        "burden": r"Burden $\log_{10}(N+1)$",
        "freq": "Fraction of cells",
        "gen": "Generation",
        "therapy": "Therapy start",
        "nadir": "Nadir",
        "relapse": "Relapse",
        "genotypes": "Genotypes (main panel; indent = nesting)",
        "key_title": "Color code",
        "key": {
            "onc": "Warm: oncogene\n(ONC1 coral, ONC2 orange, ONC3 gold)",
            "tsg": "Green: TSG loss of function\n(TSG1 green, TSG2 teal green)",
            "tsg_one": "Pale green: one TSG allele hit",
            "both": "Purple: oncogene + TSG",
            "neutral": "Gray: no functional hit (WT, passengers)",
            "res": "Spotted: carries resistance (RES)",
        },
        "small": "Other trials (same color code as the main panel)",
        "rep": "rep {r}: relapse t={g}",
        "norelapse": "rep {r}: no relapse",
        "foot": (
            "Nesting: each genotype sits inside the largest shown genotype whose"
            " mutations are a subset of its own (ties: the larger one when it first"
            " appeared). This is genotype containment, not a reconstructed cell"
            " lineage."
        ),
        "foot2": (
            "Labels: ONC = oncogene, TSGa/b = two alleles of one tumor suppressor, RES"
            " = resistance, P = passenger (genes numbered from 1)."
        ),
        "freq_title": (
            "Mutant-cell fraction per gene ({n} trials, aligned at therapy start)"
        ),
        "freq_sub": (
            "Line: median; band: 10th–90th percentile (trial variability, not a"
            " confidence interval); only generations observed in every trial."
        ),
        "panels": (
            "A  Oncogenes",
            "B  Tumor suppressors (any / biallelic)",
            "C  Resistance gene",
            "D  Passengers (neutral control)",
        ),
        "rel": "Generations from therapy start",
        "mutfrac": "Mutant-cell fraction",
        "any": "any allele",
        "bi": "biallelic loss",
        "pmean": "mean of 14 passengers",
        "freq_foot": (
            "Symlog y-axis (linear below 10⁻⁴) keeps zeros. Passengers have no fitness"
            " effect: a control for pure accumulation and hitchhiking."
        ),
    },
}


def mix(color, target, amount):
    rgb = np.array(to_rgb(color))
    end = np.ones(3) if target == "white" else np.zeros(3)
    return tuple(rgb + (end - rgb) * amount)


def color_group(genotype, config):
    """(色組, 基準色)：只看 oncogene 與 TSG；抗藥與 passenger 不影響色組。"""
    oncogenes = [
        ONCOGENE_COLORS[i]
        for i, gene in enumerate(config.driver_genes)
        if genotype & (1 << gene)
    ]
    alleles = [
        suppressor_alleles(genotype, gene, config) for gene in config.suppressor_genes
    ]
    tsg = [TSG_COLORS[i] for i, n in enumerate(alleles) if n]
    lost = suppressor_loss(genotype, config)
    if oncogenes and tsg:
        return ("both", tuple(oncogenes + tsg)), to_rgb(BOTH_COLOR)
    if oncogenes:
        # 多個 oncogene 取平均色（例如 ONC2+ONC3 介於橘與金黃之間）。
        return ("onc", tuple(oncogenes)), tuple(
            np.mean([to_rgb(c) for c in oncogenes], axis=0)
        )
    if tsg:
        base = tuple(np.mean([to_rgb(c) for c in tsg], axis=0))
        # 失去功能（達 suppressor_hits）用原色；只有單一等位基因突變用淡綠。
        if not lost:
            base = mix(base, "white", SINGLE_ALLELE_PALE)
        return ("tsg", tuple(tsg), lost), base
    return ("neutral",), to_rgb(NEUTRAL_COLOR)


def genotype_colors(shown, config):
    """依功能決定基準色；同一色組內依出現順序混白或混黑，讓相鄰色帶好分辨。"""
    used, colors = {}, {}
    for g in shown:
        group, base = color_group(g, config)
        i = used.get(group, 0)
        used[group] = i + 1
        shades = NEUTRAL_SHADES if group == ("neutral",) else SHADES
        target, amount = shades[i % len(shades)]
        colors[g] = mix(base, target, amount) if target else base
    return colors


def spot_color(color):
    """深色帶用白點、淺色帶用深色點，斑點在任何底色上都看得見。"""
    r, g, b = color[:3]
    luminance = 0.2126 * r + 0.7152 * g + 0.0722 * b
    return (1, 1, 1, 0.8) if luminance < 0.6 else (0.1, 0.13, 0.17, 0.55)


def color_key(t):
    """顏色規則的圖例：每一類取代表色。"""
    samples = (
        ("onc", ONCOGENE_COLORS[0], None),
        ("tsg", TSG_COLORS[0], None),
        ("tsg_one", mix(TSG_COLORS[0], "white", SINGLE_ALLELE_PALE), None),
        ("both", BOTH_COLOR, None),
        ("neutral", NEUTRAL_COLOR, None),
        ("res", mix(NEUTRAL_COLOR, "black", 0.3), SPOTS),
    )
    return [spot_patch(color, spots, t["key"][name]) for name, color, spots in samples]


def spot_patch(color, spots, label):
    color = to_rgb(color)
    return Patch(
        facecolor=color,
        edgecolor=spot_color(color) if spots else color,
        hatch=spots,
        linewidth=0,  # 只留斑點，不畫外框
        label=label,
    )


def resistant(genotype, config):
    return bool(genotype & (1 << config.resistance_gene))


def draw_muller(ax, result, config, label=True, k=TOP_K):
    shown, parent, children, own = muller_layout(result, k)
    bands = muller_bands(shown, children, own)
    colors = genotype_colors(shown, config)
    t = np.arange(len(own[0]))
    for g in shown:  # 深度優先：父帶先畫，子帶疊在上面
        low, high = bands[g]
        ax.fill_between(t, low, high, color=colors[g], linewidth=0.4, edgecolor="white")
        if resistant(g, config):
            ax.fill_between(
                t,
                low,
                high,
                facecolor="none",
                edgecolor=spot_color(colors[g]),
                hatch=SPOTS,
                linewidth=0,
            )
    if label:
        for g in shown:
            low, high = bands[g]
            width = own[g]
            best = int(np.argmax(width))
            if width[best] >= 0.06:
                mid = low[best] + (high[best] - low[best]) / 2
                text = genotype_label(g, config)
                edge = best > 0.85 * (len(t) - 1)
                ax.text(
                    best - (0.4 if edge else 0),
                    mid,
                    text if len(text) < 26 else text[:24] + "…",
                    ha="right" if edge else "center",
                    va="center",
                    fontsize=9.5,
                    color=INK,
                    bbox={
                        "facecolor": "white",
                        "alpha": 0.8,
                        "edgecolor": "none",
                        "pad": 1.2,
                    },
                )
    ax.set_xlim(0, len(t) - 1)
    ax.set_ylim(0, 1)
    return shown, parent, colors


def event_lines(ax, event, t, labels=True):
    for key, name, style in (
        ("therapy_generation", "therapy", "--"),
        ("relapse_generation", "relapse", ":"),
    ):
        g = event[key]
        if g is not None:
            ax.axvline(g, color=INK, ls=style, lw=1.1)
            if labels:
                ax.text(
                    g + 0.3,
                    1.02,
                    t[name],
                    fontsize=10.5,
                    color=INK,
                    transform=ax.get_xaxis_transform(),
                )


def plot_muller(config, results, events, lang):
    t = TEXT[lang]
    rep = choose_exemplar({e["replicate"]: e for e in events})
    result, event = results[rep], events[rep]
    fig = plt.figure(figsize=(16.8, 13))
    fig.text(0.05, 0.975, t["title"], fontsize=22, weight="bold", va="top")
    fig.text(
        0.05,
        0.937,
        t["subtitle"].format(r=rep, s=result.seed, k=TOP_K),
        fontsize=12.5,
        color=GRAY,
        va="top",
    )
    top = fig.add_axes((0.07, 0.74, 0.64, 0.13))
    muller = fig.add_axes((0.07, 0.38, 0.64, 0.34), sharex=top)
    totals = np.array([sum(p.values()) for p in result.populations])
    top.plot(np.log10(totals + 1), color=INK, lw=2)
    top.set_ylabel(t["burden"], fontsize=12)
    top.grid(color="#E6EBEE")
    top.tick_params(labelbottom=False)
    event_lines(top, event, t)
    if event["therapy_generation"] is not None:
        after = range(event["therapy_generation"] + 1, len(totals))
        nadir = min(after, key=lambda g: totals[g])
        top.scatter(
            [nadir], [np.log10(totals[nadir] + 1)], marker="v", color=INK, zorder=5
        )
        top.annotate(
            t["nadir"],
            (nadir, np.log10(totals[nadir] + 1)),
            xytext=(6, -4),
            textcoords="offset points",
            fontsize=10.5,
        )
    shown, parent, colors = draw_muller(muller, result, config)
    event_lines(muller, event, t, labels=False)
    muller.set_ylabel(t["freq"], fontsize=12)
    muller.set_xlabel(t["gen"], fontsize=12)
    muller.yaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{v:.0%}"))
    # 主圖：依巢狀（深度優先）逐一列出 Top-k 基因型，縮排表示層級；抗藥者加斑點。
    handles = [
        spot_patch(
            colors[g],
            SPOTS if resistant(g, config) else None,
            "    " * depth(g, parent) + genotype_label(g, config),
        )
        for g in shown
    ]
    first = fig.legend(
        handles=handles,
        title=t["genotypes"],
        loc="upper left",
        bbox_to_anchor=(0.73, 0.87),
        frameon=False,
        fontsize=11,
        title_fontsize=12,
        alignment="left",
        handlelength=2.6,
        handleheight=1.3,
    )
    fig.legend(
        handles=color_key(t),
        title=t["key_title"],
        loc="upper left",
        bbox_to_anchor=(0.73, 0.5),
        frameon=False,
        fontsize=11,
        title_fontsize=12,
        alignment="left",
        handlelength=2.6,
        handleheight=1.3,
    )
    fig.add_artist(first)
    # 小分圖：最早與最晚復發，加上復發時間的 4 個分位數代表。
    relapsing = sorted(
        (
            e
            for e in events
            if e["relapse_generation"] is not None and e["replicate"] != rep
        ),
        key=lambda e: (e["relapse_generation"], e["replicate"]),
    )
    picks = (
        [
            relapsing[int(round(q * (len(relapsing) - 1)))]
            for q in (0, 0.2, 0.4, 0.6, 0.8, 1)
        ]
        if len(relapsing) >= 6
        else relapsing
    )
    others = [e for e in events if e["relapse_generation"] is None]
    picks = (picks + others)[:6]
    fig.text(0.07, 0.315, t["small"], fontsize=13, weight="bold")
    for i, e in enumerate(picks):
        ax = fig.add_axes((0.07 + i * 0.147, 0.1, 0.13, 0.19))
        draw_muller(ax, results[e["replicate"]], config, label=False)
        event_lines(ax, e, t, labels=False)
        ax.set_title(
            (
                t["rep"].format(r=e["replicate"], g=e["relapse_generation"])
                if e["relapse_generation"] is not None
                else t["norelapse"].format(r=e["replicate"])
            ),
            fontsize=10.5,
        )
        ax.set_yticks([0, 0.5, 1], ["0%", "50%", "100%"] if i == 0 else [])
        ax.tick_params(labelsize=9.5)
    fig.text(0.07, 0.045, t["foot"], fontsize=10.5, color=GRAY)
    fig.text(0.07, 0.02, t["foot2"], fontsize=10.5, color=GRAY)
    return fig, rep, shown


def plot_gene_frequencies(rel, data, config, n_trials: int, lang: str):
    """每個基因的突變細胞比例：中位數與 10–90 百分位帶，以治療開始對齊。"""
    t = TEXT[lang]
    fig, axes = plt.subplots(2, 2, figsize=(16.8, 11), sharex=True, sharey=True)
    fig.subplots_adjust(
        left=0.08, right=0.97, top=0.84, bottom=0.13, hspace=0.22, wspace=0.08
    )
    fig.text(
        0.05,
        0.965,
        t["freq_title"].format(n=n_trials),
        fontsize=22,
        weight="bold",
        va="top",
    )
    fig.text(0.05, 0.918, t["freq_sub"], fontsize=12.5, color=GRAY, va="top")
    palette = ("#2a78d6", "#eb6834", "#1baf7a", "#4a3aa7")
    groups = [
        [(g, f"ONC{i+1}", "-") for i, g in enumerate(config.driver_genes)],
        [
            (f"{g}_{kind}", f"TSG{i+1} {t[kind]}", "-" if kind == "any" else "--")
            for i, g in enumerate(config.suppressor_genes)
            for kind in ("any", "bi")
        ],
        [(config.resistance_gene, "RES", "-")],
        [("passenger_mean", t["pmean"], "-")],
    ]
    for ax, title, group in zip(axes.flat, t["panels"], groups):
        for j, (key, label, style) in enumerate(group):
            color = palette[(j // 2) if group is groups[1] else j]
            values = data[key]
            ax.fill_between(
                rel,
                np.percentile(values, 10, 0),
                np.percentile(values, 90, 0),
                color=color,
                alpha=0.15,
                lw=0,
            )
            ax.plot(rel, np.median(values, 0), color=color, ls=style, lw=2, label=label)
        ax.axvline(0, color=INK, ls="--", lw=1)
        ax.set_yscale("symlog", linthresh=1e-4)
        ax.set_ylim(0, 1.05)
        ax.set_yticks(
            [0, 1e-4, 1e-3, 1e-2, 0.1, 1], ["0", "0.01%", "0.1%", "1%", "10%", "100%"]
        )
        ax.set_title(title, loc="left", fontsize=15, weight="bold")
        ax.grid(color="#E6EBEE")
        ax.legend(frameon=False, fontsize=11, loc="upper left")
    for ax in axes[1]:
        ax.set_xlabel(t["rel"], fontsize=12)
    for ax in axes[:, 0]:
        ax.set_ylabel(t["mutfrac"], fontsize=12)
    fig.text(0.08, 0.05, t["freq_foot"], fontsize=11, color=GRAY)
    return fig
