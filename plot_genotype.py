"""基因型演化圖：Muller 圖（Top-k 實際基因型）與每基因突變頻率軌跡。

Muller 圖以「基因型巢狀」排列：子基因型畫在已呈現的最大子集基因型內。
同一基因型可能來自不同祖先，因此這不是重建的親緣樹。
巢狀版面與頻率由 analysis.py 計算；這裡只負責畫圖。
"""

from __future__ import annotations

import matplotlib.pyplot as plt
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

# 依功能配色（CIE LCh 色相角）：oncogene 暖色、TSG 冷色、兩者皆有為紫色、無功能突變為灰；
# 抗藥另加白色斜線。同一色組內再以明度與小幅色相差區分不同基因型。
ONCOGENE_HUES = (25, 55, 85)  # ONC1 紅、ONC2 橘、ONC3 金
TSG_HUES = (250, 205)  # TSG1 藍、TSG2 青
BOTH_HUE = 315  # oncogene + TSG：紫
NEUTRAL_HUE = 250  # 灰色帶一點冷調，與暖色 oncogene 分開
GROUP_LIGHTNESS = (60, 74, 47, 84, 40, 67, 53)  # 同色組內依序使用，相鄰兩者明暗差大
NEUTRAL_LIGHTNESS = (72, 56, 84, 44, 64, 50)
HUE_OFFSETS = (0, 7, -7, 12, -12)
HATCH = "////"
# CIE XYZ（D65）→ 線性 sRGB 的標準矩陣。
XYZ_TO_LINEAR_RGB = np.array(
    [
        [3.2406, -1.5372, -0.4986],
        [-0.9689, 1.8758, 0.0415],
        [0.0557, -0.2040, 1.0570],
    ]
)
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
            "onc": "暖色：oncogene\n（ONC1 紅、ONC2 橘、ONC3 金）",
            "tsg": "冷色：TSG 失去功能\n（TSG1 藍、TSG2 青）",
            "tsg_one": "淡冷色：TSG 只有一個等位基因突變",
            "both": "紫色：oncogene + TSG",
            "neutral": "灰色：無功能突變（WT、passenger）",
            "res": "斜線：帶抗藥突變（RES）",
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
            "onc": "Warm: oncogene\n(ONC1 red, ONC2 orange, ONC3 gold)",
            "tsg": "Cool: TSG loss of function\n(TSG1 blue, TSG2 teal)",
            "tsg_one": "Pale cool: one TSG allele hit",
            "both": "Purple: oncogene + TSG",
            "neutral": "Gray: no functional hit (WT, passengers)",
            "res": "Hatched: carries resistance (RES)",
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


def lch_to_rgb(lightness, chroma, hue):
    """CIE LCh（D65）→ sRGB；同明度的顏色看起來一樣亮，避免黃色特別刺眼。"""
    h = np.radians(hue)
    a, b = chroma * np.cos(h), chroma * np.sin(h)
    fy = (lightness + 16) / 116
    fx, fz = fy + a / 500, fy - b / 200

    def inverse(f):
        return f**3 if f**3 > 0.008856 else (f - 16 / 116) / 7.787

    xyz = np.array([0.95047 * inverse(fx), inverse(fy), 1.08883 * inverse(fz)])
    linear = np.clip(XYZ_TO_LINEAR_RGB @ xyz, 0, 1)
    rgb = np.where(
        linear <= 0.0031308, 12.92 * linear, 1.055 * linear ** (1 / 2.4) - 0.055
    )
    return tuple(np.clip(rgb, 0, 1))


def mean_hue(hues):
    angles = np.radians(hues)
    return float(np.degrees(np.arctan2(np.sin(angles).mean(), np.cos(angles).mean())))


def color_group(genotype, config):
    """(色組, 色相, 彩度)：只看 oncogene 與 TSG；抗藥與 passenger 不影響色組。"""
    oncogenes = [
        ONCOGENE_HUES[i]
        for i, gene in enumerate(config.driver_genes)
        if genotype & (1 << gene)
    ]
    alleles = [
        suppressor_alleles(genotype, gene, config) for gene in config.suppressor_genes
    ]
    tsg = [TSG_HUES[i] for i, n in enumerate(alleles) if n]
    lost = suppressor_loss(genotype, config)
    if oncogenes and tsg:
        return ("both", tuple(oncogenes + tsg)), BOTH_HUE, 50
    if oncogenes:
        return ("onc", tuple(oncogenes)), mean_hue(oncogenes), 62
    if tsg:
        # 失去功能（達 suppressor_hits）用飽和冷色；只有單一等位基因突變用淡冷色。
        return ("tsg", tuple(tsg), lost), mean_hue(tsg), 50 if lost else 24
    return ("neutral",), NEUTRAL_HUE, 6


def genotype_colors(shown, config):
    """依功能決定色相；同一色組內依出現順序換明度並微調色相，讓相鄰色帶好分辨。"""
    used, colors = {}, {}
    for g in shown:
        group, hue, chroma = color_group(g, config)
        i = used.get(group, 0)
        used[group] = i + 1
        levels = NEUTRAL_LIGHTNESS if group == ("neutral",) else GROUP_LIGHTNESS
        lightness = levels[i % len(levels)]
        hue = (hue + HUE_OFFSETS[i % len(HUE_OFFSETS)]) % 360
        if 60 <= hue <= 110:
            lightness = max(lightness, 70)  # 中低明度的黃色會變成橄欖綠
        colors[g] = lch_to_rgb(lightness, chroma, hue)
    return colors


def color_key(t):
    """顏色規則的圖例：每一類取代表色。"""
    samples = (
        ("onc", lch_to_rgb(62, 62, ONCOGENE_HUES[1]), None),
        ("tsg", lch_to_rgb(60, 50, TSG_HUES[0]), None),
        ("tsg_one", lch_to_rgb(74, 24, TSG_HUES[0]), None),
        ("both", lch_to_rgb(60, 50, BOTH_HUE), None),
        ("neutral", lch_to_rgb(72, 6, NEUTRAL_HUE), None),
        ("res", lch_to_rgb(56, 6, NEUTRAL_HUE), HATCH),
    )
    return [
        Patch(
            facecolor=color,
            edgecolor="white" if hatch else color,
            hatch=hatch,
            label=t["key"][name],
        )
        for name, color, hatch in samples
    ]


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
            with plt.rc_context({"hatch.linewidth": 0.5}):
                ax.fill_between(
                    t,
                    low,
                    high,
                    facecolor="none",
                    edgecolor=(1, 1, 1, 0.55),
                    hatch=HATCH,
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
    # 主圖：依巢狀（深度優先）逐一列出 Top-k 基因型，縮排表示層級；抗藥者畫斜線。
    handles = [
        Patch(
            facecolor=colors[g],
            edgecolor=(1, 1, 1, 0.8) if resistant(g, config) else colors[g],
            hatch=HATCH if resistant(g, config) else None,
            label="    " * depth(g, parent) + genotype_label(g, config),
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
