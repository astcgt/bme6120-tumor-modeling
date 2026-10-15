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
from simulation import suppressor_alleles

# 依突變組成配色（CIELAB，感知均勻）：WT 為淺灰；每個功能基因朝自己的色向推一段固定距離，
# 五個色向在色環上平均分開，passenger 只推一小段；每多一個突變就變暗一些。
# 所以突變組成越相近，顏色越接近；差越多突變，顏色差越大。抗藥不影響顏色，只加斜線。
WT_LIGHTNESS = 84
DARKEN_PER_HIT = 11  # 每個功能突變降低的明度（passenger 與單一 TSG 等位基因按權重）
FUNCTIONAL_PUSH = 42  # 功能基因在 a*b* 平面上的推移距離
PASSENGER_PUSH = 14
PASSENGER_WEIGHT = 0.35
# CIELAB 色相角；避開中明度會變橄欖色的純黃（約 90–110°）。
ONCOGENE_DIRECTIONS = (28, 72, 160)  # ONC1 珊瑚紅、ONC2 琥珀橘、ONC3 青綠
TSG_DIRECTIONS = (250, 318)  # TSG1 藍、TSG2 紫
AMBER_LIFT = 6  # 帶 ONC2 時略提高明度，橘色才不會發灰
GOLDEN_ANGLE = 137.5  # passenger 色向依基因編號以黃金角分散，彼此不重疊
STRIPES = "///"
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
        "key_title": "顏色規則：每個突變一個色向，\n突變組成越近顏色越近、突變越多越深",
        "key_tsg": "TSG{i}（兩個等位基因）",
        "key_passenger": "passenger（只稍微偏移）",
        "key_res": "斜線：帶抗藥突變（RES）",
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
        "key_title": (
            "Color code: one direction per mutation;\nsimilar genotypes, similar"
            " colors; more hits, darker"
        ),
        "key_tsg": "TSG{i} (both alleles)",
        "key_passenger": "passenger (small shift)",
        "key_res": "Hatched: carries resistance (RES)",
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


def lab_to_rgb(lightness, a, b):
    """CIELAB（D65）→ sRGB，不在色域內就保持色相、逐步降低彩度。"""
    for scale in np.linspace(1, 0, 21):
        fy = (lightness + 16) / 116
        fx, fz = fy + a * scale / 500, fy - b * scale / 200

        def inverse(f):
            return f**3 if f**3 > 0.008856 else (f - 16 / 116) / 7.787

        xyz = np.array([0.95047 * inverse(fx), inverse(fy), 1.08883 * inverse(fz)])
        linear = XYZ_TO_LINEAR_RGB @ xyz
        if np.all(linear >= -1e-6) and np.all(linear <= 1 + 1e-6):
            break
    linear = np.clip(linear, 0, 1)
    rgb = np.where(
        linear <= 0.0031308, 12.92 * linear, 1.055 * linear ** (1 / 2.4) - 0.055
    )
    return tuple(float(x) for x in np.clip(rgb, 0, 1))


def push(angle, distance):
    radians = np.radians(angle)
    return np.array([np.cos(radians), np.sin(radians)]) * distance


def genotype_color(genotype, config):
    """WT 淺灰 + 每個突變的色向推移；突變越多越暗。抗藥位元不計入。"""
    ab, hits = np.zeros(2), 0.0
    for i, gene in enumerate(config.driver_genes):
        if genotype & (1 << gene):
            ab += push(ONCOGENE_DIRECTIONS[i], FUNCTIONAL_PUSH)
            hits += 1 - (AMBER_LIFT / DARKEN_PER_HIT if i == 1 else 0)
    for i, gene in enumerate(config.suppressor_genes):
        # 每個等位基因各推一半：單一等位基因突變時顏色介於 WT 與完全失活之間。
        alleles = suppressor_alleles(genotype, gene, config)
        ab += push(TSG_DIRECTIONS[i], FUNCTIONAL_PUSH * alleles / 2)
        hits += alleles / 2
    functional = set(config.driver_genes + config.suppressor_genes)
    for gene in range(config.n_genes):
        if gene in functional or gene == config.resistance_gene:
            continue
        if genotype & (1 << gene):
            ab += push(gene * GOLDEN_ANGLE, PASSENGER_PUSH)
            hits += PASSENGER_WEIGHT
    lightness = max(WT_LIGHTNESS - DARKEN_PER_HIT * hits, 36)
    return lab_to_rgb(lightness, *ab)


def stripe_color(color):
    """深色帶用白色斜線、淺色帶用深色斜線，在任何底色上都看得見。"""
    r, g, b = color[:3]
    luminance = 0.2126 * r + 0.7152 * g + 0.0722 * b
    return (1, 1, 1, 0.7) if luminance < 0.45 else (0.1, 0.13, 0.17, 0.32)


def swatch(color, stripes, label):
    color = to_rgb(color)
    return Patch(
        facecolor=color,
        edgecolor=stripe_color(color) if stripes else color,
        hatch=stripes,
        linewidth=0,  # 只留斜線，不畫外框
        label=label,
    )


def color_key(config, t):
    """各單一突變的顏色，加上一個抗藥（斜線）範例。"""
    examples = [
        (1 << gene, f"ONC{i + 1}") for i, gene in enumerate(config.driver_genes)
    ]
    for i, gene in enumerate(config.suppressor_genes):
        both = (1 << gene) | (1 << (config.n_genes + i))
        examples.append((both, t["key_tsg"].format(i=i + 1)))
    passenger = next(
        g
        for g in range(config.n_genes)
        if g not in config.driver_genes + config.suppressor_genes
        and g != config.resistance_gene
    )
    examples.append((1 << passenger, t["key_passenger"]))
    handles = [swatch(genotype_color(g, config), None, name) for g, name in examples]
    handles.append(swatch(genotype_color(0, config), STRIPES, t["key_res"]))
    return handles


def resistant(genotype, config):
    return bool(genotype & (1 << config.resistance_gene))


def draw_muller(ax, result, config, label=True, k=TOP_K):
    shown, parent, children, own = muller_layout(result, k)
    bands = muller_bands(shown, children, own)
    colors = {g: genotype_color(g, config) for g in shown}
    t = np.arange(len(own[0]))
    for g in shown:  # 深度優先：父帶先畫，子帶疊在上面
        low, high = bands[g]
        ax.fill_between(t, low, high, color=colors[g], linewidth=0.4, edgecolor="white")
        if resistant(g, config):
            with plt.rc_context({"hatch.linewidth": 0.7}):
                ax.fill_between(
                    t,
                    low,
                    high,
                    facecolor="none",
                    edgecolor=stripe_color(colors[g]),
                    hatch=STRIPES,
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
    # 主圖：依巢狀（深度優先）逐一列出 Top-k 基因型，縮排表示層級；抗藥者加斜線。
    handles = [
        swatch(
            colors[g],
            STRIPES if resistant(g, config) else None,
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
        handles=color_key(config, t),
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
