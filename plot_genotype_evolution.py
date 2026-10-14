"""基因型分布的演變：Muller 圖（Top-k 實際基因型）與每基因突變頻率軌跡。

執行：.venv/bin/python plot_genotype_evolution.py   （需先有 outputs/ 的治療試驗資料）
Muller 圖以「基因型巢狀」排列：子基因型畫在已呈現的最大子集基因型內。
同一基因型可能來自不同祖先，因此這不是重建的親緣樹。只用基礎模型。
"""
from __future__ import annotations

import argparse
import colorsys
import csv
import json
from pathlib import Path

from modeling_cancer import Config, Result, driver_count, suppressor_alleles, suppressor_loss
from phase_state import GRAY, INK, LANGS, ROOT, genotype_label, reproduce_trials, save, set_style
from plot_phase_portraits import choose_exemplar
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgb
from matplotlib.patches import Patch
import numpy as np


TOP_K = 10
# 功能類別決定色相，類別內以明度區分；身分另有直接標籤，不只靠顏色。
FAMILY = {"wt": "#9AA5AD", "onc": "#2a78d6", "tsg": "#1baf7a", "onc_tsg": "#4a3aa7",
          "res": "#eda100", "res_onc": "#eb6834", "res_tsg": "#e87ba4", "res_onc_tsg": "#e34948"}
TEXT = {
    "zh": {"title": "基因型組成的演變：治療篩選出預先存在的抗藥基因型",
           "subtitle": "replicate {r}（seed={s}）：每代前 {k} 個實際基因型（以峰值頻率選取）；其餘細胞併入最接近的已呈現祖先基因型。",
           "burden": r"負荷 $\log_{10}(N+1)$", "freq": "細胞比例", "gen": "世代",
           "therapy": "治療開始", "nadir": "谷底", "relapse": "復發",
           "families": {"wt": "無功能突變", "onc": "oncogene", "tsg": "TSG 失活", "onc_tsg": "oncogene+TSG",
                        "res": "抗藥", "res_onc": "抗藥+oncogene", "res_tsg": "抗藥+TSG", "res_onc_tsg": "抗藥+onc+TSG"},
           "genotypes": "基因型（主圖）", "family_title": "功能類別（色相）",
           "small": "其他試驗（顏色依功能類別，相同類別不代表相同基因型）",
           "rep": "rep {r}：復發 t={g}", "norelapse": "rep {r}：未復發",
           "foot": "巢狀規則：子基因型位於「已呈現且位元為其子集的最大基因型」內；這是基因型包含關係，不是重建的細胞譜系。",
           "foot2": "標籤：ONC=oncogene，TSGa/b=同一抑癌基因的兩個等位基因，RES=抗藥，P=passenger（基因編號由 1 起算）。",
           "freq_title": "每個基因的突變細胞比例（{n} 次試驗，以治療開始對齊）",
           "freq_sub": "線：中位數；帶：10–90 百分位（試驗間變異，非信賴區間）；只呈現所有試驗都有觀察的世代。",
           "panels": ("A  Oncogene", "B  抑癌基因（任一 / 雙等位基因）", "C  抗藥基因", "D  Passenger（中性對照）"),
           "rel": "距治療開始的世代", "mutfrac": "突變細胞比例",
           "any": "任一等位基因", "bi": "雙等位基因失活", "pmean": "14 個 passenger 平均",
           "freq_foot": "縱軸 symlog（10⁻⁴ 以下線性），保留 0；passenger 無適應度效果，用來對照純突變累積與搭便車（hitchhiking）。"},
    "en": {"title": "Genotype composition over time: therapy selects pre-existing resistant genotypes",
           "subtitle": "Replicate {r} (seed={s}): top {k} observed genotypes by peak frequency; all other cells merged into their nearest shown ancestor genotype.",
           "burden": r"Burden $\log_{10}(N+1)$", "freq": "Fraction of cells", "gen": "Generation",
           "therapy": "Therapy start", "nadir": "Nadir", "relapse": "Relapse",
           "families": {"wt": "No functional hit", "onc": "Oncogene", "tsg": "TSG loss", "onc_tsg": "Oncogene+TSG",
                        "res": "Resistant", "res_onc": "Resistant+oncogene", "res_tsg": "Resistant+TSG",
                        "res_onc_tsg": "Resistant+onc+TSG"},
           "genotypes": "Genotypes (main panel)", "family_title": "Functional class (hue)",
           "small": "Other trials (colors by functional class; same class ≠ same genotype)",
           "rep": "rep {r}: relapse t={g}", "norelapse": "rep {r}: no relapse",
           "foot": "Nesting: each genotype sits inside the largest shown genotype whose mutations are a subset of its own. This is genotype containment, not a reconstructed cell lineage.",
           "foot2": "Labels: ONC = oncogene, TSGa/b = two alleles of one tumor suppressor, RES = resistance, P = passenger (genes numbered from 1).",
           "freq_title": "Mutant-cell fraction per gene ({n} trials, aligned at therapy start)",
           "freq_sub": "Line: median; band: 10th–90th percentile (trial variability, not a confidence interval); only generations observed in every trial.",
           "panels": ("A  Oncogenes", "B  Tumor suppressors (any / biallelic)", "C  Resistance gene", "D  Passengers (neutral control)"),
           "rel": "Generations from therapy start", "mutfrac": "Mutant-cell fraction",
           "any": "any allele", "bi": "biallelic loss", "pmean": "mean of 14 passengers",
           "freq_foot": "Symlog y-axis (linear below 10⁻⁴) keeps zeros. Passengers have no fitness effect: a control for pure accumulation and hitchhiking."},
}


def family(genotype: int, config: Config) -> str:
    res = bool(genotype & (1 << config.resistance_gene))
    onc = driver_count(genotype, config) >= config.driver_hits
    tsg = suppressor_loss(genotype, config)
    name = "_".join(n for n, flag in (("onc", onc), ("tsg", tsg)) if flag) or "wt"
    return ("res_" + name if name != "wt" else "res") if res else name


def is_subset(a: int, b: int) -> bool:
    return a & b == a


def muller_layout(result: Result, k: int = TOP_K):
    """選 Top-k 基因型並建立巢狀樹；未選者併入最大的已選子集祖先，因此頻率總和恆為 1。"""
    peak = {}
    for population in result.populations:
        total = sum(population.values())
        for g, n in population.items():
            peak[g] = max(peak.get(g, 0), n / total)
    # 基因型 0（起始狀態）永遠是根，確保每個細胞都有可併入的祖先。
    shown = set(sorted(peak, key=lambda g: (-peak[g], g))[:k]) | {0}
    shown = sorted(shown, key=lambda g: (bin(g).count("1"), g))

    def ancestor(g):
        """已呈現、位元為 g 子集的最大基因型（不含 g 本身）。"""
        candidates = [a for a in shown if a != g and is_subset(a, g)]
        return max(candidates, key=lambda a: (bin(a).count("1"), -a))

    parent = {g: (None if g == 0 else ancestor(g)) for g in shown}
    children = {g: [c for c in shown if parent[c] == g] for g in shown}
    own = {g: np.zeros(len(result.populations)) for g in shown}
    for t, population in enumerate(result.populations):
        total = sum(population.values())
        for g, n in population.items():
            own[g if g in own else ancestor(g)][t] += n / total
    assert np.allclose(sum(own.values()), 1)
    return shown, parent, children, own


def muller_bands(shown, children, own):
    """遞迴：父基因型帶寬 = 自身 + 所有後代；子帶置於父帶中央，父色在上下兩側可見。"""
    subtree = {}
    def total(g):
        if g not in subtree:
            subtree[g] = own[g] + sum(total(c) for c in children[g])
        return subtree[g]
    total(0)
    bands = {}
    def place(g, bottom):
        bands[g] = (bottom, bottom + subtree[g])
        cursor = bottom + own[g] / 2
        for c in children[g]:
            place(c, cursor)
            cursor = cursor + subtree[c]
    place(0, np.zeros_like(own[0]))
    return bands


def colors_for(shown, config):
    """同一功能類別內依出現順序變化明度。"""
    counts, colors = {}, {}
    for g in shown:
        fam = family(g, config)
        i = counts.get(fam, 0); counts[fam] = i + 1
        h, l, s = colorsys.rgb_to_hls(*to_rgb(FAMILY[fam]))
        shift = (0, .16, -.12, .28, -.2, .36)[i % 6]
        colors[g] = colorsys.hls_to_rgb(h, float(np.clip(l + shift, .25, .88)), s)
    return colors


def draw_muller(ax, result, config, label=True, k=TOP_K):
    shown, parent, children, own = muller_layout(result, k)
    bands = muller_bands(shown, children, own)
    colors = colors_for(shown, config)
    t = np.arange(len(own[0]))
    order = sorted(shown, key=lambda g: (depth(g, parent), g))
    for g in order:
        low, high = bands[g]
        ax.fill_between(t, low, high, color=colors[g], linewidth=.4, edgecolor="white")
    if label:
        for g in shown:
            low, high = bands[g]
            width = own[g]
            best = int(np.argmax(width))
            if width[best] >= .06:
                mid = low[best] + (high[best] - low[best]) / 2
                text = genotype_label(g, config)
                edge = best > .85 * (len(t) - 1)
                ax.text(best - (.4 if edge else 0), mid, text if len(text) < 26 else text[:24] + "…",
                        ha="right" if edge else "center", va="center",
                        fontsize=9.5, color=INK, bbox={"facecolor": "white", "alpha": .75, "edgecolor": "none", "pad": 1.2})
    ax.set_xlim(0, len(t) - 1); ax.set_ylim(0, 1)
    return shown, own


def depth(g, parent):
    d = 0
    while parent[g] is not None:
        g = parent[g]; d += 1
    return d


def event_lines(ax, event, t, labels=True):
    for key, name, style in (("therapy_generation", "therapy", "--"), ("relapse_generation", "relapse", ":")):
        g = event[key]
        if g is not None:
            ax.axvline(g, color=INK, ls=style, lw=1.1)
            if labels:
                ax.text(g + .3, 1.02, t[name], fontsize=10.5, color=INK, transform=ax.get_xaxis_transform())


def plot_muller(config, results, events, lang):
    t = TEXT[lang]
    rep = choose_exemplar({e["replicate"]: e for e in events})
    result, event = results[rep], events[rep]
    fig = plt.figure(figsize=(16.8, 13))
    fig.text(.05, .975, t["title"], fontsize=22, weight="bold", va="top")
    fig.text(.05, .937, t["subtitle"].format(r=rep, s=result.seed, k=TOP_K), fontsize=12.5, color=GRAY, va="top")
    top = fig.add_axes((.07, .74, .64, .13))
    muller = fig.add_axes((.07, .38, .64, .34), sharex=top)
    totals = np.array([sum(p.values()) for p in result.populations])
    top.plot(np.log10(totals + 1), color=INK, lw=2)
    top.set_ylabel(t["burden"], fontsize=12); top.grid(color="#E6EBEE"); top.tick_params(labelbottom=False)
    event_lines(top, event, t)
    if event["therapy_generation"] is not None:
        after = range(event["therapy_generation"] + 1, len(totals))
        nadir = min(after, key=lambda g: totals[g])
        top.scatter([nadir], [np.log10(totals[nadir] + 1)], marker="v", color=INK, zorder=5)
        top.annotate(t["nadir"], (nadir, np.log10(totals[nadir] + 1)), xytext=(6, -4), textcoords="offset points", fontsize=10.5)
    shown, _ = draw_muller(muller, result, config)
    event_lines(muller, event, t, labels=False)
    muller.set_ylabel(t["freq"], fontsize=12); muller.set_xlabel(t["gen"], fontsize=12)
    muller.yaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{v:.0%}"))
    # 主圖：逐一列出 Top-k 基因型（色塊＋名稱），不只靠顏色辨識。
    colors = colors_for(shown, config)
    first = fig.legend(handles=[Patch(color=colors[g], label=genotype_label(g, config)) for g in shown],
                       title=t["genotypes"], loc="upper left", bbox_to_anchor=(.73, .87), frameon=False,
                       fontsize=11, title_fontsize=12, alignment="left")
    families = [f for f in FAMILY if any(family(g, config) == f for r in results for g in muller_layout(r)[0])]
    fig.legend(handles=[Patch(color=FAMILY[f], label=t["families"][f]) for f in families],
               title=t["family_title"], loc="upper left", bbox_to_anchor=(.86, .87), frameon=False,
               fontsize=11, title_fontsize=12, alignment="left")
    fig.add_artist(first)
    # 小分圖：最早與最晚復發，加上復發時間的 4 個分位數代表。
    relapsing = sorted((e for e in events if e["relapse_generation"] is not None and e["replicate"] != rep),
                       key=lambda e: (e["relapse_generation"], e["replicate"]))
    picks = [relapsing[int(round(q * (len(relapsing) - 1)))] for q in (0, .2, .4, .6, .8, 1)] if len(relapsing) >= 6 else relapsing
    others = [e for e in events if e["relapse_generation"] is None]
    picks = (picks + others)[:6]
    fig.text(.07, .315, t["small"], fontsize=13, weight="bold")
    for i, e in enumerate(picks):
        ax = fig.add_axes((.07 + i * .147, .1, .13, .19))
        draw_muller(ax, results[e["replicate"]], config, label=False)
        event_lines(ax, e, t, labels=False)
        ax.set_title(t["rep"].format(r=e["replicate"], g=e["relapse_generation"]) if e["relapse_generation"] is not None
                     else t["norelapse"].format(r=e["replicate"]), fontsize=10.5)
        ax.set_yticks([0, .5, 1], ["0%", "50%", "100%"] if i == 0 else [])
        ax.tick_params(labelsize=9.5)
    fig.text(.07, .045, t["foot"], fontsize=10.5, color=GRAY)
    fig.text(.07, .02, t["foot2"], fontsize=10.5, color=GRAY)
    return fig, rep, shown


def gene_frequencies(population, config):
    total = sum(population.values())
    freq = {}
    for gene in config.driver_genes + (config.resistance_gene,):
        freq[gene] = sum(n for g, n in population.items() if g & (1 << gene)) / total
    for gene in config.suppressor_genes:
        freq[f"{gene}_any"] = sum(n for g, n in population.items() if suppressor_alleles(g, gene, config)) / total
        freq[f"{gene}_bi"] = sum(n for g, n in population.items() if suppressor_alleles(g, gene, config) == 2) / total
    passengers = [g for g in range(config.n_genes) if g not in
                  config.driver_genes + config.suppressor_genes + (config.resistance_gene,)]
    freq["passenger_mean"] = np.mean([sum(n for g, n in population.items() if g & (1 << p)) / total for p in passengers])
    return freq


def plot_gene_frequencies(config, results, lang):
    """以治療開始對齊；只取所有試驗都存活且有觀察的相對世代。"""
    t = TEXT[lang]
    starts = [r.therapy_generation for r in results]
    lo = -min(starts); hi = min(len(r.populations) - 1 - s for r, s in zip(results, starts))
    rel = np.arange(lo, hi + 1)
    series = {}
    for r, s in zip(results, starts):
        for x in rel:
            for key, value in gene_frequencies(r.populations[s + x], config).items():
                series.setdefault(key, []).append(value)
    data = {key: np.array(v).reshape(len(results), len(rel)) for key, v in series.items()}
    fig, axes = plt.subplots(2, 2, figsize=(16.8, 11), sharex=True, sharey=True)
    fig.subplots_adjust(left=.08, right=.97, top=.84, bottom=.13, hspace=.22, wspace=.08)
    fig.text(.05, .965, t["freq_title"].format(n=len(results)), fontsize=22, weight="bold", va="top")
    fig.text(.05, .918, t["freq_sub"], fontsize=12.5, color=GRAY, va="top")
    palette = ("#2a78d6", "#eb6834", "#1baf7a", "#4a3aa7")
    groups = [[(g, f"ONC{i+1}", "-") for i, g in enumerate(config.driver_genes)],
              [(f"{g}_{kind}", f"TSG{i+1} {t[kind]}", "-" if kind == "any" else "--")
               for i, g in enumerate(config.suppressor_genes) for kind in ("any", "bi")],
              [(config.resistance_gene, "RES", "-")],
              [("passenger_mean", t["pmean"], "-")]]
    for ax, title, group in zip(axes.flat, t["panels"], groups):
        for j, (key, label, style) in enumerate(group):
            color = palette[(j // 2) if group is groups[1] else j]
            values = data[key]
            ax.fill_between(rel, np.percentile(values, 10, 0), np.percentile(values, 90, 0), color=color, alpha=.15, lw=0)
            ax.plot(rel, np.median(values, 0), color=color, ls=style, lw=2, label=label)
        ax.axvline(0, color=INK, ls="--", lw=1)
        ax.set_yscale("symlog", linthresh=1e-4); ax.set_ylim(0, 1.05)
        ax.set_yticks([0, 1e-4, 1e-3, 1e-2, .1, 1], ["0", "0.01%", "0.1%", "1%", "10%", "100%"])
        ax.set_title(title, loc="left", fontsize=15, weight="bold")
        ax.grid(color="#E6EBEE"); ax.legend(frameon=False, fontsize=11, loc="upper left")
    for ax in axes[1]:
        ax.set_xlabel(t["rel"], fontsize=12)
    for ax in axes[:, 0]:
        ax.set_ylabel(t["mutfrac"], fontsize=12)
    fig.text(.08, .05, t["freq_foot"], fontsize=11, color=GRAY)
    summary = {key: {"relative_generation": rel.tolist(), "median": np.median(v, 0).tolist(),
                     "p10": np.percentile(v, 10, 0).tolist(), "p90": np.percentile(v, 90, 0).tolist()}
               for key, v in data.items()}
    return fig, summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=ROOT / "outputs")
    parser.add_argument("--output", type=Path, default=ROOT / "outputs/genotype_evolution")
    args = parser.parse_args()
    config, results, events = reproduce_trials(args.input)
    args.output.mkdir(parents=True, exist_ok=True)
    record = {}
    for lang in LANGS:
        set_style(lang)
        fig, rep, shown = plot_muller(config, results, events, lang)
        save(fig, args.output, "01_muller_top_genotypes", lang)
        fig, summary = plot_gene_frequencies(config, results, lang)
        save(fig, args.output, "02_gene_mutation_frequencies", lang)
    _, _, _, own = muller_layout(results[rep])
    with (args.output / "muller_representative.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["generation"] + [genotype_label(g, config) for g in shown])
        for t in range(len(own[0])):
            writer.writerow([t] + [own[g][t] for g in shown])
    with (args.output / "gene_frequency_summary.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["series", "relative_generation", "median", "p10", "p90"])
        for key, s in summary.items():
            for i, x in enumerate(s["relative_generation"]):
                writer.writerow([key, x, s["median"][i], s["p10"][i], s["p90"][i]])
    record = {"representative_replicate": rep, "seed": results[rep].seed, "top_k": TOP_K,
              "shown_genotypes": [genotype_label(g, config) for g in shown],
              "csv_note": "muller_representative.csv stores each shown genotype's own fraction after merging unshown cells into their nearest shown ancestor",
              "nesting": "genotype containment (bit subset), not reconstructed lineage", "base_model_only": True}
    (args.output / "genotype_evolution.json").write_text(json.dumps(record, ensure_ascii=False, indent=2) + "\n")
    print(json.dumps(record, ensure_ascii=False, indent=1))


if __name__ == "__main__":
    main()
