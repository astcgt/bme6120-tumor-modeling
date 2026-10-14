"""理論期望場 vs 經驗漂移場：兩組相座標，分開治療關閉與持續治療。

執行：.venv/bin/python plot_phase_comparison.py   （需先有 outputs/ 的治療試驗資料）
理論：固定隱藏基因型組成的截面，箭頭端點 = φ(原始計數一步期望)。
經驗：30 次實際試驗中，起點落在同一格的一步位移平均 E[φ(X_{t+1}) − φ(X_t)]。
兩者的差異，即被投影掉的自由度（driver/TSG 組成、隨機性）造成的效果。只用基礎模型。
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

from modeling_cancer import Config
from phase_state import (GRAY, INK, LANGS, ROOT, check_expectation, expected_next, fitness, logit_rho,
                         reproduce_trials, save, set_style, trial_feature_rows)
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D
import numpy as np


THEORY, EMPIRICAL_CMAP, TRACE = "#2a78d6", "Oranges", "#C3CCD2"
MIN_TRANSITIONS, MIN_TRIALS = 5, 2
RHO_TICKS = (0.001, 0.01, 0.1, 0.5, 0.9, 0.99, 0.999, 0.9999, 0.99999, 0.999999)
TEXT = {
    "zh": {"title": {"logN": "負荷 × 抗藥組成：理論截面與實際流向", "fitness": "適應度 × 抗藥組成：演化與選擇的分工"},
           "subtitle": "上排：固定隱藏基因型的一步期望；下排：{t} 次實際試驗、起點同格的平均一步位移（≥{m} 次轉移、≥{k} 次試驗）。",
           "cols": ("未治療", "持續治療"), "rows": ("理論", "經驗"),
           "x": {"logN": r"總負荷 $\log_{10}(N+1)$", "fitness": r"平均無治療適應度 $\bar w_{\mathrm{off}}$（每細胞期望子代）"},
           "y": "抗藥比例 ρ（logit 刻度）", "S0": "S = 0 邊界", "R0": "R = 0 邊界",
           "legend": ("理論期望位移", "經驗平均位移（顏色＝轉移數）", "實際試驗路徑"),
           "count": "格內轉移數", "dx": {"logN": "Δ log₁₀ N", "fitness": "Δ w̄_off"},
           "slice": {"logN": "理論截面：敏感細胞全為基因型 0，抗藥細胞只帶抗藥位元。",
                     "fitness": "理論截面：N=治療門檻 {n:,}；oncogene 比例 f 由 w̄ 換算，與抗藥獨立，無 TSG 失活。"},
           "foot": "治療中理論 logit ρ 每代上升 {d:.2f} 個 10 倍（log₁₀ w_R/w_S）；ρ→1 是組成上的穩定方向，但 log N 無不動點。",
           "foot2": "經驗箭頭是 E[φ(X′)]−φ(X)，理論是 φ(E[X′])−φ(X)；大族群時兩者接近。計算上限不是平衡點。"},
    "en": {"title": {"logN": "Burden × resistant composition: theory slice vs observed flow",
                     "fitness": "Fitness × resistant composition: evolution versus selection"},
           "subtitle": "Top: one-step expectation on a fixed hidden-genotype slice. Bottom: mean one-step displacement of {t} trials from each cell (≥{m} transitions, ≥{k} trials).",
           "cols": ("No therapy", "Continuous therapy"), "rows": ("Theory", "Empirical"),
           "x": {"logN": r"Burden $\log_{10}(N+1)$", "fitness": r"Mean no-therapy fitness $\bar w_{\mathrm{off}}$ (offspring per cell)"},
           "y": "Resistant fraction ρ (logit scale)", "S0": "S = 0 boundary", "R0": "R = 0 boundary",
           "legend": ("Theoretical expected step", "Empirical mean step (color = transitions)", "Observed trials"),
           "count": "Transitions per cell", "dx": {"logN": "Δ log₁₀ N", "fitness": "Δ w̄_off"},
           "slice": {"logN": "Theory slice: sensitive cells are genotype 0; resistant cells carry only the resistance bit.",
                     "fitness": "Theory slice: N = therapy threshold {n:,}; oncogene fraction f from w̄, independent of resistance, no TSG loss."},
           "foot": "On therapy, theoretical logit ρ rises {d:.2f} decades per generation (log₁₀ w_R/w_S): ρ→1 is a stable compositional direction, but log N has no fixed point.",
           "foot2": "Empirical arrows are E[φ(X′)]−φ(X); theory arrows are φ(E[X′])−φ(X); they agree for large populations. The computational cap is not an equilibrium."},
}


def nice(value):
    """取 1、2、5 × 10^k 中不大於 value 的最大值，作為比例鍵。"""
    exponent = np.floor(np.log10(value))
    return max(m for m in (1, 2, 5) if m * 10 ** exponent <= value) * 10 ** exponent


def rho_ticks(lang):
    positions = [np.log10(p / (1 - p)) for p in RHO_TICKS]
    labels = [f"{100 * p:.{max(0, int(-np.floor(np.log10(min(p, 1 - p)))) - 2)}f}%" for p in RHO_TICKS]
    return positions, labels


def theory_point(pair, x, y, config, therapy, size):
    """格點 → 完整族群（隱藏組成依截面固定）→ 一步期望 → 相座標位移。超出可行域時回傳 None。"""
    odds = 10 ** y
    if pair == "logN":
        total = 10 ** x - 1
        sensitive = (total + .5 - .5 * odds) / (1 + odds)
        if sensitive < 0 or sensitive > total:
            return None
        population = {0: sensitive, 1 << config.resistance_gene: total - sensitive}
    else:
        w0, w1 = fitness(0, config, False), fitness(1 << config.driver_genes[0], config, False)
        f = (x - w0) / (w1 - w0)
        sensitive = (size + .5 - .5 * odds) / (1 + odds)
        if not 0 <= f <= 1 or sensitive < 0 or sensitive > size:
            return None
        resistant, onc, res = size - sensitive, 1 << config.driver_genes[0], 1 << config.resistance_gene
        population = {0: sensitive * (1 - f), onc: sensitive * f, res: resistant * (1 - f), res | onc: resistant * f}
    nxt = expected_next(population, config, therapy)
    new_x = np.log10(nxt["N"] + 1) if pair == "logN" else nxt["w_off"]
    return new_x - x, float(logit_rho(nxt["R"], nxt["S"])) - y


def empirical_steps(rows, pair):
    """同一試驗相鄰世代的位移；治療狀態以起點世代為準（治療從該轉移起作用）。"""
    key = "log_N" if pair == "logN" else "mean_w_off"
    by_trial = {}
    for row in rows:
        by_trial.setdefault(row["replicate"], {})[row["generation"]] = row
    steps = []
    for replicate, gens in by_trial.items():
        for g, row in gens.items():
            if g + 1 in gens:
                nxt = gens[g + 1]
                steps.append({"replicate": replicate, "generation": g, "therapy": bool(row["therapy_on"]),
                              "x": row[key], "y": row["logit_R"],
                              "dx": nxt[key] - row[key], "dy": nxt["logit_R"] - row["logit_R"]})
    return steps


def build_fields(rows, config, pair):
    steps = empirical_steps(rows, pair)
    xs = np.array([s["x"] for s in steps]); ys = np.array([s["y"] for s in steps])
    pad_x = .04 * (xs.max() - xs.min())
    x_edges = np.linspace(xs.min() - pad_x, xs.max() + pad_x, 13)
    y_edges = np.linspace(np.floor(ys.min() * 2) / 2 - .25, np.ceil(ys.max() * 2) / 2 + .25, 15)
    xc, yc = (x_edges[:-1] + x_edges[1:]) / 2, (y_edges[:-1] + y_edges[1:]) / 2
    fields = []
    for therapy in (False, True):
        for i, x in enumerate(xc):
            for j, y in enumerate(yc):
                step = theory_point(pair, x, y, config, therapy, config.therapy_start_size)
                if step is not None:
                    fields.append({"pair": pair, "kind": "theory", "therapy": therapy, "x": x, "y": y,
                                   "dx": step[0], "dy": step[1], "transitions": "", "trials": ""})
        chosen = [s for s in steps if s["therapy"] == therapy]
        ix = np.digitize([s["x"] for s in chosen], x_edges) - 1
        iy = np.digitize([s["y"] for s in chosen], y_edges) - 1
        for i in range(len(xc)):
            for j in range(len(yc)):
                cell = [s for s, a, b in zip(chosen, ix, iy) if a == i and b == j]
                if len(cell) >= MIN_TRANSITIONS and len({s["replicate"] for s in cell}) >= MIN_TRIALS:
                    fields.append({"pair": pair, "kind": "empirical", "therapy": therapy,
                                   "x": np.mean([s["x"] for s in cell]), "y": np.mean([s["y"] for s in cell]),
                                   "dx": np.mean([s["dx"] for s in cell]), "dy": np.mean([s["dy"] for s in cell]),
                                   "transitions": len(cell), "trials": len({s["replicate"] for s in cell})})
    return fields, steps, x_edges, y_edges


def column_scale(fields, therapy, x_edges, y_edges):
    """同一治療狀態的理論與經驗共用比例：95 百分位箭頭 ≈ 0.9 格；長度仍與位移成正比。"""
    span = np.array([x_edges[-1] - x_edges[0], y_edges[-1] - y_edges[0]])
    cell = np.hypot(*(span / [len(x_edges) - 1, len(y_edges) - 1] / span))
    lengths = [np.hypot(f["dx"] / span[0], f["dy"] / span[1]) for f in fields if f["therapy"] == therapy]
    return .9 * cell / np.percentile(lengths, 95)


def plot_pair(pair, rows, fields, steps, x_edges, y_edges, config, lang):
    t = TEXT[lang]
    fig, axes = plt.subplots(2, 2, figsize=(16.8, 13.2), sharex=True, sharey=True)
    fig.subplots_adjust(left=.1, right=.9, bottom=.14, top=.84, hspace=.13, wspace=.07)
    fig.text(.05, .975, t["title"][pair], fontsize=24, weight="bold", va="top")
    fig.text(.05, .935, t["subtitle"].format(t=len({r["replicate"] for r in rows}), m=MIN_TRANSITIONS, k=MIN_TRIALS),
             fontsize=12.5, color=GRAY, va="top")
    norm = Normalize(MIN_TRANSITIONS, max(f["transitions"] for f in fields if f["kind"] == "empirical"))
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
                    xs = [s["x"] for s in segment] + [segment[-1]["x"] + segment[-1]["dx"]]
                    ys = [s["y"] for s in segment] + [segment[-1]["y"] + segment[-1]["dy"]]
                    ax.plot(xs, ys, color=TRACE, lw=.7, alpha=.55, zorder=1)
            chosen = [f for f in fields if f["kind"] == kind and f["therapy"] == therapy]
            if chosen:
                X, Y, U, V = (np.array([f[key] for f in chosen]) for key in ("x", "y", "dx", "dy"))
                color = THEORY if kind == "theory" else cmap(norm([f["transitions"] for f in chosen]) * .65 + .35)
                ax.quiver(X, Y, U * k, V * k, angles="xy", scale_units="xy", scale=1, color=color,
                          width=.0042, headwidth=3.6, headlength=4.2, headaxislength=3.8, zorder=4)
            if pair == "logN":
                grid = np.linspace(x_edges[0], x_edges[-1], 100)
                total = 10 ** grid - 1
                ax.plot(grid, logit_rho(total, 0), color=GRAY, ls="--", lw=1)
                ax.plot(grid, logit_rho(0, total), color=GRAY, ls=":", lw=1)
                if r == 1 and c == 1:
                    ax.text(grid[6], logit_rho(10 ** grid[6] - 1, 0) + .3, t["S0"], color=GRAY, fontsize=10.5)
                    ax.text(grid[6], logit_rho(0, 10 ** grid[6] - 1) + .3, t["R0"], color=GRAY, fontsize=10.5)
            ax.set_xlim(x_edges[0], x_edges[-1]); ax.set_ylim(y_edges[0], y_edges[-1])
            ax.set_yticks(positions, labels)
            ax.grid(color="#E6EBEE", lw=.7); ax.set_axisbelow(True)
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
        x0, y0 = x_edges[0] + .05 * span_x, y_edges[-1] - .08 * span_y
        for axis, values in (("x", [abs(f["dx"]) for f in chosen]), ("y", [abs(f["dy"]) for f in chosen])):
            value = nice(np.percentile(values, 95))
            u, v = (value * k, 0) if axis == "x" else (0, -value * k)
            ax.quiver([x0], [y0], [u], [v], angles="xy", scale_units="xy", scale=1, color=INK, width=.005, zorder=6)
            label = (t["dx"][pair] if axis == "x" else "Δ logit ρ") + f" = {value:g}"
            where = (x0 + u + .015 * span_x, y0) if axis == "x" else (x0 + .015 * span_x, y0 + v / 2)
            ax.text(*where, label, fontsize=10.5, va="center",
                    bbox={"facecolor": "white", "edgecolor": "none", "pad": 1.5}, zorder=7)
            if axis == "x":
                y0 -= .07 * span_y
    fig.legend(handles=[Line2D([], [], color=THEORY, lw=3, label=t["legend"][0]),
                        Line2D([], [], color=cmap(.8), lw=3, label=t["legend"][1]),
                        Line2D([], [], color=TRACE, lw=2, label=t["legend"][2])],
               loc="upper left", bbox_to_anchor=(.05, .905), ncol=3, frameon=False, fontsize=12)
    bar = fig.colorbar(plt.cm.ScalarMappable(norm, cmap), ax=axes[1, :].tolist(), fraction=.02, pad=.015)
    bar.set_label(t["count"]); bar.outline.set_visible(False)
    w_s = fitness(0, config, True); w_r = fitness(1 << config.resistance_gene, config, True)
    fig.text(.1, .075, t["slice"][pair].format(n=config.therapy_start_size), fontsize=11.5, color=GRAY)
    fig.text(.1, .05, t["foot"].format(d=np.log10(w_r / w_s)), fontsize=11.5)
    fig.text(.1, .025, t["foot2"], fontsize=11.5, color=GRAY)
    return fig


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=ROOT / "outputs")
    parser.add_argument("--output", type=Path, default=ROOT / "outputs/phase_plots")
    args = parser.parse_args()
    config, results, _ = reproduce_trials(args.input)
    check_expectation(config)
    rows = trial_feature_rows(config, results)
    all_fields, record = [], {}
    names = {"logN": "05_phase_field_burden_logit", "fitness": "06_phase_field_fitness_logit"}
    built = {pair: build_fields(rows, config, pair) for pair in names}
    for pair, (fields, steps, x_edges, y_edges) in built.items():
        all_fields += fields
        record[pair] = {"x_edges": x_edges.tolist(), "y_edges": y_edges.tolist(),
                        "theory_arrows": sum(f["kind"] == "theory" for f in fields),
                        "empirical_arrows": sum(f["kind"] == "empirical" for f in fields),
                        "transitions": len(steps)}
    args.output.mkdir(parents=True, exist_ok=True)
    with (args.output / "phase_field_comparison.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(all_fields[0])); writer.writeheader(); writer.writerows(all_fields)
    record.update({
        "y_coordinate": "empirical logit log10((R+0.5)/(S+0.5)); finite at rho=0 or 1",
        "empirical_rule": f">= {MIN_TRANSITIONS} transitions from >= {MIN_TRIALS} trials per cell; arrow from mean start",
        "theory_slices": {"logN": "sensitive=genotype 0, resistant=resistance bit only",
                          "fitness": f"N={config.therapy_start_size}; oncogene fraction from w_off, independent of resistance; no TSG loss"},
        "arrow_scale": "shared within each therapy column; length proportional to displacement; 95th percentile ≈ 0.9 cell",
        "base_model_only": True})
    (args.output / "phase_field_comparison.json").write_text(json.dumps(record, ensure_ascii=False, indent=2) + "\n")
    for lang in LANGS:
        set_style(lang)
        for pair, (fields, steps, x_edges, y_edges) in built.items():
            save(plot_pair(pair, rows, fields, steps, x_edges, y_edges, config, lang), args.output, names[pair], lang)
    print(json.dumps({k: {kk: v for kk, v in r.items() if "edges" not in kk} for k, r in record.items()
                      if isinstance(r, dict) and "transitions" in r}, indent=1))


if __name__ == "__main__":
    main()
