"""挑選相座標：候選狀態量的相關合併、PCA，以及對下一代成長率的解釋力。

執行：.venv/bin/python analyze_phase_coordinates.py
重跑記錄的治療試驗（逐代比對 CSV），只分析基礎模型；輸出至 outputs/phase_coordinates/。
PCA 找的是「變異最大」的方向；能否決定流向另以下一代 ln 成長率的 R² 檢驗。
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

from phase_state import GRAY, INK, LANGS, ROOT, reproduce_trials, save, set_style, trial_feature_rows
import matplotlib.pyplot as plt
import numpy as np


FEATURES = ("log_N", "log_S", "log_R", "frac_R", "logit_R", "frac_onc", "mean_drivers",
            "frac_tsg_anyhit", "frac_tsg_loss", "frac_onc_and_tsg", "mean_passengers",
            "log_richness", "shannon", "simpson", "mean_w_off", "var_w_off", "mean_w_on",
            "mean_w_now", "therapy_on", "gens_since_therapy")
NAMES = {
    "zh": {"log_N": "總負荷 log N", "log_S": "敏感 log S", "log_R": "抗藥 log R", "frac_R": "抗藥比例 ρ",
           "logit_R": "logit ρ", "frac_onc": "oncogene 活化比例", "mean_drivers": "平均 driver 數",
           "frac_tsg_anyhit": "TSG 任一擊比例", "frac_tsg_loss": "TSG 失活比例",
           "frac_onc_and_tsg": "oncogene+TSG 比例", "mean_passengers": "平均 passenger 數",
           "log_richness": "log 基因型數", "shannon": "Shannon H", "simpson": "Simpson",
           "mean_w_off": r"$\bar w_{\mathrm{off}}$（無治療適應度）", "var_w_off": "適應度變異數",
           "mean_w_on": r"$\bar w_{\mathrm{on}}$（治療下適應度）", "mean_w_now": "當下平均適應度",
           "therapy_on": "治療中", "gens_since_therapy": "治療後代數"},
    "en": {"log_N": "Burden log N", "log_S": "Sensitive log S", "log_R": "Resistant log R",
           "frac_R": "Resistant fraction ρ", "logit_R": "logit ρ", "frac_onc": "Oncogene-active fraction",
           "mean_drivers": "Mean driver count", "frac_tsg_anyhit": "TSG any-hit fraction",
           "frac_tsg_loss": "TSG-loss fraction", "frac_onc_and_tsg": "Oncogene+TSG fraction",
           "mean_passengers": "Mean passenger count", "log_richness": "log genotype count",
           "shannon": "Shannon H", "simpson": "Simpson", "mean_w_off": r"$\bar w_{\mathrm{off}}$ (no-therapy fitness)",
           "var_w_off": "Fitness variance", "mean_w_on": r"$\bar w_{\mathrm{on}}$ (on-therapy fitness)",
           "mean_w_now": "Current mean fitness", "therapy_on": "Therapy on",
           "gens_since_therapy": "Generations on therapy"},
}
TEXT = {
    "zh": {"title": "哪些自由度值得當相座標？", "subtitle": "{n} 個（試驗×世代）狀態，{t} 次治療試驗；標準化後 PCA。",
           "segments": {"all": "全部", "pre": "治療前", "on": "治療中"},
           "A": "A  前五個主成分的解釋變異", "B": "B  主成分負荷（|負荷| 越深越重要）",
           "C": "C  解釋下一代 ln 成長率的 R²", "pcs": "PC1–{k}", "variance": "解釋變異比例",
           "pairs": {"logN_logit": r"($\log N$, logit ρ)", "w_logit": r"($\bar w_{\mathrm{off}}$, logit ρ)"},
           "foot1": "PC1 主要是隨時間單調增加的量（負荷、多樣性、passenger），變異最大卻幾乎無法預測成長率。",
           "foot2": "治療中兩個選定座標即可解釋大部分成長率；治療前成長率主要是人口學隨機噪音。灰格：該段為常數。"},
    "en": {"title": "Which degrees of freedom make good phase coordinates?",
           "subtitle": "{n} (trial × generation) states from {t} therapy trials; PCA on standardized features.",
           "segments": {"all": "All", "pre": "Pre-therapy", "on": "On therapy"},
           "A": "A  Variance explained by the first five PCs", "B": "B  PC loadings (darker = larger |loading|)",
           "C": "C  R² for next-generation ln growth", "pcs": "PC1–{k}", "variance": "Fraction of variance",
           "pairs": {"logN_logit": r"($\log N$, logit ρ)", "w_logit": r"($\bar w_{\mathrm{off}}$, logit ρ)"},
           "foot1": "PC1 is mostly quantities that rise monotonically with time (burden, diversity, passengers): largest variance, almost no growth prediction.",
           "foot2": "On therapy the two chosen coordinates explain most growth; before therapy growth is dominated by demographic noise. Gray: constant in that segment."},
}
SEGMENTS = ("all", "pre", "on")


def standardized(X):
    keep = X.std(0) > 1e-12
    Z = np.full_like(X, np.nan)
    Z[:, keep] = (X[:, keep] - X[:, keep].mean(0)) / X[:, keep].std(0)
    return Z, keep


def pca(X):
    """SVD 主成分；每個主成分的符號固定為最大負荷為正，方便跨段比較。"""
    Z, keep = standardized(X)
    _, s, Vt = np.linalg.svd(Z[:, keep], full_matrices=False)
    signs = np.sign(Vt[np.arange(len(Vt)), np.abs(Vt).argmax(1)])
    Vt = Vt * signs[:, None]
    loadings = np.full((len(Vt), X.shape[1]), np.nan)
    loadings[:, keep] = Vt
    scores = Z[:, keep] @ Vt.T
    return s ** 2 / (s ** 2).sum(), loadings, scores


def r_squared(predictors, y):
    ok = np.isfinite(y)
    A = np.c_[np.ones(ok.sum()), predictors[ok]]
    residual = y[ok] - A @ np.linalg.lstsq(A, y[ok], rcond=None)[0]
    return float(1 - residual.var() / y[ok].var())


def redundant_groups(X, threshold=.95):
    C = np.corrcoef(X[:, X.std(0) > 0], rowvar=False)
    names = [f for f, k in zip(FEATURES, X.std(0) > 0) if k]
    left, groups = list(range(len(names))), []
    while left:
        i = left.pop(0)
        group = [i] + [j for j in left if abs(C[i, j]) >= threshold]
        left = [j for j in left if j not in group]
        if len(group) > 1:
            groups.append([names[j] for j in group])
    return groups


def analyse(rows):
    X = np.array([[row[f] for f in FEATURES] for row in rows], float)
    y = np.array([row["growth_next"] for row in rows], float)
    on = X[:, FEATURES.index("therapy_on")] > 0
    masks = {"all": np.ones(len(X), bool), "pre": ~on, "on": on}
    summary = {"rows": len(X), "redundant_groups_abs_r_ge_0.95": redundant_groups(X), "segments": {}}
    for name, mask in masks.items():
        variance, loadings, scores = pca(X[mask])
        pairs = {"logN_logit": ("log_N", "logit_R"), "w_logit": ("mean_w_off", "logit_R")}
        summary["segments"][name] = {
            "n": int(mask.sum()), "explained_variance": variance[:5].tolist(),
            "loadings": {f"PC{i+1}": {f: (None if np.isnan(v) else float(v)) for f, v in zip(FEATURES, loadings[i])}
                         for i in range(3)},
            "growth_r2_first_k_pcs": [r_squared(scores[:, :k], y[mask]) for k in (1, 2, 3)],
            "growth_r2_pairs": {key: r_squared(X[mask][:, [FEATURES.index(a), FEATURES.index(b)]], y[mask])
                                for key, (a, b) in pairs.items()},
            "growth_abs_corr": {f: float(abs(np.corrcoef(X[mask][np.isfinite(y[mask]), j],
                                                        y[mask][np.isfinite(y[mask])])[0, 1]))
                                for j, f in enumerate(FEATURES) if X[mask][:, j].std() > 0},
        }
    return summary


def plot(summary, trials, lang):
    t = TEXT[lang]
    fig = plt.figure(figsize=(16.8, 11))
    fig.text(.05, .965, t["title"], fontsize=24, weight="bold", va="top")
    fig.text(.05, .918, t["subtitle"].format(n=summary["rows"], t=trials), fontsize=13.5, color=GRAY, va="top")
    colors = {"all": "#687783", "pre": "#2a78d6", "on": "#eb6834"}
    # A：解釋變異。
    ax = fig.add_axes((.13, .56, .27, .28))
    for i, seg in enumerate(SEGMENTS):
        values = summary["segments"][seg]["explained_variance"]
        ax.bar(np.arange(5) + (i - 1) * .27, values, width=.25, color=colors[seg], label=t["segments"][seg])
    ax.set_xticks(range(5), [f"PC{i+1}" for i in range(5)])
    ax.set_ylabel(t["variance"]); ax.set_ylim(0, .8)
    ax.set_title(t["A"], loc="left", fontsize=15, weight="bold", pad=12)
    ax.legend(frameon=False, fontsize=11)
    ax.grid(axis="y", color="#DDE4E8"); ax.set_axisbelow(True)
    # C：成長率 R²。
    ax = fig.add_axes((.13, .12, .27, .3))
    labels = [t["pcs"].format(k=k) for k in (1, 2, 3)] + list(t["pairs"].values())
    for i, seg in enumerate(SEGMENTS):
        s = summary["segments"][seg]
        values = s["growth_r2_first_k_pcs"] + list(s["growth_r2_pairs"].values())
        ax.barh(np.arange(5) + (i - 1) * .27, values, height=.25, color=colors[seg])
    ax.set_yticks(range(5), labels); ax.invert_yaxis(); ax.set_xlim(0, 1)
    ax.set_title(t["C"], loc="left", fontsize=15, weight="bold", pad=12)
    ax.grid(axis="x", color="#DDE4E8"); ax.set_axisbelow(True)
    # B：負荷熱圖，欄為 三段 × PC1–3。
    ax = fig.add_axes((.56, .12, .38, .72))
    matrix = np.array([[summary["segments"][seg]["loadings"][f"PC{k}"][f] for seg in SEGMENTS for k in (1, 2, 3)]
                       for f in FEATURES], dtype=float)
    shown = np.ma.masked_invalid(matrix)
    cmap = plt.get_cmap("RdBu_r").copy(); cmap.set_bad("#E6EAED")
    image = ax.imshow(shown, cmap=cmap, vmin=-.75, vmax=.75, aspect="auto")
    ax.set_xticks(range(9), [f"PC{k}" for _ in SEGMENTS for k in (1, 2, 3)], fontsize=11)
    ax.set_yticks(range(len(FEATURES)), [NAMES[lang][f] for f in FEATURES], fontsize=11)
    for j, seg in enumerate(SEGMENTS):
        ax.text(3 * j + 1, -1.2, t["segments"][seg], ha="center", fontsize=12, weight="bold", color=colors[seg])
        if j:
            ax.axvline(3 * j - .5, color="white", lw=4)
    for i in range(len(FEATURES)):
        for j in range(9):
            if np.isfinite(matrix[i, j]) and abs(matrix[i, j]) >= .3:
                ax.text(j, i, f"{matrix[i, j]:+.2f}", ha="center", va="center", fontsize=8.5,
                        color="white" if abs(matrix[i, j]) > .5 else INK)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.tick_params(length=0)
    ax.set_title(t["B"], loc="left", fontsize=15, weight="bold", pad=34)
    bar = fig.colorbar(image, ax=ax, fraction=.03, pad=.02); bar.outline.set_visible(False)
    fig.text(.06, .045, t["foot1"], fontsize=11.5, color=GRAY)
    fig.text(.06, .018, t["foot2"], fontsize=11.5, color=GRAY)
    return fig


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=ROOT / "outputs")
    parser.add_argument("--output", type=Path, default=ROOT / "outputs/phase_coordinates")
    args = parser.parse_args()
    config, results, _ = reproduce_trials(args.input)
    rows = trial_feature_rows(config, results)
    summary = analyse(rows)
    summary["note"] = ("PCA ranks variance, not dynamical relevance; growth R^2 uses next-generation ln(N_t+1/N_t). "
                       "Features are projections of the full genotype state; base model only (IL11 off).")
    args.output.mkdir(parents=True, exist_ok=True)
    with (args.output / "phase_features.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)
    (args.output / "pca_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n")
    for lang in LANGS:
        set_style(lang)
        save(plot(summary, len(results), lang), args.output, "05_phase_coordinate_pca", lang)
    print(json.dumps({s: {"variance_pc1_3": [round(v, 3) for v in summary["segments"][s]["explained_variance"][:3]],
                          "r2_pcs": [round(v, 3) for v in summary["segments"][s]["growth_r2_first_k_pcs"]],
                          "r2_pairs": {k: round(v, 3) for k, v in summary["segments"][s]["growth_r2_pairs"].items()}}
                      for s in SEGMENTS}, indent=1))


if __name__ == "__main__":
    main()
