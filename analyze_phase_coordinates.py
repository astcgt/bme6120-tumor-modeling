"""PCA 與「哪個因子決定走向」：變異最大的方向，不一定是決定下一步往哪走的方向。

執行：.venv/bin/python analyze_phase_coordinates.py
重跑記錄的治療試驗（逐代比對 CSV），只分析基礎模型；輸出至 outputs/phase_coordinates/。

做法：
1. PCA 只放生物因子（負荷、抗藥組成、oncogene、TSG、多樣性、治療時間），不放 w̄ 這類
   由因子合成的適應度，讓每個主成分都能讀成「因子的組合」。治療前、治療中分開做。
2. 每一步拆成「期望漂移」（由完整基因型組成精確算出）與「隨機噪音」。走向分兩個方向：
   負荷方向 Δln N、組成方向 Δlogit ρ。
3. 對期望漂移，比較每組因子單獨能解釋多少（R²）；另報告期望漂移占實際步伐的比例。
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

from phase_state import GRAY, INK, LANGS, ROOT, fitness, reproduce_trials, save, set_style, trial_feature_rows
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np


# 因子分組：組內彼此高度相關，組間代表不同生物機制。
GROUPS = {
    "burden": ("log_N",),
    "resistance": ("frac_R", "logit_R"),
    "oncogene": ("frac_onc", "mean_drivers"),
    "tsg": ("frac_tsg_anyhit", "frac_tsg_loss", "frac_onc_and_tsg"),
    "diversity": ("shannon", "simpson", "log_richness", "mean_passengers"),
    "time": ("gens_since_therapy",),
}
# PCA 的特徵：各組因子，再加上 log S、log R（負荷與組成的混合量，常被直接當座標）。
FEATURES = ("log_N", "log_S", "log_R", "frac_R", "logit_R", "frac_onc", "mean_drivers",
            "frac_tsg_anyhit", "frac_tsg_loss", "frac_onc_and_tsg", "shannon", "simpson",
            "log_richness", "mean_passengers", "gens_since_therapy")
GROUP_OF = {f: g for g, fs in GROUPS.items() for f in fs} | {"log_S": "mixed", "log_R": "mixed"}
GROUP_COLORS = {"burden": "#4b5563", "resistance": "#c2410c", "oncogene": "#1d4ed8", "tsg": "#7c3aed",
                "diversity": "#0f766e", "time": "#a16207", "mixed": "#64748b"}
SEGMENTS = ("pre", "on")
SEG_COLORS = {"pre": "#2a78d6", "on": "#eb6834"}
# 走向的兩個分量：(期望漂移欄, 實際一步欄)。
TARGETS = {"burden": ("growth_expected", "growth_next"), "composition": ("dlogit_expected", "dlogit_next")}

NAMES = {
    "zh": {"log_N": "總負荷 log N", "log_S": "敏感 log S", "log_R": "抗藥 log R", "frac_R": "抗藥比例 ρ",
           "logit_R": "logit ρ", "frac_onc": "oncogene 比例", "mean_drivers": "平均 driver 數",
           "frac_tsg_anyhit": "TSG 任一擊比例", "frac_tsg_loss": "TSG 失活比例",
           "frac_onc_and_tsg": "oncogene+TSG 比例", "shannon": "Shannon H", "simpson": "Simpson",
           "log_richness": "log 基因型數", "mean_passengers": "平均 passenger 數", "gens_since_therapy": "治療後代數"},
    "en": {"log_N": "Burden log N", "log_S": "Sensitive log S", "log_R": "Resistant log R",
           "frac_R": "Resistant fraction ρ", "logit_R": "logit ρ", "frac_onc": "Oncogene fraction",
           "mean_drivers": "Mean driver count", "frac_tsg_anyhit": "TSG any-hit fraction",
           "frac_tsg_loss": "TSG-loss fraction", "frac_onc_and_tsg": "Oncogene+TSG fraction",
           "shannon": "Shannon H", "simpson": "Simpson", "log_richness": "log genotype count",
           "mean_passengers": "Mean passenger count", "gens_since_therapy": "Generations on therapy"},
}
TEXT = {
    "zh": {"title": "PCA 找到變異最大的方向；決定腫瘤往哪走的是哪個因子？",
           "subtitle": "{n} 個（試驗×世代）狀態，{t} 次治療試驗。每一步拆成「期望漂移」（由完整基因型精確算出）＋「隨機噪音」。",
           "segments": {"pre": "治療前", "on": "治療中"},
           "groups": {"burden": "負荷", "resistance": "抗藥組成 ρ", "oncogene": "oncogene", "tsg": "TSG",
                      "diversity": "多樣性／passenger", "time": "治療時間"},
           "A": "A  變異最大 ≠ 決定走向", "A_bar": "解釋變異（長條）", "A_dot": "單一 PC 解釋負荷漂移的 R²（點）",
           "B": "B  主成分由哪些因子組成（負荷）",
           "C": "C  哪組因子解釋期望漂移（單獨 R²）",
           "C_cols": ("治療前\n負荷方向", "治療前\n組成方向", "治療中\n負荷方向", "治療中\n組成方向"),
           "C_signal": "漂移占實際步伐", "C_note": "框＝該方向最主導的因子組；灰＝該段為常數。\n多樣性、負荷也隨時間上升，所以 R² 偏高；只有框內因子在模型中直接改變 m。",
           "D1": "D  治療前：oncogene 決定漂移，但噪音更大",
           "D1x": "oncogene 活化比例", "D1y": "每代 ln 成長率",
           "D1_exp": "期望漂移", "D1_obs": "實際一步", "D1_theory": r"理論 $\ln[(1-f)w_0+f\,w_1]$",
           "D2": "E  治療中：抗藥比例 ρ 決定負荷漲跌",
           "D2x": "抗藥比例 ρ", "D2y": "每代期望 ln 成長率", "D2_theory": r"理論 $\ln[(1-\rho)w_S+\rho\,w_R]$",
           "D2_cbar": "oncogene 比例", "D2_star": "ρ* = {r:.0%}：由縮小轉為長大（谷底）",
           "key": "結論",
           "key_lines": ["PC1 是「隨時間一起變大」的量：變異最多。",
                         "   治療前它能預測漂移，只因 oncogene 也隨時間累積；",
                         "   治療中要靠 PC2（ρ）。",
                         "治療前：漂移完全由 oncogene 比例決定（R² {pre_onc:.2f}），",
                         "   但小族群噪音更大，實際一步只有 {pre_sig:.0%} 是漂移。",
                         "治療中：負荷漲跌由 ρ 決定（R² {on_res:.2f}），",
                         "   ρ 越過 {rstar:.0%} 即從縮小轉為長大。",
                         "組成方向兩段都由 ρ 自己決定：",
                         "   治療前是突變供應，治療中是篩選。"],
           "foot": "期望漂移：負荷 = ln(E[N′]/N)，組成 = logit(E[R′],E[S′]) − logit(R,S)，logit 為 log10((R+½)/(S+½))。組成方向只用 S > 0 的狀態。"
                   "R² 為各組因子的線性迴歸；PCA 使用標準化特徵。"},
    "en": {"title": "PCA finds the largest variance; which factor decides where the tumor goes?",
           "subtitle": "{n} (trial × generation) states from {t} therapy trials. Each step = expected drift (exact, from the full genotype state) + random noise.",
           "segments": {"pre": "Pre-therapy", "on": "On therapy"},
           "groups": {"burden": "Burden", "resistance": "Resistant fraction ρ", "oncogene": "Oncogene",
                      "tsg": "TSG", "diversity": "Diversity / passengers", "time": "Time on therapy"},
           "A": "A  Largest variance ≠ direction", "A_bar": "Variance explained (bars)",
           "A_dot": "R² of burden drift from that PC alone (dots)",
           "B": "B  What the PCs are made of (loadings)",
           "C": "C  Which factor group explains the drift (R² alone)",
           "C_cols": ("Pre-therapy\nburden", "Pre-therapy\ncomposition", "On therapy\nburden", "On therapy\ncomposition"),
           "C_signal": "Drift share of actual step", "C_note": "Box = dominant group per direction; gray = constant in segment.\nDiversity and burden also rise with time, inflating their R²; only boxed factors change m in the model.",
           "D1": "D  Pre-therapy: oncogenes set drift, noise is larger",
           "D1x": "Oncogene-active fraction", "D1y": "ln growth per generation",
           "D1_exp": "Expected drift", "D1_obs": "Actual step", "D1_theory": r"Theory $\ln[(1-f)w_0+f\,w_1]$",
           "D2": "E  On therapy: ρ decides grow vs shrink",
           "D2x": "Resistant fraction ρ", "D2y": "Expected ln growth per generation",
           "D2_theory": r"Theory $\ln[(1-\rho)w_S+\rho\,w_R]$", "D2_cbar": "Oncogene fraction",
           "D2_star": "ρ* = {r:.0%}: shrinking turns to growing (nadir)",
           "key": "Takeaways",
           "key_lines": ["PC1 = everything that grows with time (most variance).",
                         "   Tracks pre-therapy drift only via oncogene accumulation;",
                         "   on therapy, direction needs PC2 (ρ).",
                         "Pre-therapy: oncogene fraction sets drift (R² {pre_onc:.2f}),",
                         "   but noise dominates: drift = {pre_sig:.0%} of actual step.",
                         "On therapy: ρ sets growth vs decline (R² {on_res:.2f});",
                         "   crossing ρ = {rstar:.0%} turns shrinking into growing.",
                         "Composition is moved by ρ itself:",
                         "   mutation supply, then selection."],
           "foot": "Expected drift: burden = ln(E[N′]/N); composition = logit(E[R′],E[S′]) − logit(R,S), logit = log10((R+½)/(S+½)); composition uses states with S > 0. "
                   "R² from linear regression on each group; PCA on standardized features."},
}


def r_squared(X, y):
    ok = np.isfinite(y) & np.isfinite(X).all(1)
    if ok.sum() < 3 or np.all(X[ok].std(0) < 1e-12):
        return np.nan
    A = np.c_[np.ones(ok.sum()), X[ok]]
    residual = y[ok] - A @ np.linalg.lstsq(A, y[ok], rcond=None)[0]
    return float(1 - residual.var() / y[ok].var())


def pca(X):
    """標準化後 SVD；常數特徵負荷記為 NaN；每個主成分的符號固定為最大負荷為正。"""
    keep = X.std(0) > 1e-12
    Z = (X[:, keep] - X[:, keep].mean(0)) / X[:, keep].std(0)
    _, s, Vt = np.linalg.svd(Z, full_matrices=False)
    Vt = Vt * np.sign(Vt[np.arange(len(Vt)), np.abs(Vt).argmax(1)])[:, None]
    loadings = np.full((len(Vt), X.shape[1]), np.nan)
    loadings[:, keep] = Vt
    return s ** 2 / (s ** 2).sum(), loadings, Z @ Vt.T


def column(rows, key):
    return np.array([row[key] for row in rows], float)


def analyse(rows, config):
    summary = {"rows": len(rows), "features": list(FEATURES), "groups": {k: list(v) for k, v in GROUPS.items()},
               "segments": {}}
    for segment in SEGMENTS:
        sub = [r for r in rows if (r["therapy_on"] > 0) == (segment == "on")]
        X = np.array([[r[f] for f in FEATURES] for r in sub], float)
        variance, loadings, scores = pca(X)
        burden = column(sub, "growth_expected")
        record = {"n": len(sub), "explained_variance": variance[:5].tolist(),
                  "loadings": {f"PC{i+1}": {f: (None if np.isnan(v) else float(v)) for f, v in zip(FEATURES, loadings[i])}
                               for i in range(3)},
                  "burden_drift_r2_single_pc": [r_squared(scores[:, [k]], burden) for k in range(5)],
                  "directions": {}}
        for direction, (expected_key, actual_key) in TARGETS.items():
            # 組成方向只看 S > 0：S = 0 後經驗 logit 只是在記錄抗藥族群的生長。
            use = [r for r in sub if direction == "burden" or r["sensitive"] > 0]
            y, actual = column(use, expected_key), column(use, actual_key)
            record["directions"][direction] = {
                "n": len(use),
                "drift_mean": float(np.mean(y)), "drift_sd": float(np.std(y)),
                "drift_share_of_actual_r2": r_squared(y[:, None], actual),
                "group_r2_alone": {g: r_squared(np.array([[r[f] for f in fs] for r in use], float), y)
                                   for g, fs in GROUPS.items()},
            }
        summary["segments"][segment] = record
    # 理論：只有基準敏感與基準抗藥細胞時，治療中 ln 成長率 = 0 的抗藥比例。
    w_s, w_r = fitness(0, config, True), fitness(1 << config.resistance_gene, config, True)
    on = [r for r in rows if r["therapy_on"] > 0]
    rho, growth = column(on, "frac_R"), column(on, "growth_expected")
    summary["theory"] = {"w0_off": fitness(0, config, False), "w1_off": fitness(1 << config.driver_genes[0], config, False),
                         "wS_on": w_s, "wR_on": w_r, "rho_star": (1 - w_s) / (w_r - w_s)}
    summary["observed_sign_change"] = {"max_rho_with_decline": float(rho[growth < 0].max()),
                                       "min_rho_with_growth": float(rho[growth > 0].min())}
    # 期望負荷漂移必須等於當下平均適應度的 ln（IL11 關閉時突變不改變總數）。
    assert np.allclose(column(rows, "growth_expected"), np.log(column(rows, "mean_w_now")), atol=1e-12)
    summary["note"] = ("PCA ranks variance, not dynamical relevance. Drift = exact one-step conditional mean from the "
                       "full genotype state; R² from linear regression of drift on each factor group. Base model only.")
    return summary


def heat_text(ax, j, i, value, dark):
    ax.text(j, i, f"{value:.2f}", ha="center", va="center", fontsize=10.5, color="white" if dark else INK)


def plot(summary, rows, trials, lang):
    t = TEXT[lang]
    seg = summary["segments"]
    fig = plt.figure(figsize=(21, 13))
    fig.text(.04, .972, t["title"], fontsize=24, weight="bold", va="top")
    fig.text(.04, .932, t["subtitle"].format(n=summary["rows"], t=trials), fontsize=13.5, color=GRAY, va="top")

    # A：解釋變異（長條）與單一 PC 對負荷漂移的 R²（點）。
    ax = fig.add_axes((.05, .56, .25, .30))
    for i, s in enumerate(SEGMENTS):
        x = np.arange(5) + (i - .5) * .36
        ax.bar(x, seg[s]["explained_variance"], width=.34, color=SEG_COLORS[s], alpha=.35)
        ax.plot(x, seg[s]["burden_drift_r2_single_pc"], "o", color=SEG_COLORS[s], markersize=9,
                markeredgecolor="white", label=t["segments"][s])
    ax.set_xticks(range(5), [f"PC{i+1}" for i in range(5)]); ax.set_ylim(0, 1.02)
    ax.set_title(t["A"], loc="left", fontsize=15, weight="bold", pad=10)
    ax.grid(axis="y", color="#DDE4E8"); ax.set_axisbelow(True)
    ax.legend(frameon=False, fontsize=11, loc="upper right")
    ax.text(.99, .78, t["A_bar"] + "\n" + t["A_dot"], transform=ax.transAxes, ha="right", va="top", fontsize=10, color=GRAY)

    # B：負荷熱圖，列依因子組上色。
    ax = fig.add_axes((.44, .56, .2, .30))
    matrix = np.array([[seg[s]["loadings"][f"PC{k}"][f] for s in SEGMENTS for k in (1, 2, 3)] for f in FEATURES], float)
    cmap = plt.get_cmap("RdBu_r").copy(); cmap.set_bad("#E6EAED")
    image = ax.imshow(np.ma.masked_invalid(matrix), cmap=cmap, vmin=-.7, vmax=.7, aspect="auto")
    ax.set_xticks(range(6), [f"PC{k}" for _ in SEGMENTS for k in (1, 2, 3)], fontsize=10.5)
    ax.set_yticks(range(len(FEATURES)), [NAMES[lang][f] for f in FEATURES], fontsize=10.5)
    for tick, f in zip(ax.get_yticklabels(), FEATURES):
        tick.set_color(GROUP_COLORS[GROUP_OF[f]])
    for j, s in enumerate(SEGMENTS):
        ax.text(3 * j + 1, -1.1, t["segments"][s], ha="center", fontsize=12, weight="bold", color=SEG_COLORS[s])
    ax.axvline(2.5, color="white", lw=4)
    for i in range(len(FEATURES)):
        for j in range(6):
            if np.isfinite(matrix[i, j]) and abs(matrix[i, j]) >= .3:
                ax.text(j, i, f"{matrix[i, j]:+.2f}", ha="center", va="center", fontsize=8.5,
                        color="white" if abs(matrix[i, j]) > .5 else INK)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.tick_params(length=0); ax.grid(False)
    ax.set_title(t["B"], loc="left", fontsize=15, weight="bold", pad=28)
    bar = fig.colorbar(image, ax=ax, fraction=.06, pad=.03); bar.outline.set_visible(False)

    # 結論框。
    ax = fig.add_axes((.71, .56, .27, .30)); ax.axis("off")
    ax.add_patch(Rectangle((0, 0), 1, 1, transform=ax.transAxes, facecolor="#F4F7F9", edgecolor="#D6DEE2"))
    ax.text(.04, .93, t["key"], fontsize=16, weight="bold", va="top")
    values = dict(pre_onc=seg["pre"]["directions"]["burden"]["group_r2_alone"]["oncogene"],
                  pre_sig=seg["pre"]["directions"]["burden"]["drift_share_of_actual_r2"],
                  on_res=seg["on"]["directions"]["burden"]["group_r2_alone"]["resistance"],
                  rstar=summary["theory"]["rho_star"])
    for i, line in enumerate(t["key_lines"]):
        ax.text(.04, .83 - i * .088, line.format(**values), fontsize=12, va="top",
                weight="bold" if not line.startswith("   ") else "normal")

    # C：因子組 × 方向的 R² 熱圖，加一列「漂移占實際步伐」。
    ax = fig.add_axes((.12, .08, .2, .36))
    cols = [(s, d) for s in SEGMENTS for d in TARGETS]
    groups = list(GROUPS)
    grid = np.array([[seg[s]["directions"][d]["group_r2_alone"][g] for (s, d) in cols] for g in groups], float)
    signal = np.array([[seg[s]["directions"][d]["drift_share_of_actual_r2"] for (s, d) in cols]], float)
    ax.imshow(np.ma.masked_invalid(np.vstack([grid, np.full((1, 4), np.nan), signal])),
              cmap=cmap.from_list("r2", ["#FFFFFF", "#1d4ed8"]).with_extremes(bad="#FFFFFF"), vmin=0, vmax=1, aspect="auto")
    for i in range(len(groups)):
        for j in range(4):
            if np.isfinite(grid[i, j]):
                heat_text(ax, j, i, grid[i, j], grid[i, j] > .55)
            else:
                ax.add_patch(Rectangle((j - .5, i - .5), 1, 1, color="#E6EAED"))
    for j in range(4):
        heat_text(ax, j, len(groups) + 1, signal[0, j], signal[0, j] > .55)
        best = int(np.nanargmax(grid[:, j]))
        ax.add_patch(Rectangle((j - .5, best - .5), 1, 1, fill=False, edgecolor=INK, lw=2.4))
    ax.set_xticks(range(4), t["C_cols"], fontsize=10.5)
    for tick, (s, _) in zip(ax.get_xticklabels(), cols):
        tick.set_color(SEG_COLORS[s])
    ax.set_yticks(list(range(len(groups))) + [len(groups) + 1],
                  [t["groups"][g] for g in groups] + [t["C_signal"]], fontsize=11)
    for tick, g in zip(ax.get_yticklabels(), groups):
        tick.set_color(GROUP_COLORS[g])
    ax.get_yticklabels()[-1].set_weight("bold")
    ax.xaxis.tick_top(); ax.tick_params(length=0); ax.grid(False)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.axvline(1.5, color="white", lw=4)
    ax.set_title(t["C"], loc="left", fontsize=15, weight="bold", pad=52, x=-.45)
    ax.text(-.45, -.04, t["C_note"], transform=ax.transAxes, fontsize=10, color=GRAY, va="top")

    theory = summary["theory"]
    # D：治療前，負荷漂移 vs oncogene 比例；灰點為實際一步，顯示噪音。
    ax = fig.add_axes((.42, .08, .23, .34))
    pre = [r for r in rows if r["therapy_on"] == 0 and np.isfinite(r["growth_next"])]
    f, expected, actual = column(pre, "frac_onc"), column(pre, "growth_expected"), column(pre, "growth_next")
    ax.scatter(f, actual, s=10, color="#AAB6BE", alpha=.45, linewidths=0, label=t["D1_obs"])
    ax.scatter(f, expected, s=14, color=SEG_COLORS["pre"], alpha=.8, linewidths=0, label=t["D1_exp"])
    grid_f = np.linspace(0, max(f.max(), .05), 100)
    ax.plot(grid_f, np.log((1 - grid_f) * theory["w0_off"] + grid_f * theory["w1_off"]), color=INK, lw=1.6,
            label=t["D1_theory"])
    ax.set(xlabel=t["D1x"], ylabel=t["D1y"])
    ax.set_title(t["D1"], loc="left", fontsize=14, weight="bold", pad=10)
    ax.grid(color="#DDE4E8"); ax.set_axisbelow(True)
    ax.legend(frameon=False, fontsize=10, loc="upper right")

    # E：治療中，負荷漂移 vs ρ，點色＝oncogene 比例；虛線為理論零成長點。
    ax = fig.add_axes((.72, .08, .22, .34))
    on = [r for r in rows if r["therapy_on"] > 0]
    rho, growth, onc = column(on, "frac_R"), column(on, "growth_expected"), column(on, "frac_onc")
    points = ax.scatter(rho, growth, c=onc, cmap="viridis", s=16, vmin=0, vmax=1, linewidths=0, alpha=.85)
    grid_r = np.linspace(0, 1, 200)
    ax.plot(grid_r, np.log((1 - grid_r) * theory["wS_on"] + grid_r * theory["wR_on"]), color=INK, lw=1.6,
            label=t["D2_theory"])
    ax.axhline(0, color=GRAY, lw=1)
    ax.axvline(theory["rho_star"], color="#c2410c", ls="--", lw=1.4, label=t["D2_star"].format(r=theory["rho_star"]))
    ax.set(xlabel=t["D2x"], ylabel=t["D2y"], xlim=(0, 1.01))
    ax.set_title(t["D2"], loc="left", fontsize=14, weight="bold", pad=10)
    ax.grid(color="#DDE4E8"); ax.set_axisbelow(True)
    ax.legend(frameon=False, fontsize=10, loc="upper left")
    bar = fig.colorbar(points, ax=ax, fraction=.05, pad=.02); bar.set_label(t["D2_cbar"], fontsize=11)
    bar.outline.set_visible(False)

    fig.text(.04, .018, t["foot"], fontsize=10.5, color=GRAY)
    return fig


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=ROOT / "outputs")
    parser.add_argument("--output", type=Path, default=ROOT / "outputs/phase_coordinates")
    args = parser.parse_args()
    config, results, _ = reproduce_trials(args.input)
    rows = trial_feature_rows(config, results)
    summary = analyse(rows, config)
    args.output.mkdir(parents=True, exist_ok=True)
    with (args.output / "phase_features.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)
    (args.output / "pca_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n")
    for lang in LANGS:
        set_style(lang)
        save(plot(summary, rows, len(results), lang), args.output, "05_phase_coordinate_pca", lang)
    print(json.dumps({s: {"variance_pc1_3": [round(v, 3) for v in summary["segments"][s]["explained_variance"][:3]],
                          "burden_r2_single_pc": [round(v, 3) for v in summary["segments"][s]["burden_drift_r2_single_pc"]],
                          **{d: {"signal": round(r["drift_share_of_actual_r2"], 3),
                                 "groups": {g: (None if np.isnan(v) else round(v, 3)) for g, v in r["group_r2_alone"].items()}}
                             for d, r in summary["segments"][s]["directions"].items()}}
                      for s in SEGMENTS} | {"theory": summary["theory"], "observed": summary["observed_sign_change"]},
                     indent=1, ensure_ascii=False))


if __name__ == "__main__":
    main()
