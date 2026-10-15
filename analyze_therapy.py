"""治療延伸分析：何時治療、何時復發、參數多敏感、間歇治療能否延後惡化。

執行：.venv/bin/python analyze_therapy.py            （全部四項，輸出至 outputs/therapy_analysis/）
      .venv/bin/python analyze_therapy.py --only 08  （只跑其中一項）

08 Goldie–Coldman：治療起始負荷 vs 治癒機率，並與分枝過程理論（突變 clone 的 Poisson 近似）比較。
09 復發時間：用治療開始時的抗藥細胞數預測每次試驗的復發時間；ρ* 決定谷底。
10 參數敏感度：抗藥細胞殺傷率 k_R × 突變率 μ 的治癒／復發相圖。
11 適應性治療：間歇 vs 持續治療，在有／無承載量競爭與抗藥代價下的惡化時間。
所有參數都是教學假設；09 重用主實驗的 30 次治療試驗，其餘使用獨立 seed 區段。
"""

from __future__ import annotations

import argparse
import json
import os
from concurrent.futures import ProcessPoolExecutor
from dataclasses import asdict, replace
from pathlib import Path

from modeling_cancer import (
    COLORS,
    Config,
    crowding_death,
    growth_probabilities,
    save_figure,
    set_plot_style,
    simulate,
    therapy_outcome,
    write_csv,
)
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.ticker import MaxNLocator
import numpy as np

ROOT = Path(__file__).resolve().parent
INK, GRAY = "#172B42", "#64748B"


# ---------------------------------------------------------------- 共用理論 -------------
def offspring_law(death: float, triple: float) -> np.ndarray:
    """單一父細胞的子代數分布 P(0), P(1), P(2), P(3)。"""
    return np.array([death, 0.0, (1 - death) * (1 - triple), (1 - death) * triple])


def pgf(law: np.ndarray, s: float) -> float:
    return float(law @ s ** np.arange(4))


def extinction(law: np.ndarray) -> float:
    """家系最終滅絕機率：f(s) = s 在 [0,1] 的最小根（由 0 迭代）。"""
    s = 0.0
    for _ in range(200_000):
        nxt = pgf(law, s)
        if abs(nxt - s) < 1e-15:
            return nxt
        s = nxt
    return s


def wilson(successes: int, n: int, z: float = 1.96) -> tuple[float, float]:
    """二項比例的 Wilson 95% 區間。"""
    if n == 0:
        return np.nan, np.nan
    p = successes / n
    centre = (p + z * z / (2 * n)) / (1 + z * z / n)
    half = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / (1 + z * z / n)
    return centre - half, centre + half


def spearman(x, y) -> float:
    rx, ry = np.argsort(np.argsort(x)), np.argsort(np.argsort(y))
    return float(np.corrcoef(rx, ry)[0, 1])


def run_trial(args) -> dict:
    """子程序：跑一次試驗，只回傳摘要與逐代 N／R／治療狀態，避免傳回整個基因型快照。"""
    config, seed = args
    result = simulate(config, seed=seed, therapy=True)
    start = result.therapy_generation
    return {
        "seed": seed,
        "stop": result.stop_reason,
        "therapy_generation": start,
        "total": [r["total"] for r in result.history],
        "resistant": [r["resistant"] for r in result.history],
        "active": result.therapy_active,
        "resistant_at_start": (
            result.history[start]["resistant"] if start is not None else None
        ),
        **therapy_outcome(result),
    }


def run_many(jobs, workers):
    with ProcessPoolExecutor(max_workers=workers) as pool:
        return list(pool.map(run_trial, jobs, chunksize=4))


# ---------------------------------------------------------------- 08 Goldie–Coldman ----
def goldie_coldman_theory(config: Config, threshold: float) -> dict:
    """突變 clone 的 Poisson 近似（Iwasa–Nowak–Michor 2006 型）。

    確定性敏感族群 S_t = N0·m^t；第 s 代出生的抗藥突變數期望 μ·m·S_{s−1}。
    該 clone 先未治療繁殖 T−s 代，再進入治療：最終存活機率 1 − f^{(T−s)}(q_on)，
    f 為未治療子代 PGF、q_on 為治療下抗藥家系滅絕機率。
    治療中敏感細胞仍會產生新抗藥突變，每個存活機率 1 − q_on。治癒 = 沒有任何抗藥家系存活 ≈ exp(−Λ)。
    只用基準細胞（忽略 oncogene／TSG），因此是近似。
    """
    d, q, mu = config.death_rate, config.triple_probability, config.mutation_rate
    off = offspring_law(d, q)
    on_r = offspring_law(1 - (1 - d) * (1 - config.resistant_therapy_kill), q)
    m = (1 - d) * (2 + q)
    m_s_on = m * (1 - config.sensitive_therapy_kill)
    q_on = extinction(on_r)
    T = int(
        np.ceil(
            np.log(max(threshold, config.initial_cells) / config.initial_cells)
            / np.log(m)
            - 1e-12
        )
    )
    S = config.initial_cells * m ** np.arange(T + 1)
    lam_pre, mean_r0, survive = 0.0, 0.0, q_on
    # 從治療開始往回推：第 s 代出生的 clone 還有 T−s 代未治療繁殖。
    for s in range(T, 0, -1):
        births = mu * m * S[s - 1]
        lam_pre += births * (1 - survive)
        mean_r0 += births * m ** (T - s)
        survive = pgf(off, survive)
    lam_on = mu * S[T] * m_s_on / (1 - m_s_on) * (1 - q_on)
    return {
        "threshold": threshold,
        "T": T,
        "N_T": float(S[T]),
        "lambda_pre": lam_pre,
        "lambda_on": lam_on,
        "p_cure": float(np.exp(-(lam_pre + lam_on))),
        "mean_R0": mean_r0,
        "q_on": q_on,
    }


def analysis_08(config, out, trials, seed, workers):
    sizes = [100, 150, 200, 300, 500, 700, 1000, 1500, 2000, 3000, 5000]
    jobs, keys = [], []
    for i, size in enumerate(sizes):
        cfg = replace(
            config, therapy_start_size=size, max_cells=4 * size, generations=300
        )
        for rep in range(trials):
            jobs.append((cfg, seed + 400_000 + 1000 * i + rep))
            keys.append(size)
    results = run_many(jobs, workers)
    rows, summary = [], []
    for size, r in zip(keys, results):
        rows.append({
            "threshold": size,
            "seed": r["seed"],
            "status": r["status"],
            "stop": r["stop"],
            "therapy_generation": r["therapy_generation"],
            "resistant_at_start": r["resistant_at_start"],
            "burden_at_start": r["total"][r["therapy_generation"]],
        })
    for size in sizes:
        sub = [r for r in rows if r["threshold"] == size]
        cures = sum(r["status"] == "extinction" for r in sub)
        lo, hi = wilson(cures, len(sub))
        r0 = np.array([r["resistant_at_start"] for r in sub], float)
        theory = goldie_coldman_theory(config, size)
        summary.append({
            "threshold": size,
            "n": len(sub),
            "cures": cures,
            "p_cure": cures / len(sub),
            "ci_low": lo,
            "ci_high": hi,
            "R0_mean": float(r0.mean()),
            "R0_median": float(np.median(r0)),
            "R0_zero_fraction": float(np.mean(r0 == 0)),
            "theory_p_cure": theory["p_cure"],
            "theory_R0_mean": theory["mean_R0"],
            "theory_T": theory["T"],
        })
    grid = np.unique(np.round(np.geomspace(100, 6000, 400)))
    curve = [goldie_coldman_theory(config, g) for g in grid]

    set_plot_style()
    fig, axes = plt.subplots(1, 2, figsize=(16, 6.3), layout="constrained")
    fig.suptitle(
        "Goldie–Coldman | treat earlier, cure more often",
        fontsize=21,
        fontweight="bold",
    )
    ax = axes[0]
    ax.step(
        grid,
        [c["p_cure"] for c in curve],
        where="post",
        color=INK,
        lw=1.8,
        label="Branching-process theory exp(−Λ)",
    )
    x = [s["threshold"] for s in summary]
    p = np.array([s["p_cure"] for s in summary])
    err = np.array([
        [s["p_cure"] - s["ci_low"] for s in summary],
        [s["ci_high"] - s["p_cure"] for s in summary],
    ])
    ax.errorbar(
        x,
        p,
        yerr=err,
        fmt="o",
        color=COLORS[1],
        ms=7,
        capsize=3,
        label=f"Simulation ({trials} trials, 95% Wilson CI)",
    )
    ax.axhline(0.5, color=GRAY, ls=":", lw=1)
    ax.set(
        xscale="log",
        ylim=(-0.03, 1.03),
        xlabel="Tumor burden when continuous therapy starts (cells)",
        ylabel="P(cure) = P(tumor goes extinct)",
        title="A  Cure probability vs size at treatment",
    )
    ax.legend(frameon=False, loc="lower left")
    ax = axes[1]
    ax.plot(
        x,
        [s["R0_mean"] for s in summary],
        "o",
        color=COLORS[1],
        label="Simulation mean",
    )
    ax.plot(
        x,
        [s["R0_median"] for s in summary],
        "s",
        mfc="white",
        color=COLORS[1],
        label="Simulation median",
    )
    ax.step(
        grid,
        [c["mean_R0"] for c in curve],
        where="post",
        color=INK,
        lw=1.8,
        label=r"Theory $E[R_0]$",
    )
    for s in summary:
        ax.annotate(
            f"{100 * s['R0_zero_fraction']:.0f}%",
            (s["threshold"], max(s["R0_mean"], 0.02)),
            xytext=(0, 9),
            textcoords="offset points",
            ha="center",
            fontsize=8.5,
            color=GRAY,
        )
    ax.set(
        xscale="log",
        yscale="symlog",
        xlabel="Tumor burden when therapy starts (cells)",
        ylabel="Resistant cells at therapy start, R₀ (symlog)",
        title="B  Resistant cells already present",
    )
    ax.set_ylim(bottom=-0.05)
    ax.text(
        0.02,
        0.97,
        "Gray % = trials with R₀ = 0",
        transform=ax.transAxes,
        va="top",
        fontsize=9.5,
        color=GRAY,
    )
    ax.legend(frameon=False, loc="lower right")
    half = next((s["threshold"] for s in curve if s["p_cure"] < 0.5), None)
    q_on = curve[0]["q_on"]
    fig.supxlabel(
        "Default model, continuous therapy from the first generation ≥ threshold. Cure"
        " = extinction; otherwise relapse. Theory: mutant clones arise as a Poisson"
        f" process (μ = {config.mutation_rate:g}); a clone born s generations before"
        f" therapy survives with 1 − f^(s)(q_on),\nq_on = {q_on:.3f} = extinction"
        " probability of one resistant lineage under therapy. Theory uses baseline"
        " cells only (no oncogene speed-up). Theory P(cure) falls below 50% at ≈"
        f" {half:,.0f} cells; the R₀ mean is dominated by rare early 'jackpot'"
        " mutations (Luria–Delbrück).",
        fontsize=10,
    )
    save_figure(fig, out, "08_goldie_coldman")
    write_csv(out / "goldie_coldman_trials.csv", rows)
    write_csv(out / "goldie_coldman_summary.csv", summary)
    return {"summary": summary, "theory_half_cure_size": half, "q_on": q_on}


# ---------------------------------------------------------------- 09 復發時間 ----------
def analysis_09(out):
    from phase_state import fitness, reproduce_trials, trial_feature_rows

    config, results, _ = reproduce_trials(ROOT / "outputs")
    w_off = fitness(0, config, False)
    m_r, m_s = fitness(1 << config.resistance_gene, config, True), fitness(
        0, config, True
    )
    rows = []
    t_grid = np.arange(0.01, 80, 0.01)
    for rep, result in enumerate(results):
        T = result.therapy_generation
        outcome = therapy_outcome(result)
        population = result.populations[T]
        n_t = result.history[T]["total"]
        resistant = {
            g: n for g, n in population.items() if g & (1 << config.resistance_gene)
        }
        sensitive = n_t - sum(resistant.values())
        r0 = sum(resistant.values())
        onc_r = (
            sum(
                n
                for g, n in resistant.items()
                if fitness(g, config, False) > w_off + 1e-9
            )
            / r0
        )
        # 基準預測：抗藥細胞全是單純 RES，每代 ×m_R，S 很快消失。
        t_simple = np.log(n_t / r0) / np.log(m_r)
        # 組成預測：每個抗藥基因型以自己的治療中適應度成長，敏感細胞以 m_S 衰減。
        weights = np.array([fitness(g, config, True) for g in resistant])
        counts = np.array(list(resistant.values()), float)
        burden = (counts[:, None] * weights[:, None] ** t_grid).sum(
            0
        ) + sensitive * m_s**t_grid
        dipped = np.cumsum(burden < n_t) > 0
        t_mix = float(t_grid[np.argmax(dipped & (burden >= n_t))])
        rows.append({
            "replicate": rep,
            "therapy_generation": T,
            "burden_at_start": n_t,
            "resistant_at_start": r0,
            "oncogene_fraction_of_resistant": onc_r,
            "observed_relapse_delay": outcome["relapse_generation"] - T,
            "predicted_delay_baseline": t_simple,
            "predicted_delay_genotype": t_mix,
        })
    obs = np.array([r["observed_relapse_delay"] for r in rows], float)
    pred1 = np.array([r["predicted_delay_baseline"] for r in rows])
    pred2 = np.array([r["predicted_delay_genotype"] for r in rows])
    onc = np.array([r["oncogene_fraction_of_resistant"] for r in rows])
    # 觀察值是整數世代（第一次 N ≥ 起始負荷），預測值向上取整後比較。
    stats = {
        "mae_baseline": float(np.mean(np.abs(np.ceil(pred1) - obs))),
        "mae_genotype": float(np.mean(np.abs(np.ceil(pred2) - obs))),
        "exact_genotype": int(np.sum(np.ceil(pred2) == obs)),
        "spearman_obs_vs_onc": spearman(onc, obs),
        "spearman_obs_vs_R0": spearman([r["resistant_at_start"] for r in rows], obs),
        "rho_star": (1 - m_s) / (m_r - m_s),
    }
    features = trial_feature_rows(config, results)
    on = [f for f in features if f["therapy_on"] > 0]
    rho = np.array([f["frac_R"] for f in on])
    growth = np.array([f["growth_expected"] for f in on])
    onc_all = np.array([f["frac_onc"] for f in on])
    stats["observed_sign_change_rho"] = [
        float(rho[growth < 0].max()),
        float(rho[growth > 0].min()),
    ]

    set_plot_style()
    fig, axes = plt.subplots(1, 3, figsize=(19, 6.3), layout="constrained")
    fig.suptitle(
        "Relapse timing is predictable from the resistant cells present at therapy"
        " start",
        fontsize=20,
        fontweight="bold",
    )
    ax = axes[0]
    lim = (obs.min() - 2, obs.max() + 2)
    jitter = np.random.default_rng(9).uniform(-0.12, 0.12, len(obs))
    ax.plot(lim, lim, color=GRAY, lw=1, ls="--")
    ax.scatter(
        pred1,
        obs + jitter,
        color=GRAY,
        s=30,
        alpha=0.7,
        label=f"Baseline: ln(N/R₀)/ln {m_r:.3f}  (MAE {stats['mae_baseline']:.2f})",
    )
    points = ax.scatter(
        pred2,
        obs + jitter,
        c=onc,
        cmap="viridis",
        vmin=0,
        vmax=max(0.5, onc.max()),
        s=46,
        edgecolor="white",
        lw=0.6,
        label=f"Genotype-aware mixture  (MAE {stats['mae_genotype']:.2f})",
    )
    ax.set(
        xlim=lim,
        ylim=lim,
        xlabel="Predicted relapse delay (generations after therapy start)",
        ylabel="Observed relapse delay (generations)",
        title="A  Prediction vs observation (30 trials)",
    )
    ax.legend(frameon=False, loc="upper left", fontsize=9.5)
    bar = fig.colorbar(points, ax=ax, pad=0.01)
    bar.set_label("Oncogene fraction among resistant cells", fontsize=10)
    ax = axes[1]
    r0 = np.array([r["resistant_at_start"] for r in rows])
    points = ax.scatter(
        r0,
        obs + jitter,
        c=onc,
        cmap="viridis",
        vmin=0,
        vmax=max(0.5, onc.max()),
        s=46,
        edgecolor="white",
        lw=0.6,
    )
    grid = np.geomspace(r0.min() * 0.8, r0.max() * 1.2, 100)
    ax.plot(
        grid,
        np.log(np.median([r["burden_at_start"] for r in rows]) / grid) / np.log(m_r),
        color=GRAY,
        ls="--",
        label="Baseline prediction at median N",
    )
    ax.set(
        xscale="log",
        xlabel="Resistant cells at therapy start, R₀",
        ylabel="Observed relapse delay (generations)",
        title="B  More pre-existing resistance → earlier relapse",
    )
    ax.text(
        0.98,
        0.97,
        f"Spearman: R₀ {stats['spearman_obs_vs_R0']:+.2f}\noncogene share"
        f" {stats['spearman_obs_vs_onc']:+.2f}",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=10,
    )
    ax.legend(frameon=False, loc="lower left", fontsize=9.5)
    ax = axes[2]
    points = ax.scatter(
        rho, growth, c=onc_all, cmap="viridis", s=12, vmin=0, vmax=1, lw=0, alpha=0.85
    )
    g = np.linspace(0, 1, 200)
    ax.plot(
        g,
        np.log((1 - g) * m_s + g * m_r),
        color=INK,
        lw=1.6,
        label=r"Theory $\ln[(1-\rho)m_S+\rho\,m_R]$",
    )
    ax.axhline(0, color=GRAY, lw=1)
    ax.axvline(
        stats["rho_star"],
        color=COLORS[1],
        ls="--",
        lw=1.5,
        label=f"ρ* = {stats['rho_star']:.0%}: decline turns to regrowth",
    )
    ax.set(
        xlabel="Resistant fraction ρ",
        ylabel="Expected ln growth per generation (on therapy)",
        title="C  The nadir happens when ρ crosses ρ*",
    )
    ax.legend(frameon=False, loc="upper left", fontsize=9.5)
    bar = fig.colorbar(points, ax=ax, pad=0.01)
    bar.set_label("Oncogene fraction (all cells)", fontsize=10)
    fig.supxlabel(
        "Main experiment: 30 trials, continuous therapy from 20,000 cells. Relapse"
        " delay = first generation after therapy start with burden ≥ burden at start."
        " Baseline prediction assumes every resistant cell grows ×m_R per"
        " generation;\nthe mixture prediction lets each resistant genotype grow at its"
        " own on-therapy mean (oncogene cells faster) and lets sensitive cells decay"
        " ×m_S. Predictions are rounded up to whole generations for the MAE. Vertical"
        " jitter ±0.12 only separates overlapping points.",
        fontsize=10,
    )
    save_figure(fig, out, "09_relapse_time_prediction")
    write_csv(out / "relapse_time_prediction.csv", rows)
    return stats


# ---------------------------------------------------------------- 10 參數敏感度 --------
def analysis_10(config, out, trials, seed, workers):
    kills = [0.0, 0.1, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6]
    mus = [1e-4, 3e-4, 1e-3, 3e-3]
    jobs, keys = [], []
    for i, k in enumerate(kills):
        for j, mu in enumerate(mus):
            cfg = replace(
                config,
                resistant_therapy_kill=k,
                mutation_rate=mu,
                max_cells=2 * config.therapy_start_size,
                generations=250,
            )
            for rep in range(trials):
                jobs.append((cfg, seed + 500_000 + 10_000 * i + 1000 * j + rep))
                keys.append((k, mu))
    results = run_many(jobs, workers)
    rows = []
    for (k, mu), r in zip(keys, results):
        # fmt: off
        status = (
            "cure" if r["stop"] == "extinction"
            else "relapse" if r["stop"] == "cell_limit"
            else "persistent"
        )
        # fmt: on
        rows.append({
            "resistant_kill": k,
            "mutation_rate": mu,
            "seed": r["seed"],
            "status": status,
            "therapy_generation": r["therapy_generation"],
            "resistant_at_start": r["resistant_at_start"],
            "relapse_delay": (
                (r["relapse_generation"] - r["therapy_generation"])
                if r["relapse_generation"] is not None
                else None
            ),
        })
    cure = np.zeros((len(mus), len(kills)))
    delay = np.full((len(mus), len(kills)), np.nan)
    persistent = np.zeros_like(cure)
    for i, k in enumerate(kills):
        for j, mu in enumerate(mus):
            sub = [
                r for r in rows if r["resistant_kill"] == k and r["mutation_rate"] == mu
            ]
            cure[j, i] = np.mean([r["status"] == "cure" for r in sub])
            persistent[j, i] = np.mean([r["status"] == "persistent" for r in sub])
            d = [r["relapse_delay"] for r in sub if r["relapse_delay"] is not None]
            if d:
                delay[j, i] = np.median(d)
    # 抗藥族群在治療下停止成長的殺傷率：(1−k)·m = 1。
    boundaries = {}
    for name, genotype in (
        ("baseline resistant", 1 << config.resistance_gene),
        ("+ oncogene", (1 << config.resistance_gene) | 1),
        (
            "+ oncogene + TSG loss",
            (1 << config.resistance_gene) | 1 | (1 << 4) | (1 << 20),
        ),
    ):
        dth, q = growth_probabilities(genotype, config)
        boundaries[name] = 1 - 1 / ((1 - dth) * (2 + q))

    set_plot_style()
    fig, axes = plt.subplots(1, 2, figsize=(18, 6.6), layout="constrained")
    fig.suptitle(
        "Parameter sensitivity | cure needs therapy that also out-kills resistant"
        " growth",
        fontsize=20,
        fontweight="bold",
    )
    cmaps = (
        LinearSegmentedColormap.from_list("cure", ["#F7F4EE", "#167D9A"]),
        LinearSegmentedColormap.from_list("delay", ["#FBE9DD", "#B4472A"]),
    )
    for ax, data, cmap, title, fmt, label in (
        (axes[0], cure, cmaps[0], "A  P(cure)", "{:.0%}", "P(cure)"),
        (
            axes[1],
            delay,
            cmaps[1],
            "B  Median relapse delay (relapsed trials)",
            "{:.0f}",
            "Generations after therapy start",
        ),
    ):
        image = ax.imshow(
            np.ma.masked_invalid(data),
            cmap=cmap,
            aspect="auto",
            origin="lower",
            vmin=0,
            vmax=1 if ax is axes[0] else np.nanmax(delay),
        )
        for j in range(len(mus)):
            for i in range(len(kills)):
                value = data[j, i]
                text = "—" if not np.isfinite(value) else fmt.format(value)
                if ax is axes[0] and persistent[j, i] > 0:
                    text += f"\n({persistent[j, i]:.0%} pers.)"
                ax.text(
                    i,
                    j,
                    text,
                    ha="center",
                    va="center",
                    fontsize=8.5,
                    color=(
                        "white"
                        if np.isfinite(value)
                        and value > (0.6 if ax is axes[0] else 0.6 * np.nanmax(delay))
                        else INK
                    ),
                )
        ax.set_xticks(range(len(kills)), [f"{k:g}" for k in kills])
        ax.set_yticks(range(len(mus)), [f"{m:g}" for m in mus])
        ax.set(
            xlabel="Therapy kill of resistant cells, k_R (default 0.04)",
            ylabel="Mutation rate μ per locus (default 0.001)",
            title=title,
        )
        ax.grid(False)
        bar = fig.colorbar(image, ax=ax, pad=0.01)
        bar.set_label(label)
        # 邊界畫在 k 軸上的內插位置。
        for (name, b), style in zip(boundaries.items(), ("-", "--", ":")):
            pos = np.interp(b, kills, range(len(kills)))
            ax.vlines(
                pos,
                -0.5,
                len(mus) - 0.5,
                color=INK,
                ls=style,
                lw=1.6,
                label=f"{name} stops growing (k = {b:.2f})",
            )
        # 上方留一條空白帶放圖例，避免遮住格子。
        ax.set_ylim(-0.5, len(mus) - 0.5 + 0.75)
        if ax is axes[0]:
            ax.legend(frameon=False, fontsize=8.5, loc="upper left", ncol=2)
    fig.supxlabel(
        f"{trials} trials per cell; continuous therapy from 20,000 cells; sensitive"
        " kill 0.90. Cure = extinction; relapse = burden reaches 40,000; persistent ="
        " still alive after 250 generations.\nVertical lines: k_R where (1 − k_R)·m ="
        " 1 for each resistant genotype. Above a line that genotype shrinks on"
        " therapy; cure still fails if a fitter resistant genotype exists.",
        fontsize=10,
    )
    save_figure(fig, out, "10_parameter_sensitivity")
    write_csv(out / "parameter_sensitivity_trials.csv", rows)
    return {
        "kills": kills,
        "mutation_rates": mus,
        "p_cure": cure.tolist(),
        "median_delay": np.where(np.isfinite(delay), delay, None).tolist(),
        "persistent": persistent.tolist(),
        "boundaries": boundaries,
    }


# ---------------------------------------------------------------- 11 適應性治療 --------
def time_to_progression(record, start_size, factor=1.2):
    """惡化：治療開始後負荷曾降到門檻以下，之後在「上一代有用藥」的情況下仍 ≥ factor × 門檻。

    要求用藥中，是因為離散世代下暫停期的一代就可能長 20–30%，恢復用藥前的那一步不算治療失敗。
    回傳治療開始後的代數；觀察期內未惡化回傳 None（右設限）。
    """
    T = record["therapy_generation"]
    if T is None:
        return None
    dipped = False
    for g in range(T + 1, len(record["total"])):
        dipped |= record["total"][g] < start_size
        if (
            dipped
            and record["active"][g - 1]
            and record["total"][g] >= factor * start_size
        ):
            return g - T
    return None


def analysis_11(config, out, trials, seed, workers):
    K = int(1.5 * config.therapy_start_size)
    settings = [
        ("No competition", None, 0.0),
        ("No competition + cost", None, 0.05),
        (f"Competition (K = {K:,})", K, 0.0),
        ("Competition + cost", K, 0.05),
    ]
    jobs, keys = [], []
    for i, (name, capacity, cost) in enumerate(settings):
        for mode in ("continuous", "adaptive"):
            cfg = replace(
                config,
                carrying_capacity=capacity,
                resistance_cost=cost,
                therapy_mode=mode,
                generations=200,
                max_cells=int(1.8 * config.therapy_start_size),
            )
            for rep in range(trials):
                # 兩種策略用同一組 seed：治療開始前的軌跡完全相同。
                jobs.append((cfg, seed + 600_000 + 1000 * i + rep))
                keys.append((name, mode))
    results = run_many(jobs, workers)
    rows, examples = [], {}
    for (name, mode), r in zip(keys, results):
        ttp = time_to_progression(r, config.therapy_start_size)
        rows.append({
            "setting": name,
            "strategy": mode,
            "seed": r["seed"],
            "stop": r["stop"],
            "therapy_generation": r["therapy_generation"],
            "time_to_progression": ttp,
            "therapy_fraction": (
                float(np.mean(r["active"][r["therapy_generation"] :]))
                if r["therapy_generation"] is not None
                else None
            ),
        })
        if name == settings[3][0] and (name, mode) not in examples:
            examples[(name, mode)] = r
    summary = []
    for name, *_ in settings:
        entry = {"setting": name}
        for mode in ("continuous", "adaptive"):
            sub = [r for r in rows if r["setting"] == name and r["strategy"] == mode]
            t = [
                r["time_to_progression"]
                for r in sub
                if r["time_to_progression"] is not None
            ]
            entry[f"{mode}_median_ttp"] = float(np.median(t)) if t else None
            entry[f"{mode}_progressed"] = len(t)
            entry[f"{mode}_n"] = len(sub)
            entry[f"{mode}_extinct"] = sum(r["stop"] == "extinction" for r in sub)
            entry[f"{mode}_median_therapy_fraction"] = float(
                np.median([r["therapy_fraction"] for r in sub])
            )
        # 配對比較（同 seed）：adaptive 惡化時間較長的試驗比例。
        pairs = [
            (a["time_to_progression"], b["time_to_progression"])
            for a, b in zip(
                [
                    r
                    for r in rows
                    if r["setting"] == name and r["strategy"] == "continuous"
                ],
                [
                    r
                    for r in rows
                    if r["setting"] == name and r["strategy"] == "adaptive"
                ],
            )
        ]
        entry["adaptive_longer_fraction"] = float(
            np.mean([
                (b if b is not None else 1e9) > (a if a is not None else 1e9)
                for a, b in pairs
            ])
        )
        summary.append(entry)

    set_plot_style()
    fig = plt.figure(figsize=(19, 7.2), layout="constrained")
    grid = fig.add_gridspec(2, 3, width_ratios=[1.15, 1.15, 1.3])
    fig.suptitle(
        "Adaptive therapy | it only helps when sensitive cells hold resistant cells"
        " back",
        fontsize=20,
        fontweight="bold",
    )
    for row, mode in enumerate(("continuous", "adaptive")):
        ax = fig.add_subplot(grid[row, :2])
        r = examples[(settings[3][0], mode)]
        t = np.arange(len(r["total"]))
        ax.plot(t, r["total"], color=INK, label="Total")
        ax.plot(t, r["resistant"], color=COLORS[1], label="Resistant")
        ax.plot(
            t,
            np.array(r["total"]) - np.array(r["resistant"]),
            color=COLORS[0],
            label="Sensitive",
        )
        for g, on in enumerate(r["active"]):
            if on:
                ax.axvspan(g, g + 1, color=COLORS[1], alpha=0.10, lw=0)
        ax.axhline(1.2 * config.therapy_start_size, color=GRAY, ls=":", lw=1)
        ax.axhline(K, color=GRAY, ls="--", lw=1)
        ttp = time_to_progression(r, config.therapy_start_size)
        title = (
            f"{'A' if row == 0 else 'B'}  {mode.title()} therapy, competition + cost"
            " (one trial; shaded = drug on)"
            + (f" — progression after {ttp} gen" if ttp else "")
        )
        end = min(
            len(t) - 1, (r["therapy_generation"] + ttp + 6) if ttp else len(t) - 1
        )
        ax.set(ylabel="Cells", title=title, xlim=(0, end))
        ax.set_yscale("symlog", linthresh=10)
        ax.set_ylim(bottom=0)
        ax.title.set_fontsize(12.5)
        if row == 0:
            ax.legend(frameon=False, ncol=3, loc="lower left", fontsize=9)
            ax.text(end, K, " K", va="center", fontsize=9, color=GRAY)
        else:
            ax.set_xlabel("Generation")
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax = fig.add_subplot(grid[:, 2])
    for i, (name, *_) in enumerate(settings):
        for k, (mode, color) in enumerate(
            (("continuous", COLORS[0]), ("adaptive", COLORS[1]))
        ):
            values = [
                r["time_to_progression"]
                for r in rows
                if r["setting"] == name
                and r["strategy"] == mode
                and r["time_to_progression"] is not None
            ]
            y = i + (k - 0.5) * 0.34
            jitter = np.random.default_rng(i * 10 + k).uniform(-0.08, 0.08, len(values))
            ax.scatter(
                values,
                y + jitter,
                color=color,
                s=18,
                alpha=0.6,
                label=mode.title() if i == 0 else None,
            )
            if values:
                ax.plot([np.median(values)] * 2, [y - 0.14, y + 0.14], color=INK, lw=3)
        entry = summary[i]
        ratio = entry["adaptive_median_ttp"] / entry["continuous_median_ttp"]
        ax.text(
            1.01,
            i,
            f"×{ratio:.2f}",
            transform=ax.get_yaxis_transform(),
            va="center",
            fontsize=11,
            weight="bold",
            color=COLORS[1] if ratio > 1.05 else GRAY,
        )
    ax.set_yticks(range(len(settings)), [s[0] for s in settings])
    ax.invert_yaxis()
    ax.set(
        xlabel="Time to progression (generations after therapy start)",
        title="C  Adaptive vs continuous, 4 settings",
    )
    ax.legend(frameon=False, loc="upper right")
    fig.supxlabel(
        f"{trials} paired trials per setting (same seeds for both strategies)."
        " Adaptive: pause when burden ≤ 50% of the burden at first treatment; resume"
        " at 20,000 cells. Progression = burden ≥ 24,000 while on drug, after first"
        " dropping below 20,000.\nCompetition: every cell gets extra death 0.21·N/K,"
        f" so a mutation-free cell has mean offspring 1 at N = K = {K:,}. Cost:"
        " resistant cells +0.05 death. Right labels: adaptive / continuous median"
        " ratio.",
        fontsize=10,
    )
    save_figure(fig, out, "11_adaptive_therapy")
    write_csv(out / "adaptive_therapy_trials.csv", rows)
    write_csv(out / "adaptive_therapy_summary.csv", summary)
    return {
        "summary": summary,
        "carrying_capacity": K,
        "crowding_at_K": crowding_death(K, replace(config, carrying_capacity=K)),
    }


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--only", choices=("08", "09", "10", "11"), nargs="*")
    parser.add_argument(
        "--trials", type=int, default=200, help="trials per threshold in 08"
    )
    parser.add_argument(
        "--grid-trials",
        type=int,
        default=30,
        help="trials per cell in 10 and per arm in 11",
    )
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--workers", type=int, default=os.cpu_count())
    parser.add_argument(
        "--output", type=Path, default=ROOT / "outputs/therapy_analysis"
    )
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    config = Config()
    summary_path = args.output / "therapy_analysis_summary.json"
    summary = json.loads(summary_path.read_text()) if summary_path.exists() else {}
    todo = args.only or ["08", "09", "10", "11"]
    if "08" in todo:
        summary["08_goldie_coldman"] = analysis_08(
            config, args.output, args.trials, args.seed, args.workers
        )
    if "09" in todo:
        summary["09_relapse_time"] = analysis_09(args.output)
    if "10" in todo:
        summary["10_parameter_sensitivity"] = analysis_10(
            config, args.output, args.grid_trials, args.seed, args.workers
        )
    if "11" in todo:
        summary["11_adaptive_therapy"] = analysis_11(
            config, args.output, args.grid_trials, args.seed, args.workers
        )
    summary["config"] = asdict(config)
    summary_path.write_text(json.dumps(summary, indent=2, default=float) + "\n")
    print(
        json.dumps(
            {k: v for k, v in summary.items() if k != "config"}, indent=1, default=float
        )[:6000]
    )


if __name__ == "__main__":
    main()
