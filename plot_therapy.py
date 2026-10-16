"""治療延伸分析圖：治療時機（Goldie–Coldman）、復發時間預測、參數敏感度、適應性治療。

資料由 experiments.py 計算；這裡只負責畫圖與存檔。
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.ticker import MaxNLocator
import numpy as np

from experiments import time_to_progression
from plot_style import COLORS, save_figure, set_plot_style
from simulation import Config

INK, GRAY = "#172B42", "#64748B"


# Goldie–Coldman：治療起始負荷 vs 治癒機率（模擬與理論），以及治療開始時已有的抗藥細胞。
def plot_goldie_coldman(data: dict, config: Config, trials: int, output: Path):
    """Goldie–Coldman：治療起始負荷 vs 治癒機率（模擬與理論），以及治療開始時已有的抗藥細胞。"""
    summary, grid, curve = data["summary"], data["grid"], data["curve"]
    half = data["half"]
    m = (1 - config.death_rate) * (2 + config.triple_probability)
    set_plot_style()
    fig, axes = plt.subplots(1, 2, figsize=(16, 6.8), layout="constrained")
    fig.suptitle(
        "Treat earlier, cure more often: the chance of cure falls with tumor size",
        fontsize=19,
        fontweight="bold",
    )
    x = np.array([s["threshold"] for s in summary])
    # A：治癒機率。階梯 = 理論；腫瘤每代 ×m 跳著長，相近門檻在同一代開始治療。
    ax = axes[0]
    ax.step(
        grid,
        [c["p_cure"] for c in curve],
        where="post",
        color=INK,
        lw=1.8,
        label="Theory: P(cure) = e^(−Λ), no fitted parameters",
    )
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
        label=f"Simulation: share of {trials} tumors cured (95% CI)",
    )
    ax.axhline(0.5, color=GRAY, ls=":", lw=1)
    if half:
        ax.axvline(half, color=GRAY, ls=":", lw=1)
        ax.annotate(
            f"50% cure at ≈ {half:,.0f} cells",
            (half, 0.5),
            xytext=(12, 14),
            textcoords="offset points",
            fontsize=11,
            color=INK,
        )
    for size, dx, dy in ((100, 14, -4), (1000, 12, 18)):
        k = int(np.where(x == size)[0][0])
        ax.annotate(
            f"start at {size:,} cells → {p[k]:.1%} cured",
            (size, p[k]),
            xytext=(dx, dy),
            textcoords="offset points",
            fontsize=10.5,
            color=COLORS[1],
            arrowprops=dict(arrowstyle="-", color=COLORS[1], lw=0.8),
        )
    ax.text(
        0.98,
        0.97,
        f"Steps: the tumor grows ×{m:.2f} per generation,\nso nearby sizes start"
        " therapy in the same generation.",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=9.5,
        color=GRAY,
    )
    ax.set(
        xscale="log",
        ylim=(-0.03, 1.03),
        xlabel="Tumor size when continuous therapy starts (cells, log scale)",
        ylabel="Probability of cure (tumor goes extinct)",
        title="A  Bigger tumor at treatment → lower chance of cure",
    )
    ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{v:.0%}"))
    ax.legend(frameon=False, loc="lower left")
    # B：原因。治療開始時已存在的抗藥細胞數；平均被少數早期 jackpot 拉高。
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
        label=r"Theory: expected $R_0$",
    )
    for s in summary:
        ax.annotate(
            f"{100 * s['R0_zero_fraction']:.0f}%",
            (s["threshold"], max(s["R0_mean"], 0.02)),
            xytext=(0, 9),
            textcoords="offset points",
            ha="center",
            fontsize=9,
            color=GRAY,
        )
    ax.set(
        xscale="log",
        yscale="symlog",
        xlabel="Tumor size when continuous therapy starts (cells, log scale)",
        ylabel="Resistant cells already present, R₀",
        title="B  Why: bigger tumors already carry resistant cells",
    )
    ax.set_ylim(-0.05, 2.5 * max(c["mean_R0"] for c in curve))
    ax.text(
        0.02,
        0.97,
        "Gray % above each point = trials with no resistant cell (R₀ = 0)\n"
        "Mean > median: a few early mutations grow into large 'jackpot' clones",
        transform=ax.transAxes,
        va="top",
        fontsize=9.5,
        color=GRAY,
    )
    ax.legend(frameon=False, loc="lower right")
    save_figure(fig, output, "goldie_coldman")


# 復發時間預測：基準預測與組成預測各自對照實際復發延遲。
def plot_relapse_prediction(data: dict, output: Path):
    """復發時間預測：基準預測與組成預測各自對照實際復發延遲。"""
    stats = data["stats"]
    obs, pred1, pred2 = data["obs"], data["pred_baseline"], data["pred_genotype"]
    onc, m_r = data["oncogene_share"], data["m_r"]
    set_plot_style()
    fig, ax = plt.subplots(figsize=(9.5, 7.2), layout="constrained")
    fig.suptitle(
        "Relapse timing predicted from the composition at therapy start",
        fontsize=18,
        fontweight="bold",
    )
    lim = (obs.min() - 2, obs.max() + 2)
    # 縱向微幅抖動只為分開重疊的整數世代點。
    jitter = np.random.default_rng(9).uniform(-0.12, 0.12, len(obs))
    ax.plot(lim, lim, color=GRAY, lw=1, ls="--", label="Perfect prediction")
    ax.scatter(
        pred1,
        obs + jitter,
        color=GRAY,
        s=34,
        alpha=0.7,
        label=f"Baseline: ln(N_T/R₀) / ln {m_r:.3f}  (MAE {stats['mae_baseline']:.2f})",
    )
    points = ax.scatter(
        pred2,
        obs + jitter,
        c=onc,
        cmap="viridis",
        vmin=0,
        vmax=max(0.5, onc.max()),
        s=52,
        edgecolor="white",
        lw=0.6,
        label=(
            f"Composition: each genotype at its own rate  (MAE"
            f" {stats['mae_genotype']:.2f}, {stats['exact_genotype']}/{len(obs)} exact)"
        ),
    )
    ax.set(
        xlim=lim,
        ylim=lim,
        xlabel="Predicted relapse delay (generations after therapy start)",
        ylabel="Observed relapse delay (generations)",
        title=f"Prediction vs observation ({len(obs)} trials)",
    )
    ax.set_aspect("equal", adjustable="box")
    ax.legend(frameon=False, loc="upper left", fontsize=10)
    bar = fig.colorbar(points, ax=ax, pad=0.01, shrink=0.85)
    bar.set_label("Oncogene fraction among resistant cells", fontsize=10)
    save_figure(fig, output, "relapse_time_prediction")


# 參數敏感度：抗藥細胞殺傷率 × 突變率的治癒機率與復發延遲熱圖，並畫出理論邊界。
def plot_sensitivity(data: dict, trials: int, output: Path):
    """參數敏感度：抗藥細胞殺傷率 × 突變率的治癒機率與復發延遲熱圖，並畫出理論邊界。"""
    kills, mus = data["kills"], data["mutation_rates"]
    cure, delay, persistent = data["cure"], data["delay"], data["persistent"]
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
    for ax, values, cmap, title, fmt, label in (
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
            np.ma.masked_invalid(values),
            cmap=cmap,
            aspect="auto",
            origin="lower",
            vmin=0,
            vmax=1 if ax is axes[0] else np.nanmax(delay),
        )
        for j in range(len(mus)):
            for i in range(len(kills)):
                value = values[j, i]
                text = "—" if not np.isfinite(value) else fmt.format(value)
                if ax is axes[0] and persistent[j, i] > 0:
                    text += f"\n({persistent[j, i]:.0%} pers.)"
                dark = np.isfinite(value) and value > (
                    0.6 if ax is axes[0] else 0.6 * np.nanmax(delay)
                )
                ax.text(
                    i,
                    j,
                    text,
                    ha="center",
                    va="center",
                    fontsize=8.5,
                    color="white" if dark else INK,
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
        for (name, b), style in zip(data["boundaries"].items(), ("-", "--", ":")):
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
    save_figure(fig, output, "parameter_sensitivity")


# 適應性治療：一次試驗的持續與間歇治療軌跡，以及四種設定下的惡化時間比較。
def plot_adaptive(data: dict, config: Config, trials: int, output: Path):
    """適應性治療：一次試驗的持續與間歇治療軌跡，以及四種設定下的惡化時間比較。"""
    rows, summary, settings, K = (
        data["rows"],
        data["summary"],
        data["settings"],
        data["K"],
    )
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
        r = data["examples"][(settings[3][0], mode)]
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
    save_figure(fig, output, "adaptive_therapy")
