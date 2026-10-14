"""相圖 01：實際軌跡疊上由完整基因型狀態算出的條件期望箭頭（單一面板）。

執行：.venv/bin/python plot_phase_portraits.py
讀取既有 CSV；只重現一條標示路徑以計算期望，不修改模擬或簡報。
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

# 先載入模型，沿用專案內的 Matplotlib cache 與 Agg backend。
from modeling_cancer import (Config, Result, growth_probabilities,
                             il11_environment, simulate, therapy_outcome)
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch
import numpy as np


ROOT = Path(__file__).resolve().parent
INK, BLUE, ORANGE, RED, GRAY = "#173042", "#236B8E", "#BC7844", "#A74F3D", "#687783"


def set_phase_style():
    """中文標籤及 Computer Modern 數學字型；SVG 保留字形路徑。"""
    available = {f.name for f in font_manager.fontManager.ttflist}
    candidates = ("Arial Unicode MS", "Microsoft JhengHei", "Noto Sans CJK TC", "PingFang TC")
    family = next((f for f in candidates if f in available), None)
    if family is None:
        raise RuntimeError("需要支援繁體中文的字型：Arial Unicode MS 或 Noto Sans CJK TC。")
    plt.rcParams.update({"font.family": family, "font.size": 13,
                         "mathtext.fontset": "cm", "axes.unicode_minus": False,
                         "axes.spines.top": False, "axes.spines.right": False,
                         "axes.labelcolor": INK, "text.color": INK,
                         "axes.edgecolor": GRAY, "xtick.color": GRAY, "ytick.color": GRAY,
                         "figure.facecolor": "white", "axes.facecolor": "white",
                         "svg.fonttype": "path"})


def read_csv(path):
    with Path(path).open(newline="") as handle:
        return list(csv.DictReader(handle))


def load_trials(directory):
    """只使用實際觀察；不在 cell_limit 後延伸路徑或連接不同試驗。"""
    directory = Path(directory)
    metadata = json.loads((directory / "run_metadata.json").read_text())
    settings = metadata["config"].copy()
    for key in ("driver_genes", "suppressor_genes"):
        settings[key] = tuple(settings[key])
    config = Config(**settings)
    config.validate()
    events = {}
    for row in read_csv(directory / "therapy_replicates.csv"):
        events[int(row["replicate"])] = {
            **row, "replicate": int(row["replicate"]), "seed": int(row["seed"]),
            **{key: int(row[key]) if row[key] else None
               for key in ("therapy_generation", "relapse_generation", "nadir")}}
    trials = {rep: [] for rep in events}
    for row in read_csv(directory / "therapy_trajectories.csv"):
        if row["extended_extinction"] == "True":
            continue
        rep = int(row["replicate"])
        state = {key: int(row[key]) for key in ("generation", "total", "sensitive", "resistant")}
        assert state["sensitive"] + state["resistant"] == state["total"]
        assert min(state["sensitive"], state["resistant"]) >= 0
        assert int(row["seed"]) == events[rep]["seed"]
        trials[rep].append(state)
    for states in trials.values():
        states.sort(key=lambda row: row["generation"])
        assert [row["generation"] for row in states] == list(range(len(states)))
    assert len(trials) == metadata["replicates"]
    return config, trials, events


def choose_exemplar(events):
    """先選最接近復發中位代數，再選谷底最接近中位數的實際試驗。"""
    eligible = [e for e in events.values() if e["relapse_generation"] is not None and e["nadir"]]
    if not eligible:
        return min(events)
    relapse = np.median([e["relapse_generation"] for e in eligible])
    nadir = np.median([e["nadir"] for e in eligible])
    return min(eligible, key=lambda e: (abs(e["relapse_generation"] - relapse),
                                       abs(np.log(e["nadir"] / nadir)), e["replicate"]))["replicate"]


def projected_offspring_mean(population: dict[int, int], config: Config,
                             therapy_active: bool) -> dict[str, float]:
    """精確的條件期望；其他基因突變不改變子代的 S/R 分組機率。

    S 子代以 μ 成為 R，既有 R 不逆轉。B_S/B_R 是突變前的期望子代數。
    這不是只依賴 (S,R) 的封閉模型：m_g 需要完整基因型與共享環境。
    """
    _, benefit = il11_environment(population, config)
    before = {"sensitive": 0.0, "resistant": 0.0}
    for genotype, count in population.items():
        death, triple = growth_probabilities(genotype, config, benefit)
        resistant = bool(genotype & (1 << config.resistance_gene))
        kill = (config.resistant_therapy_kill if resistant else config.sensitive_therapy_kill)
        mean = (1 - death) * (1 - kill if therapy_active else 1) * (2 + triple)
        before["resistant" if resistant else "sensitive"] += count * mean
    sensitive = (1 - config.mutation_rate) * before["sensitive"]
    resistant = before["resistant"] + config.mutation_rate * before["sensitive"]
    assert np.isclose(sensitive + resistant, sum(before.values()), rtol=1e-12)
    return {"expected_sensitive": sensitive, "expected_resistant": resistant,
            "expected_total": sensitive + resistant,
            "before_mutation_sensitive": before["sensitive"],
            "before_mutation_resistant": before["resistant"]}


def phase_coordinates(states, composition=False):
    total = np.array([r["total"] for r in states], dtype=float)
    if composition:
        # 滅絕時比例未定義，用 NaN 讓圖斷開，而非畫成 0% 抗藥。
        fraction = np.divide([r["resistant"] for r in states], total,
                             out=np.full(len(states), np.nan), where=total > 0)
        return np.column_stack((np.log10(total + 1), 100 * fraction))
    return np.log10(np.array([[r["sensitive"], r["resistant"]] for r in states]) + 1)


def decorate_axis(ax, title, composition=False):
    ax.set_title(title, loc="left", fontsize=17, weight="bold", pad=16)
    ax.grid(color="#DDE4E8", alpha=.6, linewidth=.7)
    ax.set_axisbelow(True)
    if composition:
        ax.set_xlabel(r"總負荷 $\log_{10}(T+1)$", labelpad=12)
        ax.set_ylabel(r"抗藥比例 $100R/T$（%）", labelpad=12)
        ax.set_ylim(-4, 107)
    else:
        ax.set_xlabel(r"敏感細胞 $\log_{10}(S+1)$", labelpad=12)
        ax.set_ylabel(r"抗藥細胞 $\log_{10}(R+1)$", labelpad=12)
        ax.set_xlim(-.16, 4.85)
        ax.set_ylim(-.18, 6.55)


def arrow(ax, a, b, color, size=12, linewidth=2, zorder=6):
    if np.isfinite([*a, *b]).all():
        ax.add_patch(FancyArrowPatch(a, b, arrowstyle="-|>", mutation_scale=size,
                                    color=color, linewidth=linewidth, shrinkA=0, shrinkB=0,
                                    zorder=zorder))


def event_generations(states, event):
    start, relapse = event["therapy_generation"], event["relapse_generation"]
    points = {"初始": 0}
    if start is not None:
        points["治療啟動"] = start
        after = states[start + 1:]
        if after:
            points["負荷谷底"] = min(after, key=lambda r: r["total"])["generation"]
    if relapse is not None:
        points["首次復發"] = relapse
    return points


def conditional_steps(result: Result):
    """每個可推進世代的原始計數期望，以及實際下一代計數。"""
    rows = []
    for g, population in enumerate(result.populations[:-1]):
        active = result.therapy_generation is not None and g >= result.therapy_generation
        current, following = result.history[g:g+2]
        mean = projected_offspring_mean(population, result.config, active)
        rows.append({"generation": g, "therapy_active": active,
                     "sensitive": current["sensitive"], "resistant": current["resistant"], **mean,
                     "expected_delta_sensitive": mean["expected_sensitive"] - current["sensitive"],
                     "expected_delta_resistant": mean["expected_resistant"] - current["resistant"],
                     "observed_next_sensitive": following["sensitive"],
                     "observed_next_resistant": following["resistant"]})
    return rows


def plot_phase_portrait(trials, events, replicate, config, result, steps):
    """單一相平面：30 條實際軌跡、突出案例，以及由完整基因型狀態算出的期望下一代箭頭。

    取代舊的 01（軌跡兩格）與 02（條件期望）：同一試驗、同一座標不再畫兩次。
    負荷 × 抗藥比例的座標改看相圖 05（logit 軸）；推導公式在圖 07。
    """
    fig = plt.figure(figsize=(15.5, 10.2))
    ax = fig.add_axes((.07, .2, .5, .58))
    side = fig.add_axes((.62, .2, .35, .58)); side.axis("off")
    fig.text(.055, .965, "相平面：負荷下降不等於抗藥消失", fontsize=25, weight="bold", va="top")
    event, states = events[replicate], trials[replicate]
    fig.text(.055, .912, f"{len(trials)} 次試驗的實際軌跡；突出案例 replicate {replicate}（seed={event['seed']}）。"
             "黑色箭頭＝由當代完整基因型組成算出的期望下一代位置。", fontsize=13.5, color=GRAY)
    start, relapse = event["therapy_generation"], event["relapse_generation"]
    colors = [BLUE if start is None or g < start else RED if relapse is not None and g >= relapse else ORANGE
              for g in range(len(states) - 1)]
    for path in trials.values():
        ax.plot(*phase_coordinates(path).T, color="#BCC7CE", alpha=.48, lw=.85, zorder=1)
    xy = phase_coordinates(states)
    for g, color in enumerate(colors):
        ax.plot(*xy[g:g+2].T, color=color, lw=2.6, zorder=3)
    ax.scatter(*xy.T, s=9, color=INK, alpha=.55, zorder=4, linewidths=0)
    # 期望箭頭：先求原始計數期望再取 log；箭頭頭端與下一個實際點的距離就是隨機偏差。
    for g in sorted({0, 5, 10, 13, 14, 15, 16, 17, 18, 19, 20, 22, 24, 28, 32, len(steps) - 1}):
        if 0 <= g < len(steps):
            row = steps[g]
            arrow(ax, xy[g], np.log10([row["expected_sensitive"] + 1, row["expected_resistant"] + 1]),
                  INK, size=11, linewidth=1.5, zorder=6)
    sensitive = np.r_[0, np.geomspace(1, config.therapy_start_size, 250)]
    resistant = np.maximum(0, config.therapy_start_size - sensitive)
    ax.plot(np.log10(sensitive + 1), np.log10(resistant + 1), color=GRAY, ls="--", lw=1, alpha=.7)
    ax.text(.65, 4.37, rf"$S+R={config.therapy_start_size:,}$（治療門檻）", color=GRAY, fontsize=11.5)
    positions = {"初始": (2.75, .56), "治療啟動": (3.64, 3.5), "負荷谷底": (1.05, 2.08), "首次復發": (1.27, 5.15)}
    points = event_generations(states, event)
    for label, g in points.items():
        marker = {"初始": "o", "治療啟動": "D", "負荷谷底": "v", "首次復發": "s"}[label]
        ax.scatter(*xy[g], s=62, marker=marker, c=INK, edgecolors="white", linewidth=.9, zorder=7)
        ax.annotate(f"{label}\n$t={g}$", xy=xy[g], xytext=positions[label], fontsize=11.5, ha="center",
                    arrowprops={"arrowstyle": "-", "color": GRAY, "lw": .9},
                    bbox={"facecolor": "white", "edgecolor": "none", "pad": 2}, zorder=8)
    ax.scatter(*xy[-1], s=50, facecolor="white", edgecolor=RED, linewidth=1.5, zorder=7)
    decorate_axis(ax, "敏感 × 抗藥細胞（log 座標）")
    handles = [Line2D([], [], color="#BCC7CE", lw=1, label=f"全部 {len(trials)} 次試驗"),
               Line2D([], [], color=BLUE, lw=3, label="治療前"),
               Line2D([], [], color=ORANGE, lw=3, label="治療至首次復發"),
               Line2D([], [], color=RED, lw=3, label="復發後"),
               Line2D([], [], color=INK, lw=1.5, marker=">", markersize=7, label="期望下一代")]
    fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(.055, .875), frameon=False, ncol=5, fontsize=12)
    # 右側：案例數字。期望與實際並列，說明大族群幾乎沒有偏差、小族群偏差變大。
    y = .97
    side.text(0, y, "讀圖", fontsize=16, weight="bold", va="top"); y -= .08
    for text in ("治療前：S、R 一起長，路徑往右上。",
                 "治療後：S 每代少約一個數量級，\n路徑幾乎水平往左。",
                 "S 歸零後：沿縱軸直直往上，\n純抗藥族群在生長。"):
        side.text(0, y, text, fontsize=12.5, va="top"); y -= .035 + .05 * (1 + text.count("\n"))
    y -= .02
    if start is not None and "負荷谷底" in points:
        nadir = states[points["負荷谷底"]]
        side.text(0, y, f"抗藥細胞在治療中從未減少", fontsize=14, weight="bold", va="top", color=RED); y -= .065
        side.text(0, y, f"治療起點 R = {states[start]['resistant']:,}\n谷底 R = {nadir['resistant']:,}"
                  f"（占 {100*nadir['resistant']/nadir['total']:.1f}%）", fontsize=12.5, va="top"); y -= .14
    nadir_g = points.get("負荷谷底")
    for g, title in ((start, "治療啟動"), (nadir_g, "負荷谷底")):
        if g is not None and g < len(steps):
            row = steps[g]
            side.text(0, y, f"{title} t={g}：(S, R) = ({row['sensitive']:,}, {row['resistant']:,})",
                      fontsize=12.5, weight="bold", va="top"); y -= .055
            side.text(0, y, f"期望下一代 ({row['expected_sensitive']:,.1f}, {row['expected_resistant']:,.1f})\n"
                      f"實際下一代 ({row['observed_next_sensitive']:,}, {row['observed_next_resistant']:,})",
                      fontsize=12, va="top", color=GRAY); y -= .12
    fig.text(.07, .115, r"箭頭端點 $=\phi(\mathrm{E}[S_{t+1}\mid s],\mathrm{E}[R_{t+1}\mid s])$，"
             r"$\phi(S,R)=(\log_{10}(S+1),\log_{10}(R+1))$；公式推導見圖 07。", fontsize=12.5)
    fig.text(.07, .07, "族群越小，實際下一步偏離期望越多（人口學隨機性）。座標含 +1，因此零細胞可顯示；空心終點是計算上限，不是平衡點。",
             fontsize=11.5, color=GRAY)
    fig.text(.07, .035, "相同 (S, R) 底下的 driver／TSG 組成可能不同，因此這是完整狀態的投影，不是封閉的二維向量場。",
             fontsize=11.5, color=GRAY)
    return fig


def save_phase_figure(fig, destination, name):
    for extension in ("png", "svg"):
        fig.savefig(destination / f"{name}.{extension}", dpi=300, facecolor="white")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=ROOT / "outputs")
    parser.add_argument("--output", type=Path, default=ROOT / "outputs/phase_plots")
    parser.add_argument("--replicate", type=int, help="指定要突出的實際試驗，預設接近中位結果")
    args = parser.parse_args()
    config, trials, events = load_trials(args.input)
    rep = choose_exemplar(events) if args.replicate is None else args.replicate
    if rep not in trials:
        raise ValueError(f"replicate {rep} 不存在；可選 {sorted(trials)}")
    result = simulate(config, events[rep]["seed"], therapy=True)
    assert len(result.history) == len(trials[rep])
    for actual, reference in zip(result.history, trials[rep]):
        assert all(actual[k] == reference[k] for k in reference), (actual, reference)
    outcome = therapy_outcome(result)
    assert outcome["nadir"] == events[rep]["nadir"]
    assert outcome["relapse_generation"] == events[rep]["relapse_generation"]
    assert result.therapy_generation == events[rep]["therapy_generation"]
    steps = conditional_steps(result)
    # 沒有敏感父代且突變不可逆時，下一代敏感細胞數的期望必須為零。
    assert all(row["expected_sensitive"] == 0 for row in steps if row["sensitive"] == 0)
    args.output.mkdir(parents=True, exist_ok=True)
    with (args.output / "conditional_phase_steps.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(steps[0]))
        writer.writeheader(); writer.writerows(steps)
    points = event_generations(result.history, events[rep])
    details = {label: {**result.history[g], "resistant_percent": 100*result.history[g]["resistant"]/result.history[g]["total"]
                       if result.history[g]["total"] else None} for label, g in points.items()}
    summary = {"source": str(args.input), "trial_count": len(trials), "observed_rows": sum(map(len, trials.values())),
               "highlight_replicate": rep, "highlight_seed": result.seed,
               "selection": "nearest median relapse generation; then nearest median nadir on log scale; replicate breaks ties",
               "stop_reason": result.stop_reason, "events": details,
               "coordinates": {"A": "(log10(S+1), log10(R+1))", "B": "(log10(T+1), 100 R/T), missing at T=0"},
               "conditional_arrows": "raw-count conditional mean from full genotype state, transformed afterward; not E[log(count+1)]",
               "limitations": ["2D projection is not generally a closed Markov state", "no interpolated autonomous vector field",
                               "finite observed trajectories do not establish attractors or stable equilibria",
                               "cell_limit is a computational stop, not carrying capacity"],
               "verification": "highlight reproduced against every observed CSV count; event timing and expectation conservation verified"}
    (args.output / "phase_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2)+"\n")
    set_phase_style()
    save_phase_figure(plot_phase_portrait(trials, events, rep, config, result, steps), args.output, "01_phase_portrait")
    print(json.dumps({"replicate": rep, "seed": result.seed, "trial_count": len(trials),
                      "verified_generations": len(result.history), "output": str(args.output)}, ensure_ascii=False))


if __name__ == "__main__":
    main()
