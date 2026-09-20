"""BME 6120：離散分枝腫瘤模型及展示圖表。

執行：.venv/bin/python modeling_cancer.py
模型以基因型計數（不是追蹤家系），所有參數均為教學假設，不是臨床估計。
每代每個細胞死亡或產生 2 / 3 個子細胞；子細胞獨立發生不可逆突變。
「兩擊」指兩個不同 driver genes 突變，並非同一基因的雙等位基因失活。
"""

from __future__ import annotations

import argparse
import csv
import json
import os
from dataclasses import asdict, dataclass, replace
from pathlib import Path

# 圖形快取留在專案中，不寫入 conda 或使用者的 matplotlib 設定目錄。
os.environ.setdefault("MPLCONFIGDIR", str(Path(__file__).resolve().parent / ".mplconfig"))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import numpy as np


@dataclass(frozen=True)
class Config:
    n_genes: int = 12
    driver_genes: tuple[int, ...] = (0, 1, 2)
    resistance_gene: int = 3
    driver_hits: int = 1
    mutation_rate: float = 0.001  # 每個子細胞、每個未突變基因、每代
    death_rate: float = 0.30
    driver_death_reduction: float = 0.08
    triple_probability: float = 0.04  # 存活後產生三個子細胞的機率
    driver_triple_bonus: float = 0.12
    initial_cells: int = 100
    generations: int = 75
    therapy_start_size: int = 20_000
    sensitive_therapy_kill: float = 0.90  # 自然死亡後的額外殺傷機率
    resistant_therapy_kill: float = 0.04
    max_cells: int = 1_000_000  # 計算停止門檻，不是環境承載量

    def validate(self):
        if not 2 <= self.n_genes <= 20:
            raise ValueError("n_genes 必須介於 2–20")
        genes = self.driver_genes + (self.resistance_gene,)
        if len(set(genes)) != len(genes) or any(g < 0 or g >= self.n_genes for g in genes):
            raise ValueError("driver 與 resistance 基因必須互異且位於基因組內")
        if not 1 <= self.driver_hits <= len(self.driver_genes):
            raise ValueError("driver_hits 超出可用 driver 數目")
        for name in ("mutation_rate", "death_rate", "driver_death_reduction",
                     "triple_probability", "driver_triple_bonus",
                     "sensitive_therapy_kill", "resistant_therapy_kill"):
            if not 0 <= getattr(self, name) <= 1:
                raise ValueError(f"{name} 必須介於 0–1")
        if self.triple_probability + self.driver_triple_bonus > 1:
            raise ValueError("產生三個子細胞的總機率不能超過 1")
        if min(self.initial_cells, self.generations, self.therapy_start_size, self.max_cells) <= 0:
            raise ValueError("細胞數、代數與門檻必須為正")


@dataclass
class Result:
    config: Config
    seed: int
    history: list[dict]
    populations: list[dict[int, int]]
    therapy_generation: int | None
    stop_reason: str


def driver_count(genotype: int, config: Config) -> int:
    return sum(bool(genotype & (1 << g)) for g in config.driver_genes)


def mutate_offspring(population: dict[int, int], config: Config,
                     rng: np.random.Generator) -> dict[int, int]:
    """逐基因二項抽樣：與逐細胞獨立突變等價，允許一代多基因突變。"""
    population = dict(population)
    for gene in range(config.n_genes):
        bit = 1 << gene
        updated = dict(population)
        for genotype, count in population.items():
            if genotype & bit:
                continue
            mutants = int(rng.binomial(count, config.mutation_rate))
            if mutants:
                updated[genotype] -= mutants
                updated[genotype | bit] = updated.get(genotype | bit, 0) + mutants
        population = {g: n for g, n in updated.items() if n > 0}
    return population


def summarize(population: dict[int, int], generation: int, config: Config) -> dict:
    total = sum(population.values())
    fractions = np.array(list(population.values()), dtype=float) / max(total, 1)
    resistant = sum(n for g, n in population.items() if g & (1 << config.resistance_gene))
    return dict(generation=generation, total=total, resistant=resistant,
                sensitive=total - resistant, richness=len(population),
                shannon=float(-np.sum(fractions * np.log(fractions))),
                simpson=float(1 - np.sum(fractions ** 2)) if total else 0.0,
                driver_fraction=sum(n for g, n in population.items()
                                    if driver_count(g, config) >= config.driver_hits) / max(total, 1))


def simulate(config: Config = Config(), seed: int = 42, therapy: bool = True,
             target_size: int | None = None) -> Result:
    """同步世代分枝過程；治療從首次達門檻後的轉移開始，持續至模擬結束。

    超過 max_cells 時保留完整該代，不截斷細胞；target_size 用於等大小比較。
    """
    config.validate()
    if target_size is not None and target_size <= 0:
        raise ValueError("target_size 必須為正")
    rng = np.random.default_rng(seed)
    population = {0: config.initial_cells}
    history, populations = [], []
    therapy_generation = None
    stop_reason = "generation_limit"
    for generation in range(config.generations + 1):
        row = summarize(population, generation, config)
        history.append(row)
        populations.append(dict(population))
        if row["total"] == 0:
            stop_reason = "extinction"
            break
        if target_size is not None and row["total"] >= target_size:
            stop_reason = "target_reached"
            break
        if row["total"] >= config.max_cells:
            stop_reason = "cell_limit"
            break
        if generation == config.generations:
            break
        if therapy and therapy_generation is None and row["total"] >= config.therapy_start_size:
            therapy_generation = generation
        offspring = {}
        for genotype, count in population.items():
            advantage = driver_count(genotype, config) >= config.driver_hits
            death = max(0, config.death_rate - config.driver_death_reduction * advantage)
            if therapy_generation is not None:
                resistant = bool(genotype & (1 << config.resistance_gene))
                kill = config.resistant_therapy_kill if resistant else config.sensitive_therapy_kill
                death = 1 - (1 - death) * (1 - kill)
            survivors = int(rng.binomial(count, 1 - death))
            triples = int(rng.binomial(survivors, config.triple_probability +
                                      config.driver_triple_bonus * advantage))
            n = 2 * survivors + triples
            if n:
                offspring[genotype] = n
        population = mutate_offspring(offspring, config, rng)
    return Result(config, seed, history, populations, therapy_generation, stop_reason)


def therapy_outcome(result: Result) -> dict:
    """復發定義：先降至治療前負荷以下，再回到治療開始時負荷。"""
    start = result.therapy_generation
    if start is None:
        return dict(status="therapy_not_started", nadir=None, relapse_generation=None)
    baseline = result.history[start]["total"]
    after = result.history[start + 1:]
    declined = False
    relapse = None
    for row in after:
        declined |= row["total"] < baseline
        if declined and row["total"] >= baseline:
            relapse = row["generation"]
            break
    status = ("relapse" if relapse is not None else "extinction" if result.stop_reason == "extinction"
              else "no_relapse_observed" if declined else "no_initial_response")
    return dict(status=status, nadir=min((r["total"] for r in after), default=baseline),
                relapse_generation=relapse)


COLORS = ["#167D9A", "#ED7953", "#6C5B9C", "#57A773", "#D6A43B", "#BC5B86"]


def set_plot_style():
    """統一簡報用視覺樣式；英文圖標籤避免依賴中文系統字型。"""
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11,
                         "axes.titlesize": 13, "axes.titleweight": "bold",
                         "axes.spines.top": False, "axes.spines.right": False,
                         "axes.labelcolor": "#334155", "text.color": "#172B42",
                         "axes.grid": True, "grid.alpha": 0.16,
                         "figure.facecolor": "#FAFBFD", "axes.facecolor": "white",
                         "savefig.facecolor": "#FAFBFD", "lines.linewidth": 2.3})


def save_figure(fig, output: Path, name: str):
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    for extension in ("png", "svg"):
        fig.savefig(output / f"{name}.{extension}", dpi=240, bbox_inches="tight")
    plt.close(fig)


def plot_dashboard(result: Result, output: Path):
    """輸出生長、基因型組成、多樣性與突變頻率四面板圖。"""
    set_plot_style()
    fig, axes = plt.subplots(2, 2, figsize=(13.5, 9), layout="constrained")
    fig.suptitle("Tumor evolution | treatment & resistance", fontsize=21, fontweight="bold")
    t = np.array([r["generation"] for r in result.history])
    ax = axes[0, 0]
    for key, label, color in (("total", "Total", "#172B42"),
                              ("sensitive", "Sensitive", COLORS[0]),
                              ("resistant", "Resistant", COLORS[1])):
        ax.plot(t, [r[key] for r in result.history], label=label, color=color)
    ax.set(yscale="symlog", ylabel="Cells (symlog; linear below 1)", title="A  Population response")
    ax.set_yscale("symlog", linthresh=1)
    if result.therapy_generation is not None:
        ax.axvspan(result.therapy_generation, t[-1], color=COLORS[1], alpha=0.09)
        ax.axvline(result.therapy_generation, color=COLORS[1], ls="--", lw=1.5, label="Therapy begins")
    ax.legend(frameon=False, fontsize=9)
    ax = axes[0, 1]
    totals = np.array([r["total"] for r in result.history])
    scores = {}
    for pop, total in zip(result.populations, totals):
        for g, n in pop.items():
            scores[g] = scores.get(g, 0) + n / max(total, 1)
    top = sorted(scores, key=scores.get, reverse=True)[:5]
    bands = np.array([[p.get(g, 0) / max(n, 1) for p, n in zip(result.populations, totals)] for g in top])
    others = np.maximum(0, (totals > 0).astype(float) - bands.sum(axis=0))
    labels = [f"G{g:03X}" + (" (R)" if g & (1 << result.config.resistance_gene) else "") for g in top]
    ax.stackplot(t, *bands, others, labels=labels + ["Other"], colors=COLORS, alpha=0.9)
    ax.set(ylim=(0, 1), ylabel="Fraction of cells", title="B  Genotype composition")
    ax.legend(loc="upper left", fontsize=8, frameon=False, ncol=3)
    ax = axes[1, 0]
    ax.plot(t, [np.exp(r["shannon"]) if r["total"] else 0 for r in result.history], color=COLORS[2])
    ax.set(ylabel="Effective genotypes: exp(Shannon H)", title="C  Genetic heterogeneity")
    ax = axes[1, 1]
    frequency = np.array([[sum(n for g, n in p.items() if g & (1 << gene)) / max(total, 1)
                           for p, total in zip(result.populations, totals)]
                          for gene in range(result.config.n_genes)])
    mesh = ax.pcolormesh(np.arange(len(t) + 1) - 0.5,
                         np.arange(result.config.n_genes + 1) - 0.5,
                         frequency, cmap="YlGnBu", vmin=0, vmax=1, shading="flat")
    labels = [f"Gene {g + 1}" + (" · driver" if g in result.config.driver_genes else
                                  " · resistance" if g == result.config.resistance_gene else "")
              for g in range(result.config.n_genes)]
    ax.set(yticks=range(result.config.n_genes), yticklabels=labels, title="D  Mutation frequencies")
    ax.grid(False)
    fig.colorbar(mesh, ax=ax, label="Mutant fraction", pad=0.02)
    for ax in axes.flat:
        ax.set_xlabel("Generation")
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    save_figure(fig, output, "01_tumor_dashboard")


def sample_at_size(result: Result, size: int, rng: np.random.Generator) -> dict | None:
    """跨門檻的該代中無放回取固定數目，避免不同生長速度的 overshoot 偏差。

    這是等大小的腫瘤樣本，不是連續時間中恰好達到該大小的完整腫瘤。
    """
    if result.history[-1]["total"] < size:
        return None
    pop = result.populations[-1]
    counts = rng.multivariate_hypergeometric(np.array(list(pop.values()), dtype=np.int64), size)
    sampled = {g: int(n) for g, n in zip(pop, counts) if n}
    return summarize(sampled, result.history[-1]["generation"], result.config)


def run_comparisons(config: Config, repeats: int, seed: int, output: Path) -> list[dict]:
    """重複模擬：驅動突變對生長，以及死亡率對等大小樣本多樣性的影響。"""
    if repeats < 1:
        raise ValueError("repeats 必須至少為 1")
    set_plot_style()
    fig, axes = plt.subplots(1, 2, figsize=(13.5, 5.4), layout="constrained")
    fig.suptitle("What changes tumor growth and diversity?", fontsize=20, fontweight="bold")
    rows = []
    growth_horizon = 18
    for condition, label in enumerate(("Neutral mutations", "One driver required", "Two drivers required")):
        cfg = replace(config, generations=growth_horizon, max_cells=10**12, driver_hits=2 if condition == 2 else 1)
        if condition == 0:
            cfg = replace(cfg, driver_death_reduction=0, driver_triple_bonus=0)
        trajectories = []
        for rep in range(repeats):
            run_seed = seed + condition * 10000 + rep
            result = simulate(cfg, seed=run_seed, therapy=False)
            values = [r["total"] for r in result.history]
            if result.stop_reason == "cell_limit":
                raise RuntimeError("Growth comparison hit safety limit; shorten growth_horizon")
            values.extend([0] * (growth_horizon + 1 - len(values)))
            trajectories.append(values)
            rows.append(dict(experiment="growth", condition=label, replicate=rep,
                             seed=run_seed, reached_target="", generation=result.history[-1]["generation"],
                             total=result.history[-1]["total"], shannon=result.history[-1]["shannon"]))
        low, median, high = np.percentile(trajectories, [10, 50, 90], axis=0)
        axes[0].plot(range(growth_horizon + 1), median, color=COLORS[condition], label=label)
        axes[0].fill_between(range(growth_horizon + 1), low, high, color=COLORS[condition], alpha=0.14)
    axes[0].set(title="A  Growth advantage", xlabel="Generation", ylabel="Cells · median & 10–90% interval")
    axes[0].set_yscale("symlog", linthresh=1)
    axes[0].legend(frameon=False, fontsize=9)
    target = 5000
    for index, death in enumerate((0.20, 0.35, 0.45)):
        values = []
        for rep in range(repeats):
            run_seed = seed + 50000 + index * 10000 + rep
            cfg = replace(config, death_rate=death, generations=100, max_cells=max(target * 4, config.max_cells))
            result = simulate(cfg, seed=run_seed, therapy=False, target_size=target)
            sample = sample_at_size(result, target, np.random.default_rng(run_seed + 1_000_000))
            if sample is not None:
                values.append(sample["shannon"])
            rows.append(dict(experiment="death_rate", condition=death, replicate=rep,
                             seed=run_seed, reached_target=sample is not None,
                             generation=result.history[-1]["generation"],
                             total=result.history[-1]["total"], shannon=sample["shannon"] if sample else ""))
        if values:
            jitter = np.random.default_rng(seed + index).uniform(-0.10, 0.10, len(values))
            axes[1].scatter(index + jitter, values, color=COLORS[index], s=34, alpha=0.7)
            axes[1].plot([index - 0.18, index + 0.18], [np.median(values)] * 2, color="#172B42", lw=3)
        axes[1].text(index, 0.98, f"{len(values)}/{repeats} reached", ha="center", va="top",
                     transform=axes[1].get_xaxis_transform(), fontsize=9)
    axes[1].set(title=f"B  Diversity in {target:,}-cell samples", xlabel="Natural death probability per generation",
                ylabel="Shannon diversity (natural log)", xticks=[0, 1, 2], xticklabels=["0.20", "0.35", "0.45"])
    axes[1].margins(y=0.25, x=0.20)
    save_figure(fig, output, "02_scenario_comparisons")
    write_csv(output / "comparisons.csv", rows)
    return rows


def write_csv(path: Path, rows: list[dict]):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--replicates", type=int, default=12)
    parser.add_argument("--output", type=Path, default=Path("outputs"))
    args = parser.parse_args()
    if args.replicates < 1:
        parser.error("--replicates must be >= 1")
    args.output.mkdir(parents=True, exist_ok=True)
    config = Config()
    result = simulate(config, args.seed)
    plot_dashboard(result, args.output)
    write_csv(args.output / "trajectory.csv", result.history)
    run_comparisons(config, args.replicates, args.seed, args.output)
    # 抗藥性為隨機事件，額外報告重複試驗結果，避免只展示成功復發案例。
    outcomes = []
    for rep in range(args.replicates):
        trial = simulate(config, args.seed + 100000 + rep)
        outcomes.append(dict(replicate=rep, seed=trial.seed,
                             therapy_generation=trial.therapy_generation,
                             stop_reason=trial.stop_reason, **therapy_outcome(trial)))
    write_csv(args.output / "therapy_replicates.csv", outcomes)
    metadata = dict(config=asdict(config), seed=args.seed, replicates=args.replicates,
                    stop_reason=result.stop_reason, therapy_generation=result.therapy_generation,
                    outcome=therapy_outcome(result), numpy_version=np.__version__,
                    matplotlib_version=matplotlib.__version__)
    (args.output / "run_metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    print(json.dumps(metadata, indent=2))
    print(f"Figures and data saved to: {args.output.resolve()}")


if __name__ == "__main__":
    main()
