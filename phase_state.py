"""相圖與基因型演化圖共用：重現治療試驗、計算狀態特徵、雙語標籤。

只使用基礎模型（driver、TSG、抗藥、治療）；IL11 必須關閉。
重現的每次試驗都與 outputs/therapy_trajectories.csv 逐代比對。
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

from modeling_cancer import (Config, Result, driver_count, growth_probabilities,
                             simulate, suppressor_alleles, suppressor_loss)
from plot_phase_portraits import projected_offspring_mean, set_phase_style
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parent
INK, GRAY, GRID = "#173042", "#687783", "#DDE4E8"
LANGS = ("zh", "en")


def set_style(lang: str):
    """中文沿用既有相圖字型；英文用 DejaVu Sans，兩者數學式皆為 Computer Modern。"""
    if lang == "zh":
        set_phase_style()
    else:
        plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 13,
                             "mathtext.fontset": "cm", "axes.unicode_minus": True,
                             "axes.spines.top": False, "axes.spines.right": False,
                             "axes.labelcolor": INK, "text.color": INK,
                             "axes.edgecolor": GRAY, "xtick.color": GRAY, "ytick.color": GRAY,
                             "figure.facecolor": "white", "axes.facecolor": "white",
                             "svg.fonttype": "path"})


def save(fig, directory: Path, name: str, lang: str):
    directory.mkdir(parents=True, exist_ok=True)
    for extension in ("png", "svg"):
        fig.savefig(directory / f"{name}_{lang}.{extension}", dpi=300, facecolor="white")
    plt.close(fig)


def gene_label(bit: int, config: Config) -> str:
    """位元 → 基因名稱；基因以 1 起算，TSG 第二等位基因位於 n_genes 之後。"""
    if bit in config.driver_genes:
        return f"ONC{config.driver_genes.index(bit) + 1}"
    if bit == config.resistance_gene:
        return "RES"
    if bit in config.suppressor_genes:
        return f"TSG{config.suppressor_genes.index(bit) + 1}a"
    if bit >= config.n_genes:
        return f"TSG{bit - config.n_genes + 1}b"
    return f"P{bit + 1}"


def genotype_label(genotype: int, config: Config) -> str:
    bits = [b for b in range(config.n_genes + len(config.suppressor_genes)) if genotype & (1 << b)]
    functional = [gene_label(b, config) for b in bits if not gene_label(b, config).startswith("P")]
    passengers = [gene_label(b, config) for b in bits if gene_label(b, config).startswith("P")]
    return "+".join(functional + passengers) or "WT"


def load_config(directory: Path) -> tuple[Config, dict]:
    metadata = json.loads((directory / "run_metadata.json").read_text())
    settings = metadata["config"].copy()
    for key in ("driver_genes", "suppressor_genes"):
        settings[key] = tuple(settings[key])
    config = Config(**settings)
    config.validate()
    if config.il11_enabled:
        raise ValueError("相圖只分析基礎模型；請用 IL11 關閉的輸出重新執行。")
    return config, metadata


def reproduce_trials(directory: Path) -> tuple[Config, list[Result], list[dict]]:
    """以記錄的 seed 重跑治療試驗，取得完整基因型快照；計數必須與 CSV 完全一致。"""
    config, metadata = load_config(directory)
    with (directory / "therapy_trajectories.csv").open(newline="") as handle:
        recorded = {}
        for row in csv.DictReader(handle):
            recorded.setdefault(int(row["replicate"]), []).append(row)
    with (directory / "therapy_replicates.csv").open(newline="") as handle:
        events = [{**row, **{k: int(row[k]) if row[k] else None
                             for k in ("replicate", "seed", "therapy_generation", "nadir", "relapse_generation")}}
                  for row in csv.DictReader(handle)]
    results = []
    for event in sorted(events, key=lambda e: e["replicate"]):
        result = simulate(config, event["seed"], therapy=True)
        rows = [r for r in recorded[event["replicate"]] if r["extended_extinction"] != "True"]
        assert len(rows) == len(result.history), event["replicate"]
        for actual, row in zip(result.history, rows):
            assert all(actual[k] == int(row[k]) for k in ("total", "sensitive", "resistant"))
        assert result.therapy_generation == event["therapy_generation"]
        results.append(result)
    assert len(results) == metadata["replicates"]
    return config, results, events


def fitness(genotype: int, config: Config, therapy: bool) -> float:
    """每個父細胞的期望子代數 w = (1-d)(1-k)(2+q)；k 只在治療時作用。"""
    death, triple = growth_probabilities(genotype, config)
    kill = 0.0
    if therapy:
        resistant = bool(genotype & (1 << config.resistance_gene))
        kill = config.resistant_therapy_kill if resistant else config.sensitive_therapy_kill
    return (1 - death) * (1 - kill) * (2 + triple)


def expected_next(population: dict[int, float], config: Config, therapy: bool) -> dict[str, float]:
    """一步原始計數期望：S、R 與子代（突變後）的平均無治療適應度 w̄_off。

    w_off 只取決於 oncogene 是否活化與 TSG 是否失活；兩者及抗藥位元的突變彼此獨立，
    因此可對每個父代基因型精確計算子代的期望 w_off，不需要展開所有子代基因型。
    """
    mu = config.mutation_rate
    gain_onc = 1 - (1 - mu) ** len(config.driver_genes)
    before = {"S": 0.0, "R": 0.0}
    weighted_w = 0.0
    for genotype, count in population.items():
        offspring = count * fitness(genotype, config, therapy)
        resistant = bool(genotype & (1 << config.resistance_gene))
        before["R" if resistant else "S"] += offspring
        onc = driver_count(genotype, config) >= config.driver_hits
        p_onc = 1.0 if onc else (gain_onc if config.driver_hits == 1 else 0.0)
        # 每個 TSG 尚缺 h 個等位基因時，一代內補齊的機率為 μ^h；不同 TSG 獨立。
        keep = 1.0
        for gene in config.suppressor_genes:
            missing = max(0, config.suppressor_hits - suppressor_alleles(genotype, gene, config))
            keep *= 1 - mu ** missing
        p_tsg = 1 - keep
        w = {(o, t): (1 - max(0, config.death_rate - config.driver_death_reduction * o
                              - config.suppressor_death_reduction * t))
             * (2 + config.triple_probability + config.driver_triple_bonus * o
                + config.suppressor_triple_bonus * t)
             for o in (0, 1) for t in (0, 1)}
        mean_w = sum(w[o, t] * (p_onc if o else 1 - p_onc) * (p_tsg if t else 1 - p_tsg)
                     for o in (0, 1) for t in (0, 1))
        weighted_w += offspring * mean_w
    total = before["S"] + before["R"]
    return {"S": (1 - mu) * before["S"], "R": before["R"] + mu * before["S"], "N": total,
            "w_off": weighted_w / total if total else np.nan}


def logit_rho(resistant, sensitive):
    """經驗 logit：log10((R+½)/(S+½))；ρ=0 或 1 仍有限，且等於 log R − log S 的平滑版本。"""
    return np.log10((np.asarray(resistant, float) + .5) / (np.asarray(sensitive, float) + .5))


def state_features(population: dict[int, int], config: Config, therapy: bool) -> dict:
    """由完整基因型快照計算所有候選相座標；IL11 不在基礎模型中，故不列入。"""
    genotypes = list(population)
    counts = np.array(list(population.values()), float)
    total = counts.sum()
    p = counts / total
    resistance = 1 << config.resistance_gene
    passengers = [g for g in range(config.n_genes) if g not in
                  config.driver_genes + config.suppressor_genes + (config.resistance_gene,)]
    res = np.array([bool(g & resistance) for g in genotypes])
    onc = np.array([driver_count(g, config) >= config.driver_hits for g in genotypes])
    tsg_loss = np.array([suppressor_loss(g, config) for g in genotypes])
    tsg_any = np.array([any(suppressor_alleles(g, k, config) for k in config.suppressor_genes) for g in genotypes])
    w_off = np.array([fitness(g, config, False) for g in genotypes])
    w_on = np.array([fitness(g, config, True) for g in genotypes])
    sensitive, resistant = counts[~res].sum(), counts[res].sum()
    mean_off = (p * w_off).sum()
    return {
        "log_N": np.log10(total + 1), "log_S": np.log10(sensitive + 1), "log_R": np.log10(resistant + 1),
        "frac_R": p[res].sum(), "logit_R": float(logit_rho(resistant, sensitive)),
        "frac_onc": p[onc].sum(),
        "mean_drivers": (p * [driver_count(g, config) for g in genotypes]).sum(),
        "frac_tsg_anyhit": p[tsg_any].sum(), "frac_tsg_loss": p[tsg_loss].sum(),
        "frac_onc_and_tsg": p[onc & tsg_loss].sum(),
        "mean_passengers": (p * [sum(bool(g & (1 << k)) for k in passengers) for g in genotypes]).sum(),
        "log_richness": np.log10(len(genotypes)),
        "shannon": float(-(p * np.log(p)).sum()), "simpson": float(1 - (p ** 2).sum()),
        "mean_w_off": mean_off, "var_w_off": (p * (w_off - mean_off) ** 2).sum(),
        "mean_w_on": (p * w_on).sum(), "mean_w_now": (p * (w_on if therapy else w_off)).sum(),
        "therapy_on": float(therapy),
    }


def trial_feature_rows(config: Config, results: list[Result]) -> list[dict]:
    """每個 (試驗, 世代) 一列；growth_next 是實際下一代 ln 成長率，作為動力學目標而非座標。"""
    rows = []
    for replicate, result in enumerate(results):
        start = result.therapy_generation
        for generation, population in enumerate(result.populations):
            if not population:
                continue
            therapy = start is not None and generation >= start
            row = {"replicate": replicate, "generation": generation,
                   **state_features(population, config, therapy),
                   "gens_since_therapy": float(generation - start) if therapy else 0.0}
            total = sum(population.values())
            following = result.populations[generation + 1] if generation + 1 < len(result.populations) else None
            row["growth_next"] = np.log(sum(following.values()) / total) if following else np.nan
            rows.append(row)
    return rows


def check_expectation(config: Config) -> None:
    """expected_next 的 S/R 必須與既有 projected_offspring_mean 完全一致。"""
    population = {0: 700, 1: 200, 1 << config.resistance_gene: 90,
                  (1 << config.resistance_gene) | (1 << config.suppressor_genes[0]): 10}
    for therapy in (False, True):
        ours = expected_next(population, config, therapy)
        reference = projected_offspring_mean(population, config, therapy)
        assert np.isclose(ours["S"], reference["expected_sensitive"], rtol=1e-12)
        assert np.isclose(ours["R"], reference["expected_resistant"], rtol=1e-12)
    # μ=0 時子代 w̄_off 必須等於父代以期望子代數加權的 fitness 平均。
    frozen = Config(**{**config.__dict__, "mutation_rate": 0.0})
    ours = expected_next(population, frozen, False)
    weights = {g: n * fitness(g, frozen, False) for g, n in population.items()}
    direct = sum(w * fitness(g, frozen, False) for g, w in weights.items()) / sum(weights.values())
    assert np.isclose(ours["w_off"], direct, rtol=1e-12)
