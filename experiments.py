"""實驗：在不同條件下重複執行模型，彙整結果，並讀寫資料檔。

每個函數只做計算或讀寫，回傳資料給 main.py；畫圖在 plot_*.py。
主實驗與對照實驗的 seed 區段：治療 +100000、TSG 對照 +0／+10000／+20000、
死亡率 +50000、IL11 +300000、治療時機 +400000、敏感度 +500000、適應性治療 +600000。
"""

from __future__ import annotations

import csv
import json
from concurrent.futures import ProcessPoolExecutor
from dataclasses import asdict, replace
from importlib.metadata import version
from pathlib import Path

import numpy as np

from analysis import (
    LINEAR_SCALE,
    LOG_SCALE,
    MIN_TRANSITIONS,
    MIN_TRIALS,
    build_fields,
    conditional_steps,
    event_generations,
    fitness,
    genotype_label,
    goldie_coldman_theory,
    il11_benefit_curve,
    muller_layout,
    spearman,
    TOP_K,
    trial_feature_rows,
    wilson,
)
from simulation import (
    Config,
    Result,
    crowding_death,
    growth_probabilities,
    simulate,
    summarize,
    therapy_outcome,
)

ROOT = Path(__file__).resolve().parent
OUTPUTS = ROOT / "outputs"


# ---------------------------------------------------------------- 資料讀寫 ----------
def write_csv(path: Path, rows: list[dict]):
    # 以第一筆資料的鍵建立欄名；UTF-8 支援中文，newline 避免 CSV 多餘空行。
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def write_json(
    path: Path, data, *, ensure_ascii=False, trailing_newline=True, default=None
):
    """縮排 2 格的 JSON；預設保留中文並在結尾換行。"""
    text = json.dumps(data, ensure_ascii=ensure_ascii, indent=2, default=default)
    path.write_text(text + ("\n" if trailing_newline else ""), encoding="utf-8")


def read_csv(path):
    with Path(path).open(newline="") as handle:
        return list(csv.DictReader(handle))


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
            **row,
            "replicate": int(row["replicate"]),
            "seed": int(row["seed"]),
            **{
                key: int(row[key]) if row[key] else None
                for key in ("therapy_generation", "relapse_generation", "nadir")
            },
        }
    trials = {rep: [] for rep in events}
    for row in read_csv(directory / "therapy_trajectories.csv"):
        if row["extended_extinction"] == "True":
            continue
        rep = int(row["replicate"])
        state = {
            key: int(row[key])
            for key in ("generation", "total", "sensitive", "resistant")
        }
        assert state["sensitive"] + state["resistant"] == state["total"]
        assert min(state["sensitive"], state["resistant"]) >= 0
        assert int(row["seed"]) == events[rep]["seed"]
        trials[rep].append(state)
    for states in trials.values():
        states.sort(key=lambda row: row["generation"])
        assert [row["generation"] for row in states] == list(range(len(states)))
    assert len(trials) == metadata["replicates"]
    return config, trials, events


def reproduce_trials(directory: Path) -> tuple[Config, list[Result], list[dict]]:
    """以記錄的 seed 重跑治療試驗，取得完整基因型快照；計數必須與 CSV 完全一致。"""
    config, metadata = load_config(directory)
    with (directory / "therapy_trajectories.csv").open(newline="") as handle:
        recorded = {}
        for row in csv.DictReader(handle):
            recorded.setdefault(int(row["replicate"]), []).append(row)
    with (directory / "therapy_replicates.csv").open(newline="") as handle:
        events = [
            {
                **row,
                **{
                    k: int(row[k]) if row[k] else None
                    for k in (
                        "replicate",
                        "seed",
                        "therapy_generation",
                        "nadir",
                        "relapse_generation",
                    )
                },
            }
            for row in csv.DictReader(handle)
        ]
    results = []
    for event in sorted(events, key=lambda e: e["replicate"]):
        result = simulate(config, event["seed"], therapy=True)
        rows = [
            r
            for r in recorded[event["replicate"]]
            if r["extended_extinction"] != "True"
        ]
        assert len(rows) == len(result.history), event["replicate"]
        for actual, row in zip(result.history, rows):
            assert all(
                actual[k] == int(row[k]) for k in ("total", "sensitive", "resistant")
            )
        assert result.therapy_generation == event["therapy_generation"]
        results.append(result)
    assert len(results) == metadata["replicates"]
    return config, results, events


# ---------------------------------------------------------------- 主實驗 ------------
def run_therapy_trials(config: Config, replicates: int, seed: int) -> list[Result]:
    """主實驗：獨立重複試驗，達門檻後持續治療；第 rep 次的 seed 為 seed + 100000 + rep。"""
    return [simulate(config, seed + 100000 + rep) for rep in range(replicates)]


def therapy_ensemble(trials: list[Result]) -> tuple[list[dict], list[dict]]:
    """依絕對世代彙整治療試驗；滅絕延續零值，計算上限後視為缺失。

    回傳 (逐代摘要, 逐試驗逐代資料)。
    """
    # 不把因計算上限停止的族群當成滅絕，也不以最後數量冒充後續觀察。
    length = max(len(t.history) for t in trials)
    keys = (
        "total",
        "sensitive",
        "resistant",
        "effective_genotypes",
        "resistant_percent",
    )
    data = {k: np.full((len(trials), length), np.nan) for k in keys}
    raw = []
    for rep, trial in enumerate(trials):
        for generation in range(length):
            if generation >= len(trial.history) and trial.stop_reason != "extinction":
                continue
            row = trial.history[min(generation, len(trial.history) - 1)]
            total = row["total"]
            values = dict(
                total=total,
                sensitive=row["sensitive"],
                resistant=row["resistant"],
                effective_genotypes=np.exp(row["shannon"]) if total else 0,
                resistant_percent=100 * row["resistant"] / total if total else np.nan,
            )
            raw.append(
                dict(
                    replicate=rep,
                    seed=trial.seed,
                    generation=generation,
                    extended_extinction=generation >= len(trial.history),
                    **values,
                )
            )
            for key, value in values.items():
                data[key][rep, generation] = value
    # 每代輸出平均、樣本標準差、中位數、百分位與有效 n；n=1 時 SD 不定義。
    summary = []
    for generation in range(length):
        record = dict(generation=generation)
        for key, matrix in data.items():
            v = matrix[:, generation]
            v = v[np.isfinite(v)]
            record.update({
                f"{key}_n": len(v),
                f"{key}_mean": float(v.mean()) if len(v) else np.nan,
                f"{key}_sd": float(v.std(ddof=1)) if len(v) > 1 else np.nan,
            })
            for name, q in (("p10", 10), ("median", 50), ("p90", 90)):
                record[f"{key}_{name}"] = (
                    float(np.percentile(v, q)) if len(v) else np.nan
                )
        summary.append(record)
    return summary, raw


def therapy_events(trials: list[Result]) -> list[dict]:
    """每次試驗的治療起點、停止原因、谷底與復發世代；保留所有結果，不只挑復發案例。"""
    return [
        dict(
            replicate=rep,
            seed=trial.seed,
            therapy_generation=trial.therapy_generation,
            stop_reason=trial.stop_reason,
            **therapy_outcome(trial),
        )
        for rep, trial in enumerate(trials)
    ]


def run_metadata(
    config: Config,
    seed: int,
    replicates: int,
    trials: list[Result],
    outcomes: list[dict],
) -> dict:
    """完整參數、種子、套件版本與結果計數，方便追溯及重現本次執行。"""
    return dict(
        config=asdict(config),
        seed=seed,
        replicates=replicates,
        therapy_seeds=[t.seed for t in trials],
        therapy_outcome_counts={
            status: sum(o["status"] == status for o in outcomes)
            for status in sorted({o["status"] for o in outcomes})
        },
        interval="median and 10–90% trial percentiles; not confidence intervals",
        numpy_version=np.__version__,
        matplotlib_version=version("matplotlib"),
    )


# ---------------------------------------------------------------- 對照實驗 ----------
# TSG 比較的三種條件；第一個是對照組（TSG 失活沒有生長優勢）。
HIT_SCENARIOS = ("No TSG advantage", "TSG one-hit", "TSG two-hit")


def run_tsg_comparison(
    config: Config, repeats: int, seed: int, growth_horizon: int = 18
) -> tuple[list[dict], list[dict]]:
    """無 TSG 優勢（對照）、one-hit、two-hit 各重複 repeats 次，不治療、觀察 18 代。

    回傳 (逐代資料, 每次試驗終點)。
    """
    # 重複次數必須為正，才有可供比較的試驗分布。
    if repeats < 1:
        raise ValueError("repeats 必須至少為 1")
    rows, hit_rows = [], []
    for condition, label in enumerate(HIT_SCENARIOS):
        # 用 replace 建立新設定而不修改原始 Config；控制組只取消 TSG 優勢；
        # 所有情境保留相同 oncogene 效果與突變率。
        cfg = replace(
            config,
            generations=growth_horizon,
            max_cells=10**12,
            suppressor_hits=1 if condition == 1 else 2,
        )
        if condition == 0:
            cfg = replace(cfg, suppressor_death_reduction=0, suppressor_triple_bonus=0)
        for rep in range(repeats):
            run_seed = seed + condition * 10000 + rep
            result = simulate(cfg, seed=run_seed, therapy=False)
            # 若碰到安全上限就停止比較並報錯，以免把截短的曲線當成完整資料。
            if result.stop_reason == "cell_limit":
                raise RuntimeError(
                    "Growth comparison hit safety limit; shorten growth_horizon"
                )
            # 提早滅絕的試驗沿用最後一列（零細胞），使各次曲線長度一致。
            for generation in range(growth_horizon + 1):
                row = result.history[min(generation, len(result.history) - 1)]
                hit_rows.append(
                    dict(
                        scenario=label,
                        replicate=rep,
                        seed=run_seed,
                        generation=generation,
                        total=row["total"],
                        suppressor_fraction=row["suppressor_fraction"],
                    )
                )
            rows.append(
                dict(
                    experiment="growth",
                    condition=label,
                    replicate=rep,
                    seed=run_seed,
                    reached_target="",
                    generation=result.history[-1]["generation"],
                    total=result.history[-1]["total"],
                    shannon=result.history[-1]["shannon"],
                    suppressor_fraction=result.history[-1]["suppressor_fraction"],
                )
            )
    return hit_rows, rows


def hit_comparison_series(rows: list[dict]) -> dict:
    """把逐代資料整理成每個情境的 (試驗 × 世代) 矩陣，並計算第 18 代中位數。"""
    last = max(r["generation"] for r in rows)
    repeats = len({r["replicate"] for r in rows})
    order = [s for s in HIT_SCENARIOS if any(r["scenario"] == s for r in rows)]
    series, summaries = {}, []
    for label in order:
        # 依試驗與世代建立矩陣，各組均使用相同時間範圍及完整重複試驗。
        subset = [r for r in rows if r["scenario"] == label]
        trials = [
            [r for r in subset if r["replicate"] == rep] for rep in range(repeats)
        ]
        total = np.array([[r["total"] for r in trial] for trial in trials])
        fraction = 100 * np.array(
            [[r["suppressor_fraction"] for r in trial] for trial in trials]
        )
        series[label] = (total, fraction)
        summaries.append(
            dict(
                scenario=label,
                generation=last,
                replicates=repeats,
                median_cells=float(np.median(total[:, -1])),
                median_tsg_percent=float(np.median(fraction[:, -1])),
            )
        )
    return {
        "last": last,
        "repeats": repeats,
        "order": order,
        "series": series,
        "summaries": summaries,
    }


def sample_at_size(result: Result, size: int, rng: np.random.Generator) -> dict | None:
    """跨門檻的該代中無放回取固定數目，避免不同生長速度的 overshoot 偏差。

    這是等大小的腫瘤樣本，不是連續時間中恰好達到該大小的完整腫瘤。
    """
    # 最後一代未達取樣數量時回傳 None，由比較函數記錄為未達目標。
    if result.history[-1]["total"] < size:
        return None
    # 多變量超幾何分布等同從整個族群無放回取樣，維持固定樣本大小。
    pop = result.populations[-1]
    counts = rng.multivariate_hypergeometric(
        np.array(list(pop.values()), dtype=np.int64), size
    )
    # 重建樣本的基因型計數，沿用同一統計函數計算多樣性。
    sampled = {g: int(n) for g, n in zip(pop, counts) if n}
    return summarize(sampled, result.history[-1]["generation"], result.config)


def run_death_rate_comparison(
    config: Config,
    repeats: int,
    seed: int,
    target: int = 5000,
    deaths: tuple[float, ...] = (0.20, 0.35, 0.45),
) -> dict:
    """死亡率比較：首次達標後固定取樣 5,000 個細胞，降低跨門檻大小差異的影響。"""
    shannon, generations, rows = {}, {}, []
    for index, death in enumerate(deaths):
        shannon[death], generations[death] = [], []
        for rep in range(repeats):
            # 使用與生長比較不同的種子區段；取樣也使用獨立的亂數產生器。
            run_seed = seed + 50000 + index * 10000 + rep
            cfg = replace(
                config,
                death_rate=death,
                generations=100,
                max_cells=max(target * 4, config.max_cells),
            )
            result = simulate(cfg, seed=run_seed, therapy=False, target_size=target)
            sample = sample_at_size(
                result, target, np.random.default_rng(run_seed + 1_000_000)
            )
            if sample is not None:
                shannon[death].append(sample["shannon"])
                generations[death].append(result.history[-1]["generation"])
            # 未達門檻者仍寫入 CSV，保留達標與否及停止世代，避免隱藏失敗試驗。
            rows.append(
                dict(
                    experiment="death_rate",
                    condition=death,
                    replicate=rep,
                    seed=run_seed,
                    reached_target=sample is not None,
                    generation=result.history[-1]["generation"],
                    total=result.history[-1]["total"],
                    shannon=sample["shannon"] if sample else "",
                    suppressor_fraction=sample["suppressor_fraction"] if sample else "",
                )
            )
    return {
        "shannon": shannon,
        "generations": generations,
        "rows": rows,
        "target": target,
        "repeats": repeats,
    }


def run_il11_comparison(
    config: Config, repeats: int, seed: int, horizon: int = 18
) -> dict:
    """相同 20 個基因座下，IL11 功能開／關的生長、生產者比例與共享效益（不治療）。

    共享死亡率降低量與非生產者適應度都是生產者比例的確定性轉換，只寫進 CSV 與註腳。
    """
    if repeats < 1 or horizon < 1:
        raise ValueError("repeats 與 horizon 必須為正")
    config.validate()
    rows, summaries, series, producer = [], [], {}, None
    keys = (
        "total",
        "il11_percent",
        "shared_death_reduction",
        "nonproducer_expected_offspring",
    )
    for enabled, label in ((False, "IL11 OFF"), (True, "IL11 ON")):
        # 固定時間窗與其餘參數；同一 replicate 用同一 seed，但分歧後不視為逐細胞配對。
        cfg = replace(
            config, il11_enabled=enabled, generations=horizon, max_cells=10**12
        )
        trials = []
        for rep in range(repeats):
            run_seed = seed + 300000 + rep
            result = simulate(cfg, seed=run_seed, therapy=False)
            if result.stop_reason == "cell_limit":
                raise RuntimeError("IL11 comparison hit safety limit; shorten horizon")
            trial = []
            for generation in range(horizon + 1):
                state = result.history[min(generation, len(result.history) - 1)]
                # 非生產者本身沒有 IL11 突變；其 fitness 為基因型 0 的期望子代數。
                benefit = state["il11_shared_benefit"]
                death, triple = growth_probabilities(0, cfg, benefit)
                row = dict(
                    condition=label,
                    replicate=rep,
                    seed=run_seed,
                    generation=generation,
                    total=state["total"],
                    il11_percent=100 * state["il11_fraction"],
                    shared_death_reduction=benefit,
                    nonproducer_expected_offspring=(1 - death) * (2 + triple),
                )
                rows.append(row)
                trial.append(row)
            trials.append(trial)
        series[label] = {
            key: np.array([[r[key] for r in trial] for trial in trials])
            for key in keys[:2]
        }
        if enabled:
            producer = series[label]["il11_percent"]
        summaries.append(
            dict(
                condition=label,
                generation=horizon,
                replicates=repeats,
                **{
                    f"median_{key}": float(np.median([t[-1][key] for t in trials]))
                    for key in keys
                },
            )
        )
    return {
        "rows": rows,
        "summaries": summaries,
        "series": series,
        "producer": producer,
        "horizon": horizon,
        "repeats": repeats,
    }


def il11_benefit_rows(config: Config) -> list[dict]:
    """效益曲線的資料表，方便自行重畫或核對圖中數值。"""
    fraction, benefit = il11_benefit_curve(config)
    return [
        dict(il11_percent=100 * f, benefit_percentage_points=100 * b)
        for f, b in zip(fraction, benefit)
    ]


# ---------------------------------------------------------------- 數學圖與相圖資料 --
def write_markov_data(output: Path, markov: dict):
    """矩陣 P、Q、M 的 CSV（列＝目前／父代，欄＝下一步／子代）與假設說明。"""
    for name, values, names in (
        ("count_transition_P", markov["P"], ["0", "1", "2", "3", "overflow"]),
        ("mutation_Q", markov["Q"], ["S", "R"]),
        ("offspring_M_off", markov["M_off"], ["S", "R"]),
        ("offspring_M_on", markov["M_on"], ["S", "R"]),
    ):
        with (output / f"{name}.csv").open("w", newline="") as handle:
            writer = csv.writer(handle)
            writer.writerow(["current / next", *names])
            writer.writerows([[name, *row] for name, row in zip(names, values)])
    write_json(
        output / "matrix_assumptions.json",
        {
            "P": (
                "One-type chain, no mutation, no therapy, d=0.30, q=0.04; first"
                " N>=4 is absorbing bookkeeping overflow."
            ),
            "Q_and_M": (
                "Only resistance is retained; no oncogene/TSG advantage or IL11;"
                " mu=0.001."
            ),
            "orientation": (
                "Rows=current/parent, columns=next/offspring; probability"
                " distributions are row vectors."
            ),
            "verification": (
                "Exact convolution probabilities, P/Q row sums, and M expected row"
                " sums checked."
            ),
        },
        ensure_ascii=True,
    )


def phase_portrait_data(config: Config, trials: dict, events: dict, replicate: int):
    """重跑代表試驗並逐代核對 CSV，計算每一代由完整基因型算出的期望下一步。"""
    if replicate not in trials:
        raise ValueError(f"replicate {replicate} 不存在；可選 {sorted(trials)}")
    result = simulate(config, events[replicate]["seed"], therapy=True)
    assert len(result.history) == len(trials[replicate])
    for actual, reference in zip(result.history, trials[replicate]):
        assert all(actual[k] == reference[k] for k in reference), (actual, reference)
    outcome = therapy_outcome(result)
    assert outcome["nadir"] == events[replicate]["nadir"]
    assert outcome["relapse_generation"] == events[replicate]["relapse_generation"]
    assert result.therapy_generation == events[replicate]["therapy_generation"]
    steps = conditional_steps(result)
    # 沒有敏感父代且突變不可逆時，下一代敏感細胞數的期望必須為零。
    assert all(row["expected_sensitive"] == 0 for row in steps if row["sensitive"] == 0)
    points = event_generations(result.history, events[replicate])
    details = {
        label: {
            **result.history[g],
            "resistant_percent": (
                100 * result.history[g]["resistant"] / result.history[g]["total"]
                if result.history[g]["total"]
                else None
            ),
        }
        for label, g in points.items()
    }
    summary = {
        "source": str(OUTPUTS),
        "trial_count": len(trials),
        "observed_rows": sum(map(len, trials.values())),
        "highlight_replicate": replicate,
        "highlight_seed": result.seed,
        "selection": (
            "nearest median relapse generation; then nearest median nadir on log scale;"
            " replicate breaks ties"
        ),
        "stop_reason": result.stop_reason,
        "events": details,
        "coordinates": {
            "A": "(log10(S+1), log10(R+1))",
            "B": "(log10(T+1), 100 R/T), missing at T=0",
        },
        "conditional_arrows": (
            "raw-count conditional mean from full genotype state, transformed"
            " afterward; not E[log(count+1)]"
        ),
        "limitations": [
            "2D projection is not generally a closed Markov state",
            "no interpolated autonomous vector field",
            (
                "finite observed trajectories do not establish attractors or stable"
                " equilibria"
            ),
            "cell_limit is a computational stop, not carrying capacity",
        ],
        "verification": (
            "highlight reproduced against every observed CSV count; event timing and"
            " expectation conservation verified"
        ),
    }
    return result, steps, summary


def write_grid_csv(fields, destination):
    rows = []
    keys = (
        "sensitive",
        "resistant",
        "expected_sensitive",
        "expected_resistant",
        "delta_sensitive",
        "delta_resistant",
        "magnitude_cells",
        "x",
        "y",
        "display_u",
        "display_v",
        "display_magnitude",
    )
    for field in fields:
        scale = LINEAR_SCALE if field["coordinates"] == "linear" else LOG_SCALE
        for index in np.ndindex(field["sensitive"].shape):
            rows.append({
                "coordinates": field["coordinates"],
                "therapy_active": field["therapy_active"],
                **{key: field[key][index].item() for key in keys},
                "quiver_scale": scale,
                "drawn_u": field["display_u"][index] / scale,
                "drawn_v": field["display_v"][index] / scale,
            })
    with Path(destination).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return rows


def vector_field_record(
    config: Config, grids: dict, fields: list[dict], rows: list
) -> dict:
    """向量場資料的假設、截面與驗證說明。"""
    return {
        "source_metadata": str(OUTPUTS / "run_metadata.json"),
        "config": asdict(config),
        "slice": {
            "sensitive_genotype": 0,
            "resistant_genotype": 1 << config.resistance_gene,
            "other_bits": "all zero in parents; offspring mutations remain enabled",
            "therapy": (
                "externally fixed alpha=0 or alpha=1; no threshold switching inside"
                " panels"
            ),
        },
        "scope": (
            "Exact one-step projected conditional drift on a specified"
            " full-genotype-state slice; not an exact autonomous reduction of the"
            " evolving tumor"
        ),
        "matrix_off": fields[0]["matrix"].tolist(),
        "matrix_on": fields[1]["matrix"].tolist(),
        "grid": {
            key: {"counts": grid.tolist(), "points_per_panel": len(grid) ** 2}
            for key, grid in grids.items()
        },
        "arrow_length": {
            "linear": "Euclidean raw-count drift / 10000 / 16, equal axis units",
            "log": (
                "Euclidean phi(raw-count mean)-phi(current) displacement / 6, equal"
                " axis units"
            ),
            "normalization": (
                "one common scalar per figure; no per-vector normalization or panel"
                " autoscaling"
            ),
            "zero": "origin has zero drift, shown as a point instead of a finite arrow",
        },
        "row_count": len(rows),
        "verification": (
            "every grid point checked against projected_offspring_mean of its explicit"
            " full population; total conservation and shared scales checked"
        ),
        "interpretation": [
            "discrete-generation expected increments, not continuous-time derivatives",
            "logarithmic arrow lengths are not absolute-cell increments",
            "other mutations can move future full states off this slice",
        ],
    }


def phase_comparison_data(config: Config, results: list[Result]) -> dict:
    """相圖 05／06：兩組座標各自的理論箭頭與經驗平均箭頭。"""
    rows = trial_feature_rows(config, results)
    built = {pair: build_fields(rows, config, pair) for pair in ("logN", "fitness")}
    all_fields, record = [], {}
    for pair, (fields, steps, x_edges, y_edges) in built.items():
        all_fields += fields
        record[pair] = {
            "x_edges": x_edges.tolist(),
            "y_edges": y_edges.tolist(),
            "theory_arrows": sum(f["kind"] == "theory" for f in fields),
            "empirical_arrows": sum(f["kind"] == "empirical" for f in fields),
            "transitions": len(steps),
        }
    record.update({
        "y_coordinate": "empirical logit log10((R+0.5)/(S+0.5)); finite at rho=0 or 1",
        "empirical_rule": (
            f">= {MIN_TRANSITIONS} transitions from >= {MIN_TRIALS} trials per cell;"
            " arrow from mean start"
        ),
        "theory_slices": {
            "logN": "sensitive=genotype 0, resistant=resistance bit only",
            "fitness": (
                f"N={config.therapy_start_size}; oncogene fraction from w_off,"
                " independent of resistance; no TSG loss"
            ),
        },
        "arrow_scale": (
            "shared within each therapy column; length proportional to displacement;"
            " 95th percentile ≈ 0.9 cell"
        ),
        "base_model_only": True,
    })
    return {"rows": rows, "built": built, "fields": all_fields, "record": record}


# ---------------------------------------------------------------- 基因型資料 --------
def write_muller_csv(path: Path, result: Result, shown: list[int], config: Config):
    """代表試驗中每個已呈現基因型自己的比例（未呈現細胞已併入最近祖先）。"""
    _, _, _, own = muller_layout(result)
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["generation"] + [genotype_label(g, config) for g in shown])
        for t in range(len(own[0])):
            writer.writerow([t] + [own[g][t] for g in shown])


def write_gene_frequency_csv(path: Path, summary: dict):
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["series", "relative_generation", "median", "p10", "p90"])
        for key, s in summary.items():
            for i, x in enumerate(s["relative_generation"]):
                writer.writerow([key, x, s["median"][i], s["p10"][i], s["p90"][i]])


def genotype_record(
    replicate: int, result: Result, shown: list[int], config: Config
) -> dict:
    return {
        "representative_replicate": replicate,
        "seed": result.seed,
        "top_k": TOP_K,
        "shown_genotypes": [genotype_label(g, config) for g in shown],
        "csv_note": (
            "muller_representative.csv stores each shown genotype's own fraction after"
            " merging unshown cells into their nearest shown ancestor"
        ),
        "nesting": "genotype containment (bit subset), not reconstructed lineage",
        "base_model_only": True,
    }


# ---------------------------------------------------------------- 治療延伸分析 ------
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


def goldie_coldman_experiment(
    config: Config, trials: int, seed: int, workers: int
) -> dict:
    """08：治療起始負荷 100–5,000 細胞，每個大小 trials 次；與分枝過程理論比較。"""
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
    half = next((s["threshold"] for s in curve if s["p_cure"] < 0.5), None)
    q_on = curve[0]["q_on"]
    return {
        "rows": rows,
        "summary": summary,
        "grid": grid,
        "curve": curve,
        "half": half,
        "q_on": q_on,
        "record": {"summary": summary, "theory_half_cure_size": half, "q_on": q_on},
    }


def relapse_prediction(config: Config, results: list[Result]) -> dict:
    """09：用治療開始那一代的抗藥基因型組成預測每次試驗的復發延遲。"""
    w_off = fitness(0, config, False)
    m_r = fitness(1 << config.resistance_gene, config, True)
    m_s = fitness(0, config, True)
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
    stats["observed_sign_change_rho"] = [
        float(rho[growth < 0].max()),
        float(rho[growth > 0].min()),
    ]
    return {
        "rows": rows,
        "stats": stats,
        "obs": obs,
        "pred_baseline": pred1,
        "pred_genotype": pred2,
        "oncogene_share": onc,
        "rho": rho,
        "growth": growth,
        "oncogene_all": np.array([f["frac_onc"] for f in on]),
        "m_r": m_r,
        "m_s": m_s,
    }


def sensitivity_experiment(
    config: Config, trials: int, seed: int, workers: int
) -> dict:
    """10：抗藥細胞殺傷率 k_R × 突變率 μ，每格 trials 次；分類為治癒、復發或持續存活。"""
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
    return {
        "rows": rows,
        "kills": kills,
        "mutation_rates": mus,
        "cure": cure,
        "delay": delay,
        "persistent": persistent,
        "boundaries": boundaries,
        "record": {
            "kills": kills,
            "mutation_rates": mus,
            "p_cure": cure.tolist(),
            "median_delay": np.where(np.isfinite(delay), delay, None).tolist(),
            "persistent": persistent.tolist(),
            "boundaries": boundaries,
        },
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


def adaptive_experiment(config: Config, trials: int, seed: int, workers: int) -> dict:
    """11：間歇（適應性）與持續治療，在有／無承載量競爭與抗藥代價下的惡化時間。"""
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
        continuous = [
            r for r in rows if r["setting"] == name and r["strategy"] == "continuous"
        ]
        adaptive = [
            r for r in rows if r["setting"] == name and r["strategy"] == "adaptive"
        ]
        pairs = [
            (a["time_to_progression"], b["time_to_progression"])
            for a, b in zip(continuous, adaptive)
        ]
        entry["adaptive_longer_fraction"] = float(
            np.mean([
                (b if b is not None else 1e9) > (a if a is not None else 1e9)
                for a, b in pairs
            ])
        )
        summary.append(entry)
    return {
        "rows": rows,
        "summary": summary,
        "examples": examples,
        "settings": settings,
        "K": K,
        "record": {
            "summary": summary,
            "carrying_capacity": K,
            "crowding_at_K": crowding_death(K, replace(config, carrying_capacity=K)),
        },
    }
