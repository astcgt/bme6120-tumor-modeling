"""Extend presentation data with control, gene-frequency and therapy-analysis panels.

Regenerated trials are verified against the recorded CSVs.

Run once after prepare_results.py. Seed/config come from recorded metadata.
Core simulation is imported unchanged; saved CSV totals must match every trial.
"""

import csv
import json
from dataclasses import replace
from pathlib import Path
import sys
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from simulation import Config, simulate


def read(name):
    return list(csv.DictReader((Path("outputs") / name).open()))


def main():
    root = Path(__file__).resolve().parent
    destination = root / "data/results.json"
    data = json.loads(destination.read_text())
    meta = json.loads(Path("outputs/run_metadata.json").read_text())
    config_values = meta["config"].copy()
    for key in ("driver_genes", "suppressor_genes"):
        config_values[key] = tuple(config_values[key])
    config = Config(**config_values)
    trials = [simulate(config, seed) for seed in meta["therapy_seeds"]]
    raw = read("therapy_trajectories.csv")
    reference = {(int(r["replicate"]), int(r["generation"])): r for r in raw}
    for rep, trial in enumerate(trials):
        for generation, row in enumerate(trial.history):
            assert int(reference[rep, generation]["total"]) == row["total"]
            assert int(reference[rep, generation]["resistant"]) == row["resistant"]
    trajectory = read("trajectory.csv")
    for key in ("p10", "median", "p90"):
        data["therapy"]["effective_genotypes_" + key] = [
            float(r["effective_genotypes_" + key]) for r in trajectory[:31]
        ]
    hits, il11 = read("one_hit_vs_two_hit.csv"), read("il11_comparison.csv")
    data["endpoints"] = {}
    for rows, category, experiment in (
        (hits, "scenario", "hits"),
        (il11, "condition", "il11"),
    ):
        data["endpoints"][experiment] = {
            label: [
                float(r["total"])
                for r in rows
                if r[category] == label and int(r["generation"]) == 18
            ]
            for label in sorted({r[category] for r in rows})
        }
    # The omitted no-TSG-advantage control belongs to 02_scenario_comparisons.
    control_config = replace(
        config,
        generations=18,
        max_cells=10**12,
        suppressor_death_reduction=0,
        suppressor_triple_bonus=0,
    )
    control = [
        simulate(control_config, meta["seed"] + rep, therapy=False)
        for rep in range(meta["replicates"])
    ]
    reference_control = [
        r
        for r in read("comparisons.csv")
        if r["experiment"] == "growth" and r["condition"] == "No TSG advantage"
    ]
    assert len(reference_control) == meta["replicates"]
    for trial, row in zip(control, reference_control):
        assert trial.history[-1]["total"] == int(row["total"])
    values = np.array([
        [t.history[min(g, len(t.history) - 1)]["total"] for g in range(19)]
        for t in control
    ])
    data["no_tsg_advantage"] = {
        "generation": list(range(19)),
        **{
            f"total_{key}": np.percentile(values, p, axis=0).tolist()
            for key, p in [("p10", 10), ("median", 50), ("p90", 90)]
        },
    }
    data["il11_curve"] = {
        "fraction_percent": [i / 10 for i in range(101)],
        "benefit_percentage_points": [
            100
            * config.il11_max_death_reduction
            * (i / 1000)
            / (config.il11_half_fraction + i / 1000)
            for i in range(101)
        ],
    }
    # 取代原本的熱圖：以治療開始對齊的基因頻率中位數（main.py --only genotype 的輸出）。
    freq_rows = read("genotype_evolution/gene_frequency_summary.csv")
    series = {
        "ONC1": "0",
        "ONC2": "1",
        "ONC3": "2",
        "RES": "3",
        "TSG1 any hit": "4_any",
        "passenger mean": "passenger_mean",
    }
    data["gene_frequency"] = {
        "relative_generation": sorted(
            {int(r["relative_generation"]) for r in freq_rows}
        ),
        **{
            name: [
                100 * float(r["median"])
                for r in sorted(
                    (r for r in freq_rows if r["series"] == key),
                    key=lambda r: int(r["relative_generation"]),
                )
            ]
            for name, key in series.items()
        },
    }
    # 治療延伸分析（main.py --only therapy）：只讀取已核對的輸出，不重跑。
    therapy = json.loads(
        Path("outputs/therapy_analysis/therapy_analysis_summary.json").read_text()
    )
    data["therapy_analysis"] = {
        "goldie_coldman": [
            {
                k: r[k]
                for k in (
                    "threshold",
                    "n",
                    "cures",
                    "p_cure",
                    "ci_low",
                    "ci_high",
                    "theory_p_cure",
                    "R0_mean",
                    "R0_median",
                )
            }
            for r in therapy["goldie_coldman"]["summary"]
        ],
        "half_cure_size": therapy["goldie_coldman"]["theory_half_cure_size"],
        "relapse": therapy["relapse_prediction"],
        "relapse_rows": [
            {
                k: float(r[k])
                for k in (
                    "resistant_at_start",
                    "oncogene_fraction_of_resistant",
                    "observed_relapse_delay",
                    "predicted_delay_baseline",
                    "predicted_delay_genotype",
                )
            }
            for r in read("therapy_analysis/relapse_time_prediction.csv")
        ],
        "sensitivity": therapy["sensitivity"],
        "adaptive": therapy["adaptive"]["summary"],
    }
    data.pop("gene_heatmaps", None)
    data.pop("availability", None)
    data["verification"] = {
        "therapy_trials_reproduced": len(trials),
        "therapy_csv_rows_verified": len(reference),
        "control_trials_reproduced": len(control),
    }
    destination.write_text(json.dumps(data, ensure_ascii=False, indent=2) + "\n")
    print(json.dumps(data["verification"], indent=2))


if __name__ == "__main__":
    main()
