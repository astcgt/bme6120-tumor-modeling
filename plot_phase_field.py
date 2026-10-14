"""規則網格向量場的資料：每個網格點的一步條件期望漂移。

執行：.venv/bin/python plot_phase_field.py   → phase_vector_grid.csv、phase_vector_assumptions.json
固定基準基因型截面，分別固定治療關閉／持續治療。圖畫在 plot_markov_model.py 的圖 07
（M 矩陣旁的線性向量場）；舊的 03 線性、04 對數獨立圖已併入或刪除。
"""
from __future__ import annotations

import argparse
import csv
from dataclasses import asdict
import json
from pathlib import Path

from modeling_cancer import Config
from plot_phase_portraits import projected_offspring_mean
import numpy as np


ROOT = Path(__file__).resolve().parent
INK, BLUE, GRAY = "#173042", "#236B8E", "#687783"
CELL_UNIT = 10_000
LINEAR_SCALE = 16.0
LOG_SCALE = 6.0


def baseline_mean_matrix(config: Config, therapy_active: bool) -> np.ndarray:
    """父代只取 g=0 和 g=抗藥位元；用完整模型的一步期望取得 M。

    所有其他位元為 0。其他突變仍可出現在子代，因此這是完整狀態
    空間的截面，而不是聲稱整段基因型演化保持在兩個基因型中。
    """
    config.validate()
    resistance = 1 << config.resistance_gene
    sensitive = projected_offspring_mean({0: 1}, config, therapy_active)
    resistant = projected_offspring_mean({resistance: 1}, config, therapy_active)
    matrix = np.array([[sensitive["expected_sensitive"], sensitive["expected_resistant"]],
                       [resistant["expected_sensitive"], resistant["expected_resistant"]]])
    assert matrix[1, 0] == 0  # 抗藥突變不可逆。
    return matrix


def compute_phase_field(config: Config, therapy_active: bool, counts,
                        coordinates="linear") -> dict:
    """每個整數網格點的 E[(S',R') | 完整基準狀態, α] − (S,R)。

    display_u/v 是座標內的一步位移，不進行逐向量正規化。
    對數版先對原始計數求期望，再轉換；不是 E[log(count+1)]。
    """
    counts = np.asarray(counts)
    if counts.ndim != 1 or len(counts) < 2 or np.any(counts < 0) or np.any(counts != np.rint(counts)):
        raise ValueError("counts 必須為至少兩個非負整數的一維陣列")
    if np.any(np.diff(counts) <= 0):
        raise ValueError("counts 必須嚴格遞增")
    if coordinates not in ("linear", "log"):
        raise ValueError("coordinates 必須是 linear 或 log")
    sensitive, resistant = np.meshgrid(counts.astype(np.int64), counts.astype(np.int64))
    matrix = baseline_mean_matrix(config, therapy_active)
    current = np.stack((sensitive, resistant), axis=-1)
    mean = current @ matrix
    delta = mean - current
    if coordinates == "linear":
        position = current / CELL_UNIT
        displacement = delta / CELL_UNIT
    else:
        position = np.log10(current + 1)
        displacement = np.log10(mean + 1) - position
    return {"sensitive": sensitive, "resistant": resistant,
            "expected_sensitive": mean[..., 0], "expected_resistant": mean[..., 1],
            "delta_sensitive": delta[..., 0], "delta_resistant": delta[..., 1],
            "magnitude_cells": np.linalg.norm(delta, axis=-1),
            "x": position[..., 0], "y": position[..., 1],
            "display_u": displacement[..., 0], "display_v": displacement[..., 1],
            "display_magnitude": np.linalg.norm(displacement, axis=-1),
            "matrix": matrix, "coordinates": coordinates, "therapy_active": therapy_active}


def verify_field(field, config):
    """逐點核對完整狀態期望；另驗證總數及零向量的語意。"""
    resistance = 1 << config.resistance_gene
    for index in np.ndindex(field["sensitive"].shape):
        s, r = int(field["sensitive"][index]), int(field["resistant"][index])
        population = {g: n for g, n in ((0, s), (resistance, r)) if n}
        mean = projected_offspring_mean(population, config, field["therapy_active"])
        np.testing.assert_allclose([field["expected_sensitive"][index], field["expected_resistant"][index]],
                                   [mean["expected_sensitive"], mean["expected_resistant"]], rtol=1e-12, atol=1e-10)
        assert np.isclose(field["expected_sensitive"][index]+field["expected_resistant"][index],
                          mean["expected_total"], rtol=1e-12)
    assert field["display_magnitude"][0, 0] == 0
    assert field["magnitude_cells"][0, 0] == 0
    assert np.all(field["expected_sensitive"][field["sensitive"] == 0] == 0)


def write_grid_csv(fields, destination):
    rows = []
    keys = ("sensitive", "resistant", "expected_sensitive", "expected_resistant",
            "delta_sensitive", "delta_resistant", "magnitude_cells", "x", "y",
            "display_u", "display_v", "display_magnitude")
    for field in fields:
        scale = LINEAR_SCALE if field["coordinates"] == "linear" else LOG_SCALE
        for index in np.ndindex(field["sensitive"].shape):
            rows.append({"coordinates": field["coordinates"], "therapy_active": field["therapy_active"],
                         **{key: field[key][index].item() for key in keys},
                         "quiver_scale": scale, "drawn_u": field["display_u"][index]/scale,
                         "drawn_v": field["display_v"][index]/scale})
    with Path(destination).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=ROOT / "outputs")
    parser.add_argument("--output", type=Path, default=ROOT / "outputs/phase_plots")
    args = parser.parse_args()
    meta = json.loads((args.input / "run_metadata.json").read_text())
    settings = meta["config"].copy()
    for key in ("driver_genes", "suppressor_genes"):
        settings[key] = tuple(settings[key])
    config = Config(**settings)
    config.validate()
    grids = {"linear": np.rint(np.linspace(0, 25_000, 19)).astype(np.int64),
             "log": np.rint(np.power(10, np.linspace(0, 5.5, 21))-1).astype(np.int64)}
    fields = [compute_phase_field(config, active, counts, coordinates)
              for coordinates, counts in grids.items() for active in (False, True)]
    # 目前範圍均在計算上限內，避免對已停止執行的狀態畫出更新箭頭。
    assert max(2*grid[-1] for grid in grids.values()) < config.max_cells
    for field in fields:
        verify_field(field, config)
    args.output.mkdir(parents=True, exist_ok=True)
    rows = write_grid_csv(fields, args.output / "phase_vector_grid.csv")
    record = {"source_metadata": str(args.input / "run_metadata.json"), "config": asdict(config),
              "slice": {"sensitive_genotype": 0, "resistant_genotype": 1 << config.resistance_gene,
                        "other_bits": "all zero in parents; offspring mutations remain enabled",
                        "therapy": "externally fixed alpha=0 or alpha=1; no threshold switching inside panels"},
              "scope": "Exact one-step projected conditional drift on a specified full-genotype-state slice; not an exact autonomous reduction of the evolving tumor",
              "matrix_off": fields[0]["matrix"].tolist(), "matrix_on": fields[1]["matrix"].tolist(),
              "grid": {key: {"counts": grid.tolist(), "points_per_panel": len(grid)**2} for key, grid in grids.items()},
              "arrow_length": {"linear": "Euclidean raw-count drift / 10000 / 16, equal axis units",
                               "log": "Euclidean phi(raw-count mean)-phi(current) displacement / 6, equal axis units",
                               "normalization": "one common scalar per figure; no per-vector normalization or panel autoscaling",
                               "zero": "origin has zero drift, shown as a point instead of a finite arrow"},
              "row_count": len(rows), "verification": "every grid point checked against projected_offspring_mean of its explicit full population; total conservation and shared scales checked",
              "interpretation": ["discrete-generation expected increments, not continuous-time derivatives",
                                 "logarithmic arrow lengths are not absolute-cell increments",
                                 "other mutations can move future full states off this slice"]}
    (args.output / "phase_vector_assumptions.json").write_text(json.dumps(record, ensure_ascii=False, indent=2)+"\n")
    print(json.dumps({"verified_grid_points": len(rows), "nonzero_arrows": len(rows)-4,
                      "output": str(args.output)}, ensure_ascii=False))


if __name__ == "__main__":
    main()
