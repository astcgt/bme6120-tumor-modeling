"""規則網格相圖：箭頭方向與長度呈現一步條件期望漂移。

執行：.venv/bin/python plot_phase_field.py
固定基準基因型截面，分別固定治療關閉／持續治療；不修改 PPT。
"""
from __future__ import annotations

import argparse
import csv
from dataclasses import asdict
import json
from pathlib import Path

from modeling_cancer import Config
from plot_phase_portraits import projected_offspring_mean, set_phase_style, save_phase_figure
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
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


def drift_formula(matrix):
    a, b, c = matrix[0, 0]-1, matrix[0, 1], matrix[1, 1]-1
    return rf"$\Delta S={a:.7g}S,\quad\Delta R={b:.7g}S+{c:.7g}R$"


def plot_phase_field(fields, config):
    """兩面板共用線性箭頭倍率及色階，保留大小差異。"""
    coordinates = fields[0]["coordinates"]
    linear = coordinates == "linear"
    scale = LINEAR_SCALE if linear else LOG_SCALE
    fig = plt.figure(figsize=(17.2, 10.3))
    axes = [fig.add_axes((.075, .28, .355, .51)), fig.add_axes((.535, .28, .355, .51))]
    color_ax = fig.add_axes((.925, .31, .012, .445))
    fig.text(.055, .966, "每個網格點都有方向，箭頭長短保留變化大小", fontsize=25, weight="bold", va="top")
    fig.text(.055, .907, "固定基準基因型截面：S 的所有位元為 0；R 僅抗藥位元為 1。治療控制在各面板固定。", fontsize=14, color=GRAY)
    magnitude_key = "magnitude_cells" if linear else "display_magnitude"
    norm = Normalize(0, max(np.max(field[magnitude_key]) for field in fields))
    quivers = []
    for i, (ax, field) in enumerate(zip(axes, fields)):
        quiver = ax.quiver(field["x"], field["y"], field["display_u"], field["display_v"], field[magnitude_key],
                           angles="xy", scale_units="xy", scale=scale, pivot="tail", cmap="viridis", norm=norm,
                           width=.0031, headwidth=3.4, headlength=4.4, headaxislength=3.8,
                           minlength=0, minshaft=1)
        quivers.append(quiver)
        # 原點的期望漂移為零；以點表示，不能偽造一支有長度的箭頭。
        ax.scatter(0, 0, s=22, color=INK, zorder=4)
        ax.set_aspect("equal", adjustable="box")
        ax.grid(color="#E1E7EB", lw=.65, alpha=.75)
        ax.set_axisbelow(True)
        matrix = field["matrix"]
        m_s, m_r = np.sum(matrix[0]), np.sum(matrix[1])
        title = (r"A  治療關閉 ($\alpha=0$)" if i == 0 else r"B  持續治療 ($\alpha=1$)")
        title += f"\n每父細胞期望子代：S={m_s:.5g}，R={m_r:.5g}"
        ax.set_title(title, loc="left", fontsize=15.5, weight="bold", pad=15)
        if linear:
            ax.set(xlim=(-.065, 2.70), ylim=(-.065, 2.70),
                   xlabel=r"敏感細胞 $S$（$\times10^4$）", ylabel=r"抗藥細胞 $R$（$\times10^4$）")
            ax.set_xticks(np.arange(0, 2.6, .5)); ax.set_yticks(np.arange(0, 2.6, .5))
        else:
            ax.set(xlim=(-.12, 5.74), ylim=(-.12, 5.74),
                   xlabel=r"敏感細胞 $\log_{10}(S+1)$", ylabel=r"抗藥細胞 $\log_{10}(R+1)$")
            ax.set_xticks(range(6)); ax.set_yticks(range(6))
        ax.tick_params(labelsize=12)
        ax.xaxis.labelpad = 10; ax.yaxis.labelpad = 10
    bar = fig.colorbar(quivers[0], cax=color_ax)
    bar.ax.tick_params(labelsize=11)
    if linear:
        bar.set_label(r"$\|\Delta\|_2$（細胞／代）", labelpad=12, fontsize=12)
        bar.ax.ticklabel_format(style="sci", scilimits=(4,4))
        for x, field in zip((.075, .535), fields):
            fig.text(x, .207, drift_formula(field["matrix"]), fontsize=16, color=BLUE)
        axes[0].quiverkey(quivers[0], X=.089, Y=.143, U=1, label="10,000 細胞／代", labelpos="E",
                          coordinates="figure", fontproperties={"size": 12}, color=BLUE)
        fig.text(.34, .137, r"箭頭長度 $\propto\sqrt{(\Delta S)^2+(\Delta R)^2}$；兩圖共用倍率，未做單支箭頭正規化。", fontsize=12.5)
        fig.text(.075, .087, r"$\Delta=\mathrm{E}[(S_{t+1},R_{t+1})\mid s,\alpha]-(S,R)$。離散世代的一步期望變化；原點零變化以圓點標記。", fontsize=12.5, color=GRAY)
    else:
        bar.set_label("對數座標位移大小／代", labelpad=12, fontsize=12)
        fig.text(.075, .207, r"$\mathbf{D}=\phi(\mathrm{E}[(S_{t+1},R_{t+1})\mid s,\alpha])-\phi(S,R),\quad\phi(S,R)=(\log_{10}(S+1),\log_{10}(R+1))$", fontsize=16, color=BLUE)
        axes[0].quiverkey(quivers[0], X=.089, Y=.143, U=.5, label="0.5 座標單位／代", labelpos="E",
                          coordinates="figure", fontproperties={"size": 12}, color=BLUE)
        fig.text(.34, .137, r"長度 $\propto\|\mathbf{D}\|_2$，色階亦表示對數位移；兩圖共用倍率，原點為零變化。", fontsize=12.5)
        fig.text(.075, .087, "先對原始計數求期望，再取對數；長度代表對數座標的位移，絕對細胞數變化請看線性版。", fontsize=12.5, color=GRAY)
    fig.text(.075, .04, f"μ={config.mutation_rate:g}；driver／TSG 未達標、IL11 無生產者。子代仍可突變；後續離開此截面時，需回到完整基因型模型。", fontsize=12, color=GRAY)
    # 不依面板單獨自動縮放，也不將向量除以自己的長度。
    assert quivers[0].scale == quivers[1].scale == scale
    return fig


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
    set_phase_style()
    for coordinates, name in (("linear", "03_phase_vector_field"), ("log", "04_phase_vector_field_log")):
        fig = plot_phase_field([f for f in fields if f["coordinates"] == coordinates], config)
        save_phase_figure(fig, args.output, name)
    print(json.dumps({"verified_grid_points": len(rows), "nonzero_arrows": len(rows)-4,
                      "output": str(args.output)}, ensure_ascii=False))


if __name__ == "__main__":
    main()
