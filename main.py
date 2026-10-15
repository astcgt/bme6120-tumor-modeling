"""執行入口：依序跑模擬、存資料、畫圖。

    .venv/bin/python main.py                          # 全部步驟
    .venv/bin/python main.py --only results           # 只跑主實驗與對照實驗
    .venv/bin/python main.py --only therapy --analyses 08 09

步驟（後面的步驟讀取 results 寫出的 CSV，所以第一次要先跑 results）：
    results   主實驗（治療與復發）與對照實驗：圖 01–03
    model     數學圖：06 單細胞規則、07 矩陣與流向、A1 附錄
    phase     相圖：01、05、06
    genotype  Muller 圖與基因頻率
    therapy   治療延伸分析：08–11

檔案分工：simulation.py 模型本身；experiments.py 重複試驗與資料讀寫；
analysis.py 理論公式；plot_*.py 畫圖。
"""

from __future__ import annotations

import argparse
import json
import os
from dataclasses import asdict
from functools import lru_cache

import plot_style  # 先設定 matplotlib 後端與快取位置
import analysis
import experiments as exp
import plot_genotype
import plot_model
import plot_phase
import plot_results
import plot_therapy
from simulation import Config

STEPS = ("results", "model", "phase", "genotype", "therapy")
ANALYSES = ("08", "09", "10", "11")


def main():
    args = parse_args()
    steps = args.only or STEPS
    if "results" in steps:
        run_results(args)
    if "model" in steps:
        run_model()
    if "phase" in steps:
        run_phase(args)
    if "genotype" in steps:
        run_genotype()
    if "therapy" in steps:
        run_therapy(args)


def parse_args():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--only", choices=STEPS, nargs="+", help="只跑指定步驟")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--replicates", type=int, default=30, help="每組重複試驗數")
    parser.add_argument(
        "--suppressor-hits",
        type=int,
        choices=(1, 2),
        default=2,
        help="主實驗的 TSG 失活門檻（one-hit 或 two-hit）",
    )
    parser.add_argument("--replicate", type=int, help="相圖 01 要突出的試驗編號")
    parser.add_argument(
        "--analyses", choices=ANALYSES, nargs="+", help="只跑部分治療分析"
    )
    parser.add_argument(
        "--therapy-trials", type=int, default=200, help="08 每個大小的試驗數"
    )
    parser.add_argument(
        "--grid-trials", type=int, default=30, help="10 每格、11 每組試驗數"
    )
    parser.add_argument("--workers", type=int, default=os.cpu_count())
    args = parser.parse_args()
    if args.replicates < 1:
        parser.error("--replicates must be >= 1")
    return args


@lru_cache(maxsize=1)
def reproduced_trials():
    """相圖、基因型與復發分析共用：依記錄的 seed 重跑治療試驗，取得完整基因型快照。"""
    return exp.reproduce_trials(exp.OUTPUTS)


def run_results(args):
    """主實驗與兩組對照實驗：圖 01–03 與對應的 CSV。"""
    out = exp.OUTPUTS
    out.mkdir(parents=True, exist_ok=True)
    config = Config(suppressor_hits=args.suppressor_hits)

    # 主實驗：30 次治療試驗。
    trials = exp.run_therapy_trials(config, args.replicates, args.seed)
    summary, raw = exp.therapy_ensemble(trials)
    plot_results.plot_therapy_ensemble(trials, summary, out)
    exp.write_csv(out / "trajectory.csv", summary)
    exp.write_csv(out / "therapy_trajectories.csv", raw)

    # 對照實驗：死亡率與異質性、TSG one-hit vs two-hit。
    tsg_rows, growth_rows = exp.run_tsg_comparison(config, args.replicates, args.seed)
    death = exp.run_death_rate_comparison(config, args.replicates, args.seed)
    plot_results.plot_death_rate(
        death["shannon"],
        death["generations"],
        config,
        args.replicates,
        death["target"],
        args.seed,
        out,
    )
    exp.write_csv(out / "comparisons.csv", growth_rows + death["rows"])
    hits = exp.hit_comparison_series(tsg_rows)
    plot_results.plot_hit_comparison(hits, config, out)
    exp.write_csv(out / "one_hit_vs_two_hit.csv", tsg_rows)
    exp.write_csv(out / "one_hit_vs_two_hit_summary.csv", hits["summaries"])

    # 每次試驗的結果與本次執行的完整設定。
    outcomes = exp.therapy_events(trials)
    exp.write_csv(out / "therapy_replicates.csv", outcomes)
    metadata = exp.run_metadata(config, args.seed, args.replicates, trials, outcomes)
    exp.write_json(
        out / "run_metadata.json", metadata, ensure_ascii=True, trailing_newline=False
    )
    print(json.dumps(metadata, indent=2))
    print(f"Figures and data saved to: {out.resolve()}")


def run_model():
    """數學圖：矩陣與向量場由 analysis.py 計算，實際路徑讀取主實驗的 CSV。"""
    out = plot_model.OUT
    out.mkdir(parents=True, exist_ok=True)
    config = Config()
    markov = analysis.markov_matrices(config)
    fields = analysis.linear_flow_fields(config)
    _, trials, events = exp.load_trials(exp.OUTPUTS)
    plot_style.set_model_style()
    plot_model.single_cell_figure(
        config, markov["P"], markov["law"], markov["phenotypes"]
    )
    plot_model.matrices_to_flow_figure(
        config,
        markov["Q"],
        markov["M_off"],
        markov["M_on"],
        fields,
        trials,
        events,
        analysis.choose_exemplar(events),
    )
    plot_model.appendix_figure(config)
    exp.write_markov_data(out, markov)
    print(f"Verified matrices and saved PNG/SVG figures to {out}")


def run_phase(args):
    """相圖 01（實際軌跡與期望下一步）、向量場資料、相圖 05／06（理論 vs 經驗）。"""
    out = exp.OUTPUTS / "phase_plots"
    out.mkdir(parents=True, exist_ok=True)
    config, trials, events = exp.load_trials(exp.OUTPUTS)

    # 相圖 01：重跑代表試驗，計算每一代的期望下一步。
    rep = analysis.choose_exemplar(events) if args.replicate is None else args.replicate
    result, steps, summary = exp.phase_portrait_data(config, trials, events, rep)
    exp.write_csv(out / "conditional_phase_steps.csv", steps)
    exp.write_json(out / "phase_summary.json", summary)
    plot_style.set_phase_style()
    figure = plot_phase.plot_phase_portrait(trials, events, rep, config, result, steps)
    plot_phase.save_phase_figure(figure, out, "01_phase_portrait")

    # 向量場資料（圖 07 使用同一組公式）。
    grids, fields = analysis.vector_field_grids(config)
    rows = exp.write_grid_csv(fields, out / "phase_vector_grid.csv")
    record = exp.vector_field_record(config, grids, fields, rows)
    exp.write_json(out / "phase_vector_assumptions.json", record)

    # 相圖 05／06：兩組座標的理論箭頭與經驗平均箭頭。
    base_config, results, _ = reproduced_trials()
    analysis.check_expectation(base_config)
    comparison = exp.phase_comparison_data(base_config, results)
    exp.write_csv(out / "phase_field_comparison.csv", comparison["fields"])
    exp.write_json(out / "phase_field_comparison.json", comparison["record"])
    names = {
        "logN": "05_phase_field_burden_logit",
        "fitness": "06_phase_field_fitness_logit",
    }
    for lang in plot_style.LANGS:
        plot_style.set_style(lang)
        for pair, (fields, steps, x_edges, y_edges) in comparison["built"].items():
            figure = plot_phase.plot_pair(
                pair,
                comparison["rows"],
                fields,
                steps,
                x_edges,
                y_edges,
                base_config,
                lang,
            )
            plot_style.save_lang(figure, out, names[pair], lang)
    print(f"Phase portraits saved to {out}")


def run_genotype():
    """Muller 圖與每基因突變頻率（雙語），以及繪圖數值的 CSV。"""
    out = exp.OUTPUTS / "genotype_evolution"
    out.mkdir(parents=True, exist_ok=True)
    config, results, events = reproduced_trials()
    rel, data, summary = analysis.gene_frequency_series(config, results)
    for lang in plot_style.LANGS:
        plot_style.set_style(lang)
        figure, rep, shown = plot_genotype.plot_muller(config, results, events, lang)
        plot_style.save_lang(figure, out, "01_muller_top_genotypes", lang)
        figure = plot_genotype.plot_gene_frequencies(
            rel, data, config, len(results), lang
        )
        plot_style.save_lang(figure, out, "02_gene_mutation_frequencies", lang)
    exp.write_muller_csv(out / "muller_representative.csv", results[rep], shown, config)
    exp.write_gene_frequency_csv(out / "gene_frequency_summary.csv", summary)
    record = exp.genotype_record(rep, results[rep], shown, config)
    exp.write_json(out / "genotype_evolution.json", record)
    print(f"Genotype figures saved to {out}")


def run_therapy(args):
    """治療延伸分析 08–11；摘要累積寫入 therapy_analysis_summary.json。"""
    out = exp.OUTPUTS / "therapy_analysis"
    out.mkdir(parents=True, exist_ok=True)
    config = Config()
    summary_path = out / "therapy_analysis_summary.json"
    summary = json.loads(summary_path.read_text()) if summary_path.exists() else {}
    todo = args.analyses or ANALYSES
    if "08" in todo:
        data = exp.goldie_coldman_experiment(
            config, args.therapy_trials, args.seed, args.workers
        )
        plot_therapy.plot_goldie_coldman(data, config, args.therapy_trials, out)
        exp.write_csv(out / "goldie_coldman_trials.csv", data["rows"])
        exp.write_csv(out / "goldie_coldman_summary.csv", data["summary"])
        summary["08_goldie_coldman"] = data["record"]
    if "09" in todo:
        base_config, results, _ = reproduced_trials()
        data = exp.relapse_prediction(base_config, results)
        plot_therapy.plot_relapse_prediction(data, out)
        exp.write_csv(out / "relapse_time_prediction.csv", data["rows"])
        summary["09_relapse_time"] = data["stats"]
    if "10" in todo:
        data = exp.sensitivity_experiment(
            config, args.grid_trials, args.seed, args.workers
        )
        plot_therapy.plot_sensitivity(data, args.grid_trials, out)
        exp.write_csv(out / "parameter_sensitivity_trials.csv", data["rows"])
        summary["10_parameter_sensitivity"] = data["record"]
    if "11" in todo:
        data = exp.adaptive_experiment(
            config, args.grid_trials, args.seed, args.workers
        )
        plot_therapy.plot_adaptive(data, config, args.grid_trials, out)
        exp.write_csv(out / "adaptive_therapy_trials.csv", data["rows"])
        exp.write_csv(out / "adaptive_therapy_summary.csv", data["summary"])
        summary["11_adaptive_therapy"] = data["record"]
    summary["config"] = asdict(config)
    exp.write_json(summary_path, summary, ensure_ascii=True, default=float)
    print(f"Therapy analyses saved to {out}")


if __name__ == "__main__":
    main()
