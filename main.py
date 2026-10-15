"""執行入口：依序跑模擬、存資料、畫圖。

    .venv/bin/python main.py                          # 全部步驟
    .venv/bin/python main.py --only results           # 只跑主實驗與對照實驗
    .venv/bin/python main.py --only therapy --analyses goldie_coldman relapse_prediction

步驟（後面的步驟讀取 results 寫出的 CSV，所以第一次要先跑 results）：
    results   主實驗（治療與復發）與對照實驗：腫瘤總覽、死亡率與異質性、one-hit vs two-hit
    model     數學圖：單細胞規則、矩陣與流向、附錄（轉移核與程式流程）
    phase     相圖：實際軌跡、負荷 × 抗藥比例、適應度 × 抗藥比例
    genotype  Muller 圖與基因頻率
    therapy   治療延伸分析：Goldie–Coldman、復發預測、參數敏感度、適應性治療

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
ANALYSES = ("goldie_coldman", "relapse_prediction", "sensitivity", "adaptive")


# 解析命令列參數，依指定步驟執行實驗、數學圖、相圖、基因型圖與治療分析。
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


# 解析執行步驟、種子、試驗數與分析選項，並檢查主實驗重複次數。
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
    parser.add_argument("--replicate", type=int, help="實際軌跡相圖要突出的試驗編號")
    parser.add_argument(
        "--analyses", choices=ANALYSES, nargs="+", help="只跑部分治療分析"
    )
    parser.add_argument(
        "--therapy-trials",
        type=int,
        default=200,
        help="goldie_coldman 每個大小的試驗數",
    )
    parser.add_argument(
        "--grid-trials",
        type=int,
        default=30,
        help="sensitivity 每格、adaptive 每組試驗數",
    )
    parser.add_argument("--workers", type=int, default=os.cpu_count())
    args = parser.parse_args()
    if args.replicates < 1:
        parser.error("--replicates must be >= 1")
    return args


# 依輸出紀錄重現完整治療試驗並快取，供相圖、基因型與復發分析共用。
@lru_cache(maxsize=1)
def reproduced_trials():
    """相圖、基因型與復發分析共用：依記錄的 seed 重跑治療試驗，取得完整基因型快照。"""
    return exp.reproduce_trials(exp.OUTPUTS)


# 執行治療、TSG 與死亡率實驗，輸出比較圖、CSV 及可重現的執行紀錄。
def run_results(args):
    """主實驗與兩組對照實驗：腫瘤總覽、死亡率與異質性、one-hit vs two-hit 三張圖與對應的 CSV。"""
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


# 計算模型矩陣與線性向量場，輸出數學圖、流程附錄及矩陣資料。
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


# 重現代表軌跡並建立相圖、向量場及雙語理論與經驗比較圖，輸出相關資料。
def run_phase(args):
    """實際軌跡相圖（與期望下一步）、向量場資料、兩張理論 vs 經驗相場圖。"""
    out = exp.OUTPUTS / "phase_plots"
    out.mkdir(parents=True, exist_ok=True)
    config, trials, events = exp.load_trials(exp.OUTPUTS)

    # 實際軌跡相圖：重跑代表試驗，計算每一代的期望下一步。
    rep = analysis.choose_exemplar(events) if args.replicate is None else args.replicate
    result, steps, summary = exp.phase_portrait_data(config, trials, events, rep)
    exp.write_csv(out / "conditional_phase_steps.csv", steps)
    exp.write_json(out / "phase_summary.json", summary)
    plot_style.set_phase_style()
    figure = plot_phase.plot_phase_portrait(trials, events, rep, config, result, steps)
    plot_phase.save_phase_figure(figure, out, "phase_portrait")

    # 向量場資料（矩陣與流向圖使用同一組公式）。
    grids, fields = analysis.vector_field_grids(config)
    rows = exp.write_grid_csv(fields, out / "phase_vector_grid.csv")
    record = exp.vector_field_record(config, grids, fields, rows)
    exp.write_json(out / "phase_vector_assumptions.json", record)

    # 理論 vs 經驗相場：兩組座標的理論箭頭與經驗平均箭頭。
    base_config, results, _ = reproduced_trials()
    analysis.check_expectation(base_config)
    comparison = exp.phase_comparison_data(base_config, results)
    exp.write_csv(out / "phase_field_comparison.csv", comparison["fields"])
    exp.write_json(out / "phase_field_comparison.json", comparison["record"])
    names = {
        "logN": "phase_field_burden_logit",
        "fitness": "phase_field_fitness_logit",
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


# 重現治療試驗，輸出雙語 Muller 圖、基因頻率圖及對應資料。
def run_genotype():
    """Muller 圖與每基因突變頻率（雙語），以及繪圖數值的 CSV。"""
    out = exp.OUTPUTS / "genotype_evolution"
    out.mkdir(parents=True, exist_ok=True)
    config, results, events = reproduced_trials()
    rel, data, summary = analysis.gene_frequency_series(config, results)
    for lang in plot_style.LANGS:
        plot_style.set_style(lang)
        figure, rep, shown = plot_genotype.plot_muller(config, results, events, lang)
        plot_style.save_lang(figure, out, "muller_top_genotypes", lang)
        figure = plot_genotype.plot_gene_frequencies(
            rel, data, config, len(results), lang
        )
        plot_style.save_lang(figure, out, "gene_mutation_frequencies", lang)
    exp.write_muller_csv(out / "muller_representative.csv", results[rep], shown, config)
    exp.write_gene_frequency_csv(out / "gene_frequency_summary.csv", summary)
    record = exp.genotype_record(rep, results[rep], shown, config)
    exp.write_json(out / "genotype_evolution.json", record)
    print(f"Genotype figures saved to {out}")


# 執行選定的治療延伸分析，輸出圖表與 CSV 並更新共用摘要。
def run_therapy(args):
    """治療延伸分析（ANALYSES）；摘要依分析名稱累積寫入 therapy_analysis_summary.json。"""
    out = exp.OUTPUTS / "therapy_analysis"
    out.mkdir(parents=True, exist_ok=True)
    config = Config()
    summary_path = out / "therapy_analysis_summary.json"
    summary = json.loads(summary_path.read_text()) if summary_path.exists() else {}
    todo = args.analyses or ANALYSES
    if "goldie_coldman" in todo:
        data = exp.goldie_coldman_experiment(
            config, args.therapy_trials, args.seed, args.workers
        )
        plot_therapy.plot_goldie_coldman(data, config, args.therapy_trials, out)
        exp.write_csv(out / "goldie_coldman_trials.csv", data["rows"])
        exp.write_csv(out / "goldie_coldman_summary.csv", data["summary"])
        summary["goldie_coldman"] = data["record"]
    if "relapse_prediction" in todo:
        base_config, results, _ = reproduced_trials()
        data = exp.relapse_prediction(base_config, results)
        plot_therapy.plot_relapse_prediction(data, out)
        exp.write_csv(out / "relapse_time_prediction.csv", data["rows"])
        summary["relapse_prediction"] = data["stats"]
    if "sensitivity" in todo:
        data = exp.sensitivity_experiment(
            config, args.grid_trials, args.seed, args.workers
        )
        plot_therapy.plot_sensitivity(data, args.grid_trials, out)
        exp.write_csv(out / "parameter_sensitivity_trials.csv", data["rows"])
        summary["sensitivity"] = data["record"]
    if "adaptive" in todo:
        data = exp.adaptive_experiment(
            config, args.grid_trials, args.seed, args.workers
        )
        plot_therapy.plot_adaptive(data, config, args.grid_trials, out)
        exp.write_csv(out / "adaptive_therapy_trials.csv", data["rows"])
        exp.write_csv(out / "adaptive_therapy_summary.csv", data["summary"])
        summary["adaptive"] = data["record"]
    summary["config"] = asdict(config)
    exp.write_json(summary_path, summary, ensure_ascii=True, default=float)
    print(f"Therapy analyses saved to {out}")


if __name__ == "__main__":
    main()
