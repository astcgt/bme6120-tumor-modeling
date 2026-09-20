"""BME 6120：離散分枝腫瘤模型及展示圖表。

執行：.venv/bin/python modeling_cancer.py
模型以基因型計數（不是追蹤家系），所有參數均為教學假設，不是臨床估計。
每代每個細胞死亡或產生 2 / 3 個子細胞；子細胞獨立發生不可逆突變。
抑癌基因以兩個等位基因建模；one-hit／two-hit 指同一抑癌基因的失活門檻。
"""

# 延後解析型別註記，讓函數與資料類別的型別宣告更容易互相引用。
from __future__ import annotations

# 標準函式庫：解析命令列、讀寫 CSV／JSON、處理路徑與建立參數資料類別。
import argparse
import csv
import json
import os
from dataclasses import asdict, dataclass, replace
from pathlib import Path

# 圖形快取留在專案中，不寫入 conda 或使用者的 matplotlib 設定目錄。
os.environ.setdefault("MPLCONFIGDIR", str(Path(__file__).resolve().parent / ".mplconfig"))
import matplotlib
# 使用不需要視窗的繪圖後端，適合直接輸出簡報圖片。
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
# NumPy 負責隨機抽樣、陣列運算與統計摘要。
import numpy as np


# 將所有模型參數集中管理；frozen=True 防止執行中意外修改設定。
@dataclass(frozen=True)
class Config:
    # 20 個基因座：一般基因用一個位元；每個抑癌基因另加一個等位基因位元。
    n_genes: int = 20
    driver_genes: tuple[int, ...] = (0, 1, 2)
    resistance_gene: int = 3
    driver_hits: int = 1  # 保留：不同 oncogene 的突變數門檻
    suppressor_genes: tuple[int, ...] = (4, 5)
    suppressor_hits: int = 2  # 同一抑癌基因需要失活的等位基因數
    suppressor_death_reduction: float = 0.08
    suppressor_triple_bonus: float = 0.12
    # 突變與生長設定：driver 達標後同時影響自然死亡率和三子細胞機率。
    mutation_rate: float = 0.001  # 每個子細胞、每個未突變位元、每代；TSG 為每等位基因
    death_rate: float = 0.30
    driver_death_reduction: float = 0.08
    triple_probability: float = 0.04  # 存活後產生三個子細胞的機率
    driver_triple_bonus: float = 0.12
    # 模擬規模及治療設定：控制初始負荷、觀察時間、啟動門檻與額外殺傷率。
    initial_cells: int = 100
    generations: int = 75
    therapy_start_size: int = 20_000
    sensitive_therapy_kill: float = 0.90  # 自然死亡後的額外殺傷機率
    resistant_therapy_kill: float = 0.04
    max_cells: int = 1_000_000  # 計算停止門檻，不是環境承載量

    def validate(self):
        # 先檢查基因組大小、基因索引是否有效，以及 driver 門檻是否可達。
        if not 2 <= self.n_genes <= 20:
            raise ValueError("n_genes 必須介於 2–20")
        genes = self.driver_genes + self.suppressor_genes + (self.resistance_gene,)
        if len(set(genes)) != len(genes) or any(g < 0 or g >= self.n_genes for g in genes):
            raise ValueError("oncogene、suppressor 與 resistance 基因必須互異且位於基因組內")
        if not 1 <= self.driver_hits <= len(self.driver_genes):
            raise ValueError("driver_hits 超出可用 driver 數目")
        if not self.suppressor_genes or self.suppressor_hits not in (1, 2):
            raise ValueError("需要至少一個抑癌基因，suppressor_hits 必須為 1 或 2")
        # 逐一檢查機率參數，避免傳入二項分布時出現無效機率。
        for name in ("mutation_rate", "death_rate", "driver_death_reduction",
                     "triple_probability", "driver_triple_bonus",
                     "suppressor_death_reduction", "suppressor_triple_bonus",
                     "sensitive_therapy_kill", "resistant_therapy_kill"):
            if not 0 <= getattr(self, name) <= 1:
                raise ValueError(f"{name} 必須介於 0–1")
        # 額外檢查加總後的機率與模擬規模；單個參數有效不代表組合必然有效。
        if self.triple_probability + self.driver_triple_bonus + self.suppressor_triple_bonus > 1:
            raise ValueError("產生三個子細胞的總機率不能超過 1")
        if min(self.initial_cells, self.generations, self.therapy_start_size, self.max_cells) <= 0:
            raise ValueError("細胞數、代數與門檻必須為正")


# 封裝一次試驗的設定、逐代統計、完整基因型計數及停止／治療資訊。
@dataclass
class Result:
    config: Config
    seed: int
    history: list[dict]
    populations: list[dict[int, int]]
    therapy_generation: int | None
    stop_reason: str


def driver_count(genotype: int, config: Config) -> int:
    # 1 << g 選取第 g 個基因；位元 AND 判斷是否突變，再加總 driver 數。
    return sum(bool(genotype & (1 << g)) for g in config.driver_genes)


def suppressor_alleles(genotype: int, gene: int, config: Config) -> int:
    """回傳同一抑癌基因失活的等位基因數（0、1 或 2）。"""
    # 第一個位元沿用基因座索引；第二個放在 n_genes 後，不增加基因座數。
    second = config.n_genes + config.suppressor_genes.index(gene)
    return int(bool(genotype & (1 << gene))) + int(bool(genotype & (1 << second)))


def suppressor_loss(genotype: int, config: Config) -> bool:
    # 必須是同一基因達到門檻；兩個不同基因各一擊不等於同一基因兩擊。
    return any(suppressor_alleles(genotype, gene, config) >= config.suppressor_hits
               for gene in config.suppressor_genes)


def growth_probabilities(genotype: int, config: Config) -> tuple[float, float]:
    """分開計算 oncogene 活化及抑癌功能喪失的生長優勢，可相加。"""
    oncogene = driver_count(genotype, config) >= config.driver_hits
    suppressor = suppressor_loss(genotype, config)
    # 每一類優勢最多給一次；多個抑癌基因達標不再重複累加。
    death = max(0, config.death_rate - config.driver_death_reduction * oncogene
                - config.suppressor_death_reduction * suppressor)
    triple = (config.triple_probability + config.driver_triple_bonus * oncogene
              + config.suppressor_triple_bonus * suppressor)
    return death, triple


def mutate_offspring(population: dict[int, int], config: Config,
                     rng: np.random.Generator) -> dict[int, int]:
    """逐基因二項抽樣：與逐細胞獨立突變等價，允許一代多基因突變。"""
    # 複製輸入字典，避免突變抽樣直接改動呼叫端保存的族群。
    population = dict(population)
    # 依序處理一般基因與 TSG 的兩個等位基因；允許同一代兩個等位基因都突變。
    for gene in range(config.n_genes + len(config.suppressor_genes)):
        bit = 1 << gene
        # 將讀取與寫入分開，避免同一輪新增的基因型被重複處理。
        updated = dict(population)
        for genotype, count in population.items():
            # 已突變的基因保持不變；只替尚未突變的細胞抽樣。
            if genotype & bit:
                continue
            # 二項分布一次抽出此基因型的突變細胞數，省去逐細胞迴圈。
            mutants = int(rng.binomial(count, config.mutation_rate))
            if mutants:
                # 從舊基因型移出突變細胞，以位元 OR 設定新突變並合併相同基因型。
                updated[genotype] -= mutants
                updated[genotype | bit] = updated.get(genotype | bit, 0) + mutants
        # 刪除零細胞的項目；突變只重新分類細胞，不改變細胞總數。
        population = {g: n for g, n in updated.items() if n > 0}
    return population


def summarize(population: dict[int, int], generation: int, config: Config) -> dict:
    # 計算總負荷及各基因型比例；分母至少為 1，讓滅絕族群也能安全處理。
    total = sum(population.values())
    fractions = np.array(list(population.values()), dtype=float) / max(total, 1)
    # 抗藥基因位元決定敏感／抗藥分類，與是否已達 driver 門檻分開判斷。
    resistant = sum(n for g, n in population.items() if g & (1 << config.resistance_gene))
    # 彙整多樣性：richness 是基因型數；Shannon 用自然對數；Simpson 為 1−Σp²。
    # driver_fraction 只計 oncogene 達標；suppressor_fraction 計 TSG 達標。
    # suppressor_biallelic_fraction 計至少一個 TSG 雙等位基因失活；滅絕時皆為零。
    return dict(generation=generation, total=total, resistant=resistant,
                sensitive=total - resistant, richness=len(population),
                shannon=float(-np.sum(fractions * np.log(fractions))),
                simpson=float(1 - np.sum(fractions ** 2)) if total else 0.0,
                driver_fraction=sum(n for g, n in population.items()
                                    if driver_count(g, config) >= config.driver_hits) / max(total, 1),
                suppressor_fraction=sum(n for g, n in population.items()
                                        if suppressor_loss(g, config)) / max(total, 1),
                suppressor_biallelic_fraction=sum(n for g, n in population.items()
                    if any(suppressor_alleles(g, gene, config) == 2
                           for gene in config.suppressor_genes)) / max(total, 1))


def simulate(config: Config = Config(), seed: int = 42, therapy: bool = True,
             target_size: int | None = None) -> Result:
    """同步世代分枝過程；治療從首次達門檻後的轉移開始，持續至模擬結束。

    超過 max_cells 時保留完整該代，不截斷細胞；target_size 用於等大小比較。
    """
    # 先驗證設定與目標大小，再建立一次模擬專用的亂數產生器。
    config.validate()
    if target_size is not None and target_size <= 0:
        raise ValueError("target_size 必須為正")
    # 固定 seed 可重現同一試驗；基因型 0 代表所有基因都尚未突變。
    rng = np.random.default_rng(seed)
    population = {0: config.initial_cells}
    # history 存摘要，populations 存基因型快照；尚未治療時代數記為 None。
    history, populations = [], []
    therapy_generation = None
    stop_reason = "generation_limit"
    # 包含第 0 代初始狀態；每輪先記錄當前族群，再決定是否推進下一代。
    for generation in range(config.generations + 1):
        row = summarize(population, generation, config)
        history.append(row)
        populations.append(dict(population))
        # 依序處理滅絕、達取樣目標、達計算上限與觀察時間結束，保留停止原因。
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
        # 首次達到治療門檻就記下當前世代；接下來的轉移持續套用治療。
        if therapy and therapy_generation is None and row["total"] >= config.therapy_start_size:
            therapy_generation = generation
        # 依基因型批次處理所有父細胞；新世代完全由 offspring 取代。
        offspring = {}
        for genotype, count in population.items():
            # 將 oncogene 與抑癌基因的效果合併，再套用治療殺傷。
            death, triple = growth_probabilities(genotype, config)
            # 治療殺傷作用於自然存活者，因此用兩個存活機率相乘計算總死亡率。
            if therapy_generation is not None:
                resistant = bool(genotype & (1 << config.resistance_gene))
                kill = config.resistant_therapy_kill if resistant else config.sensitive_therapy_kill
                death = 1 - (1 - death) * (1 - kill)
            # 先抽出存活父細胞，再抽出其中產生三個子細胞者；其餘各產生兩個。
            survivors = int(rng.binomial(count, 1 - death))
            triples = int(rng.binomial(survivors, triple))
            # 每個存活者先貢獻兩個子細胞，三子細胞事件再各加一個。
            n = 2 * survivors + triples
            if n:
                offspring[genotype] = n
        # 所有分裂完成後才讓子代突變；新突變的生存效果從下一次轉移開始。
        population = mutate_offspring(offspring, config, rng)
    # 把逐代資料與停止資訊一起回傳，供繪圖和治療結果判讀使用。
    return Result(config, seed, history, populations, therapy_generation, stop_reason)


def therapy_outcome(result: Result) -> dict:
    """復發定義：先降至治療前負荷以下，再回到治療開始時負荷。"""
    # 未啟動治療的試驗單獨標記，不當作治療成功或失敗。
    start = result.therapy_generation
    if start is None:
        return dict(status="therapy_not_started", nadir=None, relapse_generation=None)
    # 以啟動治療時的負荷為基準，只搜尋治療開始後的世代。
    baseline = result.history[start]["total"]
    after = result.history[start + 1:]
    declined = False
    relapse = None
    # 記住是否曾低於基準；一旦之後回升至基準，就記錄第一次復發代數。
    for row in after:
        declined |= row["total"] < baseline
        if declined and row["total"] >= baseline:
            relapse = row["generation"]
            break
    # 區分復發、滅絕、尚未觀察到復發及無初始反應；nadir 是觀察期最低負荷。
    status = ("relapse" if relapse is not None else "extinction" if result.stop_reason == "extinction"
              else "no_relapse_observed" if declined else "no_initial_response")
    return dict(status=status, nadir=min((r["total"] for r in after), default=baseline),
                relapse_generation=relapse)


# 集中定義展示圖的配色，讓各情境與面板使用一致的視覺樣式。
COLORS = ["#167D9A", "#ED7953", "#6C5B9C", "#57A773", "#D6A43B", "#BC5B86"]


def set_plot_style():
    """統一簡報用視覺樣式；英文圖標籤避免依賴中文系統字型。"""
    # 統一字型、線寬、背景與淡色網格，並移除上方及右側邊框。
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11,
                         "axes.titlesize": 13, "axes.titleweight": "bold",
                         "axes.spines.top": False, "axes.spines.right": False,
                         "axes.labelcolor": "#334155", "text.color": "#172B42",
                         "axes.grid": True, "grid.alpha": 0.16,
                         "figure.facecolor": "#FAFBFD", "axes.facecolor": "white",
                         "savefig.facecolor": "#FAFBFD", "lines.linewidth": 2.3})


def save_figure(fig, output: Path, name: str):
    # 建立輸出目錄；同時支援呼叫端傳入字串或 Path 路徑。
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    # 同時輸出高解析度點陣 PNG 與可縮放向量 SVG，緊密裁切多餘留白。
    for extension in ("png", "svg"):
        fig.savefig(output / f"{name}.{extension}", dpi=240, bbox_inches="tight")
    # 儲存後釋放圖形資源，避免重複繪圖時累積記憶體。
    plt.close(fig)


def plot_dashboard(result: Result, output: Path):
    """輸出生長、基因型組成、多樣性與突變頻率四面板圖。"""
    # 建立 2×2 面板，使用自動排版協調標題、座標標籤與色條的空間。
    set_plot_style()
    fig, axes = plt.subplots(2, 2, figsize=(15, 11), layout="constrained")
    fig.suptitle("Tumor evolution | treatment & resistance", fontsize=21, fontweight="bold")
    t = np.array([r["generation"] for r in result.history])
    # 面板 A：並列總細胞、敏感細胞與抗藥細胞的逐代數量。
    ax = axes[0, 0]
    for key, label, color in (("total", "Total", "#172B42"),
                              ("sensitive", "Sensitive", COLORS[0]),
                              ("resistant", "Resistant", COLORS[1])):
        ax.plot(t, [r[key] for r in result.history], label=label, color=color)
    # symlog 同時呈現零細胞及跨數量級的增長，1 以下保留線性尺度。
    ax.set(yscale="symlog", ylabel="Cells (symlog; linear below 1)", title="A  Population response")
    ax.set_yscale("symlog", linthresh=1)
    # 以垂直虛線標示治療起點，淡色背景表示持續接受治療的區間。
    if result.therapy_generation is not None:
        ax.axvspan(result.therapy_generation, t[-1], color=COLORS[1], alpha=0.09)
        ax.axvline(result.therapy_generation, color=COLORS[1], ls="--", lw=1.5, label="Therapy begins")
    ax.legend(frameon=False, fontsize=9)
    # 面板 B：累加各世代的基因型比例，選出整段歷程中主要的五種基因型。
    ax = axes[0, 1]
    totals = np.array([r["total"] for r in result.history])
    scores = {}
    for pop, total in zip(result.populations, totals):
        for g, n in pop.items():
            scores[g] = scores.get(g, 0) + n / max(total, 1)
    top = sorted(scores, key=scores.get, reverse=True)[:5]
    # 將主要基因型轉為逐代比例，其餘合併為 Other；滅絕時總比例為零。
    bands = np.array([[p.get(g, 0) / max(n, 1) for p, n in zip(result.populations, totals)] for g in top])
    others = np.maximum(0, (totals > 0).astype(float) - bands.sum(axis=0))
    # 以十六進位顯示基因型 bitmask；含抗藥基因者附上 R，方便辨識。
    labels = [f"G{g:06X}" + (" (R)" if g & (1 << result.config.resistance_gene) else "") for g in top]
    ax.stackplot(t, *bands, others, labels=labels + ["Other"], colors=COLORS, alpha=0.9)
    ax.set(ylim=(0, 1), ylabel="Fraction of cells", title="B  Genotype composition")
    ax.legend(loc="upper left", fontsize=8, frameon=False, ncol=2)
    # 面板 C：exp(Shannon H) 將熵轉為有效基因型數，滅絕時另外定義為零。
    ax = axes[1, 0]
    ax.plot(t, [np.exp(r["shannon"]) if r["total"] else 0 for r in result.history], color=COLORS[2])
    ax.set(ylabel="Effective genotypes: exp(Shannon H)", title="C  Genetic heterogeneity")
    # 面板 D：TSG 顯示任一等位基因失活的細胞比例，不代表已達功能喪失門檻。
    ax = axes[1, 1]
    frequency = np.array([[sum(n for g, n in p.items() if (suppressor_alleles(g, gene, result.config) > 0
                               if gene in result.config.suppressor_genes else bool(g & (1 << gene)))) / max(total, 1)
                           for p, total in zip(result.populations, totals)]
                          for gene in range(result.config.n_genes)])
    # 用網格邊界繪製熱圖，固定色階 0–1，便於比較不同基因的突變頻率。
    mesh = ax.pcolormesh(np.arange(len(t) + 1) - 0.5,
                         np.arange(result.config.n_genes + 1) - 0.5,
                         frequency, cmap="YlGnBu", vmin=0, vmax=1, shading="flat")
    # 顯示時將零起算索引改成 Gene 1 起算，並標出 driver 與 resistance。
    labels = [f"Gene {g + 1}" + (" · driver" if g in result.config.driver_genes else
                                  " · resistance" if g == result.config.resistance_gene else
                                  " · TSG (any hit)" if g in result.config.suppressor_genes else "")
              for g in range(result.config.n_genes)]
    ax.set(yticks=range(result.config.n_genes), yticklabels=labels, title="D  Mutation frequencies")
    ax.tick_params(axis="y", labelsize=8)
    ax.grid(False)
    fig.colorbar(mesh, ax=ax, label="Mutant fraction", pad=0.02)
    # 所有面板共用世代標籤及整數刻度，最後一次輸出整張展示圖。
    for ax in axes.flat:
        ax.set_xlabel("Generation")
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    save_figure(fig, output, "01_tumor_dashboard")


def sample_at_size(result: Result, size: int, rng: np.random.Generator) -> dict | None:
    """跨門檻的該代中無放回取固定數目，避免不同生長速度的 overshoot 偏差。

    這是等大小的腫瘤樣本，不是連續時間中恰好達到該大小的完整腫瘤。
    """
    # 最後一代未達取樣數量時回傳 None，由比較函數記錄為未達目標。
    if result.history[-1]["total"] < size:
        return None
    # 多變量超幾何分布等同從整個族群無放回取樣，維持固定樣本大小。
    pop = result.populations[-1]
    counts = rng.multivariate_hypergeometric(np.array(list(pop.values()), dtype=np.int64), size)
    # 重建樣本的基因型計數，沿用同一統計函數計算多樣性。
    sampled = {g: int(n) for g, n in zip(pop, counts) if n}
    return summarize(sampled, result.history[-1]["generation"], result.config)


def run_comparisons(config: Config, repeats: int, seed: int, output: Path) -> list[dict]:
    """重複模擬：驅動突變對生長，以及死亡率對等大小樣本多樣性的影響。"""
    # 重複次數必須為正，才有可供比較的試驗分布。
    if repeats < 1:
        raise ValueError("repeats 必須至少為 1")
    # 建立左右兩個比較面板，rows 收集每次試驗的結果供 CSV 匯出。
    set_plot_style()
    fig, axes = plt.subplots(1, 2, figsize=(13.5, 5.4), layout="constrained")
    fig.suptitle("What changes tumor growth and diversity?", fontsize=20, fontweight="bold")
    rows = []
    # 生長比較固定觀察 18 代，固定 oncogene 設定，比較無 TSG 優勢、TSG 一擊與兩擊條件。
    growth_horizon = 18
    # 保留同批試驗的逐代資料，供 one-hit／two-hit 專屬比較圖重用。
    hit_rows = []
    for condition, label in enumerate(("No TSG advantage", "TSG one-hit", "TSG two-hit")):
        # 用 replace 建立新設定而不修改原始 Config；控制組只取消 TSG 優勢；所有情境保留相同 oncogene 效果與突變率。
        cfg = replace(config, generations=growth_horizon, max_cells=10**12, suppressor_hits=1 if condition == 1 else 2)
        if condition == 0:
            cfg = replace(cfg, suppressor_death_reduction=0, suppressor_triple_bonus=0)
        # 收集每次試驗的完整生長曲線；不同條件與重複試驗使用不同種子。
        trajectories = []
        for rep in range(repeats):
            run_seed = seed + condition * 10000 + rep
            result = simulate(cfg, seed=run_seed, therapy=False)
            # 提早滅絕後補零；未滅絕的截短試驗仍由下方安全上限檢查拒絕。
            if condition in (1, 2):
                for generation in range(growth_horizon + 1):
                    row = result.history[min(generation, len(result.history) - 1)]
                    hit_rows.append(dict(scenario=label, replicate=rep, seed=run_seed,
                                         generation=generation, total=row["total"],
                                         suppressor_fraction=row["suppressor_fraction"]))
            values = [r["total"] for r in result.history]
            # 若碰到安全上限就停止比較並報錯，以免把截短的曲線當成完整資料。
            if result.stop_reason == "cell_limit":
                raise RuntimeError("Growth comparison hit safety limit; shorten growth_horizon")
            # 提早滅絕的試驗以零補齊剩餘世代，使各次曲線長度一致。
            values.extend([0] * (growth_horizon + 1 - len(values)))
            trajectories.append(values)
            rows.append(dict(experiment="growth", condition=label, replicate=rep,
                             seed=run_seed, reached_target="", generation=result.history[-1]["generation"],
                             total=result.history[-1]["total"], shannon=result.history[-1]["shannon"],
                             suppressor_fraction=result.history[-1]["suppressor_fraction"]))
        # 逐代計算中位數與 10–90 百分位；陰影是試驗間變異，不是信賴區間。
        low, median, high = np.percentile(trajectories, [10, 50, 90], axis=0)
        axes[0].plot(range(growth_horizon + 1), median, color=COLORS[condition], label=label)
        axes[0].fill_between(range(growth_horizon + 1), low, high, color=COLORS[condition], alpha=0.14)
    # 完成生長面板的對數尺度與圖例，讓三種情境易於比較。
    axes[0].set(title="A  Tumor suppressor hit requirement", xlabel="Generation", ylabel="Cells · median & 10–90% interval")
    axes[0].set_yscale("symlog", linthresh=1)
    axes[0].legend(frameon=False, fontsize=9)
    # 死亡率比較：首次達標後固定取樣 5,000 個細胞，降低跨門檻大小差異的影響。
    target = 5000
    for index, death in enumerate((0.20, 0.35, 0.45)):
        # 每種死亡率各做多次試驗，只有達標樣本的 Shannon H 會加入散點圖。
        values = []
        for rep in range(repeats):
            # 使用與生長比較不同的種子區段；取樣也使用獨立的亂數產生器。
            run_seed = seed + 50000 + index * 10000 + rep
            cfg = replace(config, death_rate=death, generations=100, max_cells=max(target * 4, config.max_cells))
            result = simulate(cfg, seed=run_seed, therapy=False, target_size=target)
            sample = sample_at_size(result, target, np.random.default_rng(run_seed + 1_000_000))
            if sample is not None:
                values.append(sample["shannon"])
            # 未達門檻者仍寫入 CSV，保留達標與否及停止世代，避免隱藏失敗試驗。
            rows.append(dict(experiment="death_rate", condition=death, replicate=rep,
                             seed=run_seed, reached_target=sample is not None,
                             generation=result.history[-1]["generation"],
                             total=result.history[-1]["total"], shannon=sample["shannon"] if sample else "",
                             suppressor_fraction=sample["suppressor_fraction"] if sample else ""))
        # 散點左右微幅偏移只改善重疊顯示；黑色橫線標示該組中位數。
        if values:
            jitter = np.random.default_rng(seed + index).uniform(-0.10, 0.10, len(values))
            axes[1].scatter(index + jitter, values, color=COLORS[index], s=34, alpha=0.7)
            axes[1].plot([index - 0.18, index + 0.18], [np.median(values)] * 2, color="#172B42", lw=3)
        # 標示達到目標的試驗數，提醒多樣性統計是以達標為條件。
        axes[1].text(index, 0.98, f"{len(values)}/{repeats} reached", ha="center", va="top",
                     transform=axes[1].get_xaxis_transform(), fontsize=9)
    # 設定死亡率面板的軸標籤與留白，再儲存圖表及逐次試驗資料。
    axes[1].set(title=f"B  Diversity in {target:,}-cell samples", xlabel="Natural death probability per generation",
                ylabel="Shannon diversity (natural log)", xticks=[0, 1, 2], xticklabels=["0.20", "0.35", "0.45"])
    axes[1].margins(y=0.25, x=0.20)
    save_figure(fig, output, "02_scenario_comparisons")
    write_csv(output / "comparisons.csv", rows)
    plot_hit_comparison(hit_rows, config, output)
    return rows


def plot_hit_comparison(rows: list[dict], config: Config, output: Path):
    """用同批重複試驗直接比較生長、TSG 達標比例及固定世代的負荷。"""
    # 圖 A 採線性尺度凸顯絕對數量差異；圖 B 顯示機制差異；圖 C 保留個別試驗。
    set_plot_style()
    output = Path(output)
    fig, axes = plt.subplots(1, 3, figsize=(17, 5.8), layout="constrained")
    last = max(r["generation"] for r in rows)
    repeats = len({r["replicate"] for r in rows})
    fig.suptitle("One-hit vs two-hit | tumor suppressor loss", fontsize=21, fontweight="bold")
    fig.supxlabel(f"{config.n_genes} genes · no therapy · identical oncogene effects & mutation rates · "
                  f"{repeats} trials/scenario\n"
                  "Lines: median; bands: 10–90% of trials (not confidence intervals). "
                  "One-hit: either allele; two-hit: both alleles of the same TSG.", fontsize=10)
    summaries = []
    for index, (label, color) in enumerate((("TSG one-hit", COLORS[1]), ("TSG two-hit", COLORS[2]))):
        # 依試驗與世代建立矩陣，兩組均使用相同時間範圍及完整重複試驗。
        subset = [r for r in rows if r["scenario"] == label]
        series = [[r for r in subset if r["replicate"] == rep] for rep in range(repeats)]
        total = np.array([[r["total"] for r in trial] for trial in series])
        fraction = 100 * np.array([[r["suppressor_fraction"] for r in trial] for trial in series])
        for ax, values in zip(axes[:2], (total, fraction)):
            low, median, high = np.percentile(values, [10, 50, 90], axis=0)
            ax.plot(range(last + 1), median, color=color, label=label,
                    ls="-" if index == 0 else "--")
            ax.fill_between(range(last + 1), low, high, color=color, alpha=0.15)
        # 固定終點比較避免拿不同世代的負荷相比；散點左右位移僅為減少遮蔽。
        endpoint = total[:, -1]
        jitter = np.random.default_rng(900 + index).uniform(-0.10, 0.10, repeats)
        axes[2].scatter(index + jitter, endpoint, color=color, alpha=0.7, s=40)
        median = float(np.median(endpoint))
        axes[2].plot([index - 0.2, index + 0.2], [median] * 2, color=color, lw=3)
        axes[2].text(index, 0.98, f"Median: {median:,.0f} cells", ha="center", va="top",
                     transform=axes[2].get_xaxis_transform(), fontsize=10)
        summaries.append(dict(scenario=label, generation=last, replicates=repeats,
                              median_cells=median,
                              median_tsg_percent=float(np.median(fraction[:, -1]))))
    # TSG 比例差距可能跨數量級；symlog 可同時保留零值和低頻 two-hit 細胞。
    axes[0].set(title="A  Tumor growth (linear scale)", xlabel="Generation", ylabel="Total cells")
    axes[1].set(title="B  Cells meeting the TSG threshold", xlabel="Generation", ylabel="Cells with TSG loss (%)")
    axes[1].set_yscale("symlog", linthresh=0.01)
    axes[1].text(0.03, 0.97, "Symlog; linear below 0.01%", transform=axes[1].transAxes,
                 va="top", fontsize=9)
    axes[2].set(title=f"C  Tumor burden at generation {last}", ylabel="Total cells",
                xticks=[0, 1], xticklabels=["One-hit", "Two-hit"], xlim=(-0.5, 1.5))
    axes[2].margins(y=0.25)
    for ax in axes[:2]:
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))
        ax.legend(frameon=False, loc="lower right", fontsize=9)
    save_figure(fig, output, "03_one_hit_vs_two_hit")
    # 匯出逐代及終點摘要，圖中的差異可直接追溯到原始重複試驗。
    write_csv(output / "one_hit_vs_two_hit.csv", rows)
    write_csv(output / "one_hit_vs_two_hit_summary.csv", summaries)


def write_csv(path: Path, rows: list[dict]):
    # 以第一筆資料的鍵建立欄名；UTF-8 支援中文，newline 避免 CSV 多餘空行。
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main():
    # 命令列入口：允許指定亂數種子、重複次數及輸出資料夾。
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--suppressor-hits", type=int, choices=(1, 2), default=2,
                        help="TSG allele-loss threshold for the main treatment scenario")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--replicates", type=int, default=12)
    parser.add_argument("--output", type=Path, default=Path("outputs"))
    # 解析並檢查命令列輸入，建立後續輸出的資料夾。
    args = parser.parse_args()
    if args.replicates < 1:
        parser.error("--replicates must be >= 1")
    args.output.mkdir(parents=True, exist_ok=True)
    # 執行預設治療案例，輸出四面板圖與逐代統計，再執行情境比較。
    config = Config(suppressor_hits=args.suppressor_hits)
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
    # 保存完整參數、種子、套件版本與主要結果，方便追溯及重現本次執行。
    metadata = dict(config=asdict(config), seed=args.seed, replicates=args.replicates,
                    stop_reason=result.stop_reason, therapy_generation=result.therapy_generation,
                    outcome=therapy_outcome(result), numpy_version=np.__version__,
                    matplotlib_version=matplotlib.__version__)
    (args.output / "run_metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    # 在終端顯示結果摘要和輸出位置，方便找到生成的圖表與資料。
    print(json.dumps(metadata, indent=2))
    print(f"Figures and data saved to: {args.output.resolve()}")


# 只在直接執行本檔案時啟動 main；被其他程式 import 時不自動跑模擬。
if __name__ == "__main__":
    main()
