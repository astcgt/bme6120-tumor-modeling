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
from dataclasses import asdict, dataclass, field, replace
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
    # IL11 開關：啟用時將第 7 個基因座由 passenger 改為 IL11；關閉時保留中性位點。
    il11_enabled: bool = False
    il11_gene: int = 6
    il11_max_death_reduction: float = 0.06
    il11_half_fraction: float = 0.02  # IL11+ 占 2% 時達最大共享效果的一半
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
    # 延伸分析用（預設關閉，主結果不受影響）：
    # adaptive：負荷降到首次治療負荷的 adaptive_pause_fraction 以下就暫停，回到治療門檻再恢復。
    therapy_mode: str = "continuous"
    adaptive_pause_fraction: float = 0.5
    resistance_cost: float = 0.0  # 抗藥細胞額外的自然死亡機率（抗藥代價）
    # 承載量 K：每個細胞額外死亡 c·N/K，c 使無突變細胞在 N=K 時平均子代剛好為 1。None＝無競爭。
    carrying_capacity: int | None = None

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
        # 保留 20 個基因座及相同突變機會，避免有無 IL11 比較混入基因組大小差異。
        if not isinstance(self.il11_enabled, bool):
            raise ValueError("il11_enabled 必須為 bool")
        if self.il11_gene in genes or not 0 <= self.il11_gene < self.n_genes:
            raise ValueError("IL11 必須是獨立且有效的基因座")
        if not 0 < self.il11_half_fraction <= 1:
            raise ValueError("il11_half_fraction 必須介於 0（不含）與 1")
        # 逐一檢查機率參數，避免傳入二項分布時出現無效機率。
        for name in ("mutation_rate", "death_rate", "driver_death_reduction",
                     "triple_probability", "driver_triple_bonus",
                     "suppressor_death_reduction", "suppressor_triple_bonus", "il11_max_death_reduction",
                     "sensitive_therapy_kill", "resistant_therapy_kill", "resistance_cost"):
            if not 0 <= getattr(self, name) <= 1:
                raise ValueError(f"{name} 必須介於 0–1")
        # 額外檢查加總後的機率與模擬規模；單個參數有效不代表組合必然有效。
        if self.triple_probability + self.driver_triple_bonus + self.suppressor_triple_bonus > 1:
            raise ValueError("產生三個子細胞的總機率不能超過 1")
        if min(self.initial_cells, self.generations, self.therapy_start_size, self.max_cells) <= 0:
            raise ValueError("細胞數、代數與門檻必須為正")
        if self.therapy_mode not in ("continuous", "adaptive"):
            raise ValueError("therapy_mode 必須是 continuous 或 adaptive")
        if not 0 < self.adaptive_pause_fraction < 1:
            raise ValueError("adaptive_pause_fraction 必須介於 0 與 1 之間")
        if self.carrying_capacity is not None and self.carrying_capacity <= 0:
            raise ValueError("carrying_capacity 必須為正或 None")


# 封裝一次試驗的設定、逐代統計、完整基因型計數及停止／治療資訊。
@dataclass
class Result:
    config: Config
    seed: int
    history: list[dict]
    populations: list[dict[int, int]]
    therapy_generation: int | None
    stop_reason: str
    # 每個已推進世代是否施加治療；continuous 時從 therapy_generation 起全為 True。
    therapy_active: list[bool] = field(default_factory=list)


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


def growth_probabilities(genotype: int, config: Config, shared_benefit: float = 0.0) -> tuple[float, float]:
    """分開計算 oncogene 活化及抑癌功能喪失的生長優勢，可相加。"""
    oncogene = driver_count(genotype, config) >= config.driver_hits
    suppressor = suppressor_loss(genotype, config)
    # 每一類優勢最多給一次；多個抑癌基因達標不再重複累加。
    death = max(0, config.death_rate - config.driver_death_reduction * oncogene
                - config.suppressor_death_reduction * suppressor - shared_benefit)
    # 抗藥代價只在明確設定時加上，預設為 0 時完全不改變計算。
    if config.resistance_cost and genotype & (1 << config.resistance_gene):
        death = min(1.0, death + config.resistance_cost)
    triple = (config.triple_probability + config.driver_triple_bonus * oncogene
              + config.suppressor_triple_bonus * suppressor)
    return death, triple


def il11_environment(population: dict[int, int], config: Config) -> tuple[float, float]:
    """由當前 IL11+ 細胞比例計算下一次轉移的共享死亡率降低量。"""
    # 關閉時此位點只是 passenger，不產生 IL11，也不提供共享優勢。
    total = sum(population.values())
    if not config.il11_enabled or total == 0:
        return 0.0, 0.0
    fraction = sum(n for g, n in population.items() if g & (1 << config.il11_gene)) / total
    # 飽和函數是教學假設；均勻混合環境中所有細胞（含非生產者）皆受益。
    benefit = config.il11_max_death_reduction * fraction / (config.il11_half_fraction + fraction)
    return fraction, benefit


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
    il11_fraction, shared_benefit = il11_environment(population, config)
    # 彙整多樣性：richness 是基因型數；Shannon 用自然對數；Simpson 為 1−Σp²。
    # driver_fraction 只計 oncogene 達標；suppressor_fraction 計 TSG 達標。
    # suppressor_biallelic_fraction 計至少一個 TSG 雙等位基因失活；滅絕時皆為零。
    return dict(generation=generation, total=total, resistant=resistant,
                il11_fraction=il11_fraction, il11_shared_benefit=shared_benefit,
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


def crowding_death(total: int, config: Config) -> float:
    """承載量競爭的額外死亡機率；係數讓無突變細胞在 N = K 時期望子代恰為 1。"""
    if config.carrying_capacity is None:
        return 0.0
    coefficient = (1 - config.death_rate) - 1 / (2 + config.triple_probability)
    return coefficient * total / config.carrying_capacity


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
    therapy_on = False
    active = []
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
            therapy_on = True
        # adaptive：降到首次治療負荷的一定比例以下就暫停，回到治療門檻再恢復。
        elif therapy_generation is not None and config.therapy_mode == "adaptive":
            reference = history[therapy_generation]["total"]
            if therapy_on and row["total"] <= config.adaptive_pause_fraction * reference:
                therapy_on = False
            elif not therapy_on and row["total"] >= config.therapy_start_size:
                therapy_on = True
        active.append(therapy_on)
        crowding = crowding_death(row["total"], config)
        # 依基因型批次處理所有父細胞；新世代完全由 offspring 取代。
        # 每代只計算一次共享環境；本代新突變會從下一次轉移開始提供 IL11。
        _, shared_benefit = il11_environment(population, config)
        offspring = {}
        for genotype, count in population.items():
            # 將 oncogene 與抑癌基因的效果合併，再套用治療殺傷。
            death, triple = growth_probabilities(genotype, config, shared_benefit)
            if crowding:
                death = min(1.0, death + crowding)
            # 治療殺傷作用於自然存活者，因此用兩個存活機率相乘計算總死亡率。
            if therapy_on:
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
    return Result(config, seed, history, populations, therapy_generation, stop_reason, active)


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
# TSG 比較的三種條件；第一個是對照組（TSG 失活沒有生長優勢）。
HIT_SCENARIOS = ("No TSG advantage", "TSG one-hit", "TSG two-hit")
HIT_SHORT = {"No TSG advantage": "Control\n(no TSG advantage)", "TSG one-hit": "One-hit", "TSG two-hit": "Two-hit"}


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
                                  " · TSG (any hit)" if g in result.config.suppressor_genes else
                                  " · IL11" if result.config.il11_enabled and g == result.config.il11_gene else "")
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
    """重複模擬：TSG 規則對生長（畫在圖 03），以及死亡率對等大小樣本多樣性的影響（圖 02）。"""
    # 重複次數必須為正，才有可供比較的試驗分布。
    if repeats < 1:
        raise ValueError("repeats 必須至少為 1")
    rows = []
    # 生長比較固定觀察 18 代，固定 oncogene 設定，比較無 TSG 優勢、TSG 一擊與兩擊條件。
    growth_horizon = 18
    # 三種條件的逐代資料都交給圖 03；無 TSG 優勢是對照組。
    hit_rows = []
    for condition, label in enumerate(HIT_SCENARIOS):
        # 用 replace 建立新設定而不修改原始 Config；控制組只取消 TSG 優勢；所有情境保留相同 oncogene 效果與突變率。
        cfg = replace(config, generations=growth_horizon, max_cells=10**12, suppressor_hits=1 if condition == 1 else 2)
        if condition == 0:
            cfg = replace(cfg, suppressor_death_reduction=0, suppressor_triple_bonus=0)
        for rep in range(repeats):
            run_seed = seed + condition * 10000 + rep
            result = simulate(cfg, seed=run_seed, therapy=False)
            # 若碰到安全上限就停止比較並報錯，以免把截短的曲線當成完整資料。
            if result.stop_reason == "cell_limit":
                raise RuntimeError("Growth comparison hit safety limit; shorten growth_horizon")
            # 提早滅絕的試驗沿用最後一列（零細胞），使各次曲線長度一致。
            for generation in range(growth_horizon + 1):
                row = result.history[min(generation, len(result.history) - 1)]
                hit_rows.append(dict(scenario=label, replicate=rep, seed=run_seed,
                                     generation=generation, total=row["total"],
                                     suppressor_fraction=row["suppressor_fraction"]))
            rows.append(dict(experiment="growth", condition=label, replicate=rep,
                             seed=run_seed, reached_target="", generation=result.history[-1]["generation"],
                             total=result.history[-1]["total"], shannon=result.history[-1]["shannon"],
                             suppressor_fraction=result.history[-1]["suppressor_fraction"]))
    # 死亡率比較：首次達標後固定取樣 5,000 個細胞，降低跨門檻大小差異的影響。
    target, deaths = 5000, (0.20, 0.35, 0.45)
    shannon, generations = {}, {}
    for index, death in enumerate(deaths):
        shannon[death], generations[death] = [], []
        for rep in range(repeats):
            # 使用與生長比較不同的種子區段；取樣也使用獨立的亂數產生器。
            run_seed = seed + 50000 + index * 10000 + rep
            cfg = replace(config, death_rate=death, generations=100, max_cells=max(target * 4, config.max_cells))
            result = simulate(cfg, seed=run_seed, therapy=False, target_size=target)
            sample = sample_at_size(result, target, np.random.default_rng(run_seed + 1_000_000))
            if sample is not None:
                shannon[death].append(sample["shannon"])
                generations[death].append(result.history[-1]["generation"])
            # 未達門檻者仍寫入 CSV，保留達標與否及停止世代，避免隱藏失敗試驗。
            rows.append(dict(experiment="death_rate", condition=death, replicate=rep,
                             seed=run_seed, reached_target=sample is not None,
                             generation=result.history[-1]["generation"],
                             total=result.history[-1]["total"], shannon=sample["shannon"] if sample else "",
                             suppressor_fraction=sample["suppressor_fraction"] if sample else ""))
    plot_death_rate(shannon, generations, config, repeats, target, seed, output)
    write_csv(output / "comparisons.csv", rows)
    plot_hit_comparison(hit_rows, config, output)
    return rows


def plot_death_rate(shannon: dict, generations: dict, config: Config, repeats: int, target: int,
                    seed: int, output: Path):
    """同樣 5,000 個細胞：死亡率越高越異質；B 顯示機制是需要更多代（更多分裂）才長到同樣大小。"""
    set_plot_style()
    fig, axes = plt.subplots(1, 2, figsize=(14, 5.8), layout="constrained")
    fig.suptitle("Higher cell turnover → more heterogeneity at the same tumor size", fontsize=20, fontweight="bold")
    deaths = list(shannon)
    for ax, values, key in ((axes[0], shannon, "shannon"), (axes[1], generations, "generations")):
        for index, death in enumerate(deaths):
            data = values[death]
            if not data:
                continue
            # 散點左右微幅偏移只改善重疊顯示；黑色橫線標示該組中位數。
            jitter = np.random.default_rng(seed + index).uniform(-0.10, 0.10, len(data))
            ax.scatter(index + jitter, data, color=COLORS[index], s=34, alpha=0.7)
            median = float(np.median(data))
            ax.plot([index - 0.18, index + 0.18], [median] * 2, color="#172B42", lw=3)
            ax.annotate(f"{median:.2f}" if key == "shannon" else f"{median:.0f} gen", (index + 0.22, median),
                        va="center", fontsize=10.5, color="#172B42")
            # 標示達到目標的試驗數，提醒統計是以達標為條件。
            ax.text(index, 0.98, f"{len(data)}/{repeats} reached", ha="center", va="top",
                    transform=ax.get_xaxis_transform(), fontsize=9)
    # 無突變細胞的每代淨倍增 (1-d)(2+q)：死亡率越高，越接近 1，達標需要越多代。
    labels = [f"{d:.2f}\n×{(1 - d) * (2 + config.triple_probability):.2f}/gen" for d in deaths]
    axes[0].set(title=f"A  Diversity in {target:,}-cell samples", ylabel="Shannon diversity (natural log)")
    axes[1].set(title=f"B  Generations needed to reach {target:,} cells", ylabel="Generation first ≥ target")
    axes[1].set_ylim(bottom=0)
    for ax in axes:
        ax.set(xticks=range(len(deaths)), xticklabels=labels,
               xlabel="Natural death probability (net growth of a mutation-free cell)")
        ax.margins(y=0.2, x=0.20)
        # 留出頂端空間給「達標數」文字，避免壓到散點。
        low, high = ax.get_ylim(); ax.set_ylim(low, high + 0.12 * (high - low))
    fig.supxlabel(f"No therapy · {repeats} trials per death rate · each trial sampled once at its first generation ≥ {target:,} cells "
                  f"(without replacement).\nMore generations → more divisions before reaching the same size → more mutations "
                  "accumulated → higher diversity.", fontsize=10)
    save_figure(fig, output, "02_death_rate_heterogeneity")


def plot_hit_comparison(rows: list[dict], config: Config, output: Path):
    """用同批重複試驗直接比較生長、TSG 達標比例及固定世代的負荷；含無 TSG 優勢對照組。"""
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
                  "One-hit: either allele; two-hit: both alleles of the same TSG. "
                  "Control: TSG loss gives no advantage (two-hit counting).", fontsize=10)
    summaries = []
    styles = {HIT_SCENARIOS[0]: ("#64748B", ":"), HIT_SCENARIOS[1]: (COLORS[1], "-"), HIT_SCENARIOS[2]: (COLORS[2], "--")}
    order = [s for s in HIT_SCENARIOS if any(r["scenario"] == s for r in rows)]
    for index, label in enumerate(order):
        color, linestyle = styles[label]
        # 依試驗與世代建立矩陣，各組均使用相同時間範圍及完整重複試驗。
        subset = [r for r in rows if r["scenario"] == label]
        series = [[r for r in subset if r["replicate"] == rep] for rep in range(repeats)]
        total = np.array([[r["total"] for r in trial] for trial in series])
        fraction = 100 * np.array([[r["suppressor_fraction"] for r in trial] for trial in series])
        for ax, values in zip(axes[:2], (total, fraction)):
            low, median, high = np.percentile(values, [10, 50, 90], axis=0)
            ax.plot(range(last + 1), median, color=color, label=label, ls=linestyle)
            ax.fill_between(range(last + 1), low, high, color=color, alpha=0.10 if index == 0 else 0.15)
        # 固定終點比較避免拿不同世代的負荷相比；散點左右位移僅為減少遮蔽。
        endpoint = total[:, -1]
        jitter = np.random.default_rng(900 + index).uniform(-0.10, 0.10, repeats)
        axes[2].scatter(index + jitter, endpoint, color=color, alpha=0.7, s=40)
        median = float(np.median(endpoint))
        axes[2].plot([index - 0.2, index + 0.2], [median] * 2, color=color, lw=3)
        axes[2].text(index, 0.98, f"Median\n{median:,.0f}", ha="center", va="top",
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
                xticks=range(len(order)), xticklabels=[HIT_SHORT[s] for s in order], xlim=(-0.5, len(order) - 0.5))
    axes[2].margins(y=0.25)
    for ax in axes[:2]:
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))
        ax.legend(frameon=False, loc="lower right", fontsize=9)
    save_figure(fig, output, "03_one_hit_vs_two_hit")
    # 匯出逐代及終點摘要，圖中的差異可直接追溯到原始重複試驗。
    write_csv(output / "one_hit_vs_two_hit.csv", rows)
    write_csv(output / "one_hit_vs_two_hit_summary.csv", summaries)


def run_il11_comparison(config: Config, repeats: int, seed: int, output: Path, horizon: int = 18):
    """比較相同 20 個基因座下，IL11 功能有無的生長；並把實際生產者比例標在效益曲線上。

    C 直接畫模型公式 b(f) = Bmax f/(K+f)（教學假設，不是實驗擬合）；共享死亡率降低量與
    非生產者適應度都是 B 的確定性轉換，所以只寫進 CSV 與註腳，不另畫面板。
    """
    if repeats < 1 or horizon < 1:
        raise ValueError("repeats 與 horizon 必須為正")
    config.validate()
    set_plot_style()
    fig, axes = plt.subplots(1, 3, figsize=(18.5, 6.2), layout="constrained")
    fig.suptitle("IL11+ producers | a few cells speed up the whole tumor", fontsize=21, fontweight="bold")
    rows, summaries, producer = [], [], None
    for enabled, label, color in ((False, "IL11 OFF", COLORS[0]), (True, "IL11 ON", COLORS[1])):
        # 固定時間窗與其餘參數；同一 replicate 用同一 seed，但分歧後不視為逐細胞配對。
        cfg = replace(config, il11_enabled=enabled, generations=horizon, max_cells=10**12)
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
                row = dict(condition=label, replicate=rep, seed=run_seed, generation=generation,
                           total=state["total"], il11_percent=100 * state["il11_fraction"],
                           shared_death_reduction=benefit,
                           nonproducer_expected_offspring=(1 - death) * (2 + triple))
                rows.append(row)
                trial.append(row)
            trials.append(trial)
        keys = ("total", "il11_percent", "shared_death_reduction", "nonproducer_expected_offspring")
        for ax, key in zip(axes[:2], keys[:2]):
            values = np.array([[r[key] for r in trial] for trial in trials])
            low, med, high = np.percentile(values, [10, 50, 90], axis=0)
            ax.plot(range(horizon + 1), med, label=label, color=color)
            ax.fill_between(range(horizon + 1), low, high, color=color, alpha=0.15)
        if enabled:
            producer = np.array([[r["il11_percent"] for r in trial] for trial in trials])
        summaries.append(dict(condition=label, generation=horizon, replicates=repeats,
                              **{f"median_{key}": float(np.median([t[-1][key] for t in trials])) for key in keys}))
    for ax, title, ylabel in zip(axes[:2], ("A  Whole-tumor growth", "B  IL11+ producer cells"),
                                 ("Total cells", "IL11+ cells (%)")):
        ax.set(title=title, xlabel="Generation", ylabel=ylabel)
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))
        ax.legend(frameon=False, loc="upper left")
    off, on = summaries
    axes[0].annotate(f"Gen {horizon}: {on['median_total']:,.0f} vs {off['median_total']:,.0f}\n"
                     f"(+{100 * (on['median_total'] / off['median_total'] - 1):.0f}% with IL11)",
                     xy=(horizon, on["median_total"]), xytext=(0.30, 0.55), textcoords="axes fraction",
                     arrowprops=dict(arrowstyle="->", color=COLORS[1]), fontsize=10.5, color=COLORS[1])
    # C：公式曲線（放大 0–10%），並把 IL11 ON 試驗實際出現的生產者比例範圍疊上去。
    ax = axes[2]
    maximum, half = config.il11_max_death_reduction, config.il11_half_fraction
    fraction = np.linspace(0, 1, 1001)
    benefit = maximum * fraction / (half + fraction)
    zoom = min(1.0, max(0.10, 3 * half))
    lo, hi = np.percentile(producer[:, 1:], [10, 90])
    ax.axvspan(lo, hi, color=COLORS[1], alpha=0.13, lw=0,
               label=f"Observed producers, gen 1–{horizon} (10–90%)")
    ax.plot(100 * fraction, 100 * benefit, color=COLORS[1], label=r"$b(f)=B_{max}\,f/(K+f)$")
    ax.axhline(100 * maximum, color="#64748B", ls="--", lw=1.4, label=f"Cap {100 * maximum:g} pp")
    ax.plot([100 * half, 100 * half, 0], [0, 50 * maximum, 50 * maximum], color=COLORS[2], ls=":", lw=1.5)
    ax.scatter([100 * half], [50 * maximum], color=COLORS[2], s=55, zorder=5)
    ax.annotate(f"Half-saturation K = {100 * half:g}%", xy=(100 * half, 50 * maximum), xytext=(4.2, 1.6),
                arrowprops=dict(arrowstyle="->", color=COLORS[2]), fontsize=10.5, color=COLORS[2])
    end = float(np.median(producer[:, -1]))
    ax.scatter([end], [100 * maximum * end / 100 / (half + end / 100)], marker="D", color="#172B42", s=45, zorder=6,
               label=f"Median at gen {horizon}: {end:.1f}%")
    ax.set(xlim=(0, 100 * zoom), ylim=(-0.2, 115 * maximum), title="C  Where the tumor sits on the benefit curve",
           xlabel="IL11+ cells (% of tumor)", ylabel="Shared reduction in death probability (pp)")
    ax.legend(frameon=False, loc="lower right", fontsize=9)
    fig.supxlabel(f"No therapy · {repeats} trials/condition · median and 10–90% interval (not CI). OFF: neutral passenger at the same locus; "
                  f"ON: IL11 mutation creates producers and ALL cells benefit.\nAt gen {horizon}: shared death reduction "
                  f"{100 * on['median_shared_death_reduction']:.1f} pp; a mutation-free nonproducer's expected offspring "
                  f"{off['median_nonproducer_expected_offspring']:.3f} → {on['median_nonproducer_expected_offspring']:.3f}. "
                  "Benefit formula is a teaching assumption, not fitted data.", fontsize=10)
    output = Path(output)
    save_figure(fig, output, "04_il11_shared_fitness")
    write_csv(output / "il11_comparison.csv", rows)
    write_csv(output / "il11_comparison_summary.csv", summaries)
    # 保存曲線資料，方便自行重畫或核對圖中數值。
    write_csv(output / "il11_benefit_curve.csv",
              [dict(il11_percent=100 * f, benefit_percentage_points=100 * b) for f, b in zip(fraction, benefit)])
    return rows


def write_csv(path: Path, rows: list[dict]):
    # 以第一筆資料的鍵建立欄名；UTF-8 支援中文，newline 避免 CSV 多餘空行。
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def plot_therapy_ensemble(trials: list[Result], output: Path):
    """依絕對世代彙整治療試驗；滅絕延續零值，計算上限後視為缺失。"""
    config = trials[0].config
    # 不把因計算上限停止的族群當成滅絕，也不以最後數量冒充後續觀察。
    length = max(len(t.history) for t in trials)
    keys = ("total", "sensitive", "resistant", "effective_genotypes", "resistant_percent")
    data = {k: np.full((len(trials), length), np.nan) for k in keys}
    raw = []
    for rep, trial in enumerate(trials):
        for generation in range(length):
            if generation >= len(trial.history) and trial.stop_reason != "extinction":
                continue
            row = trial.history[min(generation, len(trial.history) - 1)]
            total = row["total"]
            values = dict(total=total, sensitive=row["sensitive"], resistant=row["resistant"],
                          effective_genotypes=np.exp(row["shannon"]) if total else 0,
                          resistant_percent=100 * row["resistant"] / total if total else np.nan)
            raw.append(dict(replicate=rep, seed=trial.seed, generation=generation,
                            extended_extinction=generation >= len(trial.history), **values))
            for key, value in values.items():
                data[key][rep, generation] = value
    # 每代輸出平均、樣本標準差、中位數、百分位與有效 n；n=1 時 SD 不定義。
    summary = []
    for generation in range(length):
        record = dict(generation=generation)
        for key, matrix in data.items():
            v = matrix[:, generation]; v = v[np.isfinite(v)]
            record.update({f"{key}_n": len(v), f"{key}_mean": float(v.mean()) if len(v) else np.nan,
                           f"{key}_sd": float(v.std(ddof=1)) if len(v)>1 else np.nan})
            for name, q in (("p10", 10), ("median", 50), ("p90", 90)):
                record[f"{key}_{name}"] = float(np.percentile(v, q)) if len(v) else np.nan
        summary.append(record)
    # 三格：負荷、異質性、抗藥比例。逐基因頻率改由 plot_genotype_evolution.py 的基因頻率圖呈現；
    # 可用試驗數只是記帳資訊，改寫進註腳。
    set_plot_style()
    fig, axes = plt.subplots(1, 3, figsize=(18, 6.2), layout="constrained")
    fig.suptitle(f"Tumor evolution under continuous therapy | {len(trials)} stochastic trials",
                 fontsize=21, fontweight="bold")
    x = np.arange(length)
    starts = [t.therapy_generation for t in trials if t.therapy_generation is not None]
    for key, color in zip(keys[:3], ("#172B42", COLORS[0], COLORS[1])):
        axes[0].plot(x, [r[f"{key}_median"] for r in summary], label=key.title(), color=color)
        axes[0].fill_between(x, [r[f"{key}_p10"] for r in summary],
                             [r[f"{key}_p90"] for r in summary], color=color, alpha=0.12)
    axes[0].set(title="A  Population response", ylabel="Cells (symlog)")
    axes[0].set_yscale("symlog", linthresh=1)
    for ax, key, title, ylabel in ((axes[1], "effective_genotypes", "B  Genetic heterogeneity", "Effective genotypes: exp(Shannon H)"),
                                  (axes[2], "resistant_percent", "C  Resistance in surviving tumors", "Resistant cells (%)")):
        ax.plot(x, [r[f"{key}_median"] for r in summary], color=COLORS[2])
        ax.fill_between(x, [r[f"{key}_p10"] for r in summary], [r[f"{key}_p90"] for r in summary], color=COLORS[2], alpha=0.18)
        ax.set(title=title, ylabel=ylabel)
    # 治療起點因試驗而異，以淡色區間標示最早到最晚的啟動世代。
    if starts:
        for ax in axes:
            ax.axvspan(min(starts), max(starts), color=COLORS[1], alpha=0.10, lw=0,
                       label=f"Therapy starts (gen {min(starts)}–{max(starts)})")
        # 治療前已存在的抗藥細胞：在各試驗自己的治療起點取值。
        at_start = [100 * t.history[t.therapy_generation]["resistant"] / t.history[t.therapy_generation]["total"]
                    for t in trials if t.therapy_generation is not None]
        present = sum(t.history[t.therapy_generation]["resistant"] > 0 for t in trials if t.therapy_generation is not None)
        axes[2].annotate(f"Already present at therapy start:\nmedian {np.median(at_start):.1f}% resistant\n"
                         f"({present}/{len(at_start)} trials > 0)",
                         xy=(np.median(starts), np.median(at_start)), xytext=(0.36, 0.42), textcoords="axes fraction",
                         arrowprops=dict(arrowstyle="->", color=COLORS[1]), fontsize=10.5, color=COLORS[1])
    axes[0].legend(frameon=False, fontsize=9)
    for ax in axes:
        ax.set_xlabel("Generation"); ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    # 可用試驗數改成文字：全部試驗持續觀察到哪一代，之後因計算上限而截尾。
    complete = next((g for g, r in enumerate(summary) if r["total_n"] < len(trials)), length)
    fig.supxlabel("Lines: median; bands: 10–90% of trials, not confidence intervals. Absolute generation alignment.\n"
                  f"All {len(trials)} trials observed through generation {complete - 1}; later generations include only trials "
                  f"still below the {config.max_cells:,}-cell computational cap (not tumor death).", fontsize=10)
    save_figure(fig, output, "01_tumor_dashboard")
    write_csv(output / "trajectory.csv", summary)
    write_csv(output / "therapy_trajectories.csv", raw)


def main():
    # 命令列入口：允許指定亂數種子、重複次數及輸出資料夾。
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--suppressor-hits", type=int, choices=(1, 2), default=2,
                        help="TSG allele-loss threshold for the main treatment scenario")
    parser.add_argument("--il11", action="store_true", help="Enable IL11-mediated shared fitness advantage")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--replicates", type=int, default=30)
    parser.add_argument("--output", type=Path, default=Path("outputs"))
    # 解析並檢查命令列輸入，建立後續輸出的資料夾。
    args = parser.parse_args()
    if args.replicates < 1:
        parser.error("--replicates must be >= 1")
    args.output.mkdir(parents=True, exist_ok=True)
    # 執行預設治療案例，輸出四面板圖與逐代統計，再執行情境比較。
    config = Config(suppressor_hits=args.suppressor_hits, il11_enabled=args.il11)
    trials = [simulate(config, args.seed + 100000 + rep) for rep in range(args.replicates)]
    plot_therapy_ensemble(trials, args.output)
    run_comparisons(config, args.replicates, args.seed, args.output)
    run_il11_comparison(config, args.replicates, args.seed, args.output)
    # 抗藥性為隨機事件，額外報告重複試驗結果，避免只展示成功復發案例。
    outcomes = []
    for rep, trial in enumerate(trials):
        outcomes.append(dict(replicate=rep, seed=trial.seed,
                             therapy_generation=trial.therapy_generation,
                             stop_reason=trial.stop_reason, **therapy_outcome(trial)))
    write_csv(args.output / "therapy_replicates.csv", outcomes)
    # 保存完整參數、種子、套件版本與主要結果，方便追溯及重現本次執行。
    metadata = dict(config=asdict(config), seed=args.seed, replicates=args.replicates,
                    therapy_seeds=[t.seed for t in trials],
                    therapy_outcome_counts={status: sum(o["status"] == status for o in outcomes)
                                           for status in sorted({o["status"] for o in outcomes})},
                    interval="median and 10–90% trial percentiles; not confidence intervals",
                    numpy_version=np.__version__,
                    matplotlib_version=matplotlib.__version__)
    (args.output / "run_metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    # 在終端顯示結果摘要和輸出位置，方便找到生成的圖表與資料。
    print(json.dumps(metadata, indent=2))
    print(f"Figures and data saved to: {args.output.resolve()}")


# 只在直接執行本檔案時啟動 main；被其他程式 import 時不自動跑模擬。
if __name__ == "__main__":
    main()
