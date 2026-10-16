"""腫瘤分枝模型：只包含模擬本身。

一代的規則：每個細胞死亡，或產生 2 / 3 個子細胞；子細胞每個位元獨立以 μ 發生不可逆突變。
腫瘤以「基因型 → 細胞數」的字典表示（不是追蹤每個細胞的家系）。
治療在負荷首次達到門檻後開始，敏感與抗藥細胞有不同的額外殺傷率。
所有參數均為教學假設，不是臨床估計。

閱讀順序：Config（參數）→ 基因型與表型 → growth_probabilities（表型 → 機率）
→ mutate_offspring（突變）→ simulate（一代一代更新）→ therapy_outcome（判讀復發）。
重複試驗、對照實驗與分析在 experiments.py；理論公式在 analysis.py；圖在 plot_*.py。

Tumor branching model: the simulation only.

Rule for one generation: each cell dies or makes 2 / 3 daughters; each daughter bit mutates irreversibly with probability μ.
A tumor is a dictionary genotype → cell count (lineages of individual cells are not tracked).
Therapy starts when the burden first reaches a threshold; sensitive and resistant cells get different extra kill rates.
All parameters are teaching assumptions, not clinical estimates.

Reading order: Config (parameters) → genotype and phenotype → growth_probabilities (phenotype → probabilities)
→ mutate_offspring (mutation) → simulate (generation-by-generation update) → therapy_outcome (relapse classification).
Repeated trials, controls and analyses are in experiments.py; theory formulas in analysis.py; figures in plot_*.py.
"""

# 延後解析型別註記，讓函數與資料類別的型別宣告更容易互相引用。
# Postpone evaluation of type annotations so functions and dataclasses can refer to each other's types.
from __future__ import annotations

from dataclasses import dataclass, field

# NumPy 負責隨機抽樣與統計摘要。
# NumPy handles random sampling and summary statistics.
import numpy as np


# ---------------------------------------------------------------- 參數與結果 / Parameters and results --------
# 將所有模型參數集中管理；frozen=True 防止執行中意外修改設定。
# All model parameters in one place; frozen=True prevents accidental changes during a run.
@dataclass(frozen=True)
class Config:
    # 20 個基因座：一般基因用一個位元；每個抑癌基因另加一個等位基因位元。
    # 20 loci: one bit per regular gene; each tumor suppressor gets an extra bit for its second allele.
    n_genes: int = 20
    driver_genes: tuple[int, ...] = (0, 1, 2)
    resistance_gene: int = 3
    driver_hits: int = 1  # 保留：不同 oncogene 的突變數門檻 / kept: mutation-count threshold across oncogenes
    suppressor_genes: tuple[int, ...] = (4, 5)
    suppressor_hits: int = 2  # 同一抑癌基因需要失活的等位基因數 / number of alleles that must be lost in one tumor suppressor
    suppressor_death_reduction: float = 0.08
    suppressor_triple_bonus: float = 0.12
    # 突變與生長設定：driver 達標後同時影響自然死亡率和三子細胞機率。
    # Mutation and growth settings: a driver past its threshold changes both natural death and triple division.
    mutation_rate: float = 0.001  # 每個子細胞、每個未突變位元、每代；TSG 為每等位基因 / per daughter, per unmutated bit, per generation; per allele for TSGs
    death_rate: float = 0.30
    driver_death_reduction: float = 0.08
    triple_probability: float = 0.04  # 存活後產生三個子細胞的機率 / probability that a survivor makes three daughters
    driver_triple_bonus: float = 0.12
    # 模擬規模及治療設定：控制初始負荷、觀察時間、啟動門檻與額外殺傷率。
    # Scale and therapy settings: initial burden, observation time, start threshold and extra kill rates.
    initial_cells: int = 100
    generations: int = 75
    therapy_start_size: int = 20_000
    sensitive_therapy_kill: float = 0.90  # 自然死亡後的額外殺傷機率 / extra kill probability after natural death
    resistant_therapy_kill: float = 0.04
    max_cells: int = 1_000_000  # 計算停止門檻，不是環境承載量 / computational stopping cap, not a carrying capacity
    # 延伸分析用（預設關閉，主結果不受影響）：
    # For follow-up analyses (off by default; main results are unaffected):
    # adaptive：負荷降到首次治療負荷的 adaptive_pause_fraction 以下就暫停，回到治療門檻再恢復。
    # adaptive: pause when the burden falls below adaptive_pause_fraction of the first-treatment burden, resume at the threshold.
    therapy_mode: str = "continuous"
    adaptive_pause_fraction: float = 0.5
    resistance_cost: float = 0.0  # 抗藥細胞額外的自然死亡機率（抗藥代價） / extra natural death probability of resistant cells (cost of resistance)
    # 承載量 K：每個細胞額外死亡 c·N/K，c 使無突變細胞在 N=K 時平均子代剛好為 1。None＝無競爭。
    # Carrying capacity K: each cell gets extra death c·N/K, with c set so a mutation-free cell has mean offspring 1 at N = K. None = no competition.
    carrying_capacity: int | None = None

    # 驗證基因索引、機率與模擬設定，遇到無效參數時拋出 ValueError。
    # Validate gene indices, probabilities and simulation settings; raise ValueError on invalid parameters.
    def validate(self):
        # 先檢查基因組大小、基因索引是否有效，以及 driver 門檻是否可達。
        # Check genome size, valid gene indices, and that the driver threshold can be reached.
        if not 2 <= self.n_genes <= 20:
            raise ValueError("n_genes 必須介於 2–20")
        genes = self.driver_genes + self.suppressor_genes + (self.resistance_gene,)
        if len(set(genes)) != len(genes) or any(
            g < 0 or g >= self.n_genes for g in genes
        ):
            raise ValueError(
                "oncogene、suppressor 與 resistance 基因必須互異且位於基因組內"
            )
        if not 1 <= self.driver_hits <= len(self.driver_genes):
            raise ValueError("driver_hits 超出可用 driver 數目")
        if not self.suppressor_genes or self.suppressor_hits not in (1, 2):
            raise ValueError("需要至少一個抑癌基因，suppressor_hits 必須為 1 或 2")
        # 逐一檢查機率參數，避免傳入二項分布時出現無效機率。
        # Check each probability parameter so no invalid probability reaches the binomial draws.
        for name in (
            "mutation_rate",
            "death_rate",
            "driver_death_reduction",
            "triple_probability",
            "driver_triple_bonus",
            "suppressor_death_reduction",
            "suppressor_triple_bonus",
            "sensitive_therapy_kill",
            "resistant_therapy_kill",
            "resistance_cost",
        ):
            if not 0 <= getattr(self, name) <= 1:
                raise ValueError(f"{name} 必須介於 0–1")
        # 額外檢查加總後的機率與模擬規模；單個參數有效不代表組合必然有效。
        # Also check summed probabilities and simulation scale; valid single parameters do not guarantee a valid combination.
        if (
            self.triple_probability
            + self.driver_triple_bonus
            + self.suppressor_triple_bonus
            > 1
        ):
            raise ValueError("產生三個子細胞的總機率不能超過 1")
        if (
            min(
                self.initial_cells,
                self.generations,
                self.therapy_start_size,
                self.max_cells,
            )
            <= 0
        ):
            raise ValueError("細胞數、代數與門檻必須為正")
        if self.therapy_mode not in ("continuous", "adaptive"):
            raise ValueError("therapy_mode 必須是 continuous 或 adaptive")
        if not 0 < self.adaptive_pause_fraction < 1:
            raise ValueError("adaptive_pause_fraction 必須介於 0 與 1 之間")
        if self.carrying_capacity is not None and self.carrying_capacity <= 0:
            raise ValueError("carrying_capacity 必須為正或 None")


# 封裝一次試驗的設定、逐代統計、完整基因型計數及停止／治療資訊。
# One trial: settings, per-generation statistics, full genotype counts, and stop / therapy information.
@dataclass
class Result:
    config: Config
    seed: int
    history: list[dict]
    populations: list[dict[int, int]]
    therapy_generation: int | None
    stop_reason: str
    # 每個已推進世代是否施加治療；continuous 時從 therapy_generation 起全為 True。
    # Whether therapy was applied in each advanced generation; with continuous therapy all True from therapy_generation on.
    therapy_active: list[bool] = field(default_factory=list)


# ---------------------------------------------------------------- 基因型 → 表型 / Genotype → phenotype ------
# 計算指定基因型中已突變的 oncogene 數目。
# Number of mutated oncogenes in a genotype.
def driver_count(genotype: int, config: Config) -> int:
    # 1 << g 選取第 g 個基因；位元 AND 判斷是否突變，再加總 driver 數。
    # 1 << g selects gene g; bitwise AND tests whether it is mutated, then sum the drivers.
    return sum(bool(genotype & (1 << g)) for g in config.driver_genes)


# 回傳指定抑癌基因已失活的等位基因數（0、1 或 2）。
# Number of inactivated alleles of a given tumor suppressor (0, 1 or 2).
def suppressor_alleles(genotype: int, gene: int, config: Config) -> int:
    """回傳同一抑癌基因失活的等位基因數（0、1 或 2）。

    Return the number of inactivated alleles of one tumor suppressor (0, 1 or 2).
    """
    # 第一個位元沿用基因座索引；第二個放在 n_genes 後，不增加基因座數。
    # The first allele uses the locus index; the second sits after n_genes, adding no locus.
    second = config.n_genes + config.suppressor_genes.index(gene)
    return int(bool(genotype & (1 << gene))) + int(bool(genotype & (1 << second)))


# 判斷是否至少有一個抑癌基因達到設定的失活門檻。
# Whether at least one tumor suppressor reaches the inactivation threshold.
def suppressor_loss(genotype: int, config: Config) -> bool:
    # 必須是同一基因達到門檻；兩個不同基因各一擊不等於同一基因兩擊。
    # The same gene must reach the threshold; one hit in each of two genes is not two hits in one gene.
    return any(
        suppressor_alleles(genotype, gene, config) >= config.suppressor_hits
        for gene in config.suppressor_genes
    )


# 依 oncogene、抑癌基因與抗藥代價，回傳自然死亡率及存活後產生三個子細胞的機率。
# Natural death rate and triple-division probability of survivors, from oncogenes, tumor suppressors and resistance cost.
def growth_probabilities(genotype: int, config: Config) -> tuple[float, float]:
    """分開計算 oncogene 活化及抑癌功能喪失的生長優勢，可相加。

    Growth advantages from oncogene activation and tumor-suppressor loss, computed separately; they add.
    """
    oncogene = driver_count(genotype, config) >= config.driver_hits
    suppressor = suppressor_loss(genotype, config)
    # 每一類優勢最多給一次；多個抑癌基因達標不再重複累加。
    # Each class of advantage counts at most once; several TSGs reaching threshold do not stack.
    death = max(
        0,
        config.death_rate
        - config.driver_death_reduction * oncogene
        - config.suppressor_death_reduction * suppressor,
    )
    # 抗藥代價只在明確設定時加上，預設為 0 時完全不改變計算。
    # The resistance cost is added only when set; the default 0 leaves the calculation unchanged.
    if config.resistance_cost and genotype & (1 << config.resistance_gene):
        death = min(1.0, death + config.resistance_cost)
    triple = (
        config.triple_probability
        + config.driver_triple_bonus * oncogene
        + config.suppressor_triple_bonus * suppressor
    )
    return death, triple


# 依總細胞數與承載量計算競爭造成的額外死亡機率，未設定承載量時回傳零。
# Extra death probability from competition, from total cells and carrying capacity; zero when no capacity is set.
def crowding_death(total: int, config: Config) -> float:
    """承載量競爭的額外死亡機率；係數讓無突變細胞在 N = K 時期望子代恰為 1。

    Extra death from carrying-capacity competition; the coefficient gives a mutation-free cell mean offspring exactly 1 at N = K.
    """
    if config.carrying_capacity is None:
        return 0.0
    coefficient = (1 - config.death_rate) - 1 / (2 + config.triple_probability)
    return coefficient * total / config.carrying_capacity


# ---------------- 突變與一代更新 ----
# 逐位元抽樣不可逆突變並合併相同基因型，保持子代細胞總數不變。
# Sample irreversible mutations bit by bit and merge identical genotypes, keeping the total offspring count.
def mutate_offspring(
    population: dict[int, int], config: Config, rng: np.random.Generator
) -> dict[int, int]:
    """逐基因二項抽樣：與逐細胞獨立突變等價，允許一代多基因突變。

    Binomial sampling per gene: equivalent to independent per-cell mutation; several genes may mutate in one generation.
    """
    # 複製輸入字典，避免突變抽樣直接改動呼叫端保存的族群。
    # Copy the input dictionary so sampling does not modify the caller's population.
    population = dict(population)
    # 依序處理一般基因與 TSG 的兩個等位基因；允許同一代兩個等位基因都突變。
    # Process regular genes and both TSG alleles in turn; both alleles may mutate in the same generation.
    for gene in range(config.n_genes + len(config.suppressor_genes)):
        bit = 1 << gene
        # 將讀取與寫入分開，避免同一輪新增的基因型被重複處理。
        # Separate reads from writes so genotypes created in this round are not processed twice.
        updated = dict(population)
        for genotype, count in population.items():
            # 已突變的基因保持不變；只替尚未突變的細胞抽樣。
            # Mutated genes stay mutated; sample only cells not yet mutated.
            if genotype & bit:
                continue
            # 二項分布一次抽出此基因型的突變細胞數，省去逐細胞迴圈。
            # One binomial draw gives the number of mutant cells for this genotype, with no per-cell loop.
            mutants = int(rng.binomial(count, config.mutation_rate))
            if mutants:
                # 從舊基因型移出突變細胞，以位元 OR 設定新突變並合併相同基因型。
                # Move mutants out of the old genotype, set the new bit with OR, and merge identical genotypes.
                updated[genotype] -= mutants
                updated[genotype | bit] = updated.get(genotype | bit, 0) + mutants
        # 刪除零細胞的項目；突變只重新分類細胞，不改變細胞總數。
        # Drop empty entries; mutation only relabels cells and keeps the total.
        population = {g: n for g, n in updated.items() if n > 0}
    return population


# 彙整單代族群的負荷、抗藥組成、多樣性與功能性突變比例。
# Summarize one generation: burden, resistant composition, diversity and functional-mutation fractions.
def summarize(population: dict[int, int], generation: int, config: Config) -> dict:
    # 計算總負荷及各基因型比例；分母至少為 1，讓滅絕族群也能安全處理。
    # Total burden and genotype fractions; the denominator is at least 1 so extinct populations are safe.
    total = sum(population.values())
    fractions = np.array(list(population.values()), dtype=float) / max(total, 1)
    # 抗藥基因位元決定敏感／抗藥分類，與是否已達 driver 門檻分開判斷。
    # The resistance bit decides sensitive / resistant, separately from whether a driver threshold is reached.
    resistant = sum(
        n for g, n in population.items() if g & (1 << config.resistance_gene)
    )
    # 彙整多樣性：richness 是基因型數；Shannon 用自然對數；Simpson 為 1−Σp²。
    # Diversity: richness is the number of genotypes; Shannon uses the natural log; Simpson is 1−Σp².
    # driver_fraction 只計 oncogene 達標；suppressor_fraction 計 TSG 達標。
    # driver_fraction counts oncogenes past threshold only; suppressor_fraction counts TSGs past threshold.
    # suppressor_biallelic_fraction 計至少一個 TSG 雙等位基因失活；滅絕時皆為零。
    # suppressor_biallelic_fraction counts at least one TSG with both alleles lost; all zero at extinction.
    return dict(
        generation=generation,
        total=total,
        resistant=resistant,
        sensitive=total - resistant,
        richness=len(population),
        shannon=float(-np.sum(fractions * np.log(fractions))),
        simpson=float(1 - np.sum(fractions**2)) if total else 0.0,
        driver_fraction=sum(
            n
            for g, n in population.items()
            if driver_count(g, config) >= config.driver_hits
        )
        / max(total, 1),
        suppressor_fraction=sum(
            n for g, n in population.items() if suppressor_loss(g, config)
        )
        / max(total, 1),
        suppressor_biallelic_fraction=sum(
            n
            for g, n in population.items()
            if any(
                suppressor_alleles(g, gene, config) == 2
                for gene in config.suppressor_genes
            )
        )
        / max(total, 1),
    )


# 執行含突變與可選持續或適應性治療的逐代分枝模擬，回傳族群快照、摘要及停止資訊。
# Run the generation-by-generation branching simulation with mutation and optional continuous or adaptive therapy; return snapshots, summaries and stop information.
def simulate(
    config: Config = Config(),
    seed: int = 42,
    therapy: bool = True,
    target_size: int | None = None,
) -> Result:
    """同步世代分枝過程；治療從首次達門檻後的轉移開始，持續至模擬結束。

    超過 max_cells 時保留完整該代，不截斷細胞；target_size 用於等大小比較。

    Synchronous-generation branching process; therapy starts at the first transition after the threshold and lasts to the end.

    Above max_cells the whole generation is kept, not truncated; target_size is used for equal-size comparisons.
    """
    # 先驗證設定與目標大小，再建立一次模擬專用的亂數產生器。
    # Validate the settings and target size, then create a random generator for this run.
    config.validate()
    if target_size is not None and target_size <= 0:
        raise ValueError("target_size 必須為正")
    # 固定 seed 可重現同一試驗；基因型 0 代表所有基因都尚未突變。
    # A fixed seed reproduces the trial; genotype 0 means no gene is mutated yet.
    rng = np.random.default_rng(seed)
    population = {0: config.initial_cells}
    # history 存摘要，populations 存基因型快照；尚未治療時代數記為 None。
    # history stores summaries, populations stores genotype snapshots; the therapy generation is None until therapy starts.
    history, populations = [], []
    therapy_generation = None
    therapy_on = False
    active = []
    stop_reason = "generation_limit"
    # 包含第 0 代初始狀態；每輪先記錄當前族群，再決定是否推進下一代。
    # Includes generation 0; each round records the current population, then decides whether to advance.
    for generation in range(config.generations + 1):
        row = summarize(population, generation, config)
        history.append(row)
        populations.append(dict(population))
        # 依序處理滅絕、達取樣目標、達計算上限與觀察時間結束，保留停止原因。
        # Check, in order, extinction, sampling target, computational cap and end of observation; keep the stop reason.
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
        # Record the generation when the therapy threshold is first reached; therapy applies to all later transitions.
        if (
            therapy
            and therapy_generation is None
            and row["total"] >= config.therapy_start_size
        ):
            therapy_generation = generation
            therapy_on = True
        # adaptive：降到首次治療負荷的一定比例以下就暫停，回到治療門檻再恢復。
        # adaptive: pause below a fraction of the first-treatment burden, resume at the threshold.
        elif therapy_generation is not None and config.therapy_mode == "adaptive":
            reference = history[therapy_generation]["total"]
            if (
                therapy_on
                and row["total"] <= config.adaptive_pause_fraction * reference
            ):
                therapy_on = False
            elif not therapy_on and row["total"] >= config.therapy_start_size:
                therapy_on = True
        active.append(therapy_on)
        crowding = crowding_death(row["total"], config)
        # 依基因型批次處理所有父細胞；新世代完全由 offspring 取代。
        # Process all parent cells in batches by genotype; the new generation is entirely offspring.
        offspring = {}
        for genotype, count in population.items():
            # 將 oncogene 與抑癌基因的效果合併，再套用治療殺傷。
            # Combine oncogene and tumor-suppressor effects, then apply therapy kill.
            death, triple = growth_probabilities(genotype, config)
            if crowding:
                death = min(1.0, death + crowding)
            # 治療殺傷作用於自然存活者，因此用兩個存活機率相乘計算總死亡率。
            # Therapy kills natural survivors, so total death uses the product of the two survival probabilities.
            if therapy_on:
                if genotype & (1 << config.resistance_gene):
                    kill = config.resistant_therapy_kill
                else:
                    kill = config.sensitive_therapy_kill
                death = 1 - (1 - death) * (1 - kill)
            # 先抽出存活父細胞，再抽出其中產生三個子細胞者；其餘各產生兩個。
            # Draw surviving parents, then those among them that make three daughters; the rest make two.
            survivors = int(rng.binomial(count, 1 - death))
            triples = int(rng.binomial(survivors, triple))
            # 每個存活者先貢獻兩個子細胞，三子細胞事件再各加一個。
            # Each survivor contributes two daughters; each triple event adds one more.
            n = 2 * survivors + triples
            if n:
                offspring[genotype] = n
        # 所有分裂完成後才讓子代突變；新突變的生存效果從下一次轉移開始。
        # Mutate offspring only after all divisions; new mutations affect survival from the next transition.
        population = mutate_offspring(offspring, config, rng)
    # 把逐代資料與停止資訊一起回傳，供繪圖和治療結果判讀使用。
    # Return per-generation data with stop information for plotting and outcome classification.
    return Result(
        config, seed, history, populations, therapy_generation, stop_reason, active
    )


# ---------------------------------------------------------------- 判讀治療結果 / Classify the therapy outcome ------
# 依治療後負荷是否先下降再回升至起始負荷，判定治療結果並回傳谷底與首次復發世代。
# Classify the outcome by whether the burden first falls and then returns to the starting burden; return the nadir and first relapse generation.
def therapy_outcome(result: Result) -> dict:
    """復發定義：先降至治療前負荷以下，再回到治療開始時負荷。

    Relapse: first fall below the pre-therapy burden, then return to the burden at therapy start.
    """
    # 未啟動治療的試驗單獨標記，不當作治療成功或失敗。
    # Trials where therapy never started are marked separately, not counted as success or failure.
    start = result.therapy_generation
    if start is None:
        return dict(status="therapy_not_started", nadir=None, relapse_generation=None)
    # 以啟動治療時的負荷為基準，只搜尋治療開始後的世代。
    # Use the burden at therapy start as the baseline and search only generations after it.
    baseline = result.history[start]["total"]
    after = result.history[start + 1 :]
    declined = False
    relapse = None
    # 記住是否曾低於基準；一旦之後回升至基準，就記錄第一次復發代數。
    # Remember whether the burden dipped below baseline; the first return to baseline afterwards is the relapse generation.
    for row in after:
        declined |= row["total"] < baseline
        if declined and row["total"] >= baseline:
            relapse = row["generation"]
            break
    # 區分復發、滅絕、尚未觀察到復發及無初始反應；nadir 是觀察期最低負荷。
    # Distinguish relapse, extinction, no relapse yet and no initial response; nadir is the lowest burden observed.
    # fmt: off
    status = (
        "relapse" if relapse is not None
        else "extinction" if result.stop_reason == "extinction"
        else "no_relapse_observed" if declined
        else "no_initial_response"
    )
    # fmt: on
    return dict(
        status=status,
        nadir=min((r["total"] for r in after), default=baseline),
        relapse_generation=relapse,
    )
