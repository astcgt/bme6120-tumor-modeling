# BME 6120 — Quantitative tumor modeling

[English](#english) | [繁體中文](#繁體中文)

## English

A stochastic discrete branching model for classroom discussion. Parameters are teaching assumptions, not calibrated clinical estimates. Code documentation and comments are in Chinese. Core simulation figures use English labels; the standalone teaching diagrams and phase portraits use Traditional Chinese labels.

### Run without modifying conda

Use the project virtual environment directly; activation is unnecessary:

```bash
.venv/bin/python main.py                                   # every step (about 3 minutes on 12 cores)
.venv/bin/python main.py --only results                    # main and control experiments only
.venv/bin/python main.py --only results --suppressor-hits 1
.venv/bin/python main.py --only therapy --analyses 08 09   # part of the therapy analyses
```

Steps run in order: `results` (figures 01–03), `model` (06, 07, A1), `phase` (phase portraits), `genotype` (Muller plot and gene frequencies), `therapy` (08–11). Later steps read the CSVs written by `results`, so run `results` first on a fresh checkout. Outputs always go to `outputs/`; `--seed` and `--replicates` change the main and control experiments.

### Code structure

| File | Contents |
| --- | --- |
| `simulation.py` | The model only: `Config`, genotype → phenotype, `growth_probabilities`, `mutate_offspring`, `simulate`, `therapy_outcome`. Start here when explaining the code. |
| `experiments.py` | Runs the model under different conditions (main experiment, TSG and death-rate controls, therapy analyses), summarizes trials, and reads/writes data files. |
| `analysis.py` | Formulas derived from the model rules: expected offspring and one-step drift, branching-process extinction, Goldie–Coldman theory, phase-plane vector fields, Muller layout and gene frequencies. |
| `plot_style.py`, `plot_results.py`, `plot_model.py`, `plot_phase.py`, `plot_genotype.py`, `plot_therapy.py` | Figures only, grouped by topic; each takes data computed above. |
| `main.py` | Entry point: `main()` calls one short function per step. |
| `test_simulation.py` | Unit tests for the model and the ensemble statistics. |

On a new computer, first create an isolated environment and install dependencies:

```bash
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
```

In VS Code, select this project's `.venv/bin/python` as the Python interpreter. Packages reside in `.venv`; the default Matplotlib cache is the project's `.mplconfig` directory.

### Model assumptions

- Start with 100 unmutated cells and 20 gene loci: 3 oncogenes (Python indices 0–2), 1 resistance gene (3), 2 tumor suppressor genes / TSGs (4–5), and 14 passengers (6–19). Each TSG has two independently mutable alleles; other loci use a binary mutation state. The genotype has 22 bits representing 20 genes, not 22 genes.
- Each generation replaces each parent cell with 0, 2, or 3 offspring; the parent is not retained. With death probability `d` and probability `q` of three offspring conditional on survival, expected offspring per parent are `(1-d)*(2+q)`.
- Oncogene activation and TSG functional loss each reduce death probability and increase the probability of three offspring; these two effects add, once per category. `driver_hits` retains the separate threshold for the number of mutated oncogenes. `suppressor_hits=1` grants a TSG advantage after either allele is lost; `suppressor_hits=2` (default) requires loss of both alleles of the same TSG. One hit in each of two different TSGs does not satisfy the two-hit rule. The one-hit case is a simplified haploinsufficiency scenario, not a claim about every TSG.
- Offspring inherit mutations, and each unmutated bit independently mutates with probability `mutation_rate`. Mutations are irreversible; multiple genes may mutate in one generation. The default rate of 0.001 per bit per offspring (per allele for TSGs) is chosen for a short demonstration. Both TSG alleles can mutate in one generation. Allele loss is represented by irreversible inactivating mutations; LOH, deletions, and epigenetic silencing are not separately modeled.
- Grouping cells by genotype and sampling binomial counts is equivalent to individual-cell sampling under these assumptions, but more efficient. Identical genotypes may have different ancestors, so a genotype is not necessarily a lineage clone.
- When the tumor first reaches 20,000 cells, continuous treatment starts with the next generation transition. After natural survival, sensitive cells have a 90% additional probability of death and resistant cells have a 4% probability; combined death probability is `1-(1-d)*(1-kill)`.
- Resistance mutations can arise before or during treatment; treatment does not directly increase the mutation rate. New mutations arise in offspring and affect survival from the next generation.
- Relapse means the population first falls below its treatment-start burden and then returns to that burden. No observed relapse does not imply permanent cure.
- `max_cells` is a computational stopping threshold. The full generation that crosses it is retained; it is not a carrying capacity or population truncation. The model excludes spatial structure, immunity, resource competition, and quiescent cells.

### Figures and data

The `outputs/` directory contains:

| File | Purpose |
| --- | --- |
| `01_tumor_dashboard.png / .svg` | Ensemble burden, effective genotypes and resistant fraction under therapy; resistance present at therapy start is annotated |
| `02_death_rate_heterogeneity.png / .svg` | Shannon diversity of 5,000-cell samples and generations needed to reach 5,000 cells at three death rates |
| `03_one_hit_vs_two_hit.png / .svg` | No-TSG-advantage control, one-hit and two-hit: linear growth, TSG-loss fraction, and trial burdens at generation 18 |
| `one_hit_vs_two_hit.csv / one_hit_vs_two_hit_summary.csv` | Per-generation comparison data and endpoint medians |
| `trajectory.csv` | Ensemble statistics per generation: n, mean, sample SD, median, p10, p90 |
| `therapy_trajectories.csv` | Individual treatment trajectories used in the ensemble |
| `comparisons.csv` | Replicate data, including trials that fail to reach the target size |
| `therapy_replicates.csv` | Minimum burden, relapse generation, and stopping reason across treatment trials |
| `run_metadata.json` | Parameters, random seed, and package versions |

PNG figures use 240 dpi; SVG figures scale without loss of quality. Growth curves use symlog to display zero, with an approximately logarithmic scale above 1. `suppressor_fraction` records cells meeting the selected TSG threshold; `suppressor_biallelic_fraction` records cells with both alleles lost at at least one TSG. `driver_fraction` continues to refer only to the oncogene threshold.

Growth comparisons keep oncogene effects, mutation rates, and loci identical, varying only the TSG advantage: disabled (control), one-hit, or two-hit; all three appear in figure 03. By default, growth comparisons show the median and 10th–90th percentile interval across 30 trials. This describes variation across trials, not a confidence interval. Death-rate comparisons sample 5,000 cells without replacement from the first generation to cross 5,000 cells and calculate Shannon H; horizontal bars show medians. These are equal-sized **samples**, conditional on reaching the threshold, rather than whole tumors at an exactly matched size. Figure 02 panel B plots the generation at which each trial first reached 5,000 cells (medians 8 / 14 / 26.5 at death 0.20 / 0.35 / 0.45): more generations mean more divisions and more accumulated mutations at the same size. Trials that do not reach the threshold are excluded but retained in the CSV; the plot reports how many trials reached it.

### Direct one-hit versus two-hit comparison

![One-hit versus two-hit comparison](outputs/03_one_hit_vs_two_hit.png)

The figure shows the same 30 trials per scenario from the growth experiment, including the no-TSG-advantage control (dotted gray; TSG loss is counted with the two-hit rule but gives no advantage); no therapy is applied. At generation 18 the median burdens are 76,446 (control), 94,618 (one-hit) and 73,879 (two-hit), and 0.05% / 27.9% / 0.07% of cells meet the TSG threshold. Panel A uses a linear cell-count axis, panel B shows the percentage meeting each scenario's TSG-loss threshold (symlog, linear below 0.01%), and panel C shows every trial at generation 18 with median bars. The threshold differs by definition; panel B measures the resulting phenotype, not an identical mutation state in both groups. Bands show trial variation, not confidence intervals. Overlap and stochastic reversals are retained rather than selecting favorable seeds.

### Replication and uncertainty

`--replicates` controls **all stochastic experiment groups**, including the main treatment dashboard; default 30. For example, `.venv/bin/python main.py --only results --replicates 100`. The main dashboard now uses ensemble statistics instead of one trajectory. Population counts and diversity include extinct trials as zero; after a computational cell limit, later values are missing, never padded with the last population. The footnote states the last generation at which every trial is observed (generation 31 by default). Resistant proportions condition on observed living tumors; empty populations have no defined fraction. Per-gene frequencies are shown in `outputs/genotype_evolution/02_gene_mutation_frequencies_{zh,en}`.

Trials are aligned by absolute generation; treatment starts separately at each trial's size threshold, so the dashboard shades the range of start generations instead of a single line, and panel C annotates the median resistant fraction at each trial's own therapy start. `therapy_replicates.csv` retains each start time and outcome. `trajectory.csv` is now an ensemble summary; `therapy_trajectories.csv` stores individual trajectories. The standalone `plot_dashboard(result, ...)` remains available for explicitly requested single-trial exploration.

Growth comparisons show medians and 10–90% trial intervals; endpoint and diversity scatter plots show individual trials and medians. More replicates stabilize summaries but do not remove biological stochasticity; percentile bands describe trial variability, **not confidence intervals or standard errors**.

### Phase portraits

`.venv/bin/python main.py --only phase` creates one 300 dpi PNG/vector SVG figure, `outputs/phase_plots/01_phase_portrait`, without modifying the presentation. It shows sensitive vs resistant counts in `log10(count+1)` coordinates. Gray curves retain every observed trial; one actual trial is highlighted, selected first by proximity to the median relapse generation, then to the median nadir, and colored by pre-therapy, therapy-to-relapse and post-relapse. Treatment onset, minimum burden and relapse are marked. Black arrows end at the one-generation conditional mean calculated from the highlighted trial's full genotype snapshot, including the exact S-to-R mutation contribution; the distance from each arrowhead to the next observed point is that generation's stochastic deviation. A side column lists expected and observed next states at therapy start and at the nadir. The reproduced trial is checked against its recorded CSV counts. Arrows end at the log-transformed raw-count expectation, which differs from the expectation of log-transformed counts; the formulas are in figure 07.

This is a projection of the full genotype/therapy state, not an autonomous two-dimensional vector field: the same S/R counts can have different driver or TSG composition. Computational endpoints are not equilibria. `conditional_phase_steps.csv` retains the computed means and observed next states; `phase_summary.json` records selection, events and interpretation limits. Use `.venv/bin/python main.py --only phase --replicate 0` to highlight another actual trial (this overwrites the phase-portrait files). Chinese labels require a supported font such as Arial Unicode MS or Noto Sans CJK TC; SVG files preserve glyph paths. The burden × resistant-fraction view is in `05_phase_field_burden_logit` below, which uses a logit axis instead of a linear percentage.

The same step computes the grid vector-field data whose formulas figure 07 uses: `phase_vector_grid.csv` and `phase_vector_assumptions.json`. Each grid point specifies a full parent population with genotype 0 for S and only the resistance bit for R, and the one-step projected mean is checked against the full model at all 1,604 points (linear and log grids). Other offspring mutations remain enabled and can move later states off this slice. Treatment is fixed within each panel rather than switched at the size threshold. The former standalone linear and log vector-field figures were merged into figure 07 or removed.

### Theory versus empirical phase fields

The phase step also draws two coordinate pairs, each shown as theory versus empirical drift with therapy off and on:

- `05_phase_field_burden_logit_{zh,en}`: x = log10(N+1), y = resistant fraction on an empirical-logit axis, log10((R+½)/(S+½)), labeled in percent. It stays finite when ρ is 0 or 1, and the dashed and dotted lines mark S = 0 and R = 0.
- `06_phase_field_fitness_logit_{zh,en}`: x = mean no-therapy fitness w̄_off, y = the same logit axis. It separates evolution (oncogene accumulation) from selection (resistance).

Theory arrows use a fixed hidden-genotype slice: sensitive cells are genotype 0 and resistant cells carry only the resistance bit. For the fitness plane, the oncogene fraction is derived from w̄, independent of resistance, at N equal to the therapy threshold. Arrows end at the transformed one-step raw-count expectation, which is checked against `projected_offspring_mean`. Empirical arrows average the observed one-step displacement of all transitions starting in a cell, drawn only when a cell has ≥5 transitions from ≥2 trials; color shows the transition count. Within each therapy column, theory and empirical arrows share one length scale with separate x and y keys. Under therapy, theoretical logit ρ rises about 0.98 decades per generation: ρ → 1 is a stable compositional direction, but log N has no fixed point. Differences between rows show the effect of projected-out degrees of freedom and stochasticity. `phase_field_comparison.csv/.json` record every arrow, bin edge and assumption.

### Therapy analyses

`.venv/bin/python main.py --only therapy` (about 2 minutes on 12 cores; `--analyses 08 09` runs a subset) writes four figures with CSV/JSON data to `outputs/therapy_analysis/`. Analysis 09 reuses the 30 recorded main-experiment trials; the others use separate seed blocks. All parameters are teaching assumptions.

- `08_goldie_coldman`: cure probability versus burden at the start of continuous therapy (100–5,000 cells, 200 trials each, Wilson 95% intervals), against a branching-process theory with no fitted parameters. Resistant mutant clones arise as a Poisson process; a clone born s generations before therapy survives with 1 − f^(s)(q_on), where f is the untreated offspring generating function and q_on = 0.479 is the extinction probability of one resistant lineage under therapy. Cure ≈ exp(−Λ). Theory and simulation agree; 50% cure falls near 416 cells (99.5% at 100, 10% at 1,000, 0% at 5,000). Panel B shows the resistant cells already present at therapy start, whose mean exceeds the median (Luria–Delbrück jackpots).
- `09_relapse_time_prediction`: each main trial's relapse delay predicted from its genotypes at therapy start. A baseline prediction ln(N/R₀)/ln 1.371 is off by 2.0 generations on average; letting each resistant genotype grow at its own on-therapy mean (oncogene cells ×1.62) predicts the exact generation in 24/30 trials (mean error 0.2). Relapse is earlier with more resistant cells (Spearman −0.70) and with a larger oncogene share among them (−0.51). Panel C shows expected on-therapy growth versus resistant fraction ρ; it crosses zero at ρ* = (1 − m_S)/(m_R − m_S) ≈ 70%, the nadir.
- `10_parameter_sensitivity`: cure probability and median relapse delay over resistant-cell kill k_R (0–0.6) × mutation rate μ (10⁻⁴–3×10⁻³), 30 trials per cell. Plain resistant cells stop growing at k_R = 0.30, but cure requires about 0.41, where oncogene-carrying resistant cells also stop growing; higher μ makes such cells more likely to pre-exist.
- `11_adaptive_therapy`: intermittent (pause at 50% of the treatment-start burden, resume at 20,000) versus continuous therapy, 30 paired trials per setting. Progression = burden ≥ 24,000 while on drug, after first dropping below 20,000; the on-drug condition prevents a single off-drug growth step from counting as treatment failure. Without competition adaptive therapy is worse (median ratio ×0.85; ×0.91 with a resistance cost). With carrying-capacity competition (K = 30,000) it delays progression ×1.14 (14 → 16 generations, longer in 30/30 pairs; 17.5 → 20 with cost) while using drug about half the time.

The adaptive analysis uses three `Config` options that are off by default and leave all other outputs unchanged: `therapy_mode="adaptive"` with `adaptive_pause_fraction`, `resistance_cost` (extra natural death for resistant cells), and `carrying_capacity` (extra death 0.21·N/K per cell, so a mutation-free cell has mean offspring 1 at N = K). The PCA coordinate analysis used in an earlier version was removed because its leading component only tracked time and size.

### Genotype distribution over time

`.venv/bin/python main.py --only genotype` writes to `outputs/genotype_evolution/`:

- `01_muller_top_genotypes_{zh,en}`: a Muller plot of the representative trial (selected by the same rule as the phase portraits) beneath its burden curve, plus six other trials spanning early to late relapse. For each trial, the 10 genotypes with the highest peak frequency are shown. Each sits inside the largest shown genotype whose mutations are a subset of its own, and unshown cells merge into their nearest shown ancestor, so fractions always sum to 1. This is genotype containment, not a reconstructed lineage. Hue encodes functional class and lightness separates genotypes within a class. The main panel lists every genotype by name.
- `02_gene_mutation_frequencies_{zh,en}`: mutant-cell fractions per gene across all 30 trials, aligned at therapy start, as median and 10th–90th percentile bands. Panels cover oncogenes, TSGs (any allele and biallelic loss), resistance, and the passenger mean as a neutral control. Only relative generations observed in every trial are shown, on a symlog y-axis.

`muller_representative.csv`, `gene_frequency_summary.csv` and `genotype_evolution.json` hold the plotted values.

### Mathematical diagrams

`.venv/bin/python main.py --only model` regenerates the publication-style 300 dpi PNG and vector SVG figures in `outputs/mathematical_model/`:

- `06_single_cell_rules`: the offspring distribution of a mutation-free cell, the same rule as an exact five-state matrix (mutations and therapy disabled; the overflow state is imposed bookkeeping absorption, not biological saturation), and a table of d, q, offspring probabilities, mean offspring m with and without therapy, and the probability that one cell's lineage eventually dies out (smallest root of s = f(s)). A baseline resistant lineage survives therapy with probability 0.52; sensitive lineages under therapy always die out.
- `07_matrices_to_flow`: the reduced sensitive/resistant mutation matrix Q and expected-offspring matrices M, the conditional-mean formulas E[S′ | s] and E[R′ | s], and beside them the linear S–R vector fields with therapy off and on (13×13 grid up to 30,000 cells, shared length scale). The 30 main-experiment trajectories are overlaid: the pre-therapy segment on the untreated panel and the segment from therapy start on the treated panel, with the representative trial highlighted. The actual paths follow the arrows even though the arrows use baseline cells only.
- `A1_markov_kernel_and_controller` (appendix): the full genotype-count Markov kernel written with probability-generating functions, and the execution state machine.

P and Q are row-stochastic; M contains offspring means and is not a probability transition matrix. Rows are current/parent states and columns are next/offspring states. CSV matrices and a JSON assumptions record are exported alongside the figures; analytic checks run during generation.

### Teaching presentation

`outputs/presentations/tumor_markov_code_results_zh_tw_v10.pptx` is a 54-slide Traditional Chinese teaching deck. It opens with treatment response, explains the implementation and its mathematics, then examines population composition, diversity, therapy analyses, control experiments and future tests. Mathematical derivations appear in the relevant content pages, with no appendices.

The code lessons use excerpts from `simulation.py` (and `experiments.py` for repeated trials), with variable meanings, numerical examples and detailed speaker notes. They cover Config, sparse genotype counts, bitmasks, TSG alleles, IL11, update timing, Binomial reproduction, mutation conservation, P/Q/M, relapse identification, ensemble statistics, sampling and tests. Excerpts are selected by line number; when `simulation.py` or `experiments.py` changes, the ranges in `presentations/teach_model.mjs` must be re-checked (v10 re-mapped every range to the current source). Equations use LaTeX math syntax rendered with Computer Modern fonts at high resolution; they are image objects rather than native PowerPoint equation objects. Text, matrices, diagrams, tables and charts remain editable.

v10 follows the same merge rules as the figures: the trial-availability chart became a footnote, the two gene-frequency heatmaps became one aligned gene-frequency chart with a passenger control, the full-range IL11 curve was dropped, and the TSG fraction and endpoint slides include the no-TSG-advantage control. Four slides summarize the therapy analyses (Goldie–Coldman cure curve, relapse-time prediction, the k_R × μ cure table, and adaptive versus continuous therapy). Shared treatment curves show generations 0–30, where all 30 trials are observed and alive. Log transformations are labeled explicitly.

The IL11 module was later removed from the code. The existing decks (v3–v10) still contain the IL11 slides, and the deck build scripts (`presentations/prepare_results.py`, `prepare_extended_results.py`, `teach_model.mjs`, `tell_results.mjs`) still reference the deleted IL11 outputs and old `simulation.py` line numbers, so they must be updated before building a new deck.

The data preparation reproduces the treatment populations and the no-TSG-advantage control, checks recorded CSV counts, and reads the verified outputs of the `genotype` and `therapy` steps. Rebuild the data and deck with:

```bash
.venv/bin/python presentations/prepare_results.py
.venv/bin/python presentations/prepare_extended_results.py
PRESENTATION_BUILD_DIR=.presentation-build/code-v11 PRESENTATION_FILENAME=tumor_markov_code_results_zh_tw_v11.pptx "$HOME/.cache/codex-runtimes/codex-primary-runtime/dependencies/node/bin/node" presentations/build_tumor_markov.mjs
```

The builder uses the installed bundled Artifact Tool and validators. Equation rendering uses the project-local Python without an external TeX installation. Private candidates, equation assets, renders and validation records stay in the ignored `.presentation-build/` directory. No conda environment is changed. The previous v3 and v9 decks are retained separately.

### Discussion guide

`outputs/discussion_guide/` holds a 20-minute cue sheet for the one-on-one discussion, as a web page (`discussion_guide.html`, with downscaled figures in `figs/`) and a Word file (`腫瘤模擬_20分鐘討論提示卡.docx`). A read-first code walkthrough opens it: a layer strip from `main.py` down to the plotting modules, UML-style class cards (`Config` and `Result` as dataclasses, each module as a card with its responsibility, attributes and public methods), the step-by-step flow of one `main.py` run, and talking points for architecture questions. The rest is not a script: each timed segment lists the figure to open, where to point, the code lines to show, key English sentences with Traditional Chinese glosses, concept anchors, the numbers to cite and a transition line. Likely questions and a number cheat sheet follow. The web page has a talk timer that highlights the current segment.

Both files come from one script, so they stay consistent. Code excerpts are read from the source files by line range; recheck the ranges in `SEGMENTS` after editing `simulation.py` or `analysis.py`.

```bash
.venv/bin/python presentations/build_discussion_guide.py
```

### Customize and reuse

Code style: all Python files are formatted with black (settings in `pyproject.toml`; long expressions are split one argument or condition per line). Install it into the project environment with `.venv/bin/python -m pip install black`, then run `.venv/bin/python -m black --unstable *.py presentations/*.py`. Each plotting style function first restores Matplotlib's initial settings, so running every step in one process gives the same figures as running them separately. A few chained conditionals are kept in a one-condition-per-line layout with `# fmt: off`. Formatting moves lines, so recheck the code-excerpt line ranges in `presentations/teach_model.mjs` before rebuilding the deck.

Edit parameters in `Config`, or call the model from another Python program:

```python
from pathlib import Path
from simulation import Config, simulate
from plot_results import plot_dashboard

config = Config(suppressor_hits=1, mutation_rate=0.0005)
result = simulate(config, seed=7, therapy=True)
plot_dashboard(result, Path("my_figures"))
```

Run the allele-state, mutation-conservation, reproducibility and ensemble-statistics tests with `.venv/bin/python -m unittest -v test_simulation.py`.

Experiment functions in `experiments.py` return plain data (lists, dicts, arrays), and the `plot_*` functions take that data, so either side can be reused alone. The Agg backend saves figures directly without opening a GUI.

Discussion topics include how drivers change expected offspring per generation; why requiring loss of both TSG alleles delays selective advantage; how pre-existing resistant cells become dominant during treatment; and how death rate, time to target size, and mutation accumulation relate at a fixed sample size. Interpret trends using the actual outputs rather than treating one random trajectory as a universal result.

## 繁體中文

這是一個用於課堂討論的隨機離散分枝模型。參數是教學假設，未用臨床資料校準。程式及註解以中文說明；核心模擬圖表使用英文標籤，獨立教學圖及相圖使用繁體中文標籤。

### 執行（不修改 conda 環境）

直接指定專案虛擬環境，不必 activate：

```bash
.venv/bin/python main.py                                   # 全部步驟（12 核心約 3 分鐘）
.venv/bin/python main.py --only results                    # 只跑主實驗與對照實驗
.venv/bin/python main.py --only results --suppressor-hits 1
.venv/bin/python main.py --only therapy --analyses 08 09   # 只跑部分治療分析
```

步驟依序為：`results`（圖 01–03）、`model`（06、07、A1）、`phase`（相圖）、`genotype`（Muller 圖與基因頻率）、`therapy`（08–11）。後面的步驟讀取 `results` 寫出的 CSV，所以剛取得專案時要先跑 `results`。輸出一律寫到 `outputs/`；`--seed` 與 `--replicates` 改變主實驗與對照實驗。

### 程式架構

| 檔案 | 內容 |
| --- | --- |
| `simulation.py` | 只有模型：`Config`、基因型 → 表型、`growth_probabilities`、`mutate_offspring`、`simulate`、`therapy_outcome`。解釋程式時從這裡開始。 |
| `experiments.py` | 在不同條件下重複執行模型（主實驗、TSG／死亡率對照、治療分析），彙整試驗結果並讀寫資料檔。 |
| `analysis.py` | 由模型規則推出的公式：期望子代與一步漂移、分枝過程滅絕機率、Goldie–Coldman 理論、相平面向量場、Muller 版面與基因頻率。 |
| `plot_style.py`、`plot_results.py`、`plot_model.py`、`plot_phase.py`、`plot_genotype.py`、`plot_therapy.py` | 只畫圖，依主題分檔；各自使用上面算好的資料。 |
| `main.py` | 執行入口：`main()` 依序呼叫每個步驟的短函數。 |
| `test_simulation.py` | 模型與多次試驗統計的單元測試。 |

如果在新電腦上使用，先以你的 Python 建立獨立環境並安裝依賴：

```bash
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
```

VS Code 可將 Python interpreter 選為此資料夾的 `.venv/bin/python`。套件位於 `.venv`，Matplotlib 快取預設位於本專案的 `.mplconfig`。

### 模型假設

- 初始有 100 個未突變細胞、20 個基因座：3 個 oncogenes（Python 索引 0–2）、1 個抗藥基因（3）、2 個抑癌基因 TSG（4–5）及 14 個 passengers（6–19）。每個 TSG 各有兩個可獨立突變的等位基因，其餘基因採二元突變狀態。基因型使用 22 個位元表示 20 個基因，並非 22 個基因。
- 每一代父細胞由 0、2 或 3 個子細胞取代；不另保留父細胞。死亡機率為 `d`，存活後三子細胞機率為 `q`，因此每個父細胞的期望子代數為 `(1-d)*(2+q)`。
- Oncogene 活化與 TSG 功能喪失各自降低死亡率並提高三子細胞機率，兩類效果可相加，每類最多計算一次。`driver_hits` 保留為不同 oncogenes 的突變數門檻。`suppressor_hits=1` 表示同一 TSG 任一等位基因失活即產生優勢；`suppressor_hits=2`（預設）則需同一 TSG 的兩個等位基因都失活。不同 TSG 各一擊不符合 two-hit 條件。One-hit 是簡化的單倍劑量不足情境，不代表所有抑癌基因皆如此。
- 子細胞繼承突變，每個尚未突變位元再以 `mutation_rate` 獨立突變。無回復突變，可在同一代累積多個突變。預設每位元每子細胞突變率 0.001（TSG 為每等位基因）是為短時間演示而設定。TSG 的兩個等位基因可在同一代都突變。失活以不可逆突變表示，未另外區分雜合性缺失（LOH）、缺失或表觀遺傳沉默。
- 按基因型計數並使用二項分布抽樣，與此假設下逐細胞抽樣等價，但更有效率。同基因型可來自不同祖先，因此圖中的 genotype 並不等同於譜系 clone。
- 腫瘤首次達到 20,000 個細胞時，下一次世代轉移開始持續治療。自然存活後，敏感細胞有 90% 額外死亡機率，抗藥細胞有 4%；總死亡率為 `1-(1-d)*(1-kill)`。
- 抗藥突變可在治療前或治療期間出現；治療不會直接提高突變率。新突變在子代產生時發生，其生存優勢於下一代生效。
- 復發定義：治療後先下降至治療開始負荷以下，再回升至該負荷。未觀察到復發不代表永久治癒。
- `max_cells` 是計算停止門檻，保留跨過門檻的完整一代，不是承載量或截斷取樣。模型沒有空間、免疫、資源競爭或靜止細胞。

### 展示圖表與資料

`outputs/` 內含：

| 檔案 | 用途 |
| --- | --- |
| `01_tumor_dashboard.png / .svg` | 治療下多次試驗的負荷、有效基因型數與抗藥比例；標示治療開始時已存在的抗藥細胞 |
| `02_death_rate_heterogeneity.png / .svg` | 三種死亡率下 5,000 細胞樣本的 Shannon 多樣性，以及長到 5,000 細胞所需代數 |
| `03_one_hit_vs_two_hit.png / .svg` | 無 TSG 優勢對照組、one-hit、two-hit：線性生長、TSG 達標比例與第 18 代各次試驗負荷 |
| `one_hit_vs_two_hit.csv / one_hit_vs_two_hit_summary.csv` | 逐代比較資料與終點中位數 |
| `trajectory.csv` | 多次試驗逐代統計：n、平均、樣本 SD、中位數、p10、p90 |
| `therapy_trajectories.csv` | 彙整圖使用的各次治療試驗逐代資料 |
| `comparisons.csv` | 重複試驗數據，包含未達到指定大小的試驗 |
| `therapy_replicates.csv` | 多次治療試驗的最低負荷、復發代數、停止原因 |
| `run_metadata.json` | 參數、隨機種子與套件版本 |

PNG 為 240 dpi；SVG 可無損縮放。生長曲線使用 symlog，0 可顯示，1 以上近似對數尺度。`suppressor_fraction` 記錄符合所選 TSG 門檻的細胞比例；`suppressor_biallelic_fraction` 記錄至少一個 TSG 雙等位基因失活的細胞比例。`driver_fraction` 仍只表示達到 oncogene 門檻的比例。

生長比較固定 oncogene 效果、突變率與基因座，只改變 TSG 優勢條件：取消優勢（對照組）、一擊或兩擊，三者都畫在圖 03。預設展示 30 次試驗的中位數及 10–90 百分位區間，這是試驗間變異，並不是信賴區間。死亡率比較使用各試驗首次跨過 5,000 細胞的世代，無放回取樣 5,000 細胞計算 Shannon H；圖上的橫線是中位數。此結果是在達門檻條件下的等大小**樣本**比較，並非精確命中該大小的完整腫瘤。圖 02 的 B 面板畫出每次試驗首次達到 5,000 細胞的世代（死亡率 0.20／0.35／0.45 時中位數 8／14／26.5 代）：代數越多，分裂越多，同樣大小時累積的突變也越多。未達門檻者不納入，但會保留於 CSV 並顯示成功試驗數。

### One-hit 與 two-hit 直接比較

![One-hit 與 two-hit 比較](outputs/03_one_hit_vs_two_hit.png)

此圖呈現生長試驗每種情境的同一批 30 次試驗，並加入無 TSG 優勢的對照組（灰色點線；TSG 失活以 two-hit 規則計數，但沒有生長優勢），不施加治療。第 18 代負荷中位數為 76,446（對照）、94,618（one-hit）、73,879（two-hit），達到 TSG 門檻的細胞分別占 0.05%／27.9%／0.07%。A 面板以線性座標顯示細胞數，B 面板顯示各情境達到 TSG 失活門檻的細胞百分比（symlog，0.01% 以下為線性），C 面板顯示第 18 代每次試驗的負荷與中位數橫線。兩情境的門檻定義不同，B 比較的是功能表型，而非完全相同的突變狀態。陰影表示試驗間變異，不是信賴區間；保留重疊與隨機反轉，不挑選有利種子。

### 重複試驗與不確定性

`--replicates` 統一控制**所有隨機實驗組**，包含主治療圖，預設 30 次。例如：`.venv/bin/python main.py --only results --replicates 100`。主圖改用多次試驗彙整，不再只展示一條軌跡。細胞數與多樣性將滅絕試驗延續為零；達計算上限後的後續值視為缺失，不使用最後數量補值；註腳寫出所有試驗都有觀察的最後世代（預設第 31 代）。抗藥比例只針對仍存活且有觀察值的腫瘤，空族群不定義比例。逐基因突變頻率見 `outputs/genotype_evolution/02_gene_mutation_frequencies_{zh,en}`。

試驗依絕對世代對齊，治療各自在達到大小門檻後開始，因此主圖以淡色區間標出治療起點的範圍，而不是單一條線；C 面板標出各試驗在自己治療起點的抗藥比例中位數。`therapy_replicates.csv` 保留每次起點與結果。`trajectory.csv` 現為多次試驗摘要，`therapy_trajectories.csv` 保留各次軌跡。仍可自行呼叫 `plot_dashboard(result, ...)` 查看明確指定的單次案例。

生長比較顯示中位數及 10–90% 試驗區間；終點負荷與多樣性散點保留個別試驗及中位數。增加次數能穩定摘要估計，但不會消除生物隨機性；百分位帶代表試驗間變異，**不是信賴區間或標準誤**。

### 相平面圖

執行 `.venv/bin/python main.py --only phase`，即可在 `outputs/phase_plots/` 生成一張 300 dpi PNG／向量 SVG：`01_phase_portrait`，不修改簡報。圖中以 `log10(count+1)` 座標呈現敏感／抗藥細胞數。灰線保留所有實際觀察軌跡；突出一條實際試驗，先選最接近復發中位代數，再選谷底最接近中位數的案例，並依治療前、治療至復發、復發後上色，標示治療啟動、最低負荷與復發。黑色箭頭的端點是由突出案例完整基因型快照算出的下一代條件期望，包含敏感細胞突變成抗藥細胞的精確貢獻；箭頭頭端到下一個實際點的距離就是該代的隨機偏差。右欄列出治療啟動與谷底時的期望與實際下一代。重現案例會逐代核對記錄的 CSV 計數。箭頭端點是原始計數期望的對數轉換，並非對數計數的期望；公式見圖 07。

這是完整基因型／治療狀態的投影，並非封閉的二維自治向量場：相同 S／R 數目仍可能有不同 driver 或 TSG 組成。計算終點不是平衡點。`conditional_phase_steps.csv` 保留期望及實際下一代計數；`phase_summary.json` 記錄案例選法、事件及解讀限制。可用 `.venv/bin/python main.py --only phase --replicate 0` 突出另一條實際試驗（會覆寫相圖檔案）。中文標籤需要 Arial Unicode MS 或 Noto Sans CJK TC 等支援字型；SVG 保留字形路徑。負荷 × 抗藥比例的視角改看下方的 `05_phase_field_burden_logit`，它用 logit 軸取代線性百分比。

同一步驟也計算網格向量場資料（圖 07 使用相同公式）：`phase_vector_grid.csv` 與 `phase_vector_assumptions.json`。每個網格點都指定完整父代族群：S 的基因型為 0，R 只有抗藥位元為 1；全部 1,604 點（線性與對數網格）的一步投影期望均與完整模型核對。子代仍可發生其他突變，後續狀態可能離開此截面。治療在每面板固定，不依大小門檻切換。原本獨立的線性與對數向量場圖已併入圖 07 或刪除。

### 理論與經驗相場

相圖步驟也畫出兩組座標，每組皆以「理論 vs 經驗」漂移場呈現，並分開治療關閉與持續治療：

- `05_phase_field_burden_logit_{zh,en}`：x = log10(N+1)，y = 以經驗 logit log10((R+½)/(S+½)) 呈現的抗藥比例，刻度標為百分比。ρ 為 0 或 1 時仍有限值，虛線與點線分別標示 S = 0 與 R = 0。
- `06_phase_field_fitness_logit_{zh,en}`：x = 無治療平均適應度 w̄_off，y 軸同上。這張圖把演化（oncogene 累積）和選擇（抗藥）分開。

理論箭頭使用固定的隱藏基因型截面：敏感細胞為基因型 0，抗藥細胞只帶抗藥位元。適應度平面中，oncogene 比例由 w̄ 換算，與抗藥獨立，N 取治療門檻。箭頭端點是一步原始計數期望經座標轉換後的位置，並與 `projected_offspring_mean` 核對。經驗箭頭是起點落在同一格的所有一步轉移之平均位移，僅在該格有 ≥5 次轉移且來自 ≥2 次試驗時繪製；顏色表示轉移數。同一治療欄中，理論與經驗箭頭共用長度比例，x、y 各有比例鍵。治療中理論 logit ρ 每代約上升 0.98 個 10 倍：ρ → 1 是組成上的穩定方向，但 log N 沒有不動點。上下排的差異反映被投影掉的自由度與隨機性。`phase_field_comparison.csv/.json` 記錄每支箭頭、分格邊界與假設。

### 治療延伸分析

執行 `.venv/bin/python main.py --only therapy`（12 核心約 2 分鐘；`--analyses 08 09` 只跑部分），會在 `outputs/therapy_analysis/` 輸出四張圖與 CSV／JSON 資料。分析 09 重用主實驗記錄的 30 次試驗，其餘使用獨立的 seed 區段。所有參數都是教學假設。

- `08_goldie_coldman`：開始持續治療時的負荷（100–5,000 細胞，每個大小 200 次，Wilson 95% 區間）與治癒機率，並與沒有擬合參數的分枝過程理論比較。抗藥突變 clone 以 Poisson 過程出現；在治療前 s 代出生的 clone 最終存活機率為 1 − f^(s)(q_on)，f 是未治療子代生成函數，q_on = 0.479 是一個抗藥家系在治療下的滅絕機率。治癒 ≈ exp(−Λ)。理論與模擬吻合，50% 治癒點約在 416 個細胞（100 個細胞時 99.5%，1,000 時 10%，5,000 時 0%）。B 面板是治療開始時已存在的抗藥細胞，平均數大於中位數（Luria–Delbrück 的 jackpot）。
- `09_relapse_time_prediction`：用治療開始時的基因型組成預測主實驗每次試驗的復發延遲。基準預測 ln(N/R₀)/ln 1.371 平均誤差 2.0 代；讓每種抗藥基因型以自己的治療中平均子代數成長（oncogene 細胞 ×1.62），30 次中 24 次預測到正確世代（平均誤差 0.2）。抗藥細胞越多（Spearman −0.70）、其中 oncogene 比例越高（−0.51），復發越早。C 面板是治療中期望成長率對抗藥比例 ρ，在 ρ* = (1 − m_S)/(m_R − m_S) ≈ 70% 穿過零，也就是谷底。
- `10_parameter_sensitivity`：抗藥細胞殺傷率 k_R（0–0.6）× 突變率 μ（10⁻⁴–3×10⁻³）的治癒機率與復發延遲中位數，每格 30 次。單純抗藥細胞在 k_R = 0.30 停止成長，但治癒需要約 0.41，也就是帶 oncogene 的抗藥細胞也停止成長的殺傷率；μ 越高，這種細胞越可能事先存在。
- `11_adaptive_therapy`：間歇治療（負荷降到治療起點的 50% 就暫停，回到 20,000 再恢復）與持續治療，每種設定 30 對試驗。惡化 = 負荷曾低於 20,000 後，在用藥中仍 ≥ 24,000；要求「用藥中」是為了不把暫停期的一步成長算成治療失敗。沒有競爭時間歇治療較差（中位比 ×0.85；加抗藥代價 ×0.91）。有承載量競爭（K = 30,000）時延後惡化 ×1.14（14 → 16 代，30/30 對較長；加代價 17.5 → 20 代），且只用約一半的藥量。

適應性治療分析使用三個預設關閉的 `Config` 選項，不影響其他輸出：`therapy_mode="adaptive"` 與 `adaptive_pause_fraction`、`resistance_cost`（抗藥細胞額外自然死亡率），以及 `carrying_capacity`（每個細胞額外死亡 0.21·N/K，使無突變細胞在 N = K 時平均子代為 1）。前一版的 PCA 相座標分析已刪除，因為它的主成分只反映時間與大小。

### 基因型分布的演變

執行 `.venv/bin/python main.py --only genotype`，輸出至 `outputs/genotype_evolution/`：

- `01_muller_top_genotypes_{zh,en}`：代表試驗（與相平面圖相同的選取規則）的 Muller 圖，上方附負荷曲線；另有 6 個由早到晚復發的其他試驗。每個試驗呈現峰值頻率最高的 10 個基因型。每個基因型位於「已呈現且突變為其子集的最大基因型」內，未呈現的細胞併入最接近的已呈現祖先，因此比例總和恆為 1。這是基因型包含關係，不是重建的譜系。色相表示功能類別，同類別內以明度區分基因型；主圖逐一列出基因型名稱。
- `02_gene_mutation_frequencies_{zh,en}`：30 次試驗每個基因的突變細胞比例，以治療開始對齊，呈現中位數與 10–90 百分位帶。面板涵蓋 oncogene、TSG（任一等位基因與雙等位基因失活）、抗藥基因，以及作為中性對照的 passenger 平均。只呈現所有試驗都有觀察的相對世代，縱軸為 symlog。

`muller_representative.csv`、`gene_frequency_summary.csv` 與 `genotype_evolution.json` 保存繪圖數值。

### 數學模型圖

執行 `.venv/bin/python main.py --only model`，可在 `outputs/mathematical_model/` 重新生成論文風格的 300 dpi PNG 與向量 SVG：

- `06_single_cell_rules`：未突變細胞的子代分布、同一規則寫成的精確五狀態矩陣（停用突變與治療；溢出狀態是人為設定的記帳吸收狀態，不是生物飽和），以及 d、q、子代機率、治療前後平均子代數 m 和單一細胞家系最終滅絕機率（s = f(s) 的最小根）的表格。基準抗藥細胞的家系在治療下有 0.52 的機率存活；敏感細胞家系在治療下必定滅絕。
- `07_matrices_to_flow`：簡化的敏感／抗藥突變矩陣 Q、期望子代矩陣 M、條件期望公式 E[S′ | s] 與 E[R′ | s]，旁邊是治療關閉與持續治療的線性 S–R 向量場（13×13 網格至 30,000 細胞，共用長度比例）。疊上主實驗 30 次試驗的實際路徑：未治療面板畫治療前那段，治療面板畫治療開始之後那段，並突出代表試驗。雖然箭頭只用基準細胞計算，實際路徑仍順著箭頭走。
- `A1_markov_kernel_and_controller`（附錄）：以機率生成函數寫出的完整基因型計數 Markov 轉移核，以及程式執行狀態機。

P 與 Q 每列加總為 1；M 記錄期望子代數，不是機率轉移矩陣。列代表當前／父代，欄代表下一步／子代。圖旁另存 CSV 矩陣及 JSON 假設紀錄，生成時會執行解析數值檢查。

### 教學簡報

`outputs/presentations/tumor_markov_code_results_zh_tw_v10.pptx` 是 54 張繁體中文教學投影片。以治療反應開場，講解程式與數學，再依序檢查族群組成、多樣性、治療延伸分析、控制實驗與後續驗證方向。數學推導整合於相應主內容頁，沒有附錄。

程式教學採用 `simulation.py` 的節錄（重複試驗部分取自 `experiments.py`），搭配變數意義、數值例子及詳細講者備註，涵蓋 Config、稀疏基因型計數、bitmask、TSG 等位基因、IL11、更新時序、二項繁殖抽樣、突變守恆、P／Q／M、復發判讀、多次試驗統計、取樣及測試。節錄依行號選取；`simulation.py` 或 `experiments.py` 改動後必須重新檢查 `presentations/teach_model.mjs` 的行號範圍（v10 已全部對應到目前的原始碼）。公式以 LaTeX 數學語法及 Computer Modern 字型高解析度排版，是圖片物件，並非 PowerPoint 原生公式物件。文字、矩陣、狀態圖、表格與圖表均可編輯。

v10 依照與圖表相同的合併原則：可用試驗數圖改為註腳，兩張基因頻率熱圖改為一張以治療開始對齊、含 passenger 對照的基因頻率圖，刪除 IL11 全範圍曲線，TSG 達標比例與終點投影片加入無 TSG 優勢對照組。另有四張治療延伸分析投影片（Goldie–Coldman 治癒曲線、復發時間預測、k_R × μ 治癒表、適應性與持續治療比較）。主要治療曲線顯示 0–30 代的共同完整窗口，每代皆有 30 次存活且觀察中的試驗。所有對數轉換均明確標示。

IL11 模組之後已從程式中刪除。既有簡報（v3–v10）仍保留 IL11 投影片；簡報建置腳本（`presentations/prepare_results.py`、`prepare_extended_results.py`、`teach_model.mjs`、`tell_results.mjs`）仍引用已刪除的 IL11 輸出與舊的 `simulation.py` 行號，重建新簡報前必須先更新。

資料準備會重現治療族群與無 TSG 優勢對照組，核對 CSV 計數，並讀取 `genotype` 與 `therapy` 步驟已核對的輸出。重建資料與簡報：

```bash
.venv/bin/python presentations/prepare_results.py
.venv/bin/python presentations/prepare_extended_results.py
PRESENTATION_BUILD_DIR=.presentation-build/code-v11 PRESENTATION_FILENAME=tumor_markov_code_results_zh_tw_v11.pptx "$HOME/.cache/codex-runtimes/codex-primary-runtime/dependencies/node/bin/node" presentations/build_tumor_markov.mjs
```

繪製程式使用已安裝的 bundled Artifact Tool 與檢查工具。公式使用專案 Python 排版，不需外部 TeX 安裝。候選檔、公式素材、預覽與檢查紀錄放在 Git 忽略的 `.presentation-build/`，不修改 conda。前版 v3 與 v9 簡報另行保留。

### 討論提示卡

`outputs/discussion_guide/` 是一對一討論用的 20 分鐘提示卡，有網頁（`discussion_guide.html`，縮小的圖在 `figs/`）與 Word（`腫瘤模擬_20分鐘討論提示卡.docx`）兩種版本。開頭是「講前先讀」的程式導覽：從 `main.py` 到繪圖模組的分層、UML 風格的類別卡（`Config` 與 `Result` 為 dataclass，每個模組一張卡，列出職責、屬性與對外方法）、執行一次 `main.py` 的逐步流程，以及被問到程式架構時可用的關鍵句。其餘內容不是逐字稿：每個時段列出要打開的圖、指圖順序、要對照的程式行、可直接說出口的英文關鍵句（附繁體中文意思）、概念錨點、要引用的數字與銜接句；最後附老師可能問的問題與數字小抄。網頁有計時器，會標出目前該講的段落。

兩個檔案由同一個腳本產生，內容一致。程式節錄依行號從原始檔讀取；修改 `simulation.py` 或 `analysis.py` 後，須重新檢查 `SEGMENTS` 內的行號。

```bash
.venv/bin/python presentations/build_discussion_guide.py
```

### 修改與重用

程式風格：所有 Python 檔以 black 排版（設定在 `pyproject.toml`；長運算式拆成每個參數或條件一行）。以 `.venv/bin/python -m pip install black` 安裝到專案環境後，執行 `.venv/bin/python -m black --unstable *.py presentations/*.py`。每個繪圖樣式函數都先還原 Matplotlib 初始設定，所以在同一個行程跑完所有步驟，圖與分開執行時完全相同。少數連續條件式以 `# fmt: off` 保留「一個條件一行」的寫法。排版會移動行號，重建簡報前須重新檢查 `presentations/teach_model.mjs` 的程式節錄行號。

修改 `Config` 中的參數，或從其他 Python 程式呼叫：

```python
from pathlib import Path
from simulation import Config, simulate
from plot_results import plot_dashboard

config = Config(suppressor_hits=1, mutation_rate=0.0005)
result = simulate(config, seed=7, therapy=True)
plot_dashboard(result, Path("my_figures"))
```

執行 `.venv/bin/python -m unittest -v test_simulation.py` 可驗證等位基因狀態、突變守恆、可重現性與多次試驗統計。

`experiments.py` 的實驗函數回傳一般資料（list、dict、陣列），`plot_*` 函數接收這些資料，兩邊都可以單獨重用。使用 Agg backend，執行後直接存檔，不彈出 GUI。

討論時可著重：driver 如何改變每代期望子代數；TSG 雙等位基因失活條件為何延遲選擇優勢；治療前已存在的抗藥細胞如何在治療後占優；相同樣本大小下，死亡率、達到大小所需代數和突變累積的關係。請以實際輸出數據解釋趨勢，避免將單一隨機案例當成普遍定律。
