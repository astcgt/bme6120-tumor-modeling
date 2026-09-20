# BME 6120 — Quantitative tumor modeling

這是一個用於課堂討論的隨機離散分枝模型。參數是教學假設，未用臨床資料校準。程式及註解以中文說明，圖表使用英文標籤，方便簡報且不依賴中文字型。

## 執行（不修改 conda 環境）

直接指定專案虛擬環境，不必 activate：

```bash
.venv/bin/python modeling_cancer.py
.venv/bin/python modeling_cancer.py --seed 123 --replicates 30 --output outputs_seed123
```

如果在新電腦上使用，先以你的 Python 建立獨立環境並安裝依賴：

```bash
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
```

VS Code 可將 Python interpreter 選為此資料夾的 `.venv/bin/python`。套件位於 `.venv`，Matplotlib 快取預設位於本專案的 `.mplconfig`。

## 模型假設

- 初始有 100 個未突變細胞、12 個二元基因。Python 索引 0–2 是 driver，3 是抗藥基因，其餘是 passenger。
- 每一代父細胞由 0、2 或 3 個子細胞取代；不另保留父細胞。死亡機率為 `d`，存活後三子細胞機率為 `q`，因此每個父細胞的期望子代數為 `(1-d)*(2+q)`。
- 達到 driver 門檻會降低死亡率並提高產生三個子細胞的機率。多出的 driver 不再疊加優勢。「兩擊」是兩個不同 driver genes 突變，不是同一抑癌基因的兩個等位基因失活。
- 子細胞繼承突變，每個尚未突變基因再以 `mutation_rate` 獨立突變。無回復突變，可在同一代累积多個突變。預設每基因每子細胞突變率 0.001 是為短時間演示而設定。
- 按基因型計數並使用二項分布抽樣，與此假設下逐細胞抽樣等價，但更有效率。同基因型可來自不同祖先，因此圖中的 genotype 並不等同於譜系 clone。
- 腫瘤首次達到 20,000 個細胞時，下一次世代轉移開始持續治療。自然存活後，敏感細胞有 90% 額外死亡機率，抗藥細胞有 4%；總死亡率為 `1-(1-d)*(1-kill)`。
- 抗藥突變可在治療前或治療期間出現；治療不會直接提高突變率。新突變在子代產生時發生，其生存優勢於下一代生效。
- 復發定義：治療後先下降至治療開始負荷以下，再回升至該負荷。未觀察到復發不代表永久治癒。
- `max_cells` 是計算停止門檻，保留跨過門檻的完整一代，不是承載量或截斷取樣。模型沒有空間、免疫、資源競爭或靜止細胞。

## 展示圖表與資料

`outputs/` 內含：

| 檔案 | 用途 |
| --- | --- |
| `01_tumor_dashboard.png / .svg` | 細胞數、敏感／抗藥族群、基因型組成、多樣性、基因突變頻率 |
| `02_scenario_comparisons.png / .svg` | 無選擇優勢／一擊／兩擊的生長，以及不同死亡率下的異質性 |
| `trajectory.csv` | 單次治療模擬逐代統計 |
| `comparisons.csv` | 重複試驗數據，包含未達到指定大小的試驗 |
| `therapy_replicates.csv` | 多次治療試驗的最低負荷、復發代數、停止原因 |
| `run_metadata.json` | 參數、隨機種子與套件版本 |

PNG 為 240 dpi；SVG 可無損縮放。生長曲線使用 symlog，0 可顯示，1 以上近似對數尺度。組成图中的 Gxxx 是基因型 bitmask 的十六進位表示，R 表示含抗藥突變；主要基因型依跨世代累積比例選出。

生長比較展示 12 次試驗的中位數及 10–90 百分位區間，這是試驗間變異，並不是信賴區間。死亡率比較使用各試驗首次跨過 5,000 細胞的世代，無放回取樣 5,000 細胞計算 Shannon H；圖上的橫線是中位數。此結果是在達門檻條件下的等大小**樣本**比較，並非精確命中該大小的完整腫瘤。未達門檻者不納入 H，但會保留於 CSV 並顯示成功試驗數。

## 修改與重用

修改 `Config` 中的參數，或從其他 Python 程式呼叫：

```python
from pathlib import Path
from modeling_cancer import Config, simulate, plot_dashboard

config = Config(driver_hits=2, mutation_rate=0.0005)
result = simulate(config, seed=7, therapy=True)
plot_dashboard(result, Path("my_figures"))
```

繪圖函數 `plot_dashboard`、`run_comparisons`、`set_plot_style`、`save_figure` 可獨立重用。使用 Agg backend，執行後直接存檔，不彈出 GUI。

討論時可著重：driver 如何改變每代期望子代數；兩擊條件為何延遲選擇優勢；治療前已存在的抗藥細胞如何在治療後占優；相同樣本大小下，死亡率、達到大小所需代數和突變累積的關係。請以實際輸出數據解釋趨勢，避免將單一隨機案例當成普遍定律。
