"""理論與衍生量：由模型規則直接推出的公式，以及從模擬結果計算的統計量。

這裡不執行重複試驗（在 experiments.py），也不畫圖（在 plot_*.py）。
內容：統計工具、適應度與一步期望、分枝過程理論、相圖向量場、基因型組成。
"""

from __future__ import annotations

import numpy as np

from simulation import (
    Config,
    Result,
    driver_count,
    growth_probabilities,
    suppressor_alleles,
    suppressor_loss,
)


# ---------------------------------------------------------------- 統計工具 ----------
def wilson(successes: int, n: int, z: float = 1.96) -> tuple[float, float]:
    """二項比例的 Wilson 95% 區間。"""
    if n == 0:
        return np.nan, np.nan
    p = successes / n
    centre = (p + z * z / (2 * n)) / (1 + z * z / n)
    half = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / (1 + z * z / n)
    return centre - half, centre + half


def spearman(x, y) -> float:
    rx, ry = np.argsort(np.argsort(x)), np.argsort(np.argsort(y))
    return float(np.corrcoef(rx, ry)[0, 1])


# ---------------------------------------------------------------- 適應度與一步期望 --
def fitness(genotype: int, config: Config, therapy: bool) -> float:
    """每個父細胞的期望子代數 w = (1-d)(1-k)(2+q)；k 只在治療時作用。"""
    death, triple = growth_probabilities(genotype, config)
    kill = 0.0
    if therapy:
        resistant = bool(genotype & (1 << config.resistance_gene))
        kill = (
            config.resistant_therapy_kill
            if resistant
            else config.sensitive_therapy_kill
        )
    return (1 - death) * (1 - kill) * (2 + triple)


def expected_next(
    population: dict[int, float], config: Config, therapy: bool
) -> dict[str, float]:
    """一步原始計數期望：S、R 與子代（突變後）的平均無治療適應度 w̄_off。

    w_off 只取決於 oncogene 是否活化與 TSG 是否失活；兩者及抗藥位元的突變彼此獨立，
    因此可對每個父代基因型精確計算子代的期望 w_off，不需要展開所有子代基因型。
    """
    mu = config.mutation_rate
    gain_onc = 1 - (1 - mu) ** len(config.driver_genes)
    before = {"S": 0.0, "R": 0.0}
    weighted_w = 0.0
    for genotype, count in population.items():
        offspring = count * fitness(genotype, config, therapy)
        resistant = bool(genotype & (1 << config.resistance_gene))
        before["R" if resistant else "S"] += offspring
        onc = driver_count(genotype, config) >= config.driver_hits
        p_onc = 1.0 if onc else (gain_onc if config.driver_hits == 1 else 0.0)
        # 每個 TSG 尚缺 h 個等位基因時，一代內補齊的機率為 μ^h；不同 TSG 獨立。
        keep = 1.0
        for gene in config.suppressor_genes:
            missing = max(
                0, config.suppressor_hits - suppressor_alleles(genotype, gene, config)
            )
            keep *= 1 - mu**missing
        p_tsg = 1 - keep
        w = {
            (o, t): (
                (
                    1
                    - max(
                        0,
                        config.death_rate
                        - config.driver_death_reduction * o
                        - config.suppressor_death_reduction * t,
                    )
                )
                * (
                    2
                    + config.triple_probability
                    + config.driver_triple_bonus * o
                    + config.suppressor_triple_bonus * t
                )
            )
            for o in (0, 1)
            for t in (0, 1)
        }
        mean_w = sum(
            w[o, t] * (p_onc if o else 1 - p_onc) * (p_tsg if t else 1 - p_tsg)
            for o in (0, 1)
            for t in (0, 1)
        )
        weighted_w += offspring * mean_w
    total = before["S"] + before["R"]
    return {
        "S": (1 - mu) * before["S"],
        "R": before["R"] + mu * before["S"],
        "N": total,
        "w_off": weighted_w / total if total else np.nan,
    }


def projected_offspring_mean(
    population: dict[int, int], config: Config, therapy_active: bool
) -> dict[str, float]:
    """精確的條件期望；其他基因突變不改變子代的 S/R 分組機率。

    S 子代以 μ 成為 R，既有 R 不逆轉。B_S/B_R 是突變前的期望子代數。
    這不是只依賴 (S,R) 的封閉模型：m_g 需要完整基因型。
    """
    before = {"sensitive": 0.0, "resistant": 0.0}
    for genotype, count in population.items():
        death, triple = growth_probabilities(genotype, config)
        resistant = bool(genotype & (1 << config.resistance_gene))
        kill = (
            config.resistant_therapy_kill
            if resistant
            else config.sensitive_therapy_kill
        )
        mean = (1 - death) * (1 - kill if therapy_active else 1) * (2 + triple)
        before["resistant" if resistant else "sensitive"] += count * mean
    sensitive = (1 - config.mutation_rate) * before["sensitive"]
    resistant = before["resistant"] + config.mutation_rate * before["sensitive"]
    assert np.isclose(sensitive + resistant, sum(before.values()), rtol=1e-12)
    return {
        "expected_sensitive": sensitive,
        "expected_resistant": resistant,
        "expected_total": sensitive + resistant,
        "before_mutation_sensitive": before["sensitive"],
        "before_mutation_resistant": before["resistant"],
    }


def logit_rho(resistant, sensitive):
    """經驗 logit：log10((R+½)/(S+½))；ρ=0 或 1 仍有限，且等於 log R − log S 的平滑版本。"""
    return np.log10(
        (np.asarray(resistant, float) + 0.5) / (np.asarray(sensitive, float) + 0.5)
    )


def check_expectation(config: Config) -> None:
    """expected_next 的 S/R 必須與既有 projected_offspring_mean 完全一致。"""
    population = {
        0: 700,
        1: 200,
        1 << config.resistance_gene: 90,
        (1 << config.resistance_gene) | (1 << config.suppressor_genes[0]): 10,
    }
    for therapy in (False, True):
        ours = expected_next(population, config, therapy)
        reference = projected_offspring_mean(population, config, therapy)
        assert np.isclose(ours["S"], reference["expected_sensitive"], rtol=1e-12)
        assert np.isclose(ours["R"], reference["expected_resistant"], rtol=1e-12)
    # μ=0 時子代 w̄_off 必須等於父代以期望子代數加權的 fitness 平均。
    frozen = Config(**{**config.__dict__, "mutation_rate": 0.0})
    ours = expected_next(population, frozen, False)
    weights = {g: n * fitness(g, frozen, False) for g, n in population.items()}
    direct = sum(w * fitness(g, frozen, False) for g, w in weights.items()) / sum(
        weights.values()
    )
    assert np.isclose(ours["w_off"], direct, rtol=1e-12)


# ---------------------------------------------------------------- 分枝過程理論 ------
# ---------------------------------------------------------------- 共用理論 -------------
def offspring_law(death: float, triple: float) -> np.ndarray:
    """單一父細胞的子代數分布 P(0), P(1), P(2), P(3)。"""
    return np.array([death, 0.0, (1 - death) * (1 - triple), (1 - death) * triple])


def pgf(law: np.ndarray, s: float) -> float:
    return float(law @ s ** np.arange(4))


def extinction(law: np.ndarray) -> float:
    """家系最終滅絕機率：f(s) = s 在 [0,1] 的最小根（由 0 迭代）。"""
    s = 0.0
    for _ in range(200_000):
        nxt = pgf(law, s)
        if abs(nxt - s) < 1e-15:
            return nxt
        s = nxt
    return s


def lineage_extinction(death, triple):
    """單一細胞家系最終滅絕機率：f(s)=δ+(1−δ)[(1−q)s²+qs³] 在 [0,1] 的最小不動點。"""
    s = 0.0
    for _ in range(100_000):
        s_next = death + (1 - death) * ((1 - triple) * s**2 + triple * s**3)
        if abs(s_next - s) < 1e-15:
            break
        s = s_next
    return s_next


def finite_count_matrix(death, triple, threshold=4):
    """Exact stopped one-type chain: 0,...,threshold-1, and overflow dagger.

    No mutations or therapy. Overflow is absorbing by construction, rather
    than a claim that grouping all large populations preserves Markovity.
    """
    offspring = np.array([death, 0.0, (1 - death) * (1 - triple), (1 - death) * triple])
    matrix = np.zeros((threshold + 1, threshold + 1))
    for n in range(threshold):
        law = np.array([1.0])
        for _ in range(n):
            law = np.convolve(law, offspring)
        for j, prob in enumerate(law):
            matrix[n, min(j, threshold)] += prob
    matrix[threshold, threshold] = 1
    return matrix, offspring


def markov_matrices(config: Config) -> dict:
    """單細胞規則圖與矩陣流向圖的數學物件：有限計數鏈 P、突變矩陣 Q、期望子代矩陣 M 與代表細胞表格。

    P 停用突變與治療；Q、M 只保留抗藥位元（敏感 S／抗藥 R 兩類）。
    """
    p, law = finite_count_matrix(config.death_rate, config.triple_probability)
    q = np.array([[1 - config.mutation_rate, config.mutation_rate], [0, 1.0]])
    m0 = (1 - config.death_rate) * (2 + config.triple_probability)
    off = np.diag([m0, m0]) @ q
    on = (
        np.diag([
            m0 * (1 - config.sensitive_therapy_kill),
            m0 * (1 - config.resistant_therapy_kill),
        ])
        @ q
    )
    # 解析檢查：有限分布的精確值，以及 P／Q 列加總為 1、M 列加總為平均子代數。
    np.testing.assert_allclose(p.sum(axis=1), 1, atol=1e-14)
    np.testing.assert_allclose(p[2], [0.09, 0, 0.4032, 0.0168, 0.49])
    np.testing.assert_allclose(p[3], [0.027, 0, 0.18144, 0.00756, 0.784])
    np.testing.assert_allclose(q.sum(axis=1), 1)
    np.testing.assert_allclose(off.sum(axis=1), [1.428, 1.428])
    np.testing.assert_allclose(on.sum(axis=1), [0.1428, 1.37088])
    phenotypes = []
    for name, genotype, resistant in (
        ("Baseline sensitive", 0, False),
        ("Oncogene sensitive", 1, False),
        ("Oncogene + TSG loss", 1 | (1 << 4) | (1 << 20), False),
        ("Baseline resistant", 1 << 3, True),
    ):
        d, triple = growth_probabilities(genotype, config)
        mean = (1 - d) * (2 + triple)
        if resistant:
            kill = config.resistant_therapy_kill
        else:
            kill = config.sensitive_therapy_kill
        on_death = 1 - (1 - d) * (1 - kill)
        phenotypes.append([
            name,
            f"{d:.2f}",
            f"{triple:.2f}",
            f"{d:.3f}",
            f"{(1-d)*(1-triple):.3f}",
            f"{(1-d)*triple:.3f}",
            f"{mean:.4f}",
            f"{mean*(1-kill):.4f}",
            f"{lineage_extinction(d,triple):.3f}",
            f"{lineage_extinction(on_death,triple):.3f}",
        ])
    # 子代平均小於 1 時家系必然滅絕；大於 1 時滅絕機率小於 1。
    assert phenotypes[0][9] == "1.000" and float(phenotypes[3][9]) < 1
    extinct = lineage_extinction(config.death_rate, config.triple_probability)
    np.testing.assert_allclose(extinct, law @ extinct ** np.arange(4), rtol=1e-12)
    return {
        "P": p,
        "law": law,
        "Q": q,
        "M_off": off,
        "M_on": on,
        "phenotypes": phenotypes,
    }


# ---------------------------------------------------------------- Goldie–Coldman ------
def goldie_coldman_theory(config: Config, threshold: float) -> dict:
    """突變 clone 的 Poisson 近似（Iwasa–Nowak–Michor 2006 型）。

    確定性敏感族群 S_t = N0·m^t；第 s 代出生的抗藥突變數期望 μ·m·S_{s−1}。
    該 clone 先未治療繁殖 T−s 代，再進入治療：最終存活機率 1 − f^{(T−s)}(q_on)，
    f 為未治療子代 PGF、q_on 為治療下抗藥家系滅絕機率。
    治療中敏感細胞仍會產生新抗藥突變，每個存活機率 1 − q_on。治癒 = 沒有任何抗藥家系存活 ≈ exp(−Λ)。
    只用基準細胞（忽略 oncogene／TSG），因此是近似。
    """
    d, q, mu = config.death_rate, config.triple_probability, config.mutation_rate
    off = offspring_law(d, q)
    on_r = offspring_law(1 - (1 - d) * (1 - config.resistant_therapy_kill), q)
    m = (1 - d) * (2 + q)
    m_s_on = m * (1 - config.sensitive_therapy_kill)
    q_on = extinction(on_r)
    T = int(
        np.ceil(
            np.log(max(threshold, config.initial_cells) / config.initial_cells)
            / np.log(m)
            - 1e-12
        )
    )
    S = config.initial_cells * m ** np.arange(T + 1)
    lam_pre, mean_r0, survive = 0.0, 0.0, q_on
    # 從治療開始往回推：第 s 代出生的 clone 還有 T−s 代未治療繁殖。
    for s in range(T, 0, -1):
        births = mu * m * S[s - 1]
        lam_pre += births * (1 - survive)
        mean_r0 += births * m ** (T - s)
        survive = pgf(off, survive)
    lam_on = mu * S[T] * m_s_on / (1 - m_s_on) * (1 - q_on)
    return {
        "threshold": threshold,
        "T": T,
        "N_T": float(S[T]),
        "lambda_pre": lam_pre,
        "lambda_on": lam_on,
        "p_cure": float(np.exp(-(lam_pre + lam_on))),
        "mean_R0": mean_r0,
        "q_on": q_on,
    }


# ---------------------------------------------------------------- 相圖 --------------
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
    matrix = np.array([
        [sensitive["expected_sensitive"], sensitive["expected_resistant"]],
        [resistant["expected_sensitive"], resistant["expected_resistant"]],
    ])
    assert matrix[1, 0] == 0  # 抗藥突變不可逆。
    return matrix


def compute_phase_field(
    config: Config, therapy_active: bool, counts, coordinates="linear"
) -> dict:
    """每個整數網格點的 E[(S',R') | 完整基準狀態, α] − (S,R)。

    display_u/v 是座標內的一步位移，不進行逐向量正規化。
    對數版先對原始計數求期望，再轉換；不是 E[log(count+1)]。
    """
    counts = np.asarray(counts)
    if (
        counts.ndim != 1
        or len(counts) < 2
        or np.any(counts < 0)
        or np.any(counts != np.rint(counts))
    ):
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
    return {
        "sensitive": sensitive,
        "resistant": resistant,
        "expected_sensitive": mean[..., 0],
        "expected_resistant": mean[..., 1],
        "delta_sensitive": delta[..., 0],
        "delta_resistant": delta[..., 1],
        "magnitude_cells": np.linalg.norm(delta, axis=-1),
        "x": position[..., 0],
        "y": position[..., 1],
        "display_u": displacement[..., 0],
        "display_v": displacement[..., 1],
        "display_magnitude": np.linalg.norm(displacement, axis=-1),
        "matrix": matrix,
        "coordinates": coordinates,
        "therapy_active": therapy_active,
    }


def verify_field(field, config):
    """逐點核對完整狀態期望；另驗證總數及零向量的語意。"""
    resistance = 1 << config.resistance_gene
    for index in np.ndindex(field["sensitive"].shape):
        s, r = int(field["sensitive"][index]), int(field["resistant"][index])
        population = {g: n for g, n in ((0, s), (resistance, r)) if n}
        mean = projected_offspring_mean(population, config, field["therapy_active"])
        np.testing.assert_allclose(
            [field["expected_sensitive"][index], field["expected_resistant"][index]],
            [mean["expected_sensitive"], mean["expected_resistant"]],
            rtol=1e-12,
            atol=1e-10,
        )
        assert np.isclose(
            field["expected_sensitive"][index] + field["expected_resistant"][index],
            mean["expected_total"],
            rtol=1e-12,
        )
    assert field["display_magnitude"][0, 0] == 0
    assert field["magnitude_cells"][0, 0] == 0
    assert np.all(field["expected_sensitive"][field["sensitive"] == 0] == 0)


def linear_flow_fields(config: Config) -> list[dict]:
    """矩陣流向圖的線性向量場：0–30,000 細胞的 13×13 網格，治療關閉與持續治療各一張，逐點核對。"""
    counts = np.rint(np.linspace(0, 30_000, 13)).astype(np.int64)
    fields = [
        compute_phase_field(config, active, counts, "linear")
        for active in (False, True)
    ]
    for field in fields:
        verify_field(field, config)
    return fields


def vector_field_grids(config: Config) -> tuple[dict, list[dict]]:
    """向量場資料檔用的網格：線性 19×19、對數 21×21，各分治療關閉與持續治療。"""
    grids = {
        "linear": np.rint(np.linspace(0, 25_000, 19)).astype(np.int64),
        "log": np.rint(np.power(10, np.linspace(0, 5.5, 21)) - 1).astype(np.int64),
    }
    fields = [
        compute_phase_field(config, active, counts, coordinates)
        for coordinates, counts in grids.items()
        for active in (False, True)
    ]
    # 目前範圍均在計算上限內，避免對已停止執行的狀態畫出更新箭頭。
    assert max(2 * grid[-1] for grid in grids.values()) < config.max_cells
    for field in fields:
        verify_field(field, config)
    return grids, fields


def state_features(population: dict[int, int], config: Config, therapy: bool) -> dict:
    """由完整基因型快照計算所有候選相座標。"""
    genotypes = list(population)
    counts = np.array(list(population.values()), float)
    total = counts.sum()
    p = counts / total
    resistance = 1 << config.resistance_gene
    passengers = [
        g
        for g in range(config.n_genes)
        if g
        not in config.driver_genes + config.suppressor_genes + (config.resistance_gene,)
    ]
    res = np.array([bool(g & resistance) for g in genotypes])
    onc = np.array([driver_count(g, config) >= config.driver_hits for g in genotypes])
    tsg_loss = np.array([suppressor_loss(g, config) for g in genotypes])
    tsg_any = np.array([
        any(suppressor_alleles(g, k, config) for k in config.suppressor_genes)
        for g in genotypes
    ])
    w_off = np.array([fitness(g, config, False) for g in genotypes])
    w_on = np.array([fitness(g, config, True) for g in genotypes])
    sensitive, resistant = counts[~res].sum(), counts[res].sum()
    mean_off = (p * w_off).sum()
    return {
        "log_N": np.log10(total + 1),
        "log_S": np.log10(sensitive + 1),
        "log_R": np.log10(resistant + 1),
        "frac_R": p[res].sum(),
        "logit_R": float(logit_rho(resistant, sensitive)),
        "frac_onc": p[onc].sum(),
        "mean_drivers": (p * [driver_count(g, config) for g in genotypes]).sum(),
        "frac_tsg_anyhit": p[tsg_any].sum(),
        "frac_tsg_loss": p[tsg_loss].sum(),
        "frac_onc_and_tsg": p[onc & tsg_loss].sum(),
        "mean_passengers": (
            (p * [sum(bool(g & (1 << k)) for k in passengers) for g in genotypes]).sum()
        ),
        "log_richness": np.log10(len(genotypes)),
        "shannon": float(-(p * np.log(p)).sum()),
        "simpson": float(1 - (p**2).sum()),
        "mean_w_off": mean_off,
        "var_w_off": (p * (w_off - mean_off) ** 2).sum(),
        "mean_w_on": (p * w_on).sum(),
        "mean_w_now": (p * (w_on if therapy else w_off)).sum(),
        "therapy_on": float(therapy),
    }


def trial_feature_rows(config: Config, results: list[Result]) -> list[dict]:
    """每個 (試驗, 世代) 一列，供相圖與復發分析使用。

    growth/dlogit 的 _next 是實際一步，_expected 是期望漂移；兩者是動力學目標，不是座標。
    """
    rows = []
    for replicate, result in enumerate(results):
        start = result.therapy_generation
        for generation, population in enumerate(result.populations):
            if not population:
                continue
            therapy = start is not None and generation >= start
            row = {
                "replicate": replicate,
                "generation": generation,
                **state_features(population, config, therapy),
                "gens_since_therapy": float(generation - start) if therapy else 0.0,
            }
            total = sum(population.values())
            resistance = 1 << config.resistance_gene
            sensitive = sum(n for g, n in population.items() if not g & resistance)
            following = (
                result.populations[generation + 1]
                if generation + 1 < len(result.populations)
                else None
            )
            # 一步拆成「期望漂移」與「隨機噪音」：期望由完整基因型組成精確計算。
            # 負荷方向用 ln 成長率；組成方向用相圖的經驗 logit（log10 勝算）。
            expected = expected_next(population, config, therapy)
            here = float(logit_rho(total - sensitive, sensitive))
            row["sensitive"] = sensitive
            row["growth_expected"] = np.log(expected["N"] / total)
            row["dlogit_expected"] = (
                float(logit_rho(expected["R"], expected["S"])) - here
            )
            row["growth_next"] = (
                np.log(sum(following.values()) / total) if following else np.nan
            )
            if following:
                next_sensitive = sum(
                    n for g, n in following.items() if not g & resistance
                )
                row["dlogit_next"] = (
                    float(
                        logit_rho(
                            sum(following.values()) - next_sensitive, next_sensitive
                        )
                    )
                    - here
                )
            else:
                row["dlogit_next"] = np.nan
            rows.append(row)
    return rows


def conditional_steps(result: Result):
    """每個可推進世代的原始計數期望，以及實際下一代計數。"""
    rows = []
    for g, population in enumerate(result.populations[:-1]):
        active = (
            result.therapy_generation is not None and g >= result.therapy_generation
        )
        current, following = result.history[g : g + 2]
        mean = projected_offspring_mean(population, result.config, active)
        rows.append({
            "generation": g,
            "therapy_active": active,
            "sensitive": current["sensitive"],
            "resistant": current["resistant"],
            **mean,
            "expected_delta_sensitive": (
                mean["expected_sensitive"] - current["sensitive"]
            ),
            "expected_delta_resistant": (
                mean["expected_resistant"] - current["resistant"]
            ),
            "observed_next_sensitive": following["sensitive"],
            "observed_next_resistant": following["resistant"],
        })
    return rows


def event_generations(states, event):
    start, relapse = event["therapy_generation"], event["relapse_generation"]
    points = {"初始": 0}
    if start is not None:
        points["治療啟動"] = start
        after = states[start + 1 :]
        if after:
            points["負荷谷底"] = min(after, key=lambda r: r["total"])["generation"]
    if relapse is not None:
        points["首次復發"] = relapse
    return points


def choose_exemplar(events):
    """先選最接近復發中位代數，再選谷底最接近中位數的實際試驗。"""
    eligible = [
        e for e in events.values() if e["relapse_generation"] is not None and e["nadir"]
    ]
    if not eligible:
        return min(events)
    relapse = np.median([e["relapse_generation"] for e in eligible])
    nadir = np.median([e["nadir"] for e in eligible])
    return min(
        eligible,
        key=lambda e: (
            abs(e["relapse_generation"] - relapse),
            abs(np.log(e["nadir"] / nadir)),
            e["replicate"],
        ),
    )["replicate"]


MIN_TRANSITIONS, MIN_TRIALS = 5, 2


def theory_point(pair, x, y, config, therapy, size):
    """格點 → 完整族群（隱藏組成依截面固定）→ 一步期望 → 相座標位移。超出可行域時回傳 None。"""
    odds = 10**y
    if pair == "logN":
        total = 10**x - 1
        sensitive = (total + 0.5 - 0.5 * odds) / (1 + odds)
        if sensitive < 0 or sensitive > total:
            return None
        population = {0: sensitive, 1 << config.resistance_gene: total - sensitive}
    else:
        w0, w1 = fitness(0, config, False), fitness(
            1 << config.driver_genes[0], config, False
        )
        f = (x - w0) / (w1 - w0)
        sensitive = (size + 0.5 - 0.5 * odds) / (1 + odds)
        if not 0 <= f <= 1 or sensitive < 0 or sensitive > size:
            return None
        resistant, onc, res = (
            size - sensitive,
            1 << config.driver_genes[0],
            1 << config.resistance_gene,
        )
        population = {
            0: sensitive * (1 - f),
            onc: sensitive * f,
            res: resistant * (1 - f),
            res | onc: resistant * f,
        }
    nxt = expected_next(population, config, therapy)
    new_x = np.log10(nxt["N"] + 1) if pair == "logN" else nxt["w_off"]
    return new_x - x, float(logit_rho(nxt["R"], nxt["S"])) - y


def empirical_steps(rows, pair):
    """同一試驗相鄰世代的位移；治療狀態以起點世代為準（治療從該轉移起作用）。"""
    key = "log_N" if pair == "logN" else "mean_w_off"
    by_trial = {}
    for row in rows:
        by_trial.setdefault(row["replicate"], {})[row["generation"]] = row
    steps = []
    for replicate, gens in by_trial.items():
        for g, row in gens.items():
            if g + 1 in gens:
                nxt = gens[g + 1]
                steps.append({
                    "replicate": replicate,
                    "generation": g,
                    "therapy": bool(row["therapy_on"]),
                    "x": row[key],
                    "y": row["logit_R"],
                    "dx": nxt[key] - row[key],
                    "dy": nxt["logit_R"] - row["logit_R"],
                })
    return steps


def build_fields(rows, config, pair):
    steps = empirical_steps(rows, pair)
    xs = np.array([s["x"] for s in steps])
    ys = np.array([s["y"] for s in steps])
    pad_x = 0.04 * (xs.max() - xs.min())
    x_edges = np.linspace(xs.min() - pad_x, xs.max() + pad_x, 13)
    y_edges = np.linspace(
        np.floor(ys.min() * 2) / 2 - 0.25, np.ceil(ys.max() * 2) / 2 + 0.25, 15
    )
    xc, yc = (x_edges[:-1] + x_edges[1:]) / 2, (y_edges[:-1] + y_edges[1:]) / 2
    fields = []
    for therapy in (False, True):
        for i, x in enumerate(xc):
            for j, y in enumerate(yc):
                step = theory_point(
                    pair, x, y, config, therapy, config.therapy_start_size
                )
                if step is not None:
                    fields.append({
                        "pair": pair,
                        "kind": "theory",
                        "therapy": therapy,
                        "x": x,
                        "y": y,
                        "dx": step[0],
                        "dy": step[1],
                        "transitions": "",
                        "trials": "",
                    })
        chosen = [s for s in steps if s["therapy"] == therapy]
        ix = np.digitize([s["x"] for s in chosen], x_edges) - 1
        iy = np.digitize([s["y"] for s in chosen], y_edges) - 1
        for i in range(len(xc)):
            for j in range(len(yc)):
                cell = [s for s, a, b in zip(chosen, ix, iy) if a == i and b == j]
                if (
                    len(cell) >= MIN_TRANSITIONS
                    and len({s["replicate"] for s in cell}) >= MIN_TRIALS
                ):
                    fields.append({
                        "pair": pair,
                        "kind": "empirical",
                        "therapy": therapy,
                        "x": np.mean([s["x"] for s in cell]),
                        "y": np.mean([s["y"] for s in cell]),
                        "dx": np.mean([s["dx"] for s in cell]),
                        "dy": np.mean([s["dy"] for s in cell]),
                        "transitions": len(cell),
                        "trials": len({s["replicate"] for s in cell}),
                    })
    return fields, steps, x_edges, y_edges


# ---------------------------------------------------------------- 基因型組成 --------
def gene_label(bit: int, config: Config) -> str:
    """位元 → 基因名稱；基因以 1 起算，TSG 第二等位基因位於 n_genes 之後。"""
    if bit in config.driver_genes:
        return f"ONC{config.driver_genes.index(bit) + 1}"
    if bit == config.resistance_gene:
        return "RES"
    if bit in config.suppressor_genes:
        return f"TSG{config.suppressor_genes.index(bit) + 1}a"
    if bit >= config.n_genes:
        return f"TSG{bit - config.n_genes + 1}b"
    return f"P{bit + 1}"


def genotype_label(genotype: int, config: Config) -> str:
    bits = [
        b
        for b in range(config.n_genes + len(config.suppressor_genes))
        if genotype & (1 << b)
    ]
    functional = [
        gene_label(b, config) for b in bits if not gene_label(b, config).startswith("P")
    ]
    passengers = [
        gene_label(b, config) for b in bits if gene_label(b, config).startswith("P")
    ]
    return "+".join(functional + passengers) or "WT"


TOP_K = 10


def family(genotype: int, config: Config) -> str:
    res = bool(genotype & (1 << config.resistance_gene))
    onc = driver_count(genotype, config) >= config.driver_hits
    tsg = suppressor_loss(genotype, config)
    name = "_".join(n for n, flag in (("onc", onc), ("tsg", tsg)) if flag) or "wt"
    return ("res_" + name if name != "wt" else "res") if res else name


def is_subset(a: int, b: int) -> bool:
    return a & b == a


def muller_layout(result: Result, k: int = TOP_K):
    """選 Top-k 基因型並建立巢狀樹；未選者併入最大的已選子集祖先，因此頻率總和恆為 1。"""
    peak = {}
    for population in result.populations:
        total = sum(population.values())
        for g, n in population.items():
            peak[g] = max(peak.get(g, 0), n / total)
    # 基因型 0（起始狀態）永遠是根，確保每個細胞都有可併入的祖先。
    shown = set(sorted(peak, key=lambda g: (-peak[g], g))[:k]) | {0}
    shown = sorted(shown, key=lambda g: (bin(g).count("1"), g))

    def ancestor(g):
        """已呈現、位元為 g 子集的最大基因型（不含 g 本身）。"""
        candidates = [a for a in shown if a != g and is_subset(a, g)]
        return max(candidates, key=lambda a: (bin(a).count("1"), -a))

    parent = {g: None if g == 0 else ancestor(g) for g in shown}
    children = {g: [c for c in shown if parent[c] == g] for g in shown}
    own = {g: np.zeros(len(result.populations)) for g in shown}
    for t, population in enumerate(result.populations):
        total = sum(population.values())
        for g, n in population.items():
            own[g if g in own else ancestor(g)][t] += n / total
    assert np.allclose(sum(own.values()), 1)
    return shown, parent, children, own


def muller_bands(shown, children, own):
    """遞迴：父基因型帶寬 = 自身 + 所有後代；子帶置於父帶中央，父色在上下兩側可見。"""
    subtree = {}

    def total(g):
        if g not in subtree:
            subtree[g] = own[g] + sum(total(c) for c in children[g])
        return subtree[g]

    total(0)
    bands = {}

    def place(g, bottom):
        bands[g] = (bottom, bottom + subtree[g])
        cursor = bottom + own[g] / 2
        for c in children[g]:
            place(c, cursor)
            cursor = cursor + subtree[c]

    place(0, np.zeros_like(own[0]))
    return bands


def depth(g, parent):
    d = 0
    while parent[g] is not None:
        g = parent[g]
        d += 1
    return d


def gene_frequencies(population, config):
    total = sum(population.values())
    freq = {}
    for gene in config.driver_genes + (config.resistance_gene,):
        freq[gene] = sum(n for g, n in population.items() if g & (1 << gene)) / total
    for gene in config.suppressor_genes:
        freq[f"{gene}_any"] = (
            sum(n for g, n in population.items() if suppressor_alleles(g, gene, config))
            / total
        )
        freq[f"{gene}_bi"] = (
            sum(
                n
                for g, n in population.items()
                if suppressor_alleles(g, gene, config) == 2
            )
            / total
        )
    passengers = [
        g
        for g in range(config.n_genes)
        if g
        not in config.driver_genes + config.suppressor_genes + (config.resistance_gene,)
    ]
    freq["passenger_mean"] = np.mean([
        sum(n for g, n in population.items() if g & (1 << p)) / total
        for p in passengers
    ])
    return freq


def gene_frequency_series(config: Config, results: list[Result]):
    """以治療開始對齊；只取所有試驗都存活且有觀察的相對世代。

    回傳相對世代、每個序列的 (試驗 × 世代) 矩陣，以及中位數與 10／90 百分位摘要。
    """
    starts = [r.therapy_generation for r in results]
    lo = -min(starts)
    hi = min(len(r.populations) - 1 - s for r, s in zip(results, starts))
    rel = np.arange(lo, hi + 1)
    series = {}
    for r, s in zip(results, starts):
        for x in rel:
            for key, value in gene_frequencies(r.populations[s + x], config).items():
                series.setdefault(key, []).append(value)
    data = {
        key: np.array(v).reshape(len(results), len(rel)) for key, v in series.items()
    }
    summary = {
        key: {
            "relative_generation": rel.tolist(),
            "median": np.median(v, 0).tolist(),
            "p10": np.percentile(v, 10, 0).tolist(),
            "p90": np.percentile(v, 90, 0).tolist(),
        }
        for key, v in data.items()
    }
    return rel, data, summary
