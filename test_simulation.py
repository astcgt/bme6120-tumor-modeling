"""驗證模型：抑癌基因等位基因、突變守恆、生長機率，以及多次試驗統計。

Model tests: tumor-suppressor alleles, mutation conservation, growth probabilities and multi-trial statistics.
"""

import unittest
from dataclasses import replace
import numpy as np
from experiments import therapy_ensemble
from simulation import (
    Config,
    suppressor_alleles,
    suppressor_loss,
    growth_probabilities,
    mutate_offspring,
    simulate,
)


class SuppressorTests(unittest.TestCase):
    # 驗證多次試驗統計會區分計算上限後缺失、滅絕後零值及仍在觀察的族群。
    # Multi-trial statistics must distinguish missing after the cap, zero after extinction, and populations still observed.
    def test_ensemble_censoring_and_extinction(self):
        from simulation import Result, summarize

        c = Config()

        # 上限後缺失、滅絕後零值、仍觀察中的正常族群應分開處理。
        # Missing after the cap, zero after extinction and normal observed populations are handled separately.
        # 依指定族群快照與停止原因建立測試用 Result，並自動產生逐代摘要。
        # Build a test Result from population snapshots and a stop reason, generating the per-generation summaries.
        def make(pops, reason):
            return Result(
                c,
                1,
                [summarize(p, g, c) for g, p in enumerate(pops)],
                pops,
                None,
                reason,
            )

        trials = [
            make([{0: 100}], "cell_limit"),
            make([{0: 100}, {}], "extinction"),
            make([{0: 100}, {0: 150}, {0: 200}], "generation_limit"),
        ]
        summary, _ = therapy_ensemble(trials)
        self.assertEqual(summary[2]["total_n"], 2)
        self.assertEqual(summary[2]["total_mean"], 100)
        self.assertEqual(summary[2]["resistant_percent_n"], 1)
        self.assertTrue(np.isnan(summary[2]["resistant_percent_sd"]))

    # 驗證 TSG 雙擊必須發生於同一基因，並檢查 one-hit 設定與生長優勢。
    # A TSG two-hit must be in the same gene; also check the one-hit setting and the growth advantage.
    def test_same_gene_two_hits(self):
        c = Config()
        first = 1 << c.suppressor_genes[0]
        second = 1 << c.n_genes
        other_gene = 1 << c.suppressor_genes[1]
        self.assertEqual(c.n_genes, 20)
        self.assertEqual(suppressor_alleles(first | second, 4, c), 2)
        self.assertFalse(suppressor_loss(first, c))
        self.assertFalse(suppressor_loss(first | other_gene, c))
        self.assertTrue(suppressor_loss(first | second, c))
        self.assertTrue(suppressor_loss(second, replace(c, suppressor_hits=1)))
        self.assertEqual(growth_probabilities(first, c), growth_probabilities(0, c))
        self.assertLess(
            growth_probabilities(first | second, c)[0], growth_probabilities(0, c)[0]
        )

    # 驗證突變保留細胞總數、零與必然突變邊界，以及 TSG 等位基因的獨立突變比例。
    # Mutation keeps the total; check the zero and certain-mutation limits and the independent mutation rate of TSG alleles.
    def test_mutation_conservation_and_independent_alleles(self):
        c = Config(mutation_rate=1)
        full = (1 << (c.n_genes + len(c.suppressor_genes))) - 1
        self.assertEqual(
            mutate_offspring({0: 100}, c, np.random.default_rng(0)), {full: 100}
        )
        c = replace(c, mutation_rate=0)
        self.assertEqual(
            mutate_offspring({0: 100, 16: 20}, c, np.random.default_rng(0)),
            {0: 100, 16: 20},
        )
        # 只留一個 TSG 的兩個等位基因可突變；檢查 0/1/2 擊的二項機率。
        # Only the two alleles of one TSG can mutate; check the binomial probabilities of 0/1/2 hits.
        c = replace(c, mutation_rate=0.2)
        base = full ^ (1 << 4) ^ (1 << c.n_genes)
        pop = mutate_offspring({base: 100000}, c, np.random.default_rng(8))
        counts = [
            sum(n for g, n in pop.items() if suppressor_alleles(g, 4, c) == k)
            for k in range(3)
        ]
        np.testing.assert_allclose(
            np.array(counts) / 100000, [0.64, 0.32, 0.04], atol=0.006
        )
        self.assertEqual(sum(pop.values()), 100000)

    # 驗證固定種子的可重現性、族群摘要一致性及無效參數的錯誤處理。
    # Fixed seeds are reproducible, population summaries agree, and invalid parameters raise errors.
    def test_reproducibility_and_probability_validation(self):
        c = Config(generations=4)
        a = simulate(c, seed=9)
        self.assertEqual(a.history, simulate(c, seed=9).history)
        for row, pop in zip(a.history, a.populations):
            self.assertEqual(row["total"], sum(pop.values()))
            self.assertLessEqual(
                row["suppressor_biallelic_fraction"], row["suppressor_fraction"]
            )
        for invalid in (
            replace(c, suppressor_hits=3),
            replace(c, suppressor_genes=(0,)),
            replace(c, suppressor_triple_bonus=0.99),
        ):
            with self.assertRaises(ValueError):
                invalid.validate()


if __name__ == "__main__":
    unittest.main()
