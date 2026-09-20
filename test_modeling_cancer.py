"""驗證抑癌基因等位基因、突變守恆與生長機率。"""
import unittest
from dataclasses import replace
import numpy as np
from modeling_cancer import (Config, suppressor_alleles, suppressor_loss,
                             growth_probabilities, mutate_offspring, simulate)


class SuppressorTests(unittest.TestCase):
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
        self.assertLess(growth_probabilities(first | second, c)[0], growth_probabilities(0, c)[0])

    def test_mutation_conservation_and_independent_alleles(self):
        c = Config(mutation_rate=1)
        full = (1 << (c.n_genes + len(c.suppressor_genes))) - 1
        self.assertEqual(mutate_offspring({0: 100}, c, np.random.default_rng(0)), {full: 100})
        c = replace(c, mutation_rate=0)
        self.assertEqual(mutate_offspring({0: 100, 16: 20}, c, np.random.default_rng(0)), {0: 100, 16: 20})
        # 只留一個 TSG 的兩個等位基因可突變；檢查 0/1/2 擊的二項機率。
        c = replace(c, mutation_rate=0.2)
        base = full ^ (1 << 4) ^ (1 << c.n_genes)
        pop = mutate_offspring({base: 100000}, c, np.random.default_rng(8))
        counts = [sum(n for g, n in pop.items() if suppressor_alleles(g, 4, c) == k) for k in range(3)]
        np.testing.assert_allclose(np.array(counts) / 100000, [0.64, 0.32, 0.04], atol=0.006)
        self.assertEqual(sum(pop.values()), 100000)

    def test_reproducibility_and_probability_validation(self):
        c = Config(generations=4)
        a = simulate(c, seed=9)
        self.assertEqual(a.history, simulate(c, seed=9).history)
        for row, pop in zip(a.history, a.populations):
            self.assertEqual(row['total'], sum(pop.values()))
            self.assertLessEqual(row['suppressor_biallelic_fraction'], row['suppressor_fraction'])
        for invalid in (replace(c, suppressor_hits=3), replace(c, suppressor_genes=(0,)),
                        replace(c, suppressor_triple_bonus=0.99)):
            with self.assertRaises(ValueError):
                invalid.validate()


if __name__ == '__main__':
    unittest.main()
