"""Extend presentation data with omitted panels; verify regenerated trials.

Run once after prepare_results.py. Seed/config come from recorded metadata.
Core simulation is imported unchanged; saved CSV totals must match every trial.
"""
import csv
import json
from dataclasses import replace
from pathlib import Path
import sys
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from modeling_cancer import Config, simulate, suppressor_alleles


def read(name):
    return list(csv.DictReader((Path('outputs') / name).open()))


def main():
    root = Path(__file__).resolve().parent
    destination = root / 'data/results.json'
    data = json.loads(destination.read_text())
    meta = json.loads(Path('outputs/run_metadata.json').read_text())
    config_values = meta['config'].copy()
    for key in ('driver_genes', 'suppressor_genes'):
        config_values[key] = tuple(config_values[key])
    config = Config(**config_values)
    trials = [simulate(config, seed) for seed in meta['therapy_seeds']]
    raw = read('therapy_trajectories.csv')
    reference = {(int(r['replicate']), int(r['generation'])): r for r in raw}
    length = max(len(t.history) for t in trials)
    frequencies = np.full((len(trials), config.n_genes, length), np.nan)
    for rep, trial in enumerate(trials):
        for generation, (row, pop) in enumerate(zip(trial.history, trial.populations)):
            assert int(reference[rep, generation]['total']) == row['total']
            assert int(reference[rep, generation]['resistant']) == row['resistant']
            if row['total']:
                for gene in range(config.n_genes):
                    frequencies[rep, gene, generation] = sum(n for g, n in pop.items()
                        if (suppressor_alleles(g, gene, config) > 0 if gene in config.suppressor_genes
                            else bool(g & (1 << gene)))) / row['total']
    # Keep the shared complete-observation window for interpretable heatmaps.
    freq = frequencies[:, :, :31]
    assert np.isfinite(freq).all()
    data['gene_heatmaps'] = {'generation': list(range(31)), 'mean': freq.mean(axis=0).tolist(),
                             'sd': freq.std(axis=0, ddof=1).tolist(),
                             'gene_labels': ['G'+str(g+1)+(' driver' if g in config.driver_genes else
                               ' resistance' if g == config.resistance_gene else ' TSG any hit' if g in config.suppressor_genes else '')
                               for g in range(config.n_genes)]}
    trajectory = read('trajectory.csv')
    data['availability'] = {'generation': [int(r['generation']) for r in trajectory],
                            'observed_n': [int(r['total_n']) for r in trajectory],
                            'alive_n': [int(r['resistant_percent_n']) for r in trajectory],
                            'therapy_started': [sum(t.therapy_generation is not None and t.therapy_generation <= g
                                                  for t in trials) for g in range(length)]}
    for key in ('p10','median','p90'):
        data['therapy']['effective_genotypes_'+key] = [float(r['effective_genotypes_'+key]) for r in trajectory[:31]]
    hits, il11 = read('one_hit_vs_two_hit.csv'), read('il11_comparison.csv')
    data['endpoints'] = {}
    for rows, category, experiment in ((hits,'scenario','hits'),(il11,'condition','il11')):
        data['endpoints'][experiment] = {label: [float(r['total']) for r in rows if r[category] == label and int(r['generation']) == 18]
                                         for label in sorted({r[category] for r in rows})}
    # The omitted no-TSG-advantage control belongs to 02_scenario_comparisons.
    control_config = replace(config, generations=18, max_cells=10**12,
                             suppressor_death_reduction=0, suppressor_triple_bonus=0)
    control = [simulate(control_config, meta['seed']+rep, therapy=False) for rep in range(meta['replicates'])]
    reference_control = [r for r in read('comparisons.csv') if r['experiment']=='growth' and r['condition']=='No TSG advantage']
    assert len(reference_control) == meta['replicates']
    for trial, row in zip(control, reference_control):
        assert trial.history[-1]['total'] == int(row['total'])
    values = np.array([[t.history[min(g,len(t.history)-1)]['total'] for g in range(19)] for t in control])
    data['no_tsg_advantage'] = {'generation':list(range(19)), **{f'total_{key}':np.percentile(values,p,axis=0).tolist()
                                                            for key,p in [('p10',10),('median',50),('p90',90)]}}
    data['il11_curve'] = {'fraction_percent': [i/10 for i in range(101)],
                          'benefit_percentage_points': [100*config.il11_max_death_reduction*(i/1000)/(config.il11_half_fraction+i/1000)
                                                        for i in range(101)]}
    data['verification'] = {'therapy_trials_reproduced': len(trials), 'therapy_csv_rows_verified': len(reference),
                            'control_trials_reproduced':len(control), 'heatmap_window':'0-30; all 30 alive and observed'}
    destination.write_text(json.dumps(data,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(data['verification'],indent=2))


if __name__ == '__main__':
    main()
