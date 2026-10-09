"""Prepare presentation chart data from recorded simulation outputs.

Run with the project-local Python. No simulation is rerun by this helper.
"""
import argparse
import csv
import json
from collections import Counter
from pathlib import Path
import numpy as np


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', type=Path, default=Path('outputs'))
    args = parser.parse_args()
    source = args.input

    def read(name):
        with (source / name).open(newline='') as handle:
            return list(csv.DictReader(handle))

    meta = json.loads((source / 'run_metadata.json').read_text())
    n = meta['replicates']
    outcomes = read('therapy_replicates.csv')
    assert len(outcomes) == n
    assert Counter(r['status'] for r in outcomes) == meta['therapy_outcome_counts']
    trajectory = read('trajectory.csv')
    early = [r for r in trajectory if int(r['generation']) <= 30]
    assert len(early) == 31
    assert all(int(r['total_n']) == n and int(r['resistant_percent_n']) == n for r in early)
    therapy = {'generation': [int(r['generation']) for r in early]}
    for key in ('total', 'sensitive', 'resistant', 'resistant_percent'):
        for stat in ('p10', 'median', 'p90'):
            therapy[f'{key}_{stat}'] = [float(r[f'{key}_{stat}']) for r in early]
    therapy['outcome_counts'] = dict(Counter(r['status'] for r in outcomes))
    for key in ('therapy_generation', 'nadir', 'relapse_generation'):
        values = [float(r[key]) for r in outcomes if r[key]]
        therapy[f'{key}_quantiles'] = np.percentile(values, [10, 50, 90]).tolist()
    relapse = Counter(int(r['relapse_generation']) for r in outcomes if r['relapse_generation'])
    generations = list(range(min(relapse), max(relapse)+1))
    therapy['relapse_histogram'] = {'generation': generations, 'count': [relapse.get(g, 0) for g in generations]}

    def group_curves(rows, condition_key, fields):
        curves = {}
        for condition in sorted({r[condition_key] for r in rows}):
            rr = [r for r in rows if r[condition_key] == condition]
            assert len({r['replicate'] for r in rr}) == n
            generations = sorted({int(r['generation']) for r in rr})
            group = {'generation': generations}
            for field in fields:
                for stat, p in [('p10', 10), ('median', 50), ('p90', 90)]:
                    values = []
                    for g in generations:
                        a = [float(r[field]) for r in rr if int(r['generation']) == g]
                        assert len(a) == n
                        values.append(float(np.percentile(a, p)))
                    group[f'{field}_{stat}'] = values
            curves[condition] = group
        return curves

    hits = group_curves(read('one_hit_vs_two_hit.csv'), 'scenario', ['total', 'suppressor_fraction'])
    il11 = group_curves(read('il11_comparison.csv'), 'condition', ['total', 'il11_percent', 'shared_death_reduction'])
    # Independently reconcile plotted endpoint medians with recorded summaries.
    for r in read('one_hit_vs_two_hit_summary.csv'):
        assert np.isclose(hits[r['scenario']]['total_median'][-1], float(r['median_cells']))
    for r in read('il11_comparison_summary.csv'):
        assert np.isclose(il11[r['condition']]['total_median'][-1], float(r['median_total']))
    death = []
    comparisons = read('comparisons.csv')
    for condition in sorted({r['condition'] for r in comparisons if r['experiment'] == 'death_rate'}):
        rr = [r for r in comparisons if r['experiment'] == 'death_rate' and r['condition'] == condition]
        assert len(rr) == n
        observed = [r for r in rr if r['shannon']]
        values = [float(r['shannon']) for r in observed]
        death.append({'death': float(condition), 'reached': len(observed), 'trials': n,
                      'shannon': values, 'median_shannon': float(np.median(values)),
                      'median_generation_reached': float(np.median([float(r['generation']) for r in observed]))})
    result = {'replicates': n, 'seed': meta['seed'], 'config': meta['config'], 'therapy': therapy,
              'hits': hits, 'il11': il11, 'death_rate': death,
              'one_hit_median_ratio': hits['TSG one-hit']['total_median'][-1] / hits['TSG two-hit']['total_median'][-1],
              'il11_median_ratio': il11['IL11 ON']['total_median'][-1] / il11['IL11 OFF']['total_median'][-1],
              'source_files': ['run_metadata.json','trajectory.csv','therapy_replicates.csv',
                               'one_hit_vs_two_hit.csv','il11_comparison.csv','comparisons.csv']}
    destination = Path(__file__).resolve().parent / 'data' / 'results.json'
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(json.dumps(result, ensure_ascii=False, indent=2) + '\n')
    print(f'Verified {n} trials per condition and saved {destination}')


if __name__ == '__main__':
    main()
