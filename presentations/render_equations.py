"""Render LaTeX math syntax with Matplotlib Computer Modern math fonts.

No external TeX installation or conda changes are required. The presentation
stores the exact source in speaker notes and regenerates 4x-resolution assets.
"""
import json
import sys
from pathlib import Path
import os
os.environ.setdefault('MPLCONFIGDIR', str(Path(__file__).resolve().parents[1] / '.mplconfig'))
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from matplotlib.mathtext import MathTextParser


def main():
    specification = Path(sys.argv[1])
    formulas = json.loads(specification.read_text())
    plt.rcParams.update({'mathtext.fontset': 'cm', 'font.family': 'serif'})
    parser = MathTextParser('path')
    for item in formulas:
        item['render_latex'] = item['latex'].replace(r'\bigl(', r'\left(').replace(r'\bigr)', r'\right)')
        parser.parse('$' + item['render_latex'] + '$')
    for item in formulas:
        fig = plt.figure(figsize=(16, 2), dpi=400)
        fig.text(0.01, 0.35, '$' + item['render_latex'] + '$', fontsize=32,
                 color='#236B8E', va='center')
        destination = Path(item['file'])
        destination.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(destination, transparent=True, bbox_inches='tight', pad_inches=0.035)
        plt.close(fig)
    print(f'Rendered {len(formulas)} Computer Modern equations')


if __name__ == '__main__':
    main()
