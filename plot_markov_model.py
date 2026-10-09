"""Generate publication-style mathematical diagrams using the project's defaults.

Run: .venv/bin/python plot_markov_model.py
The finite count chain and two-type matrices are explicitly reduced examples.
"""
from pathlib import Path
import csv
import json

from modeling_cancer import Config, growth_probabilities
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import numpy as np


OUT = Path(__file__).resolve().parent / "outputs" / "mathematical_model"
INK, BLUE, RED, GRAY = "#182634", "#21658B", "#A54C38", "#66717A"


def finite_count_matrix(death, triple, threshold=4):
    """Exact stopped one-type chain: 0,...,threshold-1, and overflow dagger.

    No mutations or therapy. Overflow is absorbing by construction, rather
    than a claim that grouping all large populations preserves Markovity.
    """
    offspring = np.array([death, 0., (1-death)*(1-triple), (1-death)*triple])
    matrix = np.zeros((threshold+1, threshold+1))
    for n in range(threshold):
        law = np.array([1.])
        for _ in range(n):
            law = np.convolve(law, offspring)
        for j, prob in enumerate(law):
            matrix[n, min(j, threshold)] += prob
    matrix[threshold, threshold] = 1
    return matrix, offspring


def panel(ax, letter, title):
    ax.set(xlim=(0, 1), ylim=(0, 1))
    ax.axis("off")
    ax.text(0, 1, letter, fontsize=17, weight="bold", va="top", color=BLUE)
    ax.text(.055, .995, title, fontsize=13, weight="bold", va="top", color=INK)


def label(ax, x, y, text, size=12, **kw):
    return ax.text(x, y, text, fontsize=size, color=kw.pop("color", INK), va="top", **kw)


def node(ax, xy, text, width=.25, height=.115, color=BLUE):
    x, y = xy
    patch = FancyBboxPatch((x-width/2, y-height/2), width, height,
                          boxstyle="round,pad=0.008,rounding_size=0.015",
                          facecolor="#F3F6F8", edgecolor=color, linewidth=1.1)
    ax.add_patch(patch)
    ax.text(x, y, text, ha="center", va="center", fontsize=11, color=INK)


def arrow(ax, a, b, text=None, at=None, rad=0, color=GRAY):
    ax.add_patch(FancyArrowPatch(a, b, arrowstyle="-|>", mutation_scale=12,
                                linewidth=1.2, color=color,
                                connectionstyle=f"arc3,rad={rad}"))
    if text:
        x, y = at if at else ((a[0]+b[0])/2, (a[1]+b[1])/2)
        ax.text(x, y, text, fontsize=11, ha="center", va="center", color=INK,
                bbox=dict(facecolor="white", edgecolor="none", pad=2))


def matrix_display(ax, values, x, y, width, height, rows=None, cols=None, fmt=".6g"):
    """Draw an actual numeric matrix, with conventional row-to-column orientation."""
    nr, nc = values.shape
    cw, ch = width/nc, height/nr
    for i in range(nr):
        for j in range(nc):
            ax.text(x+(j+.5)*cw, y-(i+.5)*ch, format(values[i,j], fmt),
                    ha="center", va="center", fontsize=12,
                    color=BLUE if values[i,j] else GRAY)
    for bx, inward in ((x, .015), (x+width, -.015)):
        ax.plot([bx+inward, bx, bx, bx+inward], [y, y, y-height, y-height],
                color=INK, lw=1.3)
    if rows:
        for i, name in enumerate(rows):
            ax.text(x-.035, y-(i+.5)*ch, name, ha="right", va="center", fontsize=12)
    if cols:
        for j, name in enumerate(cols):
            ax.text(x+(j+.5)*cw, y+.04, name, ha="center", va="bottom", fontsize=12)


def save(fig, name):
    for suffix in ("png", "svg"):
        fig.savefig(OUT / f"{name}.{suffix}", dpi=300, facecolor="white")
    plt.close(fig)


def markov_figure(config, matrix, offspring):
    fig = plt.figure(figsize=(15, 10.8))
    fig.text(.04, .965, "Markov formulation of stochastic tumor evolution", fontsize=22, weight="bold", color=INK)
    fig.text(.04, .93, "Exact population transition kernel and a fully specified finite numerical example", fontsize=12, color=GRAY)
    ax = fig.add_axes([.04, .52, .92, .36])
    panel(ax, "A", "Full genotype-count process")
    label(ax, .015, .84, r"$S_t=(\mathbf{N}_t,A_t),\quad \mathbf{N}_t\in\mathbb{N}_0^{2^{22}},\quad A_t\in\{0,1\}$", 17)
    label(ax, .015, .65, r"$\alpha=A_t\vee\mathbf{1}\{\mathrm{therapy\ enabled},\ \sum_g N_{g,t}\geq20{,}000\}$", 15)
    label(ax, .015, .47, r"$b(\mathbf{N})=0.06\,f/(0.02+f),\qquad f=\sum_{g:\,g_6=1}N_g/\sum_gN_g$", 15)
    label(ax, .015, .29, r"$\delta_g=1-(1-d_g(\mathbf{N}))(1-\alpha k_g),\qquad u_g(\mathbf{z})=\sum_hQ_{gh}z_h$", 15)
    label(ax, .015, .10, "IL11 OFF: b = 0. Empty population: b = 0. Fixed Config; kernel applies before stopping.", 11, color=GRAY)
    label(ax, .57, .82, r"$F_g(\mathbf{z}\mid S)=\delta_g+(1-\delta_g)\left[(1-q_g)u_g^2+q_gu_g^3\right]$", 15)
    label(ax, .57, .61, r"$P_{S,S'}=\mathbf{1}\{A'=\alpha\}\,[\mathbf{z}^{\mathbf{N}'}]\prod_gF_g(\mathbf{z}\mid S)^{N_g}$", 16)
    label(ax, .57, .39, r"$Q_{gh}=\mu^{|h|-|g|}(1-\mu)^{22-|h|}\quad\mathrm{if}\ g\subseteq h;\quad 0\ \mathrm{otherwise}$", 14)
    label(ax, .57, .19, "Coefficient extraction gives the exact next-population probability.\nQ is inherited, irreversible, independent-bit mutation; mu = 0.001.", 11)

    ax = fig.add_axes([.04, .11, .40, .36])
    panel(ax, "B", "Single-parent offspring distribution")
    label(ax, .015, .86, r"$d=0.30,\quad q=0.04,\quad \mu=0,\quad\mathrm{no\ therapy}$", 14)
    node(ax, (.16, .54), "1 parent", .23, .15)
    for y, n, p in zip((.73, .50, .27), (0, 2, 3), offspring[[0,2,3]]):
        node(ax, (.80, y), f"{n} offspring", .29, .14, RED if n == 0 else BLUE)
        arrow(ax, (.285, .54), (.645, y), f"{p:.3f}", (.48, (y+.54)/2+.035))
    label(ax, .015, .12, r"$f(z)=0.300+0.672z^2+0.028z^3$", 16)
    label(ax, .015, -.005, r"$E[X]=f'(1)=1.428$", 15)

    ax = fig.add_axes([.49, .11, .47, .36])
    panel(ax, "C", "Exact finite transition matrix")
    label(ax, .015, .86, r"$P_{ij}=[z^j]f(z)^i;\qquad \boldsymbol{\pi}_{t+1}=\boldsymbol{\pi}_tP$", 15)
    labels = ["0", "1", "2", "3", r"$\dagger$"]
    matrix_display(ax, matrix, .12, .65, .83, .47, labels, labels)
    label(ax, .015, .08, r"$\dagger$: first entry into $N\geq4$; stop and remain in this bookkeeping state.", 11)
    label(ax, .015, -.005, "Rows: current state. Columns: next state. Each row sums to 1.", 11, color=GRAY)
    fig.text(.04, .035, "Figure 1 | Panel A describes the full model. Panels B-C remove mutations, genotype advantages and therapy; the threshold 4 is illustrative.\n"
             "Extinction (0) is biological absorption. Overflow (dagger) is imposed stopping, not a carrying capacity or a valid unrestricted aggregation.", fontsize=10, color=GRAY)
    save(fig, "06_markov_chain_formulation")


def state_figure(config, qmatrix, before, during, fitness):
    fig = plt.figure(figsize=(15, 10.8))
    fig.text(.04, .965, "Control states, mutation probabilities and offspring means", fontsize=21, weight="bold", color=INK)
    fig.text(.04, .93, "The finite execution controller is distinct from the population-state Markov kernel", fontsize=12, color=GRAY)
    ax = fig.add_axes([.04, .14, .43, .74])
    panel(ax, "A", "Execution state machine")
    node(ax, (.50, .87), "Initialize\nN = 100, A = 0", .45, .09)
    node(ax, (.50, .73), "Record current state\nCheck stopping conditions", .62, .10)
    arrow(ax, (.50, .822), (.50, .785))
    node(ax, (.50, .54), "Continue\nNo stop condition holds", .53, .09)
    arrow(ax, (.50, .675), (.50, .591), "otherwise", (.50, .637))
    node(ax, (.50, .32), "Update treatment\nalpha = A OR (N >= 20,000)", .68, .11)
    arrow(ax, (.50, .491), (.50, .380))
    node(ax, (.50, .14), "Draw offspring + mutations\nN <- N', A <- alpha; t <- t+1", .68, .105)
    arrow(ax, (.50, .26), (.50, .197))
    arrow(ax, (.151, .14), (.05, .14))
    ax.plot([.05, .05], [.14, .73], color=GRAY, lw=1.2)
    arrow(ax, (.05, .73), (.18, .73))
    node(ax, (.88, .62), "STOP", .18, .055, RED)
    arrow(ax, (.814, .73), (.88, .655), rad=-.12, color=RED)
    label(ax, .09, .445, "A = 0: untreated  |  A = 1: persistent therapy", 11,
          bbox=dict(facecolor="white", edgecolor="none", pad=2))
    label(ax, .02, .035, "Stop priority: extinction > target_size > max_cells > final generation.\nTreatment update shown for therapy=True; otherwise alpha = A.\nN denotes total burden here; transition draws use the full genotype vector.", 10, color=GRAY)

    ax = fig.add_axes([.53, .50, .43, .38])
    panel(ax, "B", "Two-type matrices: sensitive (S) and resistant (R)")
    label(ax, .015, .85, "Reduced example: resistance bit only; other growth effects and IL11 off.", 10, color=GRAY)
    label(ax, .015, .74, r"$Q=$", 17)
    matrix_display(ax, qmatrix, .14, .74, .40, .19, ["S", "R"], ["S", "R"])
    label(ax, .61, .70, r"$\mu=0.001$"+"\n"+r"$\sum_hQ_{gh}=1$", 14)
    label(ax, .015, .48, r"$M_{gh}=m_gQ_{gh},\quad m_g=(1-\delta_g)(2+q_g)$", 15)
    label(ax, .015, .35, r"$M_{\mathrm{off}}=$", 14)
    matrix_display(ax, before, .23, .34, .72, .13, fmt=".7g")
    label(ax, .015, .16, r"$M_{\mathrm{on}}=$", 14)
    matrix_display(ax, during, .23, .16, .72, .13, fmt=".7g")

    ax = fig.add_axes([.53, .13, .43, .32])
    panel(ax, "C", "Numerical fitness under default parameters")
    table = ax.table(cellText=fitness, colLabels=["Parent phenotype", "d", "q", "m: OFF", "m: ON"],
                     colWidths=[.43,.12,.12,.16,.17], bbox=[0,.20,1,.65], cellLoc="center")
    table.auto_set_font_size(False)
    table.set_fontsize(11)
    for (r,c), cell in table.get_celld().items():
        cell.set_edgecolor("#CED5DA")
        cell.set_linewidth(.5)
        cell.set_text_props(color=INK)
        if r == 0:
            cell.set_facecolor("#E9EFF3")
            cell.set_text_props(weight="bold")
    label(ax, .015, .13, r"$E[\mathbf{N}_{t+1}\mid\mathbf{N}_t,A_t]=\mathbf{N}_tM(\mathbf{N}_t,\alpha)$", 15)
    label(ax, .015, .01, "m > 1: expected expansion; m < 1: expected decline.\nM is NOT a transition-probability matrix; its row sums equal m.", 11)
    fig.text(.04, .035, "Figure 2 | Q and M use row-parent / column-offspring ordering (S, R). Panel C gives conditional offspring means before mutation.\n"
             "Natural death d excludes therapy; ON uses kS = 0.90, kR = 0.04. IL11 is off. Relapse is a trajectory classification, not a controller state.", fontsize=10, color=GRAY)
    save(fig, "07_state_machine_and_matrices")


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.family":"DejaVu Sans", "mathtext.fontset":"stix",
                         "svg.fonttype":"none", "axes.unicode_minus":False})
    c = Config()
    p, law = finite_count_matrix(c.death_rate, c.triple_probability)
    q = np.array([[1-c.mutation_rate, c.mutation_rate], [0, 1.]])
    m0 = (1-c.death_rate)*(2+c.triple_probability)
    off = np.diag([m0, m0]) @ q
    on = np.diag([m0*(1-c.sensitive_therapy_kill), m0*(1-c.resistant_therapy_kill)]) @ q
    # Meaningful analytic checks: the exact finite law and row/mean conventions.
    np.testing.assert_allclose(p.sum(axis=1), 1, atol=1e-14)
    np.testing.assert_allclose(p[2], [.09,0,.4032,.0168,.49])
    np.testing.assert_allclose(p[3], [.027,0,.18144,.00756,.784])
    np.testing.assert_allclose(q.sum(axis=1), 1)
    np.testing.assert_allclose(off.sum(axis=1), [1.428,1.428])
    np.testing.assert_allclose(on.sum(axis=1), [.1428,1.37088])
    fitness = []
    phenotypes = [("Baseline sensitive",0,False), ("Oncogene sensitive",1,False),
                  ("Oncogene + TSG loss",1|(1<<4)|(1<<20),False),
                  ("Baseline resistant",1<<3,True)]
    for name, genotype, resistant in phenotypes:
        d, triple = growth_probabilities(genotype,c)
        mean = (1-d)*(2+triple)
        kill = c.resistant_therapy_kill if resistant else c.sensitive_therapy_kill
        fitness.append([name,f"{d:.2f}",f"{triple:.2f}",f"{mean:.4f}",f"{mean*(1-kill):.5f}"])
    markov_figure(c,p,law)
    state_figure(c,q,off,on,fitness)
    for name, values, names in (("count_transition_P",p,["0","1","2","3","overflow"]),
                               ("mutation_Q",q,["S","R"]), ("offspring_M_off",off,["S","R"]),
                               ("offspring_M_on",on,["S","R"])):
        with (OUT/f"{name}.csv").open("w",newline="") as handle:
            writer=csv.writer(handle)
            writer.writerow(["current / next",*names])
            writer.writerows([[name,*row] for name,row in zip(names,values)])
    (OUT/"matrix_assumptions.json").write_text(json.dumps({
        "P":"One-type chain, no mutation, no therapy, d=0.30, q=0.04; first N>=4 is absorbing bookkeeping overflow.",
        "Q_and_M":"Only resistance is retained; no oncogene/TSG advantage or IL11; mu=0.001.",
        "orientation":"Rows=current/parent, columns=next/offspring; probability distributions are row vectors.",
        "verification":"Exact convolution probabilities, P/Q row sums, and M expected row sums checked."},indent=2)+"\n")
    print(f"Verified matrices and saved PNG/SVG figures to {OUT}")


if __name__ == "__main__":
    main()
