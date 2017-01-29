"""Create a plot of ORFs in a DNA sequence.

Usage:
        python gene_plot.py
"""

import seaborn as sns
from matplotlib import pyplot as plt
from matplotlib.collections import BrokenBarHCollection

from gene_finder import find_all_ORFs_both_strands, longest_ORF_noncoding
from load import load_seq


# from matplotlib.patches import Rectangle

def create_gene_figure(output_path, dna, orfs):
    """Create a matplotlib figure that shows ORFs against a DNA sequence."""
    plt.figure(figsize=(10, 1))
    plt.yticks([], [])
    ax = plt.gca()
    ax.set_xlim(0, len(dna))
    ax.set_ylim(0, 2)
    ax.patch.set_facecolor('white')
    ax.grid(False)
    coll = BrokenBarHCollection([(0, len(dna))], (0, 1), facecolors=sns.color_palette('pastel'))
    ax.add_collection(coll)
    xspans = [(dna.index(orf), len(orf)) for orf in orfs]
    coll = BrokenBarHCollection(xspans, (1, 1), facecolors=sns.color_palette('pastel')[1:])
    ax.add_collection(coll)
    plt.savefig(output_path, bbox_inches='tight')


if __name__ == "__main__":
    image_path = 'genes.png'
    data_file = 'data/X73525.fa'

    print('loading sequence from', data_file)
    dna = load_seq(data_file)
    print("finding longest noncoding ORF...")
    threshold = longest_ORF_noncoding(dna, 1500)
    print("longest noncoding ORF length:", threshold)
    print("searching for ORFs...")
    orfs = [orf for orf in find_all_ORFs_both_strands(dna) if len(orf) > threshold]
    create_gene_figure(image_path, dna, orfs)
    print('wrote', image_path)
