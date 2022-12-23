import matplotlib
import matplotlib.pyplot as plt

from seqSNP import seqSnp21
from getSNP import getSnpInfo
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas

plt.ioff
matplotlib.rcParams["toolbar"] = "None"


def set_colors(seq):
    def decide_color(base):
        base = base.upper()
        if base == "A":
            return "coral"  # "tab:red"
        elif base == "T":
            return "skyblue"  # "tab:blue"
        elif base == "G":
            return "yellowgreen"  # "tab:green"
        elif base == "C":
            return "orange"  # "tab:orange"
        else:
            return "tab:gray"

    resSeq = []
    for base in seq:
        resSeq = resSeq + [decide_color(base)]

    return tuple(resSeq)


# input : seqs(list[2]), gene(str)
def display_gene(seqs, gene):
    seqlen = len(seqs[0])
    boxPos = []

    for i in range(0, seqlen):
        boxPos = boxPos + [(i * 10, 10)]

    fig, ax = plt.subplots()
    ax.broken_barh(
        boxPos, (2 + 4, 2), facecolors=set_colors(seqs[1]), linewidth=0.2, edgecolor="w"
    )
    ax.broken_barh(
        boxPos, (4 + 4, 2), facecolors=set_colors(seqs[0]), linewidth=0.2, edgecolor="w"
    )

    i = 0
    for base in seqs[1]:
        ax.text(i * 10 + 0.2, 2 + 0.5 + 4, base, {"color": "k", "fontsize": 14})
        i = i + 1

    i = 0
    for base in seqs[0]:
        ax.text(i * 10 + 0.2, 4 + 0.5 + 4, base, {"color": "k", "fontsize": 14})
        i = i + 1

    ax.text(
        seqlen * 5 - 2,
        11,
        "Sequence of " + gene,
        {"color": "black", "fontsize": 14},
        weight="bold",
        family="monospace",
        verticalalignment="bottom",
        horizontalalignment="center",
    )

    ax.set_ylim(6, 12)
    ax.set_xlim(0, seqlen * 10)

    ax.set_yticks([3 + 4, 5 + 4], labels=["MN", "ACS"])
    ax.grid(False)
    ax.set_xticks([])

    fig.set_size_inches(seqlen / 5, 1, forward=True)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.spines["left"].set_visible(False)

    return FigureCanvas(fig)


# seq1 = "ATGCTTGAGCTTAACTGCTGAACTTGAG"
# seq2 = "AGGCTAAGTATAGCTTTTGTGCGCTTGC"
# seq1 = "ATGCTTGAGCTTAACTGC"
# seq2 = "AGGCTAAGTATAGCTTTT"
# seqs = [seq1, seq2]
# display_gene(seqs, "aaa1")


def display_snp(snp_id):
    info = getSnpInfo(int(snp_id))
    seqs = seqSnp21(info[0], info[1], info[2], info[3])
    return display_gene(seqs, info[4])
