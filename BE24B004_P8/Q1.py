import matplotlib.pyplot as plt
import numpy as np

values = {
    'A': 13.85, 'D': 11.61, 'C': 15.37, 'E': 11.38, 'F': 13.93,
    'G': 13.34, 'H': 13.82, 'I': 15.28, 'K': 11.58, 'L': 14.13,
    'M': 13.86, 'N': 13.02, 'P': 12.35, 'Q': 12.61, 'R': 13.10,
    'S': 13.39, 'T': 12.70, 'V': 14.56, 'W': 15.48, 'Y': 13.88
}
def read_fasta(filepath):
    with open(filepath, "r") as f:
        sequences = []
        seq = ""
        for line in f:
            if line.startswith(">"):
                if seq:
                    sequences.append(seq)
                    seq = ""
            else:
                seq += line.strip()
        if seq:
            sequences.append(seq)
    return sequences

sequences = read_fasta("Smaller_Q1.fasta" if False else "Q1.fasta")

fig, ax = plt.subplots(nrows=len(sequences), dpi=150, figsize=(9, 2.5 * len(sequences)))

for count, sequence in enumerate(sequences):

    hydrophobicity_profile = [values[residue] for residue in sequence]
    average = np.mean(hydrophobicity_profile)
    deviation = [0 if hydrophobicity < average else 1 for hydrophobicity in hydrophobicity_profile]
    
    helices = set()
    sheets = set()
    
    for i in range(len(sequence) - 4 + 1):
        if deviation[i] == deviation[i + 1] and deviation[i] != deviation[i + 2] and deviation[i + 2] == deviation[i + 3]:
            helices.update({i, i + 1, i + 2, i + 3})
        elif deviation[i] != deviation[i + 1] and deviation[i] == deviation[i + 2] and deviation[i] != deviation[i + 3]:
            sheets.update({i, i + 1, i + 2, i + 3})

    ax[count].plot(range(len(sequence)), hydrophobicity_profile, c="black", marker="o", markersize=3)
    ax[count].set_title(f"\nHydrophobicity profile of sequence {count + 1}", fontsize=10)
    ax[count].set_xlabel("Residue number", fontsize=9)
    ax[count].set_ylabel("Hydrophobicity value", fontsize=9)
    ax[count].axhline(average, c="black", ls="--", lw=0.75)

    for position in helices:
        ax[count].axvspan(position - 0.5, position + 0.5, alpha=0.5, color='orange')
    for position in sheets:
        ax[count].axvspan(position - 0.5, position + 0.5, alpha=0.5, color='green')

fig.tight_layout()
plt.show()