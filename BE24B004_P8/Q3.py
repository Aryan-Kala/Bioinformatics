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

sequences = read_fasta("Q2.fasta")

window_size1 = 9
half_width1 = window_size1 // 2
window_size2 = 19
half_width2 = window_size2 // 2

for sequence in sequences:
    
    hydrophobicity_profile1 = [np.mean([values[residue] for residue in sequence[i - half_width1:i + half_width1 + 1]]) 
                              for i in range(half_width1, len(sequence) - half_width1)]
    
    average1 = np.mean(hydrophobicity_profile1)
    deviation1 = [0 if hydro < average1 else 1 for hydro in hydrophobicity_profile1]
    
    transmembrane1 = set()
    
    for i in range(len(hydrophobicity_profile1) - 4 + 1):
        if deviation1[i] == 1 and deviation1[i] == deviation1[i + 1] and deviation1[i + 1] == deviation1[i + 2] and deviation1[i + 2] == deviation1[i + 3]:
            transmembrane1.update({i, i + 1, i + 2, i + 3})
            
    plt.figure(dpi=200, figsize=(10, 5))
    plt.plot(range(len(hydrophobicity_profile1)), hydrophobicity_profile1, c="black", marker="o", markersize=3)
    plt.title(f"Hydrophobicity profile of sequence with window length {window_size1}")
    plt.xlabel("Residue number")
    plt.ylabel("Hydrophobicity value")
    plt.axhline(average1, c="black", ls="--", lw=0.75)
    
    for position in transmembrane1:
        plt.axvspan(position - 0.5, position + 0.5, alpha=0.5, color='blue', lw=0)

    hydrophobicity_profile2 = [np.mean([values[residue] for residue in sequence[i - half_width2:i + half_width2 + 1]]) 
                              for i in range(half_width2, len(sequence) - half_width2)]
    
    average2 = np.mean(hydrophobicity_profile2)
    deviation2 = [0 if hydro < average2 else 1 for hydro in hydrophobicity_profile2]
    
    transmembrane2 = set()
    
    for i in range(len(hydrophobicity_profile2) - 4 + 1):
        if deviation2[i] == 1 and deviation2[i] == deviation2[i + 1] and deviation2[i + 1] == deviation2[i + 2] and deviation2[i + 2] == deviation2[i + 3]:
            transmembrane2.update({i, i + 1, i + 2, i + 3})
    

    plt.figure(dpi=200, figsize=(10, 5))
    plt.plot(range(len(hydrophobicity_profile2)), hydrophobicity_profile2, c="black", marker="o", markersize=3)
    plt.title(f"Hydrophobicity profile of sequence with window length {window_size2}")
    plt.xlabel("Residue number")
    plt.ylabel("Hydrophobicity value")
    plt.axhline(average2, c="black", ls="--", lw=0.75)
    
    for position in transmembrane2:
        plt.axvspan(position - 0.5, position + 0.5, alpha=0.5, color='blue', lw=0)
    
    plt.show()