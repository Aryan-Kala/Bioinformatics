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

sequences = read_fasta("Q1.fasta")

count = 1

for sequence in sequences:
    hydrophobicity_profile = [values[residue] for residue in sequence]
    average = np.mean(hydrophobicity_profile)
    deviation = [0 if hydrophobicity < average else 1 for hydrophobicity in hydrophobicity_profile]
    
    alphas = [i for i in range(len(sequence) - 8 + 1) if ((deviation[i:i + 8] == [0, 0, 1, 1, 0, 0, 1, 1]) or (deviation[i:i + 8] == [1, 1, 0, 0, 1, 1, 0, 0]))]
    betas = [i for i in range(len(sequence) - 6 + 1) if ((deviation[i:i + 6] == [0, 1, 0, 1, 0, 1]) or (deviation[i:i + 6] == [1, 0, 1, 0, 1, 0]))]
    
    alpha_amphipathicities = []
    
    for alpha_start in alphas:
        a1 = 0
        a2 = 0
        a3 = 0
        a4 = 0
        for j in range(0, 8, 4):
            a1 += values[sequence[alpha_start + j]]
            a2 += values[sequence[alpha_start + 1 + j]]
            a3 += values[sequence[alpha_start + 2 + j]]
            a4 += values[sequence[alpha_start + 3 + j]]
        a1 = a1 / 2
        a2 = a2 / 2
        a3 = a3 / 2
        a4 = a4 / 2
        if deviation[alpha_start] == deviation[alpha_start + 1]:
            alpha_amphipathicities.append(abs((a1 + a2) - (a3 + a4)))
        else:
            alpha_amphipathicities.append(abs((a1 + a4) - (a2 + a3)))
            
    beta_amphipathicities = []
    
    for beta_start in betas:
        b1 = 0
        b2 = 0
        for j in range(0, 6, 2):
            b1 += values[sequence[beta_start + j]]
            b2 += values[sequence[beta_start + 1 + j]]
        b1 = b1 / 3
        b2 = b2 / 3
        if deviation[beta_start] == 0:
            beta_amphipathicities.append(b2 - b1)
        else:
            beta_amphipathicities.append(b1 - b2)
            
    print(f"Sequence {count}")
    if len(alphas) != 0:
        for i in range(len(alphas)):
            print(f"Helix of length 8 found at position {alphas[i]} with amphipathicity {round(alpha_amphipathicities[i], 3)}")
    else:
        print("No helices of length 8 were found")
        
    if len(betas) != 0:
        for i in range(len(betas)):
            print(f"Sheet of length 6 found at position {betas[i]} with amphipathicity {round(beta_amphipathicities[i], 3)}")
    else:
        print("No sheets of length 6 were found")
    print("")
    
    count += 1