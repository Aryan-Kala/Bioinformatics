# Provided molecular weights
amino_acid_weights = {
    "A": 85, "C": 115, "D": 130, "E": 145, "F": 160, 
    "G": 70, "H": 150, "I": 125, "K": 145, "L": 125, 
    "M": 143, "N": 130, "P": 110, "Q": 140, "R": 170, 
    "S": 100, "T": 115, "V": 110, "W": 200, "Y": 175
}

# The sequences
seq1 = 'RATPTRWPVGCFNRPWTKWSYDEALDGIKAAGYAWTGLLTASKPSLHHATATPEYLAALKQKSRHAA'
seq2 = 'AAAVMMGLAAIGAAIGIGILGGKFLEGAARQPDLIPLLRTQFFIVMGLVDAIPMIAVGLGLYVMFAVA'
seq3 = 'AADVSAAVGATGQSGMTYRLGLSWDWDKSWWQTSTGRLTGYWDAGYTYWEGGDEGAGKHSLSFAPVFVYEFAGDSIKPFIEAGIGVAAFSGTRVGDQNLGSSLNFEDRIGAGLKFANGQSVGVRAIHYSNAGLKQPNDGIESYSLFYKIPI'

def compute_mw(seq):
    """Calculates both the direct sum and the water-loss adjusted molecular weight."""
    raw_mw = 0
    for aa in seq.upper():
        if aa in amino_acid_weights:
            raw_mw += amino_acid_weights[aa]
            
    adjusted_mw = raw_mw - (18 * (len(seq) - 1))
    
    return raw_mw, adjusted_mw

seqs = [seq1, seq2, seq3]

for i, seq in enumerate(seqs, start=1):
    raw_weight, adjusted_weight = compute_mw(seq)
    print(f"Sequence {i}:")
    print(f"  - Direct Sum Weight: {raw_weight}")
    print(f"  - Adjusted Weight (water loss): {adjusted_weight}\n")