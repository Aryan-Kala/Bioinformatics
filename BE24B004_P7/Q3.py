seq1 = 'RATPTRWPVGCFNRPWTKWSYDEALDGIKAAGYAWTGLLTASKPSLHHATATPEYLAALKQKSRHAA'
seq2 = 'AAAVMMGLAAIGAAIGIGILGGKFLEGAARQPDLIPLLRTQFFIVMGLVDAIPMIAVGLGLYVMFAVA'
seq3 = 'AADVSAAVGATGQSGMTYRLGLSWDWDKSWWQTSTGRLTGYWDAGYTYWEGGDEGAGKHSLSFAPVFVYEFAGDSIKPFIEAGIGVAAFSGTRVGDQNLGSSLNFEDRIGAGLKFANGQSVGVRAIHYSNAGLKQPNDGIESYSLFYKIPI'


group_data = {
    "A": (8.47, 8.95), "D": (5.97, 5.91), "C": (1.39, 0.47), "E": (6.32, 4.78),
    "T": (5.79, 6.54), "F": (3.91, 3.68), "G": (7.82, 8.54), "H": (2.26, 1.25),
    "I": (5.71, 4.77), "V": (7.02, 6.76), "K": (5.76, 4.93), "L": (8.48, 8.78),
    "M": (2.21, 1.56), "N": (4.54, 5.74), "W": (1.44, 1.24), "P": (4.63, 3.74),
    "Q": (3.82, 4.75), "R": (4.93, 5.24), "S": (5.94, 8.05), "Y": (3.58, 4.13)
}

def get_composition(seq):
    """Calculates the percentage composition of each amino acid."""
    length = len(seq)
    comp = {}
    for aa in seq.upper():
        comp[aa] = comp.get(aa, 0) + 1
    for aa in comp:
        comp[aa] = round((comp[aa] / length) * 100, 2)
    return comp

def identify_group(seq, seq_name):
    """Calculates the absolute difference and identifies the closest group."""
    seq_comp = get_composition(seq)
    
    diff_A = 0
    diff_B = 0
    
    
    for aa, (val_A, val_B) in group_data.items():    
        seq_val = seq_comp.get(aa, 0)
        
        diff_A += abs(seq_val - val_A)
        diff_B += abs(seq_val - val_B)
        
    group = "Group A" if diff_A < diff_B else "Group B"
    
    print(f"{seq_name}: Belongs to {group}")
    print(f"  -> Difference from Group A: {diff_A:.2f}")
    print(f"  -> Difference from Group B: {diff_B:.2f}\n")

identify_group(seq1, "Sequence 1")
identify_group(seq2, "Sequence 2")
identify_group(seq3, "Sequence 3")