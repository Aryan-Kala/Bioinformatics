import pandas as pd

seq1 = 'RATPTRWPVGCFNRPWTKWSYDEALDGIKAAGYAWTGLLTASKPSLHHATATPEYLAALKQKSRHAA'
seq2 = 'AAAVMMGLAAIGAAIGIGILGGKFLEGAARQPDLIPLLRTQFFIVMGLVDAIPMIAVGLGLYVMFAVA'
seq3 = 'AADVSAAVGATGQSGMTYRLGLSWDWDKSWWQTSTGRLTGYWDAGYTYWEGGDEGAGKHSLSFAPVFVYEFAGDSIKPFIEAGIGVAAFSGTRVGDQNLGSSLNFEDRIGAGLKFANGQSVGVRAIHYSNAGLKQPNDGIESYSLFYKIPI'

def compute_composition(seq):

    length = len(seq)
    composition = {}
    
    for aa in seq.upper():
        composition[aa] = composition.get(aa, 0) + 1
        
    for aa in composition:
        composition[aa] = round((composition[aa] / length) * 100, 2)
        
    return composition


comp1 = compute_composition(seq1)
comp2 = compute_composition(seq2)
comp3 = compute_composition(seq3)


amino_acids = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", 
               "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]

table_data = {
    "Amino Acid": amino_acids,
    "Sequence 1 (%)": [comp1.get(aa, 0.0) for aa in amino_acids],
    "Sequence 2 (%)": [comp2.get(aa, 0.0) for aa in amino_acids],
    "Sequence 3 (%)": [comp3.get(aa, 0.0) for aa in amino_acids]
}

df = pd.DataFrame(table_data)


print("Amino Acid Composition Table:\n")
print(df.to_string(index=False))
