import pandas as pd

residues = ["G", "A", "V", "L", "I", "P", "F", "Y", "W", "S", "T", "C", "M", "N", "Q", "D", "E", "K", "R", "H"]

def composition(sequence):
    comp = {residue: 100 * sequence.count(residue) / len(sequence) for residue in residues}
    return comp

def hamming(sequence1, sequence2):
    comp1 = composition(sequence1)
    comp2 = composition(sequence2)
    
    ham_dist = 0
    
    for residue in residues:
        ham_dist += abs(comp1[residue] - comp2[residue])
        
    return ham_dist

def euclidean(sequence1, sequence2):
    comp1 = composition(sequence1)
    comp2 = composition(sequence2)
    
    euc_dist = 0
    
    for residue in residues:
        euc_dist += (comp1[residue] - comp2[residue]) ** 2
        
    euc_dist = euc_dist ** 0.5
    
    return euc_dist

seq1 = "AMENLNMDLLYMAAAVMMGLAAIGAAIGIGILGGKFLEGAARQPDLIPLLRTQFFIVMGLVDAIPMIAVGLGLYVMFAVA"
seq2 = "AADVSAAVGATGQSGMTYRLGLSWDWDKSWWQTSTGRLTGYWDAGYTYWEGGDEGAGKHSLSFAPVFVYEFAGDSIKPFIEAGIGVAAFSGTRVGDQNLGSSLNFEDRIGAGLKFANGQSVGVRAIHYSNAGLKQPNDGIESYSLFYKIPI"
seq3 = "MALLPAAPGAPARATPTRWPVGCFNRPWTKWSYDEALDGIKAAGYAWTGLLTASKPSLHHATATPEYLAALKQKSRHAA"

distances = pd.DataFrame({"Hamming Distance": [hamming(seq1, seq2), hamming(seq1, seq3), hamming(seq2, seq3)],
                          "Euclidean Distance": [euclidean(seq1, seq2), euclidean(seq1, seq3), euclidean(seq2, seq3)]})
distances.index = ["1 and 2", "1 and 3", "2 and 3"]

print(distances)
print("")
print(f"Based on Hamming distance, sequences {distances[distances.columns[0]].idxmin()} are the closest.")
print(f"Based on Euclidean distance, sequences {distances[distances.columns[1]].idxmin()} are the closest.")