import pandas as pd

seq1 = 'RATPTRWPVGCFNRPWTKWSYDEALDGIKAAGYAWTGLLTASKPSLHHATATPEYLAALKQKSRHAA'
seq2 = 'AAAVMMGLAAIGAAIGIGILGGKFLEGAARQPDLIPLLRTQFFIVMGLVDAIPMIAVGLGLYVMFAVA'
seq3 = 'AADVSAAVGATGQSGMTYRLGLSWDWDKSWWQTSTGRLTGYWDAGYTYWEGGDEGAGKHSLSFAPVFVYEFAGDSIKPFIEAGIGVAAFSGTRVGDQNLGSSLNFEDRIGAGLKFANGQSVGVRAIHYSNAGLKQPNDGIESYSLFYKIPI'

amino_acids = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", 
               "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]

def compute_pair_preferences(seq, seq_name):
    N = len(seq)
    
    ni_counts = {aa: seq.count(aa) for aa in amino_acids}
    
    nij_counts = {(a1, a2): 0 for a1 in amino_acids for a2 in amino_acids}
    for k in range(N - 1):
        pair = (seq[k], seq[k+1])
        if pair in nij_counts:
            nij_counts[pair] += 1
            
    pref_a = {aa: {} for aa in amino_acids}
    pref_b = {aa: {} for aa in amino_acids}
    pref_c = {aa: {} for aa in amino_acids}
    
    list_a, list_b, list_c = [], [], []
    
    for i in amino_acids:
        for j in amino_acids:
            nij = nij_counts[(i, j)]
            ni = ni_counts[i]
            nj = ni_counts[j]
            pair_name = f"{i}-{j}"
            
            val_a = round((nij * 100) / (ni + nj), 2) if (ni + nj) > 0 else 0.0
            pref_a[i][j] = val_a
            list_a.append((pair_name, val_a))
            
            val_b = round((nij * 100) / (N - 1), 2) if (N - 1) > 0 else 0.0
            pref_b[i][j] = val_b
            list_b.append((pair_name, val_b))
            
            val_c = round((nij * 100) / (ni * nj), 2) if (ni * nj) > 0 else 0.0
            pref_c[i][j] = val_c
            list_c.append((pair_name, val_c))

    top10_a = sorted(list_a, key=lambda x: x[1], reverse=True)[:10]
    top10_b = sorted(list_b, key=lambda x: x[1], reverse=True)[:10]
    top10_c = sorted(list_c, key=lambda x: x[1], reverse=True)[:10]
    
    df_a = pd.DataFrame(pref_a).T
    df_b = pd.DataFrame(pref_b).T
    df_c = pd.DataFrame(pref_c).T
    
    print(f"-- {seq_name} --\n")
    print("Preference (a) Table ")
    print(df_a.to_string())
    print("\nTop 10 Pairs (a):", top10_a, "\n")
    
    print(" Preference (b) Table ")
    print(df_b.to_string())
    print("\nTop 10 Pairs (b):", top10_b, "\n")
    
    print(" Preference (c) Table ")
    print(df_c.to_string())
    print("\nTop 10 Pairs (c):", top10_c, "\n")
    print("----\n")

seqs = [(seq1, "Sequence 1"), (seq2, "Sequence 2"), (seq3, "Sequence 3")]

for seq, name in seqs:
    compute_pair_preferences(seq, name)