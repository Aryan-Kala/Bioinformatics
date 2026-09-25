
seq1 = 'RATPTRWPVGCFNRPWTKWSYDEALDGIKAAGYAWTGLLTASKPSLHHATATPEYLAALKQKSRHAA'
seq2 = 'AAAVMMGLAAIGAAIGIGILGGKFLEGAARQPDLIPLLRTQFFIVMGLVDAIPMIAVGLGLYVMFAVA'
seq3 = 'AADVSAAVGATGQSGMTYRLGLSWDWDKSWWQTSTGRLTGYWDAGYTYWEGGDEGAGKHSLSFAPVFVYEFAGDSIKPFIEAGIGVAAFSGTRVGDQNLGSSLNFEDRIGAGLKFANGQSVGVRAIHYSNAGLKQPNDGIESYSLFYKIPI'


Hgm_data = {"A": 0.25, "C": 0.38, "D": -0.90, "E": -0.74, "F": 1.19, "G": 0.00, "H": -0.40, "I": 1.38, "K": -1.50, "L": 1.06, "M": 0.64, "N": -0.78, "P": -0.55, "Q": -0.85, "R": -1.76, "S": -0.18, "T": -0.05, "V": 1.08, "W": 0.81, "Y": 0.26}
Ca_data = {"A": 32.0, "C": 34.0, "D": 42.0, "E": 48.0, "F": 62.0, "G": 22.0, "H": 54.0, "I": 56.0, "K": 60.0, "L": 58.0, "M": 60.0, "N": 44.0, "P": 38.0, "Q": 52.0, "R": 68.0, "S": 36.0, "T": 42.0, "V": 48.0, "W": 72.0, "Y": 66.0}
Et_data = {"A": -1.8, "C": -2.1, "D": -2.3, "E": -2.6, "F": -3.8, "G": -1.2, "H": -2.9, "I": -3.5, "K": -3.1, "L": -3.6, "M": -3.4, "N": -2.4, "P": -2.0, "Q": -2.8, "R": -3.5, "S": -1.9, "T": -2.2, "V": -3.0, "W": -4.2, "Y": -3.9}

def calc_properties(seq, name):
    total_Hgm, total_Ca, total_Et = 0, 0, 0
    N = len(seq)
    
    for aa in seq.upper():
        if aa in Hgm_data:
            total_Hgm += Hgm_data[aa]
            total_Ca += Ca_data[aa]
            total_Et += Et_data[aa]
            
    avg_Hgm = round(total_Hgm / N, 3)
    avg_Ca = round(total_Ca / N, 3)
    avg_Et = round(total_Et / N, 3)
    
    print(f"-- {name} --")
    print(f"Average Hgm (Hydrophobicity)  : {avg_Hgm}")
    print(f"Average Ca (Helical Contact)  : {avg_Ca}")
    print(f"Average Et (Total Non-bonded) : {avg_Et}")
    
    return avg_Hgm, avg_Ca, avg_Et


res1 = calc_properties(seq1, "Sequence 1")
res2 = calc_properties(seq2, "Sequence 2")
res3 = calc_properties(seq3, "Sequence 3")
