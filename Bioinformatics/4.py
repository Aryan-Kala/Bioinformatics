import numpy as np

def needleman_wunsch(seq1, seq2, match=2, mismatch=-1, gap=-2):
    n, m = len(seq1), len(seq2)
    
    # Scoring matrix
    M = np.zeros((n + 1, m + 1), dtype=int)
    M[:, 0] = np.arange(n + 1) * gap
    M[0, :] = np.arange(m + 1) * gap
    
    # Fill matrix
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            s = match if seq1[i-1] == seq2[j-1] else mismatch
            M[i, j] = max(M[i-1, j-1] + s,
                          M[i-1, j] + gap,
                          M[i, j-1] + gap)
                          
    # Traceback
    a, b = "", ""
    i, j = n, m
    
    while i or j:
        if i and j:
            s = match if seq1[i-1] == seq2[j-1] else mismatch
            if M[i, j] == M[i-1, j-1] + s:
                a, b = seq1[i-1] + a, seq2[j-1] + b
                i, j = i - 1, j - 1
                continue
                
        if i and M[i, j] == M[i-1, j] + gap:
            a, b = seq1[i-1] + a, "-" + b
            i -= 1
        else:
            a, b = "-" + a, seq2[j-1] + b
            j -= 1
            
    return M, a, b

seq1 = "ACCGTCCG"
seq2 = "ACAGTCGAACG"

M, a1, a2 = needleman_wunsch(seq1, seq2)

print("Scoring Matrix:\n", M)
print("\nAlignment:")
print(a2)
print(a1)
print("Score:", M[-1, -1])