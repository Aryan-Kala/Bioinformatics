import numpy as np

def smith_waterman(seq1, seq2, match=2, mismatch=-1, gap=-2):
    n, m = len(seq1), len(seq2)
    
    # Scoring matrix (local alignment: borders start at 0, no negative scores)
    M = np.zeros((n + 1, m + 1), dtype=int)
    
    max_score = 0
    max_pos = (0, 0)
    
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            s = match if seq1[i-1] == seq2[j-1] else mismatch
            M[i, j] = max(0,
                          M[i-1, j-1] + s,      # diagonal (match/mismatch)
                          M[i-1, j] + gap,      # gap in seq2
                          M[i, j-1] + gap)      # gap in seq1
            if M[i, j] > max_score:
                max_score = M[i, j]
                max_pos = (i, j)
                
    # Traceback from the highest-scoring cell, stop when a cell hits 0
    a, b = "", ""
    i, j = max_pos
    
    while M[i, j] != 0:
        s = match if seq1[i-1] == seq2[j-1] else mismatch
        if M[i, j] == M[i-1, j-1] + s:
            a, b = seq1[i-1] + a, seq2[j-1] + b
            i, j = i - 1, j - 1
        elif M[i, j] == M[i-1, j] + gap:
            a, b = seq1[i-1] + a, "-" + b
            i -= 1
        else:
            a, b = "-" + a, seq2[j-1] + b
            j -= 1
            
    return M, a, b, max_score

seq1 = "ACGTATCGCGTATA"
seq2 = "GATGCGTATCG"

M, a1, a2, score = smith_waterman(seq1, seq2)

print("Scoring Matrix:\n", M)
print("\nBest local alignment:")
print(a1)
print(a2)
print("Score:", score)