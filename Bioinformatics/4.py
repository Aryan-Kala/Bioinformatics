def needleman_wunsch(seq1, seq2, match=2, mismatch=-1, gap=-2):
    n = len(seq1)
    m = len(seq2)
    
    # Initialize scoring matrix with zeros
    F = [[0 for j in range(n + 1)] for i in range(m + 1)]
    
    # Initialize gap penalties for first row and column
    for i in range(m + 1):
        F[i][0] = i * gap
    for j in range(n + 1):
        F[0][j] = j * gap
        
    # Fill the scoring matrix
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            match_score = F[i-1][j-1] + (match if seq2[i-1] == seq1[j-1] else mismatch)
            delete_score = F[i-1][j] + gap
            insert_score = F[i][j-1] + gap
            F[i][j] = max(match_score, delete_score, insert_score)
            
    # Backtrack to find the optimal alignment
    alignment1 = ""
    alignment2 = ""
    i = m
    j = n
    
    while i > 0 and j > 0:
        score_current = F[i][j]
        score_diag = F[i-1][j-1]
        score_up = F[i-1][j]
        score_left = F[i][j-1]
        
        match_val = match if seq2[i-1] == seq1[j-1] else mismatch
        
        if score_current == score_diag + match_val:
            alignment1 += seq1[j-1]
            alignment2 += seq2[i-1]
            i -= 1
            j -= 1
        elif score_current == score_up + gap:
            alignment1 += '-'
            alignment2 += seq2[i-1]
            i -= 1
        else:
            alignment1 += seq1[j-1]
            alignment2 += '-'
            j -= 1
            
    while i > 0:
        alignment1 += '-'
        alignment2 += seq2[i-1]
        i -= 1
    while j > 0:
        alignment1 += seq1[j-1]
        alignment2 += '-'
        j -= 1
        
    return alignment1[::-1], alignment2[::-1]

seq1 = "ACAGTCGAACG"
seq2 = "ACCGTCCG"
aln1, aln2 = needleman_wunsch(seq1, seq2)

print("Sequence 1:", aln1)
print("Sequence 2:", aln2)