def smith_waterman(seq1, seq2, match=2, mismatch=-1, gap=-2):
    n = len(seq1)
    m = len(seq2)
    
    F = [[0 for j in range(n + 1)] for i in range(m + 1)]
    
    max_score = 0
    max_i = 0
    max_j = 0
    
    # Fill the scoring matrix with 0 to prevent negative scores
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            match_score = F[i-1][j-1] + (match if seq2[i-1] == seq1[j-1] else mismatch)
            delete_score = F[i-1][j] + gap
            insert_score = F[i][j-1] + gap
            
            F[i][j] = max(0, match_score, delete_score, insert_score)
            
            if F[i][j] > max_score:
                max_score = F[i][j]
                max_i = i
                max_j = j
                
    # Backtrack from the maximum score cell until reaching 0
    alignment1 = ""
    alignment2 = ""
    i = max_i
    j = max_j
    
    while i > 0 and j > 0 and F[i][j] > 0:
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
        elif score_current == score_left + gap:
            alignment1 += seq1[j-1]
            alignment2 += '-'
            j -= 1
        else:
            break
            
    return alignment1[::-1], alignment2[::-1], max_score

seq1 = "ACGTATCGCGTATA"
seq2 = "GATGCGTATCG"
aln1, aln2, score = smith_waterman(seq1, seq2)

print("Max Local Score:", score)
print("Sequence 1:", aln1)
print("Sequence 2:", aln2)