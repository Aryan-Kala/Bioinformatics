seq1 = "AATCTATA"
seq2 = "AAG--ATA"

match_score = 1
mismatch_score = 0
gap_origination = -2
gap_length = -1

score = 0
in_gap = False

for a, b in zip(seq1, seq2):
    if a == b:
        score += match_score
        in_gap = False
    elif a == '-' or b == '-':
        if not in_gap:
            score += gap_origination
            in_gap = True
        score += gap_length
    else:
        score += mismatch_score
        in_gap = False

print(f"Alignment Score: {score}")