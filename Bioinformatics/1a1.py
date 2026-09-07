# Two protein sequences
seq1 = "MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH"
seq2 = "MVHWTAEEKQLITGLWGKVNVAECGAEALARLLIVYPWTQRFFASFGNLSPTAILGPNMVRAHGKKVLTSFGDAVKNLDNIKNTFSQLSELHCDKLHVDPENFRLLGDILIIVLAAHFSKDFTPECQAAWQKLRVVVAHALARKYH"

similar_seq = []
p1 = 0
sim = ""

# Limit to the shorter sequence length to prevent out-of-range errors
max_len = min(len(seq1), len(seq2))

while p1 < max_len:
    if seq1[p1] == seq2[p1]:
        sim += seq1[p1]
        p1 += 1
    else:
        if sim != "":
            similar_seq.append(sim)
            sim = ""
        p1 += 1

if sim != "":
    similar_seq.append(sim)

print(similar_seq)