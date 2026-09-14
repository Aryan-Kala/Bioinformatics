Human_Beta = 'MVHLTPEEKTAVNALWGKVNVDAVGGEALGRLLVVYPWTQRFFESFGDLSSPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFSQLSELHCDKLHVDPENFRLLGNVLVCVLARNFGKEFTPQMQAAYQKVVAGVANALAHKYH'
Chicken_Beta = 'MVHWTAEEKQLITGLWGKVNVAECGAEALARLLIVYPWTQRFFASFGNLSSPTAILGNPMVRAHGKKVLTSFGDAVKNLDNIKNTFSQLSELHCDKLHVDPENFRLLGDILIIVLAAHFSKDFTPECQAAWQKLVRVVAHALARKYH'

blosum62_positives = {('A', 'G'), ('A', 'S'), ('R', 'K'), ('R', 'Q'), ('R', 'H'), ('N', 'D'), ('N', 'H'), ('N', 'S'), ('N', 'E'), ('N', 'K'), ('D', 'E'), ('Q', 'E'), ('Q', 'K'), ('Q', 'R'), ('Q', 'M'), ('H', 'Y'), ('I', 'V'), ('I', 'L'), ('I', 'M'), ('L', 'V'), ('L', 'M'), ('L', 'F'), ('K', 'E'), ('M', 'V'), ('F', 'Y'), ('F', 'W'), ('S', 'T'), ('W', 'Y')}

alignment_len = len(Human_Beta)
identities = 0

def Protein_Alignment(query_seq, target_seq):
    query_ungapped_len = len(query_seq.replace('-', ''))
    similarities = 0
    gaps = 0
    query_aligned_residues = 0
    global identities
    
    for i, j in zip(query_seq.upper(), target_seq.upper()):
        if i == '-' or j == '-':
            gaps += 1
        elif i == j:
            identities += 1
            similarities += 1
            query_aligned_residues += 1
        elif (i, j) in blosum62_positives or (j, i) in blosum62_positives:
            similarities += 1
            
    ide = (identities / alignment_len) * 100
    pos = (similarities / alignment_len) * 100
    cov = (query_aligned_residues / query_ungapped_len) * 100
    gap = (gaps / alignment_len) * 100
    return (ide, pos, cov, gap)

ide, pos, cov, gap = Protein_Alignment(Human_Beta, Chicken_Beta)
print("Sequence Identity:", ide)
print("Sequence Similarity:", pos)
print("Query Coverage:", cov)
print("Gap Percentage:", gap)