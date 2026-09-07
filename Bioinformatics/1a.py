import matplotlib.pyplot as plt

# Two protein sequences
seq1 = "MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH"
seq2 = "MVHWTAEEKQLITGLWGKVNVAECGAEALARLLIVYPWTQRFFASFGNLSPTAILGPNMVRAHGKKVLTSFGDAVKNLDNIKNTFSQLSELHCDKLHVDPENFRLLGDILIIVLAAHFSKDFTPECQAAWQKLRVVVAHALARKYH"

plt.figure(figsize=(10, 10))

for i in range(len(seq1)):
    for j in range(len(seq2)):
        if seq1[i] == seq2[j]:
            plt.plot(i + 1, j + 1, 'k.')

plt.xlabel("Human Heamoglobin")
plt.ylabel("Chicken Heamoglobin")
plt.title("Dot Plot of Two Protein Sequences")

plt.xlim(0, len(seq1) + 1)
plt.ylim(0, len(seq2) + 1)

plt.grid(True)

plt.show()