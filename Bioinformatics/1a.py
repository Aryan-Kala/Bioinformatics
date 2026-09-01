import matplotlib.pyplot as plt

human_hbb = "MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH"
chicken_hbb = "MVHWTAEEKQLITGLWGKVNVAECGAEALARLLIVYPWTQRFFASFGNLSSPTAILGNPMVRAHGKKVLTSFGDAVKNLDNIKNTFAQLSELHCDKLHVDPENFRLLGDILIIVLAAHFSKDFTPECQAAWQKLVRVVAHALARKYH"

x_coords = []
y_coords = []

for i, res1 in enumerate(human_hbb):
    for j, res2 in enumerate(chicken_hbb):
        if res1 == res2:
            x_coords.append(i + 1)
            y_coords.append(j + 1)

# Generate Plot
plt.figure(figsize=(7, 7))
plt.scatter(x_coords, y_coords, color='black', s=6)
plt.xlabel("Human Hb β Chain (Residues 1-147)")
plt.ylabel("Chicken Hb β Chain (Residues 1-147)")
plt.title("Dot Plot: Human vs. Chicken Hemoglobin β Chain")
plt.grid(True, linestyle='--', alpha=0.5)
plt.show()