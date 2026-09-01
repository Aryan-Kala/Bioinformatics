import matplotlib.pyplot as plt

human_20 = "MVHLTPEEKSAVTALWGKVN"
chicken_20 = "MVHWTAEEKQLITGLWGKVN"

x_coords = []
y_coords = []

for i, r1 in enumerate(human_20):
    for j, r2 in enumerate(chicken_20):
        if r1 == r2:
            x_coords.append(i + 1)
            y_coords.append(j + 1)

plt.figure(figsize=(6, 6))
plt.scatter(x_coords, y_coords, color='black', s=80)
plt.xlim(0.5, 20.5)
plt.ylim(0.5, 20.5)
plt.xticks(range(1, 21))
plt.yticks(range(1, 21))
plt.grid(True, linestyle='--', alpha=0.6)
plt.xlabel("Human Hb β Chain (1-20)")
plt.ylabel("Chicken Hb β Chain (1-20)")
plt.title("Dot Plot Verification: Residues 1-20")
plt.show()