with open("q2.fasta", "r") as f:
    pdbs = list(set([line[1:5] for line in f if line.startswith(">")]))

with open("pisces_list.txt", "w") as f:
    f.write("\n".join(pdbs))
print("Done! ")