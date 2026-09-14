Human_Beta = "MVHLTPEEKTAVNALWGKVNVDAVGGEALGRLLVVYPWTQRFFESFGDLSSPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFSQLSELHCDKLHVDPENFRLLGNVLVCVLARNFGKEFTPQMQAAYQKVVAGVANALAHKYH"

Chicken_Beta = "MVHWTAEEKQLITGLWGKVNVAECGAEALARLLIVYPWTQRFFASFGNLSSPTAILGNPMVRAHGKKVLTSFGDAVKNLDNIKNTFSQLSELHCDKLHVDPENFRLLGDILIIVLAAHFSKDFTPECQAAWQKLVRVVAHALARKYH"

min_len = 5
Con_seq = []

for i in range(len(Human_Beta)):
    for j in range(len(Chicken_Beta)):

        if (
            Human_Beta[i] == Chicken_Beta[j]
            and (i == 0 or j == 0 or Human_Beta[i - 1] != Chicken_Beta[j - 1])
        ):
            k = 0

            while (
                i + k < len(Human_Beta)
                and j + k < len(Chicken_Beta)
                and Human_Beta[i + k] == Chicken_Beta[j + k]
            ):
                k += 1

            if k >= min_len:
                Con_seq.append({
                    "segment": Human_Beta[i:i + k],
                    "human_pos": i + 1,
                    "chicken_pos": j + 1,
                    "Seq_len": k
                })

print("Found Conserved Segments:")

if not Con_seq:
    print("No conserved segments found.")
else:
    for segment in Con_seq:
        print(
            "Segment:", segment["segment"],
            "| Human Pos:", segment["human_pos"],
            "| Chicken Pos:", segment["chicken_pos"],
            "| Length:", segment["Seq_len"]
        )