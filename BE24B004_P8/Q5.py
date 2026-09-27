f = open("Q4.fasta", mode="r")
content = f.readlines()

sequences = {content[i][1::].strip(): content[i + 1].strip() for i in range(0, len(content) - 1, 2)}
f.close()

for header, sequence in sequences.items():
    
    pattern1_length = 6
    pattern2_length = 14
    
    # PATTERN 1 
    for i in range(len(sequence) - pattern1_length + 1):
        sub_sequence = sequence[i:i + pattern1_length]
        if (sub_sequence[0] == "S" or sub_sequence[0] == "V") and (sub_sequence[1] == "T") and \
           (sub_sequence[2] == "V" or sub_sequence[2] == "T") and \
           (sub_sequence[3] == "D" or sub_sequence[3] == "E" or sub_sequence[3] == "R" or sub_sequence[3] == "K") and \
           (sub_sequence[4] == "D" or sub_sequence[4] == "E" or sub_sequence[4] == "R" or sub_sequence[4] == "K") and \
           (sub_sequence[5] != "I" and sub_sequence[5] != "L"):
            
            print(f"Match for Pattern 1 found in {header} at position {i + 1}")
            
    #  PATTERN 2 
    for i in range(len(sequence) - pattern2_length + 1):
        sub_sequence = sequence[i:i + pattern2_length]
        if (sub_sequence[0] in ["F", "I", "L", "V"]) and (sub_sequence[1] == "Q") and \
           (sub_sequence[5] not in ["R", "K"]) and (sub_sequence[6] == "G") and \
           (sub_sequence[10] in ["R", "K"]) and \
           (sub_sequence[13] in ["F", "I", "L", "V", "W", "Y"]):
            
            print(f"Match for Pattern 2 found in {header} at position {i + 1}")