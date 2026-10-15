import math

sequence = (
    "LGASGIAAFAFGSTAILIILFNMAAEVHFDPLQFFRQFFWLGLYPPKAQY"
    "GMGIPPLHDGGWWLMAGLFMTLSLGSWWIRVYSRARALGLGTHIAWNFAA"
    "AIFFVLCIGCIHPTLVGSWSEGVPFGIWPHIDWLTAFSIRYGNFYYCPWH"
    "GFSIGFAYGCGLLFAAHGATILAVARFGGDREIEQITDRGTAVERAALFW"
)

sec_str = (
    "XHHHHHHHHHHHHHHHHHHHHHHHXXXXXXXXXXXXXXXXXXXXXXXXXX"
    "XXXXXXXXXXHHHHHHHHHHHHHHHHHHHHHHHHHHHHXXHHHHHHHHHH"
    "HHHHHHHHHXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXHHH"
    "HHHHHHHHHHHHHHHHHHHHHHHXXXXXXXXXXXXXXXXXXXXXXXXXXX"
)

sec_str = sec_str.replace('Η', 'H').replace('Χ', 'X')

total_N = len(sequence)
total_H = sec_str.count('H')

print(f"{'Amino Acid':<12} | {'Propensity (P)':<16} | {'Info Content (I)'}")
print("-" * 50)

for aa in sorted(set(sequence)):
    total_aa = sequence.count(aa)
    
    aa_in_H = sum(1 for i in range(total_N) if sequence[i] == aa and sec_str[i] == 'H')
    
    if total_aa > 0 and total_H > 0:
        propensity = (aa_in_H / total_H) / (total_aa / total_N)
    else:
        propensity = 0
        
    info_content = math.log2(propensity) if propensity > 0 else 0
    
    print(f"{aa:<12} | {propensity:<16.4f} | {info_content:.4f}")