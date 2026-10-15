import numpy as np

sequence = "KVFGRCELAAAMKRHGLDNYRGYSLGNWVCAAKFESNFNTQATNRNTDGSTDYGILQINSRWWCNDGRTPGSRNLCNIPCSALLSSDITASVNCAKKIVSDGNGMNAWVAWRNRCKGTDVQAWIRGCRL"

helix_params = {res: 1 for res in ["E","A","L","H","M","Q","W","V","F"]}
helix_params.update({res: 0.5 for res in ["K","I"]})
helix_params.update({res: 0 for res in ["D","T","S","R","C"]})
helix_params.update({res: -1 for res in ["N","Y","P","G"]})

sheet_params = {res: 1 for res in ["M","V","I","C","Y","F","Q","L","T","W"]}
sheet_params.update({res: 0.5 for res in ["A"]})
sheet_params.update({res: 0 for res in ["R","G","D"]})
sheet_params.update({res: -1 for res in ["K","S","H","N","P","E"]})

helix_propensities = {
    'E': 1.53, 'A': 1.45, 'L': 1.34, 'H': 1.24, 'M': 1.20, 'Q': 1.17, 'W': 1.14, 'V': 1.14, 
    'F': 1.12, 'K': 1.07, 'I': 1.00, 'D': 0.98, 'T': 0.82, 'S': 0.79, 'R': 0.79, 'C': 0.77, 
    'N': 0.73, 'Y': 0.61, 'P': 0.59, 'G': 0.53
}

sheet_propensities = {
    'M': 1.67, 'V': 1.65, 'I': 1.60, 'C': 1.30, 'Y': 1.29, 'F': 1.28, 'Q': 1.23, 'L': 1.22, 
    'T': 1.20, 'W': 1.19, 'A': 0.97, 'R': 0.90, 'G': 0.81, 'D': 0.80, 'K': 0.74, 'S': 0.72, 
    'H': 0.71, 'N': 0.65, 'P': 0.62, 'E': 0.26
}

def find_segments(seq, init_params, prop_params, window_size, init_threshold, ext_window):
    segments = []
    i = 0
    while i <= len(seq) - window_size:
        window_score = sum([init_params[res] for res in seq[i:i+window_size]])
        
        if window_score >= init_threshold:
            start = i
            end = i + window_size
            
            # Extend right
            while end < len(seq):
                last_k_score = sum([prop_params[res] for res in seq[end-ext_window+1:end+1]])
                if last_k_score < init_threshold:
                    break
                end += 1
            
            # Extend left
            while start > 0:
                first_k_score = sum([prop_params[res] for res in seq[start-1:start-1+ext_window]])
                if first_k_score < init_threshold:
                    break
                start -= 1
                
            segments.append([start, end])
            i = end  # Jump past the identified segment
        else:
            i += 1
            
    return segments

helices = find_segments(sequence, helix_params, helix_propensities, 6, 4, 4)
sheets = find_segments(sequence, sheet_params, sheet_propensities, 5, 3, 3)

final_helices = []
final_sheets = []

for h_start, h_end in helices:
    h_set = set(range(h_start, h_end))
    overlap_found = False
    
    for s_idx, (s_start, s_end) in enumerate(sheets):
        s_set = set(range(s_start, s_end))
        common = h_set.intersection(s_set)
        
        if common:
            overlap_found = True
            min_c, max_c = min(common), max(common) + 1
            overlap_seq = sequence[min_c:max_c]
            
            h_score = sum([helix_propensities[res] for res in overlap_seq])
            s_score = sum([sheet_propensities[res] for res in overlap_seq])
            
            if h_score > s_score:
                final_helices.append([h_start, h_end])
                # Truncate sheet by adjusting its start or end
                if s_start < min_c:
                    sheets[s_idx] = [s_start, min_c]
                elif s_end > max_c:
                    sheets[s_idx] = [max_c, s_end]
                else:
                    sheets[s_idx] = [0, 0] # Eliminated entirely
            else:
                final_sheets.append([s_start, s_end])
                # Truncate helix
                if h_start < min_c:
                    final_helices.append([h_start, min_c])
                elif h_end > max_c:
                    final_helices.append([max_c, h_end])
            break 
            
    if not overlap_found:
        final_helices.append([h_start, h_end])

for s_start, s_end in sheets:
    if s_end - s_start > 0 and [s_start, s_end] not in final_sheets:
        final_sheets.append([s_start, s_end])

print("Final Helices:")
for start, end in sorted(final_helices):
    seq_str = sequence[start:end]
    if len(seq_str) > 0:
        prop_sum = sum([helix_propensities[res] for res in seq_str])
        print(f"{seq_str} at index {start+1} with propensity {prop_sum:.3f}")

print("\nFinal Sheets:")
for start, end in sorted(final_sheets):
    seq_str = sequence[start:end]
    if len(seq_str) > 0:
        prop_sum = sum([sheet_propensities[res] for res in seq_str])
        print(f"{seq_str} at index {start+1} with propensity {prop_sum:.3f}")
