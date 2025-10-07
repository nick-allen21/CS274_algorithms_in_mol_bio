# run as a scratch script
import re
from difflib import SequenceMatcher

def load_secstr(p):
    with open(p) as f:
        lines = [l.rstrip("\n") for l in f]
    sec = next(l for l in lines if l.startswith(">secstr")+1)  # next line after >secstr
    # If your editor doesn’t support that, do:
    # idx = lines.index(">secstr"); sec = lines[idx+1]
    return sec

h = load_secstr("/Users/nickallen/Documents/GitHub/-CS274-Algorithms-in-Molecular-Biology/Project1/quiz_input/quiz4_struct_human.txt")
b = load_secstr("/Users/nickallen/Documents/GitHub/-CS274-Algorithms-in-Molecular-Biology/Project1/quiz_input/quiz4_struct_bacteria.txt")

# strip spaces – annotations use spaces as separators
h, b = h.replace(" ", ""), b.replace(" ", "")

sm = SequenceMatcher(None, h, b)
largest = (0, 0, 0)  # (h_start, b_start, length) where human has a gap (bacterial insertion)
for tag, i1, i2, j1, j2 in sm.get_opcodes():
    if tag == "insert":  # inserted in b vs h
        if (j2 - j1) > largest[2]:
            largest = (i1, j1, j2 - j1)

i1, j1, L = largest
print(f"Largest bacterial insertion: length {L}, at bacterial secstr positions {j1}..{j1+L-1}, aligned after human secstr index {i1-1}")