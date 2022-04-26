import os
import sys

def get_script_path():
    return os.path.dirname(os.path.realpath(sys.argv[0]))

rows = []

corr_dir = os.path.join(get_script_path(), "correlations")
for filepath in os.listdir(corr_dir):
    with open(os.path.join(corr_dir, f"{filepath}"), "r") as f:
        for line in f:
            T, L, t, C = line.split()
            rows.append({
                "L":L,
                "T":T,
                "t":t,
                "C":C,
            })

with open(os.path.join(get_script_path(), f"aggregated_correlations.txt"), "w") as f:
    f.write("L T t C\n")
    for row in rows:
        L = row["L"]
        T = row["T"]
        t = row["t"]
        C = row["C"]
        f.write(f"{L} {T} {t} {C}\n")
