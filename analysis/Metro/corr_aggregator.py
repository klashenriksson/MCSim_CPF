import os

rows = []

for filepath in os.listdir("correlations/"):
    with open(f"correlations/{filepath}", "r") as f:
        for line in f:
            T, L, t, C = line.split()
            rows.append({
                "L":L,
                "T":T,
                "t":t,
                "C":C,
            })

with open(f"aggregated_correlations.txt", "w") as f:
    f.write("L T t C\n")
    for row in rows:
        L = row["L"]
        T = row["T"]
        t = row["t"]
        C = row["C"]
        f.write(f"{L} {T} {t} {C}\n")
