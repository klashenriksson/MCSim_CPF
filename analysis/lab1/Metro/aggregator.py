import os

rows = []

for filepath in os.listdir("results/"):
    with open(f"results/{filepath}", "r") as f:
        for line in f:
            T, L, e, Cv, m = line.split()
            rows.append({
                "L":L,
                "T":T,
                "e":e,
                "Cv":Cv,
                "m":m
            })

with open(f"aggregated_results.txt", "w") as f:
    f.write("L T e Cv m\n")
    for row in rows:
        L = row["L"]
        T = row["T"]
        e = row["e"]
        Cv = row["Cv"]
        m = row["m"]
        f.write(f"{L} {T} {e} {Cv} {m}\n")
