from asyncio import subprocess
import json
import subprocess

with open("sim_confs.json", "r") as config_file:
    data = json.load(config_file)
    for config in data["configs"]:
        Ls = config["L"]
        Ts = config["T"]

        for L in Ls:
            for T in Ts:
                params = [
                    f"ising.exe",
                    f"L={L}",
                    f"T={T}",
                    "run"
                ]

                print(f"Executing {params}")
                subprocess.Popen(params, stdout=subprocess.DEVNULL)