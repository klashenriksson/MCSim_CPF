from asyncio import subprocess
import json
import subprocess
import os
import sys
import platform

def get_script_path():
    return os.path.dirname(os.path.realpath(sys.argv[0]))

with open("sim_confs.json", "r") as config_file:
    data = json.load(config_file)
    filepath = os.path.join(get_script_path(), "ising.exe" if platform.system() == "Windows" else "ising")
    for config in data["configs"]:
        Ls = config["L"]
        Ts = config["T"]

        for L in Ls:
            for T in Ts:
                params = [
                    filepath,
                    "nblock=64",
                    f"L={L}",
                    f"T={T}",
                    "run"
                ]

                print(f"Executing {params}")
                subprocess.Popen(params, stdout=subprocess.DEVNULL, cwd=get_script_path())