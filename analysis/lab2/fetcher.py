import json
import subprocess
import os
import sys
import platform

def get_script_path():
    return os.path.dirname(os.path.realpath(sys.argv[0]))

with open(sys.argv[1], "r") as config_file:
    data = json.load(config_file)
    filepath = os.path.join(get_script_path(), "ising.exe" if platform.system() == "Windows" else "ising")
    for config in data["configs"]:
        Ls = config["L"]
        Ts = config["T"]
        nblock = config["nblock"]

        for L in Ls:
            params = [
                filepath,
                f"nblock={nblock}",
                f"L={L}"
            ]

            for T in Ts:
                params.append(f"T={T}")
                params.append("run")

            print(f"Executing {params}")
            subprocess.Popen(params, stdout=subprocess.DEVNULL, cwd=get_script_path())