import json
import subprocess
import os
import sys
import platform

def get_script_path():
    return os.path.dirname(os.path.realpath(sys.argv[0]))

with open(sys.argv[1], "r") as config_file:
    data = json.load(config_file)
    filepath = os.path.join(get_script_path(), "rw.exe" if platform.system() == "Windows" else "rw")
    for config in data["configs"]:
        Ns = config["N"]
        nblock = config["nblock"]

        for N in Ns:
            params = [
                filepath,
                f"nblock={nblock}",
                f"N={N}",
                "run"
            ]
            print(f"Executing {params}")
            subprocess.Popen(params, stdout=subprocess.DEVNULL, cwd=get_script_path())