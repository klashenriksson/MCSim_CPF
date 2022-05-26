import json
import subprocess
import os
import sys
import platform
import numpy as np

def get_script_path():
    return os.path.dirname(os.path.realpath(sys.argv[0]))

rw_name = sys.argv[1]
with open(sys.argv[2], "r") as config_file:
    data = json.load(config_file)
    filepath = os.path.join(get_script_path(), rw_name)
    n_chunks = 10
    for config in data["configs"]:
        nblock = config["nblock"]
        nsamp = config["nsamp"]
        Ls = config["L"]
        Ps = config["p"]

        for L in Ls:
            params = [
                filepath,
                f"nblock={nblock}",
                f"nsamp={nsamp}",
                f"nblock={nblock}",
            ]

            for p in Ps:
                params.append(f"L={L}")
                params.append(f"p={p}")
                params.append("run")

            print(f"Executing {params}")
            subprocess.Popen(params, stdout=subprocess.DEVNULL, cwd=get_script_path())