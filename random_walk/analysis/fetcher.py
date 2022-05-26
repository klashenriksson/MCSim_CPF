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

    run_sequential = data["sequential"] if "sequential" in data else False
    for config in data["configs"]:
        nblock = config["nblock"]
        nsamp = config["nsamp"]
        Ns = []
        if "N_interval" in config:
            N_interval = config["N_interval"]
            Ns = range(N_interval[0], N_interval[1]+1)
        else:
            Ns = config["N"]


        if run_sequential:
            params = [
                filepath,
                f"nblock={nblock}",
                f"nsamp={nsamp}",
            ]
            
            for N in Ns:
                params.append(f"N={N}")
                params.append("run")

            print(f"Executing {params}")
            subprocess.Popen(params, stdout=subprocess.DEVNULL, cwd=get_script_path())

        else:
            n = int(np.ceil(len(Ns)/n_chunks))
            chunks = [Ns[i:i+n] for i in range(0, len(Ns), n)]
            processes = []
            for chunk in chunks:
                params = [
                    filepath,
                    f"nblock={nblock}",
                    f"nsamp={nsamp}",
                ]

                for N in chunk:
                    params.append(f"N={N}")
                    params.append("run")

                print(f"Executing {params}")
                process = subprocess.Popen(params, stdout=subprocess.DEVNULL, cwd=get_script_path())
                processes.append(process)
            
            for process in processes:
                process.wait()