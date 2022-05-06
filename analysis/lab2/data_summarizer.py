import os
import sys
import subprocess
import platform

def get_script_path():
    return os.path.dirname(os.path.realpath(sys.argv[0]))

args = sys.argv[1:]
files = {arg: [] for arg in args}
summary_filepath = os.path.join(get_script_path(), "summary.exe" if platform.system() == "Windows" else "summary")

if not os.path.isdir(os.path.join(get_script_path(), "data_summary")):
    os.mkdir(os.path.join(get_script_path(), "data_summary"))

for filepath in os.listdir("data/"):
    for arg in args:
        if arg in filepath:
            fullpath = os.path.join("data", filepath)
            files[arg].append(fullpath)

for arg in files.keys():
    with open(os.path.join("data_summary", f"{arg}.txt"), "w+") as f:
        process = subprocess.Popen([summary_filepath], text=True, stdin=subprocess.PIPE, stdout=f, stderr=subprocess.PIPE, cwd=get_script_path())
        input = '\n'.join([str(path) + "\n" for path in files[arg]])
        process.communicate(input=input)