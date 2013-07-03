import subprocess
import pandas as pd

subprocess.call("python ghan_test.py", shell=True)

g = pd.read_csv("fort.100", header=0)
print g