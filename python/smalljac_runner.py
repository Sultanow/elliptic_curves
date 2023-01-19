import numpy as np
import subprocess
import sys
import pandas as pd

with open('C:/Users/esultano/git/elliptic_curves/data/candidates_small.txt') as f:
    for line in f.readlines():
        line = line.strip()
        curve = eval(line)
        # call external process (smalljac) to obtain the sums and append them to the line
        # result = subprocess.run(
        #    [sys.executable, "-c", "C:/Users/esultano/git/elliptic_curves/python/smalljac_emulator.bat"], capture_output=True, text=True
        #)
        #sums = result.stdout
        cmd1 = [r'C:/Users/esultano/git/elliptic_curves/python/smalljac_emulator.bat']
        p = subprocess.Popen(cmd1,stdin=subprocess.PIPE,stdout=subprocess.PIPE)
        output = p.communicate(input='x'.encode())[0]
        sums = eval(output.decode('ascii').strip())
        curve+=sums
        print(curve)