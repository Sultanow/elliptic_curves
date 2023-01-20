import subprocess
import sys
import csv

def main() -> int:
    curvesfile = 'C:/Users/esultano/git/elliptic_curves/data/candidates_small.txt'
    numcurves = sum(1 for _ in open(curvesfile))
    result = [None]*numcurves

    with open(curvesfile) as f:
        for i, line in enumerate(f.readlines()):
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
            #result.insert(i, curve+sums)
            result[i] = curve+sums

    with open(curvesfile + '.result.csv', 'w+', newline='') as result_csv:
        csvWriter = csv.writer(result_csv, delimiter=',')
        csvWriter.writerows(result)
        return 0

if __name__ == '__main__':
    sys.exit(main())
