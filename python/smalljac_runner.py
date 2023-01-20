import subprocess
import sys
import csv
import timeit

def main() -> int:
    #curvesfile = 'C:/Users/esultano/git/elliptic_curves/data/candidates_small.txt'
    curvesfile = './candidates_small.txt'
    numcurves = sum(1 for _ in open(curvesfile))
    result = [None]*numcurves

    start = timeit.default_timer()
    with open(curvesfile) as f:
        for i, line in enumerate(f.readlines()):
            line = line.strip()
            curve = eval(line)
            del curve[0]
            del curve[-1]
            # call external process (smalljac) to obtain the sums and append them to the line
            cmd_line = "./smalljac/lpdata " + "testfile " + '"' + str(curve) + '"' + " 2e17"
            processoutput = subprocess.run(
                cmd_line, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE
            )           
            sums = eval(processoutput.stdout)
            
            # cmd1 = [r'C:/Users/esultano/git/elliptic_curves/python/smalljac_emulator.bat']
            # p = subprocess.Popen(cmd1,stdin=subprocess.PIPE,stdout=subprocess.PIPE)
            # output = p.communicate(input='x'.encode())[0]
            # sums = eval(output.decode('ascii').strip())
            
            #result.insert(i, curve+sums)
            result[i] = curve+sums
    end = timeit.default_timer()
    elapsed_time = round((t_1 - t_0) / 60, 3)
    print(f"Elapsed time: {elapsed_time} mins")
    
    with open(curvesfile + '.result.csv', 'w+', newline='') as result_csv:
        csvWriter = csv.writer(result_csv, delimiter=',')
        csvWriter.writerows(result)
        return 0

if __name__ == '__main__':
    sys.exit(main())
