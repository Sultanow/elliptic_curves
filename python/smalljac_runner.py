from multiprocessing import Pool
import subprocess
import sys
import csv
import timeit

def runtask(curvesfile: str) -> str:
    numcurves = sum(1 for _ in open(curvesfile))
    result = [None]*numcurves
    j = 0
    starttime = timeit.default_timer()
    with open(curvesfile) as f:
        # we parse lines looking like this: [999999999847,[1,1,1,-8578,-313130],-1]
        for i, line in enumerate(f.readlines()):
            line = line.strip()
            curve = eval(line)
            conductor = curve[0]
            coeffs = curve[1]
            # discriminant = curve[2]

            # call external process (smalljac) to obtain the sums and append them to the line
            cmd_line = "./smalljac/lpdata " + "testfile " + '"' + str(coeffs) + '"' + " 2e17"
            processoutput = subprocess.run(
                cmd_line, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE
            )
            retval = processoutput.stdout
            sum3 = 0
            try:
                sum3 = float(retval)
            except Exception as ex:
                print("skipping " + line + " (sum3 was " + retval +")")
            
            if sum3 > 2.9:
                result[j] = [conductor]+coeffs+[sum3]
                j+=1
    endtime = timeit.default_timer()
    elapsed_time = round((endtime - starttime) / 60, 3)
    print(f"Elapsed time: {elapsed_time} mins")
    
    if j > 0:
        result = result[0:j]
        resultfile = curvesfile + '.result.csv'
        with open(resultfile, 'w+', newline='') as result_csv:
            csvWriter = csv.writer(result_csv, delimiter=',')
            csvWriter.writerows(result)    
        return resultfile

def main() -> int:
    #curvesfile = 'C:/Users/esultano/git/elliptic_curves/data/candidates_small.txt'
    
    numer_of_processes = 36
    curvesfiles = ["./subfile-"+str(i).zfill(2) for i in range(numer_of_processes)]

    with Pool(numer_of_processes) as p:
        results = p.map(runtask, curvesfiles)
    print(results)
    

if __name__ == '__main__':
    sys.exit(main())
