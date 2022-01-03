# GPU-optimized Python Program
# generating 6-tuple (s,t,u,t+u,t+u-s,t-s) of squares

from math import sqrt
import pandas as pd
import numpy as np

from numba import jit
import time
import sys

@jit('void(uint64)')
def findIntegerSolutions(limit: np.uint64):
    for x in np.arange(2, limit, dtype=np.uint64):
        y = x**6-4*x**2+4
        sqr = np.uint64(sqrt(y))
        if sqr*sqr == y:
            print([x,sqr])
                                
def main() -> int:
    limit = 100000000
    start = time.time()
    findIntegerSolutions(limit)
    end = time.time()
    print("Time elapsed: {0}".format(end - start))
    
    return 0

if __name__ == '__main__':
    sys.exit(main())