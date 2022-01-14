import numpy as np
from numba import jit
import gmpy2
from gmpy2 import mpz, xmpz

import time
import sys

@jit('void(uint64)')
def findIntegerSolutions(limit: np.uint64):
    for x in np.arange(0, limit+1, dtype=np.uint64):
        y = np.uint64(x**6-4*x**2+4)
        sqr = np.uint64(np.sqrt(y))
        if np.uint64(sqr*sqr) == y:
            print([x,sqr,y])

@jit('void(uint64)')
def findIntegerSolutions4(limit: np.uint64):
    for x in np.arange(1, limit+1, dtype=np.uint64):
        y = np.uint64(7**x+x**4+47)
        sqr = np.uint64(np.sqrt(y))
        if np.uint64(sqr*sqr) == y:
            print([x,sqr,y])

@jit('void(uint64)')
def findIntegerSolutionsGmpy2(limit: np.uint64):
    for x in np.arange(0, limit+1, dtype=np.uint64):
        y = mpz(x**6-4*x**2+4)
        if gmpy2.is_square(y):
            print([x,gmpy2.sqrt(y),y])

def main() -> int:

    print(7**426768035+426768035**4+47)

    limit = 1000000000
    start = time.time()
    #findIntegerSolutions4(limit)
    end = time.time()
    print("Time elapsed: {0}".format(end - start))
    return 0

if __name__ == '__main__':
    sys.exit(main())