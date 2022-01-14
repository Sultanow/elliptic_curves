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
def findIntegerSolutionsGmpy2(limit: np.uint64):
    for x in np.arange(0, limit+1, dtype=np.uint64):
        x_ = mpz(int(x))**2
        y = x_**3-mpz(4)*x_+mpz(4)
        if gmpy2.is_square(y):
            print([x,gmpy2.sqrt(y),y])

def main() -> int:
    limit = 1000000000
    start = time.time()
    findIntegerSolutionsGmpy2(limit)
    end = time.time()
    print("Time elapsed: {0}".format(end - start))
    return 0

if __name__ == '__main__':
    sys.exit(main())