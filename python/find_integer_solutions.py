from numba.cuda import target
import numpy as np
from numba import jit, cuda
import gmpy2
from gmpy2 import mpz, xmpz

import time
import sys

@jit('void(uint64)')
def find_integer_solutions(limit: np.uint64):
    for x in np.arange(0, limit+1, dtype=np.uint64):
        y = np.uint64(x**6-4*x**2+4)
        sqr = np.uint64(np.sqrt(y))
        if np.uint64(sqr*sqr) == y:
            print([x,sqr,y])

@jit('void(uint64)')
def find_integer_solutions_gmpy2(limit: np.uint64):
    for x in np.arange(0, limit+1, dtype=np.uint64):
        x_ = mpz(int(x))**2
        y = x_**3-mpz(4)*x_+mpz(4)
        if gmpy2.is_square(y):
            print([x,gmpy2.sqrt(y),y])

def main() -> int:
    print(cuda.gpus)

    limit = 1000000000
    start = time.time()
    find_integer_solutions(limit)
    end = time.time()
    print("Time elapsed: {0}".format(end - start))
    return 0

if __name__ == '__main__':
    sys.exit(main())
