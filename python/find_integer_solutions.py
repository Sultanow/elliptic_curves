from numba.cuda import target
import numpy as np
from numba import jit, cuda
import gmpy2
from gmpy2 import mpz

import time
import sys

@jit('void(uint64)')
def find_integer_solutions(limit: np.uint64):
    for x in np.arange(0, limit+1, dtype=np.uint64):
        y = np.uint64(x**6-4*x**2+4)
        sqr = np.uint64(np.sqrt(y))
        if np.uint64(sqr*sqr) == y:
            print([x,sqr,y])

#dauert ewig
@jit('void(uint64)', forceobj = True)
def find_integer_solutions_gmpy2(limit: np.uint64):
    for x in np.arange(0, limit+1, dtype=np.uint64):
        x = mpz(int(x))
        y = mpz(x**6-4*x**2+4)
        if gmpy2.is_square(y):
            print([x,gmpy2.sqrt(y),y])

def find_integer_solutions_gmpy2_for_mp(begin, end, *, cache = {}):
    def call(begin, end):
        for x in range(begin, end):
            x = mpz(int(x))
            y = mpz(x**6-4*x**2+4)
            if gmpy2.is_square(y):
                print([x,gmpy2.sqrt(y),y])
    if 'compiled' not in cache:
        cache['compiled'] = jit('void(uint64, uint64)', forceobj = True)(call)
    cache['compiled'](begin, end)

def find_integer_solutions_gmpy2_multiprocessing(limit: np.uint64):
    import multiprocessing as mp
    limit = int(limit)
    nthreads = mp.cpu_count()
    print(f'Using {nthreads} cores.')
    block = max(1, (limit + 1 + nthreads * 8 - 1) // (nthreads * 8))
    with mp.Pool(nthreads) as pool:
        pool.starmap(find_integer_solutions_gmpy2_for_mp, [(i, min(i + block, limit + 1)) for i in range(0, limit + 1, block)])

#dauert 5min
@jit('void(uint64)', forceobj = True)
def find_integer_solutions_gmpy2_optimized(limit: np.uint64):
    for x in np.arange(0, limit+1, dtype=np.uint64):
        x_ = mpz(int(x))**2
        y = x_**3-mpz(4)*x_+mpz(4)
        if gmpy2.is_square(y):
            print([x,gmpy2.sqrt(y),y])

def main() -> int:
    #print(cuda.gpus)

    limit = 100_000_000
    start = time.time()
    find_integer_solutions_gmpy2_multiprocessing(limit)
    end = time.time()
    print("Time elapsed: {0}".format(end - start))
    return 0

if __name__ == '__main__':
    sys.exit(main())