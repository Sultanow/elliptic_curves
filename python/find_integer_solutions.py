import numpy as np
from numba import jit
import sys

@jit('void(uint64)')
def findIntegerSolutions(limit: np.uint64):
    for x in np.arange(0, limit, dtype=np.uint64):
        y = np.uint64(x**6-4*x**2+4)
        sqr = np.uint64(np.sqrt(y))
        if np.uint64(sqr*sqr) == y:
            print([x,sqr,y])

def main() -> int:
    limit = 100000000
    findIntegerSolutions(limit)
    return 0

if __name__ == '__main__':
    sys.exit(main())