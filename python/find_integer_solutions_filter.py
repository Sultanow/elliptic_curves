# Filter-based solution by Artry (published at Stack Overflow):
# https://stackoverflow.com/questions/70664185/preventing-overflow-of-large-integers-in-gpu-optimized-methods-such-as-gmpy2-a/70810780#70810780
import numpy as np, numba

@numba.njit('u8[:](u8[:], u8, u8, u1[:])', cache = True, parallel = True)
def do_filt(x, i, K, filt):
    x += i; x %= K
    x2 = x
    x2 *= x2;     x2 %= K
    x6 = x2 * x2; x6 %= K
    x6 *= x2;     x6 %= K
    x6 += np.uint64(4 * K + 4)
    x2 <<= np.uint64(2)
    x6 -= x2; x6 %= K
    y = x6
    #del x2
    filt_y = filt[y]
    filt_y_i = np.flatnonzero(filt_y).astype(np.uint64)
    return filt_y_i

def main():
    import math
    gmpy2 = None
    import gmpy2
    
    Int = lambda x: (int(x) if gmpy2 is None else gmpy2.mpz(x))
    IsSquare = lambda x: gmpy2.is_square(x)
    Sqrt = lambda x: Int(gmpy2.sqrt(x))
    
    Ks = [2 * 2 * 3 * 5 * 7 * 11 * 13 * 17 * 19,    23 * 29 * 31 * 37 * 41]
    filts = []
    for i, K in enumerate(Ks):
        a = np.arange(K, dtype = np.uint64)
        a *= a
        a %= K
        filts.append((K, np.zeros((K,), dtype = np.uint8)))
        filts[-1][1][a] = 1
        print(f'filter {i} ratio', round(len(np.flatnonzero(filts[-1][1])) / K, 4))
    
    limit = 1 << 30
    block = 1 << 26
    
    for i in range(0, limit, block):
        print(f'i block {i // block:>3} (2^{math.log2(i + 1):>6.03f})')
        x = np.arange(0, min(block, limit - i), dtype = np.uint64)
        
        for ifilt, (K, filt) in enumerate(filts):
            len_before = len(x)
            x = do_filt(x, i, K, filt)
            print(f'squares filtered by filter {ifilt}:', round(len(x) / len_before, 4))
        
        x_to_check = x
        print(f'remain to check {len(x_to_check)}')
        
        sq_x = []
        for x0 in x_to_check:
            x = Int(i + x0)
            y = x ** 6 - 4 * x ** 2 + 4
            if not IsSquare(y):
                continue
            yr = Sqrt(y)
            assert yr * yr == y
            sq_x.append((int(x), int(yr)))
        print('squares found', len(sq_x))
        print(sq_x)
        
        del x

if __name__ == '__main__':
    main()