import sys, math, threading, multiprocessing, multiprocessing.pool, traceback
# A^3 * y^2 = B^2 * x^3 - A^2 * B^2 * p * q * x

p = 23
q = 29
lim = 2 ** 14
blockAB = 2 ** 3
stripe_begin, stripe_end = 0, 2 ** 7 // blockAB
pq = p * q

def BinarySearch(f, a, b):
    assert a <= b
    if a >= b:
        return b
    begin, end = a, b
    while a + 1 < b:
        m = (a + b - 1) // 2
        if f(m):
            b = m + 1
        else:
            a = m + 1
    assert a + 1 == b, 'bin: ' + str((begin, end))
    if not f(a):
        return b
    return a

def ExpSearch(f, a, b):
    begin, end = a, b
    for bits in range(0, 1 << 10):
        l = 1 << bits
        b0 = a + l - 1
        if b0 >= b:
            b0 = b
            break
        if f(b0):
            b0 += 1
            break
    else:
        assert False
    try:
        return BinarySearch(f, a, b0)
    except Exception:
        assert False, 'exp: ' + str((begin, end)) + '\n' + traceback.format_exc()

def PrintBase(*pargs, get_lines = False, shared, print_lock, **nargs):
    if 'print_lines' not in shared:
        shared['print_lines'] = 0
    if get_lines:
        return shared['print_lines']
    shared['print_lines'] += nargs.get('end', '\n').count('\n')
    with print_lock:
        print(*pargs, **nargs, flush = True)

def PrintResBase(force = False, *, all_results, shared):
    if 'print_res_lines' not in shared:
        shared['print_res_lines'] = 0
    if not force and Print(get_lines = True) - shared['print_res_lines'] < 25:
        return
    Print('Res: [' + ', '.join([f'(x = {x}/{A}, y = {y}/{B})' for x, A, y, B in sorted(all_results.keys())]) + ']')
    shared['print_res_lines'] = Print(get_lines = True)

def ProcessAB(Abegin, Aend, Bbegin, Bend, *, all_results, shared, print_lock):
    PrintRes = lambda *pargs, **nargs: PrintResBase(*pargs, all_results = all_results, shared = shared, **nargs)
    Print(f'AxB = [{str(Abegin).rjust(4)}, {str(Aend).rjust(4)}) x [{str(Bbegin).rjust(4)}, {str(Bend).rjust(4)})')
    for A in range(Abegin, Aend):
        if A == 0:
            continue
        A2 = A ** 2
        A3 = A ** 3
        for B in range(Bbegin, Bend):
            if B == 0:
                continue
            B2 = B ** 2
            x_end = BinarySearch(lambda x: A3 * (lim - 1) ** 2 - B2 * x ** 3 + A2 * B2 * pq * x < 0, 0, lim)
            x_begin = BinarySearch(lambda x: - B2 * x ** 3 + A2 * B2 * pq * x <= 0, -lim, 0)
            assert x_begin <= x_end
            ylast = 1
            for x in range(x_begin, x_end):
                if x == 0:
                    continue
                if math.gcd(A, x) != 1:
                    continue
                x3 = x ** 3
                right = B2 * x3 - A2 * B2 * pq * x
                f = lambda y: A3 * y ** 2 - right
                y = ExpSearch(lambda y: f(y) >= 0, ylast if f(ylast) < 0 else 1, lim)
                if f(y) == 0 and math.gcd(y, B) == 1:
                    Print(f'x = {x}/{A}, y = {y}/{B}')
                    resb = sorted(all_results.keys())
                    all_results[(x // math.gcd(x, A), A // math.gcd(x, A), y // math.gcd(y, B), B // math.gcd(y, B))] = None
                    if sorted(all_results.keys()) != resb:
                        PrintRes(force = True)
                ylast = y
    PrintRes()

def Run(args):
    global Print
    pargs, nargs = args
    Print = lambda *pa, **na: PrintBase(*pa, shared = nargs['shared'], print_lock = nargs['print_lock'], **na)
    
    try:
        ProcessAB(*pargs, **nargs)
        return {}
    except:
        return {'error': traceback.format_exc()}

def Main():
    global Print

    nproc = multiprocessing.cpu_count()
    with multiprocessing.Manager() as manager, multiprocessing.Pool(nproc) as pool:
        print_lock = manager.Lock()
        shared = manager.dict()
        all_results = manager.dict()
        Print = lambda *pa, **na: PrintBase(*pa, shared = shared, print_lock = print_lock, **na)
        PrintRes = lambda *pargs, **nargs: PrintResBase(*pargs, all_results = all_results, shared = shared, **nargs)
        
        try:
            Print(f'p = {p}, q = {q}')
            for istripe in range(stripe_begin, stripe_end):
                Print('Stripe', istripe)
                
                tasks = []
                for iAB in range(0, istripe + 1):
                    tasks.append(((blockAB * iAB, blockAB * (iAB + 1), blockAB * istripe, blockAB * (istripe + 1)), dict(all_results = all_results, shared = shared, print_lock = print_lock)))
                    if iAB < istripe:
                        tasks.append(((blockAB * istripe, blockAB * (istripe + 1), blockAB * iAB, blockAB * (iAB + 1)), dict(all_results = all_results, shared = shared, print_lock = print_lock)))
                
                for r in pool.imap_unordered(Run, tasks):
                    if 'error' in r:
                        assert False, '\n' + r['error']
        finally:
            pool.terminate()
            pool.join()
            PrintRes(force = True)

if __name__ == '__main__':
    Main()