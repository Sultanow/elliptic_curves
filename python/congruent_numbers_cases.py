# https://www.overleaf.com/read/prsvxqtdwthm
import multiprocessing as mp, traceback, math, threading

limA = 1000
limB = 60

def Factor(n):
    import sympy
    fs = []
    for k, v in sympy.factorint(n).items():
        fs.extend([k] * v)
    return sorted(fs)

def PQs():
    for n in range(1 << 20):
        fs = Factor(n)
        if len(fs) != 2:
            continue
        if fs[0] == fs[1]:
            continue
        if fs[0] == 2 or fs[1] == 2:
            continue
        yield (fs[0], fs[1])

def Task(ipq, *, ab = [(a, b) for a in range(-limA + 1, limA) for b in range(1, limB) if a != 0 and math.gcd(a, b) == 1]):
    try:
        i, (p, q) = ipq
        cases = {}
        for a, b in ab:
            def Case(i, cond):
                if cond:
                    cases.setdefault(i, []).append((a, b))
            Case(1, q == a ** 2 + q * b ** 4)
            Case(2, p * q == a ** 2 + b ** 4)
            Case(3, p * q * 4 * b ** 4 == a ** 2 + 1)
            Case(4, p == q * b ** 4 - a ** 2)
            Case(5, p * 4 * b ** 3 == 2 * a ** 2 + q * b)
            Case(6, p * q == b ** 4 - a ** 2)
        for v in cases.values():
            v.sort()
        cases = dict(sorted(cases.items(), key = lambda e: e[0]))
        return (i, p, q), cases
    except:
        return {'error': traceback.format_exc()}

def Log(*pargs, _state = {'lock': threading.Lock(), 'file': None}, **nargs):
    with _state['lock']:
        if _state['file'] is None:
            _state['file'] = open('ec_pq_find_rational_xy_conditions_log.txt', 'a', encoding = 'utf-8')
        print(*pargs, **nargs, file = _state['file'], flush = True)

def Main():
    import sympy
    with mp.Pool(mp.cpu_count()) as pool:
        cnt, scnt, k, kc, avgA, avgB, maxA, maxB, primeB, primeOddB = [0] * 10
        unprimeBs = {}
        for e in pool.imap(Task, enumerate(PQs())):
            if type(e) is dict:
                assert False, e['error']
            (i, p, q), cases = e
            unprimeB_logged = False
            for _, v in cases.items():
                for a, b in v:
                    maxA = max(abs(a), maxA)
                    maxB = max(b, maxB)
                    avgA += abs(a)
                    avgB += b
                    primeB += int(sympy.isprime(b))
                    t = b
                    while t > 0 and (t & 1) == 0:
                        t >>= 1
                    if t == 1 or sympy.isprime(t):
                        primeOddB += 1
                    else:
                        fst = tuple(Factor(t))
                        unprimeBs[fst] = unprimeBs.get(fst, 0) + 1
                        if not unprimeB_logged:
                            Log('un-prime-B', b, f'({Factor(b)})', f'p={p},q={q}', f' cases {i}:', cases)
                            unprimeB_logged = True
                    k += abs(a) / b
                    kc += 1
            cnt += 1
            if cases:
                scnt += 1
            print(f'{i:>5}: {p:>5}, {q:>5}', '   ', *((cases, f'     solved {round(scnt / cnt * 100, 1)}%', 'avg(a/b)', round(k / max(1, kc), 1), 'maxA', maxA, 'avgA', round(avgA / kc, 1), 'maxB', maxB, 'avgB', round(avgB / kc, 1), 'non-primeB', kc - primeB, 'prime-B', primeB, f'({round(primeB / kc * 100, 1)}%)', 'prime-odd-B', primeOddB, f'({round(primeOddB / kc * 100, 1)}%)', 'un-prime-Bs', dict(sorted(unprimeBs.items(), key = lambda e: e[0]))) if cases else ()), flush = True)

if __name__ == '__main__':
    Main()