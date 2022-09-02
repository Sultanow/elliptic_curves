# https://www.overleaf.com/read/prsvxqtdwthm
import multiprocessing as mp, traceback, math

limA = 1000
limB = 60

def PQs():
    import sympy
    for n in range(1 << 20):
        fs = []
        for k, v in sympy.factorint(n).items():
            fs.extend([k] * v)
        fs = sorted(fs)
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

def Main():
    with mp.Pool(mp.cpu_count()) as pool:
        k, kc, avgA, avgB, maxA, maxB = 0, 0, 0, 0, 0, 0
        for e in pool.imap(Task, enumerate(PQs())):
            if type(e) is dict:
                assert False, e['error']
            (i, p, q), cases = e
            for _, v in cases.items():
                for a, b in v:
                    maxA = max(abs(a), maxA)
                    maxB = max(b, maxB)
                    avgA += abs(a)
                    avgB += b
                    k += abs(a) / b
                    kc += 1
            print(f'{i:>5}: {p:>5}, {q:>5}', '   ', *((cases, '     avg(a/b)', round(k / max(1, kc), 1), 'maxA', maxA, 'avgA', round(avgA / kc, 1), 'maxB', maxB, 'avgB', round(avgB / kc, 1)) if cases else ()))

if __name__ == '__main__':
    Main()