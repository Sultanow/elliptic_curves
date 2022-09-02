# https://www.overleaf.com/read/prsvxqtdwthm
import multiprocessing as mp, traceback, math

lim = 200

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

def Task(ipq, *, ab = [(a, b) for a in range(-lim + 1, lim) for b in range(1, lim) if a != 0 and math.gcd(a, b) == 1]):
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
        for e in pool.imap(Task, enumerate(PQs())):
            if type(e) is dict:
                assert False, e['error']
            (i, p, q), cases = e
            print(f'{i:>3}: {p:>3}, {q:>3}')
            if cases:
                print(cases)

if __name__ == '__main__':
    Main()