import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
from matplotlib import colors
from sympy import sieve, prime
import itertools
import pandas as pd


matrix_size = 500
matrix = np.zeros((matrix_size, matrix_size))

pq_table = []
with open('../../data/elliptic_curves.csv') as f:
    for line in f.readlines():
        line = line.strip()
        p, q, cases = eval(line)
        cases = [(c, a, b) for c, a, b in cases if c in [2]]
        if p >= 461 and p <= 509:
            idx_p = sieve.search(p)[0]-2
            idx_q = sieve.search(q)[0]-2
            if idx_p < matrix_size and idx_q < matrix_size:
                if len(cases) > 0:
                    matrix[idx_p,idx_q] = 1
                else:
                    pq_table.append([p,q])
                    matrix[idx_p,idx_q] = 0.2

pd.DataFrame(pq_table).to_csv('pq_cases_case_26_c=2.csv', header=False, index=False)

max_prime = prime(matrix_size+2)
axis_labels = list(enumerate(sieve.primerange(3, max_prime)))
axis_labels = axis_labels[::len(axis_labels) // 9][:-1] + [axis_labels[-1]]
ticks = [e[0] for e in axis_labels]
ticklabels = [e[1] for e in axis_labels]

fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(111)
matrix_plot = ax.imshow(matrix.transpose(), origin='lower', interpolation='nearest')
#fig.colorbar(matrix_plot)
ax.set_xticks(ticks); ax.set_xticklabels(ticklabels)
ax.set_yticks(ticks); ax.set_yticklabels(ticklabels)

ax.set_xlabel('p')
ax.set_ylabel('q')

plt.savefig('pattern_case_26_c=2.png')
plt.show()






matrix_size = 500
matrix = np.zeros((matrix_size, matrix_size))

with open('../../data/elliptic_curves.csv') as f:
    for line in f.readlines():
        line = line.strip()
        p, q, cases = eval(line)
        # [1, 2] means: leave only case 0 and 1 (i.e. Case-1 and Case-2)
        cases = [(c, a, b) for c, a, b in cases if c in [1, 2]]
        if len(cases) == 0:
            continue
        idx_p = sieve.search(p)[0]-2
        idx_q = sieve.search(q)[0]-2
        if idx_p < matrix_size and idx_q < matrix_size:
            cases_set = set()
            for case in cases:
                cases_set.add(case[0])
            cases_list = list(cases_set)
            val = sum(i*i for i in cases_list)
            matrix[idx_p,idx_q] = val

max_prime = prime(matrix_size+2)
axis_labels = list(enumerate(sieve.primerange(3, max_prime)))
axis_labels = axis_labels[::len(axis_labels) // 9][:-1] + [axis_labels[-1]]
ticks = [e[0] for e in axis_labels]
ticklabels = [e[1] for e in axis_labels]

fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(111)
matrix_plot = ax.imshow(matrix.transpose(), origin='lower', interpolation='nearest')
#fig.colorbar(matrix_plot)
ax.set_xticks(ticks); ax.set_xticklabels(ticklabels)
ax.set_yticks(ticks); ax.set_yticklabels(ticklabels)

ax.set_xlabel('p')
ax.set_ylabel('q')

plt.savefig('pattern_case_26_c=1_c=2.png')
plt.show()