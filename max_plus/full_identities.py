r"""
Full identities
"""

from itertools import izip

from .sage_import import *

from .max_plus_symbolic import symbolic_max_plus_matrices
from .max_plus_int import (random_integer_max_plus_matrix,
                           integer_max_plus_matrix_identity)

def full_identities_iterator(d, n, num_mats=10):
    a = [random_integer_max_plus_matrix(2, -1000000, 1000000, 0.01) \
            for _ in range(num_mats)]
    b = [random_integer_max_plus_matrix(2, -1000000, 1000000, 0.01) \
            for _ in range(num_mats)]
    m = [a[:]]

    ms = symbolic_max_plus_matrices(2,2)

    i = [0]
    d = {}
    while True:
        assert len(i) == len(m)
        assert len(m[-1]) == num_mats
        for _ in xrange(n-len(i)):
            i.append(0)
            m.append(tuple(mat * e for mat,e in izip(m[-1],a)))

        h = hash(m[-1])
        if h in d:
            m1 = prod(ms[k] for k in i)
            for ii in d[h]:
                m2 = prod(ms[k] for k in ii)
                if m1 == m2:
                    yield (tuple(i),ii)
            d[h].append(tuple(i))
        else:
            d[h] = [tuple(i)]

        while i[-1] == 1:
            i.pop()
            m.pop()
        if len(i) == 1:
            break
        i.pop()
        i.append(1)
        m.pop()
        m.append(tuple(mat * e for mat,e in izip(m[-1],b)))
