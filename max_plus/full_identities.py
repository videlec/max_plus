r"""
Identities of full matrices
"""

from itertools import izip

from .sage_import import *

from .max_plus_symbolic import symbolic_max_plus_matrices
from .max_plus_int import (random_integer_max_plus_matrix,
                           integer_max_plus_matrix_identity,
                           is_relation)

def full_identities_iterator(d, n, num_mats1=10, num_mats2=500,
        symbolic_check=True):
    r"""
    Return an iterator over identities of length ``n`` over full matrices of
    size ``d``.

    INPUT:

    - ``d`` -- dimension

    - ``n`` -- length

    - ``num_mats1``, ``num_mats2`` -- number of integer matrices used in the
      checks

    - ``symbolic_check`` -- whether a check with symbolic matrices is performed

    EXAMPLES::

        sage: from max_plus import *
        sage: list(full_identities_iterator(2,13,500))
        []
        sage: next(full_identities_iterator(2,17,500))    # long time
        ((0, 0, 1, 1, 1, 0, 1, 0, 1, 0, 0, 0, 1, 1, 1, 0, 0), 
         (0, 0, 1, 1, 1, 0, 0, 0, 1, 0, 1, 0, 1, 1, 1, 0, 0))
    """
    a = [random_integer_max_plus_matrix(2, -1000000, 1000000, 0.01) \
            for _ in range(num_mats1)]
    b = [random_integer_max_plus_matrix(2, -1000000, 1000000, 0.01) \
            for _ in range(num_mats1)]
    elts = [(random_integer_max_plus_matrix(2, -1000000, 1000000, 0.01),
             random_integer_max_plus_matrix(2, -1000000, 1000000, 0.01))
            for _ in range(num_mats2)]
    ms = symbolic_max_plus_matrices(2,2)

    i = [0]
    m = [a[:]]
    d = {}
    while True:
        for _ in xrange(n-len(i)):
            i.append(0)
            m.append(tuple(mat * e for mat,e in izip(m[-1],a)))

        h = hash(m[-1])
        ti = tuple(i)
        if h in d:
            for ii in d[h]:
                m1 = None
                if is_relation(ti, ii, elts, True):
                    if symbolic_check:
                        if m1 is None:
                            m1 = prod(ms[k] for k in ti)
                        m2 = prod(ms[k] for k in ii)
                        if m1 == m2:
                            yield (ti,ii)
                    else:
                        yield (ti,ii)
            d[h].append(ti)
        else:
            d[h] = [ti]

        while i[-1] == 1:
            i.pop()
            m.pop()
        if len(i) == 1:
            break
        i.pop()
        i.append(1)
        m.pop()
        m.append(tuple(mat * e for mat,e in izip(m[-1],b)))
