r"""
Full identities
"""

from .max_plus_symbolic import symbolic_max_plus_matrices

def full_candidates(d, n, verbose=False):
    a,b = symbolic_max_plus_matrices(d,2)
    m = [a]
    i = [0]
    d = {}
    while True:
        for _ in xrange(n-len(i)):
            i.append(0)
            m.append(m[-1] * a)

        h = hash(m[-1])
        if h in d:
            if verbose:
                print "POTENTIAL IDENTITY: {} = {}".format(d[h][0], i)
            d[h].append(i[:])
        else:
            d[h] = [i[:]]

        while i[-1] == 1:
            i.pop()
            m.pop()
        if len(i) == 1:
            break
        i.pop()
        m.pop()
        i.append(1)
        m.append(m[-1] * b)
    return [v for v in d.itervalues() if len(v) > 1]
