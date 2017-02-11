r"""
Explore identities for full 3x3 matrices by substituting 2x2 full identities
into 3x3 vv identities

For each pair (u,v) of 2x2 full identities a file
../candidates_full/killed_{u}_{v} is created that contains the list of
pairs (s,t) of 3x3 vv identities for which the program has found a
counterexample.

To run a parallel test, just use the function ``test_all`` as in::

    sage: %runfile RechercheI3_vincent.py  # not tested
    sage: test_all(2**20, 1000)            # not tested
"""
from max_plus import *
from max_plus.max_plus_int import IntegerMatrixProduct, IntegerMaxPlusMatrix

def t2s(t): return ''.join(map(str,t))
def s2t(s): return tuple(map(int,s))

def read_identities(filename):
    r"""
    Return the list of identities read from the file with name ``filename``
    """
    try:
        f = open(filename, 'r')
    except IOError:
        return []

    identities = []
    for line in f.readlines():
        identities.append(tuple(map(s2t, line.split())))
    f.close()
    return identities

def print_candidates(k1=17, k2=20, filename=None):
    r"""
    Print the candidate identities.

    If an argument is provided, the result is written into that file. Otherwise,
    the output is printed on the screen.

    Each line corresponds to a quadruple ``s t u v`` in this order.
    """
    Id2 = []
    for l in range(k1,k2):
        Id2.extend(read_identities('../full_identities/2_' + str(l)))
    Id3 = []
    for l in range(22,24):
        Id3.extend(read_identities('../vv_identities/3_' + str(l)))

    if filename is None:
        import sys
        output = sys.stdout
    else:
        output = open(filename, 'w')

    for u,v in Id2:
        killedf = '../candidates_full/killed_{}_{}'.format(t2s(u), t2s(v))
        killed = set(read_identities(killedf))
        for s,t in set(Id3).difference(killed):
            output.write(t2s(s) + ' ' + t2s(t) + ' ' + t2s(u) + ' ' + t2s(v) + '\n')
    if filename is not None:
        output.close()

def is_relation_stuv(s, t, u, v, K=2**20, num=100, p=0.01):
    Pu = IntegerMatrixProduct(u)
    Pv = IntegerMatrixProduct(v)

    while True:
        a = random_integer_max_plus_matrix(3,-K,K,p)**6
        b = random_integer_max_plus_matrix(3,-K,K,p)**6
        U = Pu(a,b)
        V = Pv(a,b)

        if U != V and not is_relation(s, t, [(U,V)], False):
            return False

def test_uv(u, v, Id3, K=2**20, num=100, p=0.01):
    r"""
    Check pairs obtained from a given ``(u,v)`` of 2x2 full identities.

    The file ``killed_{u}_{v}`` is updated accordingly.

    INPUT:

    - ``u``, ``v`` -- the 2x2 identity

    - ``Id3`` -- list of vv identities

    - ``K`` -- bound for matrix entries

    - ``num`` -- number of matrices to test

    - ``p`` -- probability for minus infinity to appear as coefficient
    """
    Pu = IntegerMatrixProduct(u)
    Pv = IntegerMatrixProduct(v)
    
    filename = '../candidates_full/killed_{}_{}'.format(t2s(u), t2s(v))
    killed = set(read_identities(filename))

    Id3_to_test = set(Id3).difference(killed)
    if not Id3_to_test:
        return u,v,-1

    eltsUV = []
    while len(eltsUV) < num:
        a = random_integer_max_plus_matrix(3,-K,K,p)**6
        b = random_integer_max_plus_matrix(3,-K,K,p)**6
        U = Pu(a,b)
        V = Pv(a,b)
        if U != V:
            eltsUV.append((U,V))

    nkilled = 0
    for s,t in Id3_to_test:
        if not is_relation(s, t, eltsUV, False):
            nkilled += 1
            killed.add((s,t))

    killed = sorted(killed)
    f = open(filename, 'w')
    for s,t in killed:
        f.write(t2s(s) + ' ' + t2s(t) + '\n')
    f.close()
    return u,v,nkilled

def test_all(k1, k2, K, num, ncpus=None, verbose=False, logfile=None):
    r"""
    Run through all 2x2 full identities in parallel

    EXAMPLES::

        sage: import sys                                 # not tested
        sage: %runfile RechercheI3_vincent.py            # not tested
        sage: test_all(2**10, 1000, logfile=sys.stdout) # not tested
    """
    import multiprocessing as mp
    from max_plus.misc import parallel_unfold

    close = False
    if isinstance(logfile, str):
        close = True
        logfile = open(logfile, 'w')
    if ncpus is None:
        ncpus = mp.cpu_count()
    pool = mp.Pool(ncpus)

    Id2 = []
    for l in range(k1,k2):
        Id2.extend(read_identities('../full_identities/2_' + str(l)))

    Id3 = []
    for l in range(22,24):
        Id3.extend(read_identities('../vv_identities/3_' + str(l)))

    tasks = ((verbose, test_uv, u, v, Id3, K, num) for u,v in Id2)
    total_killed = 0
    for u,v,nkilled in pool.imap_unordered(parallel_unfold, tasks):
        if logfile is not None:
            logfile.write("{} identities killed with u={} v={}\n".format(nkilled,t2s(u),t2s(v)))
            logfile.flush()
        total_killed += nkilled
    pool.close()
    pool.join()
    print "{} identities killed".format(total_killed)
    if close:
        logfile.close()
    return total_killed
