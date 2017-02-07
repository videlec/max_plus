r"""
Explore identities for full 3x3 matrices by substituting 2x2 full identities
into 3x3 vv identities

For each pair (u,v) of 2x2 full identities a file
../candidates_full/killed_{u}_{v} is created that contains the list of
pairs (s,t) of 3x3 vv identities for which the program has found a
counterexample.

To run a parallel check, just use the function ``check_all`` as in::

    sage: %runfile RechercheI3_vincent.py  # not tested
    sage: check_all(2**20, 1000)           # not tested
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
        left, right = line.split()
        identities.append((s2t(left), s2t(right)))
    f.close()
    return identities

Id2 = []
for l in range(18,20):
    Id2.extend(read_identities('../full_identities/2_' + str(l)))

Id3 = []
for l in range(22,24):
    Id3.extend(read_identities('../vv_identities/3_' + str(l)))

def print_candidates(filename=None):
    r"""
    Print the candidate identities.

    If an argument is provided, the result is written into that file. Otherwise,
    the output is printed on the screen.

    Each line corresponds to a quadruple ``s t u v`` in this order.
    """
    if filename is None:
        import sys
        output = sys.stdout
    else:
        output = open(filename, 'w')
    n = 0
    for u,v in Id2:
        killedf = '../candidates_full/killed_{}_{}'.format(t2s(u), t2s(v))
        killed = set(read_identities(killedf))
        for s,t in set(Id3).difference(killed):
            output.write(t2s(s) + ' ' + t2s(t) + ' ' + t2s(u) + ' ' + t2s(v) + '\n')
            n += 1
    if filename is not None:
        output.close()

def check(u, v, K=2**20, num=100, p=0.01):
    r"""
    Check pairs obtained from a given ``(u,v)`` of 2x2 full identities.

    The file ``killed_{u}_{v}`` is updated accordingly.

    INPUT:

    - ``u``, ``v`` -- the 2x2 identity

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
        return

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
            break

    killed = sorted(killed)
    f = open(filename, 'w')
    for s,t in killed:
        f.write(t2s(s) + ' ' + t2s(t) + '\n')
    f.close()
    return u,v,nkilled

def check_all(K, num, ncpus=None, verbose=False, logfile=None):
    r"""
    Run through all 2x2 full identities in parallel

    EXAMPLES::

        sage: import sys                                 # not tested
        sage: %runfile RechercheI3_vincent.py            # not tested
        sage: check_all(2**10, 1000, logfile=sys.stdout) # not tested
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

    tasks = ((verbose, check, u, v, K, num) for u,v in Id2)
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
