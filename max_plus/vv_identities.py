from sage_import import *

import itertools
from time import time
from combinat import occurrences, runs, has_all_subwords
from convex_hull import ppl_polytope

def is_vv_identity(left, right, d, W=None, prefix=(), status=False):
    r"""
    Check whether ``(left, right)`` is a `B^{vv}_d` identity.

    INPUT:

    - ``left``, ``right`` -- the two words of the identity

    - ``d`` -- the dimension (positive integer)

    - ``W`` -- an optional set of finite words

    - ``prefix`` -- an optional prefix to do a partial check

    - ``status`` -- an optional boolean. If ``True`` outputs a pair of made of a
      boolean and a string.

    EXAMPLES::

        sage: from max_plus import *

        sage: p,s = vincent_sv_prefix_suffix(4)
        sage: u = (p+'x'+s).replace('x', 'ab').replace('y', 'ba')
        sage: v = (p+'y'+s).replace('x', 'ab').replace('y', 'ba')
        sage: is_vv_identity(u, v, 4)
        True

        sage: u = (p+'x'+s).replace('x','a').replace('y','b')
        sage: v = (p+'y'+s).replace('x','a').replace('y','b')
        sage: is_vv_identity(u, v, 4)
        False

        sage: F = FiniteWords([0,1])
        sage: t = WordMorphism({0: [0,1], 1: [1,0]}, domain=F)
        sage: p = F((0,1,1,0,0))
        sage: s = F((1,1,0,0,1))
        sage: u = t(p) + t(0) + t(s)
        sage: v = t(p) + t(1) + t(s)
        sage: is_vv_identity(u, v, 3, W=F)
        True
        sage: ans, info = is_vv_identity(u, v, 3, W=F, status=True)
        sage: print info
        u = 00
        num ext. occ.: 45
        num int. occ.: 10
        num faces    : 13
        num verts    : 16
        polytope computation in ...secs
        containance test in ...secs
        <BLANKLINE>
        u = 01
        num ext. occ.: 50
        num int. occ.: 11
        num faces    : 12
        num verts    : 15
        polytope computation in ...secs
        containance test in ...secs
        <BLANKLINE>
        u = 10
        num ext. occ.: 50
        num int. occ.: 10
        num faces    : 12
        num verts    : 15
        polytope computation in ...secs
        containance test in ...secs
        <BLANKLINE>
        u = 11
        num ext. occ.: 45
        num int. occ.: 10
        num faces    : 13
        num verts    : 16
        polytope computation in ...secs
        containance test in ...secs

        sage: all(is_vv_identity(t(i[0]), t(i[1]), 3, W=F) for i in sv_identities(11, 3))
        True
    """
    if W is None:
        alphabet = sorted(set(left).union(right))
        W = FiniteWords(alphabet)
    else:
        alphabet = W.alphabet()

    if len(alphabet) != 2:
        raise ValueError("must be on two letters")

    a,b = alphabet

    left = W(left)
    right = W(right)
    n = len(left)

    if left == right:
        return True

    if n != len(right) or \
       left.count(a) != right.count(a) or \
       left[0] != right[0] or \
       left[1] != right[1] or \
       not has_all_subwords(left, d-1) or \
       not has_all_subwords(right, d-1) or \
       left.factor_set(d-1) != right.factor_set(d-1):
        return False

    # we need 2(d-1) runs in common prefix/suffix
    lruns = runs(left)
    rruns = runs(right)
    if len(lruns) < 4*(d-1) or len(rruns) < 4*(d-1) or \
        lruns[:2*(d-1)] != rruns[:2*(d-1)] or \
        lruns[-2*(d-1):] != rruns[-2*(d-1):]:
        return False

    # maximal run should appear in same first and last positions
    ilmax = max(lruns)
    irmax = max(rruns)
    if ilmax != irmax:
        return False
    il = 0
    while lruns[il] != ilmax: il += 1
    ir = 0
    while rruns[ir] != ilmax: ir += 1
    if sum(lruns[:il]) != sum(rruns[:ir]):
        return False

    il = len(lruns) - 1
    while lruns[il] != ilmax: il -= 1
    ir = len(rruns) - 1
    while rruns[ir] != irmax: ir -= 1
    if sum(lruns[il:]) != sum(rruns[ir:]):
        return False

    if status:
        output = ''

    # now iterate through all words `u` with given prefix
    pref = W(prefix)
    for q in W.iterate_by_length(d-1-len(pref)):
        u = pref + q

        # compute the polytope of occurrences and check that it contains the
        # "middle occurrences"
        # TODO: count the number of points in extremal_mid_occurrences and add
        # it to status
        # similarly, if there is a bad descent word, add it
        left_occ = set(tuple(left[:i].count(a) for i in occ) + \
                       tuple(left[:i].count(b) for i in occ) \
                       for occ in occurrences(left, u))
        right_occ = set(tuple(right[:i].count(a) for i in occ) + \
                       tuple(right[:i].count(b) for i in occ) \
                       for occ in occurrences(right, u))

        inter = left_occ.intersection(right_occ)
        if not inter:
            if status:
                output += 'u = {}\n'.format(join(map(str,u)))
                output += 'no int. occ.'
            continue

        union = left_occ.difference(right_occ)
        if status:
            t0 = time()
        P = ppl_polytope(inter)
        if status:
            output += 'u = {}\n'.format(''.join(map(str,u)))
            output += 'num ext. occ.: {}\n'.format(len(inter))
            output += 'num int. occ.: {}\n'.format(len(union))
            output += 'num faces    : {}\n'.format(len(P.minimized_constraints()))
            output += 'num verts    : {}\n'.format(len(P.minimized_generators()))
            output += 'polytope computation in {}secs\n'.format(time()-t0)
            t0 = time()
        for o in union:
            pt = C_Polyhedron(point(Linear_Expression(o,0)))
            if not P.contains(pt):
                if status:
                    output += 'bad occurrence {}\n'.format(o, ''.join(left[i] for i in o))
                return (False,output) if status else False
        if status:
            output += 'containance test in {}secs\n'.format(time()-t0)
            output += '\n'
    return (True,output) if status else True

def is_vv_identity_parallel(left, right, d, W=None, prefix_length=None, ncpus=None, verbose=False, logfile=None):
    r"""
    Check identity using parallelization features

    INPT:

    - ``left``, ``right`` -- the identity to check

    - ``d`` -- the dimension of the space

    - ``prefix_length`` -- length of the prefix

    - ``ncpus`` -- (optional, default the number of cpus on the computer) number
      of cpus to use

    - ``verbose`` -- (optional, default ``False``) whether some additional
      information about the workers will be printed

    - ``logfile`` -- can be a string that specifies a filename for writing
      information about each occurrence polytope or sys.stdout
    """
    import multiprocessing as mp
    from misc import parallel_unfold

    close = False
    if isinstance(logfile, str):
        close = True
        logfile = open(logfile, 'w')
    if ncpus is None:
        ncpus = mp.cpu_count()
    pool = mp.Pool(ncpus)

    if prefix_length is None:
        from math import log
        prefix_length = int(round(log(ncpus) / log(2)))

    if W is None:
        alphabet = sorted(set(left).union(right))
        W = FiniteWords(alphabet)
    else:
        alphabet = W.alphabet()

    tasks =((verbose,is_vv_identity,left,right,d,W,prefix,True) \
            for prefix in W.iterate_by_length(prefix_length))
    for ans,status in pool.imap_unordered(parallel_unfold, tasks):
        if logfile is not None:
            logfile.write(status)
            logfile.flush()
        if ans is False:
            break
    pool.terminate()
    pool.join()
    if close:
        logfile.close()
    return ans

