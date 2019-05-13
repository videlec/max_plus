r"""
Identities for matrices in `B^{sv}_d`.
"""

from __future__ import print_function, division, absolute_import

from time import time
import itertools
import math
import sys


from .sage_import import *
from .combinat import (extremal_occurrences,
                      prefix_suffix_all_subwords,
                      iterate_over_holes,
                      has_all_subwords)
from .perm_lex_order import PermLexOrder
from .convex_hull import ppl_polytope

##############################
# Global immutable constants #
##############################

t01 = (0,1)
W01 = FiniteWords(t01)

def fill_sv_with_random_lex_samples(u, v, m, W=None, n_max=5):
    r"""
    Try to fill position of ``v`` with random occurrences of ``m`` in ``u``

    EXAMPLES::

        sage: from max_plus.sv_identities import fill_sv_with_random_lex_samples

        sage: u = (0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 1)
        sage: v = [None] * len(u)
        sage: m = (0, 0)
        sage: fill_sv_with_random_lex_samples(u, v, [0,0])
        sage: v  # random
        [0, None, None, 0, 0, None, None, None, 0, 0, None]
        sage: fill_sv_with_random_lex_samples(u, v, [0,1])
        sage: v  # random
        [0, 1, None, 0, 0, None, None, None, 0, 0, 1]
        sage: fill_sv_with_random_lex_samples(u, v, [1,0])
        sage: v  # random
        [0, 1, 1, 0, 0, None, None, 1, 0, 0, 1]
        sage: fill_sv_with_random_lex_samples(u, v, [1,1])
        sage: v  # random
        [0, 1, 1, 0, 0, None, 1, 1, 0, 0, 1]
        sage: v[5] is None
        True

    TESTS::

        sage: for u in ((0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 0),
        ....:           (1, 0, 0, 1, 1, 0, 1, 0, 1, 1, 0),
        ....:           (1, 1, 0, 0, 1, 0, 1, 0, 1, 1, 0),
        ....:           (0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 0),
        ....:           (1, 0, 0, 1, 1, 1, 1, 0, 1, 1, 0),
        ....:           (1, 1, 0, 0, 1, 1, 1, 0, 1, 1, 0)):
        ....:     for _ in range(100):
        ....:         v = [None] * 11
        ....:         fill_sv_with_random_lex_samples(u, v, [0,0], n_max=50)
        ....:         fill_sv_with_random_lex_samples(u, v, [0,1], n_max=50)
        ....:         fill_sv_with_random_lex_samples(u, v, [1,0], n_max=50)
        ....:         fill_sv_with_random_lex_samples(u, v, [1,1], n_max=50)
        ....:         assert v[5] is None, "u = {}".format(u)
    """
    k = len(m)
    if W is None:
        alphabet = sorted(set(u))
        W = FiniteWords(alphabet)
    else:
        alphabet = W.alphabet()

    # this is a list of tuple of positions
    # we need to map (i0, i1, i2, i3) to
    # (i0, i1-i0, i2-i1, i3-i2, -i3)
    occ = extremal_occurrences(u, m, W)
    if not occ:
        raise ValueError("u = {}, v = {}, m = {}".format(u, v, m))
    for i in range(len(occ)):
        o = occ[i]
        assert len(o) == k, "o = {} m = {}".format(o, m)
        occ[i] = (o[0],) + tuple([o[j+1]-o[j] for j in range(k-1)]) + o
        assert len(occ[i]) == 2*k

    p = PermLexOrder(n=k)
    n = 0
    while n < n_max:
        o_min_max = p.min_max(occ)
        for pos in o_min_max:
            for i in pos[k:]:
                if v[i] is None:
                    n = 0
                    v[i] = u[i]
        n += 1
        p.randomize()

def is_sv_identity(left, right, d, W=None, prefix=(), check_common_factors=True, status=False):
    r"""
    Check if ``(left, right)`` is a B^{sv} identity in dimension ``d``.

    This method go through all subwords of length ``d-1`` and for each of them
    see whether some polytopes coincide.
    THIS IS WRONG: we should go through all words of length up to d-1.

    INPUT:

    - ``left``, ``right`` -- the identity to test

    - ``d`` -- dimension

    - ``prefix`` -- an optional prefix (mostly used for parallelization, see the
      function ``is_sv_identity_parallel`` below).

    - ``check_common_factors`` -- (default is ``False``) whether to skip the
      computation of the polytope for common factors of ``p`` and ``s``

    - ``status`` -- if ``True``, then instead of returning a boolean, returns a
      pair ``(boolean, status_string)`` where ``status_string`` gives some
      details about the computation.

    OUTPUT:

    Either a boolean or a pair ``(boolean, status_string)`` if ``status=True``.

    EXAMPLES::

        sage: from max_plus import *

        sage: p = 'xxyyx'
        sage: s = 'xxyxy'
        sage: u = p + 'x' + s
        sage: v = p + 'y' + s
        sage: is_sv_identity(u, v, 3)
        True
        sage: is_sv_identity(u, v, 4)
        False
        sage: ans, info = is_sv_identity(u, v, 3, status=True)
        sage: print(info)   # only one factor tested!
        u = yy
        num ext. occ.: 3
        num int. occ.: 3
        num faces    : 3
        num verts    : 3
        polytope computation in ...secs
        containance test in ...secs

        sage: ans, info = is_sv_identity(u, v, 3, check_common_factors=False, status=True)
        sage: print(info)   # all factors are tested
        u = xx
        num ext. occ.: 5
        num int. occ.: 1
        num faces    : 4
        num verts    : 4
        ...
        u = xy
        num ext. occ.: 4
        num int. occ.: 0
        num faces    : 4
        num verts    : 4
        ...
        u = yx
        num ext. occ.: 4
        num int. occ.: 0
        num faces    : 4
        num verts    : 4
        ...
        u = yy
        num ext. occ.: 3
        num int. occ.: 3
        num faces    : 3
        num verts    : 3
        ...

        sage: p,s = vincent_sv_prefix_suffix(4)
        sage: u = p+'x'+s
        sage: v = p+'y'+s
        sage: is_sv_identity(u, v, 4) and is_sv_identity(u, v, 4, check_common_factors=False)
        True
        sage: is_sv_identity(u, v, 5) or is_sv_identity(u, v, 5, check_common_factors=False)
        False

        sage: p,s = vincent_sv_prefix_suffix(5)
        sage: u = p+'x'+s
        sage: v = p+'y'+s
        sage: is_sv_identity(u, v, 5)
        True
        sage: is_sv_identity(u, v, 6)
        False

        sage: x,y = symbolic_max_plus_matrices_band(5, 2, 's', 'v')
        sage: p = x*y*y*x*x*x*y*y*y*y*x*x*x*x  # long time - 1 sec
        sage: s = y*y*y*y*x*x*x*x*y*y*y*x*x*y  # long time - 1 sec
        sage: p*x*s == p*y*s                   # long time - 6 secs
        True

        sage: p = 'xxxyyxyxx'
        sage: s = 'xxxyxxyxy'
        sage: u = p + 'x' + s
        sage: v = p + 'y' + s
        sage: is_sv_identity(u, v, 4)
        True
        sage: is_sv_identity(u, v, 5)
        False

        sage: for u in [(0, 1, 1, 0, 0, 0, 0, 1, 0, 1, 0),
        ....:           (1, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0),
        ....:           (1, 1, 0, 0, 1, 0, 0, 0, 1, 0, 1),
        ....:           (1, 1, 0, 0, 1, 0, 0, 1, 1, 0, 1)]:
        ....:     v = u[:5] + (1,) + u[-5:]
        ....:     test1 = is_sv_identity(u,v,3,check_common_factors=True)
        ....:     test2 = is_sv_identity(u,v,3,check_common_factors=False)
        ....:     if test1 != test2:
        ....:         raise RuntimeError("u = {}, test1 = {}, test2 = {}".format(
        ....:                        u,test1,test2))

        sage: u = (0, 1, 1, 0, 0, 0, 0, 1, 0, 1, 0)
        sage: v = u[:5] + (1,) + u[-5:]
        sage: ans,info = is_sv_identity(u, v, 3, status=True)
        sage: u = (0, 1, 1, 0, 0, 0, 0, 0, 1, 1, 0)
        sage: v = u[:5] + (1,) + u[-5:]
        sage: ans,info = is_sv_identity(u, v, 3, status=True)

        sage: for u,v in [((0,0),(0,0)), ((0,),(0,0)), ((0,0,0),(0,1,0)), ((0,),(1,1)),
        ....:     ((0,1,1,0),(1,0,0,1)), ((0,1,0,1,0),(0,1,1,1,0)),
        ....:     ((0,1,1,1,0,0,1,0,1,0,1,0,1,0,1,0),(0,1,1,1,0,0,1,0,1,1,0,0,1,0,1,0))]:
        ....:     ans,info = is_sv_identity(u,v,3,status=True)
        ....:     assert isinstance(ans, bool) and isinstance(info,str), 'u = {}  v = {}'.format(u,v)
    """
    if W is None:
        alphabet = sorted(set().union(left,right))
        W = FiniteWords(alphabet)
    else:
        alphabet = W.alphabet()

    left = W(left)
    right = W(right)

    if len(alphabet) <= 1:
        if len(left) == len(right):
            return (True, 'equal') if status else True
        else:
            return (False, 'one letter different lengths') if status else False
    if len(alphabet) != 2:
        raise ValueError("must be on two letters, got {}".format(alphabet))

    if left == right:
        return (True, 'equal') if status else True

    if len(left) != len(right) or \
       not has_all_subwords(left, d-1) or \
       not has_all_subwords(right, d-1):
       return (False, 'subwords') if status else False

    n = len(left)

    # compute common prefix and suffix
    i = 0
    while left[i] == right[i]:
        i += 1
    p = left[:i]
    i = n-1
    while left[i] == right[i]:
        i -= 1
    s = left[i+1:]
    assert left.has_prefix(p) and right.has_prefix(p)
    assert left.has_suffix(s) and right.has_suffix(s)
    if not has_all_subwords(p, d-1) or not has_all_subwords(s, d-1):
        return (False, 'subwords') if status else False

    if left.factor_set(d-1) != right.factor_set(d-1):
        return (False, 'factor set') if status else False

    facts = p.factor_set(d-1).intersection(s.factor_set(d-1))

    if status:
        output = ''

    # now iterate through all words `u` with given prefix
    pref = W(prefix)
    for q in W.iterate_by_length(d-1-len(pref)):
        u = pref + q
        if check_common_factors and u in facts:
            # ignore `u` that have a factor occurrence in both p and s
            continue

        # compute the polytope of occurrences and check that it contains the
        # "middle occurrences"
        # TODO: count the number of points in extremal_mid_occurrences and add
        # it to status
        # similarly, if there is a bad descent word, add it
        left_occ = set(extremal_occurrences(left, u))
        right_occ = set(extremal_occurrences(right, u))

        inter = left_occ.intersection(right_occ)
        if not inter:
            raise RuntimeError("THIS SHOULD NOT HAPPEN!")

        sym_diff = left_occ.symmetric_difference(right_occ)
        t0 = time()
        P = ppl_polytope(inter)
        t0 = time() - t0
        if status:
            output += 'u = {}\n'.format(''.join(map(str,u)))
            output += 'num ext. occ.: {}\n'.format(len(inter))
            output += 'num int. occ.: {}\n'.format(len(sym_diff))
            output += 'num faces    : {}\n'.format(len(P.minimized_constraints()))
            output += 'num verts    : {}\n'.format(len(P.minimized_generators()))
            output += 'polytope computation in {}secs\n'.format(t0)
        t0 = time()
        for o in sym_diff:
            pt = C_Polyhedron(point(Linear_Expression(o,0)))
            if not P.contains(pt):
                t0 = time() - t0
                if status:
                    output += 'bad occurrence {} in {}secs\n'.format(o,
                            ''.join(str(left[i]) for i in o), t0)
                return (False,output) if status else False
        t0 = time() - t0
        if status:
            output += 'containance test in {}secs\n'.format(t0)
            output += '\n'
    return (True,output) if status else True

def fill_sv(u, d, alphabet=None, verbose=False):
    r"""
    Iterator through the candidates compatible with ``u`` for a (s,v)-relation
    in `B^{sv}_d`

    Return a pair ``(v, is_full)``.

    EXAMPLES::

        sage: from max_plus.sv_identities import fill_sv

        sage: u = (0,1,1,0,0,0,1,1,0,0,1)
        sage: v, is_trivial = fill_sv(u, 3, (0,1))
        sage: v    # random
        [0, 1, 1, 0, None, None, None, 1, 0, 0, 1]
        sage: is_trivial
        False
        sage: v[5] is None
        True

        sage: u =[0, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 0, 0, 1]
        sage: v, is_trivial = fill_sv(u, 4, (0,1))
        sage: v  # random
        [0, 1, 1, 0, 0, 0, 1, 1, 1, None, 0, 0, 0, 1, 1, 1, 0, 0, 1]
        sage: is_trivial
        False
        sage: v[9] is None
        True

    TESTS::

        sage: for u in ((0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 0),
        ....:           (0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 0),
        ....:           (1, 1, 0, 0, 1, 0, 1, 0, 1, 1, 0),
        ....:           (1, 1, 0, 0, 1, 1, 1, 0, 1, 1, 0),
        ....:           (1, 0, 0, 1, 1, 0, 1, 0, 1, 1, 0),
        ....:           (1, 0, 0, 1, 1, 1, 1, 0, 1, 1, 0)):
        ....:     assert fill_sv(u, 3, (0,1), verbose=True)[0][5] is None
    """
    n = len(u)
    if alphabet is None:
        alphabet = set(u)
    is_trivial,v = prefix_suffix_all_subwords(u, d-1)
    if is_trivial:
        if verbose: print("trivial")
        return v,True
    for m in itertools.product(alphabet, repeat=d-1):
        fill_sv_with_random_lex_samples(u, v, m, W01)
        if all(k is not None for k in v):
            return v,True

    return v,False

def is_sv_identity_parallel(left, right, d, W=None, prefix_length=None, ncpus=None, verbose=False,
        check_common_factors=True, logfile=None):
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

    - ``check_common_factors`` -- (optional, default ``True``) whether to also
      test factors that are common to both ``p`` and ``s``

    - ``logfile`` -- can be a string that specifies a filename for writing
      information about each occurrence polytope or sys.stdout

    EXAMPLES::

        sage: from max_plus import *

        sage: p,s = vincent_sv_prefix_suffix(5)
        sage: u = p+'x'+s
        sage: v = p+'y'+s
        sage: is_sv_identity_parallel(u, v, 5, 3)  # not tested (fail in doctest)
        True

        sage: p,s = vincent_sv_prefix_suffix(6)
        sage: u = p+'x'+s
        sage: v = p+'y'+s
        sage: is_sv_identity_parallel(u, v, 6) # not tested
        True

    The following will write its output in a file named 'logfile.txt'::

        sage: is_sv_identity_parallel(u, v, 6, prefix_length=4, logfile='logfile.txt') # not tested
        True

    And for live information you can either turn on the ``verbose`` option (for
    a per job information) or set the ``logfile`` option to ``sys.stdout`` (for
    a per polytope information)::

        sage: is_sv_identity_parallel(u, v, 6, prefix_length=4, verbose=True)  # not tested
        True
        sage: import sys
        sage: is_sv_identity_parallel(u, v, 6, prefix_length=4, logfile=sys.stdout)  # not tested
        True
    """
    import multiprocessing as mp
    from .misc import parallel_unfold

    close = False
    if isinstance(logfile, str):
        close = True
        logfile = open(logfile, 'w')
    if ncpus is None:
        ncpus = mp.cpu_count()
    pool = mp.Pool(ncpus)

    if W is None:
        alphabet = sorted(set(left).union(right))
        W = FiniteWords(alphabet)
    else:
        alphabet = W.alphabet()

    left = W(left)
    right = W(right)

    if len(left) != len(right):
        return False

    if prefix_length is None:
        from math import log
        prefix_length = min(1 + int(log(ncpus) / log(2)), d-1)

    tasks = ((verbose,is_sv_identity,left,right,d,W,prefix,check_common_factors,True) \
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

def vincent_sv_prefix_suffix(d):
    r"""
    Return the ``p``,``s`` from the conjecture

    EXAMPLES::

        sage: from max_plus import *

        sage: vincent_sv_prefix_suffix(5)
        ('xyyxxxyyyyxxxx', 'yyyyxxxxyyyxxy')
    """
    p = ''
    s = ''
    pletter = 'x'
    sletter = 'y'
    for i in range(1,d):
        p = p + pletter * i
        s = sletter * i + s
        pletter,sletter = sletter,pletter
    p = p + pletter * (d-1)
    s = sletter * (d-1) + s
    return p,s

def sv_candidates(n, d, u_start=None, u_stop=None, nb_mats=10000):
    r"""
    Iterator through the candidates for identities.

    There is some randomness involved in the generation. Two runs of this
    function might be different!

    EXAMPLES::

        sage: from max_plus.sv_identities import (
        ....:     sv_candidates,
        ....:     is_sv_identity)

        sage: def print_identity(i):
        ....:     print(''.join(map(str,i[0])) + ' ' + ''.join(map(str,i[1])))

        sage: for i in sv_candidates(5, 2):
        ....:     if is_sv_identity(i[0], i[1], 2):
        ....:         print_identity(i)
        01110 01010
        10101 10001
        10110 10010

        sage: for i in sv_candidates(11, 3):
        ....:     if is_sv_identity(i[0], i[1], 3):
        ....:         print_identity(i)
        00110101100 00110001100
        01011111010 01011011010
        01100100110 01100000110
        01100101100 01100001100
        10011100110 10011000110
        10011101100 10011001100
        10011111001 10011011001
        10011111010 10011011010
        10100100101 10100000101
        10100100110 10100000110
        10100101100 10100001100
        11001100110 11001000110
        11001101100 11001001100
        11001110011 11001010011
        11001111001 11001011001
        11001111010 11001011010

        sage: for i in sv_candidates(11, 3,
        ....:     u_start = (1,)*2 + (0,)*9,
        ....:     u_stop  = (1,)*2 + (0,0) + (1,)*7):
        ....:     if is_sv_identity(i[0], i[1], 3):
        ....:         print_identity(i)
        11001100110 11001000110
        11001101100 11001001100
        11001110011 11001010011
        11001111001 11001011001
        11001111010 11001011010

        sage: l = [[i for i in sv_candidates(12,3) \
        ....:      if is_sv_identity(i[0], i[1], 3, check_common_factors=True)] \
        ....:      for _ in range(5)]                                               # long time
        sage: l.extend([[i for i in sv_candidates(12,3) \
        ....:      if is_sv_identity(i[0], i[1], 3, check_common_factors=False)] \
        ....:      for _ in range(5)])                                              # long time
        sage: for i in range(len(l)):                                               # long time
        ....:     for j in range(i):
        ....:         if l[i] != l[j]:
        ....:             print(set(l[i]).symmetric_difference(l[j]))
    """
    from .max_plus_int import (random_integer_max_plus_matrices_band,
            filter_upper_relation)

    from .word import product_start_stop

    a = int(math.sqrt(sys.maxint)/2)
    elements = [random_integer_max_plus_matrices_band(d, -a, a, ord('s'), ord('v')) for _ in range(nb_mats)]

    if u_start is None:
        u_start = (0,)*n
    if u_stop is None:
        u_stop = (1,)*n

    for u in product_start_stop(u_start, u_stop):
        if u[::-1] > u:
            continue
        if not has_all_subwords(u, d-1):
            continue
        v, no_hole = fill_sv(u, d, alphabet=t01)
        if no_hole:
            continue
        holes = [i for i in range(len(v)) if v[i] is None]
        tu = tuple(u)
        for i in filter_upper_relation(
                   ((tu,v) for v in iterate_over_holes(u, v, holes, t01) if has_all_subwords(v, d-1)),
                   n, d, elements):
            yield i

def sv_identities_iterator(n, d, u_start=None, u_stop=None, nb_mats=1000):
    r"""
    Iterator through the identities in `B^{sv}_d` of length ``n`` (up to
    symmetry).

    EXAMPLES::

        sage: from max_plus import *

        sage: for i in sv_identities_iterator(6, 2):
        ....:     print(i)
        ((0, 1, 1, 0, 1, 0), (0, 1, 0, 0, 1, 0))
        ((0, 1, 1, 0, 1, 0), (0, 1, 0, 1, 1, 0))
        ((0, 1, 1, 1, 0, 0), (0, 1, 0, 1, 0, 0))
        ((0, 1, 1, 1, 1, 0), (0, 1, 0, 0, 1, 0))
        ((0, 1, 1, 1, 1, 0), (0, 1, 0, 1, 1, 0))
        ((0, 1, 1, 1, 1, 0), (0, 1, 1, 0, 1, 0))
        ((1, 0, 0, 1, 1, 0), (1, 0, 0, 0, 1, 0))
        ((1, 0, 1, 0, 0, 1), (1, 0, 0, 0, 0, 1))
        ((1, 0, 1, 0, 0, 1), (1, 0, 0, 1, 0, 1))
        ((1, 0, 1, 0, 1, 0), (1, 0, 0, 0, 1, 0))
        ((1, 0, 1, 0, 1, 0), (1, 0, 0, 1, 1, 0))
        ((1, 0, 1, 1, 0, 0), (1, 0, 0, 1, 0, 0))
        ((1, 0, 1, 1, 0, 1), (1, 0, 0, 0, 0, 1))
        ((1, 0, 1, 1, 0, 1), (1, 0, 0, 1, 0, 1))
        ((1, 0, 1, 1, 0, 1), (1, 0, 1, 0, 0, 1))
        ((1, 0, 1, 1, 1, 0), (1, 0, 0, 0, 1, 0))
        ((1, 0, 1, 1, 1, 0), (1, 0, 0, 1, 1, 0))
        ((1, 0, 1, 1, 1, 0), (1, 0, 1, 0, 1, 0))
        ((1, 1, 0, 1, 0, 1), (1, 1, 0, 0, 0, 1))
        ((1, 1, 0, 1, 1, 0), (1, 1, 0, 0, 1, 0))
    """
    for i in sv_candidates(n, d, u_start=u_start, u_stop=u_stop):
        if is_sv_identity(i[0], i[1], d):
            yield i
