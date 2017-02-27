r"""
Compute relations obtained from the Fibonacci word
"""

from itertools import product
from sage.combinat.words.all import words, FiniteWords
from sage.misc.cachefunc import cached_function

fib = words.FibonacciWord([0,1])
W = FiniteWords([0,1])

wp = fib[:10000]
wp1 = W([0,1]) * wp
wp2 = W([1,0]) * wp

@cached_function
def glue_factors(u1, u2):
    r"""
    Return the factor of minimal length of the Fibonacci word that contains
    ``u1`` as a prefix, ``u2`` as a suffix that do not overlap.

    EXAMPLES::

        sage: from max_plus.fibonacci import glue_factors

        sage: W = Words([0,1])

        sage: glue_factors(W([1]), W([1]))
        (word: 101,)
        sage: glue_factors(W([0,1]), W([0]))
        (word: 010,)
        sage: glue_factors(W([0,1,0,0]), W([0]))
        (word: 010010,)
        sage: glue_factors(W([0,1,0]), W([0]))
        (word: 0100,)

    Sometimes the answer is not unique::

        sage: glue_factors(W([0,0,1,0]), W([0,1,0,0]))
        (word: 0010100100, word: 0010010100)

    TESTS::

        sage: assert glue_factors(W([0]), W([0])) == (W([0,0]),)
        sage: assert glue_factors(W([1]), W([1])) == (W([1,0,1]),)
        sage: assert glue_factors(W([0,1]), W([0])) == (W([0,1,0]),)
        sage: assert glue_factors(W([0,1,0,0]), W([0])) == (W([0,1,0,0,1,0]),)
        sage: assert glue_factors(W([0,1,0]), W([0])) == (W([0,1,0,0]),)

        sage: wp = words.FibonacciWord([0,1])[:1000]
        sage: F = set().union(*[wp.factor_set(k) for k in range(1,7)])
        sage: for u1 in F:
        ....:     for u2 in F:
        ....:         for w in glue_factors(u1, u2):
        ....:             assert w.has_prefix(u1)
        ....:             assert w.has_suffix(u2)
        ....:             assert w.length() >= u1.length() + u2.length()
    """
    i = wp.find(u2)     # automatic left extension of u2
    m2 = u2.length()
    if i == -1:
        raise RuntimeError("u2 not found in Fibonacci")
    prefixes = [i]
    ans = []
    while prefixes:
        #print prefixes
        i = prefixes.pop()
        assert wp[i:i+m2] == u2, i
        j = wp[:i].rfind(u1)

        if j != -1:
            #print "from i ={} found j = {}... {}".format(i,j,wp[j:i+m2])
            ans.append(wp[j:i+m2])
        else:
            if u1[0] == 0:
                i1 = 0; i2 = 1
            else:
                i1 = 1; i2 = 0

            assert wp1[i1] ==wp2[i2] == u1[0]
                
            if u1.length() < i+2-i1 and wp1[i1:].has_prefix(u1):
                ans.append(wp1[i1:i+2+m2])
            else:
                k = wp.find(wp1[:i+2+m2])
                if k == -1:
                    raise RuntimeError
                prefixes.append(k + i + 2)

            if u1.length() < i+2-i2 and wp2[i2:].has_prefix(u1):
                ans.append(wp2[i2:i+2+m2])
            else:
                k = wp.find(wp2[:i+2+m2])
                if k == -1:
                    raise RuntimeError
                prefixes.append(k+i+2)

    m = min(x.length() for x in ans)
    return tuple(x for x in ans if x.length() == m)
    

def tree_write(tree, words):
    r"""
    complete the tree
    """
    for left, right, parent in tree:
        assert 0 <= left < len(words)
        assert 0 <= right < len(words)
        assert 0 <= parent < len(words)
        assert words[left] is not None
        assert words[right] is not None
        assert words[parent] is None
        words[parent] = []
        for l in words[left]:
            for r in words[right]:
                words[parent].extend(glue_factors(l, r))

@cached_function
def trees_bis(m):
    r"""
    List of complete binary trees with ``m`` leaves (and hence 2m-1 vertices)
    """
    if m == 1:
        return ((),)
    else:
        return tuple((l,r) for i in range(1,m) for l in trees_bis(i) for r in trees_bis(m-i))

@cached_function
def num_leaves(t):
    if not t:
        return 1
    else:
        return num_leaves(t[0]) + num_leaves(t[1])


@cached_function
def trees(m):
    r"""
    List of trees with ``m`` leaves (and hence 2m-1 vertices)

    EXAMPLES::

        sage: from max_plus.fibonacci import trees
        sage: for t in trees(2): print t
        (((1, 2, 0),), (1, 2))
        sage: for t in trees(3): print t
        (((3, 4, 2), (1, 2, 0)), (1, 3, 4))
        (((2, 3, 1), (1, 4, 0)), (2, 3, 4))
    """
    if m == 1:
        return (((), (0,)),)
    else:
        ans = []
        lroot = 1
        for i in range(1,m):
            rroot = 2*i
            for l in _trees(i, root=lroot):
                for r in _trees(m-i, root=rroot):
                    ans.append((l[0]+r[0]+((lroot,rroot,0),), l[1]+r[1]))
        return tuple(ans)

def _trees(m, root):
    r"""
    Trees with shifted root
    """
    return tuple((tuple((root+i,root+j,root+k) for (i,j,k) in t[0]),
                  tuple(root+i for i in t[1])) for t in trees(m))

# list of trees with 2 leaves
#T2 = [
#    (((1,2,0),), (1,2))
#]
# list of trees with 3 leaves
#T3 = [
#    (((3,4,2), (1,2,0)), (1,3,4)),
#    (((3,4,1),(1,2,0)), (3,4,2))
#]
# list of trees with 4 leaves
#T4 = [
#    (((3,4,1),(5,6,2),(1,2,0)), (3,4,5,6)),
#    (((5,6,3),(3,4,1),(1,2,0)), (5,6,4,2)),
#    (((5,6,4),(3,4,1),(1,2,0)), (3,5,6,2)),
#    (((5,6,3),(3,4,2),(1,2,0)), (1,5,6,4)),
#    (((5,6,4),(3,4,2),(1,2,0)), (1,3,5,6))
#]
#T = [None, None, T2, T3, T4]

@cached_function
def minimal_lexico_factors(btree, leaves):
    r"""
    Return the list of minimal lexicographic occurrences of ``leaves`` along ``btree``.
    """
    assert num_leaves(btree) == len(leaves), (btree, len(leaves))

    if len(leaves) == 1:
        return leaves[0]

    left, right = btree
    nl = num_leaves(left)
    nr = num_leaves(right)
    assert nl + nr == len(leaves)
    ul = minimal_lexico_factors(left, leaves[:nl])
    ur = minimal_lexico_factors(right, leaves[nl:])
    ans = set()
    for l in ul:
        for r in ur:
            ans.update(glue_factors(l,r))
    return tuple(sorted(ans))

# this is method is very slow but we might be able to make it faster
# For example, if u is a factor of Fibonacci we know the answer
# if u is a concatenation of two factors as well...
@cached_function
def minimal_factors(u):
    r"""
    Return the set of factors of the Fibonacci word that contain
    all minimal occurrences of ``u`` (for the product order).

    EXAMPLES::

        sage: from max_plus.fibonacci import minimal_factors
        sage: W = FiniteWords([0,1])
        sage: minimal_factors(W([1,1,1]))
        [word: 100101, word: 101001]
        sage: minimal_factors(Word([0,0,0,0]))
        [word: 00100, word: 010010]
    """
    if len(u) == 1:
        return [u]
    elif len(u) == 2:
        return glue_factors(u[0:1], u[1:2])
    else:
        ans = set()
        for i in range(1,len(u)):
            l = minimal_factors(u[:i])
            r = minimal_factors(u[i:])
            for l in minimal_factors(u[:i]):
                for r in minimal_factors(u[i:]):
                    ans.update(glue_factors(l,r))
        return sorted(ans)

#    m = u.length()
#    labels = tuple((u[i:i+1],) for i in range(m))
#    for tree in trees_bis(m):
#        ans.update(minimal_lexico_factors(tree, labels))
#    return sorted(ans)

def sv_prefix_length(m, prefix=None):
    r"""
    Return the length of the smallest prefix of the Fibonacci word that contains
    all occurrences of subwords of length ``m`` for the product order.

    EXAMPLES::

        sage: from max_plus.fibonacci import sv_prefix_length
        sage: sv_prefix_length(1)
        2
        sage: sv_prefix_length(2)
        7
        sage: sv_prefix_length(3)
        10
        sage: sv_prefix_length(4)
        20
        sage: sv_prefix_length(5)
        23
        sage: sv_prefix_length(6)
        28
        sage: sv_prefix_length(7)
        31
        sage: sv_prefix_length(8)
        54
        sage: sv_prefix_length(9)    # not tested # too long
        57
        sage: sv_prefix_length(10)   # not tested # too long
        62
        sage: sv_prefix_length(11)   # not tested # too long
        65
        sage: sv_prefix_length(12)   # not tested # too long
        75
        sage: sv_prefix_length(13)   # not tested # too long
        78
        sage: sv_prefix_length(14)   # not tested # too long
        91
        sage: sv_prefix_length(15)   # not tested # too long
        96
        sage: sv_prefix_length(16)   # not tested # too long
        143
        sage: sv_prefix_length(17)   # not tested # too long
        146
        sage: sv_prefix_length(18)   # not tested # too long
        151

        sage: l = [2, 7, 10, 20, 23, 28, 31, 54, 57, 62, 65, 75, 78, 91, 96, 143, 146, 151]
        sage: [l[i+1]-l[i] for i in range(len(l)-1)]
        [5, 3, 10, 3, 5, 3, 23, 3, 5, 3, 10, 3, 13, 5, 47, 3, 5]
    """
    if prefix is None:
        prefix = W()
    ans = set()
    for u in W.iterate_by_length(m - len(prefix)):
        u = prefix + u
        print u
        ans.update(minimal_factors(u))
    return max(wp.find(u) + u.length() for u in ans)

def sv_prefix_length_parallel(m, ncpus=None, prefix_length=None, verbose=False):
    r"""
    The main problem with this naive parallel version is that the cache
    are not shared and many computations are done several times.

    (moreover once the pool is killed computations should be lost)

    However, it seems that the fork make a copy of the current cache...
    """
    import multiprocessing as mp
    from misc import parallel_unfold

    if ncpus is None:
        ncpus = mp.cpu_count()
    pool = mp.Pool(ncpus)

    if prefix_length is None:
        from math import log
        prefix_length = min(1 + int(log(ncpus) / log(2)), m)

    tasks = ((verbose, sv_prefix_length, m, prefix) for prefix in W.iterate_by_length(prefix_length))
    final_ans = 0
    for ans in pool.imap_unordered(parallel_unfold, tasks):
        final_ans = max(ans, final_ans)
    pool.close()
    pool.join()
    return final_ans
