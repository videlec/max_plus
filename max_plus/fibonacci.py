r"""
Compute relations obtained from the Fibonacci word
"""

from itertools import product
from sage.combinat.words.all import words, FiniteWords

fib = words.FibonacciWord([0,1])
W = FiniteWords([0,1])

wp = fib[:10000]
wp1 = W([0,1]) * wp
wp2 = W([1,0]) * wp

def fibo_glue_factor(u1, u2):
    r"""
    Return the factor of minimal length of the Fibonacci word that contains
    ``u1`` as a prefix, ``u2`` as a suffix that do not overlap.

    EXAMPLES::

        sage: fibo_glue_factor(W([1]), W([1]))
        word: 101
        sage: fibo_glue_factor(W([0,1]), W([0]))
        word: 010
        sage: fibo_glue_factor(W([0,1,0,0]), W([0]))
        word: 010010
        sage: fibo_glue_factor(W([0,1,0]), W([0]))
        word: 0100

        sage: assert fibo_glue_factor(W([0]), W([0])) == W([0,0])
        sage: assert fibo_glue_factor(W([1]), W([1])) == W([1,0,1])
        sage: assert fibo_glue_factor(W([0,1]), W([0])) == W([0,1,0])
        sage: assert fibo_glue_factor(W([0,1,0,0]), W([0])) == W([0,1,0,0,1,0])
        sage: assert fibo_glue_factor(W([0,1,0]), W([0])) == W([0,1,0,0])
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
                prefixes.append(k + i +2)

    m = min(x.length() for x in ans)
    ans = [x for x in ans if x.length() == m]
    # unicity seems to be true
    if len(ans) != 1:
        raise RuntimeError("no unicity for u1={} and u2={}\ngot {}".format(u1,u2,ans))
    return ans[0]
    

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
        words[parent] = fibo_glue_factor(words[left], words[right])
        #print "({}, {}) -> {}".format(words[left], words[right], words[parent])
        
# Each list below correspond to a list of trees (whose vertices are arbitrarily
# labeled). A tree is encoded as a pair made of
#  * a list of triples (left, right, parent)
#  * the list of leaves

# list of trees with 2 leaves
T2 = [
    (((1,2,0),), (1,2))
]
# list of trees with 3 leaves
T3 = [
    (((3,4,2), (1,2,0)), (1,3,4)),
    (((3,4,1),(1,2,0)), (3,4,2))
]
# list of trees with 4 leaves
T4 = [
    (((3,4,1),(5,6,2),(1,2,0)), (3,4,5,6)),
    (((5,6,3),(3,4,1),(1,2,0)), (5,6,4,2)),
    (((5,6,4),(3,4,1),(1,2,0)), (3,5,6,2)),
    (((5,6,3),(3,4,2),(1,2,0)), (1,5,6,4)),
    (((5,6,4),(3,4,2),(1,2,0)), (1,3,5,6))
]
T = [None, None, T2, T3, T4]

def fibo_minimal_factors(m):
    r"""
    EXAMPLES::

        sage: fibo_minimal_factors(2)
        [word: 00, word: 01, word: 10, word: 101]
        sage: fibo_minimal_factors(3)
        [word: 001,
         word: 0010,
         word: 010,
         word: 0100,
         word: 0101,
         word: 100101,
         word: 100,
         word: 1010,
         word: 101,
         word: 101001]
    """
    ans = set()
    for labels in product([W([0]), W([1])], repeat=m):
        for tree,leaves in T[m]:
            words = [None] * (2*m-1)
            for i,j in zip(leaves,labels):
                words[i] = j
            tree_write(tree, words)
            ans.add(words[0])
    return sorted(ans)

def fibo_prefix(m):
    r"""
    Return the length of the smallest prefix of the Fibonacci word that contains
    all occurrences of subwords of length ``m`` for the product order.

    ONLY WORKS FOR m=2,3,4 FOR NOW

    EXAMPLES::

        sage: fibo_prefix(2)
        7
        sage: fibo_prefix(3)
        10
        sage: fibo_prefix(4)
        20
    """
    return max(wp.find(u) + u.length() for u in fibo_minimal_factors(m))

