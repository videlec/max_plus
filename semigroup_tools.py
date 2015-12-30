r"""
Generic semigroup tools
"""

import itertools
from sage.misc.misc_c import prod

indent = 0

def occurrences(w, u):
    r"""
    Construct the polytope of all occurrences of ``u`` in ``w`` as a subword.
    """
#    global indent
#    print " " * indent + "NEW CALL"
#    indent += 2
    if len(u) == 0:
        yield ()
    else:
        for pos,letter in enumerate(w):
            if letter == u[0]:
#                print " " * indent + "found letter {} at {}".format(u[0], pos)
                for i in occurrences(w[pos+1:], u[1:]):
                    yield (pos,) + tuple(j+pos+1 for j in i)
#    indent -= 2
#    print " " * indent + "CLOSE"

def products_n(a, b, n, one):
    r"""
    Run through all products of the elements ``a`` and ``b`` of given length.

    EXAMPLES::

        sage: m1 = matrix([[1,1],[0,1]])
        sage: m2 = matrix([[1,0],[1,1]])
        sage: one = matrix([[1,0],[0,1]])
        sage: g = products_n(m1, m2, 3, one)
        sage: g.next()
        (
                   [1 3]
        [0, 0, 0], [0 1]
        )
        sage: g.next()
        (
                   [3 2]
        [0, 0, 1], [1 1]
        )
        sage: g.next()
        (
                   [2 3]
        [0, 1, 0], [1 2]
        )
        sage: g.next()
        (
                   [3 1]
        [0, 1, 1], [2 1]
        )
        sage: g.next()
        (
                   [1 2]
        [1, 0, 0], [1 3]
        )
        sage: g.next()
        (
                   [2 1]
        [1, 0, 1], [3 2]
        )
        sage: g.next()
        (
                   [1 1]
        [1, 1, 0], [2 3]
        )
        sage: g.next()
        (
                   [1 0]
        [1, 1, 1], [3 1]
        )
    """
    i = []
    branch = [one]

    while True:
        while len(i) < n:
            i.append(0)
            branch.append(branch[-1] * a)
        yield i, branch[-1]

        while i and i[-1] == 1:
            i.pop()
            branch.pop()
        if not i:
            return
        i[-1] = 1
        branch[-1] = branch[-2] * b

def products_p(a, b, na, nb, one):
    r"""
    Run through all products of the elements ``a`` and ``b`` containing exactly
    ``na`` times ``a`` and ``nb`` times ``b``.

    At each step, the iterator returns a pair that consists in a list of `0` and
    `1` which is the "code" of the product (`0` stands for `a` and `1` stands
    for `b`). The second item is the element of the semigroup.

    INPUT:

    - ``a``, ``b`` -- two elements on which we can make products

    - ``na``, ``nb`` -- the number of ``a`` and the number of ``b``

    - ``one`` -- the identity in the semigroup

    EXAMPLES::

        sage: m1 = matrix(2,2, [[1,1],[0,1]])
        sage: m2 = matrix(2,2, [[1,0],[1,1]])
        sage: for p in products(m1,m2,2,2,identity_matrix(2)):
        ....:     print p[0], p[1].list()
        [0, 0, 1, 1] [5, 2, 2, 1]
        [0, 1, 0, 1] [5, 3, 3, 2]
        [0, 1, 1, 0] [3, 4, 2, 3]
        [1, 0, 0, 1] [3, 2, 4, 3]
        [1, 0, 1, 0] [2, 3, 3, 5]
        [1, 1, 0, 0] [1, 2, 2, 5]
    """
    i = []
    branch = [one]
    na = int(na)
    nb = int(nb)
    n = na+nb

    while True:
        for _ in range(na):
            i.append(0)
            branch.append(branch[-1] * a)
        for _ in range(nb):
            i.append(1)
            branch.append(branch[-1] * a)
        na = 0
        nb = 0
        yield i, branch[-1]

        while i and (i[-1] == 1 or nb == 0):
            if i.pop() == 0:
                na += 1
            else:
                nb += 1
            branch.pop()
        if not i:
            return
        na += 1
        nb -= 1
        i[-1] = 1
        branch[-1] = branch[-2] * b

def products_all_subwords(a, b, n, m, one):
    r"""
    An iterator through all words on a,b of length n that contains all words of
    length m.
    """
    branch = [one]
    raise NotImplementedError


def minimal_all_subwords_prefix(n, k):
    r"""
    Iterator through all words of length up to that contains all subwords and
    such that it is not the case for any of its prefix

    k=2

    00111110
    00101
    """
    if k == 0:
        yield []
    else:
        for w in minimal_all_subwords_prefix(n-2,k-1):
            for m in xrange(1,n-len(w)):
                yield w + [0]*m + [1]
                yield w + [1]*m + [0]

def minimal_all_subwords_suffix(length=None, max_length=None, k=None):
    r"""
    Iterator through all words of length up to that contains all subwords and
    such that it is not the case for any of its prefix

    k=2

    00111110
    00101
    """
    if k == 0:
        yield []
    else:
        for w in minimal_all_subwords_suffix(max_length=n-2, k=k-1):
            if length is not None:
                yield [0] + [1]*(n-1-len(w)) + w
                yield [1] + [0]*(n-1-len(w)) + w
            else:
                for m in range(1, max_length-len(w)-1):
                    yield [0] + [1]*m + w
                    yield [1] + [0]*m + w

def constraints_from_subwords(i, d):
    r"""
    Given a word ``i`` and a dimension ``d`` return the word ``j`` with holes (i.e.
    a word on ``0``, ``1``, ``None``) so that

    - any first/last occurrence of `(d-1)`-subword in ``i`` are also in ``j`` at
      the same position

    EXAMPLES::

        sage: w = (0,0,1,0,1,0,0,1,0,0,1)
        sage: p,s = constraints_from_subwords(w, 3)
        sage: w[:p]
        (0, 0, 1, 0, 1)
        sage: w[s:]
        (1, 0, 0, 1)
        sage: constraints_from_subwords((0,1,0,0), 3)
        (-1, -1)
    """
    if len(i) < d-1:
        return (-1,-1)

    # subwords
    k = d-1
    zero = False
    one  = False
    kk   = 0
    ppos = 0

    while ppos < len(i) and kk < k:
        if i[ppos] == 0:
            zero = True
        else:
            one = True
        if zero and one:
            zero = one = False
            kk += 1
        ppos += 1

    spos = len(i)-1
    kk = 0
    zero = False
    one  = False

    while spos >= 0 and kk < k:
        if i[spos] == 0:
            zero = True
        else:
            one = True
        if zero and one:
            zero = one = False
            kk += 1
        spos -= 1
    spos += 1

    if spos < ppos:
        return -1,-1
    else:
        return ppos,spos

    # factors
    #w = ''.join(str(letter) for letter in i)
    #F = set(w[pos:pos+d-1] for pos in range(len(i)-d+1))
    #for f in F:
    #    pos = w.find(f)
    #    ii[pos:pos+d-1] = i[pos:pos+d-1]
    #    pos = w.rfind(f)
    #    ii[pos:pos+d-1] = i[pos:pos+d-1]

