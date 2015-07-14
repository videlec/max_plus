r"""
Generic semigroup tools
"""

import itertools
from sage.misc.misc_c import prod


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

def is_relation(r1, r2, elements):
    r"""
    Test if the relation ``r1 = r2`` is valid on ``pairs``.

    INPUT:

    - ``r1``, ``r2`` -- two lists of the non-negative integers `{0,1,2,...}`
      that encode a word in the free group

    - ``elements`` -- a list of elements that will be tested for the relation.
      If the relation involves number in `{0, ..., n-1}` then each element must
      be at least a `n`-tuple.

    EXAMPLES:

    Upper triangular matrices in dimension two commute::

        sage: m1 = matrix(2, [1,1,0,1])
        sage: m2 = matrix(2, [1,-1,0,1])
        sage: is_relation([0,1], [1,0], [(m1,m2)])
        True

    But not in dimension three::

        sage: m1 = matrix(3, [1,1,1,0,1,1,0,0,1])
        sage: m2 = matrix(3, [1,2,3,0,1,4,0,0,1])
        sage: m3 = matrix(3, [1,5,1,0,1,3,0,0,1])
        sage: is_relation([0,1], [1,0], [(m1,m2), (m1,m3), (m2,m3)])
        False

    But their commutator do::

        sage: r1 = [0,1,2,3,1,0,3,2]
        sage: r2 = [1,0,3,2,0,1,2,3]
        sage: elements = [(m1,m2,~m1,~m2), (m1,m3,~m1,~m3), (m2,m3,~m2,~m3)]
        sage: is_relation(r1, r2, elements)
        True
    """
    for p in elements:
        if prod(p[x] for x in r1) != prod(p[x] for x in r2):
            return False
    return True
