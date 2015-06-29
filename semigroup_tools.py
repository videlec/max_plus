r"""
Generic semigroup tools
"""
from sage.misc.misc_c import prod

def products(a, b, na, nb, one):
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

def is_relation(i1, i2, elements):
    r"""
    Test if the relation ``i1 = i2`` is valid on ``elements``.

    INPUT:

    - ``i1``, ``i2`` -- two lists of `0` and `1` that encode a word in the free
      group on two generators

    - ``elements`` -- a set of element of a common semigroup

    EXAMPLES::

        sage: m1 = matrix(2, [1,1,0,1])
        sage: m2 = matrix(2, [1,-1,0,1])
        sage: is_relation([0,1], [1,0], [m1,m2])
        True
    """
    for a in elements:
        for b in elements:
            if a is b:
                continue
            p1 = prod(a if x == 0 else b for x in i1)
            p2 = prod(a if x == 0 else b for x in i2)
            if p1 != p2:
                return False
    return True


