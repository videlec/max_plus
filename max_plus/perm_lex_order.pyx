r"""
Generalization of the lexicographic order on tuples.
"""
from libc.stdlib cimport malloc, free, rand, RAND_MAX

cdef class PermLexOrder:
    r"""
    A permutation lexicographic order

    This class implements an order on tuples with a fixed permutation of indices
    and for each coordinate whether a lex or revlex is used.

    EXAMPLES::

        sage: from max_plus.perm_lex_order import PermLexOrder

        sage: p = PermLexOrder(n=5)
        sage: p
        Perm lexicographic order 0+ 1+ 2+ 3+ 4+
        sage: p = PermLexOrder([2,0,1,3], [1,0,0,1])
        sage: p
        Perm lexicographic order 2+ 0- 1- 3+
    """
    cdef int * _p
    cdef int * _o
    cdef int n

    def __init__(self, perm=None, min_o_max=None, n=None):
        r"""
        INPUT:

        - ``perm`` -- a permutation on `\{0, 1, ..., n\}`

        - ``min_o_max`` -- a list of ``n`` booleans

        - ``n`` -- an optional integer

        TESTS::

            sage: from max_plus.perm_lex_order import PermLexOrder
            sage: PermLexOrder([3,0,1], [1,0,1])  # indirect doctest
            Perm lexicographic order 3+ 0- 1+
            sage: PermLexOrder(n=4)
            Perm lexicographic order 0+ 1+ 2+ 3+
        """
        if perm is None:
            self.n = int(n)
            perm = range(n)
        else:
            self.n = len(perm)
        if min_o_max == None:
            min_o_max = [True] * self.n

        self._p = <int *> malloc(2*self.n*sizeof(int))
        self._o = self._p + self.n

        cdef int i
        for i in range(self.n):
            self._p[i] = perm[i]
            self._o[i] = 1 if min_o_max[i] else 0

    def __dealloc__(self):
        free(self._p)

    def perm(self):
        r"""
        TESTS::

            sage: from max_plus.perm_lex_order import PermLexOrder
            sage: PermLexOrder(n=4).perm()
            [0, 1, 2, 3]
            sage: PermLexOrder([3,0,1,2], [1,0,1,0]).perm()
            [3, 0, 1, 2]
        """
        return [self._p[i] for i in range(self.n)]

    def o(self):
        r"""
        TESTS::

            sage: from max_plus.perm_lex_order import PermLexOrder
            sage: PermLexOrder(n=4).o()
            [1, 1, 1, 1]
            sage: PermLexOrder([3,2,0,1], [1,0,1,0]).o()
            [1, 0, 1, 0]
        """
        return [self._o[i] for i in range(self.n)]

    def __repr__(self):
        r"""
        TESTS::

            sage: from max_plus.perm_lex_order import PermLexOrder
            sage: PermLexOrder(n=3)   # indirect doctest
            Perm lexicographic order 0+ 1+ 2+
        """
        return "Perm lexicographic order {}".format(
                  ' '.join('{}+'.format(i) if j else '{}-'.format(i) \
                           for (i,j) in zip(self.perm(), self.o())))

    def set_lex(self):
        r"""
        Set this order to the lexicographic one

        TESTS::

            sage: p = PermLexOrder([2,3,1,0], [0,0,0,0])
            sage: p
            Perm lexicographic order 2- 3- 1- 0-
            sage: p.set_lex()
            sage: p
            Perm lexicographic order 0- 1- 2- 3-
        """
        cdef int i
        for i in range(self.n):
            self._p[i] = i
            self._o[i] = 0

    def randomize(self):
        r"""
        Randomize this order.

        EXAMPLES::

            sage: from max_plus.perm_lex_order import PermLexOrder

            sage: p = PermLexOrder(n=5)
            sage: p
            Perm lexicographic order 0+ 1+ 2+ 3+ 4+
            sage: p.randomize()
            sage: p  # random
            Perm lexicographic order 3- 2+ 1- 0- 4+
        """
        cdef int i,j,k
        for i in range(self.n-1):
            j = i + (rand() % (self.n-1-i))
            k = self._p[i]
            self._p[i] = self._p[j]
            self._p[j] = k

            self._o[i] = rand() % 2

    def lt(self, tuple x, tuple y):
        r"""
        Check whether ``x`` is lesser than ``y``.

        EXAMPLES::

            sage: from max_plus.perm_lex_order import PermLexOrder
            sage: p = PermLexOrder([3,0,1,2], [1,1,1,1])
            sage: p.lt((0,1,3,2), (0,3,3,2))
            True
            sage: p.lt((1,1,1,0), (1,1,0,1))
            True
            sage: p.lt((2,2,2,2), (2,2,2,1))
            False

            sage: p = PermLexOrder([1,3,0,2], [1,0,0,1])
            sage: p.lt((2,0,2,2), (2,2,2,2))
            False
            sage: p.lt((2,2,2,0), (2,2,2,2))
            True
        """
        cdef int i,j
        for i in range(self.n):
            j = self._p[i]
            if x[j] < y[j]:
                return bool(self._o[j])
            elif x[j] > y[j]:
                return not self._o[j]

        return False

    def cmp(self, x, y):
        r"""
        EXAMPLES::

            sage: from max_plus.perm_lex_order import PermLexOrder
            sage: from itertools import product

            sage: l = list(product(range(2),repeat=3))

            sage: p = PermLexOrder([0,1,2], [1,1,0])
            sage: l.sort(cmp=p.cmp)
            sage: l
            [(0, 0, 1),
             (0, 0, 0),
             (0, 1, 1),
             (0, 1, 0),
             (1, 0, 1),
             (1, 0, 0),
             (1, 1, 1),
             (1, 1, 0)]
            sage: p.min_max(l)
            ((0, 0, 1), (1, 1, 0))

            sage: p = PermLexOrder([1,2,0], [0,1,0])
            sage: l.sort(cmp=p.cmp)
            sage: l
            [(1, 0, 1),
             (0, 0, 1),
             (1, 0, 0),
             (0, 0, 0),
             (1, 1, 1),
             (0, 1, 1),
             (1, 1, 0),
             (0, 1, 0)]
            sage: p.min_max(l)
            ((1, 0, 1), (0, 1, 0))

            sage: p = PermLexOrder([1,0,2],[1,0,1])
            sage: l.sort(cmp=p.cmp)
            sage: l
            [(0, 1, 0),
             (0, 1, 1),
             (1, 1, 0),
             (1, 1, 1),
             (0, 0, 0),
             (0, 0, 1),
             (1, 0, 0),
             (1, 0, 1)]
            sage: p.min_max(l)
            ((0, 1, 0), (1, 0, 1))

        TESTS::

            sage: V = FreeModule(ZZ,5)
            sage: p = PermLexOrder(n=5)
            sage: elts = [tuple(V.random_element()) for _ in range(50)]
            sage: for _ in range(100):
            ....:     p.randomize()
            ....:     p_min, p_max = p.min_max(elts)
            ....:     elts.sort(cmp=p.cmp)
            ....:     assert p_min == elts[0] and p_max == elts[-1]
        """
        return int(self.lt(y,x)) - int(self.lt(x,y))

    def min_max(self, l):
        r"""
        Return the minimum and the maximum of a list.

        EXAMPLES::

            sage: from max_plus.perm_lex_order import PermLexOrder

            sage: p01 = PermLexOrder([0,1], [1,1])
            sage: p10 = PermLexOrder([1,0], [1,1])
            sage: p01.min_max([(0,3),(1,1),(0,2),(1,2)])
            ((0, 2), (1, 2))
            sage: p10.min_max([(0,3),(1,1),(0,2),(1,2)])
            ((1, 1), (0, 3))

            sage: from itertools import product
            sage: p = PermLexOrder([1,0], [0,1])
            sage: p.min_max(product(range(3), repeat=2))
            ((2, 0), (0, 2))

            sage: p = PermLexOrder([2,1,0], [0,1,0])
            sage: p.min_max([(2,0,1),(1,2,3),(3,0,3),(4,2,0)])
            ((3, 0, 3), (4, 2, 0))
            sage: p.min_max([(1,3,3),(2,0,2),(2,2,2),(0,1,2)])
            ((1, 3, 3), (2, 2, 2))
        """
        it = iter(l)
        try:
            p_min = p_max = it.next()
        except StopIteration:
            raise ValueError
        for q in it:
            if self.lt(q, p_min):
                p_min = q
            if self.lt(p_max, q):
                p_max = q

        return p_min, p_max

