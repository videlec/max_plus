r"""
Integer max plus matrices
"""
from libc.stdlib cimport malloc, calloc, free, rand, RAND_MAX
# note: on my computer RAND_MAX is 2**31 - 1
from libc.limits cimport LONG_MIN, LONG_MAX
from libc.math cimport round
from libc.string cimport memcpy

from cpython.object cimport Py_EQ, Py_NE

from sage.ext.memory_allocator cimport MemoryAllocator

def minus_infinity():
    return <long>LONG_MIN

cdef inline long randlong(long min_coeff, long max_coeff):
    return min_coeff + ((<long>rand()) % (max_coeff - min_coeff + 1))

cdef int_max_plus_mat_set_identity(size_t n, long ** data):
    cdef size_t i
    for i in range(n * n):
        data[0][i] = LONG_MIN
    for i in range(n):
        data[i][i] = 0

cdef int_max_plus_mat_set(size_t n, long ** res, long ** m):
    cdef size_t i
    for i in range(n*n):
        res[0][i] = m[0][i]

def integer_max_plus_matrix_identity(size_t dim):
    r"""
    EXAMPLES::

        sage: from max_plus.max_plus_int import integer_max_plus_matrix_identity
        sage: integer_max_plus_matrix_identity(2)
        [   0 -oo ]
        [ -oo   0 ]
        sage: integer_max_plus_matrix_identity(3)
        [   0 -oo -oo ]
        [ -oo   0 -oo ]
        [ -oo -oo   0 ]
    """
    cdef IntegerMaxPlusMatrix ans = new_integer_max_plus_matrix(dim, dim)
    int_max_plus_mat_set_identity(dim, ans.data)
    return ans

def random_integer_max_plus_matrix(size_t dim, long min_coeff, long max_coeff,
        double minus_infinity_proba):
    r"""
    EXAMPLES::

        sage: from max_plus.max_plus_int import random_integer_max_plus_matrix
        sage: m = random_integer_max_plus_matrix(4, -2, 5, 0)
        sage: m  # random
        [  4  3 -2 2 ]
        [ -1 -1  2 4 ]
        [ -2  0 -1 1 ]
        [  1  4  4 3 ]
        sage: all(-2 <= x <= 5 for x in m.list())
        True
    """
    cdef IntegerMaxPlusMatrix ans = new_integer_max_plus_matrix(dim, dim)
    cdef size_t i
    cdef long mip = <long> round(minus_infinity_proba * RAND_MAX)
    for i in range(dim*dim):
        if mip and rand() < mip:
            ans.data[0][i] = LONG_MIN
        else:
            ans.data[0][i] = randlong(min_coeff, max_coeff)
    return ans

def random_integer_max_plus_matrices_band(size_t dim, long min_coeff, long
        max_coeff, char diag, char surdiag):
    r"""
    Return two triangular matrices.
    """
    cdef IntegerMaxPlusMatrix ans1 = new_integer_max_plus_matrix(dim, dim)
    cdef IntegerMaxPlusMatrix ans2 = new_integer_max_plus_matrix(dim, dim)
    cdef long ** m1 = ans1.data
    cdef long ** m2 = ans2.data
    cdef size_t i,j

    for i in range(dim):
        for j in range(dim):
            if i != j and j != i+1:
                m1[i][j] = m2[i][j] = LONG_MIN

    if diag == 'z':
        for i in range(dim):
            m1[i][i] = m2[i][i] = 0
    elif diag == 's':
        for i in range(dim):
            m1[i][i] = m2[i][i] = randlong(min_coeff, max_coeff)
    elif diag == 'v':
        for i in range(dim):
            m1[i][i] = randlong(min_coeff, max_coeff)
            m2[i][i] = randlong(min_coeff, max_coeff)
    else:
        raise ValueError

    if surdiag == 'z':
        for i in range(dim-1):
            m1[i][i+1] = m2[i][i+1] = 0
    elif surdiag == 's':
        for i in range(dim-1):
            m1[i][i+1] = m2[i][i+1] = randlong(min_coeff, max_coeff)
    elif surdiag == 'v':
        for i in range(dim-1):
            m1[i][i+1] = randlong(min_coeff, max_coeff)
            m2[i][i+1] = randlong(min_coeff, max_coeff)
    else:
        raise ValueError

    return ans1,ans2




def is_relation(tuple t1, tuple t2, list elements, bint upper=False):
    r"""
    Test if the relation ``t1 = t2`` is valid on the list of pairs ``elements``.

    INPUT:

    - ``t1``, ``t2`` -- two lists of 0 and 1 that encode a word in the free
      group

    - ``elements`` -- a list of pairs of integer max plus matrices

    - ``upper`` -- whether matrices are upper triangular (multiplication is
      faster)

    EXAMPLES::

        sage: from max_plus.max_plus_int import (
        ....:    random_integer_max_plus_matrix,
        ....:    random_integer_max_plus_matrices_band,
        ....:    is_relation)

        sage: elts = [random_integer_max_plus_matrices_band(3, 0, 10000, ord('s'), ord('v'))
        ....:            for _ in range(10)]
        sage: u = (0,0,1,0,1,1)
        sage: is_relation(u + (0,) + u, u + (1,) + u, elts, True)
        True
        sage: is_relation(u + (0,) + u, u + (1,) + u, elts, False)
        True
        sage: u = (0,1)
        sage: is_relation(u + (0,) + u, u + (1,) + u, elts, True)
        False

        sage: u = (0,0,1,1,1,0,0,0,1,0,1,0,1,1,1,0,0)
        sage: v = (0,0,1,1,1,0,1,0,1,0,0,0,1,1,1,0,0)
        sage: elts = [(random_integer_max_plus_matrix(2,-1000,1000,0.01), random_integer_max_plus_matrix(2,-1000,1000,0.01)) \
        ....:     for _ in range(100)]
        sage: is_relation(u, v, elts, False)
        True
    """
    cdef int n1, n2   # lengths of the relation
    cdef size_t dim   # matrix dimension
    cdef int * r1     # C copy of t1
    cdef int * r2     # C copy of t2
    cdef size_t i

    cdef IntegerMaxPlusMatrix m11,m12,m21,m22,e0,e1

    n1 = len(t1)
    n2 = len(t2)
    dim = (<IntegerMaxPlusMatrix> elements[0][0]).nrows

    r1 = <int *> malloc((n1+n2) * sizeof(int))
    r2 = r1 + n1

    for i in range(n1):
        r1[i] = t1[i]
        assert 0 <= r1[i] <= 1
    for i in range(n2):
        r2[i] = t2[i]
        assert 0 <= r2[i] <= 1

    m11 = new_integer_max_plus_matrix(dim, dim)
    m12 = new_integer_max_plus_matrix(dim, dim)
    m21 = new_integer_max_plus_matrix(dim, dim)
    m22 = new_integer_max_plus_matrix(dim, dim)
    for e0,e1 in elements:
        int_max_plus_mat_set_identity(dim, m11.data)
        int_max_plus_mat_set_identity(dim, m21.data)

        if upper:
            for i in range(n1):
                if r1[i] == 0:
                    int_max_plus_mat_prod_upper(dim, m12.data, m11.data, e0.data)
                else:
                    int_max_plus_mat_prod_upper(dim, m12.data, m11.data, e1.data)
                m11,m12 = m12,m11

            for i in range(n2):
                if r2[i] == 0:
                    int_max_plus_mat_prod_upper(dim, m22.data, m21.data, e0.data)
                else:
                    int_max_plus_mat_prod_upper(dim, m22.data, m21.data, e1.data)
                m21,m22 = m22,m21

        else:
            for i in range(n1):
                if r1[i] == 0:
                    int_max_plus_mat_prod(m12.data, m11.data, dim, dim, e0.data, dim, dim)
                else:
                    int_max_plus_mat_prod(m12.data, m11.data, dim, dim, e1.data, dim, dim)
                m11,m12 = m12,m11

            for i in range(n2):
                if r2[i] == 0:
                    int_max_plus_mat_prod(m22.data, m21.data, dim, dim, e0.data, dim, dim)
                else:
                    int_max_plus_mat_prod(m22.data, m21.data, dim, dim, e1.data, dim, dim)
                m21,m22 = m22,m21

        if m11 != m21:
            free(r1)
            return False
    free(r1)
    return True

cdef inline void prod_mat_upper(IntegerMaxPlusMatrix m11, IntegerMaxPlusMatrix m12,
                   IntegerMaxPlusMatrix m21, IntegerMaxPlusMatrix m22,
                   IntegerMaxPlusMatrix e0, IntegerMaxPlusMatrix e1,
                   int * r, int n):

    cdef size_t dim = m11.n

    int_max_plus_mat_set_identity(dim, m11.data)
    int_max_plus_mat_set_identity(dim, m21.data)

    for i in range(n):
        if r[i] == 0:
            int_max_plus_mat_prod_upper(dim, m12.data, m11.data, e0.data)
        else:
            int_max_plus_mat_prod_upper(dim, m12.data, m11.data, e1.data)
        m11,m12 = m12,m11

def filter_upper_relation(iterator,
        int n, int dim, list elements):
    r"""
    Filter ``iterator`` with candidates discarded by numerical computations.

    INPUT:

    - ``iterator`` -- iterator over the middle parts

    - ``n`` -- length of the relation

    - ``dim`` -- dimension of the matrices

    - ``nb_mats`` -- the number of pair matrices to use in the test

    - ``min_coeff``, ``max_coeff`` -- bound for coefficients of matrices

    EXAMPLES::

        sage: from max_plus.max_plus_int import filter_upper_relation,  random_integer_max_plus_matrices_band


        sage: from itertools import product
        sage: def my_iterator():
        ....:     p = (0, 1, 1, 0, 0)
        ....:     s = (1, 1, 1, 0, 0)
        ....:     for m1 in product((0,1), repeat=9):
        ....:         for m2 in product((0,1), repeat=9):
        ....:             if m1 == m2:
        ....:                 break
        ....:             yield (p + m1 + s, p + m2 + s)
        sage: elements = [random_integer_max_plus_matrices_band(4, -2**10, 2**10, ord('s'), ord('v')) for _ in range(100)]
        sage: it = filter_upper_relation(my_iterator(), 19, 4, elements)
        sage: it.next()  # random
        ((0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 0, 0),
         (0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0))

    TESTS::

        sage: from max_plus.sv_identities import is_sv_identity
        sage: def my_iterator():
        ....:     W = [(0,0,1,1,0),(0,1,1,0,0),(1,0,0,1,1),(1,1,0,0,1)]
        ....:     for p in W:
        ....:         for s in W:
        ....:             yield (p+(0,)+s, p+(1,)+s)
        ....:     p = (0,1,0,1,1)
        ....:     S = [(1,0,0,1,1), (1,1,0,0,1), (1,1,0,1,0)]
        ....:     for s in S:
        ....:         yield (p+(0,)+s, p+(1,)+s)

        sage: for _ in range(1000):
        ....:     elements = [random_integer_max_plus_matrices_band(3, -2**10, 2**10, ord('s'), ord('v'))]
        ....:     s1 = set(my_iterator())
        ....:     s2 = set(filter_upper_relation(my_iterator(), 11, 3, elements))
        ....:     assert s1 == s2, "\nm1={}\n\nm2={}".format(
        ....:                elements[0][0].list(),
        ....:                elements[0][1].list())
    """
    cdef MemoryAllocator mem = MemoryAllocator()
    cdef int * r1     # C copy of left part of relations
    cdef int * r2     # C copy of right part of relations

    cdef int i
    cdef tuple t1,t2

    cdef list p_mats
    cdef list s_mats

    cdef IntegerMaxPlusMatrix m1,m2,mt,e0,e1

    assert n > 0

    r1 = <int *> mem.allocarray(2*n,  sizeof(int))
    r2 = r1 + n

    m1 = new_integer_max_plus_matrix(dim, dim)
    m2 = new_integer_max_plus_matrix(dim, dim)
    mt = new_integer_max_plus_matrix(dim, dim)

    int_max_plus_mat_set_identity(dim, m1.data)
    int_max_plus_mat_set_identity(dim, m2.data)

    for t in iterator:
        t1,t2 = t
        assert len(t1) == len(t2) == n

        for i in range(n):
            r1[i] = t1[i]
            assert 0 <= r1[i] <= 1
        for i in range(n):
            r2[i] = t2[i]
            assert 0 <= r2[i] <= 1

        for e0,e1 in elements:
            # product for t1
            if r1[0] == 0:
                int_max_plus_mat_set(dim, m1.data, e0.data)
            else:
                int_max_plus_mat_set(dim, m1.data, e1.data)
            for i in range(1, n):
                m1,mt = mt,m1
                if r1[i] == 0:
                    int_max_plus_mat_prod_upper(dim, m1.data, mt.data, e0.data)
                else:
                    int_max_plus_mat_prod_upper(dim, m1.data, mt.data, e1.data)

            # product for t2
            if r2[0] == 0:
                int_max_plus_mat_set(dim, m2.data, e0.data)
            else:
                int_max_plus_mat_set(dim, m2.data, e1.data)
            for i in range(1, n):
                m2,mt = mt,m2
                if r2[i] == 0:
                    int_max_plus_mat_prod_upper(dim, m2.data, mt.data, e0.data)
                else:
                    int_max_plus_mat_prod_upper(dim, m2.data, mt.data, e1.data)

            if m1 != m2:
                break

        else:
            yield t

cdef IntegerMaxPlusMatrix new_integer_max_plus_matrix(size_t nrows, size_t ncols):
    cdef IntegerMaxPlusMatrix ans = IntegerMaxPlusMatrix.__new__(IntegerMaxPlusMatrix)
    ans.nrows = nrows
    ans.ncols = ncols
    ans.data = <long **> malloc(ans.nrows*sizeof(long *))
    cdef size_t i
    ans.data[0] = <long *> malloc(ans.nrows*ans.ncols*sizeof(long))
    for i in range(ans.nrows-1):
        ans.data[i+1] = ans.data[i] + ans.ncols
    return ans

cdef int_max_plus_mat_prod(long ** ans, long ** a, size_t a_nrows, size_t a_ncols,
                                        long ** b, size_t b_nrows, size_t b_ncols):
    cdef size_t i, j, k
    cdef long c, c2
    for i in range(a_nrows):
        for j in range(b_ncols):
            c = LONG_MIN
            for k in range(a_ncols):
                if a[i][k] != LONG_MIN and b[k][j] != LONG_MIN:
                    c2 = a[i][k] + b[k][j]
                    if c2 > c:
                        c = c2
            ans[i][j] = c

cdef int_max_plus_mat_prod_upper(size_t n, long ** ans, long **a, long **b):
    r"""
    product of upper triangular matrices
    """
    cdef size_t i, j, k
    cdef long c, c2
    for i in range(n):
        for j in range(i):
            ans[i][j] = LONG_MIN
        for j in range(i,n):
            c = LONG_MIN
            for k in range(i,j+1):
                if a[i][k] != LONG_MIN and b[k][j] != LONG_MIN:
                    c2 = a[i][k] + b[k][j]
                    if c2 > c:
                        c = c2
            ans[i][j] = c

cdef class IntegerMaxPlusMatrix:
    r"""
    Max-plus matrix with integer coefficients.
    """
    cdef size_t nrows, ncols    # number of rows, columns
    cdef long ** data           # entries

    def __init__(self, nrows, ncols, data, upper=False):
        r"""
        EXAMPLES::

            sage: from max_plus.max_plus_int import IntegerMaxPlusMatrix
            sage: IntegerMaxPlusMatrix(0,2,[])
            []
            sage: IntegerMaxPlusMatrix(0,2,[0,1,-3,-1])
            Traceback (most recent call last):
            ...
            ValueError: invalid input data
            sage: IntegerMaxPlusMatrix(0,2,[[0,1,-3],[-1]])
            Traceback (most recent call last):
            ...
            ValueError: invalid input data
        """
        cdef int one_list
        self.nrows = <size_t?> nrows
        self.ncols = <size_t?> ncols
        cdef size_t i,j
        if self.nrows and self.ncols:
            if not isinstance(data, (tuple,list)) or not data:
                raise ValueError("data must be a non-empty list")
            if isinstance(data[0], (tuple,list)):
                if any(not isinstance(x, (tuple,list)) or len(x) != self.ncols for x in data):
                    raise ValueError("invalid input data")
                data = [int(data[i][j]) for i in range(self.nrows) for j in range(self.ncols)]
            else:
                if len(data) != self.nrows*self.ncols:
                    raise ValueError("invalid input data")
                data = map(int,data)

            self.data = <long **> malloc(self.nrows*sizeof(long *))
            self.data[0] = <long *> malloc(self.nrows*self.ncols*sizeof(long))
            for i in range(self.nrows-1):
                self.data[i+1] = self.data[i] + self.ncols
            for i in range(self.nrows*self.ncols):
                self.data[0][i] = data[i]

        elif data:
            raise ValueError("invalid input data")
        else:
            self.data = NULL

    def barvinok_rank(self):
        r"""
        Compute the Barvinok rank of a matrix.

        .. NOTE::

            this is currently very slow because of the startup time of polymake!

        EXAMPLES::

            sage: from max_plus.max_plus_int import IntegerMaxPlusMatrix

            sage: m = IntegerMaxPlusMatrix(3, 3, [1,2,3,4,5,6,7,8,9])
            sage: m.barvinok_rank() # not tested
            1

            sage: m = IntegerMaxPlusMatrix(3, 3, [2,3,-2,1,3,5,4,5,1])
            sage: m.barvinok_rank() # not tested
            2

            sage: m = IntegerMaxPlusMatrix(3, 3, [-2,-5,-3,3,-1,-3,-1,-2,5])
            sage: m.barvinok_rank() # not tested
            3
        """
        from max_plus.rank import encode_barvinok, toric_rank
        matrix = [[self.data[i][j] for j in range(self.ncols)] for i in range(self.nrows)]
        coords = encode_barvinok(matrix)
        return toric_rank(coords)

    def _latex_(self):
        r"""
        EXAMPLES::

            sage: from max_plus.max_plus_int import IntegerMaxPlusMatrix, minus_infinity
            sage: mi = minus_infinity()
            sage: m = IntegerMaxPlusMatrix(3, 2, [0,1,mi,3,1,mi])
            sage: latex(m)
            \left(\begin{array}{rr}
            0 & 1 \\
            -\infty & 3 \\
            1 & -\infty
            \end{array}\right)
        """
        s = []
        cdef Py_ssize_t i
        for i in range(self.nrows):
            line = []
            for j in range(self.ncols):
                if self.data[i][j] == LONG_MIN:
                    line.append('-\\infty')
                else:
                    line.append(str(self.data[i][j]))
            s.append(' & '.join(line))

        start = "\\left(\\begin{{array}}{{{}}}\n".format('r'*self.ncols)
        end = '\n\\end{array}\\right)'

        return start + ' \\\\\n'.join(s) + end

    def __hash__(self):
        r"""
        TESTS::

            sage: from max_plus.max_plus_int import IntegerMaxPlusMatrix
            sage: m1 = IntegerMaxPlusMatrix(2, 2, [0,1,2,3])
            sage: m2 = IntegerMaxPlusMatrix(2, 2, [1,0,2,3])
            sage: m3 = IntegerMaxPlusMatrix(2, 2, [3,2,1,0])
            sage: hash(m1) != hash(m3) and hash(m1) != hash(m2) and hash(m2) != hash(m3)
            True

            sage: m1 = IntegerMaxPlusMatrix(2, 2, [0,1,2,3])
            sage: m2 = IntegerMaxPlusMatrix(2, 2, [0,1,2,3])
            sage: m3 = IntegerMaxPlusMatrix(2, 2, [0,1,2,3])
            sage: hash(m1) == hash(m2) == hash(m3)
            True
        """
        cdef int i
        cdef long mult = 1000003L
        cdef long x,y
        x = 0x345678L
        for i in range(self.nrows * self.ncols):
            y = self.data[0][i]
            x = (x^y) * mult
            mult += <long>(82520L + self.nrows + self.ncols)
        return x+97531L

    def __getitem__(self, key):
        cdef Py_ssize_t i,j
        i,j = key
        if i < 0 or i >= self.nrows or j < 0 or j >= self.ncols:
            raise ValueError("index out of range")
        return self.data[i][j]

    def transpose(self):
        r"""
        Return the transpose matrix

        EXAMPLES::

            sage: from max_plus.max_plus_int import IntegerMaxPlusMatrix
            sage: IntegerMaxPlusMatrix(3, 3, range(9)).transpose()
            [ 0 3 6 ]
            [ 1 4 7 ]
            [ 2 5 8 ]
            sage: IntegerMaxPlusMatrix(2, 3, range(6)).transpose()
            [ 0 3 ]
            [ 1 4 ]
            [ 2 5 ]
            sage: IntegerMaxPlusMatrix(3, 2, range(6)).transpose()
            [ 0 2 4 ]
            [ 1 3 5 ]
        """
        cdef IntegerMaxPlusMatrix ans = new_integer_max_plus_matrix(self.ncols, self.nrows)
        cdef size_t i, j

        for i in range(self.nrows):
            for j in range(self.ncols):
                ans.data[j][i] = self.data[i][j]
        return ans

    def list(self):
        r"""
        Return the list of coefficients.

        The output is a list of ``nrows * ncols`` integers.

        EXAMPLES::

            sage: from max_plus.max_plus_int import IntegerMaxPlusMatrix
            sage: IntegerMaxPlusMatrix(2, 3, [[0,1,3],[1,1,1]]).list()
            [0, 1, 3, 1, 1, 1]
        """
        cdef int i
        return [self.data[0][i] for i in range(self.nrows*self.ncols)]

    def __dealloc__(self):
        if self.data != NULL:
            free(self.data[0])
            free(self.data)

    def __repr__(self):
        r"""
        EXAMPLES::

            sage: from max_plus.max_plus_int import IntegerMaxPlusMatrix, minus_infinity
            sage: mo = minus_infinity()

            sage: IntegerMaxPlusMatrix(2, 3, [1,mo,2,4,5,mo])
            [ 1 -oo   2 ]
            [ 4   5 -oo ]
            sage: IntegerMaxPlusMatrix(3, 1, [1,2,mo])
            [   1 ]
            [   2 ]
            [ -oo ]
        """
        if self.ncols == 0 or self.nrows == 0:
            return "[]"

        col_sizes = []
        for j in range(self.ncols):
            c = 0
            for i in range(self.nrows):
                if self.data[i][j] == LONG_MIN:
                    c = max(c, 3)
                else:
                    c = max(c, len("{:}".format(self.data[i][j])))
            col_sizes.append(c)

        l = []
        for i in range(self.nrows):
            s = "[ "
            for j in range(self.ncols):
                c = "-oo" if self.data[i][j] == LONG_MIN else self.data[i][j]
                s += "{:>{width}} ".format(c, width=col_sizes[j])
            s += "]"
            l.append(s)

        return "\n".join(l)

    def infinite_coefficients(self):
        r"""
        Return the list of pairs `(i,j)` so that the corresponding entry is
        -infinity.

        EXAMPLES::

            sage: from max_plus.max_plus_int import IntegerMaxPlusMatrix, minus_infinity
            sage: mo = minus_infinity()

            sage: IntegerMaxPlusMatrix(2, 3, [1,mo,2,4,5,mo]).infinite_coefficients()
            [(0, 1), (1, 2)]
            sage: IntegerMaxPlusMatrix(3, 1, [1,2,mo]).infinite_coefficients()
            [(2, 0)]
        """
        return [(i,j) for i in range(self.nrows) for j in range(self.ncols) if self.data[i][j] == LONG_MIN]

    def __richcmp__(_self, _other, op):
        r"""
        TESTS::

            sage: from max_plus.max_plus_int import IntegerMaxPlusMatrix

            sage: IntegerMaxPlusMatrix(2,3,[0]*6) == IntegerMaxPlusMatrix(2,3,[0]*6)
            True
            sage: IntegerMaxPlusMatrix(2,3,[0]*6) != IntegerMaxPlusMatrix(2,3,[0]*6)
            False
            sage: IntegerMaxPlusMatrix(2,3,[0]*6) == IntegerMaxPlusMatrix(3,2,[0]*6)
            False
            sage: IntegerMaxPlusMatrix(2,3,[0]*6) != IntegerMaxPlusMatrix(3,2,[0]*6)
            True
            sage: IntegerMaxPlusMatrix(2,3,[0]*6) != IntegerMaxPlusMatrix(2,3,[1]*6)
            True
        """
        cdef IntegerMaxPlusMatrix self = <IntegerMaxPlusMatrix> _self
        cdef IntegerMaxPlusMatrix other = <IntegerMaxPlusMatrix?> _other
        if op != Py_EQ and op != Py_NE:
            raise TypeError

        if self.nrows != other.nrows or self.ncols != other.ncols:
            return op == Py_NE

        cdef size_t i,j
        for i in range(self.nrows):
            for j in range(self.ncols):
                if self.data[i][j] != other.data[i][j]:
                    return op == Py_NE
        return op == Py_EQ

    def __mul__(_self, _other):
        r"""
        Multiplication of matrices

        EXAMPLES::

            sage: from max_plus.max_plus_int import IntegerMaxPlusMatrix, minus_infinity
            sage: mo = minus_infinity()

            sage: m1 = IntegerMaxPlusMatrix(2, 1, [0,1])
            sage: m2 = IntegerMaxPlusMatrix(2, 2, [1,3,mo,2])
            sage: m3 = IntegerMaxPlusMatrix(1, 2, [1,1])

            sage: m2 * m1
            [ 4 ]
            [ 3 ]
            sage: m3 * m2
            [ 2 4 ]
            sage: m1 * m3
            [ 1 1 ]
            [ 2 2 ]
            sage: m3 * m1
            [ 2 ]

            sage: m1 * m2
            Traceback (most recent call last):
            ...
            ValueError: can not multiply 2x1 matrix with 2x2 matrix
            sage: m2 * m3
            Traceback (most recent call last):
            ...
            ValueError: can not multiply 2x2 matrix with 1x2 matrix
        """
        cdef IntegerMaxPlusMatrix self = <IntegerMaxPlusMatrix?> _self
        cdef IntegerMaxPlusMatrix other = <IntegerMaxPlusMatrix?> _other

        if self.ncols != other.nrows:
            raise ValueError("can not multiply {}x{} matrix with {}x{} matrix".format(
                self.nrows, self.ncols, other.nrows, other.ncols))

        cdef IntegerMaxPlusMatrix ans = new_integer_max_plus_matrix(self.nrows, other.ncols)
        int_max_plus_mat_prod(ans.data, self.data, self.nrows, self.ncols,
                                 other.data, other.nrows, other.ncols)
        return ans

    def __pow__(_self, n, mod):
        r"""
        EXAMPLES::

            sage: from max_plus.max_plus_int import IntegerMaxPlusMatrix
            sage: m = IntegerMaxPlusMatrix(2, 2, [1,2,3,4])
            sage: m**0
            [   0 -oo ]
            [ -oo   0 ]
            sage: m**1
            [ 1 2 ]
            [ 3 4 ]
            sage: m**2
            [ 5 6 ]
            [ 7 8 ]
            sage: m**3
            [  9 10 ]
            [ 11 12 ]
            sage: m**4
            [ 13 14 ]
            [ 15 16 ]


        TESTS::

            sage: m = IntegerMaxPlusMatrix(2, 1, [3,2])
            sage: m**2
            Traceback (most recent call last):
            ...
            ValueError: can not take powers of non square matrix

            sage: from max_plus.max_plus_int import (random_integer_max_plus_matrix,
            ....:          integer_max_plus_matrix_identity)
            sage: for _ in range(100):
            ....:     dim = randint(1,10)
            ....:     k = randint(0,30)
            ....:     m = random_integer_max_plus_matrix(dim, -100, 100, 0.1)
            ....:     p = integer_max_plus_matrix_identity(dim)
            ....:     for i in range(k): p = p*m
            ....:     assert p == m**k
        """
        cdef IntegerMaxPlusMatrix self = _self
        if self.nrows != self.ncols:
            raise ValueError("can not take powers of non square matrix")

        cdef unsigned long m = <unsigned long?> n
        cdef int dim = self.nrows
        cdef IntegerMaxPlusMatrix power, apow, tmp

        power = new_integer_max_plus_matrix(dim, dim)
        int_max_plus_mat_set_identity(dim, power.data)
        if not m:
            return power

        apow = new_integer_max_plus_matrix(dim, dim)
        tmp = new_integer_max_plus_matrix(dim, dim)

        memcpy(apow.data[0], self.data[0], dim*dim*sizeof(long))
        while m&1 == 0:
            int_max_plus_mat_prod(tmp.data, apow.data, dim, dim, apow.data, dim, dim)
            apow,tmp = tmp,apow
            m = m >> 1
        memcpy(power.data[0], apow.data[0], dim*dim*sizeof(long))
        m = m >> 1

        while m:
            int_max_plus_mat_prod(tmp.data, apow.data, dim, dim, apow.data, dim, dim)
            apow,tmp = tmp,apow
            if m&1:
                int_max_plus_mat_prod(tmp.data, power.data, dim, dim, apow.data, dim, dim)
                power,tmp = tmp,power
            m = m >> 1

        return power


cdef class IntegerMatrixProduct(object):
    r"""
    Fast integer max-plus matrix products

    EXAMPLES::

        sage: from max_plus.max_plus_int import IntegerMatrixProduct, IntegerMaxPlusMatrix
        sage: t = (0,1,0,0,1,0)
        sage: p = IntegerMatrixProduct(t)
        sage: p
        product along 010010
        sage: a1 = IntegerMaxPlusMatrix(2, 2, [1,2,3,4])
        sage: a2 = IntegerMaxPlusMatrix(2, 2, [1,0,1,1])
        sage: a3 = IntegerMaxPlusMatrix(2, 2, [1,3,0,2])
        sage: p(a1, a2) == prod(a1 if i == 0 else a2 for i in t)
        True
        sage: p(a1, a3) == prod(a1 if i == 0 else a3 for i in t)
        True
        sage: p(a2, a3) == prod(a2 if i == 0 else a3 for i in t)
        True

    The class can also be initialized with strings::

        sage: p = IntegerMatrixProduct('01101010')
        sage: p(a1, a2)
        [ 17 18 ]
        [ 19 20 ]

    TESTS::

        sage: p = IntegerMatrixProduct(3)
        Traceback (most recent call last):
        ...
        ValueError: can only be initialized from tuple, list or string
    """
    cdef int * t
    cdef int n

    def __cinit__(self, t):
        if not isinstance(t, (tuple, list, str)):
            raise ValueError("can only be initialized from tuple, list or string")
        self.n = len(t)
        if not t:
            self.t = NULL
        else:
            self.t = <int *> calloc(len(t), sizeof(int))

    def __dealloc__(self):
        if self.n:
            free(self.t)

    def __init__(self, t):
        if isinstance(t, (tuple, list)):
            for i,j in enumerate(t):
                self.t[i] = j
        elif isinstance(t, str):
            for i,j in enumerate(t):
                self.t[i] = int(j)

    def __repr__(self):
        cdef int i
        return "product along " + "".join("1" if self.t[i] else "0" for i in range(self.n))

    def __call__(self, IntegerMaxPlusMatrix a, IntegerMaxPlusMatrix b):
        if not (a.nrows == a.ncols == b.nrows == b.ncols):
            raise ValueError("should be square matrices of the same size")

        cdef int dim = a.nrows
        cdef IntegerMaxPlusMatrix ans1 = new_integer_max_plus_matrix(dim, dim)
        cdef IntegerMaxPlusMatrix ans2 = new_integer_max_plus_matrix(dim, dim)

        int_max_plus_mat_set_identity(dim, ans1.data)

        for i in range(self.n):
            if self.t[i] == 0:
                int_max_plus_mat_prod(ans2.data, ans1.data, dim, dim, a.data, dim, dim)
            else:
                int_max_plus_mat_prod(ans2.data, ans1.data, dim, dim, b.data, dim, dim)
            ans1,ans2 = ans2,ans1

        return ans1
