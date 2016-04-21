r"""
Integer max plus matrices.

These are integer matrices with *positive* coefficients (-1 is reserved for
-infinity).
"""
from libc.stdlib cimport malloc, free, rand, RAND_MAX
# note: on my computer RAND_MAX is 2**31 - 1
from libc.limits cimport LONG_MIN, LONG_MAX
from libc.math cimport round

from cpython.object cimport Py_EQ, Py_NE

from sage.structure.element cimport generic_power_c
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
    cdef IntegerMaxPlusMatrix ans = new_integer_max_plus_matrix(dim)
    int_max_plus_mat_set_identity(ans.n, ans.data)
    return ans

def random_integer_max_plus_matrix(size_t dim, long min_coeff, long max_coeff,
        double minus_infinity_proba):
    cdef IntegerMaxPlusMatrix ans = new_integer_max_plus_matrix(dim)
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
    cdef IntegerMaxPlusMatrix ans1 = new_integer_max_plus_matrix(dim)
    cdef IntegerMaxPlusMatrix ans2 = new_integer_max_plus_matrix(dim)
    cdef long ** m1 = ans1.data
    cdef long ** m2 = ans2.data
    ans1.upper = ans2.upper = 1
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

        sage: from max_plus import *

        sage: elts = [random_integer_max_plus_matrices_band(3, 0, 10000, ord('s'), ord('v'))
        ....:            for _ in range(10)]
        sage: u = (0,0,1,0,1,1)
        sage: is_relation(u + (0,) + u, u + (1,) + u, elts, True)
        True
        sage: u = (0,1)
        sage: is_relation(u + (0,) + u, u + (1,) + u, elts, True)
        False
    """
    cdef int n1, n2   # lengths of the relation
    cdef size_t dim   # matrix dimension
    cdef int * r1     # C copy of t1
    cdef int * r2     # C copy of t2
    cdef size_t i

    cdef IntegerMaxPlusMatrix m11,m12,m21,m22,e0,e1

    n1 = len(t1)
    n2 = len(t2)
    dim = (<IntegerMaxPlusMatrix> elements[0][0]).n

    r1 = <int *> malloc((n1+n2) * sizeof(int))
    r2 = r1 + n1

    for i in range(n1):
        r1[i] = t1[i]
        assert 0 <= r1[i] <= 1
    for i in range(n2):
        r2[i] = t2[i]
        assert 0 <= r2[i] <= 1

    m11 = new_integer_max_plus_matrix(dim)
    m12 = new_integer_max_plus_matrix(dim)
    m21 = new_integer_max_plus_matrix(dim)
    m22 = new_integer_max_plus_matrix(dim)
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
                    int_max_plus_mat_prod(dim, m12.data, m11.data, e0.data)
                else:
                    int_max_plus_mat_prod(dim, m12.data, m11.data, e1.data)
                m11,m12 = m12,m11

            for i in range(n2):
                if r2[i] == 0:
                    int_max_plus_mat_prod(dim, m22.data, m21.data, e0.data)
                else:
                    int_max_plus_mat_prod(dim, m22.data, m21.data, e1.data)
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

def filter_sv_relation(iterator,
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

        sage: from max_plus.max_plus_int import filter_sv_relation, random_integer_max_plus_matrices_band

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
        sage: it = filter_sv_relation(my_iterator(), 19, 4, elements)
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
        ....:     s2 = set(filter_sv_relation(my_iterator(), 11, 3, elements))
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

    m1 = new_integer_max_plus_matrix(dim)
    m2 = new_integer_max_plus_matrix(dim)
    mt = new_integer_max_plus_matrix(dim)

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

cdef IntegerMaxPlusMatrix new_integer_max_plus_matrix(size_t dim):
    cdef IntegerMaxPlusMatrix ans = IntegerMaxPlusMatrix.__new__(IntegerMaxPlusMatrix)
    ans.n = dim
    ans.data = <long **> malloc(ans.n*sizeof(long *))
    cdef size_t i
    ans.data[0] = <long *> malloc(ans.n*ans.n*sizeof(long))
    for i in range(ans.n-1):
        ans.data[i+1] = ans.data[i] + ans.n
    return ans

cdef int_max_plus_mat_prod(size_t n, long ** ans, long ** a, long ** b):
    cdef size_t i, j, k
    cdef long c, c2
    for i in range(n):
        for j in range(n):
            c = LONG_MIN
            for k in range(n):
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
    cdef size_t n
    cdef long ** data
    cdef int upper

    def __init__(self, n, data, upper=False):
        self.n = n
        self.data = <long **> malloc(n*sizeof(long *))
        cdef size_t i
        if self.n:
            self.data[0] = <long *> malloc(n*n*sizeof(long))
            for i in range(self.n-1):
                self.data[i+1] = self.data[i] + self.n

        if len(data) == n:
            if any(len(x) != n for x in data):
                raise ValueError
            for i in range(n):
                for j in range(n):
                    self.data[i][j] = data[i][j]
        elif len(data) == n*n:
            for i in range(n):
                for j in range(n):
                    self.data[i][j] = data[n*i + j]

        self.upper = upper

    def __getitem__(self, key):
        cdef Py_ssize_t i,j
        i,j = key
        if i < 0 or i >= self.n or j < 0 or j >= self.n:
            raise ValueError("index out of range")
        return self.data[i][j]

    def list(self):
        cdef int i
        return [self.data[0][i] for i in range(self.n*self.n)]

    def is_upper(self):
        return self.upper == 1

    def __dealloc__(self):
        if self.data != NULL:
            free(self.data[0])
            free(self.data)

    def __repr__(self):
        col_sizes = []
        for j in range(self.n):
            c = 0
            for i in range(self.n):
                if self.data[i][j] == LONG_MIN:
                    c = max(c, 3)
                else:
                    c = max(c, len("{:}".format(self.data[i][j])))
            col_sizes.append(c)

        l = []
        for i in range(self.n):
            s = "[ "
            for j in range(self.n):
                c = "-oo" if self.data[i][j] == LONG_MIN else self.data[i][j]
                s += "{:>{width}} ".format(c, width=col_sizes[j])
            s += "]"
            l.append(s)

        return "\n".join(l)

    def infinite_coefficients(self):
        return [(i,j) for i in range(self.n) for j in range(self.n) if self.data[i][j] == LONG_MIN]

    def __richcmp__(_self, _other, op):
        cdef IntegerMaxPlusMatrix self = <IntegerMaxPlusMatrix> _self
        cdef IntegerMaxPlusMatrix other = <IntegerMaxPlusMatrix?> _other
        if op != Py_EQ and op != Py_NE:
            raise TypeError

        if self.n != other.n:
            return False

        cdef size_t i,j
        if self.upper and other.upper:
            for i in range(self.n):
                for j in range(i, self.n):
                    if self.data[i][j] != other.data[i][j]:
                        return op == Py_NE
        else:
            for i in range(self.n):
                for j in range(self.n):
                    if self.data[i][j] != other.data[i][j]:
                        return op == Py_NE
        return op == Py_EQ

    def __mul__(_self, _other):
        r"""
        Multiplication of matrices
        """
        cdef IntegerMaxPlusMatrix self = <IntegerMaxPlusMatrix?> _self
        cdef IntegerMaxPlusMatrix other = <IntegerMaxPlusMatrix?> _other
        cdef IntegerMaxPlusMatrix ans = new_integer_max_plus_matrix(self.n)
        if self.upper and other.upper:
            int_max_plus_mat_prod_upper(self.n, ans.data, self.data, other.data)
            ans.upper = 1
        else:
            int_max_plus_mat_prod(self.n, ans.data, self.data, other.data)
            ans.upper = 0
        return ans

    def __pow__(self, n, mod):
        return generic_power_c(self, n, integer_max_plus_matrix_identity(n))

