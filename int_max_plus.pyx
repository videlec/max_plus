r"""
Integer max plus matrices.

These are integer matrices with *positive* coefficients (-1 is reserved for
-infinity).
"""
from libc.stdlib cimport malloc, free, rand, RAND_MAX
# note: on my computer RAND_MAX is 2**31 - 1
from libc.limits cimport LONG_MIN, LONG_MAX

from cpython.object cimport Py_EQ, Py_NE

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

def integer_max_plus_matrix_identity(size_t dim):
    cdef IntegerMaxPlusMatrix ans = new_integer_max_plus_matrix(dim)
    int_max_plus_mat_set_identity(ans.n, ans.data)
    return ans

def random_integer_max_plus_matrix(size_t dim, long min_coeff, long max_coeff):
    cdef IntegerMaxPlusMatrix ans = new_integer_max_plus_matrix(dim)
    cdef size_t i
    for i in range(dim*dim):
        ans.data[0][i] = randlong(min_coeff, max_coeff)
    return ans

def random_integer_max_plus_matrices_band(size_t dim, long min_coeff, long
        max_coeff, char diag, char surdiag):
    r"""
    Return two triangular matrices with the same diagonal.
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
    Test if the relation ``r1 = r2`` is valid on the list of pairs ``elements``.

    INPUT:

    - ``r1``, ``r2`` -- two lists of 0 and 1 that encode a word in the free
      group

    - ``elements`` -- a list of pairs of integer max plus matrices

    - ``upper`` -- whether matrices are upper triangular
    """
    cdef int n = len(t1)
    cdef size_t dim = (<IntegerMaxPlusMatrix> elements[0][0]).n
    cdef int * r1 = <int *> malloc(2 * n * sizeof(int))
    cdef int * r2 = r1 + n
    cdef size_t i

    for i in range(n):
        r1[i] = t1[i]
        r2[i] = t2[i]

    cdef IntegerMaxPlusMatrix m11,m12,m21,m22,e0,e1
    m11 = new_integer_max_plus_matrix(dim)
    m12 = new_integer_max_plus_matrix(dim)
    m21 = new_integer_max_plus_matrix(dim)
    m22 = new_integer_max_plus_matrix(dim)
    for e0,e1 in elements:
        int_max_plus_mat_set_identity(dim, m11.data)
        int_max_plus_mat_set_identity(dim, m21.data)

        if upper:
            for i in range(n):
                if r1[i] == 0:
                    int_max_plus_mat_prod_upper(dim, m12.data, m11.data, e0.data)
                else:
                    int_max_plus_mat_prod_upper(dim, m12.data, m11.data, e1.data)
                if r2[i] == 0:
                    int_max_plus_mat_prod_upper(dim, m22.data, m21.data, e0.data)
                else:
                    int_max_plus_mat_prod_upper(dim, m22.data, m21.data, e1.data)
                m11,m12 = m12,m11
                m21,m22 = m22,m21

        else:
            for i in range(n):
                if r1[i] == 0:
                    int_max_plus_mat_prod(dim, m12.data, m11.data, e0.data)
                else:
                    int_max_plus_mat_prod(dim, m12.data, m11.data, e1.data)
                if r2[i] == 0:
                    int_max_plus_mat_prod(dim, m22.data, m21.data, e0.data)
                else:
                    int_max_plus_mat_prod(dim, m22.data, m21.data, e1.data)
                m11,m12 = m12,m11
                m21,m22 = m22,m21

        if m11 != m21:
            free(r1)
            return False
    free(r1)
    return True

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
    product of 4x4 upper triangular matrices
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



############################
# some efficient iterators #
############################
 
def product_p(int n0, int n1):
    r"""
    Run through all tuples with ``n0`` zeros and ``n1`` ones.

    INPUT:

    - ``n0``, ``n1`` -- the number of ``0`` and ``1``

    EXAMPLES::

        sage: for t in product_p(2,3): print t
        [0, 0, 1, 1, 1]
        [0, 1, 0, 1, 1]
        [0, 1, 1, 0, 1]
        [0, 1, 1, 1, 0]
        [1, 0, 0, 1, 1]
        [1, 0, 1, 0, 1]
        [1, 0, 1, 1, 0]
        [1, 1, 0, 0, 1]
        [1, 1, 0, 1, 0]
        [1, 1, 1, 0, 0]
    """
    cdef list i = []
    cdef int j

    while True:
        for j in range(n0):
            i.append(0)
        for j in range(n1):
            i.append(1)
        n0 = 0
        n1 = 0
        yield i

        while i and (i[-1] == 1 or n1 == 0):
            if i.pop() == 0:
                n0 += 1
            else:
                n1 += 1
        if not i:
            return
        n0 += 1
        n1 -= 1
        i[-1] = 1
