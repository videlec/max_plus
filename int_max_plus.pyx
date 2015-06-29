r"""
Integer max plus matrices.

These are integer matrices with *positive* coefficients (-1 is reserved for
-infinity).
"""
from libc.stdlib cimport malloc, free
from libc.limits cimport LONG_MIN, LONG_MAX

from cpython.object cimport Py_EQ, Py_NE

def minus_infinity():
    return <long>LONG_MIN

def integer_max_plus_matrix_identity(dim):
    cdef IntegerMaxPlusMatrix ans = new_integer_max_plus_matrix(dim)
    cdef size_t i
    for i in range(ans.n * ans.n):
        ans.data[0][i] = LONG_MIN
    for i in range(ans.n):
        ans.data[i][i] = 0
    return ans

def random_integer_max_plus_matrix(dim, min_coeff, max_coeff):
    from random import randint
    return IntegerMaxPlusMatrix(dim, [randint(min_coeff, max_coeff) for _ in range(dim*dim)])

def random_integer_max_plus_matrix_tri(dim, min_coeff, max_coeff):
    from random import randint
    cdef IntegerMaxPlusMatrix ans = new_integer_max_plus_matrix(dim)
    cdef size_t i
    for i in range(ans.n):
        for j in range(i):
            ans.data[i][j] = LONG_MIN
        for j in range(i, ans.n):
            ans.data[i][j] = randint(min_coeff, max_coeff)
    return ans

def random_integer_max_plus_matrix_tri_sym_diag(dim, min_coeff, max_coeff):
    raise NotImplementedError

cdef IntegerMaxPlusMatrix new_integer_max_plus_matrix(size_t dim):
    cdef IntegerMaxPlusMatrix ans = IntegerMaxPlusMatrix.__new__(IntegerMaxPlusMatrix)
    ans.n = dim
    ans.data = <long **> malloc(ans.n*sizeof(long *))
    cdef size_t i
    ans.data[0] = <long *> malloc(ans.n*ans.n*sizeof(long))
    for i in range(ans.n-1):
        ans.data[i+1] = ans.data[i] + ans.n
    return ans

cdef class IntegerMaxPlusMatrix:
    cdef size_t n
    cdef long ** data

    def __init__(self, n, data):
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
        cdef size_t i,j
        if op != Py_EQ and op != Py_NE:
            raise TypeError

        if self.n != other.n:
            return False
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
        cdef size_t n = self.n
        cdef size_t i, j

        cdef IntegerMaxPlusMatrix ans = new_integer_max_plus_matrix(n)

        cdef long c, c2
        for i in range(n):
            for j in range(n):
                c = LONG_MIN
                for k in range(n):
                    if self.data[i][k] != LONG_MIN and other.data[k][j] != LONG_MIN:
                        c2 = self.data[i][k] + other.data[k][j]
                        if c2 > c:
                            c = c2
                ans.data[i][j] = c

        return ans

