"""
Some experiments on symbolic max plus matrices

We consider semigroup relations for various max-plus matrix semigroups:

 * upper triangular matrices with identical diagonals ('tri_sim_diag')
 * upper triangular matrices ('tri')
 * all matrices ('all')

If u = v is a relation then

pus = pvs is also a relation...

So we should look for primitive relations.


1. For band 2x2 matrices same diag (or triangular same diag) the relations are exactly

   p u s = p v s

   where p and s contains both x and y

Conjecture:

1. band matrices

    relations are of the form

    p u s = p v s

    with |p| >= dim and |s| >= dim

    (just use the fact that the first occurrences of factor of length d-1 appear
    at the same positions)

2. triangular matrices:

    relations are of the form

    p u s = p v s

    with |p| >= dim+1 and |s| >= dim+1

https://www.inf.ethz.ch/personal/fukudak/polyfaq/node23.html

EXAMPLES::

    sage: x,y=symbolic_max_plus_matrices(2,2,'ppl')
    sage: x.convex_hull_engine()
    'ppl'
    sage: x,y=symbolic_max_plus_matrices(2,2,'PALP')
    sage: x.convex_hull_engine()
    'PALP'
    sage: x,y=symbolic_max_plus_matrices(2,2,'cdd')
    sage: x.convex_hull_engine()
    'cdd'
"""

from itertools import product
import time
from subprocess import Popen, PIPE

from sage.misc.misc import SAGE_TMP
from sage.misc.temporary_file import tmp_filename
from sage.misc.misc import SAGE_TMP

from sage.geometry.polyhedron.parent import Polyhedra
from sage.modules.free_module import FreeModule
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ

##########################################
# Main function helper to build matrices #
##########################################

def symbolic_max_plus_identity(d, nvar, ch=None):
    r"""
    Return the ``d x d`` identity matrices on ``nvar`` variables.

    EXAMPLES::

        sage: m = symbolic_max_plus_identity(2, 4)
        sage: print m
        [   0  -oo ]
        [ -oo    0 ]
        sage: m.num_vars()
        4
    """
    d = int(d)
    nvar = int(nvar)
    V = FreeModule(ZZ, nvar)
    e = ()
    zero = (V([0]*nvar),)

    data = [[zero if i == j else e for j in range(d)] for i in range(d)]
    return SymbolicMaxPlusMatrix(d, nvar, data, ch)

def symbolic_max_plus_matrices(d, n, ch=None, sym=True):
    r"""
    Return ``n`` independent symbolic matrices in dimension ``d``.

    INPUT:

    - ``d`` -- the dimension

    - ``n`` -- the number of matrices

    - ``ch`` -- the convex hull engine

    - ``sym`` -- if ``False`` use dense matrices (i.e. store all coefficients)
      and if ``True`` (the default) uses matrices that stores only two coefficients.

    TESTS:

    Test a relation for 2x2 matrices::

        sage: A,B = symbolic_max_plus_matrices(2,2)
        sage: U = A*B*B*B*B*A*A
        sage: V = A*A*B*B*B*B*A
        sage: L = U*A*A*B*B*V
        sage: R = U*B*B*A*A*V
        sage: L == R
        True

    Check the validity of the above computation with dense matrices::

        sage: C,D = symbolic_max_plus_matrices(2,2,sym=False)
        sage: U2 = C*D*D*D*D*C*C
        sage: V2 = C*C*D*D*D*D*C
        sage: U == U2 and V == V2
        True
        sage: L2 = U2*C*C*D*D*V2
        sage: R2 = U2*D*D*C*C*V2
        sage: L == L2 and R == R2
        True
        sage: L2 == R2
        True
    """
    d = int(d)
    n = int(n)
    nvar = n * d * d
    V = FreeModule(ZZ, nvar)
    z = [0]*nvar
    matrices = []

    if sym:
        for i in range(n):
            z[i*d*d] = 1
            diag = (V(z),)
            z[i*d*d] = 0

            z[i*d*d+1] = 1
            nondiag = (V(z),)
            z[i*d*d+1] = 0

            matrices.append(SymbolicSymmetricMaxPlusMatrix(d, n, diag, nondiag, ch))
    else:
        o = 0
        for i in range(n):
            mat = []
            for j in range(d):
                row = []
                for k in range(d):
                    z[o] = 1
                    o += 1
                    row.append((V(z),))
                    z[o-1] = 0
                mat.append(row)
            matrices.append(SymbolicMaxPlusMatrix(d, nvar, mat, ch))

    return matrices

def symbolic_max_plus_matrices_band(d, n,
        diag='v', surdiag='v', ch=None):
    r"""
    INPUT:

    - ``d`` -- dimension

    - ``n`` -- number of matrices

    - ``diag`` -- either 'zero' (for zero), 'same' (for same, i.e. equal
      in each matrix) or 'var' (i.e. independent in each matrices). Shortcuts
      'c', 'z', 'v' can also be used.

    - ``surdiag`` -- either ``'zero'`` or ``'same'`` or ``'var'``.

    - ``ch`` -- convex hull engine to use

    EXAMPLES:

    For 'zv' the relations are the pairs `(u,v)` so that `Subwords_{d-1}(u) =
    Subwords_{d-1}(v)`::

        sage: x,y = symbolic_max_plus_matrices_band(3, 2, 'z', 'v')
        sage: x*y*x == x*x*x*y*x*x
        True
        sage: x*y*x*y == x*y*y*x
        True
        sage: x*y*x*y == y*x*y*x
        True
        sage: x*y*x*y == x*x*y*y
        False

        sage: x,y = symbolic_max_plus_matrices_upper(3, 2, 'z', 'v')
        sage: x*y*x == x*x*x*y*x*x
        True
        sage: x*y*x*y == x*y*y*x
        True
        sage: x*y*x*y == y*x*y*x
        True
        sage: x*y*x*y == x*x*y*y
        False

    And for 'vc'? seems to be much more complicated.

    TESTS:

    Non interesting cases::

        sage: x,y = symbolic_max_plus_matrices_band(2, 2, 'z', 'z')
        sage: assert x == y
        sage: x,y = symbolic_max_plus_matrices_band(2, 2, 'z', 's')
        sage: assert x == y
        sage: x,y = symbolic_max_plus_matrices_band(2, 2, 's', 'z')
        sage: assert y == x
        sage: x,y = symbolic_max_plus_matrices_band(2, 2, 's', 's')
        sage: assert y == y
    """
    assert diag in ('z', 's', 'v', 'zero', 'same', 'var')
    assert surdiag in ('z', 's', 'v', 'zero', 'same', 'var')
    diag = diag[0]
    surdiag = surdiag[0]

    d = int(d)
    n = int(n)
    assert d > 0 and n > 0

    nvar = 0
    if diag == 's':
        nvar += d
    elif diag == 'v':
        nvar += n*d
    if surdiag == 's':
        nvar += d-1
    elif surdiag == 'v':
        nvar += n*(d-1)

    V = FreeModule(ZZ, nvar)
    B = iter(V.basis())
    e = ()
    f = (V.zero(),)
    mat_init = [[e]*d for _ in range(d)]

    if diag == 'z':
        for k in range(d):
            mat_init[k][k] = f
    elif diag == 's':
        for k in range(d):
            mat_init[k][k] = (next(B),)

    if surdiag == 'z':
        for k in range(d-1):
            mat_init[k][k+1] = f
    elif surdiag == 's':
        for k in range(d-1):
            mat_init[k][k+1] = (next(B),)

    matrices = []
    for i in range(n):
        mat = [row[:] for row in mat_init]

        if diag == 'v':
            for k in range(d):
                mat[k][k] = (next(B),)
        if surdiag == 'v':
            for k in range(d-1):
                mat[k][k+1] = (next(B),)

        matrices.append(SymbolicMaxPlusMatrix(d, nvar, mat, ch))

    return matrices

def symbolic_max_plus_matrices_upper(d, n, 
        diag='v', surdiag='v', ch=None):
    r"""
    EXAMPLES::

        sage: x,y = symbolic_max_plus_matrices_upper(3,2)
        sage: x
        A 3x3 symbolic max plus matrix on 12 variables
        sage: print x
        [  x0   x3  x4 ]
        [ -oo   x1  x5 ]
        [ -oo  -oo  x2 ]
        sage: print y
        [  x6   x9  x10 ]
        [ -oo   x7  x11 ]
        [ -oo  -oo   x8 ]

        sage: x,y = symbolic_max_plus_matrices_upper(3,2,'z','v')
        sage: print x
        [   0   x0  x1 ]
        [ -oo    0  x2 ]
        [ -oo  -oo   0 ]
        sage: print y
        [   0   x3  x4 ]
        [ -oo    0  x5 ]
        [ -oo  -oo   0 ]

        sage: x,y = symbolic_max_plus_matrices_upper(3,2,'v','s')
        sage: print x
        [  x3   x0  x1 ]
        [ -oo   x4  x2 ]
        [ -oo  -oo  x5 ]
        sage: print y
        [  x6   x0  x1 ]
        [ -oo   x7  x2 ]
        [ -oo  -oo  x8 ]

    TESTS:

    Trivial cases::

        sage: x,y = symbolic_max_plus_matrices_upper(2, 2, 'z', 'z')
        sage: assert x == y
        sage: x,y = symbolic_max_plus_matrices_upper(2, 2, 'z', 's')
        sage: assert x == y
        sage: x,y = symbolic_max_plus_matrices_upper(2, 2, 's', 'z')
        sage: assert y == x
        sage: x,y = symbolic_max_plus_matrices_upper(2, 2, 's', 's')
        sage: assert y == y
    """
    assert diag in ('z', 's', 'v', 'zero', 'same', 'var')
    assert surdiag in ('z', 's', 'v', 'zero', 'same', 'var')
    diag = diag[0]
    surdiag = surdiag[0]

    d = int(d)
    n = int(n)

    nvar = 0
    if diag == 's':
        nvar += d
    elif diag == 'v':
        nvar += n*d

    if surdiag == 's':
        nvar += d*(d-1)//2
    elif surdiag == 'v':
        nvar += n*d*(d-1)//2

    V = FreeModule(ZZ, nvar)
    B = iter(V.basis())
    e = ()
    f = (V.zero(),)
    mat_init = [[e]*d for _ in range(d)]

    if diag == 'z':
        for k in range(d):
            mat_init[k][k] = f
    elif diag == 's':
        for k in range(d):
            mat_init[k][k] = (next(B),)

    if surdiag == 'z':
        for k1 in range(d):
            for k2 in range(k1+1,d):
                mat_init[k1][k2] = f
    elif surdiag == 's':
        for k1 in range(d):
            for k2 in range(k1+1,d):
                mat_init[k1][k2] = (next(B),)

    matrices = []
    for i in range(n):
        mat = [row[:] for row in mat_init]

        if diag == 'v':
            for k in range(d):
                mat[k][k] = (next(B),)
        if surdiag == 'v':
            for k1 in range(d):
                for k2 in range(k1+1,d):
                    mat[k1][k2] = (next(B),)

        matrices.append(SymbolicMaxPlusMatrix(d, nvar, mat, ch))

    return matrices

###########################
# helper for pretty print #
###########################

def pretty_relation_string(r1, r2, letters='mn'):
    r"""
    A nice string to represent a relation

    EXAMPLES::

        sage: pretty_relation_string([0,1,1,0],[1,1,0,0], 'xy')
        '(xyy)x = (yyx)x'
        sage: pretty_relation_string([0,0,1,0],[0,1,0,0], 'xy')
        'x(xy)x = x(yx)x'
    """
    p = 0
    while r1[p] == r2[p]: p += 1
    s = -1
    while r1[s] == r2[s]: s -= 1
    s += 1
    return '{}({}){} = {}({}){}'.format(
       ''.join(letters[i] for i in r1[:p]),
       ''.join(letters[i] for i in r1[p:s]),
       ''.join(letters[i] for i in r1[s:]),
       ''.join(letters[i] for i in r2[:p]),
       ''.join(letters[i] for i in r2[p:s]),
       ''.join(letters[i] for i in r2[s:]))

def pretty_print_poly(p):
    if not p:
        return '-oo'
    elif len(p) == 1:
        if p[0].is_zero():
            return '0'
        else:
            v = p[0]
            return '+'.join('{}x{}'.format('' if j == 1 else j,i) for i,j in enumerate(v) if j)
            
    else:
        s = []
        for v in p:
            s.append('+'.join('{}x{}'.format('' if j == 1 else j, i) for i,j in enumerate(v) if j))
        return 'max(' + ', '.join(s) + ')'

##################
# Matrix classes #
##################

class SymbolicMaxPlusMatrix(SageObject):
    r"""
    A symbolic max plus matrix.

    EXAMPLES::

        sage: M1 = symbolic_max_plus_identity(2,2)
        sage: M2 = M1 * M1
        sage: M2 == M1
        True

        sage: x1,y1,z1 = symbolic_max_plus_matrices_band(3,3,'s','v','ppl')
        sage: x2,y2,z2 = symbolic_max_plus_matrices_band(3,3,'s','v','cdd')
        sage: x3,y3,z3 = symbolic_max_plus_matrices_band(3,3,'s','v','PALP')
        sage: (x1*y1)*x1 == x1*(y1*x1) == (x2*y2)*x2 == x2*(y2*x2) == (x3*y3)*x3 == x3*(y3*x3)
        True
        sage: (y1*x1)*y1*x1 == y1*(x1*y1)*x1 == (y1*x1)*(y1*x1)
        True
        sage: (y2*x2)*y2*x2 == y2*(x2*y2)*x2 == (y2*x2)*(y2*x2)
        True
        sage: (y3*x3)*y3*x3 == y3*(x3*y3)*x3 == (y3*x3)*(y3*x3)
        True
        sage: x1*z1*x1*y1*y1*x1*z1 == x2*z2*x2*y2*y2*x2*z2 == x3*z3*x3*y3*y3*x3*z3
        True

        sage: x,y = symbolic_max_plus_matrices_band(2,2)
        sage: (x*y)*x == x*(y*x)
        True
        sage: A = x*y*y*x
        sage: A * x * y * A == A * y * x * A
        True
    """
    def __init__(self, d, nvars, data, ch=None):
        r"""
        INPUT:

        - ``d`` -- dimension

        - ``nvars` -- number of variabes

        - ``data`` -- the polyhedron given as a tuple of tuples of tuples

        - ``ch`` -- the convex hull engine to use (could be one of
          'ppl', 'cdd' or 'PALP')
        """
        self._d = d = int(d)
        self._nvars = nvars = int(nvars)
        if len(data) != d or any(len(x) != d for x in data):
            raise ValueError
        self._data = tuple(tuple(x) for x in data)

        self.convex_hull = get_convex_hull_engine(self._nvars, ch)

    def convex_hull_engine(self):
        r"""
        Return a string that describes the convex hull engine.

        EXAMPLES::

            sage: x,y = symbolic_max_plus_matrices_band(4, 2, ch='ppl')
            sage: x.convex_hull_engine()
            'ppl'

            sage: x,y = symbolic_max_plus_matrices_band(4, 2, ch='cdd')
            sage: x.convex_hull_engine()
            'cdd'

            sage: x,y = symbolic_max_plus_matrices_band(4, 2, ch='PALP') 
            sage: x.convex_hull_engine()
            'PALP'
        """
        return self.convex_hull._name

    def dim(self):
        r"""
        Return the dimension of this matrix.
        """
        return self._d

    def num_vars(self):
        r"""
        Return the number of variables used in this matrix.
        """
        return self._nvars

    def __getitem__(self, data):
        r"""
        To obtain an item of the matrix
        """
        i,j = data
        return self._data[i][j]

    # the method below are generic

    def __mul__(self, other):
        r"""
        Multiplication of matrices
        """
        assert self.dim() == other.dim() and self.num_vars() == other.num_vars()
        d = self.dim()
        new_data = [[None]*d for _ in range(d)]
        for i in range(d):
            for j in range(d):
                vertices = set()
                for k in range(d):
                    l = [x+y for x in self[i,k] for y in other[k,j]]
                    for v in l: v.set_immutable()
                    vertices.update(l)
                new_data[i][j] = self.convex_hull(vertices)

        return SymbolicMaxPlusMatrix(d, self._nvars, new_data, self.convex_hull)

    def equal_coefficients(self, other):
        r"""
        Return the list of equal coefficients between self and other.
        """
        d = self._d
        return [(i,j) for i in range(d) for j in range(d) \
                if self[i][j] == other[i][j]]

    def __eq__(self, other):
        r"""
        Equality test

        TESTS::

            sage: A,B = symbolic_max_plus_matrices(2,2,sym=True)
            sage: C,D = symbolic_max_plus_matrices(2,2,sym=False)
            sage: A == C
            True
            sage: A == A
            True
            sage: C == C
            True
            sage: A == B
            False
            sage: C == D
            False
        """
        if not isinstance(self, SymbolicMaxPlusMatrix) or \
           not isinstance(other, SymbolicMaxPlusMatrix):
            raise TypeError("can not compare {} with {}".format(type(self), type(other)))

        if self._nvars != other._nvars or self._d != other._d:
            return False

        return all(self[i,j] == other[i,j] for i in range(self._d) for j in range(self._d))

    def __ne__(self, other):
        r"""
        Difference test

        TESTS::

            sage: A,B = symbolic_max_plus_matrices(2,2,sym=True)
            sage: C,D = symbolic_max_plus_matrices(2,2,sym=False)
            sage: A != C
            False
            sage: A != A
            False
            sage: C != C
            False
            sage: A != B
            True
            sage: C != D
            True
        """
        if not isinstance(self, SymbolicMaxPlusMatrix) or \
           not isinstance(other, SymbolicMaxPlusMatrix):
            raise TypeError("can not compare {} with {}".format(type(self), type(other)))

        if self._nvars != other._nvars or self._d != other._d:
            return True

        return any(self[i,j] != other[i,j] for i in range(self._d) for j in range(self._d))

    def _repr_(self):
        r"""
        String when the object is printed
        """
        return 'A {}x{} symbolic max plus matrix on {} variables'.format(
                  self.dim(), self.dim(), self.num_vars())

    def __str__(self):
        str_data = []
        d = self.dim()
        for i in range(d):
            str_data.append([pretty_print_poly(self[i,j]) for j in range(d)])
        col_sizes = [max(len(str_data[i][j]) for i in range(d)) for j in range(d)]
        for i in range(d):
            str_data[i] = '[ ' + \
                '  '.join('{data:>{width}}'.format(data=data, width=width) for
                         data,width in zip(str_data[i], col_sizes)) + \
                          ' ]'
        return '\n'.join(str_data)

    str = __str__

    def eval(self, p):
        r"""
        Evaluates this symbolic matrix at the integer point ``p``.

        INPUT:

        - ``p`` - a vector whose length matches the number of variables of this
          matrix

        EXAMPLES::

            sage: x,y = symbolic_max_plus_matrices_band(2, 2, 'v', 'v')
            sage: v = (1,2,3,-4,5,6)
            sage: xv = x.eval(v)   # not tested
            sage: xv               # not tested
            [   1 3 ]
            [ -oo 2 ]
            sage: yv = y.eval(v)   # not tested
            sage: yv               # not tested
            [  -4 6 ]
            [ -oo 5 ]

            sage: (x*y*y*x*x).eval(v) == (xv*yv*yv*xv*xv)  # not tested
            True
        """
        F = FreeModule(ZZ, self._nvars)
        p = F(p)
        mat = []
        d = self.dim()
        for i in range(d):
            row = []
            for j in range(d):
                pts = self[i,j]
                row.append(minus_infinity() if not pts else max(p.dot_product(v) for v in pts))
            mat.append(row)
        return IntegerMaxPlusMatrix(self._d, mat) 

def vertex_swap(d, n, l, i1, i2, j1, j2):
    r"""
    Permutations (i1,j1) -> (i2,j2)
    """
    if i1 == i2 and j1 == j2:
        return l
    if i1 == j1:
        # (i1,i1) -> (i2,i2)
        assert i2 == j2
        def swap(v):
            swap2(d, n, v, i1, i2)
    elif i1 == i2:
        # (i,j1) -> (i,j2)
        def swap(v):
            swap2(d, n, v, j1, j2)
    elif j1 == j2:
        # (i1,j) -> (i2,j)
        def swap(v):
            swap2(d, n, v, i1, i2)
    elif i1 == j2 and i2 == j1:
        # (i1,j1) -> (j1,i1)
        def swap(v):
            swap2(d, n, v, i1, j1)
    elif i1 == j2:
        # (i1,j1) -> (i2,i1)
        def swap(v):
            swap3(d, n, v, j1, i1, i2)
    elif i2 == j1:
        # (i1,j1) -> (j1,j2)
        def swap(v):
            swap3(d, n, v, i1, j1, j2)
    else:
        # (i1,j1) -> (i2,j2)
        def swap(v):
            swap2(d, n, v, i1, i2)
            swap2(d, n, v, j1, j2)
    ll = []
    for v in l:
        v = v.__copy__()
        swap(v)
        v.set_immutable()
        ll.append(v)
    ll.sort()
    return tuple(ll)

def swap2(d, n, v, i, j):
    r"""
    Transposition (i,j)
    """
    for a in range(n):
        for k in range(d):
            if k == i or k == j:
                continue
            x = a*d*d + d*i + k
            y = a*d*d + d*j + k
            v[x], v[y] = v[y], v[x]

            x = a*d*d + d*k + i
            y = a*d*d + d*k + j
            v[x], v[y] = v[y], v[x]

        x = a*d*d + d*i + i
        y = a*d*d + d*j + j
        v[x], v[y] = v[y], v[x]

        x = a*d*d + d*j + i
        y = a*d*d + d*i + j
        v[x], v[y] = v[y], v[x]

def swap3(d, n, v, i, j, k):
    r"""
    permutation (i,j,k)
    """
    for a in range(n):
        for m in range(d):
            if m == i or m == j or m == k:
                continue
            x = a*d*d + d*i + m
            y = a*d*d + d*j + m
            z = a*d*d + d*k + m
            v[x],v[y],v[z] = v[z],v[x],v[y]

            x = a*d*d + d*m + i
            y = a*d*d + d*m + j
            z = a*d*d + d*m + k
            v[x],v[y],v[z] = v[z],v[x],v[y]

        x = a*d*d + d*i + i
        y = a*d*d + d*j + j
        z = a*d*d + d*k + k
        v[x],v[y],v[z] = v[z],v[x],v[y]

        x = a*d*d + d*i + j
        y = a*d*d + d*j + k
        z = a*d*d + d*k + i
        v[x],v[y],v[z] = v[z],v[x],v[y]

        x = a*d*d + d*i + k
        y = a*d*d + d*j + i
        z = a*d*d + d*k + j
        v[x],v[y],v[z] = v[z],v[x],v[y]

class SymbolicSymmetricMaxPlusMatrix(SymbolicMaxPlusMatrix):
    r"""
    Matrices that are symmetric under permutations.

    This class only concerns relations on the full matrix monoid. In practice,
    we should win a factor d^2/2 where d is the dimension.
    
    Such matrices are completely determined by the entries (0,0) and (0,1) (that
    are stored in the attributes ``self._diag`` and ``self._nondiag``).

    .. TODO::

        see whether the convex hull computation can be simplified using the many
        symmetries
    """
    def __init__(self, d, n, diag, nondiag, ch=None):
        self._diag = diag        # diagonal term (0,0)
        self._nondiag = nondiag  # nondiagonal term (0,1)
        self._d = d              # dimension
        self._n = n              # number of matrices
        self._nvars = n*d*d      # number of variables

        self.convex_hull = get_convex_hull_engine(self._nvars, ch)

    def __eq__(self, other):
        r"""
        TESTS::

            sage: x1,y1 = symbolic_max_plus_matrices(5, 2, sym=True)
            sage: x2,y2 = symbolic_max_plus_matrices(5, 2, sym=False)
            sage: x1 == x2 and x1 == x1
            True
            sage: x1 == y1 or x1 == y2
            False
        """
        if type(self) is type(other):
            return self._d == other._d and \
                   self._n == other._n and \
                   self._diag == other._diag and \
                   self._nondiag == other._nondiag
        else:
            return SymbolicMaxPlusMatrix.__eq__(self, other)

    def __ne__(self, other):
        r"""
        TESTS::

            sage: x1,y1 = symbolic_max_plus_matrices(5, 2, sym=True)
            sage: x2,y2 = symbolic_max_plus_matrices(5, 2, sym=False)
            sage: x1 != x2 or x1 != x1
            False
            sage: x1 != y1 and x1 != y2
            True
        """
        if type(self) is type(other):
            return self._d != other._d or \
                   self._n != other._n or \
                   self._diag != other._diag or \
                   self._nondiag != other._nondiag
        else:
            return SymbolicMaxPlusMatrix.__ne__(self, other)

    def __getitem__(self, data):
        r"""
        TESTS::

            sage: x,y = symbolic_max_plus_matrices(3, 2, sym=True)
            sage: lx = [x[i,j] for i in range(3) for j in range(3)]
            sage: ly = [y[i,j] for i in range(3) for j in range(3)]
            sage: lx + ly == [(v,) for v in (ZZ**18).basis()]
            True
        """
        i,j = data
        i = int(i); j = int(j)
        if i < 0 or i >= self._d or j < 0 or j >= self._d:
            raise IndexError("matrix index out of range")

        if i == j:
            return vertex_swap(self._d, self._n, self._diag, 0, i, 0, j)
        else:
            return vertex_swap(self._d, self._n, self._nondiag, 0, i, 1, j)

    def __mul__(self, other):
        r"""
        TESTS::

            sage: x1,y1 = symbolic_max_plus_matrices(5, 2, sym=True)
            sage: x2,y2 = symbolic_max_plus_matrices(5, 2, sym=False)
            sage: assert x1 == x2 and y1 == y2
            sage: assert x1*y1 == x2*y2 and y1*x1 == y2*x2
            sage: assert x1*x1 == x2*x2
            sage: assert x1*y1*x1 == x2*y2*x2
            sage: assert x1*x1*x1 == x2*x2*x2

            sage: x1,y1 = symbolic_max_plus_matrices(4, 2, sym=True)
            sage: x2,y2 = symbolic_max_plus_matrices(4, 2, sym=False)
            sage: x1*x1 == x1*x2 == x2*x1 == x2*x2
            True
        """
        if type(self) is not type(other):
            return SymbolicMaxPlusMatrix.__mul__(self, other)

        assert self._n == other._n and self._d == other._d 

        data = []
        for j in range(2):
            vertices = set()
            for k in range(self._d):
                l = [x+y for x in self[0,k] for y in other[k,j]]
                for v in l: v.set_immutable()
                vertices.update(l)
            data.append(self.convex_hull(vertices))

        return SymbolicSymmetricMaxPlusMatrix(self._d, self._n, data[0], data[1], self.convex_hull)

#######################
# convex hull engines #
#######################

def get_convex_hull_engine(nvar, convex_hull=None):
    r"""
    Call to various library for convex hull computation

    EXAMPLES::

        sage: CH1 = get_convex_hull_engine(3, 'ppl')
        sage: CH2 = get_convex_hull_engine(3, 'cdd')
        sage: CH3 = get_convex_hull_engine(3, 'PALP')

        sage: F = ZZ**3
        sage: pts = [F.random_element() for _ in range(20)]
        sage: CH1(pts) == CH2(pts) == CH3(pts)
        True

    By far the fastest is CH1. See also:

    - http://miconvexhull.codeplex.com/
    - http://www.boost.org/doc/libs/1_47_0/libs/geometry/doc/html/geometry/reference/algorithms/convex_hull.html
    """
    if isinstance(convex_hull, ConvexHull):
        return convex_hull
    elif convex_hull is None or convex_hull == 'ppl':
        return ConvexHullPolyhedra(nvar, 'ppl')
    elif convex_hull == 'cdd':
        return ConvexHullPolyhedra(nvar, 'cdd')
    elif convex_hull == 'PALP':
        return ConvexHullPalp(nvar)
    else:
        raise ValueError("convex_hull must either be 'ppl', 'cdd' or 'PALP'")

class ConvexHull(object):
    r"""
    Class for convex hull computation.

    to be implemented in subclasses:

    - ``__init__(self, dimension)``: initialize the class to work in dimension
      ``dim``
    - ``__call__(self, pts)``: return the convex hull of the list ``pts``
    """
    pass

class ConvexHullPolyhedra(ConvexHull):
    r"""
    EXAMPLES::

        sage: C = ConvexHullPolyhedra(2)
        sage: C([(0,1),(2,3),(1,2),(3,4),(0,0),(0,0),(1,1),(1,1)])
        ((0, 0), (0, 1), (1, 1), (3, 4))
    """
    _name = None

    def __init__(self, dim, backend=None):
        from sage.geometry.polyhedron.parent import Polyhedra
        self._parent = Polyhedra(ZZ, int(dim), backend)
        self._name = backend

    def __eq__(self, other):
        return type(self) is type(other) and self._parent == other._parent

    def __ne__(self, other):
        return type(self) is not type(other) or self._parent != other._parent

    def __call__(self, pts):
        if pts:
            pts = [p.vector() for p in self._parent([pts,[],[]],None).vertex_generator()]
            for a in pts: a.set_immutable()
            pts.sort()
        return tuple(pts)

class ConvexHullPalp(ConvexHull):
    r"""
    Compute convex hull using PALP.

    Note: the points must be lattice points (ie integer coordinates) and
    generate the ambient vector space.
    """
    _name = 'PALP'

    def __init__(self, dim):
        from sage.modules.free_module import FreeModule
        self._free_module = FreeModule(ZZ, int(dim))

    def __eq__(self, other):
        return type(self) is type(other) and self._free_module == other._free_module

    def __ne__(self, other):
        return type(self) is not type(other) or self._free_module != other._free_module

    def __call__(self, pts):
        filename = tmp_filename()

        pts = list(pts)
        n = len(pts)
        if n <= 2:
            return tuple(pts)
        d = len(pts[0])
        assert d == self._free_module.rank()

        # PALP only works with full dimension polyhedra!!
        ppts = [x-pts[0] for x in pts]
        U = self._free_module.submodule(ppts)
        d2 = U.rank()
        if d2 != d:
            # not full dim
            # we compute two matrices
            #  M1: small space -> big space  (makes decomposition)
            #  M2: big space -> small space  (i.e. basis of the module)
            # warning: matrices act on row vectors, i.e. left action
            from sage.matrix.constructor import matrix
            from sage.modules.free_module import FreeModule
            V2 = FreeModule(ZZ,d2)
            M1 = U.matrix()
            assert M1.nrows() == d2
            assert M1.ncols() == d
            M2 = matrix(QQ,d)
            M2[:d2,:] = M1
            i = d2
            U = U.change_ring(QQ)
            F = self._free_module.change_ring(QQ)
            for b in F.basis():
                if b not in U:
                    M2.set_row(i, b)
                    U = F.submodule(U.basis() + [b])
                    i += 1
            assert i == self._free_module.rank()
            M2 = (~M2)[:,:d2]
            assert M2.nrows() == d
            assert M2.ncols() == d2
            assert (M1*M2).is_one()
            pts2 = [p * M2 for p in ppts]
        else:
            pts2 = pts
            d2 = d

        with open(filename, "w") as output:
            output.write("{} {}\n".format(n,d2))
            for p in pts2:
                output.write(" ".join(map(str,p)))
                output.write("\n")

        args = ['poly.x', '-v', filename]
        try:
            palp_proc = Popen(args,
                    stdin=PIPE, stdout=PIPE,
                    stderr=None, cwd=str(SAGE_TMP))
        except OSError:
            raise RuntimeError("Problem with calling PALP")
        ans, err = palp_proc.communicate()
        ret_code = palp_proc.poll()

        if ret_code:
            raise RuntimeError("PALP return code is {} from input {}".format(ret_code, pts))
        ans = ans.split('\n')
        dd,nn = ans[0].split(' ')[:2]
        dd = int(dd)
        nn = int(nn)
        if dd > nn: dd,nn = nn,dd
        if d2 != int(dd):
            raise RuntimeError("dimension changed... have d={} but PALP answered dd={} and nn={}".format(d2,dd,nn))
        n2 = int(nn)
        coords = []
        for i in xrange(1,d2+1):
            coords.append(map(ZZ,ans[i].split()))
        new_pts = zip(*coords)

        if d2 != d:
            new_pts = [pts[0] + V2(p)*M1 for p in new_pts]

        t = [self._free_module(p) for p in new_pts]
        for v in t:
            v.set_immutable()
        t.sort()
        return tuple(t)
