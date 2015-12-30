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

from sage.geometry.polyhedron.parent import Polyhedra
from sage.modules.free_module import FreeModule
from sage.rings.integer_ring import ZZ

from semigroup_tools import (products_n, products_p,
        constraints_from_subwords, minimal_all_subwords_prefix,
        minimal_all_subwords_suffix)

##########################################
# Main function helper to build matrices #
##########################################

def symbolic_max_plus_identity(d, nvar,ch=None):
    r"""
    Return the ``d x d`` identity matrices on ``nvar`` variables.
    """
    d = int(d)
    nvar = int(nvar)
    V = FreeModule(ZZ, nvar)
    e = ()
    zero = (V([0]*nvar),)

    data = [[zero if i == j else e for j in range(d)] for i in range(d)]
    return SymbolicMaxPlusMatrix(d, nvar, data, ch)

def symbolic_max_plus_matrices(d, n, ch=None, sym=False):
    r"""
    Return ``n`` independent symbolic matrices in dimension ``d``.

    INPUT:

    - ``d`` -- the dimension

    - ``n`` -- the number of matrices

    - ``ch`` -- the convex hull engine

    - ``sym`` -- if ``False`` use dense matrices (i.e. store all coefficients)
      and if ``True`` uses matrices that store only two coefficients.
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

    - ``diag`` -- either 'const' (for constant), 'sim' (for similar, i.e. equal
      in each matrix) or 'var' (i.e. independent in each matrices). Shortcuts
      'c', 's', 'v' can also be used.

    - ``ch`` -- convex hull engine to use

    EXAMPLES:

    For 'cv' the relations are the pairs `(u,v)` so that `Subwords_{d-1}(u) =
    Subwords_{d-1}(v)`::


        sage: x,y = symbolic_max_plus_matrices_band(3, 2, 'c', 'v')
        sage: x*y*x == x*x*x*y*x*x
        True
        sage: x*y*x*y == x*y*y*x
        True
        sage: x*y*x*y == y*x*y*x
        True
        sage: x*y*x*y == x*x*y*y
        False

    For 'vc'?::

    TESTS:

    Non interesting cases::

        sage: x,y = symbolic_max_plus_matrices_band(2, 2, 'c', 'c')
        sage: assert x == y
        sage: x,y = symbolic_max_plus_matrices_band(2, 2, 'c', 's')
        sage: assert x == y
        sage: x,y = symbolic_max_plus_matrices_band(2, 2, 's', 'c')
        sage: assert y == x
        sage: x,y = symbolic_max_plus_matrices_band(2, 2, 's', 's')
        sage: assert y == y
    """
    assert diag in ('c', 's', 'v', 'const', 'sim', 'var')
    assert surdiag in ('c', 's', 'v', 'const', 'sim', 'var')
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

    if diag == 'c':
        for k in range(d):
            mat_init[k][k] = f
    elif diag == 's':
        for k in range(d):
            mat_init[k][k] = (next(B),)

    if surdiag == 'c':
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

    sage: x,y = symbolic_max_plus_matrices_upper(3,2,'c','v')
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
    """
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

    if diag == 'c':
        for k in range(d):
            mat_init[k][k] = f
    elif diag == 's':
        for k in range(d):
            mat_init[k][k] = (next(B),)

    if surdiag == 'c':
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

        sage: x1,y1,z1 = symbolic_max_plus_matrices_band_sim_diag(3,3,'ppl')
        sage: x2,y2,z2 = symbolic_max_plus_matrices_band_sim_diag(3,3,'cdd')
        sage: x3,y3,z3 = symbolic_max_plus_matrices_band_sim_diag(3,3,'PALP')
        sage: (x1*y1)*x1 == x1*(y1*x1) == (x2*y2)*x2 == x2*(y2*x2) == (x3*y3)*x3 == x3*(y3*x3)
        True
        sage: 
        sage: (y1*x1)*y1*x1 == y1*(x1*y1)*x1 == (y1*x1)*(y1*x1)
        True
        sage: (y2*x2)*y2*x2 == y2*(x2*y2)*x2 == (y2*x2)*(y2*x2)
        True
        sage: (y3*x3)*y3*x3 == y3*(x3*y3)*x3 == (y3*x3)*(y3*x3)
        True
        sage: x1*z1*x1*y1*y1*x1*z1 == x2*z2*x2*y2*y2*x2*z2 == x3*z3*x3*y3*y3*x3*z3
        True

        sage: x,y = symbolic_max_plus_matrices_band_sim_diag(2,2)
        sage:(x*y)*x == x*(y*x)
        True
        sage: A = x*y*y*x
        sage: A * x * y * A == A * y * x * A
        True
    """
    def __init__(self, n, nvar, data, ch=None):
        r"""
        INPUT:

        - ``n`` -- dimension

        - ``nvar`` -- number of variabes

        - ``data`` -- the polyhedron given as a tuple of tuples of tuples

        - ``ch`` -- the convex hull engine to use (could be one of
          'ppl', 'cdd' or 'PALP')
        """
        self._n = n = int(n)
        self._nvar = nvar = int(nvar)
        if len(data) != n or any(len(x) != n for x in data):
            raise ValueError
        self._data = tuple(tuple(x) for x in data)

        import convex_hull
        self.convex_hull = convex_hull.get_convex_hull_engine(self._nvar, ch)

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

    def num_variables(self):
        r"""
        Return the number of variables used in this matrix.
        """
        return self._nvar

    def equal_coefficients(self, other):
        r"""
        Return the list of equal coefficients between self and other.
        """
        n = self._n
        return [(i,j) for i in range(n) for j in range(n) \
                if self._data[i][j] == other._data[i][j]]

    def __eq__(self, other):
        r"""
        Equality test
        """
        if type(self) is type(other):
            return self._nvar == other._nvar and self._data == other._data
        else:
            raise TypeError("can not compare {} with {}".format(type(self), type(other)))

    def __ne__(self, other):
        r"""
        Difference test
        """
        if type(self) is type(other):
            return self._data != other._data or self._data != other._data
        else:
            raise TypeError("can not compare {} with {}".format(type(self), type(other)))

    def __mul__(self, other):
        r"""
        Multiplication of matrices
        """
        n = self._n
        assert self._n == other._n and self._nvar == other._nvar and self.convex_hull == other.convex_hull
        new_data = [[None]*n for _ in range(n)]
        for i in range(n):
            for j in range(n):
                vertices = set()
                for k in range(n):
                    l = [x+y for x in self._data[i][k] for y in other._data[k][j]]
                    for v in l: v.set_immutable()
                    vertices.update(l)
                new_data[i][j] = self.convex_hull(vertices)

        return SymbolicMaxPlusMatrix(self._n, self._nvar, new_data, self.convex_hull)

    def _repr_(self):
        r"""
        String when the object is printed
        """
        return 'A {}x{} symbolic max plus matrix on {} variables'.format(self._n, self._n, self._nvar)

    def __str__(self):
        str_data = []
        for i in range(self._n):
            str_data.append([pretty_print_poly(self._data[i][j]) for j in range(self._n)])
        col_sizes = [max(len(str_data[i][j]) for i in range(self._n)) for j in range(self._n)]
        for i in range(self._n):
            str_data[i] = '[ ' + \
                '  '.join('{data:>{width}}'.format(data=data, width=width) for
                         data,width in zip(str_data[i], col_sizes)) + \
                          ' ]'
        return '\n'.join(str_data)

    str = __str__

    def __getitem__(self, data):
        r"""
        To obtain an item of the matrix
        """
        i,j = data
        return self._data[i][j]

    def eval(self, p):
        r"""
        Evaluates this symbolic matrix at the integer point ``p``.

        INPUT:

        - ``p`` - a vector whose length matches the number of variables of this
          matrix

        EXAMPLES::

            sage: x,y = symbolic_max_plus_matrices_band(2, 2, 'v', 'v')
            sage: v = (1,2,3,-4,5,6)
            sage: xv = x.eval(v)
            sage: xv
            [   1 3 ]
            [ -oo 2 ]
            sage: yv = y.eval(v)
            sage: yv
            [  -4 6 ]
            [ -oo 5 ]

            sage: (x*y*y*x*x).eval(v) == (xv*yv*yv*xv*xv)
            True
        """
        F = FreeModule(ZZ,self._nvar)
        p = F(p)
        mat = []
        for i in range(self._n):
            row = []
            for j in range(self._n):
                pts = self._data[i][j]
                row.append(minus_infinity() if not pts else max(p.dot_product(v) for v in pts))
            mat.append(row)
        return IntegerMaxPlusMatrix(self._n, mat) 

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

class SymbolicSymmetricMaxPlusMatrix(SageObject):
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
        self._d = d  # dimension
        self._n = n  # number of matrices

        import convex_hull
        self.convex_hull = convex_hull.get_convex_hull_engine(self._d * self._d * self._n, ch)

    def __repr__(self):
        return 'Symmetric {}x{} max plus matrix:\n diagonal: {}\n non-diag: {}'.format(
                self._d,self._d,
                pretty_print_poly(self._diag),
                pretty_print_poly(self._nondiag))

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
            sage: mat_eq = lambda m1,m2: all(m1[i,j] == m2[i,j] for i in range(5) for j in range(5))
            sage: assert mat_eq(x1,x2)
            sage: assert mat_eq(y1,y2)
            sage: assert mat_eq(x1*y1, x2*y2)
            sage: assert mat_eq(y1*x1, y2*x2)
            sage: assert mat_eq(x1*x1, x2*x2)
            sage: assert mat_eq(x1*y1*x1, x2*y2*x2)
            sage: assert mat_eq(x1*x1*x1, x2*x2*x2)
        """
        assert type(self) is type(other) and \
               self._n == other._n and \
               self._d == other._d and \
               self.convex_hull == other.convex_hull

        data = []
        for j in range(2):
            vertices = set()
            for k in range(self._d):
                l = [x+y for x in self[0,k] for y in other[k,j]]
                for v in l: v.set_immutable()
                vertices.update(l)
            data.append(self.convex_hull(vertices))

        return SymbolicSymmetricMaxPlusMatrix(self._d, self._n, data[0], data[1], self.convex_hull)

######################
# Experimental stuff #
######################
def prod_symbolic(w, mats):
    r"""
    Recursive computation that tries to minize the number of products

    DO NOT USE... IT DOES NOT SEEMS FASTER!
    """
    assert max(w) <= len(mats)
    if len(w) <= 3 or 3*len(mats) >= 2*len(w):
        return prod(mats[i] for i in w)

    ww = [(w[i],w[i+1]) for i in range(0,len(w)-1,2)]
    P2 = list(set(ww))
    P2_index = {p:i for i,p in enumerate(P2)}
    mmats = [mats[p[0]] * mats[p[1]] for p in P2]

    ww = [P2_index[u] for u in ww]
    if len(w)%2:
        mmats.append(mats[w[-1]])
        ww.append(len(mmats)-1)
    return prod_symbolic(ww, mmats)

def relations_band(dim, start=3, num_mat=5000, limit_coeff=1073741824, filename=None):
    r"""
    List the relations for the band matrices (with identical diagonals) in a
    given dimension.

    INPUT:

    - ``dim`` -- the dimension

    - ``start`` -- the length of product we consider first (default to ``1``)

    - ``num_mat`` -- the number of integer matrices used to check the relation.
      Note that if all the integer matrices satisfy the relations, a (costly)
      symbolic check is performed to guarantee that the relation is satisfied by
      any pair of matrices.

    .. NOTE::

        All the time in spent in symbolic verification. The only way to make it
        fast is to have a huge ``num_mat`` in order to eliminate the maximum of
        relations without symbolic computations.

    EXAMPLES::

        sage: relations_band(3)
        # band similar diagonal relations for n = 11
        n = 11
        xxyyx(y)xxyxy = xxyyx(x)xxyxy
        xxyyx(y)xxyyx = xxyyx(x)xxyyx
        xxyyx(y)xyyxx = xxyyx(x)xyyxx
        xxyyx(y)yxxyy = xxyyx(x)yxxyy
        xxyyx(y)yyxxy = xxyyx(x)yyxxy
        xyxyy(y)yxxyy = xyxyy(x)yxxyy
        xyxyy(y)yyxxy = xyxyy(x)yyxxy
        xyxyy(y)yyxyx = xyxyy(x)yyxyx
        xyyxx(y)xxyxy = xyyxx(x)xxyxy
        xyyxx(y)xxyyx = xyyxx(x)xxyyx
        xyyxx(y)xyyxx = xyyxx(x)xyyxx
        xyyxx(y)yxxyy = xyyxx(x)yxxyy
        xyyxx(y)yyxxy = xyyxx(x)yyxxy
        # band similar diagonal relations for n = 12
        n = 12
        xxxyyx(y)xxyxy = xxxyyx(x)xxyxy
        xxxyyx(y)xxyyx = xxxyyx(x)xxyyx
        xxxyyx(y)xyyxx = xxxyyx(x)xyyxx
        xxxyyx(y)yxxyy = xxxyyx(x)yxxyy
        xxxyyx(y)yyxxy = xxxyyx(x)yyxxy
        xxyxyx(y)xxyxy = xxyxyx(x)xxyxy
        xxyxyx(y)xxyyx = xxyxyx(x)xxyyx
        xxyxyx(y)xyyxx = xxyxyx(x)xyyxx
        xxyxyy(y)xxyyx = xxyxyy(x)xxyyx
        xxyxyy(y)xyyxx = xxyxyy(x)xyyxx
        xxyxyy(y)yxxyy = xxyxyy(x)yxxyy
        xxyxyy(y)yyxxy = xxyxyy(x)yyxxy
        xxyxyy(y)yyxyx = xxyxyy(x)yyxyx
        xxyyxx(y)xxyxy = xxyyxx(x)xxyxy
        xxyyxx(y)xxyyx = xxyyxx(x)xxyyx
        xxyyxx(y)xyyxx = xxyyxx(x)xyyxx
        xxyyxx(y)yxxyy = xxyyxx(x)yxxyy
        xxyyxx(y)yyxxy = xxyyxx(x)yyxxy
        xxyyx(y)xxxyxy = xxyyx(x)xxxyxy
        xxyyx(yx)xxyxy = xxyyx(xy)xxyxy
        xxyyx(y)xxxyyx = xxyyx(x)xxxyyx
        xxyyx(yx)xxyyx = xxyyx(xy)xxyyx
        xxyyx(y)xxyxxy = xxyyx(x)xxyxxy
        xxyyx(y)xxyxyx = xxyyx(x)xxyxyx
        xxyyx(y)xxyxyy = xxyyx(x)xxyxyy
        xxyyx(y)xxyyxx = xxyyx(x)xxyyxx
        xxyyx(y)xxyyxy = xxyyx(x)xxyyxy
        xxyyx(yx)xyyxx = xxyyx(xy)xyyxx
        xxyyx(y)xxyyyx = xxyyx(x)xxyyyx
        xxyyx(y)xyxxyy = xxyyx(x)xyxxyy
        xxyyx(yx)yxxyy = xxyyx(xy)yxxyy
        xxyyx(y)xyxyxx = xxyyx(x)xyxyxx
        xxyyx(y)xyyxxx = xxyyx(x)xyyxxx
        xxyyx(y)xyyxxy = xxyyx(x)xyyxxy
        xxyyx(yx)yyxxy = xxyyx(xy)yyxxy
        xxyyx(y)xyyyxx = xxyyx(x)xyyyxx
        xxyyx(y)yxxxyy = xxyyx(x)yxxxyy
        xxyyx(yy)xxyyx = xxyyx(xx)xxyyx
        xxyyx(y)yxxyyx = xxyyx(x)yxxyyx
        xxyyx(y)yxxyyy = xxyyx(x)yxxyyy
        xxyyxy(y)xxyyx = xxyyxy(x)xxyyx
        xxyyx(yy)xyyxx = xxyyx(xx)xyyxx
        xxyyx(y)yxyyxx = xxyyx(x)yxyyxx
        xxyyxy(y)xyyxx = xxyyxy(x)xyyxx
        xxyyx(y)yyxxxy = xxyyx(x)yyxxxy
        xxyyx(yy)yxxyy = xxyyx(xx)yxxyy
        xxyyx(y)yyxxyx = xxyyx(x)yyxxyx
        xxyyx(y)yyxxyy = xxyyx(x)yyxxyy
        xxyyxy(y)yxxyy = xxyyxy(x)yxxyy
        ...
    """
    # the symbolic max-plus matrices
    ab = symbolic_max_plus_matrices_band(dim, 2, diag='s', surdiag='v')

    # the integer max-plus matrices
    pairs = [random_integer_max_plus_matrices_band(dim, -limit_coeff, limit_coeff) for _ in range(num_mat)]
    one_int = integer_max_plus_matrix_identity(dim)

    # n     : length of the product
    # ii1   : first word
    # ii2   : second word
    n = start
    while True:
        eliminated_from_int = 0
        eliminated_from_symb = 0
        relations = 0

        if filename is None:
            from sys import stdout
            f = stdout
        else:
            f = open(filename.format(dim,n), 'w')

        # we restrict to products that start and ends with the same letters
        # up to symmetry we can restrict to
        #   X ... X = X ... X
        #   X ... Y = X ... Y

        f.write("# band similar diagonal relations for n = {}\n".format(n))
        print "n = {}".format(n)
        for i1 in product((0,1), repeat=n-2):
            for pref,suff in ((0,),(0,)),((0,),(1,)):
                ii1 = pref + i1 + suff
                assert len(ii1) == n
                p,s = constraints_from_subwords(ii1, dim)
                if p != -1:
                    i1p = ii1[:p]; i1s = ii1[s:]
                    for i2 in product((0,1), repeat=s-p):
                        ii2 = i1p + i2 + i1s
                        assert len(ii2) == n
                        if ii1 == ii2:
                            break
                        if not is_relation(ii1, ii2, pairs):
                            eliminated_from_int += 1
                            continue
                        if prod(ab[x] for x in ii1) != prod(ab[x] for x in ii2):
                            eliminated_from_symb += 1
                            continue

                        f.write(pretty_relation_string(ii1,ii2,'xy'))
                        f.write('\n')
                        f.flush()
                        relations += 1
                        del i2  # makes itertools faster!
            del i1  # makes itertools faster!

        if filename is not None:
            f.close()
        n += 1

        print "  int elimination  : {}".format(eliminated_from_int)
        print "  symb elimination : {}".format(eliminated_from_symb)
        print "  relations        : {}".format(relations)

def relations_tri_sim_diag(dim, start=3, num_mat=20, filename=None):
    r"""
    List the relations for the upper triangular matrices in a given dimension.

    INPUT:

    - ``dim`` -- the dimension

    - ``start`` -- the length of product we consider first (default to ``1``)

    - ``num_mat`` -- the number of integer matrices used to check the relation.
      Note that if all the integer matrices satisfy the relations, a (costly)
      symbolic check is performed to guarantee that the relation is satisfied by
      any pair of matrices.

    EXAMPLES::

        sage: relations_tri(2)
        # triangular relations for n = 1
        # triangular relations for n = 2
        # triangular relations for n = 3
        # triangular relations for n = 4
        # triangular relations for n = 5
        # triangular relations for n = 6
        # triangular relations for n = 7
        # triangular relations for n = 8
        # triangular relations for n = 9
        # triangular relations for n = 10
        #  relations (5,5)
        xyyx(xy)xyyx = xyyx(yx)xyyx
        xyyx(xy)yxxy = xyyx(yx)yxxy
        xyyx(yx)xyyx = xyyx(xy)xyyx
        xyyx(yx)yxxy = xyyx(xy)yxxy
        # triangular relations for n = 11
        #  relations (4,7)
        xyyyyx(xy)yxy = xyyyyx(yx)yxy
        xyyyyx(yx)yxy = xyyyyx(xy)yxy
        ...
    """
    # the symbolic max-plus matrices
    ab = symbolic_max_plus_matrices_tri_sim_diag(dim, 2)
    abc = symbolic_max_plus_matrices_tri_sim_diag(dim, 3)
    one_ab = symbolic_max_plus_identity(dim, ab[0].num_variables())
    one_abc = symbolic_max_plus_identity(dim, abc[0].num_variables())


    # the integer max-plus matrices
    pairs = [random_integer_max_plus_matrices_tri_sim_diag(dim, -2**30, 2**30) for _ in range(num_mat)]
    A,B = AB = pairs[0]

    one_int = integer_max_plus_matrix_identity(dim)

    # n     : length of the product
    # i1,m1AA,m1AB : data for the first word
    # i2,m2AA,m2AB : data for the second word
    n = start
    while True:
        if filename is None:
            from sys import stdout
            f = stdout
        else:
            f = open(filename.format(dim,n), 'w')

        f.write("# triangular similar diagonal relations for n = {}\n".format(n))
        relations = []

        # we first want to consider the prefix/suffix such that we got a
        # relation whatever is in between (of the same length)
        # how do we identify these relations?

        for p in minimal_all_subwords_prefix(n - 2*(dim-1), dim-1):
            m_prefix = prod(AB[i] for i in p)
            for s in minimal_all_subwords_suffix(n - len(p) - 1, dim-1):
                m_suffix = prod(AB[i] for i in s)
                for i1,m1 in products_n(A,B,n-len(p)-len(s)):
                    for i2,m2 in products_n(A,B,n-len(p)-len(s)):
                        if i1 == i2:
                            break
#
#                # here we first test equality between m1 and m2
#                # then we test the relations on all matrices in mats
#                # then we check formally using symbolic matrices
#                ii1 = [0] + i1 + [0]
#                ii2 = [0] + i1[:ppos] + i2 + i1[spos:] + [0]
#                print "{} =?= {}".format(ii1,ii2)
#                if m1AA == m2AA and \
#                   is_relation(ii1, ii2, pairs) and \
#                   prod(ab[x] for x in ii1) == prod(ab[x] for x in ii2):
#                       f.write(pretty_relation_string(ii1,ii2,'xy'))
#                       f.write('\n')
#                       f.flush()
#                ii1 = [0] + i1 + [1]
#                ii2 = [0] + i1[:ppos] + i2 + i1[spos:] + [0]
#                if m1AB == m2AB and \
#                   is_relation(ii1, ii2, pairs) and \
#                   prod(ab[x] for x in ii1) == prod(ab[x] for x in ii2):
#                       f.write(pretty_relation_string(ii1,ii2,'xy'))
#                       f.write('\n')
#                       f.flush()
#
#        if filename is not None:
#            f.close()

        n += 1

def relations_tri(dim, start=1, num_mat=10, filename=None):
    r"""
    List the relations for the upper triangular matrices in a given dimension.

    INPUT:

    - ``dim`` -- the dimension

    - ``start`` -- the length of product we consider first (default to ``1``)

    - ``num_mat`` -- the number of integer matrices used to check the relation.
      Note that if all the integer matrices satisfy the relations, a (costly)
      symbolic check is performed to guarantee that the relation is satisfied by
      any pair of matrices.

    EXAMPLES::

        sage: relations_tri(2)
        # triangular relations for n = 1
        # triangular relations for n = 2
        # triangular relations for n = 3
        # triangular relations for n = 4
        # triangular relations for n = 5
        # triangular relations for n = 6
        # triangular relations for n = 7
        # triangular relations for n = 8
        # triangular relations for n = 9
        # triangular relations for n = 10
        #  relations (5,5)
        xyyx(xy)xyyx = xyyx(yx)xyyx
        xyyx(xy)yxxy = xyyx(yx)yxxy
        xyyx(yx)xyyx = xyyx(xy)xyyx
        xyyx(yx)yxxy = xyyx(xy)yxxy
        # triangular relations for n = 11
        #  relations (4,7)
        xyyyyx(xy)yxy = xyyyyx(yx)yxy
        xyyyyx(yx)yxy = xyyyyx(xy)yxy
        ...
    """
    # the symbolic max-plus matrices
    ab = symbolic_max_plus_matrices_tri(dim, 2)
    one = symbolic_max_plus_identity(dim, ab[0].num_variables())

    # the integer max-plus matrices
    mats = [random_integer_max_plus_matrix_tri(dim, -50*dim*dim, 50*dim*dim) for _ in range(num_mat)]
    pairs = [(m1,m2) for m1 in mats for m2 in mats if m1 is not m2]
    one_int = integer_max_plus_matrix_identity(dim)

    # n     : length of the product
    # i1,m1 : data for the first word
    # i2,m2 : data for the second word
    n = start
    while True:
        if filename is None:
            from sys import stdout
            f = stdout
        else:
            f = open(filename.format(dim,n), 'w')

        # TODO (see also in band relations): we can look only at the pair of
        # words
        f.write("# triangular relations for n = {}\n".format(n))
        for k in range(1,n):
            relations = []
            for i1,m1 in products_p(mats[0], mats[1], k-1, n-k, one_int):
                m1 = mats[0] * m1
                i1 = [0] + i1

                for i2,m2 in products_p(mats[0], mats[1], k-1, n-k, one_int):
                    m2 = mats[0] * m2
                    i2 = [0] + i2

                    if i1 == i2:
                        break

                    # here we first test equality between m1 and m2
                    # then we test the relations on all matrices in mats
                    # then we check formally using symbolic matrices
                    if m1 == m2 and \
                       is_relation(i1, i2, pairs) and \
                       prod(ab[x] for x in i1) == prod(ab[x] for x in i2):
                           relations.append(pretty_relation_string(i1,i2,'ab'))

#                for i2,m2 in products(mats[0], mats[1], k, n-k-1, one_int):
#                    m2 = mats[1] * m2
#                    i2 = [1] + i2
#
#                    # here we first test equality between m1 and m2
#                    # then we test the relations on all matrices in mats
#                    # then we check formally using symbolic matrices
#                    if m1 == m2 and \
#                       is_relation(i1, i2, pairs) and \
#                       prod(ab[x] for x in i1) == prod(ab[x] for x in i2):
#                           relations.append(pretty_relation_string(i1,i2))

            if relations:
                f.write("#  relations ({},{})\n".format(k,n-k))
                for r in relations:
                    f.write(r)
                    f.write('\n')
                f.flush()

        if filename is not None:
            f.close()
        n += 1

def relations_band_vc(dim, start=1, num_mat=500, filename=None):
    # the symbolic max-plus matrices
    ab = symbolic_max_plus_matrices_tri(dim, 2)
    one = symbolic_max_plus_identity(dim, ab[0].num_variables())

    # the integer max-plus matrices
    mats = [random_integer_max_plus_matrix_tri(dim, -50*dim*dim, 50*dim*dim) for _ in range(num_mat)]
    pairs = [(m1,m2) for m1 in mats for m2 in mats if m1 is not m2]
    one_int = integer_max_plus_matrix_identity(dim)

    # n     : length of the product
    # i1,m1 : data for the first word
    # i2,m2 : data for the second word
    n = start
    while True:
        if filename is None:
            from sys import stdout
            f = stdout
        else:
            f = open(filename.format(dim,n), 'w')

        f.write("# Band^cv_{} relations of length n = {}\n".format(dim,n))

        for n0 in range(1,n-1):
            eliminated_from_int = 0
            eliminated_from_symb = 0
            nb_identities = 0

            relations = []
            n1 = n - n0
            for i1,m1 in products_p(pairs[0][0], pairs[0][1], n0, n1, one_int):
                p,s = constraints_from_subwords(i1, dim)
                if p == -1:
                    continue

                mp = prod(pairs[0][k] for k in i1[:p])
                ms = prod(pairs[1][k] for k in i1[s:])

                nn1 = sum(i1[p:s])
                nn0 = s-p-nn1
                for ii2,mm2 in products_p(pairs[0][0], pairs[0][1], nn0, nn1, one_int):
                    i2 = i1[:p] + ii2 + i1[s:]
                    assert len(i1) == len(i2) and sum(i1) == sum(i2)
                    if i1 == i2:
                        break

                    if not is_relation(tuple(i1), tuple(i2), pairs):
                        eliminated_from_int += 1
                    elif prod(ab[x] for x in i1) != prod(ab[x] for x in i2):
                        eliminated_from_symb += 1
                    else:
                        nb_identities += 1
                        relations.append(pretty_relation_string(i1,i2,'xy'))

            print "({},{})".format(n0,n1)
            print "  int elimination  : {}".format(eliminated_from_int)
            print "  symb elimination : {}".format(eliminated_from_symb)
            print "  identities       : {}".format(nb_identities)

            if relations:
                f.write("#  relations ({},{})\n".format(n0,n1))
                for r in relations:
                    f.write(r)
                    f.write('\n')
                    f.flush()

        if filename is not None:
            f.close()
        n += 1

def conj_min_relations_band_vc(dim, start=1, num_mat=500, filename=None):
    r"""
    Look at relations of the form

    M * x*y * M == M * y*x * M

    for matrices in B^vc where M contains the same number of x and y.
    """
    # the symbolic max-plus matrices
    ab = symbolic_max_plus_matrices_band(dim, 2, diag='v', surdiag='c')
    one = symbolic_max_plus_identity(dim, ab[0].num_variables())
    xy = ab[0]*ab[1]
    yx = ab[1]*ab[0]

    # the integer max-plus matrices
    pairs = [random_integer_max_plus_matrices_band(dim, -50*dim*dim, 50*dim*dim,
        ord('v'), ord('c')) for _ in range(num_mat)]
    one_int = integer_max_plus_matrix_identity(dim)

    # n  : length of m
    # i  : word of m
    # i1 : word of mxym
    # i2 : word of myxm
    n = 1
    while True:
        if filename is None:
            from sys import stdout
            f = stdout
        else:
            f = open(filename.format(dim,n), 'w')

        header = "# Band^vc_{} relations with n = {}".format(dim,n)
        f.write(header)
        f.write("\n")
        if filename:
            print header

        eliminated_from_int = 0
        eliminated_from_symb = 0
        nb_identities = 0
        for i in product_p(n-1, n):
            i = (0,) + tuple(i)
            i1 = i + (0,1) + i
            i2 = i + (1,0) + i

            if not is_relation(i1, i2, pairs, upper=True):
                eliminated_from_int += 1
                continue

            print "  potential relation:", pretty_relation_string(i1,i2,'xy')
            sys.stdout.flush()

            t0 = time.clock()
            m = prod(ab[x] for x in i)
            t1 = time.clock()
            print "    m    computed in {} seconds".format(t1-t0)
            print "    nb pts in convex hull: {}".format([len(m[i,j]) for i in range(dim) for j in range(i,dim)])
            t0 = time.clock()
            mxym = m * xy * m
            t1 = time.clock()
            print "    mxym computed in {} seconds".format(t1-t0)
            print "    nb pts in convex hull: {}".format([len(mxym[i,j]) for i in range(dim) for j in range(i,dim)])
            t0 = time.clock()
            myxm = m * yx * m
            t1 = time.clock()
            print "    myxm computed in {} seconds".format(t1-t0)
            print "    nb pts in convex hull: {}".format([len(myxm[i,j]) for i in range(dim) for j in range(i,dim)])
            if mxym != myxm:
                print "    ... not a relation\n"
                sys.stdout.flush()
                eliminated_from_symb += 1
                continue

            print "    NEW RELATION\n"
            nb_identities += 1
            f.write(pretty_relation_string(i1,i2,'xy'))
            f.write('\n')
            f.flush()
            print

        print "  int elimination  : {}".format(eliminated_from_int)
        print "  symb elimination : {}".format(eliminated_from_symb)
        print "  identities       : {}".format(nb_identities)
        print

        if filename is not None:
            f.close()
        if nb_identities:
            return
        n += 1
        if n == 9:
            break

# TODO
#  relations in B^sv?
#  d-relations for B^cv are {(u,v): Subwords_{d-1}(u) = Subwords_{d-1}(v)}

def filter_relations(r, mats, ab, one):
    ans = []
    for r1,r2 in r:
        if is_relation(r1, r2, mats) and \
            prod(ab[x] for x in r1) == prod(ab[x] for x in r2):
            ans.append(pretty_relation_string(r1,r2))
    return ans

def brute_force_fibo(dim):
    dim = int(dim)
    w = list(words.FibonacciWord([0,1])[:1000])
    n = 1

    mats = random_integer_max_plus_matrices_tri_sim_diag(dim, -2**20, 2**20)

    while True:
        u = prod(mats[w[i]] for i in range(n))
        v = prod(mats[w[i]] for i in range(n-1,-1,-1))

        if v*mats[0]*mats[1]*u != v*mats[1]*mats[0]*u:
            n += 1
            print n
        else:
            mats = random_integer_max_plus_matrices_tri_sim_diag(dim, -2**20, 2**20)



