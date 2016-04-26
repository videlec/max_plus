r"""
Identities in the monoid of max-plus matrices

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

For diagonal matrices, the polyhedra of a product of length l is always
contained in the hyperplane x_0 + ... + x_{n-1} = l.

EXAMPLES::

    sage: from max_plus import *

    sage: x,y=symbolic_max_plus_matrices(2,2,ch='ppl')
    sage: x.convex_hull_engine()
    'ppl'
    sage: x,y=symbolic_max_plus_matrices(2,2,ch='PALP')
    sage: x.convex_hull_engine()
    'PALP'
    sage: x,y=symbolic_max_plus_matrices(2,2,ch='cdd')
    sage: x.convex_hull_engine()
    'cdd'

For 2x2 matrices, it seems that the Newton polytopes of the entries always
belong to a same given subspace of rank 5 (where we have 8 variables)::

    sage: x,y = symbolic_max_plus_matrices(2,2)
    sage: print x
    [ x0  x1 ]
    [ x2  x3 ]
    sage: print y
    [ x4  x5 ]
    [ x6  x7 ]
    sage: z = x*x*x*y*x*x*y*y
    sage: V = z.get_vector_span(0,0)
    sage: print V
    Free module of degree 8 and rank 5 over Integer Ring
    Echelon basis matrix:
    [ 1  0  0 -1  0  0  0  0]
    [ 0  1  0 -1  0  0  1 -1]
    [ 0  0  1 -1  0  0 -1  1]
    [ 0  0  0  0  1  0  0 -1]
    [ 0  0  0  0  0  1  1 -2]
    sage: V == z.get_vector_span(0,1) == z.get_vector_span(1,0) == z.get_vector_span(1,1)
    True
    sage: z1 = x*y*y*y*x*x
    sage: z2 = x*x*x*x
    sage: z3 = y*y*y*y
    sage: z4 = x*y
    sage: all(z.get_vector_span(i,j).is_submodule(V) for i in (0,1) for j in (0,1) for z in (z1,z2,z3))
    True

.. TODO::

    Consider the matrices as matrices over flat polynomials with countably many
    variables. In particular, increase ``num_vars`` in multiplication if needed.

    variables: x0, x1, x2, ...
    -> needs convention for zero and -infinity

    Rectangular matrices

    Polyhedron.intersection_assign(Polyhedron)
    Polyhedron.poly_hull_assign(Polyhedron) -> compute convex hull of union
    Polyhedron.poly_difference_assign(Polyhedron)

    PPL methods:
    
    - `poly_hull_assign` method for unions
    - `intersection_assign` method for intersection
    - `relation_with` for containment
"""
from __future__ import print_function, division, absolute_import

import itertools
import multiprocessing as mp
from itertools import product
from datetime import datetime
import os
import sys
from time import time

from .sage_import import *
from .convex_hull import get_convex_hull_engine

##########################################
# Main function helper to build matrices #
##########################################

def symbolic_max_plus_identity(d, nvar, ch=None):
    r"""
    Return the ``d x d`` identity matrices on ``nvar`` variables.

    EXAMPLES::

        sage: from max_plus import *
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

def symbolic_max_plus_matrices(d, n, ch=None, typ='sym'):
    r"""
    Return ``n`` independent symbolic matrices in dimension ``d``.

    INPUT:

    - ``d`` -- the dimension

    - ``n`` -- the number of matrices

    - ``ch`` -- the convex hull engine

    - ``typ`` -- either ``full`` (store all coeffs), ``sym`` storing only two
      coeffs or ``quick`` for matrices with fast multiplication but slow
      comparison.

    TESTS::

        sage: from max_plus import *

    Test a relation for 2x2 matrices::

        sage: A,B = symbolic_max_plus_matrices(2,2)
        sage: U = A*B*B*B*B*A*A
        sage: V = A*A*B*B*B*B*A
        sage: L = U*A*A*B*B*V
        sage: R = U*B*B*A*A*V
        sage: L == R
        True

    Check the validity of the above computation with dense matrices::

        sage: C,D = symbolic_max_plus_matrices(2, 2, typ='full')
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

    Compatibility between the types ``'full'`` and ``'sym'``::

        sage: a1,b1,c1 = symbolic_max_plus_matrices(3, 3, typ='full')
        sage: a2,b2,c2 = symbolic_max_plus_matrices(3, 3, typ='sym')
        sage: assert a1 == a2 and b1 == b2 and c1 == c2
        sage: assert a1*b1 == a2*b2 == a2*b1 == a1*b2
        sage: assert a1*b1*c1*b1*a1 == a2*b2*c2*b2*a2
    """
    d = int(d)
    n = int(n)
    if d <= 0:
        raise ValueError("d (= {}) must be postive".format(d))

    nvar = n * d * d

    V = FreeModule(ZZ, nvar)
    B = ((b,) for b in V.basis())

    matrices = []

    if d == 1:
        typ = 'full'

    if typ == 'sym' or typ == 'quick':
        z = [0]*nvar
        for i in range(n):
            z[i*d*d] = 1
            diag = (V(z),)
            z[i*d*d] = 0

            z[i*d*d+1] = 1
            nondiag = (V(z),)
            z[i*d*d+1] = 0

            if typ == 'sym':
                matrices.append(SymbolicSymmetricMaxPlusMatrix(d, n, diag, nondiag, ch))
            else:
                matrices.append(QuickSymbolicSymmetricMaxPlusMatrix(d, n, diag, nondiag, ch))
    elif typ == 'full':
        for i in range(n):
            mat = []
            for j in range(d):
                mat.append([next(B) for k in range(d)])
            matrices.append(SymbolicMaxPlusMatrix(d, nvar, mat, ch))
    else:
        raise ValueError

    return matrices

def symbolic_max_plus_matrices_band(d, n,
        diag='v', surdiag='v', ch=None, typ='sym'):
    r"""
    INPUT:

    - ``d`` -- dimension

    - ``n`` -- number of matrices

    - ``diag`` -- either ``'z'`` (for zero), ``'c'`` (for constant), ``'s'``
      (for same, i.e. equal in each matrix) or ``'v'`` (i.e. independent in each
      matrices).

    - ``surdiag`` -- one of ``'z'``, ``'c'``, ``'s'``, ``'v'``

    - ``ch`` -- convex hull engine to use

    - ``typ`` -- either ``'sym'`` (default) or ``'full'``. If ``'sym'``  uses
      condensed matrix form that lower the number of convex hull computations
      (`d` instead of `d^2`)

    EXAMPLES::

        sage: from max_plus import *

    For 'zv' (or equivalently 'cv') the relations are the pairs `(u,v)` so that
    `Subwords_{d-1}(u) = Subwords_{d-1}(v)`::

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

    "Sparse" versus "dense" representations::

        sage: a1,b1 = symbolic_max_plus_matrices_band(4,2,'s','v',typ='full')
        sage: a2,b2 = symbolic_max_plus_matrices_band(4,2,'s','v',typ='sym')
        sage: a1 == a2 and b1 == b2
        True
        sage: a1*a1 == a2*a2 and a1*b1 == a2*b2
        True
        sage: b1*a1*b1*a1*a1 == b2*a2*b2*a2*a2
        True

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

        sage: for a in 'zcsv':
        ....:     for b in 'zcsv':
        ....:         x,y = symbolic_max_plus_matrices_band(3,2,a,b,typ='full')
        ....:         entries = set().union(x.list(), y.list())
        ....:         nvars = len(entries) - (2 if a == 'z' or b == 'z' else 1)
        ....:         assert x.num_vars() == y.num_vars() == nvars

    Tests of "sparse" versus "dense"::

        sage: for diag in 'zsv':
        ....:     for surdiag in 'zsv':
        ....:         a1,b1 = symbolic_max_plus_matrices_band(4,2,diag,surdiag,typ='full')
        ....:         a2,b2 = symbolic_max_plus_matrices_band(4,2,diag,surdiag,typ='sym')
        ....:         assert a1 == a2 and b1 == b2 and a1*b1*a1 == a2*b2*a2
    """
    assert isinstance(diag,str) and len(diag) == 1 and diag in 'csvz'
    assert isinstance(surdiag,str) and len(surdiag) == 1 and surdiag in 'csvz'

    if typ == 'sym' and (diag=='c' or surdiag=='c'):
        raise ValueError("typ='sym' option not compatible with diag='c' or surdiag='c'")

    d = int(d)
    n = int(n)
    assert d > 0 and n > 0

    nvar = 0
    if diag == 'c':
        nvar += 1
    elif diag == 's':
        nvar += d
    elif diag == 'v':
        nvar += n*d
    if surdiag == 'c':
        nvar += 1
    elif surdiag == 's':
        nvar += d-1
    elif surdiag == 'v':
        nvar += n*(d-1)

    V = FreeModule(ZZ, nvar)
    B = iter(V.basis())
    e = ()
    zero = (V.zero(),)
    mat_init = [[e]*d for _ in range(d)]

    if diag == 'z':
        for k in range(d):
            mat_init[k][k] = zero
    elif diag == 'c':
        f = (next(B),)
        for k in range(d):
            mat_init[k][k] = f
    elif diag == 's':
        for k in range(d):
            mat_init[k][k] = (next(B),)

    if surdiag == 'z':
        for k in range(d-1):
            mat_init[k][k+1] = zero
    elif surdiag == 'c':
        f = (next(B),)
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

        if typ == 'sym':
            matrices.append(SymbolicSymmetricMaxPlusUpperMatrix(d, nvar, tuple(mat[0]), ch))
        elif typ == 'full':
            matrices.append(SymbolicMaxPlusMatrix(d, nvar, tuple(mat), ch))
        else:
            raise ValueError("typ must be 'sym' or 'full'")

    return matrices

def symbolic_max_plus_matrices_upper(d, n,
        diag='v', surdiag='v', ch=None, typ='sym'):
    r"""
    EXAMPLES::

        sage: from max_plus import *

        sage: x,y = symbolic_max_plus_matrices_upper(3,2)
        sage: x
        A 3x3 symbolic max plus matrix on 12 variables
        sage: print x
        [  x0   x3  x5 ]
        [ -oo   x1  x4 ]
        [ -oo  -oo  x2 ]
        sage: print y
        [  x6   x9  x11 ]
        [ -oo   x7  x10 ]
        [ -oo  -oo   x8 ]

        sage: x,y = symbolic_max_plus_matrices_upper(3,2,'z','v')
        sage: print x
        [   0   x0  x2 ]
        [ -oo    0  x1 ]
        [ -oo  -oo   0 ]
        sage: print y
        [   0   x3  x5 ]
        [ -oo    0  x4 ]
        [ -oo  -oo   0 ]

        sage: x,y = symbolic_max_plus_matrices_upper(3,2,'v','s')
        sage: print x
        [  x3   x0  x2 ]
        [ -oo   x4  x1 ]
        [ -oo  -oo  x5 ]
        sage: print y
        [  x6   x0  x2 ]
        [ -oo   x7  x1 ]
        [ -oo  -oo  x8 ]

        sage: for diag in 'zsv':
        ....:     for surdiag in 'zsv':
        ....:         a1,b1 = symbolic_max_plus_matrices_upper(4,2,diag,surdiag,typ='full')
        ....:         a2,b2 = symbolic_max_plus_matrices_upper(4,2,diag,surdiag,typ='sym')
        ....:         assert a1 == a2 and b1 == b2 and a1*b1*a1 == a2*b2*a2

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
    assert isinstance(diag,str) and len(diag) == 1 and diag in 'csvz'
    assert isinstance(surdiag,str) and len(surdiag) == 1 and surdiag in 'csvz'

    d = int(d)
    n = int(n)
    assert d > 0 and n > 0

    if typ == 'sym' and (diag == 'c' or surdiag == 'c'):
        raise ValueError("typ='sym' is not compatible with diag=c or surdiag=c")

    nvar = 0
    if diag == 'c':
        nvar += 1
    elif diag == 's':
        nvar += d
    elif diag == 'v':
        nvar += n*d

    if surdiag == 'c':
        nvar += 1
    elif surdiag == 's':
        nvar += d*(d-1)//2
    elif surdiag == 'v':
        nvar += n*d*(d-1)//2

    V = FreeModule(ZZ, nvar)
    B = ((b,) for b in V.basis())

    zero = (V.zero(),)
    e = ()
    mat_init = [[e]*d for _ in range(d)]

    if diag in 'zcs':
        if diag == 'z':
            f = itertools.repeat(zero)
        elif diag == 'c':
            f = itertools.repeat(next(B))
        elif diag == 's':
            f = B
        for k in range(d):
            mat_init[k][k] = next(f)

    if surdiag in 'zcs':
        if surdiag == 'z':
            f = itertools.repeat(zero)
        elif surdiag == 'c':
            f = itertools.repeat(next(B))
        else:
            f = B
        for h in range(1,d):
            for k1 in range(d-h):
                mat_init[k1][k1+h] = next(f)

    matrices = []
    for i in range(n):
        mat = [row[:] for row in mat_init]

        if diag == 'v':
            for k in range(d):
                mat[k][k] = next(B)
        if surdiag == 'v':
            for h in range(1,d):
                for k in range(d-h):
                    mat[k][k+h] = next(B)

        if typ == 'sym':
            matrices.append(SymbolicSymmetricMaxPlusUpperMatrix(d, nvar, mat[0], ch))
        elif typ == 'full':
            matrices.append(SymbolicMaxPlusMatrix(d, nvar, mat, ch))
        else:
            raise ValueError("typ must either be 'sym' or 'full'")

    return matrices

###########################
# Helper for pretty print #
###########################

def str_linear_form(v):
    if not v:
        return '0'
    else:
        return '+'.join('{}x{}'.format('' if j == 1 else j,i) for i,j in enumerate(v) if j)

def pretty_print_poly(p):
    if not p:
        return '-oo'
    elif len(p) == 1:
        if p[0].is_zero():
            return '0'
        else:
            return str_linear_form(p[0])
    else:
        return 'max(' + ', '.join(str_linear_form(v) for v in p) + ')'

##################
# Matrix classes #
##################

class SymbolicMaxPlusMatrix(SageObject):
    r"""
    A symbolic max plus matrix.

    EXAMPLES::

        sage: from max_plus import *

        sage: M1 = symbolic_max_plus_identity(2,2)
        sage: M2 = M1 * M1
        sage: M2 == M1
        True

        sage: x1,y1,z1 = symbolic_max_plus_matrices_band(3,3,'s','v', ch='ppl')
        sage: x2,y2,z2 = symbolic_max_plus_matrices_band(3,3,'s','v', ch='cdd')
        sage: x3,y3,z3 = symbolic_max_plus_matrices_band(3,3,'s','v', ch='PALP')
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
    def __init__(self, d, nvars, data, ch=None, check=True):
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
        self._data = data

        self.convex_hull = get_convex_hull_engine(self._nvars, ch)

        if check:
            self._data = tuple(tuple(tuple(map(self.convex_hull.V,data[i][j])) \
                      for j in range(d)) \
                      for i in range(d))
            for i in range(self._d):
                for j in range(self._d):
                    row = self._data[i][j]
                    for v in row: v.set_immutable()

    def __hash__(self):
        r"""
        TESTS::

            sage: from max_plus import *
            sage: a,b = symbolic_max_plus_matrices(2,2,typ='full')
            sage: hash(a)
            Traceback (most recent call last):
            ...
            TypeError
            sage: a,b = symbolic_max_plus_matrices(1,2,typ='full')
            sage: hash(a) == hash(b)
            False
            sage: hash(a*b*a) == hash(b*a*a) == hash(a*a*b)
            True
        """
        if self._d == 1:
            return hash(self._data)
        raise TypeError

    def convex_hull_engine(self):
        r"""
        Return a string that describes the convex hull engine.

        EXAMPLES::

            sage: from max_plus import *

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
        if type(self) != type(other) and \
            not isinstance(self, SymbolicMaxPlusMatrix) and \
            not isinstance(other, SymbolicMaxPlusMatrix):
                raise TypeError("can not multiply {} with {}".format(type(self),type(other)))
        if self._d != other._d:
            raise TypeError("dimension or number of variable mismatch")
        d = self._d
        new_data = [[None]*d for _ in range(d)]
        for i in range(d):
            for j in range(d):
                vertices = set()
                for k in range(d):
                    l = [x+y for x in self[i,k] for y in other[k,j]]
                    for v in l: v.set_immutable()
                    vertices.update(l)
                new_data[i][j] = self.convex_hull(vertices)

        return SymbolicMaxPlusMatrix(d, self._nvars, tuple(new_data), self.convex_hull)

    def list(self):
        r"""
        Return the list of entries of this matrix.

        EXAMPLES::

            sage: from max_plus import *

            sage: a, = symbolic_max_plus_matrices(3,1)
            sage: a.list()
            [((1, 0, 0, 0, 0, 0, 0, 0, 0),),
             ((0, 1, 0, 0, 0, 0, 0, 0, 0),),
             ((0, 0, 1, 0, 0, 0, 0, 0, 0),),
             ((0, 0, 0, 1, 0, 0, 0, 0, 0),),
             ((0, 0, 0, 0, 1, 0, 0, 0, 0),),
             ((0, 0, 0, 0, 0, 1, 0, 0, 0),),
             ((0, 0, 0, 0, 0, 0, 1, 0, 0),),
             ((0, 0, 0, 0, 0, 0, 0, 1, 0),),
             ((0, 0, 0, 0, 0, 0, 0, 0, 1),)]
        """
        return [self[i,j] for i in range(self._d) for j in range(self._d)]

    def __pow__(self, n):
        r"""
        TESTS::

            sage: from max_plus import *

            sage: a, = symbolic_max_plus_matrices(3,1)
            sage: a^0 == symbolic_max_plus_identity(3,9)
            True
            sage: a^1 == a
            True
            sage: a^2 == a*a
            True
            sage: a^3 == a*a*a
            True
        """
        n = int(n)
        one = symbolic_max_plus_identity(self._d, self._nvars, self.convex_hull)
        return generic_power(self, n, one)

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

            sage: from max_plus import *

            sage: A,B = symbolic_max_plus_matrices(2, 2, typ='sym')
            sage: C,D = symbolic_max_plus_matrices(2, 2, typ='full')
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

            sage: from max_plus import *

            sage: A,B = symbolic_max_plus_matrices(2, 2, typ='sym')
            sage: C,D = symbolic_max_plus_matrices(2, 2, typ='full')
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
        return not self == other

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

            sage: from max_plus import *

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
        from max_plus_int import minus_infinity, IntegerMaxPlusMatrix
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

    def get_vector_span(self, i, j):
        r"""
        Return the dimension of the affine space spanned generated by each
        Newton polytopes.

        for triangular matrices, seems to stay 0 on the diagonal.
        """
        from sage.rings.infinity import Infinity
        from sage.matrix.constructor import matrix
        data = self[i,j]
        if not data:
            return None
        elif len(data) == 1:
            return FreeModule(ZZ, self._nvars).submodule([])
        else:
            return matrix([x-data[0] for x in data]).row_space()

def vertex_swap(d, n, l, i1, i2, j1, j2):
    r"""
    Permutations (i1,j1) -> (i2,j2)

    Make an exchange of rows/columns in matrix data. This is used in
    multiplication of full symbolic matrices.

    INPUT:

    - ``d`` - dimension

    - ``n`` - number of matrices

    - ``l`` - data (list of integer vectors of size `n d^2`)

    - ``i1``, ``i2``, ``j1``, ``j2`` -- matrix indices in `\{0, 1, ..., d-1\}`
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
    Transposition (i,j) in the vector v
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
    permutation (i,j,k) in the vector v
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

def vertex_cyclic_swap(nvars, l, i):
    r"""
    Perform a cyclic swap on the vertices.

    This is used in multiplication of symbolic upper matrices.

    Currently it is suboptimal but on the other hand, this cost much less than
    whatever convex hull computation.
    """
    if i == 0 or not l:
        return l
    ll = []
    F = l[0].parent()
    for v in l:
        assert not v[-i:]
        ll.append(F(tuple(v[-i:]) + tuple(v[:-i])))
    for v in ll: v.set_immutable()
    return tuple(ll)

class SymbolicSymmetricMaxPlusUpperMatrix(SymbolicMaxPlusMatrix):
    r"""
    By convention, if there is a variable x_k in position (i,j) the variable
    x_{k+1} must be in position (i+1,j+1).
    """
    def __init__(self, d, nvars, row, ch=None):
        self._d = d              # size
        self._nvars = nvars      # number of variables
        self._row = row          # the terms (0,i)

        for data in self._row:
            for v in data:
                v.set_immutable()

        self.convex_hull = get_convex_hull_engine(self._nvars, ch)

    def __hash__(self):
        r"""
        TESTS::

            sage: from max_plus import *
            sage: a,b = symbolic_max_plus_matrices_band(2,2,'s','v')
            sage: hash(a*b*a*a*b) == hash(a*b*b*a*b) 
            True
        """
        return hash(self._row)

    def __eq__(self, other):
        if type(self) is type(other):
            return self._d == other._d and \
                   self._nvars == other._nvars and \
                   self._row == other._row
        else:
            return SymbolicMaxPlusMatrix.__eq__(self, other)

    def __getitem__(self, data):
        r"""
        TESTS::

            sage: from max_plus import *

            sage: a,b = symbolic_max_plus_matrices_band(4,2,'s','v')
            sage: type(a)
            <class 'max_plus.max_plus_symbolic.SymbolicSymmetricMaxPlusUpperMatrix'>
            sage: a[0,1]
            ((0, 0, 0, 0, 1, 0, 0, 0, 0, 0),)
            sage: a[1,2]
            ((0, 0, 0, 0, 0, 1, 0, 0, 0, 0),)
            sage: (a*b*a)[1,3]
            ((0, 0, 0, 1, 0, 1, 0, 0, 0, 1),
             (0, 0, 1, 0, 0, 1, 1, 0, 0, 0),
             (0, 1, 0, 0, 0, 0, 1, 0, 1, 0))
        """
        i,j = data
        i = int(i); j = int(j)
        if i < 0 or i >= self._d or j < 0 or j >= self._d:
            raise IndexError("matrix index out of range")
        if i > j:
            return ()
        else:
            return vertex_cyclic_swap(self._nvars, self._row[j-i], i)

    def __mul__(self, other):
        if type(self) is not type(other):
            return SymbolicMaxPlusMatrix.__mul__(self, other)

        assert self._nvars == other._nvars and self._d == other._d

        # NOTE: there are some cases where the minkowski sum can be much faster
        # to compute: when the points are in disjoint subspaces as for example:
        #  C1 = [(1,1,0,0), (-1,3,0,0)]
        #  C2 = [(0,0,5,-2), (0,0,1,3), (0,0,-2,2)]
        # it is not clear whether we will gain in the most expensive computation
        # of the (0,d-1) coefficient.

        row = []
        for j in range(self._d):
            vertices = set()
            for k in range(0,j+1):
                l = [x+y for x in self[0,k] for y in other[k,j]]
                for v in l: v.set_immutable()
                vertices.update(l)
            row.append(self.convex_hull(vertices))

        return SymbolicSymmetricMaxPlusUpperMatrix(self._d, self._nvars,
                tuple(row), self.convex_hull)

class SymbolicSymmetricMaxPlusMatrix(SymbolicMaxPlusMatrix):
    r"""
    Matrices that are symmetric under permutations.

    This class only concerns relations on the full matrix monoid. In practice,
    we should win a factor d^2/2 where d is the dimension.

    Such matrices are completely determined by the entries (0,0) and (0,1) (that
    are stored in the attributes ``self._diag`` and ``self._nondiag``).

    EXAMPLES::

        sage: from max_plus import *

        sage: a,b = symbolic_max_plus_matrices(2,2)
        sage: print a
        [ x0  x1 ]
        [ x2  x3 ]
        sage: print b
        [ x4  x5 ]
        [ x6  x7 ]
        sage: print a*b
        [ max(x1+x6, x0+x4)  max(x1+x7, x0+x5) ]
        [ max(x3+x6, x2+x4)  max(x3+x7, x2+x5) ]

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

        for v in self._diag:
            v.set_immutable()
        for v in self._nondiag:
            v.set_immutable()

    def __hash__(self):
        r"""
        TESTS::

            sage: from max_plus import *
            sage: a,b = symbolic_max_plus_matrices(2,2)
            sage: hash(a) == hash(b)
            False
            sage: hash(a*b) == hash(b*a)
            False
        """
        return hash(self._diag) ^ hash(self._nondiag)

    def __eq__(self, other):
        r"""
        TESTS::

            sage: from max_plus import *

            sage: x1,y1 = symbolic_max_plus_matrices(5, 2, typ='sym')
            sage: x2,y2 = symbolic_max_plus_matrices(5, 2, typ='full')
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

    def __getitem__(self, data):
        r"""
        TESTS::

            sage: from max_plus import *

            sage: x,y = symbolic_max_plus_matrices(3, 2, typ='sym')
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

            sage: from max_plus import *

            sage: x1,y1 = symbolic_max_plus_matrices(5, 2, typ='sym')
            sage: x2,y2 = symbolic_max_plus_matrices(5, 2, typ='full')
            sage: assert x1 == x2 and y1 == y2
            sage: assert x1*y1 == x2*y2 and y1*x1 == y2*x2
            sage: assert x1*x1 == x2*x2
            sage: assert x1*y1*x1 == x2*y2*x2
            sage: assert x1*x1*x1 == x2*x2*x2

            sage: x1,y1 = symbolic_max_plus_matrices(4, 2, typ='sym')
            sage: x2,y2 = symbolic_max_plus_matrices(4, 2, typ='full')
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


class QuickSymbolicSymmetricMaxPlusMatrix(SymbolicSymmetricMaxPlusMatrix):
    def __eq__(self, other):
        if type(self) != type(other):
            raise TypeError
        for u,v in [(self._diag, other._diag),
                    (self._nondiag, other._nondiag)]:
            u = set(u)
            uv = u.intersection(v)
            out = u.symmetric_difference(v)
            for q in out:
                if not in_convex_hull_LP(uv, q):
                    return False
        return True

    def __mul__(self, other):
        if type(self) != type(other):
            raise TypeError

        assert self._n == other._n and self._d == other._d

        data = []
        for j in range(2):
            vertices = set()
            for k in range(self._d):
                l = [x+y for x in self[0,k] for y in other[k,j]]
                for v in l: v.set_immutable()
                vertices.update(l)
            data.append(self.clean(vertices))

        return QuickSymbolicSymmetricMaxPlusMatrix(self._d, self._n, data[0],
                data[1])

    def clean(self, pts):
        return clean(list(pts), self._nvars - self._d + 1)


