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

"""

import itertools
import multiprocessing as mp
from itertools import product
from datetime import datetime
import os
import sys
from subprocess import Popen, PIPE
from time import time

from sage.misc.misc import SAGE_TMP
from sage.misc.temporary_file import tmp_filename
from sage.misc.misc import SAGE_TMP

from sage.structure.sage_object import SageObject
from sage.structure.element import generic_power

try:
    from sage.combinat.words.words import FiniteWords
except ImportError:
    from sage.combinat.words.words import Words as FiniteWords
from sage.combinat.words.suffix_trees import SuffixTrie

try:
    from ppl import (Variable, C_Polyhedron, point, Generator_System,
         Linear_Expression, Constraint_System, MIP_Problem)
except ImportError:
    from sage.libs.ppl import (Variable, C_Polyhedron, point, Generator_System,
         Linear_Expression, Constraint_System, MIP_Problem)

from sage.numerical.mip import MIPSolverException
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

    Compatibility of output with ``sym=False`` and ``sym=True``::

        sage: a1,b1,c1 = symbolic_max_plus_matrices(3,3,sym=False)
        sage: a2,b2,c2 = symbolic_max_plus_matrices(3,3,sym=True)
        sage: assert a1 == a2 and b1 == b2 and c1 == c2
        sage: assert a1*b1 == a2*b2 == a2*b1 == a1*b2
        sage: assert a1*b1*c1*b1*a1 == a2*b2*c2*b2*a2
    """
    d = int(d)
    n = int(n)
    nvar = n * d * d

    V = FreeModule(ZZ, nvar)
    B = ((b,) for b in V.basis())

    matrices = []

    if sym:
        z = [0]*nvar
        for i in range(n):
            z[i*d*d] = 1
            diag = (V(z),)
            z[i*d*d] = 0

            z[i*d*d+1] = 1
            nondiag = (V(z),)
            z[i*d*d+1] = 0

            matrices.append(SymbolicSymmetricMaxPlusMatrix(d, n, diag, nondiag, ch))
    else:
        for i in range(n):
            mat = []
            for j in range(d):
                mat.append([next(B) for k in range(d)])
            matrices.append(SymbolicMaxPlusMatrix(d, nvar, mat, ch))

    return matrices

def symbolic_max_plus_matrices_band(d, n,
        diag='v', surdiag='v', ch=None, sym=True):
    r"""
    INPUT:

    - ``d`` -- dimension

    - ``n`` -- number of matrices

    - ``diag`` -- either ``'z'`` (for zero), ``'c'`` (for constant), ``'s'``
      (for same, i.e. equal in each matrix) or ``'v'`` (i.e. independent in each
      matrices).

    - ``surdiag`` -- one of ``'z'``, ``'c'``, ``'s'``, ``'v'``

    - ``ch`` -- convex hull engine to use

    - ``sym`` -- if ``True`` uses condensed matrix form that lower the number of
      convex hull computations (`d` instead of `d^2`)

    EXAMPLES:

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

        sage: a1,b1 = symbolic_max_plus_matrices_band(4,2,'s','v',sym=False)
        sage: a2,b2 = symbolic_max_plus_matrices_band(4,2,'s','v',sym=True)
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
        ....:         x,y = symbolic_max_plus_matrices_band(3,2,a,b,sym=False)
        ....:         entries = set().union(x.list(), y.list())
        ....:         nvars = len(entries) - (2 if a == 'z' or b == 'z' else 1)
        ....:         assert x.num_vars() == y.num_vars() == nvars

    Tests of "sparse" versus "dense"::

        sage: for diag in 'zsv':
        ....:     for surdiag in 'zsv':
        ....:         a1,b1 = symbolic_max_plus_matrices_band(4,2,diag,surdiag,sym=False)
        ....:         a2,b2 = symbolic_max_plus_matrices_band(4,2,diag,surdiag,sym=True)
        ....:         assert a1 == a2 and b1 == b2 and a1*b1*a1 == a2*b2*a2
    """
    assert isinstance(diag,str) and len(diag) == 1 and diag in 'csvz'
    assert isinstance(surdiag,str) and len(surdiag) == 1 and surdiag in 'csvz'

    if sym and (diag=='c' or surdiag=='c'):
        raise ValueError("sym=True option not compatible with diag='c' or surdiag='c'")

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

        if sym:
            matrices.append(SymbolicSymmetricMaxPlusUpperMatrix(d, nvar, mat[0], ch))
        else:
            matrices.append(SymbolicMaxPlusMatrix(d, nvar, mat, ch))

    return matrices

def symbolic_max_plus_matrices_upper(d, n,
        diag='v', surdiag='v', ch=None, sym=True):
    r"""
    EXAMPLES::

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
        ....:         a1,b1 = symbolic_max_plus_matrices_upper(4,2,diag,surdiag,sym=False)
        ....:         a2,b2 = symbolic_max_plus_matrices_upper(4,2,diag,surdiag,sym=True)
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

    if sym and (diag == 'c' or surdiag == 'c'):
        raise ValueError("sym=True is not compatible with diag=c or surdiag=c")

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

        if sym:
            matrices.append(SymbolicSymmetricMaxPlusUpperMatrix(d, nvar, mat[0], ch))
        else:
            matrices.append(SymbolicMaxPlusMatrix(d, nvar, mat, ch))

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

        if check:
            self._data = tuple((tuple(x) for x in data))
            for i,row in enumerate(self._data):
                if len(row) != d:
                    raise ValueError("wrong {}-th row size".format(i))
                for p in row:
                    assert isinstance(p, tuple), "got p={}".format(p)

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

        return SymbolicMaxPlusMatrix(d, self._nvars, new_data, self.convex_hull)

    def list(self):
        r"""
        Return the list of entries of this matrix.

        EXAMPLES::

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

        self.convex_hull = get_convex_hull_engine(self._nvars, ch)

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

            sage: a,b = symbolic_max_plus_matrices_band(4,2,'s','v')
            sage: type(a)
            <class '__main__.SymbolicSymmetricMaxPlusUpperMatrix'>
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

        return SymbolicSymmetricMaxPlusUpperMatrix(self._d, self._nvars, row, self.convex_hull)


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


#########################
# Combinatorial methods #
#########################
def occurrences(w, u):
    r"""
    Return the set of occurrences of ``u`` in ``w``.

    If the word ``w`` contains one or more jokers (i.e. the letter ``'*'``) then
    a pair of lists is returned. The first one is made of occurrences that does
    not pass through joker and the other one is the complement.

    EXAMPLES::

        sage: occurrences('abbabab', 'abb')
        [(0, 1, 2), (0, 1, 4), (0, 2, 4), (0, 1, 6), (0, 2, 6),
         (0, 4, 6), (3, 4, 6)]

        sage: s,t=occurrences('xyyxx*yyxxy', 'xyx')
        sage: s
        [(0, 1, 3), (0, 2, 3), (0, 1, 4), ..., (3, 7, 9), (4, 7, 9)]
        sage: t
        [(0, 1, 5), (0, 2, 5), (0, 5, 8), ..., (5, 6, 9), (5, 7, 9)]
    """
    # NOTE: this method is efficient because it parses only once the word w. On
    # the other hand there is a clear waste of memory if we just want to go
    # through the occurrences and not make the list of them (which is *not* our
    # case). Though we could eliminate some of the non-extremal occurrences but
    # that would be some work to implement (?)
    pos = [[] for _ in range(len(u))]
    pos2 = [[] for _ in range(len(u))]

    for i,letter in enumerate(w):
        for j in range(len(u)-1,0,-1):
            if letter == '*':
                pos2[j].extend(x + (i,) for x in pos[j-1])
            elif letter == u[j]:
                pos[j].extend(x + (i,) for x in pos[j-1])
                pos2[j].extend(x + (i,) for x in pos2[j-1])
        if letter == '*':
            pos2[0].append((i,))
        elif letter == u[0]:
            pos[0].append((i,))

    return pos[-1] if not pos2[0] else (pos[-1], pos2[-1])

def extremal_occurrences(w, u, verbose=False):
    r"""
    Return the set of extremal occurrences of ``u`` in ``w``.

    An occurrence is extremal if the blocks can not move inside the occurrence.
    This is a subset of all occurrences of ``u`` in ``w`` but *much* that
    defines the same convex hull.

    EXAMPLES::

        sage: extremal_occurrences('aaaaa', 'aaa')
        [(0, 1, 2), (0, 1, 4), (0, 3, 4), (2, 3, 4)]
        sage: extremal_occurrences('abababa', 'ab')
        [(0, 1), (0, 5), (4, 5)]
        sage: extremal_occurrences('aabaabaabaa', 'ab')
        [(0, 2), (1, 2), (0, 8), (7, 8)]

    Note the difference with `extremal_occurrences`::

        sage: letter_extremal_occurrences('aaaaa', 'aaa')
        [(0, 1, 2), (0, 2, 3), (1, 2, 3), (0, 1, 4), (1, 2, 4), (0, 3, 4), (2, 3, 4)]
        sage: letter_extremal_occurrences('abababa', 'ab')
        [(0, 1), (2, 3), (0, 5), (4, 5)]
        sage: letter_extremal_occurrences('aabaabaabaa', 'ab')
        [(0, 2), (1, 2), (4, 5), (0, 8), (7, 8)]

    Some more tests::

        sage: extremal_occurrences('aaabcdefghaaa', 'aa')
        [(0, 1), (2, 10), (0, 12), (11, 12)]
        sage: extremal_occurrences('bcdef', 'a')
        []
        sage: extremal_occurrences('bacadeaf', 'a')
        [(1,), (6,)]

    Some check for equality of polyhedra::

        sage: w = 'aaabbbaaabbbaaabbb'
        sage: u = 'aaab'
        sage: o0 = occurrences(w,u)
        sage: o1 = letter_extremal_occurrences(w,u)
        sage: o2 = extremal_occurrences(w,u)
        sage: len(o0), len(o1), len(o2)
        (315, 41, 25)
        sage: Polyhedron(o0) == Polyhedron(o1) == Polyhedron(o2)
        True
        sage: u = 'aabaa'
        sage: o0 = occurrences(w,u)
        sage: o1 = letter_extremal_occurrences(w,u)
        sage: o2 = extremal_occurrences(w,u)
        sage: len(o0), len(o1), len(o2)
        (270, 60, 42)
        sage: Polyhedron(o0) == Polyhedron(o1) == Polyhedron(o2)
        True

        sage: for w in ('aaaabbbb', 'abbaab', 'aabbaabb', 'abababababababab',
        ....:           'abbaaabbbaaaabbbaab',
        ....:           'aaaabbbbaaaabbbbaaaa'):
        ....:     for n in (1,2,3,4):
        ....:         for u in product('ab', repeat=n):
        ....:             P0 = Polyhedron(occurrences(w,u))
        ....:             P1 = Polyhedron(letter_extremal_occurrences(w,u))
        ....:             P2 = Polyhedron(extremal_occurrences(w,u))
        ....:             assert P0 == P1 == P2
    """
    if not u:
        return [()]

    # 1. Preprocessing of u: we compute the set of its factors inside a suffix
    # trie. To each factor of u corresponds a number in {0, 1, ..., n-1} where n
    # is the number of factors. The used variables are
    #  tf: the transition function (describe what happens when adding a letter)
    #  sl: the suffix link (describe what happens when removing the left letter)
    #  lengths: the length of factors
    alphabet = sorted(set(u))
    W = FiniteWords(alphabet)
    S = SuffixTrie(W(u, check=False))
    sl = S._suffix_link    # suffix link
    assert sl[0] == -1
    sl[0] = None           # replace -1 by None
    tf = {}                # transition function
    for (r,letter),s in S._transition_function.iteritems():
        tf[(r,letter[0])] = s
    if verbose:
        print "sl = {}".format(sl)
        print "tf = {}".format(tf)
    lengths  = [None]*len(sl)
    lengths[0] = 0      # lengths of the factors
    todo = [0]
    while todo:
        s0 = todo.pop()
        for letter in alphabet:
            t = (s0,letter)
            if t in tf:
                s1 = tf[t]
                lengths[s1] = lengths[s0] + 1
                todo.append(s1)
    if verbose:
        print "lengths = {}".format(lengths)


    # 2. run through w
    state = 0  # maximal suffix of u at the current position
    fact_occ = [-1] * len(sl)  # occurrences of factors of u
    pos = [[] for _ in range(len(u))]  # extremal subword occurrences of u
                                       # are stored as
                                       # [lb, (i0,j0,f0), (i1,j1,f1), ...]
                                       # where:
                                       #   lb: whether the last block is left
                                       #   blocked
                                       #   i0: starting point of a factor
                                       #   j0: endpoint of a factor
                                       #   f0: the index of the factor in the
                                       #       suffix trie
    for iw,letter in enumerate(w):
        if verbose:
            print "iw = {} letter={}".format(iw,letter)
            print "fact_occ = {}".format(fact_occ)
            print "pos"
            for ll in range(len(pos)):
                print "    {}: {}".format(ll,pos[ll])
        # a. compute the maximum state that can be extended as well as the new
        #    state after having read letter
        old_state = state
        new_state = tf.get((state,letter))
        while state and new_state is None:
            state = sl[state]
            new_state = tf.get((state,letter))
        if new_state is None:
            new_state = 0
        if verbose:
            print "state = {}, new_state = {}".format(state, new_state)

        # b. process subword occurrences
        for length in range(len(u)-1,0,-1):
            if verbose:
                print "  length={}".format(length)
            if letter == u[length]:
                k = 0
                while k < len(pos[length-1]):
                    x = pos[length-1][k]
                    i0, j0, s0 = x[-1]
                    j1 = x[-2][1] if len(x) > 2 else 0
                    lb = x[0]
                    if verbose: print "  x={}  lb={}  j1={}".format(x,lb,j1)
                    if j0 == iw:
                        # case when the last block can be continued
                        xx = x[:-1]
                        ss = tf[(s0,letter)]
                        xx.append((i0, j0+1, ss))
                        pos[length].append(xx)
                        xx[0] = xx[0] or fact_occ[ss] < j1 + lengths[ss]
                        k += 1
                        if verbose: print "  continue last block {}->{}".format(s0,ss)
                    elif lb or fact_occ[s0] == j0:
                        # case when the last block is either blocked on the left
                        # or on the right
                        xx = x[:]
                        ss = tf[(0,letter)]
                        xx.append((iw, iw+1, ss))
                        xx[0] = fact_occ[ss] < j0 + 1
                        pos[length].append(xx)
                        k += 1
                        if verbose: print "  new block {}".format(ss)
                    else:
                        if verbose: print "  delete x"
                        # we remove x
                        del pos[length-1][k]

        if letter == u[0]:
            ss = tf[(0,letter)]
            xx = [fact_occ[ss] == -1, (iw,iw+1,ss)]
            pos[0].append(xx)

        # update the last occurrences of factors of u and switch to the new
        # state
        state = new_state
        s = state
        while s is not None:
            fact_occ[s] = iw+1
            s = sl[s]
        if verbose: print

    return [sum((tuple(range(i,j)) for (i,j,k) in x[1:]),()) for x in pos[-1] if x[0] or fact_occ[x[-1][2]] == x[-1][1]]

def letter_extremal_occurrences(w, u):
    r"""
    Return the set of letter-extremal occurrences of ``u`` in ``w``.

    An occurrence is letter-extremal if the letters can not move inside the occurrence.
    This is a subset of all occurrences of ``u`` in ``w`` but *much* that
    defines the same convex hull.

    You should actually look at :func:`letter_extremal_occurrences` which is
    even smarter.

    EXAMPLES::

        sage: extremal_occurrences('aabb','ab')
        [(0, 2), (1, 2), (0, 3), (1, 3)]
        sage: extremal_occurrences('aaabbb','ab')
        [(0, 3), (2, 3), (0, 5), (2, 5)]

        sage: p = 'xyyxxxyyy'
        sage: s = 'xxxyyyxxy'
        sage: o1 = occurrences(p+'x'+s, 'xyx')
        sage: o1
        [(0, 1, 3), (0, 2, 3), (0, 1, 4), ...,   (11, 15, 17), (12, 15, 17)]
        sage: len(o1)
        138
        sage: o2 = extremal_occurrences(p+'x'+s, 'xyx')
        sage: o2
        [(0, 1, 3), (0, 2, 3), (5, 6, 9), ...,  (0, 15, 17), (12, 15, 17)]
        sage: len(o2)
        11
        sage: Polyhedron(o1) == Polyhedron(o2)
        True

        sage: for w in ('abbaab', 'aabbaabb', 'abbaaabbbaaaabbbaab',
        ....:           'aaaabbbbaaaabbbbaaaa'):
        ....:     for n in (1,2,3,4):
        ....:         for u in product('ab', repeat=n):
        ....:             P1 = Polyhedron(occurrences(w,u))
        ....:             P2 = Polyhedron(extremal_occurrences(w,u))
        ....:             assert P1 == P2
    """
    if len(u) == 0:
        return [()]

    pos = [[] for _ in range(len(u))]

    letters = set(w)
    n = len(w)

    # 1. find next left and next right occurrences of letters
    next_left = [None]*n
    next_right = [None]*n
    last = {letter: -1 for letter in letters}
    for i,letter in enumerate(w):
        next_left[i] = last[letter]
        last[letter] = i
    last = {letter: n for letter in letters}
    for i,letter in enumerate(reversed(w)):
        next_right[n-i-1] = last[letter]
        last[letter] = n-i-1

    # 2. run through w
    for i,letter in enumerate(w):
        for j in range(len(u)-1,0,-1):
            if letter == u[j]:
                for x in pos[j-1]:
                    if next_left[x[-1]] <= x[-2] or next_right[x[-1]] >= i:
                        pos[j].append(x + (i,))
            else:
                k = 0
                while k < len(pos[j-1]):
                    x = pos[j-1][k]
                    if next_left[x[-1]] > x[-2] and next_right[x[-1]] < i:
                        del pos[j-1][k]
                    else:
                        k += 1
        if letter == u[0]:
            pos[0].append((-1,i))

    return [x[1:] for x in pos[-1] if next_left[x[-1]] <= x[-2] or next_right[x[-1]] >= n]

def barycentric_coordinates(pts, q):
    r"""
    Return rational barycentric coordinates of ``q`` with respect to ``pts`` if
    they exist.

    If ``q`` is not in the convex hull of ``pts``, this function returns
    ``None``.

    EXAMPLES::

        sage: F = FreeModule(ZZ,3)
        sage: v0 = F((2,3,6))
        sage: v1 = F((2,3,18))
        sage: v2 = F((15,16,18))
        sage: m = F((8,9,13))
        sage: p = barycentric_coordinates([v0,v1,v2], m)  # not tested
        sage: p                                           # not tested
        [5/12, 19/156, 6/13]
        sage: p[0]*v0 + p[1]*v1 + p[2]*v2 == m            # not tested
        True

        sage: barycentric_coordinates([v0,v1], m) is None  # not tested
        True
    """
    if not pts:
        return
    N = len(pts)
    dim = len(q)
    pts = list(pts)
    assert all(len(p) == dim for p in pts)
    x = [Variable(i) for i in range(N)]
    cs = Constraint_System()
    cs.insert(sum(x[i] for i in range(N)) == 1)
    for i in range(N):
        cs.insert(x[i] >= 0)
    for i in range(dim):
        cs.insert(sum(x[j]*pts[j][i] for j in range(N)) == q[i])
    m = MIP_Problem(N, cs, 0)
    try:
        pt = m.optimizing_point()
    except ValueError:
        return
    return [QQ((pt.coefficient(x[i]),pt.divisor())) for i in range(N)]

def extremal_mid_occurrences(p, s, u):
    r"""
    Iterator through the occurrences of ``u`` in ``p*s`` that go through ``*``.

    TESTS::

        sage: for p,s in (('xyyxxxyyy','xxxyyyxxy'), ('xxxyyyxxy','xyxyxyxy')):
        ....:     for n in (1,2,3):
        ....:         for u in product('xy', repeat=n):
        ....:             w = p+'*'+s
        ....:             o1 = extremal_occurrences(w, u)
        ....:             o2 = list(extremal_mid_occurrences(p, s, u))
        ....:             o3,o4 = occurrences(w, u)
        ....:             assert Polyhedron(o3) == Polyhedron(o1)
        ....:             assert Polyhedron(o4) == Polyhedron(o2)
    """
    n = len(p)
    for i in range(len(u)):
        for o1 in extremal_occurrences(p, u[:i]):
            o1 = o1 + (n,)
            for o2 in extremal_occurrences(s, u[i+1:]):
                yield o1 + tuple(n+1+j for j in o2)

def ppl_polytope(pts):
    r"""
    Build a ppl polytope (i.e. a ``C_Polyhedron``).

    This seems to be twice faster as Sage function... there is something wrong
    in Sage class for polyhedra.
    """
    gs = Generator_System()
    for p in pts:
        gs.insert(point(Linear_Expression(p,0)))
    return C_Polyhedron(gs)

def is_sv_identity(p, s, d, prefix=(), skip_common_factors=True, status=False):
    r"""
    Check if ``(pxs, pys)`` is a B^{sv} identity in dimension ``d``.

    This method go through all subwords of length ``d-1`` and for each of them
    see whether some polytopes coincide.

    INPUT:

    - ``p``, ``s`` -- prefix and suffix for the identity

    - ``d`` -- dimension

    - ``prefix`` -- an optional prefix (mostly used for parallelization, see the
      function ``is_sv_identity_parallel`` below).

    - ``skip_common_factors`` -- (default is ``False``) whether to skip the
      computation of the polytope for common factors of ``p`` and ``s``

    - ``status`` -- if ``True``, then instead of returning a boolean, returns a
      pair ``(boolean, status_string)`` where ``status_string`` gives some
      details about the computation.

    OUTPUT:

    Either a boolean or a pair ``(boolean, status_string)`` if ``status=True``.

    EXAMPLES::

        sage: p = 'xxyyx'
        sage: s = 'xxyxy'
        sage: is_sv_identity(p, s, 3)
        True
        sage: is_sv_identity(p, s, 4)
        False
        sage: ans, info = is_sv_identity(p, s, 3, status=True)
        sage: print info    # only one factor tested!
        u = yy
        num ext. occ.: 4
        num faces    : 3
        num verts    : 3
        polytope computation in ...secs
        sage: ans, info = is_sv_identity(p, s, 3, skip_common_factors=False, status=True)
        sage: print info    # all factors are tested
        u = xx
        num ext. occ.: 6
        num faces    : 4
        num verts    : 4
        polytope computation in ...
        <BLANKLINE>
        u = xy
        num ext. occ.: 4
        num faces    : 4
        num verts    : 4
        polytope computation in ...
        <BLANKLINE>
        u = yx
        num ext. occ.: 4
        num faces    : 4
        num verts    : 4
        polytope computation in ...
        <BLANKLINE>
        u = yy
        num ext. occ.: 4
        num faces    : 3
        num verts    : 3
        polytope computation in ...
        <BLANKLINE>

        sage: p,s = vincent_sv_prefix_suffix(4)
        sage: is_sv_identity(p, s, 4)
        True
        sage: is_sv_identity(p, s, 5)
        False

        sage: p,s = vincent_sv_prefix_suffix(5)
        sage: is_sv_identity(p, s, 5)
        True
        sage: is_sv_identity(p, s, 6)
        False

        sage: x,y = symbolic_max_plus_matrices_band(5, 2, 's', 'v')
        sage: p = x*y*y*x*x*x*y*y*y*y*x*x*x*x  # ~0.5sec
        sage: s = y*y*y*y*x*x*x*x*y*y*y*x*x*y  # ~0.5sec
        sage: p*x*s == p*y*s                   # ~6secs
        True

        sage: p,s=('xxxyyxyxx', 'xxxyxxyxy')
        sage: is_sv_identity(p, s, 4)
        True
        sage: is_sv_identity(p, s, 5)
        False
    """
    if status:
        output = ''
    # compute the common factors of p and s
    if skip_common_factors:
        F = FiniteWords('xy')
        facts = set(tuple(f) for f in F(p).factor_set(d-1)).intersection(
                    (tuple(f) for f in F(s).factor_set(d-1)))

    # now iterate through all words `u` with given prefix
    pref = tuple(prefix)
    n = len(p)
    for q in product('xy', repeat=d-1-len(prefix)):
        u = prefix+q
        if skip_common_factors and u in facts:
            # ignore `u` that have a factor occurrence in both p and s
            continue

        # compute the polytope of occurrences and check that it contains the
        # "middle occurrences"
        avoid = extremal_occurrences(p+'*'+s, u)
        t0 = time()
        P = ppl_polytope(avoid)
        if status:
            output += 'u = {}\n'.format(''.join(u))
            output += 'num ext. occ.: {}\n'.format(len(avoid))
            output += 'num faces    : {}\n'.format(len(P.minimized_constraints()))
            output += 'num verts    : {}\n'.format(len(P.minimized_generators()))
            output += 'polytope computation in {}secs\n'.format(time()-t0)
            output += '\n'
        for o in extremal_mid_occurrences(p, s, u):
            pt = C_Polyhedron(point(Linear_Expression(o,0)))
            if not P.contains(pt):
                return (False,output) if status else False
    return (True,output) if status else True

def parallel_unfold(args):
    r"""
    Helper for parallel functions.
    """
    verbose = args[0]
    f = args[1]
    args = args[2:]
    if verbose:
        t = datetime.now()
        t0 = time()
        print "{}: new job at {}:{}:{}\n  {}".format(
                mp.current_process().name, t.hour, t.minute,
                t.second, args)
        sys.stdout.flush()
    ans = f(*args)
    if verbose:
        print "{}: job done in {} seconds".format(mp.current_process().name, time()-t0)
        sys.stdout.flush()
    return ans

def is_sv_identity_parallel(p, s, d, prefix_length, ncpus=None, verbose=False,
        skip_common_factors=True, logfile=None):
    r"""
    Check identity using parallelization features

    INPT:

    - ``p``, ``s`` -- the identity to check

    - ``d`` -- the dimension of the space

    - ``prefix_length`` -- length of the prefix

    - ``ncpus`` -- (optional, default the number of cpus on the computer) number
      of cpus to use

    - ``verbose`` -- (optional, default ``False``) whether some additional
      information about the workers will be printed

    - ``skip_common_factors`` -- (optional, default ``True``) whether to also
      test factors that are common to both ``p`` and ``s``

    - ``logfile`` -- can be a string that specifies a filename for writing
      information about each occurrence polytope or sys.stdout

    EXAMPLES::

        sage: p,s = vincent_sv_prefix_suffix(5)
        sage: is_sv_identity_parallel(p, s, 5, 3)  # not tested (fail in doctest)
        True

        sage: p,s = vincent_sv_prefix_suffix(6)
        sage: is_sv_identity_parallel(p, s, 6, 4) # not tested
        True

    The following will write its output in a file named 'logfile.txt'::

        sage: is_sv_identity_parallel(p, s, 6, 4, logfile='logfile.txt') # not tested

    And for live information you can either turn on the ``verbose`` option (for
    a per job information) or set the ``logfile`` option to ``sys.stdout`` (for
    a per polytope information)::

        sage: is_sv_identity_parallel(p, s, 6, 4, verbose=True)  # not tested
        sage: import sys
        sage: is_sv_identity_parallel(p, s, 6, 4, logfile=sys.stdout)  # not tested
    """
    close = False
    if isinstance(logfile, str):
        close = True
        logfile = open(logfile, 'w')
    if ncpus is None:
        ncpus = mp.cpu_count()
    pool = mp.Pool(ncpus)
    tasks = ((verbose,is_sv_identity,p,s,d,prefix,skip_common_factors,True) for prefix in product('xy', repeat=prefix_length))
    for ans,status in pool.imap_unordered(parallel_unfold, tasks):
        if logfile is not None:
            logfile.write(status)
            logfile.flush()
        if ans is False:
            break
    pool.terminate()
    pool.join()
    if close:
        logfile.close()
    return ans

def vincent_sv_prefix_suffix(d):
    r"""
    Return the ``p``,``s`` from the conjecture
    """
    p = ''
    s = ''
    pletter = 'x'
    sletter = 'y'
    for i in range(1,d):
        p = p + pletter * i
        s = sletter * i + s
        pletter,sletter = sletter,pletter
    p = p + pletter * (d-1)
    s = sletter * (d-1) + s
    return p,s



#######################
# Convex hull engines #
#######################

def get_convex_hull_engine(nvar, convex_hull=None):
    r"""
    Call to various library for convex hull computation

    EXAMPLES::

        sage: CH0 = get_convex_hull_engine(3, 'ppl_raw')
        sage: CH1 = get_convex_hull_engine(3, 'ppl')
        sage: CH2 = get_convex_hull_engine(3, 'cdd')
        sage: CH3 = get_convex_hull_engine(3, 'PALP')

        sage: F = ZZ**3
        sage: pts = [F.random_element() for _ in range(20)]
        sage: CH0(pts) == CH1(pts) == CH2(pts) == CH3(pts)
        True

        sage: CH0 = get_convex_hull_engine(4, 'ppl_raw')
        sage: CH1 = get_convex_hull_engine(4, 'ppl')
        sage: CH2 = get_convex_hull_engine(4, 'cdd')
        sage: CH3 = get_convex_hull_engine(4, 'PALP')
        sage: F = ZZ**4
        sage: pts = [F.random_element() for _ in range(30)]
        sage: CH0(pts) == CH1(pts) == CH2(pts) == CH3(pts)
        True

    By far the fastest is CH1 ('ppl') in all situtations encountered. Note that
    PALP is not reliable (errors and overflows for large input). See also:

    - https://www.inf.ethz.ch/personal/fukudak/polyfaq/node23.html
    - http://miconvexhull.codeplex.com/
    - http://www.boost.org/doc/libs/1_47_0/libs/geometry/doc/html/geometry/reference/algorithms/convex_hull.html
    - Fukada, "From the zonotope construction to the Minkowski addition of
      convex polytopes"
      http://www.sciencedirect.com/science/article/pii/S0747717104000409
    """
    if isinstance(convex_hull, ConvexHull):
        return convex_hull
    elif convex_hull is None or convex_hull == 'ppl_raw':
        return ConvexHullPPL(nvar)
    elif convex_hull == 'ppl':
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
    _name = 'none'
    def __eq__(self, other):
        return type(self) is type(other)

    def __ne__(self, other):
        return not self == ohter

    def __repr__(self):
        return "Convex hull engine {}".format(self._name)

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

    def __call__(self, pts):
        if pts:
            pts = [p.vector() for p in self._parent([pts,[],[]],None).vertex_generator()]
            for a in pts: a.set_immutable()
            pts.sort()
        return tuple(pts)


class ConvexHullPPL(ConvexHull):
    r"""
    Compute convex hull using a raw PPL interface (`sage.libs.ppl`).

    This is not faster than polyhedra from Sage.
    """
    _name = 'ppl_raw'
    def __init__(self, dim):
        self.dim = int(dim)
        self.vars = [Variable(i) for i in range(dim)]
        self.V = FreeModule(ZZ, self.dim)

    def __eq__(self, other):
        return type(self) is type(other) and self.dim == other.dim 

    def __call__(self, pts):
        if pts:
            gs = Generator_System()
            for p in pts:
                gs.insert(point(Linear_Expression(p,0)))
            poly = C_Polyhedron(gs)
            pts = [self.V(p.coefficients()) for p in poly.minimized_generators()]
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
        a = ans.split('\n')
        try:
            dd,nn = a[0].split(' ')[:2]
            dd = int(dd)
            nn = int(nn)
            if dd > nn: dd,nn = nn,dd
        except (TypeError,ValueError):
            raise RuntimeError("PALP got wrong:\n{}".format(ans))
        if d2 != int(dd):
            raise RuntimeError("dimension changed... have d={} but PALP answered dd={} and nn={}".format(d2,dd,nn))
        n2 = int(nn)
        coords = []
        for i in xrange(1,d2+1):
            coords.append(map(ZZ,a[i].split()))
        new_pts = zip(*coords)

        if d2 != d:
            new_pts = [pts[0] + V2(p)*M1 for p in new_pts]

        t = [self._free_module(p) for p in new_pts]
        for v in t:
            v.set_immutable()
        t.sort()
        return tuple(t)
