"""
Some experiments on symbolic max plus matrices

We consider semigroup relations for various max-plus matrix semigroups:

 * upper triangular matrices with identical diagonals ('tri_sym_diag')
 * upper triangular matrices ('tri')
 * all matrices ('all')

If u = v is a relation then

pus = pvs is also a relation...

So we should look for primitive relations.
"""

from sage.geometry.polyhedron.parent import Polyhedra
from semigroup_tools import products, is_relation

###################
# Relation finder #
###################

def relations_tri(dim, num_mat=10, filename=None):
    r"""
    List the relations for the upper triangular matrices in a given dimension.

    INPUT:

    - ``dim`` -- the dimension

    - ``num_mat`` -- the number of integer matrices used to check the relation.
      Note that if all the integer matrices satisfy the relations, a (costly)
      symbolic check is performed to guarantee that the relation is satisfied by
      any pair of matrices.

    EXAMPLES::

        sage: relations_tri(2)
        xyyx(yx)xyyx = xyyx(xy)xyyx
        xyyx(yx)yxxy = xyyx(xy)yxxy
        xxyyx(yx)xyyxy = xxyyx(xy)xyyxy
        xxyyx(yx)xyyyx = xxyyx(xy)xyyyx
        xxyyx(yx)yxxyy = xxyyx(xy)yxxyy
        xxyyyx(yx)xyyx = xxyyyx(xy)xyyx
        ...
    """
    # the symbolic max-plus matrices
    a,b = symbolic_max_plus_matrices_tri(dim, 2)
    one = symbolic_max_plus_identity(dim, a.num_variables())

    # the integer max-plus matrices
    mats = [random_integer_max_plus_matrix_tri(dim, -50*dim*dim, 50*dim*dim) for _ in range(num_mat)]
    one_int = integer_max_plus_matrix_identity(dim)

    # n     : length of the product
    # i1,m1 : data for the first word
    # i2,m2 : data for the second word
    n = 1
    while True:
        if filename is None:
            from sys import stdout
            f = stdout
        else:
            f = open(filename.format(dim,2*n), 'w')

        for i1,m1 in products(mats[0], mats[1], n-1, n, one_int):
            m1 = mats[0] * m1

            # look at the relations x* = x*
            for i2,m2 in products(mats[0], mats[1], n-1, n, one_int):
                if i1 == i2:
                    break
                m2 = mats[0] * m2
                # here we first test equality between m1 and m2
                # then we test the relations on all matrices in mats
                # then we check formally using symbolic matrices
                if m1 == m2 and \
                   is_relation([0]+i1, [0]+i2, mats) and \
                   a * prod(a if x == 0 else b for x in i1) == a * prod(a if x == 0 else b for x in i2):
                       p = 0
                       while i1[p] == i2[p]: p += 1
                       s = -1
                       while i1[s] == i2[s]: s -= 1
                       s += 1
                       f.write('x{}({}){} = x{}({}){}\n'.format(
                           ''.join('x' if x == 0 else 'y' for x in i1[:p]),
                           ''.join('x' if x == 0 else 'y' for x in i1[p:s]),
                           ''.join('x' if x == 0 else 'y' for x in i1[s:]),
                           ''.join('x' if x == 0 else 'y' for x in i2[:p]),
                           ''.join('x' if x == 0 else 'y' for x in i2[p:s]),
                           ''.join('x' if x == 0 else 'y' for x in i2[s:])))

                       f.flush()

            # look at the relations x* = y*
            #for i2,m2 in products(mats[0], mats[1], n, n-1, one_int):
            #    m2 = mats[1] * m2
            #    if m1 == m2 and \
            #       is_relation([0]+i1, [1]+i2, mats) and \
            #       a * prod(a if x == 0 else b for x in i1) == a * prod(a if x == 0 else b for x in i2):
            #           f.write('x{} = y{}\n'.format(
            #               ''.join('x' if x == 0 else 'y' for x in i1),
            #               ''.join('x' if x == 0 else 'y' for x in i2)))
            #           f.flush()

        if filename is not None:
            f.close()
        n += 1

##############################
# Symbolic max-plus matrices #
##############################

def symbolic_max_plus_identity(dim, N):
    r"""
    Return the ``dim x dim`` identity matrices on ``N`` variables.
    """
    P = Polyhedra(ZZ, N)
    e = P([[],[],[]], None)
    zero = P()

    data = [[zero if i == j else e for j in range(dim)] for i in range(dim)]
    return BinaryMaxPlusMatrix(dim, data)

def symbolic_max_plus_matrices(dim, nb):
    r"""
    Return ``nb`` independent symbolic matrices in dimension ``dim``
    """
    from sage.geometry.polyhedron.parent import Polyhedra

    N = nb * dim * dim
    P = Polyhedra(ZZ, N)
    e = P([[],[],[]], None)
    o = 0
    z = [0]*N
    matrices = []
    for i in range(nb):
        mat = []
        for j in range(dim):
            row = []
            for k in range(dim):
                z[o] = 1
                o += 1
                row.append(P([[z],[],[]],None))
                z[o-1] = 0
            mat.append(row)
        matrices.append(BinaryMaxPlusMatrix(dim, mat))
    return matrices

def symbolic_max_plus_matrices_tri(dim, nb):
    r"""
    Return ``nb`` independent triangular matrices in dimension ``dim``.

    INPUT:

    - ``dim`` -- dimension

    - ``nb`` -- number of matrices
    """
    from sage.geometry.polyhedron.parent import Polyhedra

    N = nb * dim*(dim+1)//2
    P = Polyhedra(ZZ, N)
    e = P([[],[],[]], None)
    o = 0
    z = [0]*N
    matrices = []
    for i in range(nb):
        mat = []
        for j in range(dim):
            row = [e] * j
            for k in range(dim-j):
                z[o] = 1
                o += 1
                row.append(P([[z], [], []], None))
                z[o-1] = 0
            mat.append(row)
        matrices.append(BinaryMaxPlusMatrix(dim, mat))
    return matrices

def symbolic_max_plus_matrices_tri_sym_diag(dim, nb):
    r"""
    Return a set of ``nb`` symbolic matrices with the same diagonal.

    INPUT:

    - ``dim`` -- the dimension

    - ``nb`` -- the number of matrices
    """
    N = nb * dim*(dim-1)//2 + dim
    P = Polyhedra(ZZ, N)
    e = P([[],[],[]], None)

    o = 0
    z = [0]*N
    diag = []
    for k in range(dim):
        z[o] = 1
        o += 1
        diag.append(P([[z], [], []], None))
        z[o-1] = 0

    matrices = []
    for i in range(nb):
        mat = []
        for j in range(dim):
            row = [e] * j
            row.append(diag[j])
            for k in range(dim-j-1):
                z[o] = 1
                o += 1
                row.append(P([[z], [], []], None))
                z[o-1] = 0
            mat.append(row)
        matrices.append(BinaryMaxPlusMatrix(dim, mat))
    return matrices


class BinaryMaxPlusMatrix(SageObject):
    r"""
    A symbolic max plus matrix.
    """
    def __init__(self, n, data):
        self._n = n
        if len(data) != n or any(len(x) != n for x in data):
            raise ValueError
        self._data = tuple(tuple(x) for x in data)

    def num_variables(self):
        r"""
        Return the number of variables used in this matrix.
        """
        return self._data[0][0].ambient_space().dimension()

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
        return (self._n == other._n and
                all(self._data[i][j] == other._data[i][j] \
                        for i in range(self._n) for j in range(self._n)))

    def __ne__(self, other):
        r"""
        Difference test
        """
        return not self.__eq__(other)

    def __mul__(self, other):
        r"""
        Multiplication of matrices
        """
        n = self._n
        assert self._n == other._n
        new_data = [[None]*n for _ in range(n)]
        for i in range(n):
            for j in range(n):
                vertices = []
                for k in range(n):
                    vertices.extend([xi+yi for xi,yi in zip(x,y)] \
                            for x in self._data[i][k].vertices_list() \
                            for y in other._data[k][j].vertices_list())
                new_data[i][j] = self._data[0][0].parent()([vertices,[],[]],None)

        return BinaryMaxPlusMatrix(self._n, new_data)

    def _repr_(self):
        r"""
        String when the object is printed
        """
        return 'A {}x{} max plus matrix'.format(self._n,self._n)

    def __getitem__(self, data):
        r"""
        To obtain an item of the matrix
        """
        i,j = data
        return self._data[i][j].vertices_list()

