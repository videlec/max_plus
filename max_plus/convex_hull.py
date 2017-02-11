r"""
"""
from sage_import import *
from subprocess import Popen, PIPE

def clean(pts, k):
    r"""
    Return a list of points with the same convex hull as in ``pts`` but with
    potentially less points.

    INPUT:

    - ``pts`` -- list or tuple of points
    - ``k``   -- number of nonredundant point set in the cleaning
                (in the case of matrices, a wise choice is nvars - dim + 1)
    """
    dim = len(pts[0])
    if len(pts) < 2*k:
        return pts

    F = FreeModule(RDF, dim)
    pts = [(p,F(p)) for p in pts]

    go = 1
    approx_ch = set()
    while go < 3:
        print "new loop, {} points".format(len(pts))

        T = set()
        n = 0
        for _ in range(k*k):
            z = F.random_element()
            m_min = m_max = z.dot_product(pts[0][0])
            p_min = p_max = pts[0][0]
            for p,pa in pts[1:]:
                m = z.dot_product(pa)
                if m > m_max:
                    p_max = p
                    m_max = m
                elif m < m_min:
                    p_min = p
                    m_min = m
            T.add(p_min)
            T.add(p_max)
            n += 1
        approx_ch.update(T)
        if len(approx_ch) > 2*k:
            T = sample(list(approx_ch), k)
        else:
            while len(T) < k:
                T.update(x[0] for x in sample(pts,k-len(T)))

        gs = Generator_System()
        for p in T:
            gs.insert(point(Linear_Expression(p,0)))
        poly_T = C_Polyhedron(gs)

        new_pts = []
        for p,pa in pts:
            if p in T:
                new_pts.append((p,pa))
                continue

            gs = Generator_System()
            gs.insert(point(Linear_Expression(p,0)))
            poly_p = C_Polyhedron(gs)

            if poly_T.contains(poly_p):
                go = 0
            else:
                new_pts.append((p,pa))

        go += 1

        pts = new_pts

    print "approx ch with {} points".format(len(approx_ch))
    return [x[0] for x in pts]

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

#######################
# Convex hull engines #
#######################

def get_convex_hull_engine(nvar, convex_hull=None):
    r"""
    Call to various library for convex hull computation

    EXAMPLES::

        sage: from max_plus.convex_hull import get_convex_hull_engine

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
    elif convex_hull == 'nothing':
        return NoConvexHull()
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

class NoConvexHull(ConvexHull):
    _name = 'nothing'
    def __call__(self, pts): return tuple(sorted(set(pts)))

class ConvexHullPolyhedra(ConvexHull):
    r"""
    EXAMPLES::

        sage: from max_plus.convex_hull import ConvexHullPolyhedra

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

def in_convex_hull_LP(pts, q, solver='ppl'):
    r"""
    Check whether ``q`` belongs to the convex hull of ``pts``.

    Note we use a slightly modified version of

    https://www.inf.ethz.ch/personal/fukudak/polyfaq/node22.html

    Namely

    maximize z^T q - z_0
    subject to
      z^T p_i - z_0 <= 0
      z^T q - z_0 <= 1
    
    In our algorithm we set z0 = 1

    EXAMPLES::

        sage: from max_plus.convex_hull import in_convex_hull_LP

        sage: F = FreeModule(ZZ,5)
        sage: for _ in range(100):
        ....:     pts = [F.random_element() for _ in range(20)]
        ....:     c = [ZZ.random_element(0,100) for _ in range(20)]
        ....:     z = sum(i*p for i,p in zip(c,pts)) / sum(c)
        ....:     assert in_convex_hull_LP(pts, z)
        ....:     m = max(p[0] for p in pts)
        ....:     assert not in_convex_hull_LP(pts, [m+1]+[0]*19)
    """
    if not pts:
        return False

    d = len(pts[0])
    M = MixedIntegerLinearProgram(solver=solver)
    z = M.new_variable(real=True)

    for p in pts:
        M.add_constraint(M.sum(z[i]*p[i] for i in range(d)) <= 1)

    M.set_objective(M.sum(z[i]*q[i] for i in range(d)))
    M.add_constraint(M.sum(z[i]*q[i] for i in range(d)) <= 2)

    return M.solve() <= 1
