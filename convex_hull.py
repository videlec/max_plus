r"""
Call to various library for convex hull computation

EXAMPLES::

    sage: CH1 = ConvexHullPolyhedra(3)
    sage: CH2 = ConvexHullPalp(3)
    sage: CH3 = ConvexHullPalp2(3)

    sage: F = ZZ**3
    sage: pts = [F.random_element() for _ in range(20)]
    sage: CH1(pts) == CH2(pts) == CH3(pts)
    True

By far the fastest is CH1 (indirectly ppl).

- http://miconvexhull.codeplex.com/
- http://www.boost.org/doc/libs/1_47_0/libs/geometry/doc/html/geometry/reference/algorithms/convex_hull.html
"""
from subprocess import Popen, PIPE

from sage.misc.misc import SAGE_TMP
from sage.misc.temporary_file import tmp_filename
from sage.misc.misc import SAGE_TMP

from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ

def get_convex_hull_engine(nvar, convex_hull=None):
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
        [(0, 0), (0, 1), (1, 1), (3, 4)]
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

class ConvexHullPalp2(ConvexHull):
    r"""
    Compute convex hull using PALP via LatticePolytope.

    Note: the points must be lattice points (ie integer coordinates) and
    generate the ambient vector space.
    """
    _name = 'PALP2'

    def __init__(self, dim):
        from sage.modules.free_module import FreeModule
        self._free_module = FreeModule(ZZ, dim)

    def __eq__(self, other):
        return type(self) is type(other) and self._free_module == other._free_module

    def __ne__(self, other):
        return type(self) is not type(other) or self._free_module != other._free_module

    def __call__(self, pts):
        from sage.geometry.lattice_polytope import LatticePolytope
        l = map(self._free_module, LatticePolytope(pts).vertices())
        l.sort()
        for v in l: v.set_immutable()
        return tuple(l)

#############
# old stuff #
#############
def convex_hull(points):
    r"""
    Return the convex hull of points in the list ``points``.

    It seems that there is a bug in qhull if the number of points is smaller
    than the dimension
    """
    if not points:
        return points

    filename = tmp_filename()
    n_pts = len(points)
    dim   = len(points[0])

    if n_pts <= 2:
        return points
    elif n_pts <= dim:
        from sage.structure.element import parent
        from sage.geometry.polyhedron.constructor import Polyhedron
        P = parent(points[0])
        return map(P,Polyhedron(points).vertices_list())

    # we first process the points so that the space they generate is full
    # dimensional

    from sage.modules.free_module import FreeModule
    from sage.rings.integer_ring import ZZ

    V = FreeModule(ZZ, dim)
    W = V.submodule([u-v for u in points for v in points])
    if W.dim() < dim:
        # base change
        print "the points belong to a smaller affine subspace"
        return

    f = open(filename, 'w')
    f.write('{}\n'.format(dim))
    f.write('{}\n'.format(n_pts))
    template = ' '.join('{}' for _ in range(dim)) + '\n'
    for p in points:
        f.write(template.format(*p))
    f.close()

    cmd = ['qhull'.format(filename), 'Fx']
    try:
        qhull_proc = Popen(cmd, stdin=open(filename), stdout=PIPE, stderr=PIPE)
    except OSError:
        raise ValueError("The qhull must be installed (type "
                "'sage -i qhull') in a console or "
                "'install_package('qhull') at a Sage prompt)!\n")

    ans, err = qhull_proc.communicate()
    ret_code = qhull_proc.poll()
    if ret_code:
        print ans,err
        raise RuntimeError("qhull returned {}".format(ret_code))

    ans = ans.splitlines()
    return [points[int(ans[i+1])] for i in range(int(ans[0]))]

def plot_random_convex_hull(n):
    pts = [(randint(-10,10),randint(-10,10)) for _ in range(n)]
    ch = convex_hull(pts)
    pts = [x for x in pts if x not in ch]
    return point(pts, color='blue') + point(ch, color='red')


