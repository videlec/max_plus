
# These routines compute the rank of a matrix with respect to classical linear
# spaces: rank 1 matrices (Barvinok rank), symmetric rank 1 matrices (symmetric
# Barvinok rank), and star tree matrices (star tree rank). In addition, it can
# compute a lower bound on the rank as the chromatic number of the deficiency
# graph for all of these notions as well as with respect to the tropical
# Grassmannian G(2,n) (tree rank).

# The following external programs are used to do all of the real work:
# For Develin's algorithm: (exact rank with respect to classical linear space)
#   Polymake: http://www.math.tu-berlin.de/polymake/
#   lp_solve: http://lpsolve.sourceforge.net/
# For the chromatic number of the deficiency graph: (--chromatic option)
#   smallk: http://www.cs.ualberta.ca/~joe/Coloring/Colorsrc/smallk.html

import os
from subprocess import Popen, PIPE
from random import randint

# dirty hack to make polymake available on several configurations
POLYMAKE_CMD = os.getenv('POLYMAKE_CMD')
if POLYMAKE_CMD is None:
    paths_to_test = ["/opt/polymake-3.0/perl/polymake",
                     "/home/merlet.g/polymake-3.0/perl/polymake"]
    for path in paths_to_test:
        if os.path.isfile(path):
            POLYMAKE_CMD = path
            break

if POLYMAKE_CMD is None:
    raise ImportError("polymake not available. Please set the environment variable POLYMAKE_CMD")

# Note: Internally, everything is in the max-plus algebra

# Utility functions
#
# Reading a matrix and making sure it has the right form.

# Returns true if the given array of arrays is a valid matrix
def is_matrix(matrix):
    ncols = len(matrix[0])
    for row in matrix[1:]:
        if len(row) != ncols:
            return False 
    return True

# Returns true if the given array of arrays is a valid symmetric matrix
def is_symmetric_matrix(matrix):
    if not is_matrix(matrix):
        return False

    n = len(matrix)
    if len(matrix[0]) != n:
        return False  # Not a square matrix

    for i in range(n):
        for j in range(i):
            if matrix[i][j] != matrix[j][i]:
                return False
    return True

# Develin's algorithm
#
# Develin's algorithm takes a matrix and creates a point configuration. For each
# entry in the matrix, the corresponding point consists of the coefficients in
# the parametrization and the value of the matrix in that entry. Then, the lower
# convex hull is computed and the minimal collection of factes which includes
# every vertex is the rank

# encodes a symmetric for Develin's algorithm with respect to symmetric rank 1
# (diag = True) or star tree matrices (diag = False)
def matrix_to_coords(matrix, diag):
    coords = []
    n = len(matrix)
    for i in range(n):
        for j in range(i):
            exp = [int(k == i or k == j) for k in range(n)]
            coords.append([matrix[i][j]] + exp)
        if diag:
            exp = [2 * int(k == i) for k in range(n)]
            coords.append([matrix[i][i]] + exp)
    return coords

# Encodes a not necessarily symmetric matrix in order to compute the Barvinok
# rank, i.e. with respect to rank 1 matrices
def encode_barvinok(matrix):
    coords = []
    n = len(matrix)
    for i in range(n):
        for j in range(n):
            erow = [int(k == i) for k in range(n)]
            ecol = [int(k == j) for k in range(n)]
            coords.append([matrix[i][j]] + erow + ecol)
    return coords

# Computes the lower convex hull of a collection of points. Returns a list of
# facets, with for each facet a list of the indices of the points in that facet.
# "Lower" is taken to mean with respect to the first coordinate.
def lower_convex_hull(coords):
    fh = open("tmp-rank", "w")
    fh.write("POINTS\n")
    for coord in coords:
        fh.write("1 " + " ".join([str(x) for x in coord]) + "\n")
    # This is a point at infinity which means that we get the lower convex hull
    fh.write("0 1" + "".join([" 0" for x in coords[0][1:]]) + "\n")
    fh.close()

    # TODO: open polymake only once!!
    p = Popen([POLYMAKE_CMD, "tmp-rank", "POINTS_IN_FACETS"], stdout=PIPE,
            stderr=PIPE)
    fh = p.stdout
    if fh.readline() != "POINTS_IN_FACETS\n":
        stderr.write("Wrong first line")
        return []

    facets = []
    line = fh.readline()
    while line and line != "\n":
        if line[0] != "{" or line[-2:] != "}\n":
            stderr.write("Bad list of points\n")
            return []

            # list of points in our facet
        facet = [int(x) for x in line[1:-2].split()]
        if max(facet) < len(coords):
            # The facet doesn't contain any of the points in our "ceiling" so
            # it's on the lower side
            facets.append(facet)

        line = fh.readline()
    p.communicate()
    fh.close()
    return facets

def lower_convex_hull_bis(coords):
    from sage.geometry.polyhedron.constructor import Polyhedron

    n = len(coords[0])
    coords = [(1,) + tuple(x) for x in coords]
    infty = (0,1) + (0,)*(n-1)
    coords.append(infty)
    coord_to_index = {v:i for i,v in enumerate(coords)}

    facets = []
    for facet in Polyhedron(coords).faces(n-1):
        if infty not in facet.as_polyhedron():
            facets.append([coord_to_index[tuple(j)] for j in facet.as_polyhedron().vertices_list()])

    return facets

# Computes the minimum number of elements of the array sets such that each
# non-negative integer less than npts is in one of these. Thus, sets is an array
# of arrays of integers in [0, npts). Returns False if there is no such
# collection because the union of the sets does not contain every integer in [0,
# npts).
def min_cover(npts, sets):
    # check if the problem is solvable
    covered = [False for i in range(npts)]
    for set in sets:
        for ndx in set:
            covered[ndx] = True
    for c in covered:
        if not c:
            return False

    # Write the Free MPS format integer programming problem
    fh = open("tmp-rank-fmps", "w")
    fh.write("NAME min cover\n")
    fh.write("ROWS\n") # constraints
    for i in range(npts):
        fh.write(" G POINT" + str(i) + "\n")
    fh.write(" N NUMSETS\n")
    fh.write("COLUMNS\n") # variables
    for i in range(len(sets)):
        fh.write("  SET%d NUMSETS 1\n" % i)
        for point in sets[i]:
            fh.write("  SET%d POINT%d 1\n" % (i, point))
    fh.write("RHS\n") # right hand side to constraints
    for i in range(npts):
        fh.write("  COVER POINT%d 1\n" % i)
    fh.write("BOUNDS\n") # bounds on variables
    for i in range(len(sets)):
        fh.write(" BV A SET%d\n" % i)
    fh.write("ENDATA\n")
    fh.close()

    # Run the solver
    if False:
        # GLPK solver
        os.system("glpsol -w tmp-rank-glpsol tmp-rank-fmps > /dev/null")
        fh = open("tmp-rank-glpsol")
        fh.readline()
        line = fh.readline()
        min = int(line.split()[1])
        fh.close()
    else:
        # lpsolve solver
        p = Popen(["lp_solve", "-fmps", "tmp-rank-fmps"], stdout=PIPE)
        fh = p.stdout
        fh.readline() # blank
        line = fh.readline()
        start = "Value of objective function: "
        if line[0:len(start)] != start:
            stderr.write("Unexpected output from lp_solve\n")
            min = 0
        else:
            min = int(line[len(start):])

        # read to the end of the output without storing anything
        while line:
            line = fh.readline()
        p.communicate()
        fh.close()
    return min

# This returns the rank of the point/classical linear space pair encoded in the
# matrix coords
def toric_rank(coords):
    return min_cover(len(coords), lower_convex_hull(coords))

# Computes the symmetric Barvinok rank (or star tree rank if star_tree_rank is
# True) of the matrix matrix.
def symmetric_rank(matrix, star_tree_rank):
    return toric_rank(matrix_to_coords(matrix, not star_tree_rank))

# Randomization and sampling routines
#
# These functions are used to generate random matrices and sample their ranks

def symmetrize(matrix):
    n = len(matrix)
    for i in range(n):
        for j in range(i):
            matrix[i][j] = matrix[j][i]

# generates a random symmetric n-by-n matrix
def random_sym_matrix(n, bound = None):
    matrix = random_matrix(n, bound)
    symmetrize(matrix)
    return matrix

# generates a random n-by-n matrix
def random_matrix(n, bound = None):
    if bound is None: bound = 2**30 - 1
    return [[randint(0, bound) for j in range(n)] for i in range(n)]

# compute the star tree rank of samples random n by n matrices and return the
# distribution of their ranks
def sample_ranks(samples, n):
    ranks = [0 for i in range(n-1)]
    for i in range(samples):
        rank = symmetric_rank(random_sym_matrix(n), True)
        ranks[rank] += 1
        stdout.write("\r%d %d" % (i, rank))
    stdout.write("\n")
    return ranks

# Deficiency graph
#
# These routines are used to compute the deficiency graph with respect to many
# different varieties and then to compute its chromatic number. A graph is
# encoded as a dictionary, in which (i,j) maps to True whenever (i,j) are an
# edge in the graph. The order doesn't matter: it is arbitrary which of (i,j)
# and (j,i) is included. The nodes are consecutive integers beginning with 0,
# and the number of nodes is always implicit.

# Translates the unordered pair (i,j) into an index
def sym_pair_ndx(i, j):
    if i < j:
        (i, j) = (j, i)
    return (i * i + i) / 2 + j

# Deficiency graph of 2x2 minors of a symmetric matrix
def minors_def_graph(matrix):
    n = len(matrix)
    edges = {}
    for i in range(n):
        for j in range(i):
            for k in range(i + 1):
                for l in range(k):
                    # two terms of the ij x kl minor
                    a = matrix[i][k] + matrix[j][l]
                    b = matrix[i][l] + matrix[j][k]
                    if a > b:
                        edges[(sym_pair_ndx(i, k), sym_pair_ndx(j, l))] = True
                    elif b > a:
                        edges[(sym_pair_ndx(i, l), sym_pair_ndx(j, k))] = True
    return (edges, n*(n+1)/2)

# Deficiency graph of off-diagonal 2x2 minors of a symmetric matrix
def minors_skew_def_graph(matrix):
    n = len(matrix)
    edges = {}
    for i in range(n):
        for j in range(i):
            for k in range(j):
                for l in range(k):
                    a = matrix[i][j] + matrix[k][l]
                    b = matrix[i][k] + matrix[j][l]
                    c = matrix[i][l] + matrix[j][k]
                    # tropical equations are a + b, a + c, b + c
                    if a > b or a > c:
                        edges[(pair_ndx(i, j), pair_ndx(k, l))] = True
                    if b > a or b > c:
                        edges[(pair_ndx(i, k), pair_ndx(j, l))] = True
                    if c > a or c > b:
                        edges[(pair_ndx(i, l), pair_ndx(j, k))] = True
    return (edges, n*(n-1)/2)

# Returns the deficiency graph of a non-nec. symmetric matrix with respect to
# the 2x2 minors
def minors_nonsym_def_graph(matrix):
    n = len(matrix)
    edges = {}
    for i in range(n):
        for j in range(i):
            for k in range(n):
                for l in range(k):
                    # ij x kl minor
                    a = matrix[i][l] + matrix[j][k]
                    b = matrix[i][k] + matrix[j][l]
                    if a > b:
                        edges[(n*i + l, n*j + k)] = True
                    elif b > a:
                        edges[(n*i + k, n*j + l)] = True
    return (edges, n*n)

# Computes the deficiency graph with respect to the 3 term Plucker relations.
# Unlike the above functions, each node is coded as a pair (i1, i2) with i1 >
# i2. Thus, an edge consists of a pair of pairs.
def tree_rank_graph(matrix):
    n = len(matrix)

    # find edges in the graph
    edges = {}
    for i in range(n):
        for j in range(i):
            for k in range(j):
                for l in range(k):
                    t1 = matrix[i][j] + matrix[k][l]
                    t2 = matrix[i][k] + matrix[j][l]
                    t3 = matrix[i][l] + matrix[j][k]
                    if t1 > t2 and t1 > t3:
                        edges[((i, j), (k, l))] = True
                    elif t2 > t1 and t2 > t3:
                        edges[((i, k), (j, l))] = True
                    elif t3 > t1 and t3 > t2:
                        edges[((i, l), (j, k))] = True

    return edges

# returns an index for the pair (i, j) with i assumed to be greater than j
def pair_ndx(i, j):
    return (i * i - i) / 2 + j

# Flattens the node names returned by tree_rank_graph, so that the node is just
# a single index instead of a pair.
def flatten_nodes(graph):
    flat_graph = {}
    for ((i, j), (k, l)) in graph.keys():
        flat_graph[(pair_ndx(i, j), pair_ndx(k, l))] = True
    return flat_graph

def flat_tree_graph(matrix):
    n = len(matrix)
    return (flatten_nodes(tree_rank_graph(matrix)), n*(n-1)/2)

# Writes a graph in DIMACS format. The nodes are 0 through num_nodes-1 and the
# edges (i, j) are the keys of graph
def write_graph(filename, num_nodes, graph):
    fh = open(filename, "w")
    fh.write("p edge %d %d\n" % (num_nodes, len(graph)))
    for i in range(num_nodes):
        fh.write("n %d 1\n" % (i + 1))
    for (i, j) in graph.keys():
        fh.write("e %d %d\n" % (i + 1, j + 1))
    fh.close()

# Runs the smallk graph coloring program and returns whether or not the graph
# contained in filename is k-colorable
def run_smallk(filename, k):
    seed = randint(0, 2**30 - 1)

    p = Popen(["smallk", filename, str(seed), str(k)], stdout=PIPE)
    fh = p.stdout
    l = fh.readline()
    while l:
        if l == "Coloring Verified\n":
            fh.close()
            return True
        elif l == "Coloring failed to find %d-coloring\n" % k:
            fh.close()
            return False
        l = fh.readline()
    stderr.write("Failed to find result in smallk output file\n")
    fh.close()
    return False

