#!/usr/bin/python

import sage.all

from optparse import OptionParser
from max_plus.rank import *

# Read a symmetric matrix from the stream fh
def read_matrix(fh):
    l = fh.readline()
    matrix = [[int(x) for x in l.split()]]
    for j in range(1, len(matrix[0])):
        l = fh.readline()
        matrix.append([int(x) for x in l.split()])
    return matrix

# the value of type is a triple of: a boolean (must the matrix be
# symmetric), a function for encoding the matrix for Develin's algorithm,
# and a function for computing the deficiency graph.
parser = OptionParser(usage="%prog [options] <matrix file>")
parser.set_defaults(type=(True,
        lambda matrix: matrix_to_coords(matrix, True), minors_def_graph))
parser.add_option("--star-tree", "-s", action="store_const", dest="type",
        const=(True, lambda matrix: matrix_to_coords(matrix, False), 
                minors_skew_def_graph),
        help="Compute the star tree rank")
parser.add_option("--barvinok", "-b", action="store_const", dest="type",
        const=(False, encode_barvinok, minors_nonsym_def_graph),
        help="Compute the Barvinok rank")
parser.add_option("--tree", "-t", action="store_const", dest="type",
        const=(True, "No algorithm for computing tree rank",
                flat_tree_graph),
        help="Compute the tree rank")
parser.add_option("--min-plus", "-m", action="store_true", dest="min_plus",
        help="Use the min-plus algebra instead of max-plus")
parser.add_option("--chromatic", "-c", action="store_true",
        dest="chromatic", help="Find deficiency graph chromatic number")
(options, args) = parser.parse_args()
if len(args) > 1:
    parser.error("Too many arguments")

# now dispatch to methods
# it just applies coords_fn

if len(args) == 1:
    fh = open(args[0], "r")
    matrix = read_matrix(fh)
    fh.close()
else:
    matrix=read_matrix(stdin)
    stdout.write("Read %dx%d matrix\n" % (len(matrix), len(matrix)))

if options.min_plus:
    matrix = [[-x for x in row] for row in matrix]

(symmetric, coords_fn, graph_fn) = options.type

# Check that the matrix is the right form
if not is_matrix(matrix):
    stderr.write("Not a matrix\n")
    exit(1)
if symmetric and not is_symmetric_matrix(matrix):
    stderr.write("Not a symmetric matrix\n")
    exit(1)

if options.chromatic:
    # compute the chromatic number of the deficiency graph
    (graph, num_edges) = graph_fn(matrix)
    #n = len(matrix)
    write_graph("tmp-graph", num_edges, graph)
    for i in range(1, 9): # smallk only supports up to 8 colors
        if run_smallk("tmp-graph", i):
            print i
            break
    else:
        stderr.write("Couldn't find any coloring\n")
else:
    # compute rank exactly using Develin's algorithm
    if isinstance(coords_fn, str):
        # not supported
        stderr.write(coords_fn + "\n")
    else:
        print toric_rank(coords_fn(matrix))
