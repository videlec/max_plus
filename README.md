# Sage library for max-plus matrix identities

## Introduction

To run the library you need to install Sage (http://sagemath.org). Then clone
this repository and just install it as a standard Python package:

    $ sage -python setup.py install

To check that it works just start Sage and try to load `max_plus` as follows

    sage: import max_plus

If you do not get error message with the above line, then you are done.

## Symbolic and integer max plus matrices

The max plus symbolic module define a certain number of functions

    symbolic_max_plus_matrices_band(d, n, diag, surdiag, ch, typ)
	symbolic_max_plus_matrices_upper(d, n, diag, surdiag, ch, typ)
    symbolic_max_plus_matrices(d, n, ch, typ)
    symbolic_max_plus_identity(d, nvar, ch)


for respectively band, upper triangular and full matrices. The arguments are

- `d` - the dimension

- `n` - the number of matrices

- `diag` - describes whether the diagonal is zero (`'z'`), constant (`'c'`),
  the same in each of the `n` matrices (`'s'`) or variable (`'v'`). This is
  optional and default to `'v'`.

- `surdiag` - describes whether the surdiagonal is constant, zero, same or
  variable.  This is also optional.

- `ch` - (optional) set the convex hull engine (could be one of 'ppl', 'cdd'
  or 'PALP'). It defaults to 'ppl' which seems to be the fastest.

- `typ` - (optional, default to `sym`) specifies the implementation: if `'sym'`
  the matrices do not store all the coefficients (only two in the case
  of full matrices and `d` in the case of triangular ones). The multiplication
  in these case should be a little bit faster (by a factor `d^2` for dense
  matrix but much less for triangular ones). This should be set to `'full'` only
  to debug the code or to play with `diag=c` or `surdiag=c` that are not
  compatible with this option.

For example you can do

	sage: from max_plus import *

    sage: x,y = symbolic_max_plus_matrices_band(3, 2, 'z', 'v')
    sage: x
    A 3x3 symbolic max plus matrix on 4 variables
    sage: print x
	[   0  x0 -oo ]
	[ -oo   0  x1 ]
	[ -oo -oo   0 ]
	sage: y
    A 3x3 symbolic max plus matrix on 4 variables
	sage: print y
	[   0  x2 -oo ]
	[ -oo   0  x3 ]
	[ -oo -oo   0 ]
    sage: p = x*y*x
    A 3x3 symbolic max plus matrix on 4 variables
    sage: print p
	[   0 max(x2, x0) max(x1+x2, x0+x3, x0+x1) ]
	[ -oo           0              max(x3, x1) ]
	[ -oo         -oo                        0 ]
	sage: x*y*x*y == y*x*y*x == x*y*y*x == y*x*x*y
	True

	sage: x2,y2,z2 = symbolic_max_plus_matrices_upper(3, 3)
    sage: x2.num_vars()
    18
	sage: print x2
	[  x0   x3  x5 ]
	[ -oo   x1  x4 ]
	[ -oo  -oo  x2 ]
	sage: print y2
	[  x6   x9  x11 ]
	[ -oo   x7  x10 ]
	[ -oo  -oo   x8 ]
	sage: print z2
	[ x12  x15  x17 ]
	[ -oo  x13  x16 ]
	[ -oo  -oo  x14 ]

And with full 3x3 matrices

    sage: x,y,z = symbolic_max_plus_matrices(3, 3)
    sage: print x
	[ x0 x1 x2 ]
	[ x3 x4 x5 ]
	[ x6 x7 x8 ]
    sage: print y
	[  x9 x10 x11 ]
	[ x12 x13 x14 ]
	[ x15 x16 x17 ]
	sage: print x*x
	[   max(x2+x6, x1+x3, 2x0) max(x2+x7, x1+x4, x0+x1) max(x2+x8, x1+x5, x0+x2) ]
	[ max(x5+x6, x3+x4, x0+x3)   max(x5+x7, 2x4, x1+x3) max(x5+x8, x4+x5, x2+x3) ]
	[ max(x6+x8, x3+x7, x0+x6) max(x7+x8, x4+x7, x1+x6)   max(2x8, x5+x7, x2+x6) ]

From a symbolic matrix you can obtain an integer max plus matrix using the
method `eval`:: 

    sage: x,y = symbolic_max_plus_matrices_band(3, 2, 'z', 'v')
    sage: xv = x.eval((1,-1,0,3))
    sage: xv
    [   0   1 -oo ]
	[ -oo   0  -1 ]
	[ -oo -oo   0 ]
    sage: yv = y.eval((1,-1,0,3))
    sage: yv
    [   0   0 -oo ]
	[ -oo   0   3 ]
	[ -oo -oo   0 ]
    sage: p = x*y*x
    sage: pv = p.eval((1,-1,0,3))
    sage: pv
	[   0   1 4 ]
	[ -oo   0 3 ]
	[ -oo -oo 0 ]
    sage: pv == xv*yv*xv
    True

Computations with integer matrices are infinitely faster. Hence, if you intend
to test some relations it is adviced to first test them on a sample of integer
matrices.

## Combinatorics

There are currently some functions to check relations in `B^{sv}_d` and `B^{vv}_d`.

- `def occurrences(w, u)`: return the occurrences of `u` in `w`

- `def extremal_occurrences(w, u)`: return the 1-extremal occurrences of `u` in
  `w` (the convex hull of this set is the same as all occurrences)

- `def letter_extremal_occurrences(w, u)`: return the extremal occurrences of
  `u` in `w` (the convex hull of this set is the same as all occurrences)

- `def extremal_mid_occurrences(p, s, u)`: return the occurrences of `u` in
  `p*s` where the occurrence uses the joker letter `*`.

To check identities you can use the following functions (note that identities
must be written with the letters 'x' and 'y'):

- `is_sv_identity(p, s, d, prefix, skip_common_factors)`: check whether `(pxs,pys)` is an identity in `B^{sv}_d`.

- `is_sv_identity_parallel(p, s, d, prefix_length, ncpus, verbose,
  skip_common_factors, logfile)`: the same as above but with parallelization.
  The argument `prefix_length` is used to chunk the subwords into different
  jobs: if it is set to `k` then there will be `2^k` jobs which correspond to
  the `2^k` possible prefixes. The argument `logfile` can be set to either
  `sys.stdout` or a string in which case details about the execution will be
  written in a file with that name.

Here are some examples

	sage: from max_plus import *

    sage: p = 'xyxxyy'
    sage: s = 'xxyyxy'
    sage: is_sv_identity(p+'x'+s, p+'y'+s, 3)
    True

    sage: p = 'xyyxxxyyyyxxxx'
    sage: s = 'yyyyxxxxyyyxxy'
    sage: %time is_sv_identity(p+'x'+s, p+'y'+s, 5)
	CPU times: user 48 ms, sys: 12 ms, total: 60 ms
	Wall time: 55.7 ms
	True

	sage: p = 'xyyxxxyyyyxxxxxyyyyy'
	sage: s = 'xxxxxyyyyyxxxxyyyxxy'
	sage: %time is_sv_identity(p+'x'+s, p+'y'+s, 6)
	CPU times: user 276 ms, sys: 8 ms, total: 284 ms
	Wall time: 259 ms
	True
	True
    sage: is_sv_identity(p+'x'+s, p+'y'+s, 7)
    False

	sage: p,s = vincent_sv_prefix_suffix(7)
	sage: is_sv_identity_parallel(p+'x'+s, p+'y'+s, 7, 3, logfile=sys.stdout)
	u = xxxxyx
	num ext. occ.: 371
	num faces    : 16
	num verts    : 58
	polytope computation in 0.0998837947845secs

	u = xxxyxx
	num ext. occ.: 289
	num faces    : 11
	num verts    : 30
	polytope computation in 0.0397710800171secs
	...
	u = yxyyyy
	num ext. occ.: 371
	num faces    : 16
	num verts    : 58
	polytope computation in 0.0290420055389secs

	True


## Contact

Any comment or remark is welcome at vincentDOTdelecroixATlabriDOTfr
