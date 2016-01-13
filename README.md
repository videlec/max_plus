# Sage library for max-plus matrix identities

## Introduction

To run the library you need to install Sage (http://sagemath.org). Then clone
this repository (or download the files `max_plus.py` and `int_max_plus.pyx`).
Go to the directory where you cloned or downloaded the file. Start Sage and run

    sage: %runfile max_plus.py

## Symbolic and integer max plus matrices

Once you are able to load the file `max_plus.py`, you can create symbolic
matrices with the following functions

    symbolic_max_plus_matrices_band(d, n, diag, surdiag, ch, sym)
	symbolic_max_plus_matrices_upper(d, n, diag, surdiag, ch, sym)
    symbolic_max_plus_matrices(d, n, ch, sym)

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

- `sym` - (optional, default to `True`) specifies the implementation: if `sym`
  is `True` the matrices do not store all the coefficients (only two in the case
  of full matrices and `d` in the case of triangular ones). The multiplication
  in these case should be a little bit faster (by a factor `d^2` for dense
  matrix but much less for triangular ones). This should be set to `False` only
  to debug the code or to play with `diag=c` or `surdiag=c` that are not
  compatible with this option.

For example you can do

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
method `eval`. For that you first need to compile `int_max_plus.pyx`

    sage: %runfile int_max_plus.pyx
    Compiling ./int_max_plus.pyx...
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

There are currently few functions to check some relations combinatorially in `B^{sv}_d`. The identities must be written with the letters 'x' and 'y'.

- `def occurrences(w, u)`: compute the position of the occurrences of `u` in `w`

- `is_sv_identity(u1, u2, d)`: check whether `(u1,u2)` is an identity in `B^{sv}_d`.

- `is_sv_identity_parallel(u1, u2, d, prefix_length)`: the same as above but
  with parallelization. The argument `prefix_length` is used to chunk the
  subwords to test into different jobs. If it is set to `k` then there will be
  `2^k` jobs which correspond to the `2^k` possible prefixes.

Here are some examples

    sage: p = 'xyxxyy'
    sage: s = 'xxyyxy'
    sage: is_sv_identity(p+'x'+s, p+'y'+s, 3)
    True

    sage: p = 'xyyxxxyyyyxxxx'
    sage: s = 'yyyyxxxxyyyxxy'
    sage: %time is_sv_identity(p+'x'+s, p+'y'+s, 5)
	CPU times: user 1.79 s, sys: 12 ms, total: 1.8 s
	Wall time: 1.76 s
	True
    sage: %time is_sv_identity_parallel(p+'x'+s, p+'y'+s, 5, 3)
	CPU times: user 8 ms, sys: 20 ms, total: 28 ms
	Wall time: 1.03 s
	True

For d=6 the computation takes around 40secs (note the `verbose` option to get
information about the ongoing computation):

	sage: p = 'xyyxxxyyyyxxxxxyyyyy'
	sage: s = 'xxxxxyyyyyxxxxyyyxxy'
	sage: is_sv_identity_parallel(p+'x'+s, p+'y'+s, 6, 3, verbose=True)
	PoolWorker-25: new job at 22:45:26
	...
	PoolWorker-26: new job at 22:45:26
	...
	PoolWorker-27: new job at 22:45:26
	...
	PoolWorker-28: new job at 22:45:26
	...
    PoolWorker-25: job done in 19.6683559418 seconds
    PoolWorker-25: new job at 22:45:45
    PoolWorker-28: job done in 22.5837759972 seconds
    PoolWorker-28: new job at 22:45:48
    ...
    PoolWorker-26: job done in 22.8013679981 seconds
    PoolWorker-26: new job at 22:45:49
    ...
    PoolWorker-27: job done in 26.1923320293 seconds
    PoolWorker-27: new job at 22:45:52
    ...
    PoolWorker-25: job done in 21.8572430611 seconds
    PoolWorker-27: job done in 17.1816589832 seconds
    PoolWorker-26: job done in 20.6833930016 seconds
    PoolWorker-28: job done in 22.0388650894 seconds
    computation with 4 cpus performed in 44.6829161644 seconds
    True

## Contact

Any comment or remark is welcome at vincentDOTdelecroixATlabriDOTfr
