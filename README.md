# Sage library for max-plus matrix identities

To run the library you need to install Sage (http://sagemath.org). Then clone
this repository (or download the files `max_plus.py` and `int_max_plus.pyx`).
Go to the directory where you cloned or downloaded the file. Start Sage and run

    sage: %runfile max_plus.py

Once that done, you can create symbolic matrices with the following functions

    symbolic_max_plus_matrices_band(d, n, diag, surdiag, i, ch)
	symbolic_max_plus_matrices_upper(d, n, diag, surdiag, i, ch)
    symbolic_max_plus_matrices(d, n, i, ch, sym)

for respectively band, upper triangular and full matrices. The arguments are

- `d` - the dimension

- `n` - the number of matrices

- `diag` - describes whether the diagonal is zero (`'z'`), constant (`'c'`),
  the same in each of the `n` matrices (`'s'`) or variable (`'v'`). This is
  optional and default to `'v'`.

- `surdiag` - describes whether the surdiagonal is constant, zero, same or
  variable.  This is also optional.

- `i` - an optional number to set a variable to zero. If set, the number of variable
  is one less and hence convex hull computations much faster.

- `ch` - (optional) set the convex hull engine (could be one of 'ppl', 'cdd'
  or 'PALP'). It defaults to 'ppl' which seems to be the fastest.

- `sym` - (optional) specifies the implementation: if `sym` is `True` uses matrices that
  store only two coefficients (and deduce the others by symmetry). If `sym` is
  `False` uses matrices that store all coefficients.

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
	[  x0   x3  x4 ]
	[ -oo   x1  x5 ]
	[ -oo  -oo  x2 ]
	sage: print y2
	[  x6   x9  x10 ]
	[ -oo   x7  x11 ]
	[ -oo  -oo   x8 ]
	sage: print z2
	[ x12  x15  x16 ]
	[ -oo  x13  x17 ]
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

Any comment or remark is welcome at vincentDOTdelecroixATlabriDOTfr
