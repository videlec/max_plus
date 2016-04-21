
from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import Cython.Compiler.Options
import os

#Cython.Compiler.Options.old_style_globals = True

from sage.env import SAGE_LOCAL, SAGE_SRC

with open("version.txt") as f:
    version = f.read()

setup(name='max_plus',
      version=version,
      description="Identities in the max plus monoid",
      author='Vincent Delecroix',
      author_email='vincent.delecroix@labri.fr',
      url='http://www.labri.fr/perso/vdelecro/',
      license = "GPL v3",
      platforms = ["any"],
      packages=['max_plus'],
      ext_modules=[
        Extension('max_plus.max_plus_int',
                sources = [os.path.join('max_plus', 'max_plus_int.pyx')],
                include_dirs = [SAGE_SRC],
            ),
        Extension('max_plus.perm_lex_order',
                sources = [os.path.join('max_plus', 'perm_lex_order.pyx')],
            ),
        Extension('max_plus.word',
                sources = [os.path.join('max_plus', 'word.pyx')],
            ),
        ],
        cmdclass = {'build_ext': build_ext}
)
