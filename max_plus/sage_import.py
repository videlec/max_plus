from sage.misc.misc import SAGE_TMP
from sage.misc.misc_c import prod
from sage.misc.temporary_file import tmp_filename

from sage.geometry.polyhedron.parent import Polyhedra

#try:
#    from ppl import (Variable, C_Polyhedron, point, Generator_System,
#         Linear_Expression, Constraint_System, MIP_Problem)
#except ImportError:
from sage.libs.ppl import (Variable, C_Polyhedron, point, Generator_System,
         Linear_Expression, Constraint_System, MIP_Problem)

from sage.modules.free_module import FreeModule

from sage.numerical.mip import MixedIntegerLinearProgram, MIPSolverException

from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ

from sage.structure.sage_object import SageObject
from sage.structure.element import generic_power

try:
    from sage.combinat.words.words import FiniteWords
except ImportError:
    from sage.combinat.words.words import Words as FiniteWords
from sage.combinat.words.suffix_trees import SuffixTrie

