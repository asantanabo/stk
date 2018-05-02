import numpy as np

from ..utilities import normalize_vector
from ..molecular import StructUnit2

mol = StructUnit2.smiles_init('NCCCN', 'amine')


def test_set_orientation2():
    mol.set_orientation2([1, 2, 3])
    assert np.allclose(next(mol.bonder_direction_vectors())[-1],
                       normalize_vector([1, 2, 3]), atol=1e-8)
