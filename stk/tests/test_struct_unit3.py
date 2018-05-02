import numpy as np

from ..utilities import normalize_vector
from ..molecular import StructUnit3


mol = StructUnit3.smiles_init('NCC(N)CN', 'amine')


def test_set_orientation2():
    mol.set_orientation2([1, 2, 3])
    assert np.allclose(mol.bonder_plane_normal(),
                       normalize_vector([1, 2, 3]), atol=1e-4)
