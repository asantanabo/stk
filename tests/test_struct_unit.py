import os
import numpy as np
import rdkit.Chem.AllChem as rdkit
import stk

if not os.path.exists('struct_unit_tests_output'):
    os.mkdir('struct_unit_tests_output')


def test_init(amine2):
    assert len(amine2.functional_group_atoms()) == 2
    assert amine2.func_grp.name == 'amine'


def test_all_bonder_distances(tmp_aldehyde3):
    shape = (3, tmp_aldehyde3.mol.GetNumAtoms())
    tmp_aldehyde3.set_position_from_matrix(np.zeros(shape))
    bonder_distances = tmp_aldehyde3.all_bonder_distances()
    for i, (id1, id2, d) in enumerate(bonder_distances):
        assert type(id1) is int
        assert type(id2) is int
        assert abs(d) < 1e-5
    assert i == 2


def test_bonder_centroids(tmp_aldehyde3):
    shape = (3, tmp_aldehyde3.mol.GetNumAtoms())
    tmp_aldehyde3.set_position_from_matrix(np.zeros(shape))

    for i, centroid in enumerate(tmp_aldehyde3.bonder_centroids()):
        assert len(centroid) == 3
        assert sum(centroid) < 1e-5
    assert i == 2


def test_bonder_centroid(tmp_aldehyde3):
    shape = (3, tmp_aldehyde3.mol.GetNumAtoms())
    tmp_aldehyde3.set_position_from_matrix(np.zeros(shape))
    centroid = tmp_aldehyde3.bonder_centroid()
    assert len(centroid) == 3
    assert sum(centroid) < 1e-5


def test_bonder_direction_vectors(aldehyde3):
    dir_vectors = aldehyde3.bonder_direction_vectors()
    for i, (id1, id2, v) in enumerate(dir_vectors):
        assert type(id1) == int
        assert type(id2) == int
        assert type(v) == np.ndarray
    assert i == 2


def test_bonder_position_matrix(tmp_aldehyde3):
    shape = (3, tmp_aldehyde3.mol.GetNumAtoms())
    tmp_aldehyde3.set_position_from_matrix(np.zeros(shape))
    position_matrix = tmp_aldehyde3.bonder_position_matrix()
    assert np.allclose(position_matrix,
                       np.zeros(position_matrix.shape),
                       atol=1e-8)


def test_centroid_centroid_dir_vector(aldehyde3):
    c1 = aldehyde3.bonder_centroid()
    c2 = aldehyde3.centroid()
    assert np.allclose(stk.normalize_vector(c2-c1),
                       aldehyde3.centroid_centroid_dir_vector(),
                       atol=1e-8)


def test_core(amine2):
    for atom in amine2.core().GetAtoms():
        assert atom.GetAtomicNum() != 1
        assert not atom.HasProp('fg')


def test_functional_group_atoms(amine2):
        func_grp_mol = rdkit.MolFromSmarts(amine2.func_grp.fg_smarts)
        fg_atoms = amine2.mol.GetSubstructMatches(func_grp_mol)
        assert fg_atoms == amine2.functional_group_atoms()


def test_is_core_atom(amine2):
    for atom in amine2.mol.GetAtoms():
        core = (False if atom.HasProp('fg') or atom.GetAtomicNum() == 1
                else True)
        assert core is amine2.is_core_atom(atom.GetIdx())


def test_json_init(amine2):
    path = os.path.join('struct_unit_tests_output', 'mol.json')
    amine2.dump(path)
    mol2 = stk.Molecule.load(path, stk.Molecule.from_dict)

    assert isinstance(amine2.file, str)
    assert mol2.optimized
    assert mol2.bonder_ids == amine2.bonder_ids
    assert mol2.energy.__class__.__name__ == 'Energy'
    assert mol2.func_grp.name == amine2.func_grp.name
    assert amine2 is not mol2
    assert amine2.atom_props == amine2.atom_props


def test_caching():
    # This test is in a try block because pytest runs with
    # molecule cache turned off. If this test fails, the finally
    # clause ensures that the cache remains off so other tests are
    # not affected by unexpected cache problems.
    try:
        stk.OPTIONS['cache'] = True
        mol = stk.StructUnit.smiles_init('NCCCN', 'amine')
        mol2 = stk.StructUnit.smiles_init('NCCCN', 'amine')
        assert mol is mol2

        mol3 = stk.StructUnit.smiles_init('NCCCN', 'aldehyde3')
        assert mol3 is not mol

        stk.OPTIONS['cache'] = False
        mol4 = stk.StructUnit.smiles_init('NCCCN', 'amine')
        stk.OPTIONS['cache'] = True
        assert mol is not mol4

    except Exception:
        raise

    finally:
        stk.OPTIONS['cache'] = False


def test_set_bonder_centroid(tmp_amine2):
    tmp_amine2.set_bonder_centroid([1, 2, 3])
    assert np.allclose(tmp_amine2.bonder_centroid(),
                       [1, 2, 3],
                       atol=1e-8)


def test_untag_atoms(tmp_amine2):
    assert any(a.HasProp('fg') for a in tmp_amine2.mol.GetAtoms())
    tmp_amine2.untag_atoms()
    assert all(not a.HasProp('fg') for a in tmp_amine2.mol.GetAtoms())
