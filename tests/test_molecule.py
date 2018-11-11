from os.path import join
import itertools as it
import numpy as np
from scipy.spatial.distance import euclidean
import stk


def test_all_atom_coords(amine2):
    """
    Test `all_atom_coords`.

    """

    natoms = amine2.mol.GetNumAtoms()
    for i, (atom, coord) in enumerate(amine2.all_atom_coords(), 1):
        assert atom < natoms
        assert type(atom) == int
        assert type(coord) == np.ndarray
        assert len(coord) == 3
    assert natoms == i


def test_atom_coords(amine2):
    """
    Tests `atom_coords`.

    """

    conf = amine2.mol.GetConformer()
    for atom in amine2.mol.GetAtoms():
        atom_id = atom.GetIdx()
        coords = amine2.atom_coords(atom_id)
        conf_coords = conf.GetAtomPosition(atom_id)
        assert np.allclose(coords, conf_coords, atol=1e-8)


def test_atom_distance(amine2):
    """
    Test `atom_distance`.

    """

    # Go through all combinations of atoms in the molecule. Calculate
    # the distance and compare it the distance calculated by the
    # method.
    conf = amine2.mol.GetConformer()
    for atom1, atom2 in it.combinations(amine2.mol.GetAtoms(), 2):
        atom1_id = atom1.GetIdx()
        atom2_id = atom2.GetIdx()
        assert (amine2.atom_distance(atom1_id, atom2_id) ==
                euclidean(conf.GetAtomPosition(atom1_id),
                          conf.GetAtomPosition(atom2_id)))


def test_atom_symbol(amine2):
    """
    Tests the `atom_symbol` method.

    """

    for atom in amine2.mol.GetAtoms():
        atom_id = atom.GetIdx()
        atom_sym = stk.periodic_table[atom.GetAtomicNum()]
        assert atom_sym == amine2.atom_symbol(atom_id)


def test_cavity_size(cc3):
    assert abs(cc3.cavity_size()-6.30569486) < 1e-4


def test_center_of_mass(amine2):
    """
    Tests `center_of_mass`.

    """

    # Calculate the center of mass.
    coord_sum = 0
    total_mass = 0
    for atom_id, coord in amine2.all_atom_coords():
        atom_mass = amine2.mol.GetAtomWithIdx(atom_id).GetMass()
        total_mass += atom_mass
        scaled_coord = np.multiply(atom_mass, coord)
        coord_sum = np.add(scaled_coord, coord_sum)

    com = np.divide(coord_sum, total_mass)
    assert np.allclose(amine2.center_of_mass(), com, atol=1e-6)


def test_centroid_functions(tmp_amine2):
    """
    Tests functions related to centroid manipulation of the molecule.

    Functions tested:
        > centroid
        > set_position

    """

    # Save the coordinates of the new centroid.
    new_centroid = tmp_amine2.centroid() + np.array([10, 20, 4])
    tmp_amine2.set_position(new_centroid)
    # Check that the centroid is at the desired position.
    assert np.allclose(new_centroid, tmp_amine2.centroid(), atol=1e-8)


def test_graph(amine2):
    """
    Tests the output of the `graph` method.

    """

    graph = amine2.graph()
    assert len(graph.nodes()) == amine2.mol.GetNumAtoms()
    assert len(graph.edges()) == amine2.mol.GetNumBonds()


def test_max_diameter(tmp_amine2):
    # Make a position matrix which sets all atoms to the origin except
    # 2 and 13. These should be placed a distance of 100 apart.
    pos_mat = [[0 for x in range(3)] for
               y in range(tmp_amine2.mol.GetNumAtoms())]
    pos_mat[1] = [0, -50, 0]
    pos_mat[3] = [0, 50, 0]
    tmp_amine2.set_position_from_matrix(np.matrix(pos_mat).T)

    d, id1, id2 = tmp_amine2.max_diameter()
    # Note that it is not exactly 100 because of the Van der Waals
    # radii of the atoms.
    assert d > 100 and d < 105
    assert id1 == 1
    assert id2 == 3


def test_position_matrix(amine2):
    """
    Test `postion_matrix`.

    """

    # Go through each atom id. For each atom id get the column in the
    # position matrix with that id as its index. Make sure that the
    # data is the same.
    pos_mat1 = amine2.position_matrix()
    conf = amine2.mol.GetConformer()

    for atom in amine2.mol.GetAtoms():
        atom_id = atom.GetIdx()
        conf_coord = np.array(conf.GetAtomPosition(atom_id))
        mat_coord = pos_mat1.T[atom_id]

        assert np.allclose(conf_coord, mat_coord, atol=1e-8)


def test_same(amine2,
              tmp_amine2,
              aldehyde2):
    """
    Tests the `same()` method.

    """

    assert amine2 is not tmp_amine2
    assert amine2.same(tmp_amine2)

    assert amine2 is not aldehyde2
    assert not amine2.same(aldehyde2)


def test_set_position_from_matrix(tmp_amine2):
    # The new position matrix just sets all atomic positions to origin.
    new_pos_mat = np.matrix([[0 for x in range(3)] for y in
                            range(tmp_amine2.mol.GetNumAtoms())])
    tmp_amine2.set_position_from_matrix(new_pos_mat.T)
    for _, atom_coord in tmp_amine2.all_atom_coords():
        assert np.allclose(atom_coord, [0, 0, 0], atol=1e-8)


def test_shift(amine2):

    s = np.array([10, -20, 5])
    mol2 = amine2.shift(s)
    conf = mol2.GetConformer()
    for atom in mol2.GetAtoms():
        atomid = atom.GetIdx()
        pos = conf.GetAtomPosition(atomid)
        should_be = amine2.atom_coords(atomid) + s
        assert np.allclose(should_be, pos, atol=1e-8)


def test_update_from_mae(tmp_amine2):
    filename = join('data', 'molecule.mae')
    tmp_amine2.update_from_mae(filename, 1)
    assert abs(tmp_amine2.max_diameter(0)[0] -
               tmp_amine2.max_diameter(1)[0]) > 1
