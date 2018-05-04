from ..molecular import (StructUnit3, StructUnit2, Periodic,
                         Honeycomb, Hexagonal, Square, Kagome,
                         NoLinkerHoneycomb)
import os
from os.path import join
import rdkit.Chem.AllChem as rdkit

test_dir = 'cof_topology_tests'
if not os.path.exists(test_dir):
    os.mkdir(test_dir)


bb1 = StructUnit2.smiles_init('Nc1ccc(N)cc1', 'amine')


def test_honeycomb():
    bb2 = StructUnit3(join('data', 'cof', 'aldehyde3f.mol'))
    cof = Periodic([bb1, bb2], Honeycomb())
    path = join(test_dir, 'honeycomb.sdf')
    cof.write(path)
    island = cof.island([3, 3, 1])
    rdkit.MolToMolFile(island, path.replace('.sdf', '_island.sdf'))


def test_hexagonal():
    bb2 = StructUnit3(join('data', 'cof', 'aldehyde6f.mol'))
    cof = Periodic([bb1, bb2], Hexagonal())
    path = os.path.join(test_dir, 'hexagonal.sdf')
    cof.write(path)
    island = cof.island([3, 3, 1])
    rdkit.MolToMolFile(island, path.replace('.sdf', '_island.sdf'))


def test_square():
    bb2 = StructUnit3(join('data', 'cof', 'aldehyde4f_1.mol'))
    cof = Periodic([bb1, bb2], Square())
    path = os.path.join(test_dir, 'square.sdf')
    cof.write(path)
    island = cof.island([3, 3, 1])
    rdkit.MolToMolFile(island, path.replace('.sdf', '_island.sdf'))


def test_kagome():
    bb2 = StructUnit3(join('data', 'cof', 'aldehyde4f_2.mol'))
    cof = Periodic([bb1, bb2], Kagome(multitopic_aligners=[1, 0, 1]))
    path = os.path.join(test_dir, 'kagome.sdf')
    cof.write(path)
    island = cof.island([3, 3, 1])
    rdkit.MolToMolFile(island, path.replace('.sdf', '_island.sdf'))


def test_boron_cof():
    bb1 = StructUnit2.smiles_init('Oc1cc2cc(O)c(O)nc2cc1O', 'diol')
    bb2 = StructUnit3(join('data', 'cof', 'boronic_acid.sdf'))
    cof = Periodic([bb1, bb2], Square())
    path = join(test_dir, 'boron.sdf')
    cof.write(path)
    island = cof.island([3, 3, 1])
    rdkit.MolToMolFile(island, path.replace('.sdf', '_island.sdf'))


def test_nolinkerhoneycomb():
    bb1 = StructUnit3.smiles_init('O=Cc1cc(C=O)nc(C=O)c1', 'aldehyde')
    bb2 = StructUnit3.smiles_init('NCc1cc(CN)nc(CN)n1', 'amine')
    cof = Periodic([bb1, bb2], NoLinkerHoneycomb())
    path = join(test_dir, 'nolinkerhoneycomb.sdf')
    cof.write(path)
    island = cof.island([3, 3, 1])
    rdkit.MolToMolFile(island, path.replace('.sdf', '_island.sdf'))
