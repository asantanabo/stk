import os

from ..molecular import StructUnit2, Polymer, Linear

if not os.path.exists('linear_topology_tests'):
    os.mkdir('linear_topology_tests')


def test_assembly():
    bb1 = StructUnit2.smiles_init()
    bb2 = StructUnit2.smiles_init()
    bb3 = StructUnit2.smiles_init()
    bb4 = StructUnit2.smiles_init()

    p1 = Polymer([bb1, bb2], Linear())
    p2 = Polymer([bb1, bb2], Linear())
    p3 = Polymer([bb3, bb4], Linear())

    path = os.path.join('linear_topology_tests', 'p1.mol')
    p1.write(path)
    p2.write(path.replace('1', '2'))
    p3.write(path.replace('1', '3'))
