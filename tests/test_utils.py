from pytest import approx
from molgif import utils
from ase.build import molecule
import os


# test_angstrom_to_axunits__:


def test_avoid_overwrite__adds_dash1_to_file():
    # use this file as test input
    this_file = os.path.abspath(__file__)

    # should create a new file path with a '-1' appended to fname
    new_file = this_file.replace('.py', '-1.py')
    assert utils.avoid_overwrite(this_file) == new_file


def test_avoid_overwrite_dir__adds_dash1_to_dir():
    # get the directory path of this file
    this_dir = os.path.dirname(os.path.abspath(__file__))

    # should create new dir path with an appended '-1'
    new_dir = this_dir + '-1'
    assert utils.avoid_overwrite_dir(this_dir) == new_dir


def test_get_bonds__finds_all_methane_bonds():
    methane = molecule('CH4')
    bonds = utils.get_bonds(methane)
    assert bonds.tolist() == [[0, 1], [0, 2], [0, 3], [0, 4]]


def test_get_bonds__0scale_return_empty_list():
    methane = molecule('CH4')
    bonds = utils.get_bonds(methane)
    assert utils.get_bonds(methane, [0] * 5) == []


# test_get_fig_bounds__


# test_opt_angle__

# test_path2atoms__



# test_pca__


# test_smart_rotate_atoms__


