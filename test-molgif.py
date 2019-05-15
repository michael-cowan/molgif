import molgif
import os
import ase.build

desk = os.path.join(os.path.expanduser('~'), 'Desktop')

for name in ['biphenyl', 'C60', 'C7NH5', 'BDA']:
    mol = ase.build.molecule(name)
    path = os.path.join(desk, '%s_nobonds' % name)

    # no bonds
    molgif.rot_gif(mol, path, auto_rotate=True, add_bonds=False, scale=0.9)
