from __future__ import division, print_function
import molgif.utils as utils
from molgif.molecule import Molecule
import os
import sys
import re
import subprocess
import platform
import numpy as np
import ase.io
from ase.data import chemical_symbols
import matplotlib
import matplotlib.cm as cm
import matplotlib.animation as anim
import matplotlib.pyplot as plt

# ensure that ImageMagick writer is found (and not Windows' convert.exe)
if platform.system().lower().startswith('windows'):
    if 'rcParams' in anim.__dir__():
        rc = anim.rcParams
    else:
        rc = matplotlib.rcParams

    if rc['animation.convert_path'].endswith(('convert', 'convert.exe')):
        rc['animation.convert_path'] = 'magick'


def rot_gif(atoms, save_path=None, img=False, vis=False, smart_rotate=False,
            colors=None, loop_time=6, fps=20, scale=0.7, draw_bonds=True,
            hide=None, alphas=None, custom_rotate=None, rot_axis='y',
            anchor=None, max_px=600, square=False, draw_legend=False,
            leg_order='size', legend_max_ms=20, optimize=False,
            transparent=False, overwrite=False, use_charges=False,
            draw_colorbar=False, cb_min=None, cb_max=None, cmap=cm.bwr_r,
            center_data=False, labels=None, label_size=None,
            bond_color='white', bond_edgecolor='black', save_frames=False):
    """
    Creates a rotating animation .gif of a molecule from
    an ase.Atoms object or a geometry file
    - The same gifs can be built dynamically through a molgif.Molecule object

    Args:
    - atoms (ase.Atoms || str): atoms to be animated - can be atoms object or
                                path to geometry file

    KArgs:
    - save_path (str): path to save gif - if None, chem formula and gif info
                       are used
                       (Default: '' -> "<chem formula>-GIF-INFO.gif")
    - img (bool): if True, saves a png image instead of a gif
                  (Default: False)
    - vis (bool): if True, display molecule as a matplotlib figure window
                  - does not save any new files
                  (Default: False)
    - smart_rotate (bool): if True, PCA is applied to coords to orient atoms
                           such that max variance is in x-axis
                           (Default: False)
    - colors (str | iterable | dict): specify atom colors with str, dict,
                                      or values which will use the cmap
                                      - 'blue': all atoms blue
                                      - ['blue', 'white', ...]: label each atom
                                      - [0, 1.1, -2.3...]: cmap used to color
                                      - {'Au': 'purple'}: use dict to color by
                                        by atom type - types not given use jmol
                                      (Default: None -> jmol colors used)
    - loop_time (int): number of seconds for atoms to complete one rotation
                       (Default: 6)
    - fps (int): frames per second of animation
                 (Default: 20)
    - scale (float): scales size of atoms: scale * ase.data.covalent_radii
                     (Default: 0.7)
    - draw_bonds (bool): if True, bonds are drawn
                        (Default: True)
    - hide (list): select atom types and/or atom indices to hide
                   - hide occurs before smart_rotate is applied
                   (Default: None)
    - alphas (int iterable dict): set transparency of atom type or atom indices
                                  (Default: None)
    - custom_rotate (list): list of ordered rotation commands to apply
                            - [angle, (x|y|z), angle, (x|y|z), ...]
                            - Always occurs AFTER smart_rotate
    - rot_axis (str): specify axis to rotate about
                      - x (top-to-bot), y (left-to-right), or
                        z (counterclockwise)
                      - can be: 'y' | 'x' | 'z',
                      - can also add a '-' to change direction!
                      e.g. '-z' causes clockwise rotation
                      (Default: 'y')
    - anchor (int | str): if given, atoms[anchor] will be set to the origin
                          so all other atoms rotate around it while it remains
                          stationary
                          - int: index of atom to anchor
                          - "center": chooses closest atom to COP
                          - chem symbol: chooses first instance
                                         (based on index order)
                                         of atom type = <chem symbol>
                          (Default: None)
    - max_px (int): sets pixel count for longest dimension in gif image
                    (Default: 600)
    - square (bool): if True, gif will be saved with square dimensions
                     (Default: False)
    - draw_legend (bool): if True, a legend specifying atom types is added
                         (Default: False)
    - leg_order (list | str): if given, use it to order the legend
                              - can also give str of single atom type
                              - 'size': largest to smallest
                              - 'size_r': smallest to largest
                              - 'alpha': (alphabetical order)
                              - 'number': atomic number
                              (Default: 'size')
    - legend_max_ms (int): scales legend such that largest atom type
                           is represented with markersize=<legend_max_ms>
                           (Default: 20pts)
    - optimize (bool): if True, creates optimized visual (~1/2 file size)
                       - experimental; still needs additional testing
                       - NOTE: gif takes much longer to make
                       (Default: False)
    - transparent (bool): if True, images are saved as transparent images
                          - Does not affect gifs
                          (Default: False)
    - overwrite (bool): if False, '-<integer>' is added to end of name
                        to avoid overwriting
                        (Default: False)
    - use_charges (bool): if True, colored by initial_charges in atoms obj
                          (Default: False)
    - draw_colorbar (bool): Draws a colorbar if values are given in <colors>
                            (Default: False)
    - cb_min (float): min value limit for colorbar
                      (Default: calc based on values)
    - cb_max (float): max value limit for colorbar
                      (Default: calc based on values)
    - cmap (ColorMap): cmap to be used if colors is specified
                       (Default: matplotlib.cm.bwr_r)
    - center_data (bool): if True, colors are centered about middle of cmap
                          - ensures (-) and (+) values are different color
                          - ex) for RdBu cmap, 0 = 'white'
                          (Default: False)
    - labels (str | iterable): type or list of labels to add to atoms
                               - 'symbols': uses chemical symbol
                               - 'symbols-noh': same as symbols, but no H
                               - 'values': uses values from colors KArg
                               - 'charges': uses initial_charges from atoms
                               - [lab1, lab2, ...]
                               (Default: None)
    - label_size (int): size of labels
    - bond_color (str): specify color of bonds
                        (Default: white)
    - bond_edgecolor (str): specify edgecolor (border) of bonds
                            (Default: black)
    - save_frames (bool): if True, folder is made - frames saved as pngs
                          - NOTE: won't make a gif
                          (Default: False)
    """
    # use utils function to handle str to geometry path
    atoms = utils.path2atoms(atoms)
    atoms = atoms.copy()

    # hide atoms if given
    if hide is not None:
        toremove = []
        for h in hide:
            if isinstance(h, int) or h.isdigit():
                toremove.append(int(h))
            elif str(h).title() in chemical_symbols:
                toremove += [i for i
                             in np.where(atoms.symbols == h.title())[0]]
            # can hide R group (C and H atoms)
            elif str(h).upper() == 'R':
                toremove += [i for i
                             in np.where(np.isin(atoms.symbols,
                                                 ['C', 'H']))[0]]
        for r in sorted(toremove, reverse=True):
            atoms.pop(r)

    # color atoms based on charge (also centers color scheme around 0)
    if use_charges:
        colors = atoms.get_initial_charges().copy()
        center_data = True
        # draw_colorbar = True

    # build figure object
    molecule = Molecule(atoms, scale=scale, colors=colors, alphas=alphas,
                        bond_color=bond_color, bond_edgecolor=bond_edgecolor,
                        labels=labels, cb_min=cb_min, cb_max=cb_max,
                        center_data=center_data, cmap=cmap, square=square,
                        rot_axis=rot_axis, draw=[])

    # anchor specific atom to origin (all other atoms will rotate around it)
    if anchor is not None:
        molecule.anchor(anchor)

    if smart_rotate:
        molecule.smart_rotate()

    if custom_rotate is not None:
        molecule.custom_rotate(custom_rotate)

    # define label size if given
    if label_size is not None:
        molecule.label_size = label_size

    # draw objects
    molecule.draw_atoms()
    if draw_bonds:
        molecule.draw_bonds()
    if labels is not None:
        molecule.draw_labels()
    if draw_legend:
        molecule.draw_legend(leg_order=leg_order, max_ms=legend_max_ms)
    if draw_colorbar:
        molecule.draw_colorbar()

    # turn on img flag if save_path ends with png or svg extension
    if save_path is not None and save_path[-3:].lower() in ['png', 'svg']:
        img = True

    # do not create gif if vis or img flag is true
    if vis or img:
        # only open matplotlib interactive plot if vis
        if vis:
            molecule.show()
        # if img, only save an img
        if img:
            # if no path given, create one with molecular formula (name)
            if save_path is None:
                save_path = molecule.name + '.png'

            # else remove gif extension and replace with png
            elif save_path.endswith('.gif'):
                save_path = save_path.strip('.gif') + '.png'
            molecule.save(save_path, overwrite, max_px, transparent, optimize)
        return

    # create and save rotation gif
    molecule.save_rot_gif(path=save_path, fps=fps, loop_time=loop_time,
                          max_px=max_px, save_frames=save_frames,
                          overwrite=overwrite, optimize_gif=optimize)
