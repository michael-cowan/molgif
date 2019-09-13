from __future__ import division, print_function
import molgif.utils as utils
from molgif.molecule import Molecule
import os
import sys
import re
import subprocess
import platform
import ase.io
import matplotlib.cm as cm
import matplotlib.animation as anim
import matplotlib.pyplot as plt

# ensure that ImageMagick writer is found (and not Windows' convert.exe)
if platform.system().lower().startswith('windows'):
    if anim.rcParams['animation.convert_path'].endswith(('convert',
                                                        'convert.exe')):
        anim.rcParams['animation.convert_path'] = 'magick'


def rot_gif(atoms, save_path=None, overwrite=False, loop_time=6, fps=20,
            scale=0.7, draw_bonds=True, smart_rotate=False,
            anchor=None, rot_axis='y', draw_legend=False, colors=None,
            center_data=False, draw_colorbar=False, cb_min=None, cb_max=None,
            cmap=cm.bwr_r, use_charges=False, max_px=600,
            leg_order=None, legend_max_ms=20, labels=None, bond_color='white',
            bond_edgecolor='k', square=False, save_frames=False,
            optimize_gif=False, transparent=False):
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
    - overwrite (bool): if False, '-<integer>' is added to end of name
                        to avoid overwriting
                        (Default: False)
    - loop_time (int): number of seconds for atoms to complete one rotation
                       (Default: 6)
    - fps (int): frames per second of animation
                 (Default: 20)
    - scale (float): scales size of atoms: scale * ase.data.covalent_radii
                     (Default: 0.9)
    - draw_bonds (bool): if True, bonds are drawn
                        (Default: True)
    - smart_rotate (bool): if True, PCA is applied to coords to orient atoms
                           such that max variance is in x-axis
                           (Default: False)
    - anchor (int): if given, atoms[anchor] will be set to the origin
                    so all other atoms rotate around it while it remains
                    stationary
                    (Default: None)
    - rot_axis (str): specify axis to rotate about
                      - x (left-to-right), y (bot-to-top), or
                        z (counter-clockwise)
                      - can be: 'y' | 'x' | 'z',
                      - can also add a '-' to change direction!
                      e.g. '-z' causes clockwise rotation
                      (Default: 'y')
    - draw_legend (bool): if True, a legend specifying atom types is added
                         (Default: False)
    - colors (str | iterable | dict): specify atom colors with str, dict,
                                      or values which will use the cmap
                                      - 'blue': all atoms blue
                                      - ['blue', 'white', ...]: label each atom
                                      - [0, 1.1, -2.3...]: cmap used to color
                                      - {'Au': 'purple'}: use dict to color by
                                        by atom type - types not given use jmol
                                      (Default: None -> jmol colors used)
    - center_data (bool): if True, colors are centered about middle of cmap
                          - ensures (-) and (+) values are different color
                          - ex) for RdBu cmap, 0 = 'white'
                          (Default: False)
    - draw_colorbar (bool): Draws a colorbar if values are given in <colors>
                            (Default: False)
    - cb_min (float): min value limit for colorbar
                      (Default: calc based on values)
    - cb_max (float): max value limit for colorbar
                      (Default: calc based on values)
    - cmap (ColorMap): cmap to be used if colors is specified
                       (Default: matplotlib.cm.bwr_r)
    - use_charges (bool): if True, colored by initial_charges in atoms obj
                          (Default: False)
    - max_px (int): sets pixel count for longest dimension in gif image
                    (Default: 600)
    - leg_order (list | str): if given, use it to order the legend
                              - can also give str of single atom type
                              - 'size': largest to smallest
                              - 'size_r': smallest to largest
                              (Default: None (alphabetical order))
    - legend_max_ms (int): scales legend such that largest atom type
                           is represented with markersize=<legend_max_ms>
                           (Default: 20pts)
    - labels (str | iterable): type or list of labels to add to atoms
                               - 'symbols': uses chemical symbol
                               - 'symbols-noh': same as symbols, but no H
                               - 'values': uses values from colors KArg
                               - 'charges': uses initial_charges from atoms
                               - [lab1, lab2, ...]
                               (Default: None)
    - bond_color (str): specify color of bonds
                        (Default: white)
    - bond_edgecolor (str): specify edgecolor (border) of bonds
                            (Default: black)
    - square (bool): if True, gif will be saved with square dimensions
                     (Default: False)
    - save_frames (bool): if True, folder is made - frames saved as pngs
                          - NOTE: won't make a gif
                          (Default: False)
    - optimize_gif (bool): if True, creates optimized gif (~1/2 file size)
                           - experimental; still needs additional testing
                           - NOTE: gif takes much longer to make
                           (Default: False)
    - transparent (bool): if True, frames are saved as transparent images
                          (Defulat: True)
    """
    # use utils function to handle str to geometry path
    atoms = utils.path2atoms(atoms)
    atoms = atoms.copy()

    # color atoms based on charge
    if use_charges:
        colors = atoms.get_initial_charges().copy()
        center_data = True
        draw_colorbar = True

    # build figure object
    molecule = Molecule(atoms, scale=scale, colors=colors,
                        bond_color=bond_color, bond_edgecolor=bond_edgecolor,
                        labels=labels, cb_min=cb_min,
                        cb_max=cb_max, center_data=center_data,
                        cmap=cmap, square=square, rot_axis=rot_axis, draw=[])

    # anchor specific atom to origin (all other atoms will rotate around it)
    if anchor is not None:
        molecule.anchor(anchor)

    if smart_rotate:
        molecule.smart_rotate()

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

    # create and save rotation gif
    molecule.save_rot_gif(path=save_path, fps=fps, loop_time=loop_time,
                          max_px=max_px, save_frames=save_frames,
                          overwrite=overwrite, optimize_gif=optimize_gif)
