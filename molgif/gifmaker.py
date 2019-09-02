from __future__ import division, print_function
# import molgif.utils as utils
# from molgif.molecule_figure import MolFig
import utils
from molecule_figure import MolFig
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


def rot_gif(atoms, save_path='', overwrite=False, loop_time=8, fps=20,
            scale=0.7, draw_bonds=True, smart_rotate=False,
            anchor=None, rot_axis='y', draw_legend=False, colors=None,
            center_data=True, draw_colorbar=False, cb_range=[None, None],
            cmap=cm.bwr_r, use_charges=False, max_px=600,
            leg_order=None, legend_max_ms=20, labels=None, bond_color='white',
            bond_edgecolor='k', square=False, save_frames=False,
            optimize_gif=False, transparent=False):
    """
    Creates a rotating animation .gif of ase.Atoms object

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
                       (Default: 8)
    - fps (int): frames per second in animation
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
                      - x (left-to-right), y (bot-to-top)
                      - can be: 'y' | 'x' | 'z',
                      - can also be '-' to change direction!
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
                          (Default: True)
    - draw_colorbar (bool): if True and colors given, a colorbar is added to gif
                       (Default: False)
    - cb_range (tuple | list): (minval, maxval) will be used as colorbar
                               range if given
                               (Default: None)
    - cmap (ColorMap): cmap to be used if colors is specified
                       (Default: matplotlib.cm.bwr_r)
    - use_charges (bool): if True, colored by initial_charges in atoms obj
                          (Default: False)
    - max_px (int): sets pixel count for longest side
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
                               - 'symbol': uses chemical symbol
                               - 'colors': uses values from colors KArg
                               - 'charge': uses initial_charges from atoms
                               - [lab1, lab2, ...]
                               (Default: None)
    - bond_color (str): specify color of bonds
                        (Default: white)
    - bond_edgecolor (str): specify edgecolor (border) of bonds
                            (Default: black)
    - square (bool): if True, gif will be saved with square dimensions
                     (Default: False)
    - save_frames (bool): if True, folder is made and frames are saved as pngs
                          - NOTE: gif will not be made if True
                          (Default: False)
    - optimize_gif (bool): if True, creates optimized gif (~1/2 file size)
                           - experimental; still needs additional testing
                           - NOTE: gif takes much longer to make
                           (Default: False)
    - transparent (bool): if True, frames are saved as transparent images
                          (Defulat: True)
    """
    # is atoms is str, read in atoms object
    if isinstance(atoms, str):
        try:
            atoms = ase.io.read(atoms)
        except:
            raise ValueError("Invalid geometry path.")

    # if directory or nothing passed in, use chemical formula as name
    if os.path.isdir(save_path) or not save_path:
        name = atoms.get_chemical_formula()
        noname = True
        save_path = os.path.join(save_path, name)
    else:
        name = os.path.basename(save_path)
        noname = False

    # remove 3-letter extensions
    if '.' in save_path and len(save_path.split('.')[-1]) == 3:
        save_path = save_path[:-4]

    # get the directory path
    dirpath = os.path.dirname(save_path)

    # total number of frames required
    frames = int(round(fps * loop_time))

    # number of digits in max frames
    ndig = len(str(frames + 1))
    dig_str = '%0{}i'.format(ndig)

    # rotation angles for atoms object
    rot = 360 / frames

    # rot_axis must be x, -x, y, -y, z, or -z
    rot_axis = rot_axis.lower()
    if not re.match('-?[xyz]', rot_axis):
        raise ValueError('Invalid rot_axis given')

    # color atoms based on charge
    if use_charges:
        if noname:
            save_path = save_path + '-charges'
        colors = atoms.get_initial_charges().copy()
        draw_colorbar = True
        draw_legend = False
        center_data = True

    # build figure object
    molecule = MolFig(atoms, scale=scale, colors=colors,
                      bond_color=bond_color, bond_edgecolor=bond_edgecolor,
                      labels=labels, cb_min=cb_range[0],
                      cb_max=cb_range[1], center_data=center_data,
                      cmap=cmap, square=square, rot_axis=rot_axis)

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

    # define how to transition to next frame
    def next_step(i):
        # print out progress
        print(' ' * 50, end='\r')
        if (i + 1) == frames:
            print(' Wrapping things up..', end='\r')
        else:
            print(' Building frame: ' + dig_str % (i + 2), end='\r')

        # rotate atoms
        molecule.rotate(rot, rot_axis)
        # hi Mike!!

    # print rotation gif info
    print('          Title: %s' % (os.path.basename(save_path)))
    print('      Loop time: %.2f s' % loop_time)
    print('            FPS: %i' % fps)
    print('   Total frames: %i' % frames)
    print(' Building frame: ' + dig_str % 1, end='\r')

    # define a gif path
    gif_path = save_path + '.gif'

    # make sure file is not overwritten if <overwrite> = False
    if not overwrite:
        gif_path = utils.avoid_overwrite(gif_path)

    # determine runtype
    runtype = ['gif', 'frames'][save_frames]

    if optimize_gif and transparent:
        print('Can only optimize gif if images are not transparent.')
        transparent = False

    # only save png frames if <save_frames>
    if save_frames or optimize_gif:
        frame_path = os.path.join(dirpath, name + '_frames')
        # create 'frames' directory to save png frames
        if not os.path.isdir(frame_path):
            os.mkdir(frame_path)
        # TODO: NEW KARG???
        for i in range(frames):
            molecule.save(os.path.join(frame_path,
                                       name + '_%03i.png' % (i + 1)),
                          transparent=transparent)
            next_step(i)
    else:
        # build frames
        animation = anim.FuncAnimation(molecule.fig,
                                       next_step,
                                       frames=frames)

        # initialize imagemagick writer
        if anim.writers.is_available('imagemagick'):
            writer = anim.ImageMagickWriter(fps=fps)
        else:
            # ImageMagick must be used
            raise ImportError("ImageMagick must be installed to create GIF")

        # save gif
        animation.save(gif_path, writer=writer, dpi=max_px / 5)

    if optimize_gif:
        inp = os.path.join(frame_path, '*.png')
        subprocess.call(['magick',
                         'convert',
                         '-delay',
                         '%.1f' % (100 / fps),
                         inp,
                         gif_path])

        # delete frames
        for f in os.listdir(frame_path):
            os.remove(os.path.join(frame_path, f))
        os.removedirs(frame_path)

    # close figure and report completion
    plt.close()
    print(' ' * 50, end='\r')
    print('saved rotation %s' % runtype)

if __name__ == '__main__':
    a = ase.build.molecule('C2H6')
    rot_gif(a, save_path='C:\\users\\mcowa\\desktop\\test.gif')
