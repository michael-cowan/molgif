from __future__ import division, print_function
import molgif.utils as utils
import os
import re
import platform
import ase
from ase.data import covalent_radii, chemical_symbols
from ase.data.colors import jmol_colors
import numpy as np
import matplotlib
from matplotlib.lines import Line2D
import matplotlib.cm as cm
import matplotlib.animation as anim
import matplotlib.pyplot as plt

# ensure that ImageMagick writer is found (and not Windows' convert.exe)
if platform.system().lower().startswith('windows'):
    if anim.rcParams['animation.convert_path'].endswith(('convert',
                                                        'convert.exe')):
        anim.rcParams['animation.convert_path'] = 'magick'


def rot_gif(atoms, save_path, loop_time=8, fps=20, scale=0.7, add_bonds=True,
            auto_rotate=False, recenter=True, anchor=None, rot_axis='y',
            add_legend=False, colors=None, center_data=True, colorbar=False,
            cmap=cm.bwr_r, use_charges=False, max_px=600, bond_width=0.25,
            direction='ccw', leg_order=None, legend_max_ms=20, labels=None):
    """
    Creates a rotating animation .gif of ase.Atoms object

    Args:
        - atoms (ase.Atoms): atoms to be animated
        - save_path (str): path to save gif

    KArgs:
        - loop_time (int): number of seconds for atoms to complete one rotation
                           (Default: 8)
        - fps (int): frames per second in animation
                     (Default: 20)
        - scale (float): scales size of atoms: scale * ase.data.covalent_radii
                         (Default: 0.9)
        - add_bonds (bool): if True, bonds are drawn
                            (Default: True)
        - auto_rotate (bool): if True, PCA is applied to coords to orient atoms
                              such that max variance is in x-axis
                              (Default: False)
        - recenter (bool): if True, atoms are centered to origin
                           based on avg. coord
                           (Default: True)
        - anchor (int): if given, atoms[anchor] will be set to the origin
                        so all other atoms rotate around it while it remains
                        stationary
                        (Default: None)
        - rot_axis (str): specify axis to rotate about
                          - x (left-to-right), y (bot-to-top)
                          - can be: 'y' | 'x' | 'z',
                          - can also be '-x' to invert rotation
                            (same as changing direction)
                          (Default: 'y')
        - add_legend (bool): if True, a legend specifying atom types is added
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
        - colorbar (bool): if True and colors given, a colorbar is added to gif
                           (Default: False)
        - cmap (ColorMap): cmap to be used if colors is specified
                           (Default: matplotlib.cm.bwr_r)
        - use_charges (bool): if True, colored by initial_charges in atoms obj
                              (Default: False)
        - max_px (int): sets pixel count for longest side
                        (Default: 600)
        - bond_width (float): width of bonds to be drawn (in Angstroms)
                              (Default: 0.2 Angstrom)
        - direction (str): direction for molecule to rotate
                           - rot_axis='y': (looking down from the top)
                           - rot_axis='x': (looking from the right)
                           - rot_axis='z': (looking into screen)
                           OPTIONS:
                           - 'ccw': counterclockwise [left-to-right]
                           - 'cw': clockwise [right-to-left]
                           (Default: 'ccw')
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
    """
    # if save_path is not a gif, give it a gif extension
    if not save_path.lower().endswith('.gif'):
        save_path += '.gif'

    # total number of frames required
    frames = fps * loop_time

    # rotation angles for atoms object
    rot = 360 / frames

    # rot_axis must be x, -x, y, -y, z, or -z
    rot_axis = rot_axis.lower()
    if not re.match('-?[xyz]', rot_axis):
        raise ValueError('Invalid rot_axis given')

    # negate rotation angle if clockwise is specified
    if direction == 'cw':
        rot *= -1
    elif direction != 'ccw':
        print('Incorrect rotation specified - using counterclockwise (ccw)')

    # number of digits in max frames
    ndig = len(str(frames + 1))
    dig_str = '%0{}i'.format(ndig)

    # print rotation gif info
    print('          Title: %s' % (os.path.basename(save_path)))
    print('      Loop time: %i s' % loop_time)
    print('            FPS: %i' % fps)
    print('   Total frames: %i' % frames)
    print(' Building frame: ' + dig_str % 1, end='\r')

    # color atoms based on charge
    if use_charges:
        colors = atoms.get_initial_charges().copy()
        colorbar = True
        add_legend = False
        center_data = True

    # align max variance to x, y, z using PCA
    if auto_rotate:
        # get coordinates of Atoms
        coords = atoms.positions.copy()

        # transform coordinates
        new_coords = utils.pca(coords)

        # set coordinates of Atoms to new transformed coords
        atoms.positions = new_coords

    # center atoms (autorotate always centers atoms)
    elif recenter:
        atoms.positions -= atoms.positions.mean(0)

    # anchor specific atom to origin (all other atoms will rotate around it)
    if anchor:
        if isinstance(anchor, int) and 0 <= anchor <= len(atoms) - 1:
            atoms.positions -= atoms[anchor].position
        else:
            print('Invalid anchor argument was given and will be ignored')

    # calculate figure size
    # calculate axis limits (include offset as buffer)
    fig_size, xlim, ylim = utils.get_fig_bounds(atoms, rot_axis=rot_axis)

    # don't allow colorbar unless values given
    block_colorbar = True

    # don't allow legend if atoms are not colored by type
    block_legend = True

    if colors is None:
        colors = [jmol_colors[i.number] for i in atoms]
        block_legend = False
    elif isinstance(colors, str):
        colors = [colors] * len(atoms)
    elif isinstance(colors, dict):
        colors_dict = colors.copy()

        not_found = []
        symbols = set(atoms.get_chemical_symbols())
        not_found = [c for c in colors_dict if c not in symbols]
        if not_found:
            print('%s do not match atom types.' % (', '.join(not_found)))

        # use combination of color_dict and jmol_colors
        colors = [colors_dict.get(i.symbol, jmol_colors[i.number])
                  for i in atoms]
        block_legend = False

    elif type(colors) in [list, np.ndarray]:
        # if values (ex charge), create Red Blue colormap
        try:
            float(colors[0])

            # if center_data, ensure mid color is 0
            if center_data:
                maxval = max(abs(colors))
                minval = -maxval
            else:
                maxval = max(temp)
                minval = min(temp)

            norm = matplotlib.colors.Normalize(vmin=minval, vmax=maxval)

            # create color map
            colors = [cmap(norm(t)) for t in colors]
            block_colorbar = False
        except:
            assert isinstance(colors[0], str)

    # all atoms must be accounted for in colors list
    assert len(colors) == len(atoms)

    # initialize plt figure and axis
    # add extra subplot if a colorbar or legend is needed
    if (colorbar and not block_colorbar) and not add_legend:
        fig, (ax, extra_ax) = plt.subplots(1, 2,
                                           gridspec_kw={'width_ratios': [30,
                                                                         1]},
                                           figsize=fig_size)
        fig.subplots_adjust(wspace=0, hspace=0)

    # make single subplot
    else:
        fig, ax = plt.subplots(figsize=fig_size)
        if colorbar:
            print('No data given for colorbar!')

    # set axis limits (include offset as buffer)
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.axis('off')

    # set aspect ratio to 1
    ax.set_aspect(1.)

    # see if labels passed in
    if labels is not None:
        if isinstance(labels, str):
            if labels.lower() in ['symbol', 'symbols']:
                labels = atoms.get_chemical_symbols()
            elif labels.lower() in ['color', 'colors']:
                labels = colors.copy()
            elif labels.lower() in ['charge', 'charges']:
                labels = atoms.get_initial_charges()
            elif labels.lower() == 'none':
                labels = None
            else:
                print('"%s" not supported for labels' % labels)
        else:
            try:
                assert len(labels) == len(atoms)
            except:
                print('invalid input of labels.')
                labels = None

    # create atom patches and add to ax
    patches = []
    annotations = []
    for i, a in enumerate(atoms):
        circ = plt.Circle((a.x, a.y),
                          radius=covalent_radii[a.number] * scale,
                          facecolor=colors[i],
                          edgecolor='k',
                          zorder=a.z)
        patches.append(circ)
        ax.add_artist(circ)

        # add element labels (excluding H)
        if labels is not None and a.symbol != 'H':
            ann = ax.annotate(labels[i], (a.x, a.y),
                              zorder=a.z + 0.001,
                              ha='center',
                              va='center',
                              fontsize=7)
        else:
            ann = None
        annotations.append(ann)

    # draw initial bonds
    if add_bonds:
        # calculate bond width info relative to axis units
        bond_width_scaled = utils.angstrom_to_axunits(bond_width, ax)
        bond_fill_scaled = utils.angstrom_to_axunits(bond_width * 0.7, ax)
        bond_info = (bond_width_scaled, bond_fill_scaled)

        radii = np.array([covalent_radii[i.number] for i in atoms])
        bonds = utils.get_bonds(atoms, radii)
        utils.draw_bonds(atoms, ax, radii, bond_info, bonds=bonds)

    # add legend of atom types
    if add_legend:
        if block_legend:
            print('Cannot add legend unless atoms are colored by type.')
        else:
            if colorbar:
                print('Cannot have colorbar and legend. '
                      'Only adding legend to gif')

            # create an ordered, unique list of atom types
            all_symbols = atoms.get_chemical_symbols()
            symbols = sorted(set(all_symbols))

            # use custom legend order if given
            if isinstance(leg_order, str):
                # order legend by size
                if leg_order in ['size', 'size_r']:
                    reverse = True if leg_order == 'size' else False
                    leg_order = sorted(
                        symbols,
                        key=lambda z: list(covalent_radii)[
                            chemical_symbols.index(z)],
                        reverse=reverse)
                else:
                    leg_order = [leg_order]
            if isinstance(leg_order, list) or isinstance(leg_order, np.ndarry):
                # any types not in leg_order will be appended in
                # alphabetical order
                leg_order = list(leg_order) + [s for s in symbols
                                               if s not in leg_order]

                symbols = sorted(symbols, key=lambda z: leg_order.index(z))

            # get an atom object of each type
            a_objs = [atoms[all_symbols.index(s)] for s in symbols]

            # calculate sizes of each atom
            sizes = np.array([utils.angstrom_to_axunits(
                              covalent_radii[a.number] * scale, ax) * 2
                              for a in a_objs])

            # normalize sizes such that largest atom
            # has size of <legend_max_ms>
            sizes = sizes * (legend_max_ms / sizes.max())

            # create legend objects
            leg = [Line2D([0], [0], marker='o', ls='',
                          markerfacecolor=colors[a.index],
                          markeredgecolor='k',
                          markersize=ms)
                   for a, ms in zip(a_objs, sizes)]

            # create legend
            ax.legend(
                leg,
                symbols,
                frameon=False,
                prop=dict(size=11, weight='bold'),
                handletextpad=np.sqrt(utils.angstrom_to_axunits(0.01, ax)),
                borderpad=0,
                borderaxespad=0,
                columnspacing=0,
                markerscale=1,
                labelspacing=np.sqrt(utils.angstrom_to_axunits(0.08, ax)),
                loc='center left',
                framealpha=0,
                # up to ten atom types per column
                ncol=(len(leg) // 10) + 1,
                bbox_to_anchor=(1, 0.5))

    # add colorbar
    elif colorbar and not block_colorbar:
        cb = matplotlib.colorbar.ColorbarBase(extra_ax, cmap=cmap,
                                              norm=norm)
        # set font size
        cb.ax.tick_params(labelsize=11)

    # call tight_layout
    fig.tight_layout()

    def next_step(i):
        # print out progress
        print(' ' * 50, end='\r')
        if (i + 1) == frames:
            print(' Wrapping things up..', end='\r')
        else:
            print(' Building frame: ' + dig_str % (i + 2), end='\r')

        # rotate atoms
        atoms.rotate(rot, v=rot_axis)

        # move atoms (set_center) and change zorder (based on z coord)
        for i, a in enumerate(atoms):
            patches[i].center = (a.x, a.y)
            patches[i].zorder = a.z

            # translates text
            if annotations[i]:
                annotations[i].x = a.x
                annotations[i].zorder = a.z + 0.001

        # redraws bonds
        if add_bonds:
            ax.lines = []
            utils.draw_bonds(atoms, ax, radii, bond_info, bonds=bonds)
            fig.canvas.draw()
        # hi Mike!!

    animation = anim.FuncAnimation(fig, next_step, frames=frames)

    # initialize imagemagick writer
    if anim.writers.is_available('imagemagick'):
        writer = anim.ImageMagickWriter(fps=fps)
    else:
        # ImageMagick must be used
        raise ImportError("ImageMagick must be installed to create GIF")

    # save gif
    animation.save(save_path, writer=writer, dpi=max_px / 5)
    plt.close()
    print(' ' * 50, end='\r')
    print('saved rotation %s' % save_path[-3:])
