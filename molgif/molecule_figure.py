from __future__ import division, print_function
# TEMPORARY IMPORTS
import utils
import ase.build
import ase.io
# END
# import molgif.utils as utils
import os
import io
import subprocess
import copy
import ase
from ase.data import covalent_radii, chemical_symbols
from ase.data.colors import jmol_colors
import numpy as np
from PIL import Image
import matplotlib.colors as mcolors
from matplotlib.lines import Line2D
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import matplotlib.animation as anim


class MolFig(object):
    def __init__(self, atoms, scale=0.7, colors=None, name=None,
                 bond_color='white', bond_edgecolor='black',
                 labels=None, cb_min=None, cb_max=None, center_data=False,
                 cmap=cm.bwr_r, square=False, rot_axis=None):

        # string ids based on params of figure
        self.fig_params = []

        # track what has been drawn
        self._drawn = []

        # figure and axis objects
        self.fig = None
        self.ax = None
        # colorbar axis object
        self.cb_ax = None

        # atom, labels, and bond plot objects (lists)
        self.atom_objs = []
        self.label_objs = []
        self.bond_objs = []

        # rotation axis (if needed)
        self.rot_axis = rot_axis

        # padding around atoms
        self.padding = None

        # position transformation matrix
        # default to no transformation
        self.transform = np.eye(3)

        # square figure (px-x = px-y)
        self.square = square

        # set atoms object
        if isinstance(atoms, ase.Atoms):
            self.atoms = atoms
            # center atoms to COP
            self.atoms.positions -= self.atoms.positions.mean(0)
        else:
            raise ValueError("must pass in valid ase.Atoms object")

        # set name
        self.name = name
        # default name is chemical formula
        if self.name is None:
            self.name = self.atoms.get_chemical_formula()

        # track atom objects (track for redrawing atoms and bonds)
        self.pos = self.atoms.positions.copy()
        self.bond_pos = self.atoms.positions.copy()

        # get number of unique atom types
        self.atom_types = len(set(atoms.symbols))

        # border width of atoms (in Angstrom)
        self.atom_borderwidth = 0.05

        # atom size scale (covalent_radii * scale)
        self.scale = scale

        # get atomic radii based on scale
        self.radii = covalent_radii[self.atoms.numbers] * self.scale

        # bond attributes
        self.bond_color = bond_color
        self.bond_edgecolor = bond_edgecolor
        self.recalc_bonds = None
        # bond width (in Angstrom)
        self.bond_width = 0.25
        # list of bonds - only create if draw_bonds is called
        self.bond_ls = []

        # colorbar attributes
        self.cmap = cmap
        self.cb_min = cb_min
        self.cb_max = cb_max
        # store values passed into colors prop
        self._cb_values = None
        # normalize distribution of colors
        self.cb_norm = None
        # center colorbar data
        self.center_data = center_data

        # labels (write text on atoms)
        self._labels = None
        self.__check_labels__(labels)

        # legend attributes (track order)
        self.legend = None
        self.leg_order = None

        # color attributes
        self._colors = None

        # check input colors and see if they can be used with legend/colorbar
        self.__check_colors__(colors)

        # initialize figure
        self.init_fig()

        # get scaled widths based on axis size
        # atom border width
        self.scaled_borderwidth = utils.angstrom_to_axunits(
                                    self.atom_borderwidth,
                                    self.ax)

        # bond width (black line) and bond fill (white line)
        self.scaled_bond_width = utils.angstrom_to_axunits(
                                    self.bond_width,
                                    self.ax)
        # width of bond border = width of atom border
        bondedge = self.bond_width - 2 * self.atom_borderwidth
        self.scaled_bond_fill = utils.angstrom_to_axunits(
                                    bondedge,
                                    self.ax)

    def __len__(self):
        return len(self.atoms)

    def __repr__(self):
        return 'MolFig(%s)' % self.atoms.get_chemical_formula()

    # colors getter
    @property
    def colors(self):
        return self._colors

    @colors.setter
    def colors(self, value):
        self.__check_colors__(value)

    # labels getter
    @property
    def labels(self):
        return self._labels

    @labels.setter
    def labels(self, value):
        self.__check_labels__(value)

    # tracks items that have been drawn
    @property
    def drawn(self):
        return self._drawn

    # drawn is read only
    @drawn.setter
    def drawn(self, value):
        print('Not allowed to edit drawn list.')

    def __check_colors__(self, value):
        # track to see if colors have changed
        old_colors = copy.deepcopy(self._colors)

        # params to block creation of colorbar and legend
        # - both depend on color assignments given
        self.block_colorbar = True
        self.block_legend = False

        # track if colors were successfully converted
        failed = True

        # resort to default colors
        if value is None:
            failed = False
        # all atoms are same color
        elif isinstance(value, str):
            try:
                self._colors = [mcolors.to_rgba(value)] * len(self)
                if len(set(self.atoms.symbols)) != 1:
                    self.block_legend = True
                failed = False
            except ValueError:
                pass
        elif isinstance(value, dict):
            colors_dict = {k.title(): value[k] for k in value.copy()}

            # find atom types that are not in atoms object
            not_found = [c for c in colors_dict
                         if c.title() not in self.atoms.symbols]
            if not_found:
                print('%s do not match atom types.' % (', '.join(not_found)))

            # use combination of color_dict and jmol_colors
            try:
                self._colors = [mcolors.to_rgba(colors_dict.get(i.symbol,
                                list(jmol_colors[i.number])))
                                for i in self.atoms]
                failed = False
            except ValueError:
                pass

        elif type(value) in [list, np.ndarray]:
            # if values (ex charge), create Red Blue colormap
            try:
                float(value[0])

                # if center_data, ensure mid color is 0
                if self.center_data:
                    maxval = max(abs(value))
                    minval = -maxval
                else:
                    if self.cb_min is None:
                        minval = min(value)
                    else:
                        minval = self.cb_min

                    if self.cb_max is None:
                        maxval = max(value)
                    else:
                        maxval = self.cb_max

                self.cb_norm = mcolors.Normalize(
                                                vmin=minval,
                                                vmax=maxval)

                # save original values
                self._cb_values = value.copy()

                # create color map
                self._colors = [self.cmap(self.cb_norm(t))
                                for t in value]
                self.block_legend = True
                self.block_colorbar = False
                failed = False
            # else see if list contains strings
            except:
                # if length is the same, try to convert items to valid colors
                if len(value) == len(self):
                    try:
                        # try to convert list to rgba values
                        self._colors = [mcolors.to_rgba(c)
                                        for c in value]
                        # not sure if unique atom type colors - block legend
                        self.block_legend = True
                        failed = False
                    except ValueError:
                        pass

        if failed:
            print("Warning: unable to determine color scheme "
                  "- resorting to default/previous colors")

        # use jmol_colors as default
        if self.colors is None:
            self._colors = [list(jmol_colors[i.number]) for i in self.atoms]

        # all atoms must be accounted for in colors list
        if len(self._colors) != len(self.atoms):
            raise ValueError("length of colors must be equal to # atoms")

        # redraw atoms if colors have changed
        if self._colors != old_colors and 'atoms' in self._drawn:
            self.draw_atoms(force=True)

    def __check_labels__(self, value):
        # track to see if labels change
        old_labels = copy.deepcopy(self._labels)

        # see if labels passed in
        if value is not None:
            if isinstance(value, str):
                value = value.lower()
                # atomic symbols
                if value == 'symbols':
                    self._labels = list(self.atoms.get_chemical_symbols())
                # atomic symbols excluding H
                elif value == 'symbols-noh':
                    self._labels = [a.symbol if a.symbol != 'H' else ''
                                    for a in self.atoms]
                # colorbar values
                elif value == 'values' and not self.block_colorbar:
                    self._labels = list(self._cb_values)
                # atomic charges
                elif value == 'charges':
                    self._labels = list(self.atoms.get_initial_charges())
                else:
                    print('"%s" not supported for labels' % value)
            elif len(value) == len(self.atoms):
                self._labels = list(value)
            else:
                print('Warning: Invalid input of labels.')
                self._labels = None

        # redraw atoms if labels have changed
        if self._labels != old_labels and 'atoms' in self._drawn:
            self.draw_atoms(force=True)

    def init_fig(self, padding=1):
        """
        Create initial figure and axis objects
        - determines figure size and axis limits
        - also makes colorbar axis if needed
        """
        if self.padding is None:
            self.padding = padding

        # calculate figure size
        # calculate axis limits (include offset as buffer)
        self.fig_size, self.xlim, self.ylim = utils.get_fig_bounds(
                                                        self.atoms,
                                                        rot_axis=self.rot_axis,
                                                        square=self.square,
                                                        padding=self.padding)

        # create figure and main axis
        if not isinstance(self.fig, plt.Figure):
            self.fig, self.ax = plt.subplots(figsize=self.fig_size)
        else:
            # set figure size
            self.fig.set_size_inches(self.fig_size, forward=True)

        # set axis limits
        self.ax.set_xlim(self.xlim)
        self.ax.set_ylim(self.ylim)

        # turn off axis ticks and lines
        self.ax.set_xticklabels([])
        self.ax.set_yticklabels([])
        self.ax.axis('off')

        # set aspect ratio to 1
        self.ax.set_aspect(1)

    def clear_figure(self):
        """
        Removes everything from figure
        same as remove('all')
        """
        self.remove('all')

    def draw(self, items=['atoms', 'bonds', 'labels'], force=False):
        """
        Draws multiple object types at once
        - Default: draws atoms, bonds, and labels (if any specified)

        KArgs:
        - items (list): list of object types to draw (as str)
                        OPTIONS:
                        - atoms
                        - bonds
                        - labels
                        - legend
                        - colorbar
        """
        for item in items:
            method = 'draw_' + item
            if method in self.__dir__():
                getattr(self, method)(force=force)
            else:
                print(item, 'can not be drawn.')
        if isinstance(self.fig, plt.Figure):
            self.fig.canvas.draw()
            self.fig.canvas.flush_events()

    def draw_atoms(self, force=False):
        """
        Draw atoms or update atom positions
        """
        # determine if drawing or updating
        update = bool(len(self.atom_objs))

        # don't redraw if positions haven't changed
        if update:
            if (self.pos == self.atoms.positions).all():
                if not force:
                    return
            else:
                self.pos = self.atoms.positions.copy()

        # create atom patches and add to ax
        for i, a in enumerate(self.atoms):
            # update atom objects
            if update:
                self.atom_objs[i].set_facecolor(self._colors[i])
                self.atom_objs[i].center = (a.x, a.y)
                self.atom_objs[i].zorder = a.z
            # create atom object patches
            else:
                circ = plt.Circle(
                            (a.x, a.y),
                            radius=self.radii[i],
                            facecolor=self._colors[i],
                            edgecolor='k',
                            linewidth=self.scaled_borderwidth,
                            zorder=a.z)
                # add patch to atom_objs list and to axis
                self.atom_objs.append(circ)
                self.ax.add_artist(circ)

        # add atoms to list of drawn items
        if 'atoms' not in self._drawn:
            self._drawn.append('atoms')

        self.fig.tight_layout()

    def draw_bonds(self, recalc_bonds=None, force=False):
        """
        Draw bonds or update bond positions
        """
        # determine if drawing or updating bonds
        update = bool(len(self.bond_objs))

        # don't redraw if positions haven't changed
        if update:
            if (self.bond_pos == self.atoms.positions).all():
                if not force:
                    return
            else:
                self.bond_pos = self.atoms.positions.copy()

        # set recalc_bonds attribute
        if isinstance(recalc_bonds, bool):
            self.recalc_bonds = recalc_bonds
        elif not isinstance(self.recalc_bonds, bool):
            self.recalc_bonds = False

        # build bond list if needed
        if not len(self.bond_ls) or self.recalc_bonds:
            self.bond_ls = utils.get_bonds(
                                    self.atoms,
                                    covalent_radii[self.atoms.numbers])

        # draw the bonds
        for i, b in enumerate(self.bond_ls):
            p1 = self.atoms[b[0]].position
            p2 = self.atoms[b[1]].position

            # vector between two atoms
            vec = p2 - p1

            # normalize the vector
            normvec = vec / np.linalg.norm(vec)

            # calc pts for black line of bond (outline)
            xy1 = p1 + self.radii[b[0]] * normvec
            xy2 = p2 - self.radii[b[1]] * normvec

            # add buffer to bond length (max buffer when zdist is minimum)
            # use exponential decay for smooth overlap b/n bonds and atoms
            zdist = abs(xy1[2] - xy2[2])
            buffer = 0.078 * np.exp(-zdist)

            xy1 += buffer * normvec
            xy2 -= buffer * normvec

            # get x - y positions for bonds
            x = [xy1[0], xy2[0]]
            y = [xy1[1], xy2[1]]

            # use avg z for zorder
            zorder = (xy1[2] + xy2[2]) / 2

            if update:
                # update bond outline
                self.bond_objs[2 * i].set_data((x, y))
                self.bond_objs[2 * i].set_zorder(zorder - 0.001)

                # update bond
                self.bond_objs[2 * i + 1].set_data((x, y))
                self.bond_objs[2 * i + 1].set_zorder(zorder)
            else:
                # draw thinner white line over black line
                # to create bordered bonds
                self.ax.plot(
                        x,
                        y,
                        zorder=zorder - 0.001,
                        color=self.bond_edgecolor,
                        lw=self.scaled_bond_width,
                        solid_capstyle='round')

                # draw bond fill
                self.ax.plot(
                        x,
                        y,
                        zorder=zorder,
                        color=self.bond_color,
                        lw=self.scaled_bond_fill,
                        solid_capstyle='round')

        # get bond lines if first time drawing bonds
        if not self.bond_objs:
            self.bond_objs = self.ax.lines

        # add bonds to list of drawn items
        if 'bonds' not in self._drawn:
            self._drawn.append('bonds')

    def draw_labels(self, labels=None, force=False):
        # can pass in new label types when calling draw_labels
        if labels is not None:
            self.__check_labels__(labels)

        # determine if drawing or updating
        update_labels = bool(len(self.label_objs))

        if self.labels is not None:
            for i, a in enumerate(self.atoms):
                if update_labels:
                    self.label_objs[i].set_text(self._labels[i])
                    self.label_objs[i].set_x(a.x)
                    self.label_objs[i].set_y(a.y)
                    self.label_objs[i].set_zorder(a.z + 0.001)
                else:
                    ann = self.ax.annotate(
                                self._labels[i],
                                (a.x, a.y),
                                zorder=a.z + 0.001,
                                ha='center',
                                va='center',
                                fontsize=15)
                    # add annotation to label_objs list
                    self.label_objs.append(ann)

        if 'labels' not in self._drawn:
            self._drawn.append('labels')

    def draw_legend(self, leg_order=None, max_ms=16, force=False):
        if self.block_legend:
            print("Warning: Cannot add legend unless atoms "
                  "are colored by type.")
            if 'legend' in self._drawn:
                self.remove_legend()
            return

        # return if legend is already drawn
        if 'legend' in self._drawn:
            # redraw legend if different leg_order
            if isinstance(leg_order, str) and leg_order != self.leg_order:
                self._drawn.remove('legend')
            elif not force:
                return

        # default leg_order is alphabetical
        if leg_order is None:
            leg_order = 'alphabetical'

        # track current legend order
        self.leg_order = leg_order

        # create an ordered, unique list of atom types
        symbols = sorted(set(self.atoms.symbols))

        # use custom legend order if given
        if isinstance(leg_order, str) and leg_order != 'alphabetical':
            # order legend by size
            if leg_order in ['size', 'size_r']:
                leg_order = sorted(
                                    symbols,
                                    key=lambda z: list(covalent_radii)[
                                        chemical_symbols.index(z)],
                                    reverse=leg_order == 'size')
            else:
                leg_order = [leg_order]
        if isinstance(leg_order, list) or isinstance(leg_order,
                                                     np.ndarray):
            # any types not in leg_order will be appended in
            # alphabetical order
            leg_order = (list(leg_order) +
                         [s for s in symbols if s not in leg_order])

            symbols = sorted(symbols, key=lambda z: leg_order.index(z))

        # get an atom object of each type
        all_symbols = list(self.atoms.symbols)
        a_objs = [self.atoms[all_symbols.index(s)] for s in symbols]

        # calculate sizes of each atom
        sizes = np.array([utils.angstrom_to_axunits(
                          covalent_radii[a.number] * self.scale, self.ax) * 2
                          for a in a_objs])

        # normalize sizes such that largest atom
        # has size of <max_ms>
        sizes = sizes * (max_ms / sizes.max())

        # create legend objects
        leg = [Line2D([0], [0], marker='o', ls='',
                      markerfacecolor=self._colors[a.index],
                      markeredgecolor='k',
                      markersize=ms)
               for a, ms in zip(a_objs, sizes)]

        # create legend
        self.legend = self.ax.legend(
            leg,
            symbols,
            frameon=False,
            prop=dict(size=11, weight='bold'),
            handletextpad=np.sqrt(utils.angstrom_to_axunits(0.01, self.ax)),
            borderpad=0,
            borderaxespad=0,
            columnspacing=0,
            markerscale=1,
            labelspacing=np.sqrt(utils.angstrom_to_axunits(0.08, self.ax)),
            loc='center left',
            framealpha=0,
            # up to ten atom types per column
            ncol=(len(leg) // 10) + 1,
            bbox_to_anchor=(1, 0.5))

        # tight layout after making legend
        self.fig.tight_layout()

        # add legend to list of drawn items
        if 'legend' not in self._drawn:
            self._drawn.append('legend')

    def draw_colorbar(self, force=False):
        if self.block_colorbar:
            print("Warning: Unable to draw colorbar - "
                  "values needed in place of colors")
            if 'colorbar' in self._drawn:
                self.remove_colorbar()
            return

        # do not redraw colorbar
        if 'colorbar' in self._drawn and not force:
            return

        # add in colorbar axis
        gs = self.fig.add_gridspec(1, 2, width_ratios=[30, 1])
        self.cb_ax = self.fig.add_subplot(gs[0, 1])

        cb = matplotlib.colorbar.ColorbarBase(
                                        self.cb_ax,
                                        cmap=self.cmap,
                                        norm=self.cb_norm)
        # set font size
        cb.ax.tick_params(labelsize=11)

        # tighten layout after drawing colorbar
        self.fig.tight_layout()

        # add colorbar to list of drawn items
        if 'colorbar' not in self._drawn:
            self._drawn.append('colorbar')

    def remove(self, items=['atoms', 'bonds', 'labels']):
        """
        Removes multiple object types at once
        - Default: removes atoms, bonds, and labels (if any specified)

        KArgs:
        - items (list): list of object types to remove (as str)
                        items='all' will remove everything
                        OPTIONS:
                        - atoms
                        - bonds
                        - labels
                        - legend
                        - colorbar
        """
        if items == 'all':
            items = self._drawn.copy()
        for item in items:
            method = 'remove_' + item
            if method in self.__dir__():
                getattr(self, method)()

    def remove_atoms(self):
        if 'atoms' in self._drawn:
            for a in self.atom_objs[::-1]:
                a.remove()
            self.atom_objs = []
            self._drawn.remove('atoms')

    def remove_bonds(self):
        if 'bonds' in self._drawn:
            for b in self.bond_objs[::-1]:
                b.remove()
            self._drawn.remove('bonds')

    def remove_labels(self):
        if 'labels' in self._drawn:
            for lab in self.label_objs[::-1]:
                lab.remove()
            self.label_objs = []
            self._drawn.remove('labels')

    def remove_legend(self):
        if self.legend is not None:
            self.legend.remove()
            self.legend = None
            self._drawn.remove('legend')
            self.fig.tight_layout()

    def remove_colorbar(self):
        if isinstance(self.cb_ax, plt.Axes):
            self.fig.delaxes(self.cb_ax)
            self.cb_ax = None
            self._drawn.remove('colorbar')
            self.fig.tight_layout()

    def smart_rotate(self):
        new_pos, self.transform = utils.pca(self.atoms.positions,
                                            return_transform=True)
        self.atoms.positions = new_pos

        # reinitialize figure to rescale
        self.init_fig()

        # redraw objects
        self.update(force=True)

    def add_param(self, value):
        if value not in self.fig_params:
            self.fig_params.append(value)

    def adjust_padding(self, padding):
        self.init_fig(padding)
        self.update(force=True)

    def anchor(self, i):
        """
        Translates atoms[<i>] to origin
        - <i> can also be cartesian coordinates

        Args:
        - i (int): index of atom that should be translated to origin
        - i (list): cartesian coordinates that should be translated to origin
        """
        try:
            if len(i) == 3:
                self.atoms.positions -= list(i)
            else:
                raise ValueError('not cartesian coordinates')
        except:
            if 0 <= i < len(self):
                self.atoms.positions -= self.atoms[i].position
                # recalculate axis boundaries
                self.init_fig()
                self.add_param('anchor')
            else:
                print('Invalid index given to anchor')

    def rotate(self, angle, rot_axis=None):
        if rot_axis is None:
            if self.rot_axis is None:
                rot_axis = 'y'
            else:
                rot_axis = self.rot_axis
        try:
            self.atoms.rotate(angle, v=rot_axis)
            self.rot_axis = rot_axis
            self.update()
        except:
            print("Unable to rotate - make sure angle "
                  "and rot_axis (if given) are valid")

    def show(self):
        if isinstance(self.fig, plt.Figure):
            self.fig.show()
        else:
            raise TypeError('init_fig must first be called')

    def save(self, path=None, overwrite=False, max_px=600, transparent=False):
        dpi = max_px / 5

        if path is None:
            path = '%s-%s.png' % (self.atoms.get_chemical_formula(),
                                  '-'.join(self.fig_params))

        if not overwrite:
            path = utils.avoid_overwrite(path)

        if path.endswith('.png'):
            ram = io.BytesIO()
            self.fig.savefig(ram, format='png', dpi=dpi,
                             transparent=transparent)
            ram.seek(0)
            im = Image.open(ram)

            # scale number of unique colors in png with number of atom types
            ncolors = self.atom_types + 15

            im2 = im.convert('P', palette=Image.ADAPTIVE, colors=ncolors)

            im2.save(path, optimize=True)
        else:
            self.fig.savefig(path, transparent=transparent, dpi=dpi)

    def update(self, force=False, recalc_bonds=None):
        if isinstance(recalc_bonds, bool):
            self.recalc_bonds = recalc_bonds

        # call all draw methods
        self.draw(self._drawn, force=force)

    def rot_gif(self, path=None, fps=20, loop_time=6,
                max_px=600, rot_axis='y', save_frames=False, overwrite=False,
                optimize_gif=False):
        if path is None:
            path = 'C:\\users\\mcowa\\desktop\\test.gif'

        # ensure path ends with .gif
        if not path.endswith('.gif'):
            path += '.gif'

        # avoid overwriting previous gif
        if not overwrite:
            path = utils.avoid_overwrite(path)

        # save frames in *_frames folder
        if save_frames or optimize_gif:
            frame_path = os.path.join(
                            os.path.dirname(path),
                            os.path.basename(path).strip('.gif')) + '_frames'

            # avoid overwriting previous frame_path
            frame_path = utils.avoid_overwrite_dir(frame_path)
            os.mkdir(frame_path)

            # total number of frames
            frames = fps * loop_time

            # angle to rotate between frames (in degrees)
            rot = 360 / frames

            # iteratively save all frames
            print(' ' * 50, end='\r')
            for i in range(frames):
                print(' Saving frame %03i' % (i + 1), end='\r')
                self.save(os.path.join(frame_path,
                                       self.name + '_%03i.png' % (i + 1)))
                self.rotate(rot)
        # else make gif with matplotlib.animation and image magick
        else:
            animation = self.rot_animation(fps=fps, loop_time=loop_time,
                                           rot_axis=rot_axis, gif=True)

            # initialize imagemagick writer
            if anim.writers.is_available('imagemagick'):
                writer = anim.ImageMagickWriter(fps=fps)
            else:
                # ImageMagick must be used
                raise ImportError("ImageMagick must be installed "
                                  "to create GIF")

            # save gif
            animation.save(path, writer=writer, dpi=max_px / 5)

        if optimize_gif:
            print(' Creating optimized gif...', end='\r')
            subprocess.call(['magick', 'convert', '-delay', str(100 / fps),
                             os.path.join(frame_path, '*.png'), path])

            # delete frames
            for f in os.listdir(frame_path):
                os.remove(os.path.join(frame_path, f))
            os.removedirs(frame_path)

    def rot_animation(self, fps=20, loop_time=6, rot_axis='y', gif=False):
        frames = fps * loop_time
        rot = 360 / frames
        dig_str = '%0{}i'.format(len(str(frames)))

        # adjust axis limits based on rot_axis
        self.rot_axis = rot_axis
        self.init_fig()
        self.update(force=True)
        self.fig.tight_layout()

        def next_step(i):
            # print out progress if building gif
            if gif:
                print(' ' * 50, end='\r')
                if (i + 1) == frames:
                    print(' Wrapping things up..', end='\r')
                else:
                    print(' Building frame: ' + dig_str % (i + 2), end='\r')

            # rotate atoms
            self.rotate(rot, rot_axis)

        # build frames
        animation = anim.FuncAnimation(self.fig,
                                       next_step,
                                       frames=frames,
                                       interval=100/fps)
        return animation

if __name__ == '__main__':
    path = 'C:\\users\\mcowa\\desktop\\autu2oct_cn.xyz'
    # a = ase.io.read(path)
    a = ase.build.molecule('C60')
    molecule = MolFig(a)
    molecule.draw()
    # animation = molecule.rot_gif(optimize_gif=True)
    plt.show()
