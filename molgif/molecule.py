from __future__ import division, print_function
import molgif.utils as utils
import os
import io
import random
import re
import subprocess
import platform
import shutil
import tempfile
import copy
import ase
from ase.data import covalent_radii, chemical_symbols
from ase.data.colors import jmol_colors
import numpy as np
from PIL import Image
import matplotlib
import matplotlib.colors as mcolors
from matplotlib.lines import Line2D
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import matplotlib.animation as anim

# ensure that ImageMagick writer is found (and not Windows' convert.exe)
if platform.system().lower().startswith('windows'):
    if 'rcParams' in anim.__dir__():
        rc = anim.rcParams
    else:
        rc = matplotlib.rcParams

    if rc['animation.convert_path'].endswith(('convert', 'convert.exe')):
        rc['animation.convert_path'] = 'magick'


class Molecule(object):
    def __init__(self, atoms, scale=0.7, name=None, square=False,
                 colors=None, labels=None, rot_axis=None,
                 alphas=None, bond_color='white', bond_edgecolor='black',
                 cb_min=None, cb_max=None, center_data=False, cmap=cm.bwr_r,
                 draw=['atoms', 'bonds']):
        """
        Molecule figure object that enables dynamic drawing, manipulating,
        and animating a matplotlib Figure

        Args:
        - atoms (ase.Atoms || str): atoms object or path to geometry file

        KArgs:
        - scale (float): atomic radii are set to <scale> * covalent radii
                         (Default: 0.7)

        - name (str): name of molecule
                      (Default: chemical formula)

        - square (bool): if True, figure will be square
                         (Default: False)

        - colors (str iterable dict): specify atom colors with str, dict,
                                      or values which will use the cmap
                                      - 'blue': all atoms blue
                                      - ['blue', 'white', ...]: label each atom
                                      - [0, 1.1, -2.3...]: cmap used to color
                                      - {'Au': 'purple'}: use dict to color by
                                        by atom type - types not given use jmol
                                      (Default: None -> jmol colors used)

        - labels (str | iterable): type or list of labels to add to atoms
                                   - 'symbols': uses chemical symbol
                                   - 'symbols-noh': same as symbols, but no H
                                   - 'values': uses values from colors KArg
                                   - 'charge': uses initial_charges from atoms
                                   - [lab1, lab2, ...]
                                   (Default: None)

        - rot_axis (str): specify axis to rotate about
                          - x (left-to-right), y (bot-to-top), or
                            z (counter-clockwise)
                          - can be: 'y' | 'x' | 'z',
                          - can also add a '-' to change direction!
                            e.g. '-z' causes clockwise rotation
                          (Default: 'y')

        - alphas (int iterable dict): set transparency of atom type
                                      or atom indices
                                      (Default: None)

        - bond_color (str): specify color of bonds
                            (Default: white)

        - bond_edgecolor (str): specify edgecolor (border) of bonds
                                (Default: black)

        - cb_min (float): min value limit for colorbar
                          (Default: calc based on values)

        - cb_max (float): max value limit for colorbar
                          (Default: calc based on values)

        - center_data (bool): if True, colors are centered about middle of cmap
                              - ensures (-) and (+) values are different color
                              - ex) for RdBu cmap, 0 = 'white'
                              (Default: False)

        - cmap (ColorMap): cmap to be used if values are given for colors
                           (Default: matplotlib.cm.bwr_r)

        - draw (list): items to immediately draw after initialization
                       (Default: ['atoms', 'bonds'])
        """
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
        self._rot_axis = None
        if isinstance(rot_axis, str) and re.match('-?[xyz]', rot_axis.lower()):
            self._rot_axis = rot_axis

        # padding around atoms
        self.padding = None

        # position transformation matrix
        # default to no transformation
        self.transform = np.eye(3)

        # square figure (px-x = px-y)
        self.square = square

        # set atoms object
        self.atoms = utils.path2atoms(atoms)
        self.atoms = self.atoms.copy()

        # center atoms to COP
        self.atoms.positions -= self.atoms.positions.mean(0)

        # set name
        if name is None:
            # default name is chemical formula
            self.name = self.atoms.get_chemical_formula()
        else:
            self.name = str(name).strip('.gif')

        # track atom objects (track for redrawing atoms and bonds)
        self.pos = self.atoms.positions.copy()
        self.bond_pos = self.atoms.positions.copy()

        # border width of atoms (in Angstrom)
        self.atom_borderwidth = 0.05

        # atom size scale (covalent_radii * scale)
        self._scale = abs(float(scale))

        # track to see when scale changes
        self._prevscale = self._scale

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

        # default label size: 15
        self._label_size = 15

        # flags used to block drawing colorbar and legend when not applicable
        self.block_colorbar = False
        self.block_legend = False

        self._check_labels(labels)

        # legend attributes (track order)
        self.legend = None
        self.leg_order = 'size'

        # alpha value of atoms
        self._alphas = [1.] * len(self.atoms)

        if alphas is not None:
            self.alphas = alphas

        # color attributes
        self._colors = [jmol_colors[a.number] for a in self.atoms]

        # track if colors have changed
        self._prevcolors = None

        # check input colors and see if they can be used with legend/colorbar
        self._check_colors(colors)

        # initialize figure
        self.init_fig()

        # draw any initial items given in <draw>
        self.draw(draw)

    def __len__(self):
        """Returns the number of atoms"""
        return len(self.atoms)

    def __repr__(self):
        return 'molgif.Molecule(%s)' % self.atoms.get_chemical_formula()

    def __getitem__(self, i):
        return self.atoms[i]

    @property
    def alphas(self):
        return self._alphas

    @alphas.setter
    def alphas(self, value):
        changed = True
        if isinstance(value, dict):
            self._alphas = self._dict_to_list(value, self._alphas)[0]
        elif ((isinstance(value, list) or isinstance(value, np.ndarray)) and
              len(value) == len(self)):
            self._alphas = list(value)
        elif isinstance(value, int) or isinstance(value, float):
            self._alphas = [value] * len(self)
        else:
            changed = False
        if changed and self.fig is not None:
            self.draw_atoms(force=True)
            if 'bonds' in self.drawn:
                self.draw_bonds()

    # colors getter
    @property
    def colors(self):
        """
        Can dynamically change color of atoms
        OPTIONS:
        - None: default jmol_colors
        - 'blue': all atoms blue
        - ['blue', 'white', ...]: label each atom
        - [0, 1.1, -2.3...]: cmap used to color
        - {'Au': 'purple'}: use dict to color by
          by atom type - types not given use jmol
        """
        return self._colors

    @colors.setter
    def colors(self, value):
        self._check_colors(value)

    # labels getter
    @property
    def labels(self):
        """
        Can dynamically change labels on atoms
        OPTIONS:
        - None: no labels
        - 'symbols': uses chemical symbol
        - 'symbols-noh': same as symbols, but no H
        - 'values': uses values from colors KArg
        - 'charges': uses initial_charges from atoms
        - [lab1, lab2, ...]
        """
        return self._labels

    @labels.setter
    def labels(self, value):
        self._check_labels(value)

    @property
    def label_size(self):
        return self._label_size

    @label_size.setter
    def label_size(self, value):
        if isinstance(value, str) and not value.isdigit():
            value = value.lower()
            if re.match('(x{0,2}-)?(small|large)', value) or value == 'medium':
                self._label_size = value
        else:
            try:
                self._label_size = int(value)
            except ValueError:
                print('Label size must be point size or relative value '
                      '(small, medium, etc.)')

    @property
    def scale(self):
        """
        """
        return self._scale

    @scale.setter
    def scale(self, value):
        try:
            value = abs(float(value))
            if value != self._scale:
                self._scale = value
                if 'atoms' in self._drawn:
                    self.draw_atoms()
        except (ValueError, TypeError):
            print('Scale must be a number')

    # tracks items that have been drawn
    @property
    def drawn(self):
        """
        List of items currently drawn on plot
        - pptions include: atoms, bonds, labels, colorbar, and legend
        """
        return self._drawn

    @property
    def rot_axis(self):
        """
        Can dynamically change rotation axis
        - figure dimensions are adjusted based on rot_axis
        OPTIONS:
        - None: draw figure based only on current view
        - x (top-to-bot)
        - y (left-to-right)
        - z (counterclockwise)
        - can also add a '-' to change direction!
          e.g. '-z' causes clockwise rotation
        """
        return self._rot_axis

    @rot_axis.setter
    def rot_axis(self, value):
        if (isinstance(value, str) and
           re.match('-?[xyz]', value.lower()) and
           value != self._rot_axis):
            self._rot_axis = value
            self.redraw()

    def init_fig(self, padding=2):
        """
        Create initial figure and axis objects
        - determines figure size and axis limits
        - also makes colorbar axis if needed

        KArgs:
        padding (float): amount of padding to give around molecule
                         - NOTE: padding value is dependent on size of system
                         (Default: 2)
        """
        if self.padding is None:
            self.padding = padding

        # calculate figure size
        # calculate axis limits (include offset as buffer)
        self.fig_size, self.xlim, self.ylim = utils.get_fig_bounds(
                                                    self.atoms,
                                                    rot_axis=self._rot_axis,
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

        # adjust axis to use entire fig space
        self.fig.subplots_adjust(left=0, right=1, bottom=0, top=1)

        # set aspect ratio to 1
        self.ax.set_aspect(1)

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

    def clear_figure(self):
        """
        Removes everything from figure
        same as remove('all')
        """
        self.remove('all')

    def draw(self, items=None, force=False):
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

        - force (bool): if True, all items are forced to redraw
                        (Default: False)
        """
        # defaults to drawn list
        if items is None:
            items = self._drawn.copy()

        # items must be a list
        if not isinstance(items, list) or items == []:
            return
        # draw elements in item list
        for item in items:
            method = 'draw_' + item
            if method in self.__dir__():
                getattr(self, method)(force=force)
            else:
                print(item, 'can not be drawn.')

    def draw_atoms(self, force=False):
        """
        Draw atoms or update atom positions

        KArgs:
        - force (bool): if True, atoms are forced to redraw
                        (Default: False)
        """

        # remove atoms if forcing to redraw
        if force:
            [i.remove() for i in self.atom_objs]
            self.fig.canvas.draw()
            self.atom_objs = []

        # determine if drawing or updating
        update = False
        if (self.pos != self.atoms.positions).any():
            update = True
            self.pos = self.atoms.positions.copy()
        if self._scale != self._prevscale:
            update = True
            self.radii = covalent_radii[self.atoms.numbers] * self._scale
            self._prevscale = self._scale
            if 'bonds' in self._drawn:
                self.draw_bonds(force=True)

        if self._colors != self._prevcolors:
            update = True
            self._prevcolors = self._colors

        # don't redraw if params haven't changed
        if self.atom_objs:
            if not (update or force):
                return
        else:
            force = True

        # create atom patches and add to ax
        for i, a in enumerate(self.atoms):
            # update atom objects
            if update and not force:
                self.atom_objs[i].set_facecolor(self._colors[i])
                self.atom_objs[i].radius = self.radii[i]
                self.atom_objs[i].center = (a.x, a.y)
                self.atom_objs[i].zorder = a.z
            # create atom object patches
            else:
                circ = plt.Circle(
                            (a.x, a.y),
                            radius=self.radii[i],
                            facecolor=self._colors[i],
                            alpha=self._alphas[i],
                            edgecolor='k',
                            linewidth=self.scaled_borderwidth,
                            zorder=a.z)
                # add patch to atom_objs list and to axis
                self.atom_objs.append(circ)
                self.ax.add_artist(circ)

        # add atoms to list of drawn items
        if 'atoms' not in self._drawn:
            self._drawn.append('atoms')

    def draw_bonds(self, bond_color=None, bond_edgecolor=None,
                   recalc_bonds=None, force=False):
        """
        Draw bonds or update bond positions

        KArgs:
        - bond_color (str): specify color of bonds
                        (Default: white)

        - bond_edgecolor (str): specify edgecolor (border) of bonds
                                (Default: black)

        - recalc_bonds (bool): if True, bonding is recalculated
                               (Default: False)

        - force (bool): if True, bonds are forced to redraw
                        (Default: False)
        """
        # determine if drawing or updating bonds
        update = bool(len(self.bond_objs))

        # if new bond colors, force a redraw
        if self.bond_color != bond_color and bond_color is not None:
            self.bond_color = bond_color
            force = True
        if (self.bond_edgecolor != bond_edgecolor and
           bond_edgecolor is not None):
            self.bond_edgecolor = bond_edgecolor
            force = True

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
            # buffer = 0.078 * np.exp(-zdist)
            buffer = 0.1 * np.exp(-zdist)

            xy1 += buffer * normvec
            xy2 -= buffer * normvec

            # get x - y positions for bonds
            x = [xy1[0], xy2[0]]
            y = [xy1[1], xy2[1]]

            # use avg z for zorder
            zorder = (xy1[2] + xy2[2]) / 2

            # set alpha value of bond to lowest alpha atom
            alpha = min(self.alphas[b[0]], self.alphas[b[1]])

            if update:
                # update bond outline
                self.bond_objs[2 * i].set_data((x, y))
                self.bond_objs[2 * i].set_zorder(zorder - 0.001)
                self.bond_objs[2 * i].set_alpha(alpha)

                # update bond
                self.bond_objs[2 * i + 1].set_data((x, y))
                self.bond_objs[2 * i + 1].set_zorder(zorder)
                self.bond_objs[2 * i + 1].set_alpha(alpha)
            else:
                # draw thinner white line over black line
                # to create bordered bonds
                self.ax.plot(
                        x,
                        y,
                        zorder=zorder - 0.001,
                        color=self.bond_edgecolor,
                        lw=self.scaled_bond_width,
                        alpha=alpha,
                        solid_capstyle='round')

                # draw bond fill
                self.ax.plot(
                        x,
                        y,
                        zorder=zorder,
                        color=self.bond_color,
                        alpha=alpha,
                        lw=self.scaled_bond_fill,
                        solid_capstyle='round')

        # get bond lines if first time drawing bonds
        if not self.bond_objs:
            self.bond_objs = self.ax.lines

        # add bonds to list of drawn items
        if 'bonds' not in self._drawn:
            self._drawn.append('bonds')

    def draw_labels(self, labels=None, label_size=None, force=False):
        """
        Draw labels or update label position

        KArgs:
        - labels (str | iterable): type or list of labels to add to atoms
                                   - 'symbols': uses chemical symbol
                                   - 'symbols-noh': same as symbols, but no H
                                   - 'values': uses values from colors KArg
                                   - 'charges': uses initial_charges from atoms
                                   - [lab1, lab2, ...]
                                   (Default: refers to self.labels)

        - force (bool): if True, labels are forced to redraw
                        (Default: False)
        """
        # can pass in new label types when calling draw_labels
        if labels is not None:
            self._check_labels(labels)

        # make sure labels exist
        if self.labels is None:
            print('No labels to draw.')
            return

        # determine if drawing or updating
        update_labels = bool(len(self.label_objs))

        if label_size is not None:
            self.label_size = label_size

        for i, a in enumerate(self.atoms):
            if update_labels:
                self.label_objs[i].set_text(self._labels[i])
                self.label_objs[i].set_x(a.x)
                self.label_objs[i].set_y(a.y)
                self.label_objs[i].set_zorder(a.z + 0.001)
                self.label_objs[i].set_size(self._label_size)
            else:
                ann = self.ax.annotate(
                            self._labels[i],
                            (a.x, a.y),
                            zorder=a.z + 0.001,
                            ha='center',
                            va='center',
                            fontsize=self._label_size)
                # add annotation to label_objs list
                self.label_objs.append(ann)

        if 'labels' not in self._drawn:
            self._drawn.append('labels')

    def draw_legend(self, leg_order=None, max_ms=16, force=False):
        """
        Draw or update legend

        KArgs:
        - leg_order (list | str): if given, use it to order the legend
                                  - can also give str of single atom type
                                  - 'size': largest to smallest
                                  - 'size_r': smallest to largest
                                  - 'alpha': alphabetical order
                                  - 'atomic': atomic number
                                  - None: atomic number
                                  (Default: size)
        - max_ms (int): scales legend such that largest atom type
                        is represented with markersize = <max_ms>
                        (Default: 16 pts)

        - force (bool): if True, legend is forced to redraw
                        (Default: False)
        """
        if self.block_legend:
            print("Warning: Cannot add legend unless atoms "
                  "are colored by type.")
            if 'legend' in self._drawn:
                self.remove_legend()
            return

        # return if legend is already drawn
        if 'legend' in self._drawn:
            # redraw legend if different leg_order
            if (isinstance(leg_order, str) and
               leg_order != self.leg_order):
                self._drawn.remove('legend')
                self.leg_order = leg_order.lower()
            elif not force:
                return

        elif leg_order is not None:
            self.leg_order = leg_order.lower()

        # create an ordered, unique list of atom types
        symbols = sorted(set(self.atoms.symbols))

        # use custom legend order if given
        if isinstance(self.leg_order, str):
            # symbols are already in alphabetical order
            if self.leg_order == 'alpha':
                pass
            # order legend by size
            elif self.leg_order in ['size', 'size_r']:
                symbols = sorted(
                            symbols,
                            key=lambda z: list(covalent_radii)[
                                chemical_symbols.index(z)],
                            reverse=self.leg_order == 'size')
            # sort by atomic number
            elif self.leg_order == 'number':
                symbols = sorted(symbols, key=chemical_symbols.index)
            else:
                raise ValueError("Unable to handle leg_order "
                                 "= %s" % self.leg_order)
        if isinstance(self.leg_order, list) or isinstance(self.leg_order,
                                                          np.ndarray):
            # any types not in leg_order will be appended in
            # alphabetical order
            leg_list = (list(self.leg_order) +
                        [s for s in symbols if s not in self.leg_order])

            symbols = sorted(symbols, key=lambda z: leg_list.index(z))

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

        # readjust figure borders so legend does not overlap
        cut = max(0.92 - 0.12 * (len(leg) // 10), 0.6)
        self.fig.subplots_adjust(right=cut)

        # create legend
        self.legend = self.fig.legend(
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
            bbox_to_anchor=(0.99, 0.5),
            loc='center right',
            framealpha=0,
            # up to nine atom types per column
            ncol=(len(leg) // 10) + 1)

        # add legend to list of drawn items
        if 'legend' not in self._drawn:
            self._drawn.append('legend')

    def draw_colorbar(self, cb_min=None, cb_max=None, center_data=None,
                      force=False):
        """
        Draw or update colorbar

        KArgs:
        - cb_min (float): max value limit for colorbar
                          (Default: refer to self.cb_min)

        - cb_max (float): max value limit for colorbar
                          (Default: refer to self.cb_max)

        - center_data (bool): if True, colors are centered about middle of cmap
                              - ensures (-) and (+) values are different color
                              - ex) for RdBu cmap, 0 = 'white'
                              (Default: refer to self.center_data)

        - force (bool): if True, colorbar is forced to redraw
                        (Default: False)
        """
        if self.block_colorbar:
            print("Warning: Unable to draw colorbar - "
                  "values needed in place of colors")
            if 'colorbar' in self._drawn:
                self.remove_colorbar()
            return

        # do not redraw colorbar
        if 'colorbar' in self._drawn and not force:
            return

        # see if new min or max color scale (or center_data)
        remake_color = False
        if cb_min is not None and cb_min != self.cb_min:
            self.cb_min = float(cb_min)
            remake_color = True
        if cb_max is not None and cb_max != self.cb_max:
            self.cb_max = float(cb_max)
            remake_color = True
        if center_data is not None and center_data != self.center_data:
            self.center_data = center_data
            remake_color = True

        # rescale colors
        if remake_color:
            print('checking colors!')
            self._check_colors(self._cb_values)

        # add in colorbar axis
        gs = self.fig.add_gridspec(1, 2, width_ratios=[30, 1])
        self.cb_ax = self.fig.add_subplot(gs[0, 1])
        mappable = cm.ScalarMappable(norm=self.cb_norm, cmap=self.cmap)
        mappable.set_array([])
        self.fig.colorbar(mappable, cax=self.cb_ax)

        # readjust figure borders so colorbar is not cut off
        self.fig.subplots_adjust(wspace=0, right=0.85, bottom=0.05, top=0.95)

        # set font size
        self.cb_ax.tick_params(labelsize=11)

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
        """
        Removes atoms from figure
        """
        if 'atoms' in self._drawn:
            for a in self.atom_objs[::-1]:
                a.remove()
            self.atom_objs = []
            self._drawn.remove('atoms')
            self.fig.canvas.draw()

    def remove_bonds(self):
        """
        Removes bonds from figure
        """
        if 'bonds' in self._drawn:
            for b in self.bond_objs[::-1]:
                b.remove()
            self._drawn.remove('bonds')

    def remove_labels(self):
        """
        Removes labels from figure
        """
        if 'labels' in self._drawn:
            for lab in self.label_objs[::-1]:
                lab.remove()
            self.label_objs = []
            self._drawn.remove('labels')

    def remove_legend(self):
        """
        Removes legend from figure
        """
        if self.legend is not None:
            self.legend.remove()
            self.legend = None
            self._drawn.remove('legend')
            self.init_fig()

    def remove_colorbar(self):
        """
        Removes colorbar from figure
        """
        if isinstance(self.cb_ax, plt.Axes):
            self.fig.delaxes(self.cb_ax)
            self.cb_ax = None
            self._drawn.remove('colorbar')
            self.init_fig()

    def smart_rotate(self, opt_angle=False):
        """
        Applies "smart" rotation to molecule
        - attempts to find the best viewing angle using PCA (love this idea)

        KArgs:
        - opt_angle: Finds angle to rotate atoms about z-axis to maximize
                     x and y distance of atoms
                     - molecule is "squared off", which may give a better
                       overall presentation for some systems
        """
        # TODO: correct transformation matrix from extra step
        # in new smart_rotate
        self.atoms = utils.smart_rotate_atoms(self.atoms, opt_angle=opt_angle)

        # new_pos, self.transform = utils.pca(self.atoms.positions.copy(),
        #                                     return_transform=True)

        # self.atoms.positions = new_pos.copy()

        # redraw figure
        self.redraw()

    def custom_rotate(self, rot_commands):
        """
        Applies an ordered list of rotation commands to molecule
        - list should have the form [angle1, x|y|z, angle2, x|y|z, ...]
        - rotations are applied sequentially and different orders will lead
          to different results

        Args:
        - rot_commands (list): list of ordered rotation commands to apply
                               - [angle, (x|y|z), angle, (x|y|z), ...]
        """
        # store current rotation axis
        rot_axis = self.rot_axis

        # apply the rotations sequentially
        for angle, axis in zip(rot_commands[::2], rot_commands[1::2]):
            self.rotate(angle, axis)

        # set rotation axis back to original
        self.rot_axis = rot_axis

    def add_param(self, value):
        if value not in self.fig_params:
            self.fig_params.append(value)

    def adjust_padding(self, padding):
        """
        Adjust padding between molecule and figure edges

        Args:
        padding (float): amount of padding to give around molecule
                         - NOTE: padding value is dependent on size of system
        """
        self.init_fig(padding)
        self.update(force=True)

    def anchor(self, i):
        """
        Translates atoms[<i>] to origin
        - <i> can also be cartesian coordinates

        Args:
        - i (str): atom element to anchor - matches first element in list
        - i == 'center': closest atom to center of position (COP)
        - i (int): index of atom that should be translated to origin
        - i (list): cartesian coordinates that should be translated to origin
        """
        try:
            if len(i) == 3:
                self.atoms.positions -= list(i)
            else:
                raise ValueError('not cartesian coordinates')
        except Exception:
            if isinstance(i, str):
                if i.title() in self.atoms.symbols:
                    i = np.where(self.atoms.symbols == i.title())[0][0]
                elif i.isdigit():
                    i = int(i)
                elif i.lower() == 'center':
                    # ensure that molecule's center of position is at origin
                    self.atoms.positions -= self.atoms.positions.mean(0)

                    # don't need to square root since order won't change
                    dists = (self.atoms.positions**2).sum(1)

                    # find closest atom to COP
                    i = np.where(dists == dists.min())[0][0]

            if 0 <= i < len(self):
                # store index of anchored atom
                self.anchored_atom = i

                # shift anchored atom to origin
                self.atoms.positions -= self.atoms[i].position

                # recalculate axis boundaries
                self.init_fig()
                self.add_param('anchor')
            else:
                print('Invalid index given to anchor')

    def center(self):
        """
        Center atoms by center-of-position (COP)
        """
        self.anchor(self.atom.positions.mean(0))

    def rotate(self, angle, rot_axis=None):
        """
        Rotates molecule and recalculates figure dimensions

        Args:
        - angle (float): angle to rotate in degrees

        - rot_axis (str): specify axis to rotate about
                          - x (left-to-right), y (bot-to-top), or
                            z (counter-clockwise)
                          - can be: 'y' | 'x' | 'z',
                          - can also add a '-' to change direction!
                            e.g. '-z' causes clockwise rotation
                          (Default: refer to self.rot_axis or 'y')
        """
        self.rot_axis = rot_axis
        if self.rot_axis is None:
            self.rot_axis = 'y'

        try:
            self.atoms.rotate(float(angle), v=self.rot_axis,
                              center=(0, 0, 0))
            self.update()
        except Exception:
            print("Unable to rotate - make sure angle "
                  "and rot_axis (if given) are valid")

    def redraw(self):
        """
        Redraws figure, including new fig and ax objects
        """
        to_redraw = self.drawn.copy()
        self.clear_figure()
        plt.close(self.fig)
        self.fig = None
        self.ax = None
        self.init_fig()
        self.draw(to_redraw)

    def show(self):
        """
        Uses plt.show() to show figure to user
        - Redraws figure once figure is closed
        """
        if isinstance(self.fig, plt.Figure):
            plt.show()
            self.redraw()
        else:
            raise TypeError('init_fig must first be called')

    def save(self, path=None, overwrite=False, max_px=600, transparent=False,
             optimize=False, quiet=False):
        """
        Save the current figure
        - default = *.png
        - can save as *.svg vector graphic
        - supports all extensions supported by plt.Figure.savefig

        KArgs:
        - path (str): path to save image
                      (Default: refer to self.name)

        - overwrite (bool): if False, '-<integer>' is added to end of name
                            to avoid overwriting
                            (Default: False)

        - max_px (int): sets pixel count for longest dimension in image
                        (Default: 600)

        - transparent (bool): if True, saved with transparent background
                              (Default: False)

        - optimize (bool): optimizes png if no colorbar is drawn
                           - NOTE: could result in some loss of quality
                           (Default: False)
        - quiet (bool): if True, supress output to console
        """
        # inform user of attempt to save
        if not quiet:
            print('Saving %s...' % path, end='\r')

        # calculate dpi based on max_px and 5" max dimension of mpl figure
        dpi = max_px / 5

        if path is None:
            path = self.name + '.png'
            # path = '%s-%s.png' % (self.name,
            #                       '-'.join(self.fig_params))

        if not overwrite:
            path = utils.avoid_overwrite(path)

        if (optimize and
           path.endswith('.png') and
           'colorbar' not in self._drawn):
            ram = io.BytesIO()
            self.fig.savefig(ram, format='png', dpi=dpi,
                             transparent=transparent)
            ram.seek(0)
            im = Image.open(ram)

            # define number of colors in png
            # scale by number of unique colors in molecule
            uniquecolors = len(np.unique(self._colors, axis=0))
            ncolors = int(15 + 2 * uniquecolors)

            im2 = im.convert('P', palette=Image.ADAPTIVE, colors=ncolors)

            im2.save(path, optimize=True)
        else:
            self.fig.savefig(path, transparent=transparent, dpi=dpi)

        # inform user of successful save
        if not quiet:
            print('Successfully saved: %s' % path)

    def update(self, force=False, recalc_bonds=None, redraw=None):
        """
        Calls "draw_" methods of current drawn items

        KArgs:
        - force (bool): if True, all items are forced to redraw
                        (Default: False)

        - recalc_bonds (bool): if True, bonding is recalculated
                               (Default: False)

        - redraw (list): if given, only items in <redraw> will have their
                         "draw_" method called
                         (Default: None -> update everything)
        """
        if isinstance(recalc_bonds, bool):
            self.recalc_bonds = recalc_bonds

        # only redraw items in <redraw> list if given
        if redraw is not None:
            to_redraw = [d for d in self._drawn if d in redraw]
        # else redraw everything
        else:
            to_redraw = self._drawn

        # call all draw methods
        self.draw(to_redraw, force=force)

    def save_rot_gif(self, path=None, fps=20, loop_time=6,
                     max_px=600, rot_axis=None, save_frames=False,
                     overwrite=False, optimize_gif=False):
        """
        Creates a smooth rotating gif of molecule

        KArgs:
        - path (str): path/directory to save gif
                      (Default: use self.name and save to CWD)

        - fps (int): frames per second of gif
                     (Default: 20)

        - loop_time (int): number of seconds for atoms to complete one rotation
                           (Default: 6)

        - max_px (int): sets pixel count for longest dimension in gif
                        (Default: 600)

        - rot_axis (str): specify axis to rotate about
                          - x (left-to-right), y (bot-to-top), or
                            z (counter-clockwise)
                          - can be: 'y' | 'x' | 'z',
                          - can also add a '-' to change direction!
                            e.g. '-z' causes clockwise rotation
                          (Default: refer to self.rot_axis or 'y')

        - save_frames (bool): if True, folder is made - frames saved as pngs
                          - NOTE: won't make a gif
                          (Default: False)

        - overwrite (bool): if False, '-<integer>' is added to end of name
                            to avoid overwriting
                            (Default: False)

        - optimize_gif (bool): if True, creates optimized gif (~1/2 file size)
                               - experimental; still needs additional testing
                               - NOTE: gif takes much longer to make
                               (Default: False)
        """
        # get original atom positions
        orig_pos = self.atoms.positions.copy()

        if not path:
            path = self.name + '.gif'
            # TODO: Implement self.fig_params for
            # more descriptive default gif name
            # path = '%s-%s.gif' % (self.name, '-'.join(self.fig_params))
        elif os.path.isdir(path):
            path = os.path.join(path, self.name + '.gif')

        # ensure path ends with .gif
        if not path.endswith('.gif'):
            path += '.gif'

        # avoid overwriting previous gif
        if not overwrite:
            path = utils.avoid_overwrite(path)

        # calc gif info
        # total number of frames
        frames = int(fps * loop_time)
        # angle to rotate between frames (in degrees)
        rot = 360 / frames
        dig_str = '%0{}i'.format(len(str(frames)))

        # cannot optimize a gif if a colorbar is drawn
        if 'colorbar' in self._drawn and optimize_gif:
            print('Unable to optimize a gif with a colorbar.')
            optimize_gif = False

        # print gif info
        print('-' * 50)
        print('          Title: %s' % (os.path.basename(path)))
        print('      Loop time: %.2f s' % loop_time)
        print('            FPS: %i' % fps)
        print('   Total frames: %i' % frames)
        print(' Building frame: ' + dig_str % 1, end='\r')

        # save frames in *_frames folder
        if save_frames or optimize_gif:
            if save_frames:
                frame_name = os.path.basename(path).strip('.gif') + '_frames'
                frame_path = os.path.join(os.path.dirname(path), frame_name)

                # avoid overwriting previous frame_path
                frame_path = utils.avoid_overwrite_dir(frame_path)
                os.mkdir(frame_path)
            # else use temp directory for optimize_gif
            else:
                frame_path = tempfile.mkdtemp()

            # iteratively save all frames
            print(' ' * 50, end='\r')
            for i in range(frames):
                print('   Saving frame: %03i' % (i + 1), end='\r')
                self.save(os.path.join(frame_path,
                                       self.name + '_%03i.png' % (i + 1)),
                          max_px=max_px, optimize=optimize_gif, quiet=True)
                self.rotate(rot)
        # else make gif with matplotlib.animation and image magick
        else:
            animation = self.build_rot_animation(
                                        fps=fps,
                                        loop_time=loop_time,
                                        rot_axis=rot_axis,
                                        verbose=True)

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

            # delete temp folder and files if only making optimized gif
            if optimize_gif and not save_frames:
                shutil.rmtree(frame_path)

        # set atom positions back to original
        self.atoms.positions = orig_pos

        # let user know the gif has been saved
        print(' ' * 50, end='\r')
        print('          Saved: %s' % path)
        print('-' * 50)

    def build_rot_animation(self, fps=20, loop_time=6, rot_axis=None,
                            verbose=False):
        """
        Builds a matplotlib animation of the molecule rotating

        KArgs:
        - fps (int): frames per second of gif
                     (Default: 20)

        - loop_time (int): number of seconds for atoms to complete one rotation
                           (Default: 6)

        - rot_axis (str): specify axis to rotate about
                          - x (left-to-right), y (bot-to-top), or
                            z (counter-clockwise)
                          - can be: 'y' | 'x' | 'z',
                          - can also add a '-' to change direction!
                            e.g. '-z' causes clockwise rotation
                          (Default: refer to self.rot_axis or 'y')

        - verbose (bool): if True, prints status of current frame being built
                          (Default: False)
        """
        # calc gif info
        frames = int(fps * loop_time)
        rot = 360 / frames
        dig_str = '%0{}i'.format(len(str(frames)))

        # adjust axis limits based on rot_axis
        # setting self.rot_axis automatically redraws
        changed = False
        if rot_axis is None:
            if self._rot_axis is None:
                self.rot_axis = 'y'
        else:
            self.rot_axis = rot_axis

        def next_step(i):
            # print out progress if <verbose>
            if verbose:
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
                                       interval=1000/fps)
        return animation

    def _check_colors(self, value):
        """
        Handles new assignment to color property
        """
        # track to see if colors have changed
        old_colors = copy.deepcopy(self._colors)

        # params to block creation of colorbar and legend
        # - both depend on color assignments given
        self.block_colorbar = True
        self.block_legend = False

        # track if colors were successfully converted
        failed = True

        if str(value).lower() == 'random':
            self._colors = list(np.random.random((len(self), 3)))
            failed = False
            self.block_legend = True
        # resort to default colors
        elif value is None:
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
        # convert dict of color specs
        elif isinstance(value, dict):
            value, allow_legend = self._dict_to_list(value, self._colors)
            self.block_legend = not allow_legend

            for i in range(len(value)):
                if isinstance(value[i], str):
                    value[i] = list(mcolors.to_rgba(value[i]))[:-1]
            self._colors = value
            failed = False
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

                self.cb_norm = mcolors.Normalize(vmin=minval,
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
            except Exception:
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
        if (isinstance(self._colors, type(None)) or
           isinstance(value, type(None))):
            self._colors = [list(jmol_colors[i.number]) for i in self.atoms]

        # all atoms must be accounted for in colors list
        if len(self._colors) != len(self.atoms):
            raise ValueError("length of colors must be equal to # atoms")

        # redraw atoms if colors have changed
        if self._colors != self._prevcolors and 'atoms' in self._drawn:
            self.update(redraw=['atoms', 'legend'])

    def _check_labels(self, value):
        """
        Handles new assignment to labels property
        """
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
                    self._labels = [round(i, 2)
                                    for i in self.atoms.get_initial_charges()]
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

    def _dict_to_list(self, prop_dict, default_ls):
        """
        Converts dict of custom properties into array of values for each atom
        - works for colors and alpha values

        Args:
        prop_dict (dict): dictionary of properties for certain atoms
                          - key = chem symbol: all atoms of that symbol
                          - key = index: just that atom index
        default_ls (list): default values for atoms in molecule

        Returns:
        (np.ndarray): 1D array of custom and default values, ordered by atom
        allow_legend (bool): if True, each atom type has a consistent value
        """
        val_ls = copy.deepcopy(default_ls)

        # track or properties are consistant for each atom type
        allow_legend = True

        # track values for each atom type if consistent
        const_val_dict = {}

        arr = np.array(val_ls)
        # see if items are iterable (i.e. arr is larger than 1D)
        nd = len(arr.shape) > 1

        for s in set(self.atoms.symbols):
            uniq = np.unique(arr[self.atoms.symbols == s], axis=0)
            if len(uniq) > 1:
                allow_legend = False
                break
            const_val_dict[s] = uniq.tolist()[0]

        # iterate over symbols first, then atom indices
        for k in sorted(prop_dict, key=lambda i: [not isinstance(i, str), i]):
            # if chemical symbol or R (equivalent to C and H)
            if isinstance(k, str) and (k.title() in chemical_symbols or
               k.upper() == 'R'):
                # find C and H atoms
                if k.upper() == 'R':
                    inds = np.where(np.isin(self.atoms.symbols, ['C', 'H']))[0]
                else:
                    inds = np.where(self.atoms.symbols == k.title())[0]
                for i in inds:
                    val_ls[i] = prop_dict[k]
            elif ((isinstance(k, int) or str(k).isdigit()) and
                  (-1 < int(k) < len(val_ls))):
                val_ls[int(k)] = prop_dict[k]
                if prop_dict[k] != const_val_dict[self.atoms.symbols[int(k)]]:
                    allow_legend = False

        return val_ls, allow_legend
