import click
from molgif.gifmaker import rot_gif
from molgif._version import __version__
import matplotlib.cm as cm


@click.command(name='molgif',
               context_settings=dict(help_option_names=['-h', '--help']))
@click.version_option(__version__)
@click.argument('atoms')
@click.option('-p', '--save-path', default=None, type=str,
              help='path to save visual - based on extension:'
                   ' .gif = animation, .png or .svg = image'
                   ' [Default: chem formula and gif info]')
@click.option('-i', '--img', is_flag=True, help='save png instead of gif')
@click.option('-v', '--vis', is_flag=True,
              help='opens visual of molecule - does not save new files')
@click.option('-r', '--smart-rotate', is_flag=True,
              help='orients atoms such that max variance is in x-axis')
@click.option('-c', '--colors', default=None, type=str,
              help='specify atom colors or colors of specific atom types\n\n'
                   '-c blue = all atoms are blue\n\n'
                   '-c C-blue-H-green = C are blue, H are green\n\n'
                   '[Default: jmol colors]')
@click.option('-t', '--loop-time', default=6, type=float, show_default=True,
              help='number of seconds for atoms to complete one rotation')
@click.option('-f', '--fps', default=20, show_default=True,
              help='define frames per second of gif')
@click.option('--scale', default=0.7, show_default=True,
              help='scales size of atoms: scale * covalent radii')
@click.option('--no-bonds', is_flag=True,
              help='removes bond from visual')
@click.option('--rot-axis', default='y', show_default=True,
              help='gif rotation axis:'
                   ' x (left-to-right), y (bot-to-top), or z (ccw);'
                   ' can also add a "-" to change direction!')
@click.option('--anchor', default=None,
              help='define atom to anchor to center such'
                   'that all other atoms rotate around it\n\n'
                   '- index: index of atom to anchor\n\n'
                   '- "center": closest to center of position is anchored\n\n'
                   '- chem-symbol: first atom type found with that symbol'
                   ' (based on index order) is anchored')
@click.option('--max-px', default=600, show_default=True,
              help='sets pixel count for longest dimension in visual')
@click.option('--square', is_flag=True,
              help='ensure visual has 1:1 aspect ratio')
@click.option('--legend', 'draw_legend', is_flag=True,
              help='adds an atom type legend to visual')
@click.option('--leg-order', default='size', show_default=True,
              help='Used to order the legend: size, size_r, alpha, number')
@click.option('--legend-max-ms', default=20, show_default=True,
              help='scales atoms in legend')
@click.option('-o', '--optimize', is_flag=True,
              help='(experimental) attempts to optimize file size of visual')
@click.option('--transparent', is_flag=True,
              help='images (not gifs) saved with transparent background')
@click.option('--overwrite', is_flag=True,
              help='enables new visual to overwrite previous files')
@click.option('--use-charges', is_flag=True,
              help='atoms are colored by initial_charges in ase atoms object')
@click.option('--colorbar', 'draw_colorbar', is_flag=True,
              help='Draws a colorbar if values are used as colors')
@click.option('--cb-min', default=None, type=float,
              help='define min value limit for colorbar')
@click.option('--cb-max', default=None, type=float,
              help='define max value limit for colorbar')
@click.option('--cmap', default='bwr_r', show_default=True,
              help='matplotlib cmap to be used if values are used for colors')
@click.option('--center-data', is_flag=True,
              help='colors are centered about middle of cmap when using values'
                   ' as colors')
@click.option('--labels', default=None, type=str,
              help='define labels to add to atoms'
                   ' - "symbols", "symbols-noh", "values" (from colors option),'
                   ' or "charges" (see --use-charges)'
                   '[Default: No labels]')
@click.option('--label-size', default=None, help='set size of labels')
@click.option('--bond-color', default='white', show_default=True,
              help='define color of bonds')
@click.option('--bond-edgecolor', default='black', show_default=True,
              help='define color of bond edges (outline)')
@click.option('--save-frames', is_flag=True,
              help='folder is made and gif frames saved as pngs')
def cli(atoms, save_path, img, vis, smart_rotate, colors, loop_time, fps,
        scale, no_bonds, rot_axis, anchor, max_px, square, draw_legend,
        leg_order, legend_max_ms, optimize, transparent, overwrite,
        use_charges, draw_colorbar, cb_min, cb_max, cmap, center_data,
        labels, label_size, bond_color, bond_edgecolor, save_frames):
    """
    molgif: create smooth gifs of rotating molecules

    \b
    - can also visualize (--vis) or quickly create images (--img)

    ATOMS: atoms to be visualized - path to geometry file
    """
    # if - is used, convert atom type colors to dict
    if colors is not None and '-' in colors:
        c = colors.split('-')
        colors = {k: v for k, v in zip(c[::2], c[1::2])}

    # if cmap str is given, ensure that it exists as a cmap
    # if not, default to bwr_r cmap (used when cmap = None)
    if cmap is not None:
        try:
            cmap = getattr(cm, cmap)
        except AttributeError:
            print('cmap given does not exist, using bwr_r')

    rot_gif(atoms=atoms,
            save_path=save_path,
            img=img,
            vis=vis,
            smart_rotate=smart_rotate,
            colors=colors,
            loop_time=loop_time,
            fps=fps,
            scale=scale,
            draw_bonds=not no_bonds,
            rot_axis=rot_axis,
            anchor=anchor,
            max_px=max_px,
            square=square,
            draw_legend=draw_legend,
            leg_order=leg_order,
            legend_max_ms=legend_max_ms,
            optimize=optimize,
            transparent=transparent,
            overwrite=overwrite,
            use_charges=use_charges,
            draw_colorbar=draw_colorbar,
            cb_min=cb_min,
            cb_max=cb_max,
            cmap=cmap,
            center_data=center_data,
            labels=labels,
            label_size=label_size,
            bond_color=bond_color,
            bond_edgecolor=bond_edgecolor,
            save_frames=save_frames)


if __name__ == '__main__':
    cli()
