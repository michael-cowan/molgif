import click
from molgif.gifmaker import rot_gif
from molgif._version import __version__
from ase.data import chemical_symbols
import matplotlib.cm as cm


@click.command(name='molgif',
               context_settings=dict(help_option_names=['-h', '--help']))
@click.version_option(__version__)
@click.argument('atoms')
@click.option('-p', '--save-path', default=None, type=str, metavar='<s>',
              help='path to save visual - based on extension:'
                   ' .gif = animation, .png or .svg = image'
                   ' [default: chem formula and gif info]')
@click.option('-i', '--img', is_flag=True, help='save png instead of gif')
@click.option('-v', '--vis', is_flag=True,
              help='opens visual of molecule - does not save new files')
@click.option('-r', '--smart-rotate', is_flag=True,
              help='orients atoms such that max variance is in x-axis')
@click.option('-c', '--colors', '--color', default=None, type=str, metavar='<s>',
              help='specify atom colors or colors of specific atom types\n\n'
                   '-c blue = all atoms are blue\n\n'
                   '-c C-blue-H-green = C are blue, H are green\n\n'
                   '[default: jmol colors]')
@click.option('-t', '--loop-time', default=6, type=float, show_default=True,
              metavar='<f>',
              help='number of seconds for atoms to complete one rotation')
@click.option('-f', '--fps', default=20, show_default=True, metavar='<i>',
              help='define frames per second of gif')
@click.option('--scale', default=0.7, show_default=True, metavar='<f>',
              help='scales size of atoms: scale * covalent radii')
@click.option('--no-bonds', is_flag=True,
              help='removes bond from visual')
@click.option('--hide', default=None, metavar='<i|s>-*',
              help='select atom types and/or atom indices to hide\n\n'
                   '--hide Cd-10-11 -> hides all Cd atoms and atoms 10 and 11')
@click.option('--alphas', '--alpha', default=None, metavar='<i|s>-<f>-*',
              help='set transparency of atom types or atom indices\n\n'
                   '--alphas Cd-0.1 -> all Cd atoms are almost see through')
@click.option('--fade', default=None, metavar='<i|s>-*',
              help='fades selected atoms (sets alpha = 0.2 and color = white')
@click.option('--rotate', default=None, type=str, metavar='<i>-<s>-*',
              help='list of ordered rotation commands to apply\n\n'
                   '- "90-x-60-z" => rot 90deg about x THEN 60deg about z\n\n'
                   '- Always occurs AFTER smart_rotate')
@click.option('--rot-axis', default='y', show_default=True, metavar='<s>',
              help='gif rotation axis:'
                   ' x (left-to-right), y (bot-to-top), or z (ccw);'
                   ' can also add a "-" to change direction!')
@click.option('--anchor', default=None, metavar='<i|s>',
              help='define atom to anchor to center such'
                   'that all other atoms rotate around it\n\n'
                   '- index: index of atom to anchor\n\n'
                   '- "center": closest to center of position is anchored\n\n'
                   '- chem-symbol: first atom type found with that symbol'
                   ' (based on index order) is anchored')
@click.option('--max-px', default=600, show_default=True, metavar='<i>',
              help='sets pixel count for longest dimension in visual')
@click.option('--square', is_flag=True,
              help='ensure visual has 1:1 aspect ratio')
@click.option('--legend', 'draw_legend', is_flag=True,
              help='adds an atom type legend to visual')
@click.option('--leg-order', default='size', show_default=True, metavar='<s>',
              help='Used to order the legend: size, size_r, alpha, number')
@click.option('--legend-max-ms', default=20, show_default=True, metavar='<i>',
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
@click.option('--cb-min', default=None, type=float, metavar='<f>',
              help='define min value limit for colorbar')
@click.option('--cb-max', default=None, type=float, metavar='<f>',
              help='define max value limit for colorbar')
@click.option('--cmap', default='bwr_r', show_default=True, metavar='<s>',
              help='matplotlib cmap to be used if values are used for colors')
@click.option('--center-data', is_flag=True,
              help='colors are centered about middle of cmap when using values'
                   ' as colors')
@click.option('--labels', default=None, type=str, metavar='<s>',
              help='define labels to add to atoms\n\n'
                   '- "symbols", "symbols-noh", "values" (from colors option),'
                   ' or "charges" (see --use-charges)\n\n'
                   '[default: No labels]')
@click.option('--label-size', default=None, type=int, metavar='<i>',
              help='set size of labels')
@click.option('--bond-color', default='white', show_default=True,
              metavar='<s>', help='define color of bonds')
@click.option('--bond-edgecolor', default='black', show_default=True,
              metavar='<s>', help='define color of bond edges (outline)')
@click.option('--save-frames', is_flag=True,
              help='folder is made and gif frames saved as pngs')
def cli(atoms, save_path, img, vis, smart_rotate, colors, loop_time, fps,
        scale, no_bonds, hide, alphas, fade, rotate, rot_axis, anchor, max_px,
        square, draw_legend, leg_order, legend_max_ms, optimize, transparent,
        overwrite, use_charges, draw_colorbar, cb_min, cb_max, cmap,
        center_data, labels, label_size, bond_color, bond_edgecolor,
        save_frames):
    """
    molgif: create smooth gifs of rotating molecules

    \b
    - can also visualize (--vis) or quickly create images (--img)

    ATOMS: atoms to be visualized - path to geometry file
    """
    # fade specific atoms (by atom type or indices)
    # sets alpha to 0.2 and color to white
    # does not change alpha or color if user sets it themselves
    if fade is not None:
        fade_ls = fade.lower().split('-')
        # convert R group to C and H atoms
        if 'r' in fade_ls:
            fade_ls.remove('r')
            fade_ls += ['c', 'h']

        for f in fade_ls:
            if colors is None:
                colors = f'{f}-white'
            elif f not in colors.lower().split('-'):
                colors += f'-{f}-white'
            if alphas is None:
                alphas = f'{f}-0.2'
            elif f not in alphas.lower().split('-'):
                alphas += f'-{f}-0.2'

    # if - is used, convert atom type colors to dict
    if colors is not None and '-' in colors:
        c = colors.split('-')

        # if chemical symbol, R, or index given, assume dict
        if any(s.title() in chemical_symbols or s.isdigit() or s.upper() == 'R'
               for s in c[::2]):
            colors = {k if not k.isdigit() else int(k): v
                      for k, v in zip(c[::2], c[1::2])}
        # else assume list of colors
        else:
            colors = c

    # check alpha values
    if alphas is not None:
        if '-' in alphas:
            alphas = alphas.split('-')
            if len(alphas) % 2:
                print('Invalid alpha values given')
            alphas = {i if not i.isdigit() else int(i): float(j)
                      for i, j in zip(alphas[::2], alphas[1::2])}
        elif alphas.replace('.', '', 1).isdigit():
            alphas = float(alphas)

    # make sure rotate option is valid
    if rotate is not None:
        # rotate param must be in the form <int>-<xyz>-<int>-<xyz>-...
        # convert to list and check rotate form
        rotate = rotate.split('-')
        if (len(rotate) % 2 == 1 or
           any(not i.isdigit() for i in rotate[::2]) or
           any(j.lower() not in 'xyz' for j in rotate[1::2])):
            # if invalid, ignore rotate param and continue
            print('invalid rotate argument. ignoring rotation commands.')
            rotate = None

    # if hide argument, ensure that format is correct
    if hide is not None:
        hide = hide.split('-')

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
            hide=hide,
            alphas=alphas,
            custom_rotate=rotate,
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
