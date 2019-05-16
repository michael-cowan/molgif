import numpy as np
import ase.neighborlist


def pca(pos, return_transform=False, tranform=None):
    """
    PCA transformation using numpy
    - calculate covariance matrix
    - find eigenvectors (transformation matrix) of
      covariance matrix
    - pos (dot) eigenvectors = tranformed positions

    Args:
    pos (np.ndarray): atom positions

    KArgs:
    return_transform (bool): if True, eigenvectors of covariance matrix is
                             returned with the transformed positions
                             (Default: False)
    transform (np.ndarray): if given, used to transform positions
                            - transforms new positions from different PCA
                            (Default: None)
    """
    # center positions about origin
    pos -= pos.mean(0)

    # if transform supplied, apply it to pos and return
    if isinstance(tranform, np.ndarray):
        return np.dot(pos, tranform)

    # covariance = (x.T * x) / (n - 1)
    # do not need to subtract off mean since data is centered
    co = np.dot(pos.T, pos) / (len(pos) - 1)

    # find eigenvalues and eigenvectors
    evals, evecs = np.linalg.eig(co)

    # sort evals
    indices = evals.argsort()[::-1]
    evals = evals[indices]
    evecs = evecs[:, indices]

    # eigenvectors = transformation matrix
    pos_pca = np.dot(pos, evecs)

    # return evecs if return_transform
    if return_transform:
        return pos_pca, evecs
    else:
        return pos_pca
    return pos_pca, evecs if return_transform else pos_pca


def get_bonds(atoms, radii, scale=1.25):
    """
    Finds bonds between atoms based on bonding radii

    Args:
    atoms (ase.Atoms): atoms object
    radii (np.ndarray): bonding radii of atoms
                        if bonding radii overlap, a bond is drawn

    KArgs:
    scale (float): scales bonding radii array
                   (Default: 1.25)
    """
    n = ase.neighborlist.NeighborList(radii * scale, skin=0,
                                      self_interaction=False)
    n.update(atoms)
    if not n.nneighbors:
        return []

    bonds = np.zeros((n.nneighbors, 2), int)
    spot1 = 0
    for atomi in range(len(atoms)):
        # get neighbors of atomi
        neighs = n.get_neighbors(atomi)[0]

        # find second cutoff in matrix
        spot2 = spot1 + len(neighs)

        # add bonds to matrix
        bonds[spot1:spot2, 0] = atomi
        bonds[spot1:spot2, 1] = neighs

        # shift down matrix
        spot1 = spot2

        # once all bonds have been found break loop
        if spot1 == n.nneighbors:
            break

    return bonds


def draw_bonds(atoms, ax, radii, bond_info):
    """
    Adds bonds to matplotlib axis

    Args:
    atoms (ase.Atoms): atoms object
    ax (matplotib axis): axis object where bonds should be added
    radii (np.ndarray): bonding radii of atoms
                        if bonding radii overlap, a bond is drawn

    KArgs:
    bond_info (tuple): bond_width - width of black outline on bonds
                       bond_fill - width of white fill
                                   (controls outline thickness)
    """
    # bond width (black line) and bond fill (white line)
    bond_width, bond_fill = bond_info

    # find indices of bonds
    bonds = get_bonds(atoms, radii)
    for b in bonds:
        p1 = atoms[b[0]].position
        p2 = atoms[b[1]].position
        x = p1[0], p2[0]
        y = p1[1], p2[1]
        zorder = min(p1[2], p2[2])

        # draw thinner white line over black line to create bordered bonds
        ax.plot(x, y, zorder=zorder - 0.001, color='k', lw=bond_width)
        ax.plot(x, y, zorder=zorder, color='white', lw=bond_fill)


def get_fig_bounds(atoms, rot_axis='y', offset=2):
    # max x and y dists for molecule
    # start by adding in offset
    minx = -offset
    maxx = offset
    miny = -offset
    maxy = offset

    if 'x' in rot_axis:
        # account for z coord in y since atoms rotate about x
        disty = np.linalg.norm(atoms.positions[:, 1:], axis=1).max()
        miny -= disty
        maxy += disty

        # x coords will not change for atoms
        minx += atoms.positions[:, 0].min()
        maxx += atoms.positions[:, 0].max()

    elif 'y' in rot_axis:
        # account for z coord in x since atoms rotate about y
        distx = np.linalg.norm(atoms.positions[:, 0:2:2], axis=1).max()
        minx -= distx
        maxx += distx

        # y coords will not change for atoms
        miny += atoms.positions[:, 1].min()
        maxy += atoms.positions[:, 1].max()

    elif 'z' in rot_axis:
        # account for x and y coords since atoms rotate about z
        distxy = np.linalg.norm(atoms.positions[:, 0:2], axis=1).max()
        minx -= distxy
        maxx += distxy

        miny -= distxy
        maxy += distxy

    # calculate axis limits
    xlim = (minx, maxx)
    ylim = (miny, maxy)

    # calculate width and height based on max atomic positions
    width = maxx - minx
    height = maxy - miny

    # calculate figure size
    # make sure max side length is 5 inches
    if height < width:
        fig_size = (5, 5 * (height / width))
    else:
        fig_size = (5 * (width / height), 5)

    return fig_size, xlim, ylim


def angstrom_to_axunits(val, ax):
    """
    Converts value in Angstrom to units of x axis
    - used to scale bond widths regardless of axis size

    Args:
    val (float): value in Angstroms
    ax (plt.axis): matplotlib axis object used to scale value
    """
    # get x range
    xr = ax.get_xlim()[1] - ax.get_xlim()[0]

    # get width
    fig = ax.get_figure()
    width = fig.bbox_inches.width * ax.get_position().width

    # convert length to points (72 points per inch)
    width *= 72

    # scale value to axis units
    return val * (width / xr)
