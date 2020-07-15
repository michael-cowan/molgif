import os
import numpy as np
import ase.neighborlist


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


def avoid_overwrite(path):
    """
    modifies a path so it won't overwrite existing file

    Args:
    path (str): path to modify

    Returns:
    (str): modified (if needed) save path

    Raises:
    ValueError: raised if an extension is not present in path
    """
    # no changes needed if path does not exist
    if not os.path.isfile(path):
        return path

    # if no extension in path, raise ValueError
    if '.' not in path:
        raise ValueError("file path must have an extension")

    # get extension
    ext = path.split('.')[-1]

    # get base name
    base = path[:-len(ext) - 1]

    # incrementally add -(integer).ext until filename doesn't exist
    j = 1
    new_path = base + '-%i.%s' % (j, ext)
    while os.path.isfile(new_path):
        new_path = new_path.replace('-%i.%s' % (j, ext),
                                    '-%i.%s' % (j + 1, ext))
        j += 1
    return new_path


def avoid_overwrite_dir(dirpath):
    """
    Avoid overwriting directories
    """
    if not os.path.isdir(dirpath):
        return dirpath
    else:
        i = 1
        dirpath = dirpath + '-1'
        while os.path.isdir(dirpath):
            dirpath = dirpath.strip('-%i' % i) + '-%i' % (i + 1)
            i += 1
        return dirpath


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
    # remove periodic boundaries
    atoms = atoms.copy()
    atoms.pbc = False

    # create neighborlist object
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


def get_fig_bounds(atoms, rot_axis=None, square=False,
                   padding=1):
    assert padding >= 0
    # max x and y dists for molecule
    # start by adding in padding
    minx = -padding
    maxx = padding
    miny = -padding
    maxy = padding

    # only buffer based on initial perspective
    if rot_axis is None:
        minx += atoms.positions[:, 0].min()
        maxx += atoms.positions[:, 0].max()

        miny += atoms.positions[:, 1].min()
        maxy += atoms.positions[:, 1].max()

    elif 'x' in rot_axis:
        # account for z coord in y since atoms rotate about x
        disty = np.linalg.norm(atoms.positions[:, 1:], axis=1).max()
        miny -= disty
        maxy += disty

        # x coords will not change for atoms
        minx += atoms.positions[:, 0].min()
        maxx += atoms.positions[:, 0].max()

    elif 'y' in rot_axis:
        # account for z coord in x since atoms rotate about y
        distx = np.linalg.norm(atoms.positions[:, 0:3:2], axis=1).max()
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
    if square:
        fig_size = (5, 5)
    elif height < width:
        fig_size = (5, 5 * (height / width))
    else:
        fig_size = (5 * (width / height), 5)

    return fig_size, xlim, ylim


def _opt_angle(atom, tol=1e-6, verbose=False):
    """
    Finds angle to rotate atoms about z-axis to:
    - maximize x and y distance of atoms
    - prioritizes x distance over y distance

    When used with smart_rotate, molecule is "squared off",
    giving a better overall presentation for some systems.
    """
    def get_dists(atom):
        return (atom.positions.max(0) - atom.positions.min(0))[:2]

    atom = atom.copy()

    # angle increment
    inc = 2

    prevang = inc
    ang = inc
    test_atom = atom.copy()
    converged = False
    prevscore = np.product(get_dists(atom))
    count = 0
    while 1:
        # rotate test_atom by current increment
        test_atom.rotate(inc, 'z')

        # calculate current score (product of x and y dists)
        score = np.product(get_dists(test_atom))

        # if score has converged (abs(change) is less than <tol>),
        # return average of previous 2 angles
        if abs(score - prevscore) < tol:
            final = round((ang + prevang) / 2, 4)
            final2 = min([final + 90, final - 90], key=abs)

            # calculate distance with solved angle
            a = atom.copy()
            a.rotate(final, 'z')
            d = get_dists(a)

            # calculate distance with solved angle + 90 degree rotation
            a2 = atom.copy()
            a2.rotate(final2, 'z')
            d2 = get_dists(a2)

            # return angle that maximizes X
            angle = final if d[0] > d[1] else final2

            # print number of iterations needed for optimization
            # at given tolerance (<tol>) and opt angle
            if verbose:
                print('TOLERANCE: %.2e' % tol)
                print(' OPT TOOK: %i iterations!' % count)
                print('OPT ANGLE: %.2f deg.' % angle)

            return angle

        # if rotated too far, reverse and half increment
        if (score > prevscore).all():
            inc *= (-1/2.)

        # print(ang)
        prevang = ang
        ang += inc
        prevscore = score
        count += 1


def path2atoms(path):
    """
    Checks to see if path can be read in as an ase.Atoms object
    - returns path back if it is already an ase.Atoms object

    Args:
    - path (str): path to geometry file

    Returns:
    - ase.Atoms: atoms object
    """
    if isinstance(path, ase.Atoms):
        return path
    try:
        atoms = ase.io.read(path)
        return atoms
    except:
        raise ValueError("path does not lead to a "
                         "supported geometry file")


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
    # can only apply pca to 2 or more atoms
    if pos.shape[0] == 1:
        if return_transform:
            return pos, np.eye(3)

        return pos

    # center positions about origin
    pos -= pos.mean(0)

    # if transform supplied, apply it to pos and return
    if isinstance(tranform, np.ndarray):
        return np.dot(pos, tranform)

    # covariance = (x.T * x) / (n - 1)
    # do not need to subtract off mean since data is centered
    # co = np.dot(pos.T, pos) / (len(pos) - 1)
    co = np.cov(pos, rowvar=False)

    # find eigenvalues and eigenvectors
    evals, evecs = np.linalg.eigh(co)

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


def smart_rotate_atoms(atoms, opt_angle=False):
    """
    Applies "smart" rotation to atoms object

    Args:
    atoms (ase.Atoms): atoms object to rotate

    KArgs:
    - opt_angle: Finds angle to rotate atoms about z-axis to maximize
                 x and y distance of atoms
                 - molecule is "squared off", which may give a better
                   overall presentation for some systems

    Returns:
    (ase.Atoms): new atoms object with transformed (rotated) coords
    """
    # can only 'smart_rotate' 2 or more atoms
    if len(atoms) < 2:
        return atoms

    new_atoms = atoms.copy()
    new_atoms.positions = pca(atoms.positions.copy())

    # use _opt_angle to fine-tune xy-plane rotation
    if opt_angle:
        rotz = _opt_angle(new_atoms)
        new_atoms.rotate(rotz, 'z')

    return new_atoms
