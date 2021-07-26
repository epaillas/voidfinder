import subprocess
from os import path
import numpy as np
from . import utilities, tesselation


def grow_spheres(
        tracers_filename, handle, box_size,
        density_threshold, ngrid=None, rvoid_max=100,
        nthreads=1, void_centres='uniform', ncentres=None,
        use_weights=False, tracers_fileformat='ascii'
):
    '''
    First step of the spherical void finder. Grows spheres
    that satisfy a certain density_threshold from a set of
    input centres.

    Parameters:  centres_filename: str
                 Name of file unformatted file containing input centres.

                 tracers_filename: str
                 Name of unformatted file file containing tracers.

                 output_filename: str
                 Name of output file that will contain the void catalogue.

                 box_size: float
                 Size of the simulation box.

                 density_threshold: float
                 Density threshold, in units of the mean density (e.g. 0.2).

                 ngrid: int
                 Number of cells along 1D for the linked lists. If None,
                 the divisor of box_size that is closest to 100 is used.

                 rvoid_max: float
                 Maxiumum radius a void can have (defaults to 100 Mpc/h).

                 nthreads: int
                 Number of threads to speed up calculations (defaults to 1).

                 void_centres: str
                 How to sample prospective void centres: "uniform" or
                 "delaunay".

                 use_weights: boolean
                 Whether to use weights during pair counting. Defaults to
                 False.

                 tracers_fileformat: str
                 'ascii' for a text csv file, or 'unformatted' for a binary
                 Fortran 90 file.
                 '''

    # check if files exist
    for fname in [tracers_filename]:
        if not path.isfile(fname):
            raise FileNotFoundError('{} does not exist.'.format(fname))

    # read tracers
    try:
        data, *_ = utilities.read_unformatted(tracers_filename)
    except Exception:
        data = np.genfromtxt(tracers_filename)

    # get prospective void centres
    centres_filename = f'{handle}_centres.unf'
    if void_centres == 'delaunay':
        vertices = tesselation.delaunay_triangulation(data, box_size)
        centres = tesselation.get_circumcentres(vertices, box_size)
    elif void_centres == 'uniform':
        if ncentres is None:
            ncentres = len(data)
        centres = tesselation.get_random_centres(ncentres, box_size)
    else:
        raise ValueError('void_centres should be either '
                         '"uniform" or "delaunay".')
    utilities.save_as_unformatted(centres, centres_filename)

    # figure out file format
    if tracers_fileformat not in ['ascii', 'unformatted']:
        raise ValueError('File format has to be either '
                         '"ascii" or "unformatted".')

    # figure out size of linked list
    if ngrid is None:
        ngrid = utilities.get_closest_divisor(box_size, 100)
    else:
        if type(ngrid) != int:
            raise ValueError('ngrid needs to be an integer.')
        if box_size % ngrid != 0:
            raise ValueError('ngrid needs to be a divisor of box_size.')

    # grow spheres around centres
    binpath = path.join(path.dirname(__file__),
                        'bin', 'grow_spheres.exe')

    voids_filename = f'{handle}_voids.dat'

    if use_weights is True:
        use_weights = 1
    else:
        use_weights = 0

    cmd = [
        binpath,
        tracers_filename,
        centres_filename,
        voids_filename,
        str(box_size),
        str(density_threshold),
        str(rvoid_max),
        str(ngrid),
        str(nthreads),
        str(use_weights),
        tracers_fileformat
    ]

    log_filename = f'{handle}_growspheres.log'.format(handle)
    log = open(log_filename, 'w+')

    return_code = subprocess.call(cmd, stdout=log, stderr=log)

    if return_code:
        raise RuntimeError('grow_spheres.exe failed. '
                           'Check log file for further information.')

    voids = np.genfromtxt(voids_filename)

    return voids


def sort_spheres(voids_filename, radius_col=3):
    '''
    Sort an input void catalogue in
    decreasing order of radius.

    Parameters:  voids_filename: str
                 Name of the void catalogue file.

                 radius_col: int
                 Index of the column containing the void
                 radius.
    '''

    voids = np.genfromtxt(voids_filename)
    voids = voids[np.argsort(voids[:, radius_col])]
    voids = voids[::-1]

    fmt = 4*'%10.3f ' + '%10i ' + '%10.3f '
    np.savetxt(voids_filename, voids, fmt=fmt)

    return voids


def overlapping_filter(
        voids_filename,
        handle,
        box_size, ngrid,
        overlap_factor=0.5
):
    '''
    Removes overlapping spheres from a void catalogue.
    If the catalogue is not sorted by decreasing void
    radius, the sort_spheres function is applied first.

    Parameters:  voids_filename: str
                 Name of the void catalogue file.

                 overlap_factor: float
                 The overlap fraction that should be allowed.
                 Defaults to 0.5.
    '''

    binpath = path.join(path.dirname(__file__),
                        'bin', 'overlapping.exe')

    output_filename = f'{handle}_voids_overlap{overlap_factor}.dat'

    cmd = [
        binpath,
        voids_filename,
        output_filename,
        str(box_size),
        str(overlap_factor),
        str(ngrid)]

    log_filename = f'{handle}_overlapping.log'.format(handle)
    log = open(log_filename, 'w+')

    return_code = subprocess.call(cmd, stdout=log, stderr=log)

    if return_code:
        raise RuntimeError('overlapping.exe failed. '
                           'Check log file for further information.')

    voids = np.genfromtxt(output_filename)

    return voids
