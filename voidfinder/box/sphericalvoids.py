import subprocess
from os import path
import numpy as np
from . import utilities, tesselation


def grow_spheres(
    tracers_filename, handle, box_size,
    density_threshold, ngrid, rvoid_max=100,
    nthreads=1, void_centres='uniform'
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

                 rvoid_max: float
                 Maxiumum radius a void can have (defaults to 100 Mpc/h).

                 nthreads: int
                 Number of threads to speed up calculations (defaults to 1).

                 void_centres: str
                 How to sample prospective void centres: "uniform" or
                 "delaunay".
                 '''

    # check if files exist
    for fname in [tracers_filename]:
        if not path.isfile(fname):
            raise FileNotFoundError('{} does not exist.'.format(fname))

    # read tracers
    data, *_ = utilities.read_unformatted(tracers_filename)

    # get prospective void centres
    centres_filename = f'{handle}_centres.unf'
    if void_centres == 'delaunay':
        vertices = tesselation.delaunay_triangulation(data, box_size)
        centres = tesselation.get_circumcentres(vertices, box_size)
    elif void_centres == 'uniform':
        ncentres = len(data)
        centres = tesselation.get_random_centres(ncentres, box_size)
    else:
        raise ValueError('void_centres should be either '
                         '"uniform" or "delaunay".')
    utilities.save_as_unformatted(centres, centres_filename)

    # grow spheres around centres
    binpath = path.join(path.dirname(__file__),
                        'bin', 'grow_spheres.exe')

    voids_filename = f'{handle}_voids.dat'

    cmd = [
        binpath,
        tracers_filename,
        centres_filename,
        voids_filename,
        str(box_size),
        str(density_threshold),
        str(rvoid_max),
        str(ngrid),
        str(nthreads)
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
