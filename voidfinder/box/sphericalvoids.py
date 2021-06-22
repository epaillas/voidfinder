import sys
import subprocess
from os import path
import numpy as np
from scipy.spatial import Delaunay

def delaunay_triangulation(
    tracers_filename, box_size
):
    '''
    Make a Delaunay triangulation over
    the cartesian positions of the tracers.
    Returns the vertices of tetrahedra.
    '''

    # add periodic images
    images = get_periodic_images(data)
    data = np.vstack([data, images])

    triangulation = Delaunay(data)
    simplices = triangulation.simplices.copy()
    vertices = data[simplices]
    print('{} vertices found.'.format(len(vertices)))
    return vertices


def grow_spheres(
    centres_filename, tracers_filename, output_filename,
    box_size, density_threshold, ngrid, rvoid_max=100
    nthreads=1
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
                 '''

    # check if files exist
    for fname in [centres_filename, tracers_filename]:
        if not path.isfile(fname):
            raise FileNotFoundError('{} does not exist.'.format(fname))
    
    binpath = path.join(path.dirname(__file__),
    'bin', 'grow_spheres.exe')

    cmd = [
        tracers_filename,
        centres_filename,
        output_filename,
        str(box_size),
        str(density_threshold),
        str(rvoid_max),
        str(nthreads)
    ]

    log_filename = '{}.log'.format(output_filename)
    log = open(log_filename, 'w+')

    return_code = subprocess.call(cmd, stdout=log, stderr=log)

    if return_code:
        raise RuntimeError('grow_spheres.exe failed. "
                           "Check log file for further information.')
        
    voids = np.genfromtxt(output_filename)

    return voids


def get_periodic_images(data, box_size):
    '''
    Find the relevant images of a
    set of points in a box that
    has boundary conditions.
    '''
    images = []
    buffer = box_size / 10  # probably an overkill

    for point in data:
        condx = ((box_size - buffer) <
                 point[0]) or (point[0] < buffer)
        condy = ((box_size - buffer) <
                 point[1]) or (point[1] < buffer)
        condz = ((box_size - buffer) <
                 point[2]) or (point[2] < buffer)

        if condx and condy and condz:
            shiftx = point[0] + \
                np.copysign(box_size, (buffer - point[0]))
            shifty = point[1] + \
                np.copysign(box_size, (buffer - point[1]))
            shiftz = point[2] + \
                np.copysign(box_size, (buffer - point[2]))
            images.append([shiftx, shifty, shiftz])
        if condx and condy:
            shiftx = point[0] + \
                np.copysign(box_size, (buffer - point[0]))
            shifty = point[1] + \
                np.copysign(box_size, (buffer - point[1]))
            shiftz = point[2]
            images.append([shiftx, shifty, shiftz])
        if condx and condz:
            shiftx = point[0] + \
                np.copysign(box_size, (buffer - point[0]))
            shifty = point[1]
            shiftz = point[2] + \
                np.copysign(box_size, (buffer - point[2]))
            images.append([shiftx, shifty, shiftz])
        if condy and condz:
            shiftx = point[0]
            shifty = point[1] + \
                np.copysign(box_size, (buffer - point[1]))
            shiftz = point[2] + \
                np.copysign(box_size, (buffer - point[2]))
            images.append([shiftx, shifty, shiftz])
        if condx:
            shiftx = point[0] + \
                np.copysign(box_size, (buffer - point[0]))
            shifty = point[1]
            shiftz = point[2]
            images.append([shiftx, shifty, shiftz])
        if condy:
            shiftx = point[0]
            shifty = point[1] + \
                np.copysign(box_size, (buffer - point[1]))
            shiftz = point[2]
            images.append([shiftx, shifty, shiftz])
        if condz:
            shiftx = point[0]
            shifty = point[1]
            shiftz = point[2] + \
                np.copysign(box_size, (buffer - point[2]))
            images.append([shiftx, shifty, shiftz])

    images = np.asarray(images)
    return images
