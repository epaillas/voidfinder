from scipy.spatial import Delaunay
import numpy as np


def circumcentres(vertices,  box_size, radius_limit=1000):
    '''
    Find the centre of the circumspheres
    associated to an input catalogue of
    tetrahedra.
    '''
    cenx, ceny, cenz, r = [], [], [], []

    sing = 0
    for tetra in vertices:
        x0, x1, x2, x3 = tetra
        A = []
        B = []
        A.append((x1 - x0).T)
        A.append((x2 - x0).T)
        A.append((x3 - x0).T)
        A = np.asarray(A)
        B = np.sum(A**2, axis=1)
        B = np.asarray(B)
        try:
            C = np.linalg.inv(A).dot(B)
        except Exception:
            sing += 1
            continue
        centre = x0 + 0.5 * C
        radius = 0.5 * np.sqrt(np.sum(C**2))
        if radius < radius_limit:
            cenx.append(centre[0])
            ceny.append(centre[1])
            cenz.append(centre[2])
            r.append(radius)

    cenx = np.asarray(cenx)
    ceny = np.asarray(ceny)
    cenz = np.asarray(cenz)

    # keep only those centres inside the box
    in_box = ((cenx >= 0) & (cenx <= box_size) &
              (ceny >= 0) & (ceny <= box_size) &
              (cenz >= 0) & (cenz <= box_size)
              )
    cenx = cenx[in_box]
    ceny = ceny[in_box]
    cenz = cenz[in_box]

    print(len(cenx))

    cenx = cenx.reshape(len(cenx), 1)
    ceny = ceny.reshape(len(ceny), 1)
    cenz = cenz.reshape(len(cenz), 1)

    centres = np.hstack([cenx, ceny, cenz])

    return centres


def delaunay_triangulation(
    data, box_size
):
    '''
    Make a Delaunay triangulation over
    the cartesian positions of the tracers.
    Returns the vertices of tetrahedra.
    '''

    # add periodic images
    images = get_periodic_images(data, box_size)
    data = np.vstack([data, images])

    triangulation = Delaunay(data)
    simplices = triangulation.simplices.copy()
    vertices = data[simplices]
    return vertices


def periodic_images(data, box_size):
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


def random_seeds(nseeds, box_size):
    '''
    Generates random centres within the simulation
    volume following a uniform distribution.

    Parameters:  nseeds: int
                 Number of random centres.

                 box_size: float
                 Size of the simulation box
    '''

    xpos = np.random.uniform(0, box_size, nseeds)
    ypos = np.random.uniform(0, box_size, nseeds)
    zpos = np.random.uniform(0, box_size, nseeds)

    centres = np.c_[xpos, ypos, zpos]

    return centres
