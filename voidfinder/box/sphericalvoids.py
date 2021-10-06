import subprocess
from os import path
import numpy as np
import yaml
from . import utilities, tesselation


class SphericalVoidFinder:
    """
    Class to identify spherical voids in a cosmological
    volume, as in https://arxiv.org/abs/1810.02864
    """
    def __init__(self, parameter_file: str):
        """
        Initialize class.
        Args:
            parameter_file: YAML file containing configuration parameters.
        """
        with open(parameter_file) as file:
            self.params = yaml.full_load(file)

    def full_pipeline(self):
        """
        Runs the full void identification pipeline in the
        fiducial order.
        """
        self.grow_spheres()
        self.recentre_spheres()
        self.sort_spheres()
        voids = self.filter_spheres()

        return voids

    def grow_spheres(self):
        """
        First step of the spherical void finder. Grows spheres
        that satisfy a certain density_threshold from a set of
        input tracers.
        """

        # check if files exist
        for fname in [self.params["tracers_filename"]]:
            if not path.isfile(fname):
                raise FileNotFoundError('{} does not exist.'.format(fname))

        # read tracers
        try:
            data, *_ = utilities.read_unformatted(
                self.params["tracers_filename"]
            )
        except Exception:
            data = np.genfromtxt(self.params["tracers_filename"])

        # get prospective void seeds
        seeds_filename = f'{self.params["handle"]}_seeds.unf'
        if self.params["seeds_method"] == 'delaunay':
            vertices = tesselation.delaunay_triangulation(
                data, self.params["box_size"]
            )
            seeds = tesselation.circumcentres(
                vertices, self.params["box_size"]
            )
        elif self.params["seeds_method"] == 'uniform':
            if self.params["nseeds"] is None:
                self.params["nseeds"] = len(data)
            seeds = tesselation.random_seeds(
                self.params["nseeds"], self.params["box_size"]
            )
        else:
            raise ValueError('self.params["seeds_method"] should'
                             'be either "uniform" or "delaunay".')
        utilities.save_as_unformatted(seeds, seeds_filename)

        # check file format
        if self.params["tracers_fileformat"] not in ['ascii', 'unformatted']:
            raise ValueError('File format has to be either '
                             '"ascii" or "unformatted".')

        # figure out size of linked list
        if self.params["ngrid"] is None:
            self.params["ngrid"] = utilities.get_closest_divisor(
                self.params["box_size"], 100
            )
        else:
            if type(self.params["ngrid"]) != int:
                raise ValueError('self.params["ngrid"] needs to'
                                 'be an integer.')
            if self.params["box_size"] % self.params["ngrid"] != 0:
                raise ValueError('self.params["ngrid"] needs to'
                                 'be a divisor of self.params["box_size"].')

        # grow spheres around centres
        binpath = path.join(path.dirname(__file__),
                            'bin', 'grow_spheres.exe')

        self.params["uncentred_voids_filename"] = \
            f'{self.params["handle"]}_uncentred_voids.dat'

        if self.params["use_weights"] is True:
            self.params["use_weights"] = 1
        else:
            self.params["use_weights"] = 0

        cmd = [
            binpath,
            self.params["tracers_filename"],
            seeds_filename,
            self.params["uncentred_voids_filename"],
            str(self.params["box_size"]),
            str(self.params["density_threshold"]),
            str(self.params["rvoid_max"]),
            str(self.params["ngrid"]),
            str(self.params["nthreads"]),
            str(self.params["use_weights"]),
            self.params["tracers_fileformat"]
        ]

        log_filename = f'{self.params["handle"]}_growspheres.log'\
                       .format(self.params["handle"])
        log = open(log_filename, 'w+')

        return_code = subprocess.call(cmd, stdout=log, stderr=log)

        if return_code:
            raise RuntimeError('grow_spheres.exe failed. '
                               'Check log file for further information.')

        voids = np.genfromtxt(self.params["uncentred_voids_filename"])

        return voids

    def recentre_spheres(self):
        """
        Second step of the spherical void finder. Recentres the
        voids by shifting their original positions and verifying
        if they can grow more.
        """

        # check if files exist
        for fname in [self.params["tracers_filename"],
            self.params["uncentred_voids_filename"]]:
            if not path.isfile(fname):
                raise FileNotFoundError('{} does not exist.'.format(fname))

        # check file format
        if self.params["tracers_fileformat"] not in ['ascii',
            'unformatted']:
            raise ValueError('File format must be either '
                             '"ascii" or "unformatted".')

        # figure out size of linked list
        if self.params["ngrid"] is None:
            self.params["ngrid"] = utilities.get_closest_divisor(
                self.params["box_size"], 100
            )
        else:
            if type(self.params["ngrid"]) != int:
                raise ValueError('self.params["ngrid"] needs to be '
                                 'an integer.')
            if self.params["box_size"] % self.params["ngrid"] != 0:
                raise ValueError('self.params["ngrid"] needs to be '
                                 'a divisor of self.params["box_size"].')

        # recentre spheres
        binpath = path.join(path.dirname(__file__),
                            'bin', 'recentring.exe')

        self.params["centred_voids_filename"] = \
            f'{self.params["handle"]}_recentred_voids.dat'

        if self.params["use_weights"] is True:
            self.params["use_weights"] = 1
        else:
            self.params["use_weights"] = 0

        cmd = [
            binpath,
            self.params["tracers_filename"],
            self.params["uncentred_voids_filename"],
            self.params["centred_voids_filename"],
            str(self.params["box_size"]),
            str(self.params["density_threshold"]),
            str(self.params["rvoid_max"]),
            str(self.params["ngrid"]),
            str(self.params["nthreads"]),
            str(self.params["use_weights"]),
            self.params["tracers_fileformat"],
            str(self.params["nshifts"])
        ]

        log_filename = f'{self.params["handle"]}_recentring.log'
        log = open(log_filename, 'w+')

        return_code = subprocess.call(cmd, stdout=log, stderr=log)

        if return_code:
            raise RuntimeError('recentring.exe failed. '
                               'Check log file for further information.')

        voids = np.genfromtxt(self.params["centred_voids_filename"])

        return voids

    def sort_spheres(self, voids_filename=None):
        """
        Sort an input void catalogue in
        decreasing order of radius.
        """

        if voids_filename is None:
            voids_filename = self.params["centred_voids_filename"]

        voids = np.genfromtxt(voids_filename)
        voids = voids[np.argsort(voids[:, 3])]
        voids = voids[::-1]

        self.params["sorted_voids_filename"] = voids_filename

        fmt = 4*'%10.3f ' + '%10i ' + '%10.3f '
        np.savetxt(self.params["sorted_voids_filename"], voids, fmt=fmt)

        return voids

    def filter_spheres(self, voids_filename=None):
        """
        Removes overlapping spheres from a void catalogue.
        If the catalogue is not sorted by decreasing void
        radius, the sort_spheres function is applied first.
        """

        if voids_filename is None:
            voids_filename = self.params["sorted_voids_filename"]

        binpath = path.join(path.dirname(__file__),
                            'bin', 'overlapping.exe')

        output_filename = f'{self.params["handle"]}_voids_'\
                          f'overlap{self.params["overlap_factor"]}.dat'

        cmd = [
            binpath,
            self.params["sorted_voids_filename"],
            output_filename,
            str(self.params["box_size"]),
            str(self.params["overlap_factor"]),
            str(self.params["ngrid"])]

        log_filename = f'{self.params["handle"]}_overlapping.log'\
                       .format(self.params["handle"])
        log = open(log_filename, 'w+')

        return_code = subprocess.call(
            cmd, stdout=log, stderr=log
        )

        if return_code:
            raise RuntimeError('overlapping.exe failed. Check log '
                               'file for further information.')

        voids = np.genfromtxt(output_filename)

        return voids
