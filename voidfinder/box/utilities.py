from scipy.io import FortranFile
import numpy as np

def read_unformatted(filename):
    '''
    Reads an unformatted Fortran 90 files as
    a numpy array.

    Parameters:  filename: str
                 Name of the unformatted file.
    '''

    fin = FortranFile(filename, 'r')
    nrows = fin.read_ints()[0]
    ncols = fin.read_ints()[0]
    data = fin.read_reals(dtype=np.float64).reshape(nrows, ncols)

    return data, nrows, ncols


def save_as_unformatted(data, filename):
    '''
    Saves a numpy array as an unformatted
    Fortran 90 file that can be handled by
    this package's numerical routines.

    Parameters:  data: ND array_like
                 Array to be saved.

                 filename: str
                 Name of the output file.
    '''
    data = np.asarray(data)

    nrows, ncols = np.shape(data)
    f = FortranFile(filename, 'w')
    nrows, ncols = np.shape(data)
    f.write_record(nrows)
    f.write_record(ncols)
    f.write_record(data)
    f.close()
