import sys
import numpy as np
from scipy.io import FortranFile
import argparse
import glob
from astropy.io import fits

def fits_to_unformatted(
  input_filename, output_filename, cosmology,
  is_random=False, equal_weights=False, zrange=None
):
  # open fits file
  with fits.open(input_filename) as hdul:
    cat = hdul[1].data

  if zrange is not None:
    zmin, zmax = zrange
    ind = (cat['Z'] > zmin) & (cat['Z'] < zmax)
    cat = cat[ind]

  # convert redshifts to distances
  dist = cosmology.ComovingDistance(cat['Z'])
  x = dist * np.cos(cat['DEC'] * np.pi / 180) * np.cos(cat['RA'] * np.pi / 180)
  y = dist * np.cos(cat['DEC'] * np.pi / 180) * np.sin(cat['RA'] * np.pi / 180)
  z = dist * np.sin(cat['DEC'] * np.pi / 180)

  if not equal_weights:
    if is_random:
      weight = cat['WEIGHT_FKP']
    else:
      weight = cat['WEIGHT_FKP'] * cat['WEIGHT_SYSTOT']
  else:
    weight = np.ones(len(cat))

  #write result to output file
  cout = np.c_[x, y, z, weight]
  nrows, ncols = np.shape(cout)
  f = FortranFile(output_filename, 'w')
  nrows, ncols = np.shape(cout)
  f.write_record(nrows)
  f.write_record(ncols)
  f.write_record(cout)
  f.close()

def ascii_to_unformatted(input_filename, output_filename,
  pos_cols=[0, 1, 2], vel_cols=None, weight_cols=None
):
  # import data
  data = np.genfromtxt(input_filename)
  
  pos = data[:, pos_cols]
  cout = pos # default catalogue with only positions

  if vel_cols is not None:
    vel = data[:, vel_cols]
    cout = np.c_[cout, vel]
  if weight_cols is not None:
    weight = data[:, weight_cols]
    cout = np.c_[cout, weight]

  f = FortranFile(output_filename, 'w')
  nrows, ncols = np.shape(cout)
  f.write_record(nrows)
  f.write_record(ncols)
  f.write_record(cout)
  f.close()

def mean_from_mocks(input_handle, output_filename,
  correlation_type='monopole'):
  mock_files = sorted(glob.glob(input_handle))
  data_list = []
  for mock_file in mock_files:
    data = np.genfromtxt(mock_file)
    data_list.append(data)

  data_list = np.asarray(data_list)
  data_mean = np.nanmean(data_list, axis=0)
  if correlation_type == 'monopole':
    data_std = np.nanstd(data_list, axis=0)[:, 1]
  elif correlation_type == 'smu':
    data_std = np.nanstd(data_list, axis=0)[:, 2]
  elif correlation_type == 'sigmapi':
    data_std = np.nanstd(data_list, axis=0)[:, 2]
  else:
    raise RuntimeError('Correlation type not recognized.')

  cout = np.c_[data_mean, data_std]
  cout = np.nan_to_num(cout)
  np.savetxt(output_filename, cout)

  return cout

def covariance_matrix(data, norm=False):
    """
    Assumes rows are observations,
    columns are variables
    """
    nobs, nbins = np.shape(data)
    mean = np.mean(data, axis=0)
    cov = np.zeros([nbins, nbins])

    for k in range(nobs):
        for i in range(nbins):
            for j in range(nbins):
                cov[i, j] += (data[k, i] - mean[i])*(data[k, j] - mean[j])

    cov /= nobs - 1
    
    if norm:
        corr = np.zeros_like(cov)
        for i in range(nbins):
            for j in range(nbins):
                corr[i, j] = cov[i, j] / np.sqrt(cov[i, i] * cov[j, j])
        return corr
    else:
        return cov


