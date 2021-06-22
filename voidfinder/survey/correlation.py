import sys
import subprocess
from os import path
import numpy as np

def tpcf(
  data_filename, random_filename, output_filename,
  dim1_min, dim1_max, dim1_nbin,
  ngrid, gridmin, gridmax,
  estimator='DP', nthreads=1
):
  '''
  Two-point cross-correlation function between
  two sets of points. If only one sample of points
  is provided, then it is equivalent to the two-point
  autocorrelation function.

  Input arguments:

  data_filename: name of text file containing first
  set of tracers, where the first three columns correspond
  to the positions in cartesian coordinates.

  random_filename: name of file containing 

  output_filename: name of text file where the output 
  is going to be stored.

  box_size: size of the simulation box.

  dim1_min: minimum pair distance for the correlation function.

  dim1_max: maximum pair distance for the correlation function.

  dim1_nbin: number of radial bins for the correlation function.

  ngrid: number of cells in which to divide the simulation box 
  to speed up calculation. It has to be a divisor of box_size 
  (e.g. 100 for 1000 Mpc/h box, or 128 for a 1024 Mpc/h box).

  estimator: the estimator that is used to calculate the 
  two-point correlaton function: 'DP' for Davis & Peebles (1983),
  'LS' for Landy & Szalay (1993).

  nthreads (optional): number of threads to use while running the 
  algorithm.
  '''

  # check if files exist
  if not path.isfile(data_filename):
    raise FileNotFoundError('{} does not exist.'.format(data_filename))

  binpath = path.join(path.dirname(__file__),
  'bin', 'tpcf.exe')

  cmd = [
    binpath, data_filename, random_filename,
    output_filename, str(dim1_min), str(dim1_max),
    str(dim1_nbin), str(ngrid), str(gridmin),
    str(gridmax), estimator, str(nthreads)
  ]

  log_filename = '{}.log'.format(output_filename)
  log = open(log_filename, 'w+')
  subprocess.call(cmd, stdout=log, stderr=log)

  # open output file
  data = np.genfromtxt(output_filename)
  r = data[:,0]
  corr = data[:,1]

  return r, corr

def tpccf(
  data_filename1, data_filename2, random_filename2,
  output_filename, dim1_min, dim1_max,
  dim1_nbin, ngrid, gridmin,
  gridmax, estimator='DP', nthreads=1,
  random_filename1=None, correlation_type='monopole',
  dim2_min=-1, dim2_max=1, dim2_nbin=100
):
  '''
  Two-point cross-correlation function between
  two sets of points. If only one sample of points
  is provided, then it is equivalent to the two-point
  autocorrelation function.

  Input arguments:

  data_filename: name of text file containing first
  set of tracers, where the first three columns correspond
  to the positions in cartesian coordinates.

  random_filename: name of file containing 

  output_filename: name of text file where the output 
  is going to be stored.

  box_size: size of the simulation box.

  dim1_min: minimum pair distance for the correlation function.

  dim1_max: maximum pair distance for the correlation function.

  dim1_nbin: number of radial bins for the correlation function.

  ngrid: number of cells in which to divide the simulation box 
  to speed up calculation. It has to be a divisor of box_size 
  (e.g. 100 for 1000 Mpc/h box, or 128 for a 1024 Mpc/h box).

  estimator: the estimator that is used to calculate the 
  two-point correlaton function: 'DP' for Davis & Peebles (1983),
  'LS' for Landy & Szalay (1993).

  nthreads (optional): number of threads to use while running the 
  algorithm.

  dim1_min: minimum value along the second bin dimension (mu or pi)

  dim2_max: maximum value along the second bin dimension (mu or pi)

  dim2_nbin: number of bins along the second dimension (mu or pi)
  '''

  # check if files exist
  if not path.isfile(data_filename1):
    raise FileNotFoundError(f'{data_filename1} does not exist.')

  if not path.isfile(data_filename2):
    raise FileNotFoundError(f'{data_filename2} does not exist.')

  if not path.isfile(random_filename2):
    raise FileNotFoundError(f'{random_filename2} does not exist.')

  if estimator == 'LS' and random_filename1 == None:
    raise RuntimeError('Lady-Szalay estimator requires a random catalogue for dataset 1.')

  if random_filename1 == None:
    random_filename1 = random_filename2
  

  binpath = path.join(path.dirname(__file__),
  'bin', f'tpccf_{correlation_type}.exe')

  print(binpath)

  if correlation_type == 'monopole':
    cmd = [
      binpath, data_filename1, data_filename2,
      random_filename1, random_filename2, output_filename,
      str(dim1_min), str(dim1_max), str(dim1_nbin),
      str(ngrid), str(gridmin), str(gridmax),
      estimator, str(nthreads)
    ]
  elif correlation_type == 'rmu':
    cmd = [
      binpath, data_filename1, data_filename2,
      random_filename1, random_filename2, output_filename,
      str(dim1_min), str(dim1_max), str(dim1_nbin),
      str(dim2_min), str(dim2_max), str(dim2_nbin),
      str(ngrid), str(gridmin), str(gridmax),
      estimator, str(nthreads)
    ]
  elif correlation_type == 'sigmapi':
    raise ValueError('sigma-pi correlation function is still not supported.')
  else:
    raise ValueError('Invalid correlation type.')

  log_filename = '{}.log'.format(output_filename)
  log = open(log_filename, 'w+')
  subprocess.call(cmd, stdout=log, stderr=log)

  # open output file
  data = np.genfromtxt(output_filename)
  r = data[:,0]
  if correlation_type == 'monopole':
    r = data[:,0]
    corr = data[:,1]
    return r, corr
  if correlation_type == 'rmu':
    r = data[:,0]
    mu = data[:,1]
    corr = data[:,2]
    return r, mu, corr

