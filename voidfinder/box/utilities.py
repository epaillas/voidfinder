import sys
import numpy as np
from scipy.io import FortranFile
import argparse
from scipy.interpolate import InterpolatedUnivariateSpline
from scipy.special import eval_legendre
from scipy.integrate import simps

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

def read_array_2d(input_filename):
  data = np.genfromtxt(input_filename)
  dim1 = np.unique(data[:,0])
  dim2 = np.unique(data[:,1])

  vary_dim2 = False
  if data[0,0] == data[1,0]:
      vary_dim2 = True

  result = np.zeros([len(dim1), len(dim2)])
  counter = 0
  if vary_dim2:
      for i in range(len(dim1)):
          for j in range(len(dim2)):
              result[i, j] = data[counter, 2]
              counter += 1
  else:
      for i in range(len(dim2)):
          for j in range(len(dim1)):
              result[j, i] = data[counter, 2]
              counter += 1
  return dim1, dim2, result

def get_multipole(ell, s, mu, xi_smu):
  multipole = np.zeros(xi_smu.shape[0])
  if mu.min() < 0:
      factor = 2
      mumin = -1
  else:
      factor = 1
      mumin=0
  for i in range(xi_smu.shape[0]):
      mufunc = InterpolatedUnivariateSpline(mu, xi_smu[i, :], k=3, ext=3)
      xaxis = np.linspace(mumin, 1, 1000)
      lmu = eval_legendre(ell, xaxis)
      yaxis = mufunc(xaxis) * (2 * ell + 1) / factor * lmu
      multipole[i] = simps(yaxis, xaxis)
  return s, multipole