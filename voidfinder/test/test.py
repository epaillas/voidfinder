from context import box
import numpy as np
import matplotlib.pyplot as plt

input_filename = 'subsampled_data.dat'
output_filename = 'subsampled_data.unf'

box.ascii_to_unformatted(input_filename, output_filename)

# tpcf
data_filename1 = 'subsampled_data.unf'
output_filename = 'tpcf.dat'
box_size = 1500
dim1_min = 0
dim1_max = 150
dim1_nbin = 50
ngrid = 100

r, corr = box.tpcf(data_filename1=data_filename1, output_filename=output_filename,
box_size=box_size, dim1_min=dim1_min, dim1_max=dim1_max,
dim1_nbin=dim1_nbin, ngrid=ngrid)

fig, ax = plt.subplots()
ax.plot(r, r**2*corr)
plt.show()


# tpcf_2d
data_filename1 = 'subsampled_data.unf'
output_filename = 'tpcf_2d.dat'
box_size = 1500
dim1_min = 0
dim1_max = 150
dim1_nbin = 50
ngrid = 100

r, corr = box.tpcf_2d(data_filename1=data_filename1, output_filename=output_filename,
box_size=box_size, dim1_min=dim1_min, dim1_max=dim1_max,
dim1_nbin=dim1_nbin, ngrid=ngrid)

fig, ax = plt.subplots()
ax.plot(r, corr)
plt.show()
