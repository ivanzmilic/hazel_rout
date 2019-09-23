import numpy as np
import matplotlib.pyplot as pl
import hazel
import h5py
import sys

filename = sys.argv[1]
x_coordinate = int(sys.argv[2])
y_coordinate = int(sys.argv[3])
x_range = int(sys.argv[4])
y_range = int(sys.argv[5])
l_left = int(sys.argv[6])
l_right = int(sys.argv[7])

#Here we read the data from the cube:
from astropy.io import fits

dataset = fits.open(filename)
stokes = dataset[0].data
ll = dataset[1].data

n_wvl = l_right - l_left
loc_continuum = np.amax(stokes[:,:,0,l_left:l_right],axis=2)

stokes /= loc_continuum[:,:,np.newaxis,np.newaxis]

# Save the data for the inversion
# First the wavelength axis
np.savetxt('10830_gregor.wavelength', ll[l_left:l_right], header='lambda')

# Now the wavelength dependent weights
f = open('10830_gregor.weights', 'w')
f.write('# WeightI WeightQ WeightU WeightV\n')
for i in range(n_wvl):
    f.write('1.0    1.0   1.0   1.0\n')
f.close()

noise = 2E-3
# And the Stokes parameters
n_pixel = x_range*y_range

# Generate mock Stokes parameters full of zeros in this case (so completely useless)
stokes_3d = stokes[x_coordinate:x_coordinate+x_range,y_coordinate:y_coordinate+y_range,:,l_left:l_right]
stokes_3d = stokes_3d.reshape(x_range*y_range,4,n_wvl)
stokes_3d = stokes_3d.transpose(0,2,1)
sigma_3d = np.zeros((n_pixel,n_wvl,4), dtype=np.float64)
sigma_3d[:,:,:] = noise
los_3d = np.zeros((n_pixel,3), dtype=np.float64)
los_3d[:,0] = 0.0
los_3d[:,1] = 0.0
los_3d[:,2] = 90.0
boundary_3d = np.zeros((n_pixel,n_wvl,4), dtype=np.float64)
boundary_3d[:,:,0] = 1.0

f = h5py.File('subset_hi_velocity.h5', 'w')
db_stokes = f.create_dataset('stokes', stokes_3d.shape, dtype=np.float64)
db_sigma = f.create_dataset('sigma', sigma_3d.shape, dtype=np.float64)
db_los = f.create_dataset('LOS', los_3d.shape, dtype=np.float64)
db_boundary = f.create_dataset('boundary', boundary_3d.shape, dtype=np.float64)
db_stokes[:] = stokes_3d
db_sigma[:] = sigma_3d
db_los[:] = los_3d
db_boundary[:] = boundary_3d
f.close()

# ######################
# Inversion:

#iterator = hazel.Iterator(use_mpi=True)
#mod = hazel.Model('conf_multi.ini', rank=iterator.get_rank())
#iterator.use_model(model=mod)
#iterator.run_all_pixels()