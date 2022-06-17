import numpy as np
from astropy.io import fits
import h5py
import sys

# Only important imput is the file name of the data in the fits format. 
# Here we will presume it looks the same as the GREGOR data
filename = sys.argv[1]



#Here we read the data from the cube:
dataset = fits.open(filename)
stokes_cube = dataset[0].data
wavelength = dataset[1].data

# Hard-set boundaries for the inversion limit (should this be specified in some other config file?)


if (local):
	loc_continuum = np.amax(stokes[:,:,0,l_left:l_right],axis=2)
	stokes /= loc_continuum[:,:,np.newaxis,np.newaxis]

else:
	qs = np.mean(np.amax(stokes[:,:,0,l_left:l_right],axis=2))
	stokes[:,:,:,:] /= qs


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

f = h5py.File(output, 'w')
db_stokes = f.create_dataset('stokes', stokes_3d.shape, dtype=np.float64)
db_sigma = f.create_dataset('sigma', sigma_3d.shape, dtype=np.float64)
db_los = f.create_dataset('LOS', los_3d.shape, dtype=np.float64)
db_boundary = f.create_dataset('boundary', boundary_3d.shape, dtype=np.float64)
db_stokes[:] = stokes_3d
db_sigma[:] = sigma_3d
db_los[:] = los_3d
db_boundary[:] = boundary_3d
f.close()