import numpy as np
from astropy.io import fits
import h5py
import sys

# Only important imput is the file name of the data in the fits format. 
# Here we will presume it looks the same as the GREGOR data

input_filename = sys.argv[1] # this is without the .fits extension, can it be done better?



#Here we read the data from the cube:
dataset = fits.open(input_filename+'.fits')
stokes_cube = dataset[0].data
wavelength = dataset[1].data

# Hard-set boundaries for the inversion limit 
# Should this be specified in some other config file?
# Maybe, but for the moment, invert everything.

NX = stokes_cube.shape[0]
NY = stokes_cube.shape[1]
NS = stokes_cube.shape[2]
NL = stokes_cube.shape[3]

print ("Dimensions of the field :", NX, NY)
print ("Total number of wavelenghts :", NL)

if (NS != 4):
	print("Wrong number of Stokes components. Check your data!")
	print("Should we 4, instead its: ",NS)
	exit();

NL_check = len(wavelength)

if (NL_check != NL):
	print("Data file and wavelenght file don't have the same length! Check you data!")
	exit();

# Here we select only wavelength range from 10828.5 to 10831.3 
#(This can also be in a conf file?)
l_left = np.argmin(np.abs(wavelength-10828.5))
l_right = np.argmin(np.abs(wavelength-10831.3)) + 1

n_wvl = l_right - l_left

print("Indices of wavelenght boundaries are: ",l_left,l_right)

# TODO: What if the wavelength range is smaller than this? Will this still work? Should, right.


# Normalization always with respect to the LOCAL continuum (meaning, in THAT pixel)
loc_continuum = np.amax(stokes_cube[:,:,0,l_left:l_right],axis=2)
stokes_cube /= loc_continuum[:,:,np.newaxis,np.newaxis]
del(loc_continuum)

# EXTRA STEP, correct the data for the limb darkning:
# These two need to be fetched from fits: 
X_solar = 0.0 # in arcsec
Y_solar = 0.0 # in arcsec
R_solar = np.sqrt(X_solar ** 2.0 + Y_solar ** 2.0)
theta = np.arcsin(R_solar / 960.0) # this is very important, also for conf file for hazel

# Take limb darkening function from hazel and put the coefficients here:

ld_correction = 1.0

# =======================================================================================

# Preparing the data is done. Let's write it up: 

# =======================================================================================

# Save the data for the inversion
# First the wavelength axis
np.savetxt(input_filename+'.wavelength', wavelength[l_left:l_right], header='lambda')


# Now the wavelength dependent weights
f = open(input_filename+'.weights', 'w')
f.write('# WeightI WeightQ WeightU WeightV\n')
for i in range(n_wvl):
    f.write('1.0    1.0   1.0   1.0\n')
f.close()

# This should be automatically extracted from the data.
# But, L1 is kinda promising that this will be in the fits!
# For the moment let's see if it is in the fits. If not, we can find it here
# as a standard deviation of the polarization in the continuum
noise = 1E-3

# And the Stokes parameters
x_start = 0
x_end = x_start + NX 
y_start = 0
y_end = y_start + NY
n_pixel = NX * NY

stokes_3d = stokes_cube[x_start:x_end+1,y_start:y_end+1,:,l_left:l_right] * ld_correction
stokes_3d = stokes_3d.reshape(n_pixel,4,n_wvl)
stokes_3d = stokes_3d.transpose(0,2,1)

# This should be changed according to the newest version of hazel
sigma_3d = np.zeros((n_pixel,n_wvl,4), dtype=np.float64)
sigma_3d[:,:,:] = noise

# And this:
los_3d = np.zeros((n_pixel,3), dtype=np.float64)
los_3d[:,0] = theta
los_3d[:,1] = 0.0
los_3d[:,2] = 0.0 # This is referent direction for Stokes Q. Ask DLNISRP / L1 how to find this!

# For example. For FIRS this is arctan(Y_solar / X_solar)

# And this too. Always (1,0,0,0)
boundary_3d = np.zeros((n_pixel,n_wvl,4), dtype=np.float64)
boundary_3d[:,:,0] = 1.0

# Write down the file:
output_filename = input_filename+'.h5'
f = h5py.File(output_filename, 'w')
db_stokes = f.create_dataset('stokes', stokes_3d.shape, dtype=np.float64)
db_sigma = f.create_dataset('sigma', sigma_3d.shape, dtype=np.float64)
db_los = f.create_dataset('LOS', los_3d.shape, dtype=np.float64)
db_boundary = f.create_dataset('boundary', boundary_3d.shape, dtype=np.float64)
db_stokes[:] = stokes_3d
db_sigma[:] = sigma_3d
db_los[:] = los_3d
db_boundary[:] = boundary_3d
f.close()

# Change the configuration file:

# ========================================================================================

config_f_o = open("conf_multi.ini", 'r')
original_config = config_f_o.readlines()

file_keys = ['Output file', 'Wavelength file', 'Wavelength weight file', 'Observations file']
file_keys_pos = np.zeros(len(file_keys))
file_keys_pos = file_keys_pos.astype('int')

keys_pre        = ['', '	', '	', '	']
keys_extensions = ['_l2.h5', '.wavelength', '.weights', '.h5']

#print (file_keys[0])
#print (file_keys[1])
#print (file_keys[2])
#print (file_keys[3])

# Find where each one of these is:

for i in range(len(file_keys)):
    for j in range(len(original_config)):
        if original_config[j].find(file_keys[i]) >=0:
            file_keys_pos[i]=j
            print (file_keys_pos[i])

            original_config[j] = keys_pre[i] + file_keys[i] + ' = ' + input_filename + keys_extensions[i]+ '\n'

config_f_output = open('conf_'+input_filename+'.ini','w')
config_f_output.writelines(original_config)
config_f_output.close()	

# ========================================================================================