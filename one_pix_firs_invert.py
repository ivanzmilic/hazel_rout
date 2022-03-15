import numpy as np
import matplotlib.pyplot as pl
import hazel
import h5py
import sys

filename = sys.argv[1]
wvlname = sys.argv[2]
x_coordinate = int(sys.argv[3])
y_coordinate = int(sys.argv[4])
l_left = int(sys.argv[5])
l_right = int(sys.argv[6])
local = int(sys.argv[7])

#Here we read the data from the cube:
from astropy.io import fits

dataset = fits.open(filename)
stokes = dataset[0].data
ll = np.loadtxt(wvlname,unpack=True,skiprows=1)
n_wvl = l_right - l_left

if (local):
	norm = np.amax(stokes[x_coordinate,y_coordinate,0,l_left:l_right])
else:
	norm = np.mean(np.amax(stokes[:,:,0,l_left:l_right],axis=(2)),axis=(0,1))

print (norm)


# Save the data for the inversion
# First the wavelength axis
np.savetxt('10830_firs.wavelength', ll[l_left:l_right], header='lambda')

# Now the wavelength dependent weights
f = open('10830_firs.weights', 'w')
f.write('# WeightI WeightQ WeightU WeightV\n')
for i in range(n_wvl):
    f.write('1.0    1.0   1.0   1.0\n')
f.close()

noise = 1E-2
# And the Stokes parameters
f = open('10830_firs_stokes.1d', 'wb')
f.write(b'# LOS theta_LOS, phi_LOS, gamma_LOS\n')
f.write(b'24.64 0.0 347.0\n')
f.write(b'\n')
f.write(b'# Boundary condition I/Ic(mu=1), Q/Ic(mu=1), U/Ic(mu=1), V/Ic(mu=1)\n')
f.write(b'1.0 0.0 0.0 0.0\n')
f.write(b'\n')
f.write(b'# SI SQ SU SV sigmaI sigmaQ sigmaU sigmaV\n')
tmp = np.vstack([stokes[x_coordinate,y_coordinate,:,l_left:l_right]/norm, noise*np.ones((4,n_wvl))])
np.savetxt(f, tmp.T)
f.close()


# ######################
# # And now we do the inversion using the appropriate configuration file
#mod = hazel.Model('conf_firs.ini', working_mode='inversion', verbose=3, randomization=5)
#mod.read_observation()
#mod.open_output()
#mod.invert()
#mod.write_output()
#mod.close_output()

iterator = hazel.Iterator(use_mpi=False)
mod = hazel.Model('conf_firs.ini', working_mode='inversion', verbose=3, rank=iterator.get_rank(), randomization=5)
iterator.use_model(model=mod)
iterator.run_all_pixels()



# Do some plots

# Open the file
f = h5py.File('output_onepixel.h5', 'r')

# Check the sizes of the output
npix,nrand,nstokes,nlambda = f['spec1']['stokes'].shape
print('(npix,nrand,nstokes,nlambda) -> {0}'.format(f['spec1']['stokes'].shape))

fig, ax = pl.subplots(nrows=2, ncols=2, figsize=(10,10))
ax = ax.flatten()
for i in range(4):
    for j in range (0,nrand):
        print (i,j)
        ax[i].plot(f['spec1']['wavelength'][:] - 10830, stokes[x_coordinate,y_coordinate,i,l_left:l_right]/norm)
        ax[i].plot(f['spec1']['wavelength'][:] - 10830, f['spec1']['stokes'][0,j,i,:])

for i in range(4):
    ax[i].set_xlabel('Wavelength - 10830[$\AA$]')
    #ax[i].set_xlim([-4,3])

print(f['ch1']['Bx'])
print(f['ch1']['By'])
print(f['ch1']['Bz'])
    
pl.tight_layout()
pl.savefig('fit'+str(x_coordinate)+'_'+str(y_coordinate)+'.png',fmt='png',bbox_inches='tight')
f.close()
