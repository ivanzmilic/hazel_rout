import h5py
import numpy as np
from astropy.io import fits
import sys

# Read the data from the hazel output and pack it up as a fits file with all the necessary information
# Also supply the information in the header

input_file = sys.argv[1] # name of the ORIGINAL fits file
input_file = input_file[:-5]


# A bit brute force but will work:
original_fits = fits.open(input_file+".fits")
NX = int(original_fits[0].header['NAXIS4'])
NY = int(original_fits[0].header['NAXIS3'])
print("NX = ", NX, "NY = ", NY)


results = h5py.File(input_file+"_l2.h5",'r')

# The output contains relevant parameters an their uncertainties + chi2 + a specific 'flag' entry that can specify additionnal info:
# Optical depth + error (unitless)
# velocity + error (km/s)
# line broadening + error (km/s)
# line damping + error (no units)
# beta + erorr (no units) : this is an ad-hoc parameter that will 99% be hard-set to 1 in final version. But let's have it here 
# 	just in case!
# Bx + error (Gauss)
# By + error (Gauss)
# Bz + error (Gauss)
# Chi2 
# flag : bad/good/ something?

atmos = np.zeros([18,NX*NY])

atmos[0] = np.copy(results['ch1']['tau'])[:,0,1,0]
atmos[1] = np.copy(results['ch1']['tau_err'])[:,0,1]
atmos[2] = np.copy(results['ch1']['v'])[:,0,1,0]
atmos[3] = np.copy(results['ch1']['v_err'])[:,0,1]
atmos[4] = np.copy(results['ch1']['deltav'])[:,0,1,0]
atmos[5] = np.copy(results['ch1']['deltav_err'])[:,0,1]
atmos[6] = np.copy(results['ch1']['a'])[:,0,1,0]
atmos[7] = np.copy(results['ch1']['a_err'])[:,0,1]
atmos[8] = np.copy(results['ch1']['a'])[:,0,1,0]
atmos[9] = np.copy(results['ch1']['a_err'])[:,0,1]
atmos[10] = np.copy(results['ch1']['Bx'])[:,0,1,0]
atmos[11] = np.copy(results['ch1']['Bx_err'])[:,0,1]
atmos[12] = np.copy(results['ch1']['By'])[:,0,1,0]
atmos[13] = np.copy(results['ch1']['By_err'])[:,0,1]
atmos[14] = np.copy(results['ch1']['Bz'])[:,0,1,0]
atmos[15] = np.copy(results['ch1']['Bz_err'])[:,0,1]
atmos[16] = np.copy(results['spec1']['chi2'])[:,0,0]
atmos[17] = 1

atmos = atmos.reshape(18,NX,NY)
# The easy part is done. Hard part is filling up the fits header

# Can I create header from scratch? Let's try the dumb way first
output = fits.PrimaryHDU(atmos)
header = output.header

#Hardcoded proxy values to be written to the header (read this from the config file?)
# Now we now how to do this
# TODO on Ivan to do for the next week's meeting:
topology = 'ch1'
height = 5.0
noise = 1E-3
NCYC = 2
Bx0 = 50.0
By0 = 50.0
Bz0 = 50.0
tau = 1.0
v = 0.0
deltav = 5.0
beta = 1.0
a = 0.001
lleft = 250
lright = 370

header['TOPO'] = (topology, 'Topology of the slab and underlying atmosphere')
header['HEIGHT'] = (height, 'Height above the solar surface')
header['NOISE'] = (noise, 'Noise given to Hazel')
header['NCYC'] = (NCYC, 'Number of cycles used in the inversion')
#header['IWEIGHT'] = (weigthts[0], 'W') # What to do with weights? There is a lot of redundant data here (see config file)
header['Bx0'] = (Bx0,'Initial value of Bx')
header['By0'] = (By0,'Initial value of By')
header['Bz0'] = (Bz0,'Initial value of Bz')
header['tau'] = (tau,'Initial value of tau')
header['v'] = (v,'Initial value of v')
header['deltav'] = (deltav,'Initial value of deltav')
header['beta'] = (beta,'Initial value of beta')
header['a'] = (a,'Initial value of a')
header['llo'] = (lleft,'Lower boundary of the wavelength grid')
header['lhi'] = (lright,'Higher bounday of the wavelength grid')

# Actually hardcoded description of the parameters:
header['P0'] = ('Bx','---')
header['P1'] = ('Bxerr','---')
header['P2'] = ('By','---')
header['P3'] = ('Byerr','---')
header['P4'] = ('Bz','---')
header['P5'] = ('Bzerr','---')
header['P6'] = ('tau','---')
header['P7'] = ('tauerr','---')
header['P8'] = ('v','---')
header['P9'] = ('verr','---')
header['P10'] = ('deltav','---')
header['P11'] = ('deltaverr','---')
header['P12'] = ('beta','---')
header['P13'] = ('betaerr','---')
header['P14'] = ('a','---')
header['P15'] = ('aerr','---')
header['P16'] = ('chi2','---')
header['P17'] = ('control','---')

print (header)

# Notes about the header
# X,Y and lambda coordinates that describe what subset of the original fits file has been taken into account for inversion

output.writeto(input_file+"_l2.fits",overwrite='True')





