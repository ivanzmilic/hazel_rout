import h5py
import matplotlib.pyplot as plt 
import numpy as np 
import sys 

#this is a simple script that outputs all the relevant physcal parameters
#for a one chromosphere run of hazel (hopefully will be extended in the future)

input_file = sys.argv[1]
output_file = sys.argv[2]
NX = int(sys.argv[3])
NY = int(sys.argv[4])

param_name = ['tau','v','deltav','Bx','By','Bz']
mapp       = ['viridis','bwr','viridis','PuOr','PuOr','PuOr']

results = h5py.File(input_file,'r')

parameters = 0

for p in range(0,len(param_name)):
	temp = results['ch1'][param_name[p]]
	temp = np.reshape(temp,[1,NX,NY,1,2,1])
	if (p==0):
		parameters = temp
	else: 
		parameters = np.append(parameters,temp,axis=0)

limits = [[0,np.max(parameters[0])],[-8,8],[0,10],[-1500,1500],[-1500,1500],[-1500,1500]]

print (parameters.shape)

for p in range(0,len(param_name)):
	plt.clf()
	plt.cla()
	plt.imshow(parameters[p,:,:,0,1,0].T,origin='Lower',cmap=mapp[p],vmin=limits[p][0],vmax=limits[p][1])
	plt.colorbar()
	plt.xlabel('x [pixel]')
	plt.ylabel('y [pixel]')
	plt.title(param_name[p])
	plt.tight_layout()
	plt.savefig(output_file+'_'+param_name[p]+'.png',bbox_inches='tight')
	plt.close('all')


