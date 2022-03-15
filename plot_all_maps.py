import matplotlib.pyplot as plt 
import numpy as np 
import sys 
import h5py

input_file = sys.argv[1]
output_file = sys.argv[2]
NX = int(sys.argv[3])
NY = int(sys.argv[4])

fit = h5py.File(input_file,'r')
ph_cube = fit['ph1']

tau = np.squeeze(ph_cube['log_tau'])
print (tau)
NT = len(np.squeeze(ph_cube['log_tau']))
print (NT)
indices = [12,22,32]

param_name = ['T','v','Bx','By','Bz']
mapp       = ['magma','bwr','PuOr','PuOr','PuOr']

parameters = 0

for p in range(0,len(param_name)):
	temp = np.squeeze(ph_cube[param_name[p]])
	temp = temp.reshape(NX,NY,2,NT)
	temp = temp[:,:,1,indices]
	
	plt.clf()
	plt.cla()
	plt.figure(figsize=[12,6])
	
	for i in range (0,len(indices)):
		plt.subplot(2,2,i+1)	
		if (p==0 or p==1):
			m = np.mean(temp[:,:,i])
			s = np.std(temp[:,:,i])
			v_max = m+4*s
			v_min = m-4*s
		else:
			v_max = 2500.0
			v_min = -2500.0

		if (p==0):
			plt.imshow(temp[:,:,i].T,cmap=mapp[p],origin='lower')
		else:
			plt.imshow(temp[:,:,i].T,cmap=mapp[p],vmin=v_min,vmax=v_max,origin='lower')
		plt.colorbar(shrink=0.8)
		plt.xlabel('x [pixel]')
		plt.ylabel('y [pixel]')
		print (tau[indices[i]])
		where = str(tau[indices[i]])
		plt.title(param_name[p]+' at $\\log \\tau = $'+str(tau[indices[i]]))
	
	#Chromosphere:
	plt.subplot(224)
	if (p==0):
		temp = fit['ch1']['tau']
	else:
		temp = fit['ch1'][param_name[p]]
	temp = np.reshape(temp,[NX,NY,1,2,1])
	if (p==0 or p==1):
		m = np.mean(temp[:,:,0,1,0])
		s = np.std(temp[:,:,0,1,0])
		v_max = m+4*s
		v_min = m-4*s
	else:
		v_max = 2500.0
		v_min = -2500.0
	plt.imshow(temp[:,:,0,1,0].T,origin='Lower',cmap=mapp[p],vmin=v_min,vmax=v_max)
	plt.colorbar(shrink=0.8)
	plt.xlabel('x [pixel]')
	plt.ylabel('y [pixel]')
	plt.title(param_name[p]+' in the chromosphere')

	plt.tight_layout()
	plt.savefig(output_file+'_'+param_name[p]+'.png',bbox_inches='tight')
	plt.close('all')


