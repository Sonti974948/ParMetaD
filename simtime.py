import h5py

import numpy as np



simtime = 0.01 # in ns

f = h5py.File('west_now.h5','r')

total_simtime = np.sum(f['summary']['n_particles'])*simtime

total_walltime = np.sum(f['summary']['walltime'])/60.0 # in hours

speed=total_simtime/total_walltime

print('total simtime = '+str(total_simtime)+' nanoseconds')

print('total walltime = '+str(total_walltime)+' minutes')

print('Speed = %.2f ns/hr = %.2f ns/day'%(speed*60,speed*60*24))

print(f['summary']['walltime'])
print(f['summary']['n_particles'])
