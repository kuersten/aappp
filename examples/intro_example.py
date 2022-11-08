#import package functions
from aappp import *

#set up a simulation for the standard Vicsek model
#we need some parameters: 
#interaction radius R
#noise strength eta
#particle speed v
#particle number N
#simulation box dimensions Lx, Ly
#specify boundary conditions in x- and y- direction bx, by (=0 for periodic boundaries)
#set some seed if you like
mysim=aappp_init(seed=21, R=1., v=1., eta=0.35, Lx=30., Ly=30., N=1000, bx=0, by=0)
#now the positions and orientations of all particles are setup at random
#if we do not like this, we need to change it
#thus, we first get positions and orientation (not really necessary)
[x, y, theta, lx, ly]=aappp_get_xythetalxly(mysim)
#now we might change all orientations, so all particles move in x-direction
for i in range(0, len(theta)):
    theta[i]=0.0
aappp_set_state(mysim, x, y, theta, lx, ly)
#nice, so now we iterate a few time steps with the Vicsek dynamics
VM_update_timesteps(mysim, 1000)
#lets measure some stuff averaged over 1000 steps
VM_measurement_timesteps(mysim, 1000)
#obtain the measurement results
[[theta, ptheta], [n, pn], polar, nematic, neighbors]=aappp_get_results(mysim)
#here:
#theta, ptheta is a (normailzed) histogram of orientations
#n, pn is a (normalized) histogram of the number of neighboring particles
#polar[0] is an instant value of the polar order parameter
#polar[1] ... polar[4] are the first 4 moments of the polar order parameter
#nematic is as polar but for the nematic order parameter
#neighbors[0] ... neighbors[3] are the first 4 moments of the number of neighbors
print(polar[1])
#we can improve the statistice of our measurement adding more data points
VM_measurement_timesteps(mysim, 1000)
[[theta, ptheta], [n, pn], polar, nematic, neighbors]=aappp_get_results(mysim)
print(polar[1])
#we might want to perform a 'new', independent measurement, so clear measurement results first (this does not reset the state of the particles)
aappp_reset_measurement(mysim)
#run new measurement
VM_measurement_timesteps(mysim, 1000)
#save my current state for later
aappp_save(mysim, 'test.dat')
#free memory of simulation
aappp_free(mysim)
#changed my mind, I want to continue now, so I load my file
mysim=aappp_load('test.dat')
#my measurement results are still there
[[theta, ptheta], [n, pn], polar, nematic, neighbors]=aappp_get_results(mysim)
print(polar[1])
#I can also continue the measurement
VM_measurement_timesteps(mysim, 1000)
[[theta, ptheta], [n, pn], polar, nematic, neighbors]=aappp_get_results(mysim)
print(polar[1])
#maybe I want to make a nice snapshot of the final state, so I should get the positions and orientations
[x, y, theta, lx, ly]=aappp_get_xythetalxly(mysim)
#do something nice now with x, y, theta
#now it's enough
aappp_free(mysim)
