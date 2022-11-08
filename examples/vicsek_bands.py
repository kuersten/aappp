import aappp

#initialize simulation
mysim=aappp.aappp_init(v=1., R=1., seed=1, N=1000000, eta=0.37, Lx=1253.3, Ly=1253.3, bx=0, by=0)
#perform 10^5 time steps of standard Vicsek model
aappp.VM_update_timesteps(mysim, 100000) #wait about 1 day
#save simulation state to a file
aappp.aappp_save(mysim, 'vicsek_bands.dat')
#observe some nice bands, there is still a defect that should disappear after longer time
#note that in most cases, the bands move not in x- or y- direction but in a somewhat diagonal direction when started from random initial conditions
#when particles are initiated e.g. moving in x- direction they usually loose this orientation initially and start to collectively move in some other (most likely) diagonal direction
#free memory
aappp.aappp_free(mysim)
