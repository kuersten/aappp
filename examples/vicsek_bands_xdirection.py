import aappp

#initialize simulation
#if we want the bands to move in x-direction we can impose reflecting boundary conditions in y-direction for some time
mysim=aappp.aappp_init(v=1., R=1., seed=1, N=1000000, eta=0.37, Lx=1253.3, Ly=1253.3, bx=0, by=1)
#perform 5*10^4 time steps
aappp.VM_update_timesteps(mysim, 50000) #wait half a day
aappp.aappp_save(mysim, 'vicsek_bands_xdirection1a.dat')
#now the bands should move in +- x direction, we can change boundary conditions to periodic
aappp.aappp_reset_parameters(mysim, by=0)
#perform another 5*10^4 time steps
aappp.VM_update_timesteps(mysim, 50000) #wait another half a day
#save simulation state to a file
aappp.aappp_save(mysim, 'vicsek_bands_xdirection1b.dat')
#now we produced bands nicely moving in x-direction
#free memory
aappp.aappp_free(mysim)
