import aappp

#initialize simulation
#if we want the bands to move in y-direction we can impose reflecting boundary conditions in x-direction for some time
mysim=aappp.aappp_init(v=1., R=1., seed=6, N=1000000, eta=0.37, Lx=1253.3, Ly=1253.3, bx=1, by=0, kn=1000)
#perform 5*10^4 time steps
aappp.VM_update_timesteps(mysim, 50000) #wait half a day
#save intermediate state
aappp.aappp_save(mysim, 'vicsek_bands_ydirection1a.dat')
#now the bands should move in +- y direction,
#but they are not, instead the bands really get reflected at the reflecting boundary
#so our trick did not work here, but nice reflections ;)
#change reflecting to periodic boundary condition in x direction anyway and see what happens
aappp.aappp_reset_parameters(mysim, bx=0)
#perform another 5*10^4 time steps
aappp.VM_update_timesteps(mysim, 50000) # wait another half a day
#save simulation state to a file
aappp.aappp_save(mysim, 'vicsek_bands_ydirection1b.dat')
#now we produced a cross sea pattern
#for the given parameter set, it should not be stable when continuing the simulation very long
#free memory
aappp.aappp_free(mysim)
