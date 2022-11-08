import aappp

mysim=aappp.aappp_init(v=0.2, eta=0.08, N=320000, Lx=2000, Ly=400, bx=0, by=1, omega=0., kn=2)
#lets use the Vicsek dynamics with metric free/ topological neighborhoods
#we want the particles to move in x- direction
#thus we apply reflecting boundaries in y-direction first
aappp.mfVM_update_timesteps(mysim, 50000) #wait 4h
aappp.aappp_save(mysim, 'mfVM1a.dat')
#now we can change boundary conditions in y direction into periodic
aappp.aappp_reset_parameters(mysim, by=0)
aappp.mfVM_update_timesteps(mysim, 50000) # wait another 4h
aappp.aappp_save(mysim, 'mfVM1b.dat')
#we see some nice band here as well
aappp.aappp_free(mysim)
