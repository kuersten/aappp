import aappp

#if we want the band structure to form itself along x- or y- direction we can use reflecting boundary conditions in one direction (here y) to avoid diagonal bands
mysim=aappp.aappp_init(N=524288, v=0.5, eta=0.065, Lx=2048., Ly=2048., bx=0, by=1, omega=0.)
aappp.NVM_update_timesteps(mysim, 100000) #wait about 10h
aappp.aappp_save(mysim, 'nematicVM2a.dat') # save result
aappp.aappp_reset_parameters(mysim, by=0) #now the band should be more or less in x- or y-direction, use periodic boundary conditions in both directions from now on
aappp.NVM_update_timesteps(mysim, 100000) #wait about 10h
aappp.aappp_save(mysim, 'nematicVM2b.dat') # watch a nice nematic band
aappp.aappp_free(mysim)
