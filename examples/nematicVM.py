import aappp

mysim=aappp.aappp_init(N=524288, v=0.5, eta=0.065, Lx=2048., Ly=2048., bx=0, by=0, omega=0.)
aappp.NVM_update_timesteps(mysim, 100000) #wait about 11h
aappp.aappp_save(mysim, 'nematicVM.dat') #there is a nematic band structure in a diagonal direction, not yet fully relaxed
aappp.aappp_free(mysim)
