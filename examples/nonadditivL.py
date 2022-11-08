import aappp

mysim=aappp.aappp_init(N=19881, Lx=19881.**0.5, Ly=19881.**0.5, R=1., gamma=1., v=1., eta=0.5**0.5, bx=0, by=0, dt=0.003)
#not only the Vicsek model can produce nice bands
#models with Brownian dynamics behave qualitatively similar
aappp.nonadditiveL_update_timesteps(mysim, 1000000) #wait 4h
aappp.aappp_save(mysim, 'nonadditiveL.dat')
#observe a nice band
aappp.aappp_free(mysim)

