FILES
src/aappp.c
src/c_functions_aappp.c
src/c_functions_aappp.h
README.txt
LICENSE.txt
MANIFEST.in

REQUIREMENTS
This software depends on:
*The GNU C Library under LGPL license 2.1 or later
*The GLib Library under LGPL license 2.1 or later
*Python3 under PSF License

INSTALLATION
Should work on Linux, macOS and other Unix-OS.
For compilation and installation type: bash install.bash

USAGE
#    import package
import aappp
#    get help 
help(aappp)
#    get help on specific function, e.g., aappp_init
help(aappp.aappp_init)

Help on module aappp:

NAME
    aappp - aappp aligning active particles py package

DESCRIPTION
    This package provides functions for agent-based simulations of various models for aligning active particles in two dimensions.
    Supported models are: Vicsek model, nematic Vicsek model, metric free (topological) Vicsek model, metric free (topological) nematic Vicsek model and overdamped Langevin models of the type:
    d/dt x_i= v*cos(phi_i),
    d/dt y_i= v*sin(phi_i),
    d/dt phi_i=sum_{j is neighbor of i}*gamma*sin((phi_j-phi_i)*order)+omega+eta*xi_i,
     where the neighborhood is either defined by an interaction distance, or metric free (topological) interactions are considered, where each particle interacts with the kn nearest neighbors.
    In each model, chirality can be added (omega) and one can simulate multiple particle species, where each species can have a different value for noise strength, velocity, chirality, coupling (=coupling matrix).
    Use function aappp_init to initialize a simulation, see >>>help(aappp_init)<<<.
    Use function ***_update_timesteps to evolve the system in time.
    *** can be VM, NVM, mfVM, mfNVM for Vicsek models and additiveL, nonadditiveL or mfL for overdamped Langevin models, see >>>help(***_update_timesteps)<<< for exact model definitions.
    Use function ***_measurement_timesteps to evolve the system in time and measure polar order (first 4 moments), nematic order (first 4 moments), number of neighbors (first four moments), histogram of the orientation and a histogram of the number of neighbors, see >>>help(***_measurement_timesteps)<<<.
    Use function aappp_get_results to obtain the results of the measurements, see >>>help(aappp_get_results)<<<.
    Use function aappp_get_xythetalxly to get positions and orientation of all particles and simulation box size, see >>>help(aappp_get_xythetalxly)<<<.
    Use function aappp_get_parameters to obtain the simulation parameters, see >>>help(aappp_get_parameters)<<<.
    Use function aappp_set_state to manipulate the positions and orientations of the particles and/or resize the simulation box, see >>>help(aappp_set_state)<<<.
    Use function aappp_reset_parameters to change some simulation parameters, see >>>help(aappp_reset_parameters)<<<.
    Use function aappp_reset_measurement to clear all measurements and possibly change the number of bins for histograms, see >>>help(aappp_reset_measurement)<<<.
    Use function aappp_save to save the state of the simulation (including measurement resuts and state of the pseudo random number generator) to a file, see >>>help(aappp_save)<<<.
    Use function aappp_load to load a previous simulation, see >>>help(aappp_load)<<<.
    Use function aappp_get_weights to obtain a list of weights used for nonadditive overdamped Langevin model, see >>>help(aappp_get_weights)<<<.
    Parallel computation is not supported by now. However, performence and memory usage should be good.
    
    EXAMPLE of use
    
    mysim=aappp.aappp_init(N=10000, eta=0.2, v=1., R=1., Lx=40., Ly=80.)
    aappp.VM_update_timesteps(mysim, 10000)
    aappp.aappp_measurement_timesteps(mysim, 1000)
    aappp.aappp_get_results(mysim)
    x, y, theta, lx, ly=aappp.aappp_get_xythetalxly(mysim)
    aappp.aappp_save(mysim, 'mydata.dat')
    aappp.aappp_free(mysim)

FUNCTIONS
    NVM_measurement_timesteps(simulation, timesteps)
        performs *timesteps* timesteps with the same dynamics as NVM_update_timesteps, see >>>help(NVM_update_timesteps)<<< and measures the same quantities as VM_measurement_timesteps, see >>>help(VM_measurement_timesteps)<<<.
    
    NVM_update_timesteps(simulation, timesteps)
        as VM_update_timesteps, see >>>help(VM_update_timesteps)<<<, but with different alignment rule 1):
        first, the new angle is determined as before:
        theta_new_i=arg(sum_{j is neighbor of i} exp(imaginary unit*theta_i))
        but then the angle is rotated by pi, if the previous orientation is closer to this rotated angle than to theta_new_i:
        theta_new_i=theta_new_i*sign(sin(theta_new_i)*sin(theta_i)+cos(theta_new_i)*cos(theta_i)).
    
    VM_measurement_timesteps(simulation, timesteps)
        follows the Vicsek model dynamics for *timesteps* timesteps as VM_update_timesteps, see >>>help(VM_update_timesteps)<<<, however, certain quantities are measured in each step.
        The measured quantities are:
        polar order parameter: p=|sum_{i}exp(imaginary unit*theta_i)|/N,
        nematic order paremeter:q=|sum_{i}exp(imaginary unit*2*theta_i)|/N,
        the number of neighbors: for each particle the number of (different) neighbors is counted.
        for polar and nematic order parameter the first to fourth moment is calculated from the time series, for the number of neighbors these four moments are also calculated averaging over both, all particles and the time series (useful to calculate Binder parameter).
        The averaging is performed over *timesteps* timesteps counting the initial configuration (before the first update), but not the last configuration (after the last update) in order to avoid to count certain states twice when calling VM_measurement_timesteps multiple times.
        For the particle orientations (angles) and the number of neighbors, additionally to the first four moments, also a histogram is recorded.
        The angles are defined in such a way that they are in the interval [-pi, pi].
        The number of neighbors can in theory be any number from 0 to N-1.
        The number of bins for the histograms can be specified in aappp_init, see >>>help(aappp_init)<<< or aappp_reset_measurement, see >>>help(aappp_reset_measurement)<<<.
        If there is more than one particle species, polar and nematic order, number of neighbors and corresponding histogram are measured for each particle species separately.
        The results of the measurement are stored in the *simulation* data structure.
        The measurement results can be accessed by calling aappp_get_results, see >>>help(aappp_get_results)<<<.
        Calling VM_measurement_timesteps mutliple times will average over all timesteps of the multiple calls (except for the very last step).
        Normalization of the moments and histograms is done only when calling aappp_get_results, see >>>help(aappp_get_results)<<<.
        Thus calling VM_measurement_timesteps(simulation, timesteps) twice introduces no rounding errors compared to a single call of VM_measurement_timesteps(simulation, 2*timesteps).
        After a call of VM_measurement_timesteps one can delete the measurement results (in order to perform a new measurement later) if desired, by calling aappp_reset_measurement, see >>>help(aappp_reset_measurement)<<<.
    
    VM_update_timesteps(simulation, timesteps)
        Performs *timesteps* timesteps of the Vicsek model dynamics with parameters specified before.
        Returns polar order parameter after penultimate time step.
        Both arguments must be given.
        The first argument *simulation* must specify a simulation state that was created earlier by the aappp_init function, see >>>help(aappp_init)<<< or loaded via aappp_load, see >>>help(aappp_load)<<<.
        The second argument *timesteps* specifies the number of timesteps that are performed. The final state after those steps is stored in the *simulation* data structure and can be accessed by the aappp_get_xythetalxly function, see >>>help(aappp_get_xythetalxly)<<<.
        
         The Vicsek dynamics that is performed consists of two parts: 1) alignment and 2) streaming. In each step 1) is performed first and 2) afterwards. In 1) for each particle i (with position x_i, y_i, orientation theta_i), all neighboring particles j are found such that (x_i-x_j)**2 + (y_i-y_j)**2<R**2. Note that this holdes true for particle i itself. Then a new angle is calculated as:
        theta_new_i=arg(sum_{j is neighbor of i} exp(imaginary unit*theta_i)).
        The new angle is disturbed by a random number:
         theta_new_i = theta_new_i + xi_i,
        where xi_i are independently drawn from a uniform distribution on [-eta*pi, eta*pi]. After updating all orientations 1) is finished. In 2) all positions are updated according to:
        x_i=x_i+v*cos(theta_new_i)
        y_i=y_i+v*sin(theta_new_i)
        and the old orientations can be forgotten:
        theta_i=theta_new_i.
    
    aappp_free(simulation)
        frees all memory that is allocated to the *simulation* data structure.
        All data get lost.
    
    aappp_get_parameters(simulation)
        returns: [v, eta, R, omega, lx, ly, gamma, dt, N, kn, order, binnum_theta, binnum_neighbors, boundary_x, boundary_y],
        where *v* is particle speed, *eta* is noise strength, *R* is interaction radius, *omega* is chirality, *lx* is simulation box size in x-direction, *ly* is simulation box size in y-direction, *gamma* is the interaction strength, *dt* is step size, *N* is particle number, *kn* is the number of distinct interacting neighbors in metric free/topological models, *order* is the interaction order for overdamped Langevin models, *binnum_theta* is the number of bins that are used for orientation histograms, *binnum_neighbors* is the number of bins used in number of neighbor histograms, *boundary_x* is a flag for boundary conditions in x-direction (0->periodic boundary conditions, 1->reflecting boundary conditions), *boundary_y* is a flag for boundary conditions in y-direction (0->periodic boundary conditions, 1->reflecting boundary conditions).
        Not all parameters are used in each model (dynamics).
        For example, *dt*, *order*, *gamma* are not used in Vicsek typ models (VM, NVM, mfVM, mfNVM), where the time evolution is a discrete map.
        *kn* is not used in metric models (VM, NVM, additiveL, nonadditive_L)
        and *R* is not used in metric free/topological models (mfVM, mfNVM, mfL).
        For periodic boundary conditions (here explained in x-direction) particles that leave x in [0, lx] are put to x->x+lx or x->x-lx such that they are back inside the simulation box.
        For calculating the distance between particles i and j, the minimum of |x_i-x_j|, |x_i-x_j-lx| and |x_i-x_j+lx| is used.
        For reflecting boundary conditions (here explained in x-direction) particles that leave x in [0, lx] are put to either -x or -x + 2*lx such that they are back inside the simulation box [0, lx],
         when this happens, the particle orientation is reflected to theta->-theta+pi.
        Particles interact with all real particles but also with image particles that are the mirror image under the above reflection (of position and orientation).
        If there is more than one particle species, then *v*, *eta*, *omega*, *N* are lists containing the corresponding quantities of all particle species and *gamma* is a list of lists (coupling matrix).
    
    aappp_get_results(simulation)
        if there is only one particle species it returns: [[ptheta, theta], [pn, n], polar, nematic, neighbors],
        where *ptheta* is the normalized histogram of particle orientations and the list *theta* is composed of the centers of the bins of the corresponding histogram.
        The list *pn* is the normalized histogram of the number of neighbor distribution and the list *n* gives the corresponding histogram bins 0, 1, 2, ..., n_max.
        Events of particles that have more than n_max neighbors are not recorded.
        Thus *pn* is not normalized to one as soon as such events are missed.
        The list *polar* contains the polar order parameter from the last time measured (after penultimate time step) in the first entry and the first four moments of the polar order parameter in the next four entries.
        The list *nematics* contains the nematic order parameter in the first entry and the first four moments of the nematic order parameter in the next four entries.
        The list *neighbors* contains the first four moments of the number of neighbor distribution.
        If there is more than one particle species, the function returns:[[[ptheta1, theta1], [ptheta2, theta2], ...], [[pn1, n1], [pn2, n2], ...], [polar1, polar2, ...], [nematic1, nematic2, ...], [neighbors1, neighbors2, ...]].
        Those quantities are analogous to the single species case, but with the same results for each species.
    
    aappp_get_weights(simulation)
        returns list of weights that are used for nonadditive overdamped Langevin dynamics nonadditiveL_update_timesteps, see >>>help(nonadditiveL_update_timesteps)<<<.
        if an empty list is returned, the standard weight function weight(n)=1/(n+1) is used.
    
    aappp_get_xythetalxly(simulation)
        returns [x, y, theta, lx, ly],
        where x, y and theta are lists of length N containing x-position, y-position and orientation of all particles.
        lx and ly give the size of the simulation box in x- and y-direction.
    
    aappp_init(v, eta, R, omega, Lx, Ly, gamma, dt, N, seed, kn, order, binnum_theta, binnum_neighbors, bx, by, weight_function, weight_vector_length)
        Returns a newly created simulation object.
        
        Arguments are optional and can be given with keyword (default values in brackets):
        Integer particle Number: N(1000)
        Float interaction radius: R(1)
        Float simulation box length in x-direction: Lx(20)
        Float simulation box length in y-direction: Ly(20)
        Float step size for Euler integration: dt(0.01)
        Integer boundary flag for x direction (if =0 ->periodic boundary, if =1 ->reflecting boundary): bx(0)
        Integer boundary flag for y direction (if =0 ->periodic boundary, if =1 ->reflecting boundary): by(0)
        Integer order of the interaction term for overdamped Langevin dynamics (the term has the form sin(order*(phi2-phi1))): order(1)
        Integer seed of the pseudo random number generator: seed(1)
        Integer number of interaction neighbors other than focus particle for metric free models: kn(6)
        Integer number of bins for angular histogram: binnum_theta(100)
        Integer number of bins for histogram of the number of neighbors: binnum_neighbors(100)
        Float speed of the particles: v(1)
        Float noise strength (which should be a number between 0 and 1 for Vicsek models and an arbitrary number for overdamped Langevin models): eta(0.5)
        Float natural rotation frequency for each particle: omega(0)
        Float coupling strength for overdamped Langevin models: gamma(0.1)
        Python function that gives a weight to interactions depending on the number of neighbors: weight_function (None)
        Ineger specifying for how many values the weight function is calcuated: weight_vector_length (100)
        
        N, v, eta, omega can be lists where each entry corresponds to one species of particles (different species behave differently).
        gamma can be a list of lists where gamma[A][B] describes, in the equation of motion of an A-particle, the strength of interaction with a B-particle.
        It is possible that e.g. N and v are lists (of the same dimension) and eta, omega ang gamma are no lists.
        Lists of length 1 are not allowed (use just the number in that case).
        If lists are used, the dimensions for those parameters must be compatible.
        Other parameters than listed before can only be single numbers.
        The weight function is only relevant for the nonadditive overdamped Langevin dynamics, see >>>help(nonadditiveL_update_timesteps)<<<.
        The weight function must take an integer as a argument and return a float.
        If no weight function is given f(n)=1/(n+1) is used.
        The weight function is only called once for the arguments 0, 1, 2, ..., weight_vector_length-1.
        The results of those calls are stored in memory and used during simulation.
        If a particle has more than weight_vector_length-1 neighbors the weight for exactly weight_vector_length-1 neighbors is used.
        particles are initalized at random with unform distribution over space and with uniform distributed orientations.
        
        EXAMPLE 1 (with only one particle species)
        mysim=aappp.aappp_init(N=100000, Lx=100., Ly=100., v=1.5, R=1.1, eta=0.23, seed=7)
        
        EXAMPLE 2 (with three particle species)
        mysim=aappp.aappp_init(N=[10000, 20000, 15000], gamma=[[0.1, 0.13, 0.15],[0.08, 0.15, 0.05],[0.11, 0.14, 0.11]], Lx=80., Ly=120., dt=0.03, v=1, eta=1.5, seed=1)
        In this example, if particle i is from species two and particle j is from species three, and i and j are neighbors, there appears the following coupling term (depending on the model, potentially with a weighting prefactor depending on the number of neighbors of particle i):
        d/dt phi_i = 0.05*sin((phi_j-phi_i)*order) + ...
    
    aappp_load(filename)
        loads and returns a simulation data structure from a file under path *filename*, that was previously saved using aappp_save, see >>>help(aappp_save)<<<.
    
    aappp_reset_measurement(simulation, binnum_theta, binnum_neighbors)
        Forgets all measured data such as moments of polar order, etc. and histograms.
        The number of bins of orientation histograms (*binnum_theta*) or number of neighbor histograms (*binnum_neighbors*) for following measurements can be set optionally by keyword.
        Returns 1 on success.
    
    aappp_reset_parameters(simulation, v, R, eta, omega, dt, kn, gamma, order, bx, by, N, weight_function, weight_vector_length)
        Is used to change (some) simulation parameters for the following time steps.
        Except for the simulation data structure (*simulation*) all arguments are otional and can be given by keyword.
        The total particle number can not be changed.
        In order to change the total particle number, a new simulation data structure has to be created using aappp_init, see >>>help(aappp_init)<<<.
        However, one can change the number of species and the number of particles in each species as long as the total particle number is not changed.
        This can be done by giving a list for *N*=[N1, N2, ...].
        Similarly, *v*, *omega* and *eta* can be lists specifying velocity, chirality and noise strength for each species.
        *gamma* can be a list of lists specifying the coupling matrix.
        Returns 1 on success.
    
    aappp_save(simulation, filename)
        saves the *simulation* data structure containing all particles positions and orientations, all parameters, the state of the pseudo random number generator and the results of measurements into a file under path *filename*.
        All data are saved in binary format from C-data types.
        Note that e.g. C doubles are not type safe.
        Thus, performing simulation on unusual devices such as mirco processors might result in data files that are unreadable on a different device such as a desktop computer.
        However, such behavior is unexpected when using usual desktop computers or high performence computing clusters only.
        In case that such problems occur it is recommended to analyze data on the same machine that performed the simulation.
        Returns 1 on success.
    
    aappp_set_state(simulation, x, y, theta, lx, ly)
        allows to change the state of the simulation.
        *x*, *y*, and *theta* are lists of x-positions, y-positions and orientations of N particles.
        the length of *x*, *y*, *theta* must be exactly N.
        It is not possible to modify the particle number of a given simulation.
        If it is desired to change the particle number one needs to create a new simulation data structure using aappp_init, see >>>help(aappp_init)<<<.
        *lx* and *ly* are the new simulation box sizes in x- and y-direction.
        All arguments must be given.
        If particle positions are outside simulation box or orientations are outside [-pi, pi], they are put inside via periodic boundary conditions.
        Returns 1 on success, -1 when failing.
    
    additiveL_measurement_timesteps(simulation, timesteps)
        performs *timesteps* timesteps as additiveL_update_timesteps, see >>>help(additiveL_update_timesteps)<<< and measures observables as VM_measurement_timesteps, see >>>help(VM_measurement_timesteps)<<<.
    
    additiveL_update_timesteps(simulation, timesteps)
        performs *timesteps* timesteps for the following overdamped Langevin dynamics:
        d/dt x_i=v*cos(phi_i)
        d/dt y_i=v*sin(phi_i)
        d/dt phi_i=gamma*sum_{j is neighbor of i}sin((phi_j-phi_i)*order) + omega + eta*xi_i
        
        two particles are considered to be neighbors if their distance is smaller than *R* as in the Vicsek model, see >>>help(VM_update_timesteps)<<<.
        If there is more than one particle species, v, omega and eta can be different for each species and gamma can depend on the species of i and j (coupling matrix), see >>>help(aappp_init)<<< for parameter definitions.
        The overdamped Langevin equations are integrated with Euler-Maruyama-scheme with step size *dt*.
        *timesteps* denotes the number of timesteps performed, such that the real time of time evolution is *timesteps*  *  *dt*.
        Returns polar order parameter before last step.
    
    mfL_measurement_timesteps(...)
        additiveL_measurement_timesteps(simulation, timesteps)
        --
        
        performs *timesteps* timesteps as mfL_update_timesteps, see >>>help(mfL_update_timesteps)<<< and measures observables as VM_measurement_timesteps, see >>>help(VM_measurement_timesteps)<<<.
    
    mfL_update_timesteps(simulation, timesteps)
        evolves the system for *timesteps* timesteps as additiveL_update_timesteps, see >>>help(additiveL_update_timesteps)<<<, but with metric free/topological definition of neighborhoods, as in mfVM_update_timesteps, see >>>help(mfVM_update_timesteps)<<<.
    
    mfNVM_measurement_timesteps(simulation, timesteps)
        performs *timesteps* timesteps with the same dynamics as mfNVM_update_timesteps, see >>>help(mfNVM_update_timesteps)<<< and measures the same quantities as VM_measurement_timesteps, see >>>help(VM_measurement_timesteps)<<<.
        Note that the number of neighbors measurements are pointless in this model because the number of neighbors is alsways the same.
        However, the number of neighbor measurements are nevertheless performed in the same way.
    
    mfNVM_update_timesteps(simulation, timesteps)
        as VM_update_timesteps, see >>>help(VM_update_timesteps)<<<, but with metric free/topological definition of neighborhoods as in mfVM_update_timesteps, see >>>help(mfVM_update_timesteps)<<<, and with a nematic interaction rule as in NVM_update_timesteps, see >>>help(NVM_update_timesteps)<<<.
    
    mfVM_measurement_timesteps(simulation, timesteps)
        performs *timesteps* timesteps with the same dynamics as mfVM_update_timesteps, see >>>help(mfVM_update_timesteps)<<< and measures the same quantities as VM_measurement_timesteps, see >>>help(VM_measurement_timesteps)<<<.
        Note that the number of neighbors measurements are pointless in this model because the number of neighbors is alsways the same.
        However, the number of neighbor measurements are nevertheless performed in the same way.
    
    mfVM_update_timesteps(simulation, timesteps)
        as VM_update_timesteps, see >>>help(VM_update_timesteps)<<<, but with a different definition of neighborhoods (called metric free or topological).
        Here particle i is considered as neighbor of particle i (as before).
        The other neighbors are the *kn* particles that are closest to particle i.
        Thus the parameter *kn* should have values 1, 2, 3, ..., see >>>help(aappp_init)<<< for initialization of the parameters.
    
    nonadditiveL_measurement_timesteps(...)
        additiveL_measurement_timesteps(simulation, timesteps)
        --
        
        performs *timesteps* timesteps as nonadditiveL_update_timesteps, see >>>help(nonadditiveL_update_timesteps)<<< and measures observables as VM_measurement_timesteps, see >>>help(VM_measurement_timesteps)<<<.
    
    nonadditiveL_update_timesteps(simulation, timesteps)
        dynamics is iterated as in additiveL_update_timesteps, see >>>help(additiveL_update_timesteps)<<<.
        here, however, the equation of motion of the orientation is:
        d/dt phi_i=gamma*weight(n)*sum_{j is neighbor of i}sin((phi_j-phi_i)*order) + omega + eta*xi_i,
        where n is the number of distinct neighbors of particle i (not counting i itself).
        if no weight function weight(n) is specified, the function weight(n)=1/(n+1) is used for all occuring values of n.
        the weight function weight(n) can be specified using aappp_init or aappp_reset_parameters, see >>>help(aappp_init)<<< or >>>help(aappp_reset_parameters)<<<.
        if a weight function is specified, also a value for *weight_vector_length* must be specified.
        in that case the weight function is used for all values of n in the range 0, 1, 2, ..., *weight_vector_length*-1
        for n>=weight_vector_length, the weight weight(*weight_vector_length*-1) is used.

FILE
    /usr/local/lib/python3.10/site-packages/aappp.cpython-310-darwin.so


None
