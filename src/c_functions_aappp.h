/*
 * =====================================================================================
 *
 *       Filename:  c_functions_aappp.h
 *
 *    Description:  
 *
 *        Version:  1.2.1
 *        Created:  12/12/2022 
 *
 *         Author:  RUEDIGER KUERSTEN 
 *
 *        License:
Copyright 2022 RUEDIGER KUERSTEN

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 *
 * =====================================================================================
 */

#ifndef C_FUNCTIONS_AAPPP //double inclusion protection
#define C_FUNCTIONS_AAPPP
#include <Python.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <glib.h>
#include <fcntl.h>
#include <complex.h>


//data structure containing degrees of freedom of each particle
struct VM_particle{
	gdouble x;	//x-position of particle
	gdouble y;	//y-position of particle
	gdouble theta;	//orientation of particle
	gdouble theta_new;	//future orientation of particle
	struct VM_particle* next_particle; //pointer to the next particle
	guint64 index;	//index numbering the particles
};

//structure for metric free models containing a list of nearest neighbors
struct mfVM_neighborhood{
	guint64 * neighbor_indexes; //list of indexes of so far known nearest neighbors
	gdouble * distances; //distances of the neighbors from the above list
	gint32 * mirror_flag;	// flag that tells if particle direction needs to be mirrored, if=0 ->no mirroring, if=1 ->mirror in x direction, if=2 ->mirror in y direction, if=3 -> mirror in x- and y-direction, always zero if periodic boundary conditions are used
	gdouble max_distance_checked; //distance in which there is for sure no neighbor that is closer than the neighbors in the above list
};

//data structure containing the simulation state including particles degrees of freedom, state of the pseudo random number generator and weights to be used with nonadditive langevin model
struct VM_state{
	struct mfVM_neighborhood nh;	//a structure that contains a list of nearest neighbors -> only needed for metric free model
	gdouble length_x;	//length of simulation box in x-direction
	gdouble length_y;	//length of simulation box in y-direction
	gdouble gauss1;	//gaussian random number 1 (they are produced pairwise by box-muller algorithm)
	gdouble gauss2;	//gaussian random number 2 (they are produced pairwise by box-muller algorithm)
	guint64 particle_number;	//number of particles
	gint32 box_muller;	//internal state variable that determines if new gaussian random variables need to be produced
	gint32 weight_vector_length;	//number of weights that are saved in advance for nonadditive Langevin model; this is the length of  weights
	gdouble * weights;	//weights for nonadditive Langevin model
	PyObject * weight_function_callback;	//python function object that provides the weight function to calculate the weights array
	struct VM_particle * particles;	//pointer to the states of the particles
	GRand * prng;	//pseudo random number generator state variable of glib
};

//data structure containing the simulation parameters
struct VM_parameters{
	gdouble interaction_radius;	// radius R: particles are neighbors if their distance is smaller than R
	gdouble mf_interaction_radius;	// effective radius for metric free models
	gdouble delta_t;	//time step for Langevin models
	gdouble sqrt_delta_t;	//square root of time step for Langevin models
	gdouble box_size_x;	//size of a box in x-dimension  
	gdouble box_size_y;	//size of a box in y-dimension
	gdouble mf_box_size_x;	//size of a box in x-dimension for metric free model
	gdouble mf_box_size_y;	//size of a box in y-dimension for metric free model
	gint32 box_number_x;	//number of boxes in one row
	gint32 box_number_y;	//number of boxes in one column -> total number of boxes=box_number_x*box_number_y
	gint32 mf_box_number_x;	//number of boxes in one row for metric free model
	gint32 mf_box_number_y;	//number of boxes in one column -> total number of boxes=mf_box_number_x*mf_box_number_y for metric free model
	gint32 mf_kn;	//for the metric free model: number of interaction partners (different from focus particle itself)
	gint32 interaction_order;	//order of interaction for langevin models: 1-> polar alignment, 2-> nematic alignment
	gint32 boundary_x;	//flag that describes boundary conditions used in x-direction: 0 -> periodic boundary, 1-> reflecting boundary
	gint32 boundary_y;	//flag that describes boundary conditions used in y-direction: 0 -> periodic boundary, 1-> reflecting boundary
	gint32 particle_species;	//number of different species that should be treated differently (different coupling, different velocities, different noise, different rotation frequencies)
	guint64 * species_particle_number;	//vector of particle numbers for each species (not needed if particle_species=1)
	gdouble * omega;	//chirality: increment of angle for each step for Vicsek models, frequency of rotation for Langevin models
	gdouble * noise_strength; // noise strength from the interval [0,1] for Vicsek models, arbitrary numbers for Langevin models
	gdouble * speed;	//velocity of each particle
	gdouble ** coupling;	//coupling constant for Langevin models
	struct VM_particle** box;	//2-d array of lists of particles (one list of particles for each box)
	struct VM_particle** mf_box;	//2-d array of lists of particles (one list of particles for each box) for metric free model
};

//data structure containing the results of measurements (averaged values of observables, histograms)
struct VM_observables{
	gdouble * polar_order;
	gdouble * polar_order_moment1;
	gdouble * polar_order_moment2;
	gdouble * polar_order_moment3;
	gdouble * polar_order_moment4;
	gdouble * nematic_order;
	gdouble * nematic_order_moment1;
	gdouble * nematic_order_moment2;
	gdouble * nematic_order_moment3;
	gdouble * nematic_order_moment4;
	gdouble * neighbors_moment1;
	gdouble * neighbors_moment2;
	gdouble * neighbors_moment3;
	gdouble * neighbors_moment4;
	guint64 measurement_steps;
	gint32 binnum_theta;
	gint32 binnum_neighbors;
	guint64 ** hist_theta;
	guint64 ** hist_neighbors;
	gint32 sf_mode_number;	//number of modes to compute structure factor
	gint32 * kx;	//x-components of modes to compute structure factore
	gint32 * ky;	//y-components of modes to compute structure factore
	double complex ** fourier_density;	//store fourier transform of density for each species
	double complex ** structure_factor;	//store structure factor for each species
};

//data structure containing the full simulation state including all previous structures
struct VM_simulation{
	struct VM_state * state;
	struct VM_parameters * parameters;
	struct VM_observables * observables;
};

gint32 min(gint32 i, gint32 j);

gdouble gauss(struct VM_simulation* simulation);

gint32 VM_get_box_index(gdouble xposition, gdouble yposition, struct VM_parameters* parameters);

gint32 get_species(guint64 index, struct VM_parameters * parameters);

void mfVM_clean_neighborhood(struct VM_simulation * simulation);

void mfVM_add_to_neighborhood(guint64 index, gdouble distance, gint32 flag, struct VM_simulation * simulation);

gint32 mfVM_get_box_index(gdouble xposition, gdouble yposition, struct VM_parameters* parameters);

void VM_append_particle_to_box(guint64 particle_index, struct VM_simulation * simulation);

void mfVM_append_particle_to_box(guint64 particle_index, struct VM_simulation * simulation);

void VM_fill_boxes(struct VM_simulation * simulation);

void mfVM_fill_boxes(struct VM_simulation * simulation);

void VM_clean_boxes(struct VM_parameters* parameters);

void mfVM_clean_boxes(struct VM_parameters* parameters);

void VM_interaction_with_neighbors_in_box(gint32 box_index, struct VM_simulation * simulation, guint64 index, gdouble * velocity_x, gdouble * velocity_y, gint32 * neighbor_n);

void VM_interaction_with_neighbors_in_box_mirrorx(gint32 box_index, struct VM_simulation * simulation, guint64 index, gdouble * velocity_x, gdouble * velocity_y, gint32 * neighbor_n);

void VM_interaction_with_neighbors_in_box_mirrory(gint32 box_index, struct VM_simulation * simulation, guint64 index, gdouble * velocity_x, gdouble * velocity_y, gint32 * neighbor_n);

void VM_interaction_with_neighbors_in_box_mirrorxy(gint32 box_index, struct VM_simulation * simulation, guint64 index, gdouble * velocity_x, gdouble * velocity_y, gint32 * neighbor_n);

void NVM_interaction_with_neighbors_in_box(gint32 box_index, struct VM_simulation * simulation, guint64 index, gdouble * velocity_x, gdouble * velocity_y, gint32 * neighbor_n);

void NVM_interaction_with_neighbors_in_box_mirrorx(gint32 box_index, struct VM_simulation * simulation, guint64 index, gdouble * velocity_x, gdouble * velocity_y, gint32 * neighbor_n);

void NVM_interaction_with_neighbors_in_box_mirrory(gint32 box_index, struct VM_simulation * simulation, guint64 index, gdouble * velocity_x, gdouble * velocity_y, gint32 * neighbor_n);

void NVM_interaction_with_neighbors_in_box_mirrorxy(gint32 box_index, struct VM_simulation * simulation, guint64 index, gdouble * velocity_x, gdouble * velocity_y, gint32 * neighbor_n);

void L_interaction_with_neighbors_in_box(gint32 box_index, struct VM_simulation * simulation, guint64 index, gdouble * delta_v, gint32 * neighbor_n);

void L_interaction_with_neighbors_in_box_mirrorx(gint32 box_index, struct VM_simulation * simulation, guint64 index, gdouble * delta_v, gint32 * neighbor_n);

void L_interaction_with_neighbors_in_box_mirrory(gint32 box_index, struct VM_simulation * simulation, guint64 index, gdouble * delta_v, gint32 * neighbor_n);

void L_interaction_with_neighbors_in_box_mirrorxy(gint32 box_index, struct VM_simulation * simulation, guint64 index, gdouble * delta_v, gint32 * neighbor_n);

gint32 VM_interaction_with_neighbors(guint64 index, struct VM_simulation * simulation);

gint32 NVM_interaction_with_neighbors(guint64 index, struct VM_simulation * simulation);

gint32 additiveL_interaction_with_neighbors(guint64 index, struct VM_simulation * simulation);

gint32 nonadditiveL_interaction_with_neighbors(guint64 index, struct VM_simulation * simulation);

gint32 mfVM_sort_in_neighbors_in_box(gint32 box_index, struct VM_simulation * simulation, guint64 index);

gint32 mfVM_sort_in_neighbors_in_box_mirrorx(gint32 box_index, struct VM_simulation * simulation, guint64 index);

gint32 mfVM_sort_in_neighbors_in_box_mirrory(gint32 box_index, struct VM_simulation * simulation, guint64 index);

gint32 mfVM_sort_in_neighbors_in_box_mirrorxy(gint32 box_index, struct VM_simulation * simulation, guint64 index);

gint32 mfVM_interaction_with_neighbors(guint64 index, struct VM_simulation * simulation);

gint32 mfNVM_interaction_with_neighbors(guint64 index, struct VM_simulation * simulation);

gint32 mfL_interaction_with_neighbors(guint64 index, struct VM_simulation * simulation);

void VM_one_step(struct VM_simulation * simulation);

void NVM_one_step(struct VM_simulation * simulation);

void mfVM_one_step(struct VM_simulation * simulation);

void mfNVM_one_step(struct VM_simulation * simulation);

void additiveL_one_step(struct VM_simulation * simulation);

void nonadditiveL_one_step(struct VM_simulation * simulation);

void mfL_one_step(struct VM_simulation * simulation);

void VM_one_step_measurement(struct VM_simulation * simulation);

void NVM_one_step_measurement(struct VM_simulation * simulation);

void mfVM_one_step_measurement(struct VM_simulation * simulation);

void mfNVM_one_step_measurement(struct VM_simulation * simulation);

void additiveL_one_step_measurement(struct VM_simulation * simulation);

void nonadditiveL_one_step_measurement(struct VM_simulation * simulation);

void mfL_one_step_measurement(struct VM_simulation * simulation);

gint32 VM_save_state(char * filename, struct VM_simulation * simulation);


#endif
