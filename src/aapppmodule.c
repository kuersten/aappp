/*
 * =====================================================================================
 *
 *       Filename:  aapppmodule.c
 *
 *    Description:  package for agent-based simulations of aligning self-propelled particles in two dimensions 
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
 *
 * =====================================================================================
 */

#define PY_SSIZE_T_CLEAN
#include "c_functions_aappp.h"

//function called from python as aappp_init
static PyObject* VM_init(PyObject* self, PyObject *args, PyObject * keywds)
{
	static char *kwlist[] = {"v", "eta", "R", "omega", "Lx", "Ly", "gamma", "dt", "N", "seed", "kn", "order", "binnum_theta", "binnum_neighbors", "bx", "by", "weight_function", "weight_vector_length", "sf_modes", NULL};
	gdouble R=1.0;
	gdouble Lx=20.0;
	gdouble Ly=20.0;
	gdouble dt=0.01;
	gint32 bx=0;
	gint32 by=0;
	gint32 order=1;
	guint32 seed=1;
	gint32 kn=6;
	gint32 binnum_theta=100;
	gint32 binnum_neighbors=100;
	gint32 k, l;
	PyObject * py_v=NULL;
	gint32 py_v_flag=0;
	PyObject * py_eta=NULL;
	gint32 py_eta_flag=0;
	PyObject * py_omega=NULL;
	gint32 py_omega_flag=0;
	PyObject * py_gamma=NULL;
	gint32 py_gamma_flag=0;
	PyObject * py_N=NULL;
	gint32 py_N_flag=0;
	PyObject * py_weight_function=NULL;
	gint32 weight_vector_length=100;
	PyObject * py_sf_modes=NULL;
	if (!PyArg_ParseTupleAndKeywords(args, keywds, "|OOdOddOdOIiiiiiiOiO", kwlist, &py_v, &py_eta, &R, &py_omega, &Lx, &Ly, &py_gamma, &dt, &py_N, &seed, &kn, &order, &binnum_theta, &binnum_neighbors, &bx, &by, &py_weight_function, &weight_vector_length, &py_sf_modes))
	{
		PyErr_SetString(PyExc_TypeError, "Error when initializing simulation, see >>>help(aappp_init)<<< for help on initialization\n");
		//dereference python objects
		if (py_v_flag==1)
			Py_DECREF(py_v);
		if (py_eta_flag==1)
			Py_DECREF(py_eta);
		if (py_omega_flag==1)
			Py_DECREF(py_omega);
		if (py_gamma_flag==1)
			Py_DECREF(py_gamma);
		if (py_N_flag==1)
			Py_DECREF(py_N);
	        return NULL;
	}
	//if python object optional arguments are not given, set them to their default values
	if (py_v==NULL)
	{
		py_v_flag=1;
		py_v=PyFloat_FromDouble(1.0);
	}
	if (py_eta==NULL)
	{
		py_eta_flag=1;
		py_eta=PyFloat_FromDouble(0.5);
	}
	if (py_omega==NULL)
	{
		py_omega_flag=1;
		py_omega=PyFloat_FromDouble(0.0);
	}
	if (py_gamma==NULL)
	{
		py_gamma_flag=1;
		py_gamma=PyFloat_FromDouble(0.1);
	}
	if (py_N==NULL)
	{
		py_N_flag=1;
		py_N=PyLong_FromLong(1000);
	}
	if (PyList_Check(py_N))
	{
		if (PyList_Size(py_N)==1)
		{
			PyErr_SetString(PyExc_TypeError, "Error when initializing simulation, dimension of N must be either >1 (for multiple species simulation) or N must not be a list (for single species simulation).\n");
			//dereference python objects
			if (py_v_flag==1)
				Py_DECREF(py_v);
			if (py_eta_flag==1)
				Py_DECREF(py_eta);
			if (py_omega_flag==1)
				Py_DECREF(py_omega);
			if (py_gamma_flag==1)
				Py_DECREF(py_gamma);
			if (py_N_flag==1)
				Py_DECREF(py_N);
			return NULL;
		}
		if (PyList_Check(py_v))
			if  ( PyList_Size(py_v) != PyList_Size(py_N))
				{
					PyErr_SetString(PyExc_TypeError, "Dimension of v must be either 1 or equal to dimension of N, see >>>help(aappp_init)<<< for help on initialization.\n");
				//dereference python objects
				if (py_v_flag==1)
					Py_DECREF(py_v);
				if (py_eta_flag==1)
					Py_DECREF(py_eta);
				if (py_omega_flag==1)
					Py_DECREF(py_omega);
				if (py_gamma_flag==1)
					Py_DECREF(py_gamma);
				if (py_N_flag==1)
					Py_DECREF(py_N);
					return NULL;
				}
		if (PyList_Check(py_eta))
			if  ( PyList_Size(py_eta) != PyList_Size(py_N))
				{
					PyErr_SetString(PyExc_TypeError, "Dimension of eta must be either 1 or equal to dimension of N, see >>>help(aappp_init)<<< for help on initialization.\n");
					//dereference python objects
					if (py_v_flag==1)
						Py_DECREF(py_v);
					if (py_eta_flag==1)
						Py_DECREF(py_eta);
					if (py_omega_flag==1)
						Py_DECREF(py_omega);
					if (py_gamma_flag==1)
						Py_DECREF(py_gamma);
					if (py_N_flag==1)
						Py_DECREF(py_N);
					return NULL;
				}
		if (PyList_Check(py_omega))
			if  ( PyList_Size(py_omega) != PyList_Size(py_N))
				{
					PyErr_SetString(PyExc_TypeError, "Dimension of omega must be either 1 or equal to dimension of N, see >>>help(aappp_init)<<< for help on initialization.\n");
					//dereference python objects
					if (py_v_flag==1)
						Py_DECREF(py_v);
					if (py_eta_flag==1)
						Py_DECREF(py_eta);
					if (py_omega_flag==1)
						Py_DECREF(py_omega);
					if (py_gamma_flag==1)
						Py_DECREF(py_gamma);
					if (py_N_flag==1)
						Py_DECREF(py_N);
					return NULL;
				}
		if (PyList_Check(py_gamma))
		{
			if  ( PyList_Size(py_gamma) != PyList_Size(py_N)) 
				{
					PyErr_SetString(PyExc_TypeError, "Dimension of gamma must be either 1 or equal to the number of species, see >>>help(aappp_init)<<< for help on initialization.\n");
					//dereference python objects
					if (py_v_flag==1)
						Py_DECREF(py_v);
					if (py_eta_flag==1)
						Py_DECREF(py_eta);
					if (py_omega_flag==1)
						Py_DECREF(py_omega);
					if (py_gamma_flag==1)
						Py_DECREF(py_gamma);
					if (py_N_flag==1)
						Py_DECREF(py_N);
					return NULL;
				}
			if (!PyList_Check(PyList_GetItem(py_gamma, 0)))
				{
					PyErr_SetString(PyExc_TypeError, "Dimension of gamma[k] must be equal to the number of species, see >>>help(aappp_init)<<< for help on initialization.\n");
					//dereference python objects
					if (py_v_flag==1)
						Py_DECREF(py_v);
					if (py_eta_flag==1)
						Py_DECREF(py_eta);
					if (py_omega_flag==1)
						Py_DECREF(py_omega);
					if (py_gamma_flag==1)
						Py_DECREF(py_gamma);
					if (py_N_flag==1)
						Py_DECREF(py_N);
					return NULL;
				}
			if (PyList_Size(PyList_GetItem(py_gamma, 0))!=PyList_Size(py_N))
				{
					PyErr_SetString(PyExc_TypeError, "Dimension of gamma[k] must be equal to the number of species, see >>>help(aappp_init)<<< for help on initialization.\n");
					//dereference python objects
					if (py_v_flag==1)
						Py_DECREF(py_v);
					if (py_eta_flag==1)
						Py_DECREF(py_eta);
					if (py_omega_flag==1)
						Py_DECREF(py_omega);
					if (py_gamma_flag==1)
						Py_DECREF(py_gamma);
					if (py_N_flag==1)
						Py_DECREF(py_N);
					return NULL;
				}
		}
	}
	struct VM_simulation * simulation=malloc(sizeof(struct VM_simulation));
	simulation->state=malloc(sizeof(struct VM_state));
	simulation->parameters=malloc(sizeof(struct VM_parameters));
	simulation->observables=malloc(sizeof(struct VM_observables));
	simulation->observables->binnum_theta=binnum_theta;
	simulation->observables->binnum_neighbors=binnum_neighbors;
	if (PyList_Check(py_N))
	{
		simulation->parameters->particle_species=PyList_Size(py_N);
		simulation->parameters->species_particle_number=malloc(simulation->parameters->particle_species*sizeof(guint64));
		simulation->parameters->speed=malloc(simulation->parameters->particle_species*sizeof(gdouble));
		simulation->parameters->noise_strength=malloc(simulation->parameters->particle_species*sizeof(gdouble));
		simulation->parameters->omega=malloc(simulation->parameters->particle_species*sizeof(gdouble));
		simulation->parameters->coupling=malloc(simulation->parameters->particle_species*sizeof(gdouble*));
		simulation->state->particle_number=0;
		for (k=0; k<simulation->parameters->particle_species; k++)
		{
			*(simulation->parameters->species_particle_number+k)=PyLong_AsUnsignedLongLong(PyList_GetItem(py_N, k));
			simulation->state->particle_number+=*(simulation->parameters->species_particle_number+k);
			if (PyList_Check(py_v))
				*(simulation->parameters->speed+k)=PyFloat_AsDouble(PyList_GetItem(py_v, k));
			else
				*(simulation->parameters->speed+k)=PyFloat_AsDouble(py_v);
			if (PyList_Check(py_eta))
				*(simulation->parameters->noise_strength+k)=PyFloat_AsDouble(PyList_GetItem(py_eta, k));
			else
				*(simulation->parameters->noise_strength+k)=PyFloat_AsDouble(py_eta);
			if (PyList_Check(py_omega))
				*(simulation->parameters->omega+k)=PyFloat_AsDouble(PyList_GetItem(py_omega, k));
			else
				*(simulation->parameters->omega+k)=PyFloat_AsDouble(py_omega);
			*(simulation->parameters->coupling+k)=malloc(simulation->parameters->particle_species*sizeof(gdouble));
			for (l=0; l<simulation->parameters->particle_species; l++)
			{
				if (PyList_Check(py_gamma))
					*(*(simulation->parameters->coupling+k)+l)=PyFloat_AsDouble(PyList_GetItem(PyList_GetItem(py_gamma, k), l));
				else
					*(*(simulation->parameters->coupling+k)+l)=PyFloat_AsDouble(py_gamma);
			}
		}
		simulation->observables->polar_order=malloc(simulation->parameters->particle_species*sizeof(gdouble));
		for (l=0; l<simulation->parameters->particle_species; l++)
			*(simulation->observables->polar_order+l)=0.0;
		simulation->observables->polar_order_moment1=malloc(simulation->parameters->particle_species*sizeof(gdouble));
		for (l=0; l<simulation->parameters->particle_species; l++)
			*(simulation->observables->polar_order_moment1+l)=0.0;
		simulation->observables->polar_order_moment2=malloc(simulation->parameters->particle_species*sizeof(gdouble));
		for (l=0; l<simulation->parameters->particle_species; l++)
			*(simulation->observables->polar_order_moment2+l)=0.0;
		simulation->observables->polar_order_moment3=malloc(simulation->parameters->particle_species*sizeof(gdouble));
		for (l=0; l<simulation->parameters->particle_species; l++)
			*(simulation->observables->polar_order_moment3+l)=0.0;
		simulation->observables->polar_order_moment4=malloc(simulation->parameters->particle_species*sizeof(gdouble));
		for (l=0; l<simulation->parameters->particle_species; l++)
			*(simulation->observables->polar_order_moment4+l)=0.0;
		simulation->observables->nematic_order=malloc(simulation->parameters->particle_species*sizeof(gdouble));
		for (l=0; l<simulation->parameters->particle_species; l++)
			*(simulation->observables->nematic_order+l)=0.0;
		simulation->observables->nematic_order_moment1=malloc(simulation->parameters->particle_species*sizeof(gdouble));
		for (l=0; l<simulation->parameters->particle_species; l++)
			*(simulation->observables->nematic_order_moment1+l)=0.0;
		simulation->observables->nematic_order_moment2=malloc(simulation->parameters->particle_species*sizeof(gdouble));
		for (l=0; l<simulation->parameters->particle_species; l++)
			*(simulation->observables->nematic_order_moment2+l)=0.0;
		simulation->observables->nematic_order_moment3=malloc(simulation->parameters->particle_species*sizeof(gdouble));
		for (l=0; l<simulation->parameters->particle_species; l++)
			*(simulation->observables->nematic_order_moment3+l)=0.0;
		simulation->observables->nematic_order_moment4=malloc(simulation->parameters->particle_species*sizeof(gdouble));
		for (l=0; l<simulation->parameters->particle_species; l++)
			*(simulation->observables->nematic_order_moment4+l)=0.0;
		simulation->observables->neighbors_moment1=malloc(simulation->parameters->particle_species*sizeof(gdouble));
		for (l=0; l<simulation->parameters->particle_species; l++)
			*(simulation->observables->neighbors_moment1+l)=0.0;
		simulation->observables->neighbors_moment2=malloc(simulation->parameters->particle_species*sizeof(gdouble));
		for (l=0; l<simulation->parameters->particle_species; l++)
			*(simulation->observables->neighbors_moment2+l)=0.0;
		simulation->observables->neighbors_moment3=malloc(simulation->parameters->particle_species*sizeof(gdouble));
		for (l=0; l<simulation->parameters->particle_species; l++)
			*(simulation->observables->neighbors_moment3+l)=0.0;
		simulation->observables->neighbors_moment4=malloc(simulation->parameters->particle_species*sizeof(gdouble));
		for (l=0; l<simulation->parameters->particle_species; l++)
			*(simulation->observables->neighbors_moment4+l)=0.0;
		simulation->observables->hist_theta=malloc(simulation->parameters->particle_species*sizeof(guint64*));
		for (l=0; l<simulation->parameters->particle_species; l++)
		{
			*(simulation->observables->hist_theta+l)=malloc(simulation->observables->binnum_theta*sizeof(guint64));
			for (k=0; k<simulation->observables->binnum_theta; k++)
				*(*(simulation->observables->hist_theta+l)+k)=0;
		}
		simulation->observables->hist_neighbors=malloc(simulation->parameters->particle_species*sizeof(guint64*));
		for (l=0; l<simulation->parameters->particle_species; l++)
		{
			*(simulation->observables->hist_neighbors+l)=malloc(simulation->observables->binnum_neighbors*sizeof(guint64));
			for (k=0; k<simulation->observables->binnum_neighbors; k++)
				*(*(simulation->observables->hist_neighbors+l)+k)=0;
		}
	}
	else
	{
		//for only one species
		simulation->parameters->particle_species=1;
		simulation->parameters->species_particle_number=NULL;
		simulation->state->particle_number=PyLong_AsUnsignedLongLong(py_N);
		simulation->observables->polar_order=malloc(sizeof(gdouble));
		*simulation->observables->polar_order=0.0;
		simulation->observables->polar_order_moment1=malloc(sizeof(gdouble));
		*simulation->observables->polar_order_moment1=0.0;
		simulation->observables->polar_order_moment2=malloc(sizeof(gdouble));
		*simulation->observables->polar_order_moment2=0.0;
		simulation->observables->polar_order_moment3=malloc(sizeof(gdouble));
		*simulation->observables->polar_order_moment3=0.0;
		simulation->observables->polar_order_moment4=malloc(sizeof(gdouble));
		*simulation->observables->polar_order_moment4=0.0;
		simulation->observables->nematic_order=malloc(sizeof(gdouble));
		*simulation->observables->nematic_order=0.0;
		simulation->observables->nematic_order_moment1=malloc(sizeof(gdouble));
		*simulation->observables->nematic_order_moment1=0.0;
		simulation->observables->nematic_order_moment2=malloc(sizeof(gdouble));
		*simulation->observables->nematic_order_moment2=0.0;
		simulation->observables->nematic_order_moment3=malloc(sizeof(gdouble));
		*simulation->observables->nematic_order_moment3=0.0;
		simulation->observables->nematic_order_moment4=malloc(sizeof(gdouble));
		*simulation->observables->nematic_order_moment4=0.0;
		simulation->observables->neighbors_moment1=malloc(sizeof(gdouble));
		*simulation->observables->neighbors_moment1=0.0;
		simulation->observables->neighbors_moment2=malloc(sizeof(gdouble));
		*simulation->observables->neighbors_moment2=0.0;
		simulation->observables->neighbors_moment3=malloc(sizeof(gdouble));
		*simulation->observables->neighbors_moment3=0.0;
		simulation->observables->neighbors_moment4=malloc(sizeof(gdouble));
		*simulation->observables->neighbors_moment4=0.0;
		simulation->observables->hist_theta=malloc(sizeof(guint64*));
		*simulation->observables->hist_theta=malloc(simulation->observables->binnum_theta*sizeof(guint64));
		for (k=0; k<simulation->observables->binnum_theta; k++)
			*(*simulation->observables->hist_theta+k)=0;
		simulation->observables->hist_neighbors=malloc(sizeof(guint64*));
		*simulation->observables->hist_neighbors=malloc(simulation->observables->binnum_neighbors*sizeof(guint64));
		for (k=0; k<simulation->observables->binnum_neighbors; k++)
			*(*simulation->observables->hist_neighbors+k)=0;
		simulation->parameters->speed=malloc(sizeof(gdouble));
		*simulation->parameters->speed=PyFloat_AsDouble(py_v);
		simulation->parameters->noise_strength=malloc(sizeof(gdouble));
		*simulation->parameters->noise_strength=PyFloat_AsDouble(py_eta);
		simulation->parameters->coupling=malloc(sizeof(gdouble *));
		*simulation->parameters->coupling=malloc(sizeof(gdouble));
		**simulation->parameters->coupling=PyFloat_AsDouble(py_gamma);
		simulation->parameters->omega=malloc(sizeof(gdouble));
		*simulation->parameters->omega=PyFloat_AsDouble(py_omega);
	}
	//init structure factor variables
	if (py_sf_modes==NULL)
	{
		simulation->observables->sf_mode_number=0;
		simulation->observables->kx=NULL;
		simulation->observables->ky=NULL;
		simulation->observables->fourier_density=NULL;
		simulation->observables->structure_factor=NULL;
	}
	else
	{
		simulation->observables->sf_mode_number=PyList_Size(py_sf_modes);
		simulation->observables->kx=malloc(simulation->observables->sf_mode_number*sizeof(gint32));
		simulation->observables->ky=malloc(simulation->observables->sf_mode_number*sizeof(gint32));
		simulation->observables->fourier_density=malloc(simulation->parameters->particle_species*sizeof(double complex*));
		simulation->observables->structure_factor=malloc(simulation->parameters->particle_species*sizeof(gdouble*));
		for (k=0; k<simulation->parameters->particle_species; k++)
		{
			*(simulation->observables->fourier_density+k)=malloc(simulation->observables->sf_mode_number*sizeof(double complex));
			*(simulation->observables->structure_factor+k)=malloc(simulation->observables->sf_mode_number*sizeof(gdouble));
		}
		for (k=0; k<simulation->observables->sf_mode_number; k++)
		{
			*(simulation->observables->kx+k)=(gint32) PyLong_AsLong(PyList_GetItem(PyList_GetItem(py_sf_modes, k), 0));
			*(simulation->observables->ky+k)=(gint32) PyLong_AsLong(PyList_GetItem(PyList_GetItem(py_sf_modes, k), 1));
			for (l=0; l<simulation->parameters->particle_species; l++)
			{
				*(*(simulation->observables->fourier_density+l)+k)=0.;
				*(*(simulation->observables->structure_factor+l)+k)=0.;
			}
		}
	}
	//dereference python objects
	if (py_v_flag==1)
		Py_DECREF(py_v);
	if (py_eta_flag==1)
		Py_DECREF(py_eta);
	if (py_omega_flag==1)
		Py_DECREF(py_omega);
	if (py_gamma_flag==1)
		Py_DECREF(py_gamma);
	if (py_N_flag==1)
		Py_DECREF(py_N);
	simulation->parameters->interaction_radius=R;
	simulation->observables->measurement_steps=0;
	simulation->state->length_x=Lx;
	simulation->state->length_y=Ly;
	simulation->state->box_muller=3;
	simulation->state->gauss1=0.0;
	simulation->state->gauss2=0.0;
	simulation->state->particles=malloc(simulation->state->particle_number*sizeof(struct VM_particle));
	simulation->state->prng=g_rand_new_with_seed(seed);
	guint64 i;
	for (i=0; i<simulation->state->particle_number; i++)
	{
		(simulation->state->particles+i)->x=g_rand_double(simulation->state->prng)*simulation->state->length_x;
		(simulation->state->particles+i)->y=g_rand_double(simulation->state->prng)*simulation->state->length_y;
		(simulation->state->particles+i)->theta=g_rand_double(simulation->state->prng)*2*M_PI -M_PI;
		(simulation->state->particles+i)->theta_new=0.0;
		(simulation->state->particles+i)->next_particle=NULL;
		(simulation->state->particles+i)->index=i;
	}
	simulation->parameters->delta_t=dt;
	simulation->parameters->sqrt_delta_t=sqrt(dt);
	simulation->parameters->mf_kn=kn;
	simulation->parameters->interaction_order=order;
	simulation->parameters->boundary_x=bx;
	simulation->parameters->boundary_y=by;
	simulation->state->nh.neighbor_indexes=malloc(simulation->parameters->mf_kn*sizeof(guint64));
	simulation->state->nh.distances=malloc(simulation->parameters->mf_kn*sizeof(gdouble));
	simulation->state->nh.mirror_flag=malloc(simulation->parameters->mf_kn*sizeof(gint32));
	mfVM_clean_neighborhood(simulation);
	simulation->parameters->mf_interaction_radius=sqrt(simulation->state->length_x*simulation->state->length_y*simulation->parameters->mf_kn/M_PI/simulation->state->particle_number);
	k=1;
	while (simulation->state->length_x/(k+1)>simulation->parameters->interaction_radius)
		k++;
	simulation->parameters->box_number_x=k;
	simulation->parameters->box_size_x=simulation->state->length_x/k;
	k=1;
	while (simulation->state->length_y/(k+1)>simulation->parameters->interaction_radius)
		k++;
	simulation->parameters->box_number_y=k;
	simulation->parameters->box_size_y=simulation->state->length_y/k;
	simulation->parameters->box=malloc(simulation->parameters->box_number_x*simulation->parameters->box_number_y*sizeof(struct VM_particle**));
	for (k=0; k<simulation->parameters->box_number_x*simulation->parameters->box_number_y; k++)
	{
		*(simulation->parameters->box+k)=NULL;
	}
	k=1;
	while (simulation->state->length_x/(k+1)>simulation->parameters->mf_interaction_radius)
		k++;
	simulation->parameters->mf_box_number_x=k;
	simulation->parameters->mf_box_size_x=simulation->state->length_x/k;
	k=1;
	while (simulation->state->length_y/(k+1)>simulation->parameters->mf_interaction_radius)
		k++;
	simulation->parameters->mf_box_number_y=k;
	simulation->parameters->mf_box_size_y=simulation->state->length_y/k;
	simulation->parameters->mf_box=malloc(simulation->parameters->mf_box_number_x*simulation->parameters->mf_box_number_y*sizeof(struct VM_particle**));
	for (k=0; k<simulation->parameters->mf_box_number_x*simulation->parameters->mf_box_number_y; k++)
	{
		*(simulation->parameters->mf_box+k)=NULL;
	}
	//set weights for nonadditive Langevin dynamics if weight function is specified
	simulation->state->weight_function_callback=NULL;	//initialize callback variable
	simulation->state->weights=NULL;	//initialize weights variable
	simulation->state->weight_vector_length=0;	//initialize weight vector length
	if (py_weight_function!=NULL)	//set weights
	{
		Py_XDECREF(simulation->state->weight_function_callback);
		Py_XINCREF(py_weight_function);
		simulation->state->weight_function_callback=py_weight_function;
		simulation->state->weight_vector_length=weight_vector_length;
		simulation->state->weights=malloc(simulation->state->weight_vector_length*sizeof(gdouble));
		for (k=0; k<simulation->state->weight_vector_length; k++)
		{
			PyObject * arglist=Py_BuildValue("(i)", k);
			PyObject * result=PyEval_CallObject(simulation->state->weight_function_callback, arglist);
			*(simulation->state->weights+k)=PyFloat_AsDouble(result);
			Py_DECREF(arglist);
			Py_DECREF(result);
		}
		Py_XDECREF(simulation->state->weight_function_callback);	//Python weight function is no longer used
		simulation->state->weight_function_callback=NULL;
	}
	return PyCapsule_New((void *) simulation, "", NULL);
}

//function called from python as aappp_reset_measurement
static PyObject* VM_reset_measurement(PyObject* self, PyObject *args, PyObject * keywds)
{
	PyObject * py_obj;
	static char *kwlist[] = {"simulation", "binnum_theta", "binnum_neighbors", NULL};
	gint32 binnum_theta=100;
	gint32 binnum_neighbors=100;
	gint32 i, j;
	if (!PyArg_ParseTupleAndKeywords(args, keywds, "O|ii", kwlist, &py_obj, &binnum_theta, &binnum_neighbors))
	{
		PyErr_SetString(PyExc_TypeError, "Error when resetting measurement, see >>>help(aappp_reset_measurement)<<< for help.\n");
	        return NULL;
	}
	struct VM_simulation * simulation= (struct VM_simulation *) PyCapsule_GetPointer(py_obj, "");
	simulation->observables->binnum_theta=binnum_theta;
	simulation->observables->binnum_neighbors=binnum_neighbors;
	for (i=0; i<simulation->parameters->particle_species; i++)
	{
		*(simulation->observables->polar_order_moment1+i)=0.0;
		*(simulation->observables->polar_order_moment2+i)=0.0;
		*(simulation->observables->polar_order_moment3+i)=0.0;
		*(simulation->observables->polar_order_moment4+i)=0.0;
		*(simulation->observables->nematic_order_moment1+i)=0.0;
		*(simulation->observables->nematic_order_moment2+i)=0.0;
		*(simulation->observables->nematic_order_moment3+i)=0.0;
		*(simulation->observables->nematic_order_moment4+i)=0.0;
		*(simulation->observables->neighbors_moment1+i)=0.0;
		*(simulation->observables->neighbors_moment2+i)=0.0;
		*(simulation->observables->neighbors_moment3+i)=0.0;
		*(simulation->observables->neighbors_moment4+i)=0.0;
	}
	simulation->observables->measurement_steps=0;
	for (i=0; i<simulation->parameters->particle_species; i++)
		free(*(simulation->observables->hist_theta+i));
	free(simulation->observables->hist_theta);
	simulation->observables->hist_theta=malloc(simulation->parameters->particle_species*sizeof(guint64*));
	for (i=0; i<simulation->parameters->particle_species; i++)
		*(simulation->observables->hist_theta+i)=malloc(simulation->observables->binnum_theta*sizeof(guint64));
	for (i=0; i<simulation->parameters->particle_species; i++)
		for (j=0; j<simulation->observables->binnum_theta; j++)
			*(*(simulation->observables->hist_theta+i)+j)=0;
	for (i=0; i<simulation->parameters->particle_species; i++)
		free(*(simulation->observables->hist_neighbors+i));
	free(simulation->observables->hist_neighbors);
	simulation->observables->hist_neighbors=malloc(simulation->parameters->particle_species*sizeof(guint64*));
	for (i=0; i<simulation->parameters->particle_species; i++)
		*(simulation->observables->hist_neighbors+i)=malloc(simulation->observables->binnum_neighbors*sizeof(guint64));
	for (i=0; i<simulation->parameters->particle_species; i++)
		for (j=0; j<simulation->observables->binnum_neighbors; j++)
			*(*(simulation->observables->hist_neighbors+i)+j)=0;
	return Py_BuildValue("i", 1);
}

//function called from python as aappp_reset_parameters
static PyObject* VM_reset_parameters(PyObject* self, PyObject *args, PyObject * keywds)
{
	PyObject * py_obj;
	static char *kwlist[] = {"simulation", "v", "R", "eta", "omega", "dt", "gamma", "kn", "order", "bx", "by", "N", "weight_function", "weight_vector_length", "sf_modes", NULL};
	gdouble R=1.0;
	gdouble dt=0.01;
	gint32 kn=6;
	gint32 order=1;
	gint32 bx=0;
	gint32 by=0;
	gint32 k, j;
	PyObject * py_v=NULL;
	gint32 py_v_flag=0;
	PyObject * py_eta=NULL;
	gint32 py_eta_flag=0;
	PyObject * py_omega=NULL;
	gint32 py_omega_flag=0;
	PyObject * py_gamma=NULL;
	gint32 py_gamma_flag=0;
	PyObject * py_N=NULL;
	gint32 py_N_flag=0;
	PyObject * py_weight_function=NULL;
	gint32 weight_vector_length=100;
	PyObject * py_sf_modes=NULL;
	//parse arguments for the first time only to get previous parameters from simulation data structure (only first argument relevant)
	if (!PyArg_ParseTupleAndKeywords(args, keywds, "O|OdOOdOiiiiOOiO", kwlist, &py_obj, &py_v, &R, &py_eta, &py_omega, &dt, &py_gamma, &kn, &order, &bx, &by, &py_N, &py_weight_function, &weight_vector_length, &py_sf_modes))
	{
		PyErr_SetString(PyExc_TypeError, "Error when resetting parameters, see >>>help(aappp_reset_parameters)<<< for help.\n");
		return NULL;
	}
	struct VM_simulation * simulation= (struct VM_simulation *) PyCapsule_GetPointer(py_obj, "");
	R=simulation->parameters->interaction_radius;
	dt=simulation->parameters->delta_t;
	kn=simulation->parameters->mf_kn;
	order=simulation->parameters->interaction_order;
	bx=simulation->parameters->boundary_x;
	by=simulation->parameters->boundary_y;
	weight_vector_length=simulation->state->weight_vector_length;
	//now parse once more to update parameters given by the function call
	if (!PyArg_ParseTupleAndKeywords(args, keywds, "O|OdOOdOiiiiOOi", kwlist, &py_obj, &py_v, &R, &py_eta, &py_omega, &dt, &py_gamma, &kn, &order, &bx, &by, &py_N, &py_weight_function, &weight_vector_length))
	{
		PyErr_SetString(PyExc_TypeError, "Error when resetting parameters, see >>>help(aappp_reset_parameters)<<< for help on resetting parameters.\n");
		return NULL;
	}
	//set python objects that have not been parsed to previous values
	if (simulation->parameters->particle_species==1)
	{
		if (py_v==NULL)
		{
			py_v_flag=1;
			py_v=PyFloat_FromDouble(*simulation->parameters->speed);
		}
		if (py_eta==NULL)
		{
			py_eta_flag=1;
			py_eta=PyFloat_FromDouble(*simulation->parameters->noise_strength);
		}
		if (py_omega==NULL)
		{
			py_omega_flag=1;
			py_omega=PyFloat_FromDouble(*simulation->parameters->omega);
		}
		if (py_gamma==NULL)
		{
			py_gamma_flag=1;
			py_gamma=PyFloat_FromDouble(**simulation->parameters->coupling);
		}
		if (py_N==NULL)
		{
			py_N_flag=1;
			py_N=PyLong_FromUnsignedLongLong(simulation->state->particle_number);
		}
	}
	else
	{
		if (py_v==NULL)
		{
			py_v_flag=1;
			py_v=PyList_New(simulation->parameters->particle_species);
			for (k=0; k<simulation->parameters->particle_species; k++)
				PyList_SET_ITEM(py_v, k, PyFloat_FromDouble(*(simulation->parameters->speed+k)));
		}
		if (py_eta==NULL)
		{
			py_eta_flag=1;
			py_eta=PyList_New(simulation->parameters->particle_species);
			for (k=0; k<simulation->parameters->particle_species; k++)
				PyList_SET_ITEM(py_eta, k, PyFloat_FromDouble(*(simulation->parameters->noise_strength+k)));
		}
		if (py_omega==NULL)
		{
			py_omega_flag=1;
			py_omega=PyList_New(simulation->parameters->particle_species);
			for (k=0; k<simulation->parameters->particle_species; k++)
				PyList_SET_ITEM(py_omega, k, PyFloat_FromDouble(*(simulation->parameters->omega+k)));
		}
		if (py_N==NULL)
		{
			py_N_flag=1;
			py_N=PyList_New(simulation->parameters->particle_species);
			for (k=0; k<simulation->parameters->particle_species; k++)
				PyList_SET_ITEM(py_N, k, PyLong_FromUnsignedLongLong(*(simulation->parameters->species_particle_number+k)));
		}
		if (py_gamma==NULL)
		{
			py_gamma_flag=1;
			py_gamma=PyList_New(simulation->parameters->particle_species);
			for (k=0; k<simulation->parameters->particle_species; k++)
			{
				PyObject * py_coup = PyList_New(simulation->parameters->particle_species);
				for (j=0; j<simulation->parameters->particle_species; j++)
					PyList_SET_ITEM(py_coup, j, PyFloat_FromDouble(*(*(simulation->parameters->coupling+k)+j)));
				PyList_SET_ITEM(py_gamma, k, py_coup);
			}
		}
	}
	//if there is only one species of particles
	if (!PyList_Check(py_N))
	{
		//check if total number of particles is unchanged
		if (PyLong_AsUnsignedLongLong(py_N)!=simulation->state->particle_number)
		{
			PyErr_SetString(PyExc_TypeError, "Error when resetting parameters: particle number can not be changed, see >>>help(aappp_reset_parameters)<<< for help on resetting parameters.\n");
			//dereference python objects
			if (py_v_flag==1)
				Py_DECREF(py_v);
			if (py_eta_flag==1)
				Py_DECREF(py_eta);
			if (py_omega_flag==1)
				Py_DECREF(py_omega);
			if (py_gamma_flag==1)
				Py_DECREF(py_gamma);
			if (py_N_flag==1)
				Py_DECREF(py_N);
			return NULL;
		}
		//(re)set velocity
		free(simulation->parameters->speed);
		simulation->parameters->speed=malloc(sizeof(gdouble));
		if (!PyList_Check(py_v))
			*simulation->parameters->speed=PyFloat_AsDouble(py_v);
		else
			*simulation->parameters->speed=PyFloat_AsDouble(PyList_GetItem(py_v, 0));
		//(re)set noise strength
		free(simulation->parameters->noise_strength);
		simulation->parameters->noise_strength=malloc(sizeof(gdouble));
		if (!PyList_Check(py_eta))
			*simulation->parameters->noise_strength=PyFloat_AsDouble(py_eta);
		else
			*simulation->parameters->noise_strength=PyFloat_AsDouble(PyList_GetItem(py_eta, 0));
		//(re)set rotation frequency (omega)
		free(simulation->parameters->omega);
		simulation->parameters->omega=malloc(sizeof(gdouble));
		if (!PyList_Check(py_omega))
			*simulation->parameters->omega=PyFloat_AsDouble(py_omega);
		else
			*simulation->parameters->omega=PyFloat_AsDouble(PyList_GetItem(py_omega, 0));
		//(re)set coupling
		for (k=0; k<simulation->parameters->particle_species; k++)
		{
			free(*(simulation->parameters->coupling+k));
		}
		free(simulation->parameters->coupling);
		simulation->parameters->coupling=malloc(sizeof(gdouble*));
		*simulation->parameters->coupling=malloc(sizeof(gdouble));
		if (!PyList_Check(py_gamma))
			**simulation->parameters->coupling=PyFloat_AsDouble(py_gamma);
		else
			**simulation->parameters->coupling=PyFloat_AsDouble(PyList_GetItem(PyList_GetItem(py_gamma, 0), 0));
		//(re)set number of species
		simulation->parameters->particle_species=1;
		// different particle species not needed
		free(simulation->parameters->species_particle_number);
		simulation->parameters->species_particle_number=NULL;
	}
	else
	{
		guint64 total_N=0;
		for (k=0; k<PyList_Size(py_N); k++)
		{
			total_N+=PyLong_AsUnsignedLongLong(PyList_GetItem(py_N, k));
		}
		//check if total number of particles is unchanged
		if (total_N!=simulation->state->particle_number)
		{
			PyErr_SetString(PyExc_TypeError, "Error when resetting parameters: total particle number can not be changed, see >>>help(aappp_reset_parameters)<<< for help on resetting parameters.\n");
			//dereference python objects
			if (py_v_flag==1)
				Py_DECREF(py_v);
			if (py_eta_flag==1)
				Py_DECREF(py_eta);
			if (py_omega_flag==1)
				Py_DECREF(py_omega);
			if (py_gamma_flag==1)
				Py_DECREF(py_gamma);
			if (py_N_flag==1)
				Py_DECREF(py_N);
			return NULL;
		}
		// set number of species
		simulation->parameters->particle_species=PyList_Size(py_N);
		//(re)set number of particles for different species
		free(simulation->parameters->species_particle_number);
		simulation->parameters->species_particle_number=malloc(simulation->parameters->particle_species*sizeof(gdouble));
		for (k=0; k<simulation->parameters->particle_species; k++)
		{
			*(simulation->parameters->species_particle_number+k)=PyLong_AsUnsignedLongLong(PyList_GetItem(py_N, k));
		}
		//(re)set velocities
		free(simulation->parameters->speed);
		simulation->parameters->speed=malloc(simulation->parameters->particle_species*sizeof(gdouble));
		if (PyList_Check(py_v)) //if velocities are a vector
		{
			if (PyList_Size(py_v)==simulation->parameters->particle_species) // if velocity vector has correct dimension
				for (k=0; k<simulation->parameters->particle_species; k++)
					*(simulation->parameters->speed+k)=PyFloat_AsDouble(PyList_GetItem(py_v, k));
			else //if velocity vector has not correct dimension, use only first entry
				for (k=0; k<simulation->parameters->particle_species; k++)
					*(simulation->parameters->speed+k)=PyFloat_AsDouble(PyList_GetItem(py_v, 0));
		}
		else //velocities are just a single number
			for (k=0; k<simulation->parameters->particle_species; k++)
				*(simulation->parameters->speed+k)=PyFloat_AsDouble(py_v);
		//(re)set noise strength
		free(simulation->parameters->noise_strength);
		simulation->parameters->noise_strength=malloc(simulation->parameters->particle_species*sizeof(gdouble));
		if (PyList_Check(py_eta)) //if noise strengths are a vector
		{
			if (PyList_Size(py_eta)==simulation->parameters->particle_species) // if noise strength vector has correct dimension
				for (k=0; k<simulation->parameters->particle_species; k++)
					*(simulation->parameters->noise_strength+k)=PyFloat_AsDouble(PyList_GetItem(py_eta, k));
			else //if noise strength vector has not correct dimension, use only first entry
				for (k=0; k<simulation->parameters->particle_species; k++)
					*(simulation->parameters->noise_strength+k)=PyFloat_AsDouble(PyList_GetItem(py_eta, 0));
		}
		else //noise strength is just a single number
			for (k=0; k<simulation->parameters->particle_species; k++)
				*(simulation->parameters->noise_strength+k)=PyFloat_AsDouble(py_eta);
		//(re)set rotation frequency (omega)
		free(simulation->parameters->omega);
		simulation->parameters->omega=malloc(simulation->parameters->particle_species*sizeof(gdouble));
		if (PyList_Check(py_omega)) //if omega are a vector
		{
			if (PyList_Size(py_omega)==simulation->parameters->particle_species) // if omega vector has correct dimension
				for (k=0; k<simulation->parameters->particle_species; k++)
					*(simulation->parameters->omega+k)=PyFloat_AsDouble(PyList_GetItem(py_omega, k));
			else //if omega vector has not correct dimension, use only first entry
				for (k=0; k<simulation->parameters->particle_species; k++)
					*(simulation->parameters->omega+k)=PyFloat_AsDouble(PyList_GetItem(py_omega, 0));
		}
		else //if omega is just a single number
			for (k=0; k<simulation->parameters->particle_species; k++)
				*(simulation->parameters->omega+k)=PyFloat_AsDouble(py_omega);
		//(re)set coupling
		for (k=0; k<simulation->parameters->particle_species; k++)
		{
			free(*(simulation->parameters->coupling+k));
		}
		free(simulation->parameters->coupling);
		simulation->parameters->coupling=malloc(simulation->parameters->particle_species*sizeof(gdouble*));
		for (k=0; k<simulation->parameters->particle_species; k++)
			*(simulation->parameters->coupling+k)=malloc(simulation->parameters->particle_species*sizeof(gdouble));
		if (PyList_Check(py_gamma)) //if gamma (coupling) is a vector
		{
			if (PyList_Size(py_gamma)==simulation->parameters->particle_species) // if gamma has correct dimension
				for (k=0; k<simulation->parameters->particle_species; k++)
					for (j=0; j<simulation->parameters->particle_species; j++)
						*(*(simulation->parameters->coupling+k)+j)=PyFloat_AsDouble(PyList_GetItem(PyList_GetItem(py_gamma, k), j));
			else //if gamma vector has not correct dimension, use only first entry
				for (k=0; k<simulation->parameters->particle_species; k++)
					for (j=0; j<simulation->parameters->particle_species; j++)
						*(*(simulation->parameters->coupling+k)+j)=PyFloat_AsDouble(PyList_GetItem(PyList_GetItem(py_gamma, 0), 0));
		}
		else //if gamma is just a single number
			for (k=0; k<simulation->parameters->particle_species; k++)
				for (j=0; j<simulation->parameters->particle_species; j++)
					*(*(simulation->parameters->coupling+k)+j)=PyFloat_AsDouble(py_gamma);
	}
	//reset structure factor variables
	free(simulation->observables->kx);
	free(simulation->observables->ky);
	for (k=0; k<simulation->parameters->particle_species; k++)
	{
		free(*(simulation->observables->fourier_density+k));
		free(*(simulation->observables->structure_factor+k));
	}
	free(simulation->observables->fourier_density);
	free(simulation->observables->structure_factor);
	if (py_sf_modes==NULL)
	{
		simulation->observables->sf_mode_number=0;
		simulation->observables->kx=NULL;
		simulation->observables->ky=NULL;
		simulation->observables->fourier_density=NULL;
		simulation->observables->structure_factor=NULL;
	}
	else
	{
		simulation->observables->sf_mode_number=PyList_Size(py_sf_modes);
		simulation->observables->kx=malloc(simulation->observables->sf_mode_number*sizeof(gint32));
		simulation->observables->ky=malloc(simulation->observables->sf_mode_number*sizeof(gint32));
		simulation->observables->fourier_density=malloc(simulation->parameters->particle_species*sizeof(double complex*));
		simulation->observables->structure_factor=malloc(simulation->parameters->particle_species*sizeof(gdouble*));
		for (k=0; k<simulation->parameters->particle_species; k++)
		{
			*(simulation->observables->fourier_density+k)=malloc(simulation->observables->sf_mode_number*sizeof(double complex));
			*(simulation->observables->structure_factor+k)=malloc(simulation->observables->sf_mode_number*sizeof(gdouble));
		}
		for (k=0; k<simulation->observables->sf_mode_number; k++)
		{
			*(simulation->observables->kx+k)=(gint32) PyLong_AsLong(PyList_GetItem(PyList_GetItem(py_sf_modes, k), 0));
			*(simulation->observables->ky+k)=(gint32) PyLong_AsLong(PyList_GetItem(PyList_GetItem(py_sf_modes, k), 1));
			for (j=0; j<simulation->parameters->particle_species; j++)
			{
				*(*(simulation->observables->fourier_density+j)+k)=0.;
				*(*(simulation->observables->structure_factor+j)+k)=0.;
			}
		}
	}
	//dereference python objects
	if (py_v_flag==1)
		Py_DECREF(py_v);
	if (py_eta_flag==1)
		Py_DECREF(py_eta);
	if (py_omega_flag==1)
		Py_DECREF(py_omega);
	if (py_gamma_flag==1)
		Py_DECREF(py_gamma);
	if (py_N_flag==1)
		Py_DECREF(py_N);
	//(re)set other parameters
	simulation->parameters->interaction_radius=R;
	simulation->parameters->delta_t=dt;
	simulation->parameters->sqrt_delta_t=sqrt(dt);
	simulation->parameters->mf_kn=kn;
	simulation->parameters->mf_interaction_radius=sqrt(simulation->state->length_x*simulation->state->length_y*simulation->parameters->mf_kn/M_PI/simulation->state->particle_number);
	simulation->parameters->interaction_order=order;
	simulation->parameters->boundary_x=bx;
	simulation->parameters->boundary_y=by;
	free(simulation->parameters->box);
	free(simulation->parameters->mf_box);
	k=1;
	while (simulation->state->length_x/(k+1)>simulation->parameters->interaction_radius)
		k++;
	simulation->parameters->box_number_x=k;
	simulation->parameters->box_size_x=simulation->state->length_x/k;
	k=1;
	while (simulation->state->length_y/(k+1)>simulation->parameters->interaction_radius)
		k++;
	simulation->parameters->box_number_y=k;
	simulation->parameters->box_size_y=simulation->state->length_y/k;
	simulation->parameters->box=malloc(simulation->parameters->box_number_x*simulation->parameters->box_number_y*sizeof(struct VM_particle**));
	for (k=0; k<simulation->parameters->box_number_x*simulation->parameters->box_number_y; k++)
	{
		*(simulation->parameters->box+k)=NULL;
	}
	k=1;
	while (simulation->state->length_x/(k+1)>simulation->parameters->mf_interaction_radius)
		k++;
	simulation->parameters->mf_box_number_x=k;
	simulation->parameters->mf_box_size_x=simulation->state->length_x/k;
	k=1;
	while (simulation->state->length_y/(k+1)>simulation->parameters->mf_interaction_radius)
		k++;
	simulation->parameters->mf_box_number_y=k;
	simulation->parameters->mf_box_size_y=simulation->state->length_y/k;
	simulation->parameters->mf_box=malloc(simulation->parameters->mf_box_number_x*simulation->parameters->mf_box_number_y*sizeof(struct VM_particle**));
	for (k=0; k<simulation->parameters->mf_box_number_x*simulation->parameters->mf_box_number_y; k++)
	{
		*(simulation->parameters->mf_box+k)=NULL;
	}
	free(simulation->state->nh.neighbor_indexes);
	simulation->state->nh.neighbor_indexes=malloc(simulation->parameters->mf_kn*sizeof(guint64));
	free(simulation->state->nh.distances);
	simulation->state->nh.distances=malloc(simulation->parameters->mf_kn*sizeof(gdouble));
	free(simulation->state->nh.mirror_flag);
	simulation->state->nh.mirror_flag=malloc(simulation->parameters->mf_kn*sizeof(gint32));
	mfVM_clean_neighborhood(simulation);
	//reset weights for nonadditive Langevin dynamics if weight function is specified
	if (py_weight_function!=NULL)	//reset weights
	{
		Py_XDECREF(simulation->state->weight_function_callback);
		Py_XINCREF(py_weight_function);
		simulation->state->weight_function_callback=py_weight_function;
		simulation->state->weight_vector_length=weight_vector_length;
		free(simulation->state->weights);
		simulation->state->weights=malloc(simulation->state->weight_vector_length*sizeof(gdouble));
		for (k=0; k<simulation->state->weight_vector_length; k++)
		{
			PyObject * arglist=Py_BuildValue("(i)", k);
			PyObject * result=PyEval_CallObject(simulation->state->weight_function_callback, arglist);
			*(simulation->state->weights+k)=PyFloat_AsDouble(result);
			Py_DECREF(arglist);
			Py_DECREF(result);
		}
		Py_XDECREF(simulation->state->weight_function_callback);	//Python weight function is no longer used
		simulation->state->weight_function_callback=NULL;
	}
	return Py_BuildValue("i", 1);
}

//function called from python as aappp_load
static PyObject* VM_load(PyObject* self, PyObject *args)
{
	char * filename;
	if (!PyArg_ParseTuple(args, "s", &filename))
	{
		PyErr_SetString(PyExc_TypeError, "Error when loading simulation state from file.\n");
		return NULL;
	}
	gint32 myfile;
	guint64 i;
	gint32 j, k, v;
	if ((myfile=open(filename, O_RDONLY))==-1)
	{
		PyErr_SetString(PyExc_TypeError, "Error when opening file for reading simulation state.\n");
		return NULL;
	}
	//check version flag in data file, this software can only handle datafiles with version flag=0 (version 1.0) or version flag=1 (version 1.1), or version flag=2 (version 1.2) or version flag=999 (version 1.2.1, side branch for structure factor)
	v=999;
	if ((read(myfile, &v, sizeof(gint32)))==-1)
	{
		PyErr_SetString(PyExc_TypeError, "Error when reading version flag from data file\n");
		return NULL;
	}
	//modified in version 1.2
	if ((v!=0) && (v!=1) && (v!=2) && (v!=999))
	{
		PyErr_SetString(PyExc_TypeError, "This is aappp version 1.2. Apparently the data file was produced by a later version and can not be read by this version.\n");
		return NULL;
	}
	//new in version 1.1
	//modified in version 1.2
	//print warning if data was produced by version 1.0
	if (v==0)
		PyErr_WarnEx(PyExc_Warning, "This is aappp version 1.2, you are loading a data file from aappp version 1.0, from v1.0 to v1.1 the measuring procedure for second to fourth moments of the number of neighbors was changed, be carefull when using those results, all other measured quantities can be used without harm, see >>>help(aappp)<<< for version information. WARNING: In versions 1.0 and 1.1 reflecting boundary conditions for the models mfL, additiveL and nonadditiveL were not correctly implemented. Periodic bc for all models, and reflecting boundary conditions for models VM, NVM, mfVM and mfNVM were ok.", 1);
	//new in version 1.2
	//print warning if data was produced by version 1.1
	if (v==1)
		PyErr_WarnEx(PyExc_Warning, "This is aappp version 1.2, you are loading a data file from aappp version 1.1. WARNING: In version 1.1 reflecting boundary conditions for the models mfL, additiveL and nonadditiveL were not correctly implemented. Periodic bc for all models, and reflecting boundary conditions for models VM, NVM, mfVM and mfNVM were ok.", 1);
	//build simulation data structure
	struct VM_simulation * simulation=malloc(sizeof(struct VM_simulation));
	//build data structure to hold state variables (including pseudo random number generator state)
	simulation->state=malloc(sizeof(struct VM_state));
	//build data structure holding parameters
	simulation->parameters=malloc(sizeof(struct VM_parameters));
	//build data structure holding observables
	simulation->observables=malloc(sizeof(struct VM_observables));
	//read in parameters
		//read particle number
	if ((read(myfile, &(simulation->state->particle_number), sizeof(guint64)))==-1)
	{
		PyErr_SetString(PyExc_TypeError, "Error when reading particle number\n");
		return NULL;
	}
		//read particles species number
	if ((read(myfile, &(simulation->parameters->particle_species), sizeof(gint32)))==-1)
	{
		PyErr_SetString(PyExc_TypeError, "Error when reading number of particle species\n");
		return NULL;
	}
		//if there is more than one species, read particle numbers for all species
	if (simulation->parameters->particle_species > 1)
	{
		//build datastructure to hold the number of particles for different species
		simulation->parameters->species_particle_number=malloc(simulation->parameters->particle_species*sizeof(guint64));
		for(k=0; k<simulation->parameters->particle_species; k++)
		{
			if ((read(myfile, simulation->parameters->species_particle_number+k, sizeof(gint32)))==-1)
			{
				PyErr_SetString(PyExc_TypeError, "Error when reading number of particles in each species\n");
				return NULL;
			}
		}
	}
	else
		simulation->parameters->species_particle_number=NULL; //we do not need this pointer
		//read interaction radius
	if ((read(myfile, &(simulation->parameters->interaction_radius), sizeof(gdouble)))==-1)
	{
		PyErr_SetString(PyExc_TypeError, "Error when reading interaction radius\n");
		return NULL;
	}
		//read step size
	if ((read(myfile, &(simulation->parameters->delta_t), sizeof(gdouble)))==-1)
	{
		PyErr_SetString(PyExc_TypeError, "Error when reading step size\n");
		return NULL;
	}
		//set sqrt(step size)
	simulation->parameters->sqrt_delta_t=sqrt(simulation->parameters->delta_t);
		//read metric free neighbor number
	if ((read(myfile, &(simulation->parameters->mf_kn), sizeof(gint32)))==-1)
	{
		PyErr_SetString(PyExc_TypeError, "Error when reading metric free neighbor number\n");
		return NULL;
	}
		//read interaction order
	if ((read(myfile, &(simulation->parameters->interaction_order), sizeof(gint32)))==-1)
	{
		PyErr_SetString(PyExc_TypeError, "Error when reading interaction order\n");
		return NULL;
	}
		//read boundary flag for x-direction
	if ((read(myfile, &(simulation->parameters->boundary_x), sizeof(gint32)))==-1)
	{
		PyErr_SetString(PyExc_TypeError, "Error when reading boundary flag for x-direction\n");
		return NULL;
	}
		//read boundary flag for y-direction
	if ((read(myfile, &(simulation->parameters->boundary_y), sizeof(gint32)))==-1)
	{
		PyErr_SetString(PyExc_TypeError, "Error when reading boundary flag for y-direction\n");
		return NULL;
	}
		//allocate memory to hold omega(s)
	simulation->parameters->omega=malloc(simulation->parameters->particle_species*sizeof(gdouble));
		//read natural frequency(ies) for chiral systems omega(s)
	for(k=0; k<simulation->parameters->particle_species; k++)
	{
		if ((read(myfile, simulation->parameters->omega+k, sizeof(gdouble)))==-1)
		{
			PyErr_SetString(PyExc_TypeError, "Error when reading omega(s)\n");
			return NULL;
		}
	}
		//allocate memory to hold eta(s)
	simulation->parameters->noise_strength=malloc(simulation->parameters->particle_species*sizeof(gdouble));
		//read noise strength(s)
	for(k=0; k<simulation->parameters->particle_species; k++)
	{
		if ((read(myfile, simulation->parameters->noise_strength+k, sizeof(gdouble)))==-1)
		{
			PyErr_SetString(PyExc_TypeError, "Error when reading eta(s)\n");
			return NULL;
		}
	}
		//allocate memory to hold velocity(ies)
	simulation->parameters->speed=malloc(simulation->parameters->particle_species*sizeof(gdouble));
		//read velocity(ies)
	for(k=0; k<simulation->parameters->particle_species; k++)
	{
		if ((read(myfile, simulation->parameters->speed+k, sizeof(gdouble)))==-1)
		{
			PyErr_SetString(PyExc_TypeError, "Error when reading v(s)\n");
			return NULL;
		}
	}
		//allocate memory to save coupling matrix
	simulation->parameters->coupling=malloc(simulation->parameters->particle_species*sizeof(gdouble*));
	for (j=0; j<simulation->parameters->particle_species; j++)
		*(simulation->parameters->coupling+j)=malloc(simulation->parameters->particle_species*sizeof(gdouble));
		//read coupling(s)
	for (j=0; j<simulation->parameters->particle_species; j++)
	{
		for (k=0; k<simulation->parameters->particle_species; k++)
		{
			if ((read(myfile, *(simulation->parameters->coupling+j)+k, sizeof(gdouble)))==-1)
			{
				PyErr_SetString(PyExc_TypeError, "Error when reading coupling matrix\n");
				return NULL;
			}
		}
	}
		//read simulation box size in x-direction Lx
	if ((read(myfile, &(simulation->state->length_x), sizeof(gdouble)))==-1)
	{
		PyErr_SetString(PyExc_TypeError, "Error when reading Lx\n");
		return NULL;
	}
		//read simulation box size in y-direction Ly
	if ((read(myfile, &(simulation->state->length_y), sizeof(gdouble)))==-1)
	{
		PyErr_SetString(PyExc_TypeError, "Error when reading Ly\n");
		return NULL;
	}
	//read internal states related to random number generation and weight function
		//read gauss1
	if ((read(myfile, &(simulation->state->gauss1), sizeof(gdouble)))==-1)
	{
		PyErr_SetString(PyExc_TypeError, "Error when reading gauss1\n");
		return NULL;
	}
		//read gauss2
	if ((read(myfile, &(simulation->state->gauss2), sizeof(gdouble)))==-1)
	{
		PyErr_SetString(PyExc_TypeError, "Error when reading gauss2\n");
		return NULL;
	}
		//read box_muller flag
	if ((read(myfile, &(simulation->state->box_muller), sizeof(gint32)))==-1)
	{
		PyErr_SetString(PyExc_TypeError, "Error when reading box muller flag\n");
		return NULL;
	}
		//allocate pseuda random number state variable 
	simulation->state->prng=malloc(624*sizeof(guint32)+sizeof(guint));
		//read pseuda random number state variable
	if ((read(myfile, simulation->state->prng, 624*sizeof(guint32)+sizeof(guint)))==-1)
	{
		PyErr_SetString(PyExc_TypeError, "Error when reading pseudo random number generator state\n");
		return NULL;
	}
		//read weight vector length
	if ((read(myfile, &(simulation->state->weight_vector_length), sizeof(gint32)))==-1)
	{
		PyErr_SetString(PyExc_TypeError, "Error when reading weight vector length\n");
		return NULL;
	}
		//read weights (if there are some)
	simulation->state->weights=NULL;
	if (simulation->state->weight_vector_length>0)
	{
		//allocate memory for weights
		simulation->state->weights=malloc(simulation->state->weight_vector_length*sizeof(gdouble));
		if ((read(myfile, simulation->state->weights, simulation->state->weight_vector_length*sizeof(gdouble)))==-1)
		{
			PyErr_SetString(PyExc_TypeError, "Error when reading weight vector\n");
			return NULL;
		}
	}
		//set weight function caller reference to NULL
	simulation->state->weight_function_callback=NULL;
		//allocate memory for particles
	simulation->state->particles=malloc(simulation->state->particle_number*sizeof(struct VM_particle));
		//read particles positions and set particles indexes
		//read x-positions
	for(i=0; i<simulation->state->particle_number; i++)
	{
		(simulation->state->particles+i)->index=i;
		(simulation->state->particles+i)->next_particle=NULL;
		if ((read(myfile, &((simulation->state->particles+i)->x), sizeof(gdouble)))==-1)
		{
			PyErr_SetString(PyExc_TypeError, "Error when reading x-positions\n");
			return NULL;
		}
	}
		//read y-positions
	for(i=0; i<simulation->state->particle_number; i++)
	{
		if ((read(myfile, &((simulation->state->particles+i)->y), sizeof(gdouble)))==-1)
		{
			PyErr_SetString(PyExc_TypeError, "Error when reading y-positions\n");
			return NULL;
		}
	}
		//read orientations
	for(i=0; i<simulation->state->particle_number; i++)
	{
		if ((read(myfile, &((simulation->state->particles+i)->theta), sizeof(gdouble)))==-1)
		{
			PyErr_SetString(PyExc_TypeError, "Error when reading orientations\n");
			return NULL;
		}
	}
	//read measurement results
		//allocate memory for polar order measurement
	simulation->observables->polar_order=malloc(simulation->parameters->particle_species*sizeof(gdouble));
		//read polar order(s)
	for (k=0; k<simulation->parameters->particle_species; k++)
	{
		if ((read(myfile, (simulation->observables->polar_order+k), sizeof(gdouble)))==-1)
		{
			PyErr_SetString(PyExc_TypeError, "Error when reading polar order(s)\n");
			return NULL;
		}
	}
		//allocate memory for polar order first moment
	simulation->observables->polar_order_moment1=malloc(simulation->parameters->particle_species*sizeof(gdouble));
		//read first moment of polar order(s)
	for (k=0; k<simulation->parameters->particle_species; k++)
	{
		if ((read(myfile, (simulation->observables->polar_order_moment1+k), sizeof(gdouble)))==-1)
		{
			PyErr_SetString(PyExc_TypeError, "Error when reading polar order first moment(s)\n");
			return NULL;
		}
	}
		//allocate memory for polar order second moment
	simulation->observables->polar_order_moment2=malloc(simulation->parameters->particle_species*sizeof(gdouble));
		//read second moment of polar order(s)
	for (k=0; k<simulation->parameters->particle_species; k++)
	{
		if ((read(myfile, (simulation->observables->polar_order_moment2+k), sizeof(gdouble)))==-1)
		{
			PyErr_SetString(PyExc_TypeError, "Error when reading polar order second moment(s)\n");
			return NULL;
		}
	}
		//allocate memory for polar order third moment
	simulation->observables->polar_order_moment3=malloc(simulation->parameters->particle_species*sizeof(gdouble));
		//read third moment of polar order(s)
	for (k=0; k<simulation->parameters->particle_species; k++)
	{
		if ((read(myfile, (simulation->observables->polar_order_moment3+k), sizeof(gdouble)))==-1)
		{
			PyErr_SetString(PyExc_TypeError, "Error when reading polar order third moment(s)\n");
			return NULL;
		}
	}
		//allocate memory for polar order fourth moment
	simulation->observables->polar_order_moment4=malloc(simulation->parameters->particle_species*sizeof(gdouble));
		//read fourth moment of polar order(s)
	for (k=0; k<simulation->parameters->particle_species; k++)
	{
		if ((read(myfile, (simulation->observables->polar_order_moment4+k), sizeof(gdouble)))==-1)
		{
			PyErr_SetString(PyExc_TypeError, "Error when reading polar oder fourth moment(s)\n");
			return NULL;
		}
	}
		//allocate memory for nematic order measurement
	simulation->observables->nematic_order=malloc(simulation->parameters->particle_species*sizeof(gdouble));
		//read nematic order(s)
	for (k=0; k<simulation->parameters->particle_species; k++)
	{
		if ((read(myfile, (simulation->observables->nematic_order+k), sizeof(gdouble)))==-1)
		{
			PyErr_SetString(PyExc_TypeError, "Error when reading nematic order(s)\n");
			return NULL;
		}
	}
		//allocate memory for nematic order first moment
	simulation->observables->nematic_order_moment1=malloc(simulation->parameters->particle_species*sizeof(gdouble));
		//read first moment of nematic order(s)
	for (k=0; k<simulation->parameters->particle_species; k++)
	{
		if ((read(myfile, (simulation->observables->nematic_order_moment1+k), sizeof(gdouble)))==-1)
		{
			PyErr_SetString(PyExc_TypeError, "Error when reading nematic order first moment(s)\n");
			return NULL;
		}
	}
		//allocate memory for nematic order second moment
	simulation->observables->nematic_order_moment2=malloc(simulation->parameters->particle_species*sizeof(gdouble));
		//read secondmoment of nematic order(s)
	for (k=0; k<simulation->parameters->particle_species; k++)
	{
		if ((read(myfile, (simulation->observables->nematic_order_moment2+k), sizeof(gdouble)))==-1)
		{
			PyErr_SetString(PyExc_TypeError, "Error when reading nematic order second moment(s)\n");
			return NULL;
		}
	}
		//allocate memory for nematic order third moment
	simulation->observables->nematic_order_moment3=malloc(simulation->parameters->particle_species*sizeof(gdouble));
		//read third moment of nematic order(s)
	for (k=0; k<simulation->parameters->particle_species; k++)
	{
		if ((read(myfile, (simulation->observables->nematic_order_moment3+k), sizeof(gdouble)))==-1)
		{
			PyErr_SetString(PyExc_TypeError, "Error when reading nematic order third moment(s)\n");
			return NULL;
		}
	}
		//allocate memory for nematic order fourth moment
	simulation->observables->nematic_order_moment4=malloc(simulation->parameters->particle_species*sizeof(gdouble));
		//read fourth moment of nematic order(s)
	for (k=0; k<simulation->parameters->particle_species; k++)
	{
		if ((read(myfile, (simulation->observables->nematic_order_moment4+k), sizeof(gdouble)))==-1)
		{
			PyErr_SetString(PyExc_TypeError, "Error when reading nematic order fourth moment(s)\n");
			return NULL;
		}
	}
		//allocate memory for neighbor number first moment
	simulation->observables->neighbors_moment1=malloc(simulation->parameters->particle_species*sizeof(gdouble));
		//read first moment of neighbor number(s)
	for (k=0; k<simulation->parameters->particle_species; k++)
	{
		if ((read(myfile, (simulation->observables->neighbors_moment1+k), sizeof(gdouble)))==-1)
		{
			PyErr_SetString(PyExc_TypeError, "Error when reading neighbors first moment(s)\n");
			return NULL;
		}
	}
		//allocate memory for neighbor number second moment
	simulation->observables->neighbors_moment2=malloc(simulation->parameters->particle_species*sizeof(gdouble));
		//read second moment of neighbor number(s)
	for (k=0; k<simulation->parameters->particle_species; k++)
	{
		if ((read(myfile, (simulation->observables->neighbors_moment2+k), sizeof(gdouble)))==-1)
		{
			PyErr_SetString(PyExc_TypeError, "Error when reading neighbors second moment(s)\n");
			return NULL;
		}
	}
		//allocate memory for neighbor number third moment
	simulation->observables->neighbors_moment3=malloc(simulation->parameters->particle_species*sizeof(gdouble));
		//read third moment of neighbor number(s)
	for (k=0; k<simulation->parameters->particle_species; k++)
	{
		if ((read(myfile, (simulation->observables->neighbors_moment3+k), sizeof(gdouble)))==-1)
		{
			PyErr_SetString(PyExc_TypeError, "Error when reading neighbors third moment(s)\n");
			return NULL;
		}
	}
		//allocate memory for neighbor number fourth moment
	simulation->observables->neighbors_moment4=malloc(simulation->parameters->particle_species*sizeof(gdouble));
		//read fourth moment of neighbor number(s)
	for (k=0; k<simulation->parameters->particle_species; k++)
	{
		if ((read(myfile, (simulation->observables->neighbors_moment4+k), sizeof(gdouble)))==-1)
		{
			PyErr_SetString(PyExc_TypeError, "Error when reading neighbors fourth moment(s)\n");
			return NULL;
		}
	}
		//read number of measurement time steps
	if ((read(myfile, &(simulation->observables->measurement_steps), sizeof(guint64)))==-1)
	{
		PyErr_SetString(PyExc_TypeError, "Error when reading number of measurement steps\n");
		return NULL;
	}
		//read number of histogram bins for orientation
	if ((read(myfile, &(simulation->observables->binnum_theta), sizeof(gint32)))==-1)
	{
		PyErr_SetString(PyExc_TypeError, "Error when reading number histogram bins for orientation\n");
		return NULL;
	}
		//read number of histogram bins for neighbor number
	if ((read(myfile, &(simulation->observables->binnum_neighbors), sizeof(gint32)))==-1)
	{
		PyErr_SetString(PyExc_TypeError, "Error when reading number histogram bins for neighbor number\n");
		return NULL;
	}
		//allocate memory for orientation histograms
	simulation->observables->hist_theta=malloc(simulation->parameters->particle_species*sizeof(guint64*));
	for (k=0; k<simulation->parameters->particle_species; k++)
		*(simulation->observables->hist_theta+k)=malloc(simulation->observables->binnum_theta*sizeof(guint64));
		//read histogram(s) for orientations
	for (k=0; k<simulation->parameters->particle_species; k++)
		for (j=0; j<simulation->observables->binnum_theta; j++)
			if ((read(myfile, (*(simulation->observables->hist_theta+k)+j), sizeof(guint64)))==-1)
			{
				PyErr_SetString(PyExc_TypeError, "Error when reading orientation histogram\n");
				return NULL;
			}
		//allocate memory for neighbor histograms
	simulation->observables->hist_neighbors=malloc(simulation->parameters->particle_species*sizeof(guint64*));
	for (k=0; k<simulation->parameters->particle_species; k++)
		*(simulation->observables->hist_neighbors+k)=malloc(simulation->observables->binnum_neighbors*sizeof(guint64));
		//read histogram(s) for neighbor numbers
	for (k=0; k<simulation->parameters->particle_species; k++)
		for (j=0; j<simulation->observables->binnum_neighbors; j++)
			if ((read(myfile, (*(simulation->observables->hist_neighbors+k)+j), sizeof(guint64)))==-1)
			{
				PyErr_SetString(PyExc_TypeError, "Error when reading neighbor histogram\n");
				return NULL;
			}
	//read/init structure factor variables
	if (v!=999)
	{
		simulation->observables->sf_mode_number=0;
		simulation->observables->kx=NULL;
		simulation->observables->ky=NULL;
		simulation->observables->fourier_density=NULL;
		simulation->observables->structure_factor=NULL;
	}
	else
	{
		if ((read(myfile, &(simulation->observables->sf_mode_number), sizeof(gint32)))==-1)
		{
			PyErr_SetString(PyExc_TypeError, "Error when reading structure factor mode number\n");
			return NULL;
		}
		simulation->observables->kx=malloc(simulation->observables->sf_mode_number*sizeof(gint32));
		simulation->observables->ky=malloc(simulation->observables->sf_mode_number*sizeof(gint32));
		for (k=0; k<simulation->observables->sf_mode_number; k++)
		{
			if ((read(myfile, (simulation->observables->kx+k), sizeof(gint32)))==-1)
			{
				PyErr_SetString(PyExc_TypeError, "Error when reading structure x-modes\n");
				return NULL;
			}
			if ((read(myfile, (simulation->observables->ky+k), sizeof(gint32)))==-1)
			{
				PyErr_SetString(PyExc_TypeError, "Error when reading structure y-modes\n");
				return NULL;
			}
		}
		simulation->observables->fourier_density=malloc(simulation->parameters->particle_species*sizeof(double complex*));
		simulation->observables->structure_factor=malloc(simulation->parameters->particle_species*sizeof(gdouble*));
		for (k=0; k<simulation->parameters->particle_species; k++)
		{
			*(simulation->observables->fourier_density+k)=malloc(simulation->observables->sf_mode_number*sizeof(double complex));
			*(simulation->observables->structure_factor+k)=malloc(simulation->observables->sf_mode_number*sizeof(gdouble));
			for (j=0; j<simulation->observables->sf_mode_number; j++)
			{
				*(*(simulation->observables->fourier_density+k)+j)=0.0;
				if ((read(myfile, (*(simulation->observables->structure_factor+k)+j), sizeof(gdouble)))==-1)
				{
					PyErr_SetString(PyExc_TypeError, "Error when reading structure factor\n");
					return NULL;
				}
			}
		}
	}
	//reading file is done
	close(myfile);
	//calculate metrix free box interaction radius
	simulation->parameters->mf_interaction_radius=sqrt(simulation->state->length_x*simulation->state->length_y*simulation->parameters->mf_kn/M_PI/simulation->state->particle_number);
	//calculate box number in x-direction
	k=1;
	while (simulation->state->length_x/(k+1)>simulation->parameters->interaction_radius)
		k++;
	simulation->parameters->box_number_x=k;
	//calculate box size in x-direction
	simulation->parameters->box_size_x=simulation->state->length_x/k;
	//calculate box number in y-direction
	k=1;
	while (simulation->state->length_y/(k+1)>simulation->parameters->interaction_radius)
		k++;
	simulation->parameters->box_number_y=k;
	//calculate box size in y-direction
	simulation->parameters->box_size_y=simulation->state->length_y/k;
	//now for metric free models
	//box number in x-direction (mf)
	k=1;
	while (simulation->state->length_x/(k+1)>simulation->parameters->mf_interaction_radius)
		k++;
	simulation->parameters->mf_box_number_x=k;
	//box size in x-direction (mf)
	simulation->parameters->mf_box_size_x=simulation->state->length_x/k;
	//box number in y-direction (mf)
	k=1;
	while (simulation->state->length_y/(k+1)>simulation->parameters->mf_interaction_radius)
		k++;
	simulation->parameters->mf_box_number_y=k;
	//box size in y-direction (mf)
	simulation->parameters->mf_box_size_y=simulation->state->length_y/k;
	//create box data structure
	simulation->parameters->box=malloc(simulation->parameters->box_number_x*simulation->parameters->box_number_y*sizeof(struct VM_particle**));
	for (k=0; k<simulation->parameters->box_number_x*simulation->parameters->box_number_y; k++)
	{
		*(simulation->parameters->box+k)=NULL;
	}
	//create box data structure for metric free models
	simulation->parameters->mf_box=malloc(simulation->parameters->mf_box_number_x*simulation->parameters->mf_box_number_y*sizeof(struct VM_particle**));
	for (k=0; k<simulation->parameters->mf_box_number_x*simulation->parameters->mf_box_number_y; k++)
	{
		*(simulation->parameters->mf_box+k)=NULL;
	}
	//create neighborhood data structures
	simulation->state->nh.neighbor_indexes=malloc(simulation->parameters->mf_kn*sizeof(guint64));
	simulation->state->nh.distances=malloc(simulation->parameters->mf_kn*sizeof(gdouble));
	simulation->state->nh.mirror_flag=malloc(simulation->parameters->mf_kn*sizeof(gint32));
	mfVM_clean_neighborhood(simulation);
	//new in version 1.1
	//if data file from version 1.0 was loaded, the number of neighbor moments have to be renormalized
	if (v==0)
	{
		//normalize ensemble averaged neighbor number
		//new in version 1.1
		if (simulation->parameters->particle_species==1)
		{
			*simulation->observables->neighbors_moment1=*simulation->observables->neighbors_moment1/simulation->state->particle_number;
			*simulation->observables->neighbors_moment2=*simulation->observables->neighbors_moment2/simulation->state->particle_number;
			*simulation->observables->neighbors_moment3=*simulation->observables->neighbors_moment3/simulation->state->particle_number;
			*simulation->observables->neighbors_moment4=*simulation->observables->neighbors_moment4/simulation->state->particle_number;
		}
		else
			for (k=0; k<simulation->parameters->particle_species; k++)
			{
				*(simulation->observables->neighbors_moment1+k)=*(simulation->observables->neighbors_moment1+k)/ *(simulation->parameters->species_particle_number+k);
				*(simulation->observables->neighbors_moment2+k)=*(simulation->observables->neighbors_moment2+k)/ *(simulation->parameters->species_particle_number+k);
				*(simulation->observables->neighbors_moment3+k)=*(simulation->observables->neighbors_moment3+k)/ *(simulation->parameters->species_particle_number+k);
				*(simulation->observables->neighbors_moment4+k)=*(simulation->observables->neighbors_moment4+k)/ *(simulation->parameters->species_particle_number+k);
			}
	}
	return PyCapsule_New((void *) simulation, "", NULL);
}

//function called from python as aappp_save
static PyObject* VM_save(PyObject* self, PyObject *args)
{
	PyObject * py_obj;
	char* filename;
	if (!PyArg_ParseTuple(args, "Os", &py_obj, &filename))
	{
		PyErr_SetString(PyExc_TypeError, "Error when saving simulation state to file, see >>>help(aappp_save)<<< for help.\n");
		return NULL;
	}
	struct VM_simulation * simulation= (struct VM_simulation *) PyCapsule_GetPointer(py_obj, "");
	if (VM_save_state(filename, simulation)==-1)
	{
		PyErr_SetString(PyExc_TypeError, "Error when saving simulation state to file, see >>>help(aappp_save)<<< for help.\n");
		return NULL;
	}
	return Py_BuildValue("i", 1);
}

//function called from python as aappp_free
static PyObject* VM_free(PyObject* self, PyObject *args)
{
	gint32 i;
	PyObject * py_obj;
	if (!PyArg_ParseTuple(args, "O", &py_obj))
	{
		PyErr_SetString(PyExc_TypeError, "Error when freeing memory from simulation data structure.\n");
		return NULL;
	}
	struct VM_simulation * simulation= (struct VM_simulation *) PyCapsule_GetPointer(py_obj, "");
	free(simulation->state->particles);
	simulation->state->particles=NULL;
	free(simulation->state->weights);
	simulation->state->weights=NULL;
	free(simulation->state->nh.neighbor_indexes);
	simulation->state->nh.neighbor_indexes=NULL;
	free(simulation->state->nh.distances);
	simulation->state->nh.distances=NULL;
	free(simulation->state->nh.mirror_flag);
	simulation->state->nh.mirror_flag=NULL;
	free(simulation->parameters->box);
	simulation->parameters->box=NULL;
	free(simulation->parameters->mf_box);
	simulation->parameters->mf_box=NULL;
	for (i=0; i<simulation->parameters->particle_species; i++)
	{
		free(*(simulation->parameters->coupling+i));
		*(simulation->parameters->coupling+i)=NULL;
		free(*(simulation->observables->hist_theta+i));
		*(simulation->observables->hist_theta+i)=NULL;
		free(*(simulation->observables->hist_neighbors+i));
		*(simulation->observables->hist_neighbors+i)=NULL;
	}
	free(simulation->parameters->speed);
	simulation->parameters->speed=NULL;
	free(simulation->parameters->noise_strength);
	simulation->parameters->noise_strength=NULL;
	free(simulation->parameters->omega);
	simulation->parameters->omega=NULL;
	free(simulation->parameters->coupling);
	simulation->parameters->coupling=NULL;
	free(simulation->parameters->species_particle_number);
	simulation->parameters->species_particle_number=NULL;
	free(simulation->observables->hist_theta);
	simulation->observables->hist_theta=NULL;
	free(simulation->observables->hist_neighbors);
	simulation->observables->hist_neighbors=NULL;
	free(simulation->observables->polar_order);
	simulation->observables->polar_order=NULL;
	free(simulation->observables->polar_order_moment1);
	simulation->observables->polar_order_moment1=NULL;
	free(simulation->observables->polar_order_moment2);
	simulation->observables->polar_order_moment2=NULL;
	free(simulation->observables->polar_order_moment3);
	simulation->observables->polar_order_moment3=NULL;
	free(simulation->observables->polar_order_moment4);
	simulation->observables->polar_order_moment4=NULL;
	free(simulation->observables->nematic_order);
	simulation->observables->nematic_order=NULL;
	free(simulation->observables->nematic_order_moment1);
	simulation->observables->nematic_order_moment1=NULL;
	free(simulation->observables->nematic_order_moment2);
	simulation->observables->nematic_order_moment2=NULL;
	free(simulation->observables->nematic_order_moment3);
	simulation->observables->nematic_order_moment3=NULL;
	free(simulation->observables->nematic_order_moment4);
	simulation->observables->nematic_order_moment4=NULL;
	free(simulation->observables->neighbors_moment1);
	simulation->observables->neighbors_moment1=NULL;
	free(simulation->observables->neighbors_moment2);
	simulation->observables->neighbors_moment2=NULL;
	free(simulation->observables->neighbors_moment3);
	simulation->observables->neighbors_moment3=NULL;
	free(simulation->observables->neighbors_moment4);
	simulation->observables->neighbors_moment4=NULL;
	free(simulation->state);
	simulation->state=NULL;
	free(simulation->parameters);
	simulation->parameters=NULL;
	free(simulation->observables);
	simulation->observables=NULL;
	free(simulation);
	simulation=NULL;
	return Py_BuildValue("i", 1);
}

//function called from python as VM_update_timesteps
static PyObject* VM_update_timesteps(PyObject* self, PyObject *args)
{
	PyObject * py_obj;
	guint64 timesteps;
	guint64 i;
	if (!PyArg_ParseTuple(args, "OL", &py_obj, &timesteps))
	{
		PyErr_SetString(PyExc_TypeError, "Error when time evolving Vicsek model, see >>>help(VM_update_timesteps)<<< for help.\n");
		return NULL;
	}
	struct VM_simulation * simulation= (struct VM_simulation *) PyCapsule_GetPointer(py_obj, "");
	for (i=0; i<timesteps; i++)
	{
		VM_one_step(simulation);
	}
	return Py_BuildValue("d", *simulation->observables->polar_order);
}

//function called from python as NVM_update_timesteps
static PyObject* NVM_update_timesteps(PyObject* self, PyObject *args)
{
	PyObject * py_obj;
	guint64 timesteps;
	guint64 i;
	if (!PyArg_ParseTuple(args, "OL", &py_obj, &timesteps))
	{
		PyErr_SetString(PyExc_TypeError, "Error when time evolving nematic Vicsek model, see >>>help(NVM_update_timesteps)<<< for help.\n");
		return NULL;
	}
	struct VM_simulation * simulation= (struct VM_simulation *) PyCapsule_GetPointer(py_obj, "");
	for (i=0; i<timesteps; i++)
	{
		NVM_one_step(simulation);
	}
	return Py_BuildValue("d", *simulation->observables->polar_order);
}

//function called from python as mfVM_update_timesteps
static PyObject* mfVM_update_timesteps(PyObject* self, PyObject *args)
{
	PyObject * py_obj;
	guint64 timesteps;
	guint64 i;
	if (!PyArg_ParseTuple(args, "OL", &py_obj, &timesteps))
	{
		PyErr_SetString(PyExc_TypeError, "Error when time evolving metric free Vicsek model, see >>>help(mfVM_update_timesteps)<<< for help.\n");
		return NULL;
	}
	struct VM_simulation * simulation= (struct VM_simulation *) PyCapsule_GetPointer(py_obj, "");
	for (i=0; i<timesteps; i++)
	{
		mfVM_one_step(simulation);
	}
	return Py_BuildValue("d", *simulation->observables->polar_order);
}

//function called from python as mfNVM_update_timesteps
static PyObject* mfNVM_update_timesteps(PyObject* self, PyObject *args)
{
	PyObject * py_obj;
	guint64 timesteps;
	guint64 i;
	if (!PyArg_ParseTuple(args, "OL", &py_obj, &timesteps))
	{
		PyErr_SetString(PyExc_TypeError, "Error when time evolving metric free nematic Vicsek model, see >>>help(mfNVM_update_timesteps)<<< for help.\n");
		return NULL;
	}
	struct VM_simulation * simulation= (struct VM_simulation *) PyCapsule_GetPointer(py_obj, "");
	for (i=0; i<timesteps; i++)
	{
		mfNVM_one_step(simulation);
	}
	return Py_BuildValue("d", *simulation->observables->polar_order);
}

//function called from python as additiveL_update_timesteps
static PyObject* additiveL_update_timesteps(PyObject* self, PyObject *args)
{
	PyObject * py_obj;
	guint64 timesteps;
	guint64 i;
	if (!PyArg_ParseTuple(args, "OL", &py_obj, &timesteps))
	{
		PyErr_SetString(PyExc_TypeError, "Error when time evolving additive overdamped Langevin model, see >>>help(additiveL_update_timesteps)<<< for help.\n");
		return NULL;
	}
	struct VM_simulation * simulation= (struct VM_simulation *) PyCapsule_GetPointer(py_obj, "");
	for (i=0; i<timesteps; i++)
	{
		additiveL_one_step(simulation);
	}
	return Py_BuildValue("d", *simulation->observables->polar_order);
}

//function called from python as nonadditiveL_update_timesteps
static PyObject* nonadditiveL_update_timesteps(PyObject* self, PyObject *args)
{
	PyObject * py_obj;
	guint64 timesteps;
	guint64 i;
	if (!PyArg_ParseTuple(args, "OL", &py_obj, &timesteps))
	{
		PyErr_SetString(PyExc_TypeError, "Error when time evolving nonadditive overdamped Langevin model, see >>>help(nonadditiveL_update_timesteps)<<< for help.\n");
		return NULL;
	}
	struct VM_simulation * simulation= (struct VM_simulation *) PyCapsule_GetPointer(py_obj, "");
	for (i=0; i<timesteps; i++)
	{
		nonadditiveL_one_step(simulation);
	}
	return Py_BuildValue("d", *simulation->observables->polar_order);
}

//function called from python as mfL_update_timesteps
static PyObject* mfL_update_timesteps(PyObject* self, PyObject *args)
{
	PyObject * py_obj;
	guint64 timesteps;
	guint64 i;
	if (!PyArg_ParseTuple(args, "OL", &py_obj, &timesteps))
	{
		PyErr_SetString(PyExc_TypeError, "Error when time evolving metric free overdamped Langevin model, see >>>help(mfL_update_timesteps)<<< for help.\n");
		return NULL;
	}
	struct VM_simulation * simulation= (struct VM_simulation *) PyCapsule_GetPointer(py_obj, "");
	for (i=0; i<timesteps; i++)
	{
		mfL_one_step(simulation);
	}
	return Py_BuildValue("d", *simulation->observables->polar_order);
}

//function called from python as VM_measurement_timesteps
static PyObject* VM_measurement_timesteps(PyObject* self, PyObject *args)
{
	PyObject * py_obj;
	guint64 timesteps;
	guint64 i;
	if (!PyArg_ParseTuple(args, "OL", &py_obj, &timesteps))
	{
		PyErr_SetString(PyExc_TypeError, "Error when time evolving and measuring Vicsek model, see >>>help(VM_measurement_timesteps)<<< for help.\n");
		return NULL;
	}
	struct VM_simulation * simulation= (struct VM_simulation *) PyCapsule_GetPointer(py_obj, "");
	for (i=0; i<timesteps; i++)
	{
		VM_one_step_measurement(simulation);
	}
	return Py_BuildValue("d", *simulation->observables->polar_order);
}

//function called from python as NVM_measurement_timesteps
static PyObject* NVM_measurement_timesteps(PyObject* self, PyObject *args)
{
	PyObject * py_obj;
	guint64 timesteps;
	guint64 i;
	if (!PyArg_ParseTuple(args, "OL", &py_obj, &timesteps))
	{
		PyErr_SetString(PyExc_TypeError, "Error when time evolving and measuring nematic Vicsek model, see >>>help(NVM_measurement_timesteps)<<< for help.\n");
		return NULL;
	}
	struct VM_simulation * simulation= (struct VM_simulation *) PyCapsule_GetPointer(py_obj, "");
	for (i=0; i<timesteps; i++)
	{
		NVM_one_step_measurement(simulation);
	}
	return Py_BuildValue("d", *simulation->observables->polar_order);
}

//function called from python as mfVM_measurement_timesteps
static PyObject* mfVM_measurement_timesteps(PyObject* self, PyObject *args)
{
	PyObject * py_obj;
	guint64 timesteps;
	guint64 i;
	if (!PyArg_ParseTuple(args, "OL", &py_obj, &timesteps))
	{
		PyErr_SetString(PyExc_TypeError, "Error when time evolving and measuring metric free Vicsek model, see >>>help(mfVM_measurement_timesteps)<<< for help.\n");
		return NULL;
	}
	struct VM_simulation * simulation= (struct VM_simulation *) PyCapsule_GetPointer(py_obj, "");
	for (i=0; i<timesteps; i++)
	{
		mfVM_one_step_measurement(simulation);
	}
	return Py_BuildValue("d", *simulation->observables->polar_order);
}

//function called from python as mfNVM_measurement_timesteps
static PyObject* mfNVM_measurement_timesteps(PyObject* self, PyObject *args)
{
	PyObject * py_obj;
	guint64 timesteps;
	guint64 i;
	if (!PyArg_ParseTuple(args, "OL", &py_obj, &timesteps))
	{
		PyErr_SetString(PyExc_TypeError, "Error when time evolving and measuring metric free nematic Vicsek model, see >>>help(mfNVM_measurement_timesteps)<<< for help.\n");
		return NULL;
	}
	struct VM_simulation * simulation= (struct VM_simulation *) PyCapsule_GetPointer(py_obj, "");
	for (i=0; i<timesteps; i++)
	{
		mfNVM_one_step_measurement(simulation);
	}
	return Py_BuildValue("d", *simulation->observables->polar_order);
}

//function called from python as additiveL_measurement_timesteps
static PyObject* additiveL_measurement_timesteps(PyObject* self, PyObject *args)
{
	PyObject * py_obj;
	guint64 timesteps;
	guint64 i;
	if (!PyArg_ParseTuple(args, "OL", &py_obj, &timesteps))
	{
		PyErr_SetString(PyExc_TypeError, "Error when time evolving and measuring additive overdamped Langevin model, see >>>help(additiveL_measurement_timesteps)<<< for help.\n");
		return NULL;
	}
	struct VM_simulation * simulation= (struct VM_simulation *) PyCapsule_GetPointer(py_obj, "");
	for (i=0; i<timesteps; i++)
	{
		additiveL_one_step_measurement(simulation);
	}
	return Py_BuildValue("d", *simulation->observables->polar_order);
}

//function called from python as nonadditiveL_measurement_timesteps
static PyObject* nonadditiveL_measurement_timesteps(PyObject* self, PyObject *args)
{
	PyObject * py_obj;
	guint64 timesteps;
	guint64 i;
	if (!PyArg_ParseTuple(args, "OL", &py_obj, &timesteps))
	{
		PyErr_SetString(PyExc_TypeError, "Error when time evolving and measuring nonadditive overdamped Langevin model, see >>>help(nonadditiveL_measurement_timesteps)<<< for help.\n");
		return NULL;
	}
	struct VM_simulation * simulation= (struct VM_simulation *) PyCapsule_GetPointer(py_obj, "");
	for (i=0; i<timesteps; i++)
	{
		nonadditiveL_one_step_measurement(simulation);
	}
	return Py_BuildValue("d", *simulation->observables->polar_order);
}

//function called from python as mfL_measurement_timesteps
static PyObject* mfL_measurement_timesteps(PyObject* self, PyObject *args)
{
	PyObject * py_obj;
	guint64 timesteps;
	guint64 i;
	if (!PyArg_ParseTuple(args, "OL", &py_obj, &timesteps))
	{
		PyErr_SetString(PyExc_TypeError, "Error when time evolving and measuring metric free overdamped Langevin model, see >>>help(mfL_measurement_timesteps)<<< for help.\n");
		return NULL;
	}
	struct VM_simulation * simulation= (struct VM_simulation *) PyCapsule_GetPointer(py_obj, "");
	for (i=0; i<timesteps; i++)
	{
		mfL_one_step_measurement(simulation);
	}
	return Py_BuildValue("d", *simulation->observables->polar_order);
}

//function called from python as aappp_get_parameters
static PyObject* VM_print_parameters(PyObject* self, PyObject *args)
{
	gint32 k, j;
	PyObject * py_obj;
	if (!PyArg_ParseTuple(args, "O", &py_obj))
	{
		PyErr_SetString(PyExc_TypeError, "Error when printing parameters\n");
		return NULL;
	}
	struct VM_simulation * simulation= (struct VM_simulation *) PyCapsule_GetPointer(py_obj, "");
	PyObject * return_value=PyList_New(15); //we return 15 numbers/objects
	PyObject * py_v;
	PyObject * py_eta;
	PyObject * py_omega;
	PyObject * py_gamma;
	PyObject * py_N;
	if (simulation->parameters->particle_species==1)
	{
		py_v=PyFloat_FromDouble(*simulation->parameters->speed);
		py_eta=PyFloat_FromDouble(*simulation->parameters->noise_strength);
		py_omega=PyFloat_FromDouble(*simulation->parameters->omega);
		py_gamma=PyFloat_FromDouble(**simulation->parameters->coupling);
		py_N=PyLong_FromUnsignedLongLong(simulation->state->particle_number);
	}
	else
	{
		py_v=PyList_New(simulation->parameters->particle_species);
		py_eta=PyList_New(simulation->parameters->particle_species);
		py_omega=PyList_New(simulation->parameters->particle_species);
		py_gamma=PyList_New(simulation->parameters->particle_species);
		py_N=PyList_New(simulation->parameters->particle_species);
		for (k=0; k<simulation->parameters->particle_species; k++)
		{
			PyList_SET_ITEM(py_v, k, PyFloat_FromDouble(*(simulation->parameters->speed+k)));
			PyList_SET_ITEM(py_eta, k, PyFloat_FromDouble(*(simulation->parameters->noise_strength+k)));
			PyList_SET_ITEM(py_omega, k, PyFloat_FromDouble(*(simulation->parameters->omega+k)));
			PyList_SET_ITEM(py_N, k, PyLong_FromUnsignedLongLong(*(simulation->parameters->species_particle_number+k)));
			PyObject * py_c=PyList_New(simulation->parameters->particle_species);
			for (j=0; j<simulation->parameters->particle_species; j++)
				PyList_SET_ITEM(py_c, j, PyFloat_FromDouble(*(*(simulation->parameters->coupling+j)+k)));
			PyList_SET_ITEM(py_gamma, k, py_c);
		}
	}
	PyList_SET_ITEM(return_value, 0, py_v);
	PyList_SET_ITEM(return_value, 1, py_eta);
	PyList_SET_ITEM(return_value, 2, PyFloat_FromDouble(simulation->parameters->interaction_radius));
	PyList_SET_ITEM(return_value, 3, py_omega);
	PyList_SET_ITEM(return_value, 4, PyFloat_FromDouble(simulation->state->length_x));
	PyList_SET_ITEM(return_value, 5, PyFloat_FromDouble(simulation->state->length_y));
	PyList_SET_ITEM(return_value, 6, py_gamma);
	PyList_SET_ITEM(return_value, 7, PyFloat_FromDouble(simulation->parameters->delta_t));
	PyList_SET_ITEM(return_value, 8, py_N);
	PyList_SET_ITEM(return_value, 9, PyLong_FromLong((long) simulation->parameters->mf_kn));
	PyList_SET_ITEM(return_value, 10, PyLong_FromLong((long) simulation->parameters->interaction_order));
	PyList_SET_ITEM(return_value, 11, PyLong_FromLong((long) simulation->observables->binnum_theta));
	PyList_SET_ITEM(return_value, 12, PyLong_FromLong((long) simulation->observables->binnum_neighbors));
	PyList_SET_ITEM(return_value, 13, PyLong_FromLong((long) simulation->parameters->boundary_x));
	PyList_SET_ITEM(return_value, 14, PyLong_FromLong((long) simulation->parameters->boundary_y));
	return return_value;
}

//function called from python as aappp_get_xythetalxly
static PyObject* VM_xythetalxly(PyObject* self, PyObject *args)
{
	guint64 i;
	PyObject * py_obj;
	if (!PyArg_ParseTuple(args, "O", &py_obj))
	{
		PyErr_SetString(PyExc_TypeError, "Error when printing positions, orientations, simulation box dimensions.\n");
		return NULL;
	}
	struct VM_simulation * simulation= (struct VM_simulation *) PyCapsule_GetPointer(py_obj, "");
	Py_ssize_t len = simulation->state->particle_number;
	PyObject * x= PyList_New(len);
	PyObject * y= PyList_New(len);
	PyObject * theta= PyList_New(len);
	for (i=0; i<simulation->state->particle_number; i++)
	{
		PyList_SET_ITEM(x, i, PyFloat_FromDouble((simulation->state->particles+i)->x));
		PyList_SET_ITEM(y, i, PyFloat_FromDouble((simulation->state->particles+i)->y));
		PyList_SET_ITEM(theta, i, PyFloat_FromDouble((simulation->state->particles+i)->theta));
	}
	PyObject *MyResult = Py_BuildValue("OOOdd", x, y, theta, simulation->state->length_x, simulation->state->length_y);
	//this is to prevent memory leaks
	Py_DECREF(x);  
	Py_DECREF(y);  
	Py_DECREF(theta);  
	return MyResult;
}

//function called from python as aappp_get_weights
static PyObject* aappp_weights(PyObject* self, PyObject *args)
{
	gint32 i;
	PyObject * py_obj;
	if (!PyArg_ParseTuple(args, "O", &py_obj))
	{
		PyErr_SetString(PyExc_TypeError, "Error when printing weights for nonadditive overdamped Langevin dynamics.\n");
		return NULL;
	}
	struct VM_simulation * simulation= (struct VM_simulation *) PyCapsule_GetPointer(py_obj, "");
	Py_ssize_t len = simulation->state->weight_vector_length;
	PyObject * weights= PyList_New(len);
	for (i=0; i<len; i++)
	{
		PyList_SET_ITEM(weights, i, PyFloat_FromDouble(*(simulation->state->weights+i) ));
	}
	return weights;
}

//function called from python as aappp_get_results
static PyObject* VM_measurement_results(PyObject* self, PyObject *args)
{
	gint32 i;
	gint32 k;
	PyObject * py_obj;
	if (!PyArg_ParseTuple(args, "O", &py_obj))
	{
		PyErr_SetString(PyExc_TypeError, "Error when printing measurement results, see >>>help(aappp_get_results)<<< for help.\n");
		return NULL;
	}
	struct VM_simulation * simulation= (struct VM_simulation *) PyCapsule_GetPointer(py_obj, "");
	guint64 measurement_steps= simulation->observables->measurement_steps;
	if (measurement_steps==0)
		measurement_steps=1;
	gdouble binsize_theta=2.0*M_PI/simulation->observables->binnum_theta;
	Py_ssize_t len1 = simulation->observables->binnum_theta;
	Py_ssize_t len2 = simulation->observables->binnum_neighbors;
	PyObject * distribution_theta=NULL;
	PyObject * distribution_neighbors=NULL;
	PyObject * polar=NULL;
	PyObject * nematic=NULL;
	PyObject * neighbors=NULL;
	//data to hold return structure factor
	PyObject * structure_factor=NULL;
	if (simulation->observables->sf_mode_number==0)
	{
		structure_factor=PyList_New(0);
	}
	else
	{
		if (simulation->parameters->particle_species==1)
		{
			structure_factor=PyList_New(simulation->observables->sf_mode_number);
			for (k=0; k<simulation->observables->sf_mode_number; k++)
			{
				PyList_SET_ITEM(structure_factor, k, PyFloat_FromDouble(*(*(simulation->observables->structure_factor)+k)/simulation->observables->measurement_steps));
			}
		}
		else
		{
			structure_factor=PyList_New(simulation->parameters->particle_species);
			for (i=0; i<simulation->parameters->particle_species; i++)
			{
				PyObject * py_sf=PyList_New(simulation->observables->sf_mode_number);
				for (k=0; k<simulation->observables->sf_mode_number; k++)
				{
					PyList_SET_ITEM(py_sf, k, PyFloat_FromDouble(*(*(simulation->observables->structure_factor+i)+k)/simulation->observables->measurement_steps));
				}
				PyList_SET_ITEM(structure_factor, i, py_sf);
			}
		}
	}
	if (simulation->parameters->particle_species==1)
	{
		PyObject * theta= PyList_New(len1);
		PyObject * ptheta= PyList_New(len1);
		PyObject * n= PyList_New(len2);
		PyObject * pn= PyList_New(len2);
		polar= PyList_New(5);
		nematic= PyList_New(5);
		neighbors= PyList_New(4);
		distribution_theta= PyList_New(2);
		distribution_neighbors= PyList_New(2);
		for (i=0; i<simulation->observables->binnum_theta; i++)
		{
			PyList_SET_ITEM(theta, i, PyFloat_FromDouble( (i+0.5)*binsize_theta -M_PI ));
			PyList_SET_ITEM(ptheta, i, PyFloat_FromDouble((1.0*(*(*simulation->observables->hist_theta+i)))/(simulation->state->particle_number*measurement_steps)/binsize_theta));
		}
		PyList_SET_ITEM(distribution_theta, 0, theta);
		PyList_SET_ITEM(distribution_theta, 1, ptheta);
		for (i=0; i<simulation->observables->binnum_neighbors; i++)
		{
			PyList_SET_ITEM(n, i, PyLong_FromUnsignedLongLong( i ));
			PyList_SET_ITEM(pn, i, PyFloat_FromDouble((1.0*(*(*simulation->observables->hist_neighbors+i)))/(simulation->state->particle_number*measurement_steps)));
		}
		PyList_SET_ITEM(distribution_neighbors, 0, n);
		PyList_SET_ITEM(distribution_neighbors, 1, pn);
		PyList_SET_ITEM(polar, 0, PyFloat_FromDouble(*simulation->observables->polar_order));
		PyList_SET_ITEM(polar, 1, PyFloat_FromDouble(*simulation->observables->polar_order_moment1/measurement_steps));
		PyList_SET_ITEM(polar, 2, PyFloat_FromDouble(*simulation->observables->polar_order_moment2/measurement_steps));
		PyList_SET_ITEM(polar, 3, PyFloat_FromDouble(*simulation->observables->polar_order_moment3/measurement_steps));
		PyList_SET_ITEM(polar, 4, PyFloat_FromDouble(*simulation->observables->polar_order_moment4/measurement_steps));
		PyList_SET_ITEM(nematic, 0, PyFloat_FromDouble(*simulation->observables->nematic_order));
		PyList_SET_ITEM(nematic, 1, PyFloat_FromDouble(*simulation->observables->nematic_order_moment1/measurement_steps));
		PyList_SET_ITEM(nematic, 2, PyFloat_FromDouble(*simulation->observables->nematic_order_moment2/measurement_steps));
		PyList_SET_ITEM(nematic, 3, PyFloat_FromDouble(*simulation->observables->nematic_order_moment3/measurement_steps));
		PyList_SET_ITEM(nematic, 4, PyFloat_FromDouble(*simulation->observables->nematic_order_moment4/measurement_steps));
		//change in version 1.1, neighbor moments are already normalized with respect to particle number
		PyList_SET_ITEM(neighbors, 0, PyFloat_FromDouble(*simulation->observables->neighbors_moment1/(measurement_steps)));
		PyList_SET_ITEM(neighbors, 1, PyFloat_FromDouble(*simulation->observables->neighbors_moment2/(measurement_steps)));
		PyList_SET_ITEM(neighbors, 2, PyFloat_FromDouble(*simulation->observables->neighbors_moment3/(measurement_steps)));
		PyList_SET_ITEM(neighbors, 3, PyFloat_FromDouble(*simulation->observables->neighbors_moment4/(measurement_steps)));
	}
	else
	{
		polar= PyList_New(simulation->parameters->particle_species);
		nematic= PyList_New(simulation->parameters->particle_species);
		neighbors= PyList_New(simulation->parameters->particle_species);
		distribution_theta= PyList_New(simulation->parameters->particle_species);
		distribution_neighbors= PyList_New(simulation->parameters->particle_species);
		for (k=0; k<simulation->parameters->particle_species; k++)
		{
			PyObject * theta= PyList_New(len1);
			PyObject * ptheta= PyList_New(len1);
			PyObject * n= PyList_New(len2);
			PyObject * pn= PyList_New(len2);
			PyObject * kpolar= PyList_New(5);
			PyObject * knematic= PyList_New(5);
			PyObject * kneighbors= PyList_New(4);
			PyObject * kdistribution_theta= PyList_New(2);
			PyObject * kdistribution_neighbors= PyList_New(2);
			for (i=0; i<simulation->observables->binnum_theta; i++)
			{
				PyList_SET_ITEM(theta, i, PyFloat_FromDouble( (i+0.5)*binsize_theta -M_PI ));
				PyList_SET_ITEM(ptheta, i, PyFloat_FromDouble((1.0*(*(*(simulation->observables->hist_theta+k)+i)))/(*(simulation->parameters->species_particle_number+k)*measurement_steps)/binsize_theta));
			}
			PyList_SET_ITEM(kdistribution_theta, 0, theta);
			PyList_SET_ITEM(kdistribution_theta, 1, ptheta);
			for (i=0; i<simulation->observables->binnum_neighbors; i++)
			{
				PyList_SET_ITEM(n, i, PyLong_FromUnsignedLongLong( i ));
				PyList_SET_ITEM(pn, i, PyFloat_FromDouble((1.0*(*(*(simulation->observables->hist_neighbors+k)+i)))/(*(simulation->parameters->species_particle_number+k)*measurement_steps)));
			}
			PyList_SET_ITEM(kdistribution_neighbors, 0, n);
			PyList_SET_ITEM(kdistribution_neighbors, 1, pn);
			PyList_SET_ITEM(kpolar, 0, PyFloat_FromDouble(*(simulation->observables->polar_order+k)));
			PyList_SET_ITEM(kpolar, 1, PyFloat_FromDouble(*(simulation->observables->polar_order_moment1+k)/measurement_steps));
			PyList_SET_ITEM(kpolar, 2, PyFloat_FromDouble(*(simulation->observables->polar_order_moment2+k)/measurement_steps));
			PyList_SET_ITEM(kpolar, 3, PyFloat_FromDouble(*(simulation->observables->polar_order_moment3+k)/measurement_steps));
			PyList_SET_ITEM(kpolar, 4, PyFloat_FromDouble(*(simulation->observables->polar_order_moment4+k)/measurement_steps));
			PyList_SET_ITEM(knematic, 0, PyFloat_FromDouble(*(simulation->observables->nematic_order+k)));
			PyList_SET_ITEM(knematic, 1, PyFloat_FromDouble(*(simulation->observables->nematic_order_moment1+k)/measurement_steps));
			PyList_SET_ITEM(knematic, 2, PyFloat_FromDouble(*(simulation->observables->nematic_order_moment2+k)/measurement_steps));
			PyList_SET_ITEM(knematic, 3, PyFloat_FromDouble(*(simulation->observables->nematic_order_moment3+k)/measurement_steps));
			PyList_SET_ITEM(knematic, 4, PyFloat_FromDouble(*(simulation->observables->nematic_order_moment4+k)/measurement_steps));
			//change in version 1.1, neighbor moments are already normalized with respect to particle number
			PyList_SET_ITEM(kneighbors, 0, PyFloat_FromDouble(*(simulation->observables->neighbors_moment1+k)/(measurement_steps)));
			PyList_SET_ITEM(kneighbors, 1, PyFloat_FromDouble(*(simulation->observables->neighbors_moment2+k)/(measurement_steps)));
			PyList_SET_ITEM(kneighbors, 2, PyFloat_FromDouble(*(simulation->observables->neighbors_moment3+k)/(measurement_steps)));
			PyList_SET_ITEM(kneighbors, 3, PyFloat_FromDouble(*(simulation->observables->neighbors_moment4+k)/(measurement_steps)));
			PyList_SET_ITEM(distribution_theta, k, kdistribution_theta);
			PyList_SET_ITEM(distribution_neighbors, k, kdistribution_neighbors);
			PyList_SET_ITEM(polar, k, kpolar);
			PyList_SET_ITEM(nematic, k, knematic);
			PyList_SET_ITEM(neighbors, k, kneighbors);
		}
	}
	PyObject * MyResult = Py_BuildValue("OOOOOO", distribution_theta, distribution_neighbors, polar, nematic, neighbors, structure_factor);
	Py_DECREF(distribution_theta);
	Py_DECREF(distribution_neighbors);
	Py_DECREF(polar);
	Py_DECREF(nematic);
	Py_DECREF(neighbors);
	Py_DECREF(structure_factor);
	return MyResult;
}

//function called from python as aappp_set_state
static PyObject* VM_set_state(PyObject* self, PyObject *args)
{
	PyObject * py_obj;
	PyObject * x;
	PyObject * y;
	PyObject *theta;
	gdouble lx, ly;
	guint64 i;
	gint32 k;
	if (!PyArg_ParseTuple(args, "OOOOdd", &py_obj, &x, &y, &theta, &lx, &ly))
		return Py_BuildValue("i", -1);
	struct VM_simulation * simulation= (struct VM_simulation *) PyCapsule_GetPointer(py_obj, "");
	if (PyList_Size(x)!=(Py_ssize_t) simulation->state->particle_number)
		return Py_BuildValue("i", -1);
	if (PyList_Size(y)!=(Py_ssize_t) simulation->state->particle_number)
		return Py_BuildValue("i", -1);
	if (PyList_Size(theta)!=(Py_ssize_t) simulation->state->particle_number)
		return Py_BuildValue("i", -1);
	simulation->state->length_x=lx;
	simulation->state->length_y=ly;
	//create new box lists
	simulation->parameters->mf_interaction_radius=sqrt(simulation->state->length_x*simulation->state->length_y*simulation->parameters->mf_kn/M_PI/simulation->state->particle_number);
	k=1;
	while (simulation->state->length_x/(k+1)>simulation->parameters->interaction_radius)
		k++;
	simulation->parameters->box_number_x=k;
	simulation->parameters->box_size_x=simulation->state->length_x/k;
	k=1;
	while (simulation->state->length_y/(k+1)>simulation->parameters->interaction_radius)
		k++;
	simulation->parameters->box_number_y=k;
	simulation->parameters->box_size_y=simulation->state->length_y/k;
	free(simulation->parameters->box);
	simulation->parameters->box=malloc(simulation->parameters->box_number_x*simulation->parameters->box_number_y*sizeof(struct VM_particle**));
	for (k=0; k<simulation->parameters->box_number_x*simulation->parameters->box_number_y; k++)
	{
		*(simulation->parameters->box+k)=NULL;
	}
	k=1;
	while (simulation->state->length_x/(k+1)>simulation->parameters->mf_interaction_radius)
		k++;
	simulation->parameters->mf_box_number_x=k;
	simulation->parameters->mf_box_size_x=simulation->state->length_x/k;
	k=1;
	while (simulation->state->length_y/(k+1)>simulation->parameters->mf_interaction_radius)
		k++;
	simulation->parameters->mf_box_number_y=k;
	simulation->parameters->mf_box_size_y=simulation->state->length_y/k;
	free(simulation->parameters->mf_box);
	simulation->parameters->mf_box=malloc(simulation->parameters->mf_box_number_x*simulation->parameters->mf_box_number_y*sizeof(struct VM_particle**));
	for (k=0; k<simulation->parameters->mf_box_number_x*simulation->parameters->mf_box_number_y; k++)
	{
		*(simulation->parameters->mf_box+k)=NULL;
	}
	//read in particles positions
	for (i=0; i<simulation->state->particle_number; i++)
	{
		(simulation->state->particles+i)->x=PyFloat_AsDouble(PyList_GetItem(x, i));
		while ((simulation->state->particles+i)->x<0.0)
			(simulation->state->particles+i)->x+=simulation->state->length_x;
		while ((simulation->state->particles+i)->x>simulation->state->length_x)
			(simulation->state->particles+i)->x-=simulation->state->length_x;
		(simulation->state->particles+i)->y=PyFloat_AsDouble(PyList_GetItem(y, i));
		while ((simulation->state->particles+i)->y<0.0)
			(simulation->state->particles+i)->y+=simulation->state->length_y;
		while ((simulation->state->particles+i)->y>simulation->state->length_y)
			(simulation->state->particles+i)->y-=simulation->state->length_y;
		(simulation->state->particles+i)->theta=PyFloat_AsDouble(PyList_GetItem(theta, i));
		while ((simulation->state->particles+i)->theta<-M_PI)
			(simulation->state->particles+i)->theta+=2.0*M_PI;
		while ((simulation->state->particles+i)->theta>M_PI)
			(simulation->state->particles+i)->theta-=2.0*M_PI;
	}
	return Py_BuildValue("i", 1);
}

static PyMethodDef aappp_methods[] = {
	{"aappp_init", (PyCFunction)VM_init, METH_VARARGS | METH_KEYWORDS, "aappp_init(v, eta, R, omega, Lx, Ly, gamma, dt, N, seed, kn, order, binnum_theta, binnum_neighbors, bx, by, weight_function, weight_vector_length)\n--\n\nReturns a newly created simulation object.\n\nArguments are optional and can be given with keyword (default values in brackets):\nInteger particle Number: N(1000)\nFloat interaction radius: R(1)\nFloat simulation box length in x-direction: Lx(20)\nFloat simulation box length in y-direction: Ly(20)\nFloat step size for Euler integration: dt(0.01)\nInteger boundary flag for x direction (if =0 ->periodic boundary, if =1 ->reflecting boundary): bx(0)\nInteger boundary flag for y direction (if =0 ->periodic boundary, if =1 ->reflecting boundary): by(0)\nInteger order of the interaction term for overdamped Langevin dynamics (the term has the form sin(order*(phi2-phi1))): order(1)\nInteger seed of the pseudo random number generator: seed(1)\nInteger number of interaction neighbors other than focus particle for metric free models: kn(6)\nInteger number of bins for angular histogram: binnum_theta(100)\nInteger number of bins for histogram of the number of neighbors: binnum_neighbors(100)\nFloat speed of the particles: v(1)\nFloat noise strength (which should be a number between 0 and 1 for Vicsek models and an arbitrary number for overdamped Langevin models): eta(0.5)\nFloat natural rotation frequency for each particle: omega(0)\nFloat coupling strength for overdamped Langevin models: gamma(0.1)\nPython function that gives a weight to interactions depending on the number of neighbors: weight_function (None)\nIneger specifying for how many values the weight function is calcuated: weight_vector_length (100)\n\nN, v, eta, omega can be lists where each entry corresponds to one species of particles (different species behave differently).\ngamma can be a list of lists where gamma[A][B] describes, in the equation of motion of an A-particle, the strength of interaction with a B-particle.\nIt is possible that e.g. N and v are lists (of the same dimension) and eta, omega ang gamma are no lists.\nLists of length 1 are not allowed (use just the number in that case).\nIf lists are used, the dimensions for those parameters must be compatible.\nOther parameters than listed before can only be single numbers.\nThe weight function is only relevant for the nonadditive overdamped Langevin dynamics, see >>>help(nonadditiveL_update_timesteps)<<<.\nThe weight function must take an integer as a argument and return a float.\nIf no weight function is given f(n)=1/(n+1) is used.\nThe weight function is only called once for the arguments 0, 1, 2, ..., weight_vector_length-1.\nThe results of those calls are stored in memory and used during simulation.\nIf a particle has more than weight_vector_length-1 neighbors the weight for exactly weight_vector_length-1 neighbors is used.\nparticles are initalized at random with unform distribution over space and with uniform distributed orientations.\n\nEXAMPLE 1 (with only one particle species)\nmysim=aappp.aappp_init(N=100000, Lx=100., Ly=100., v=1.5, R=1.1, eta=0.23, seed=7)\n\nEXAMPLE 2 (with three particle species)\nmysim=aappp.aappp_init(N=[10000, 20000, 15000], gamma=[[0.1, 0.13, 0.15],[0.08, 0.15, 0.05],[0.11, 0.14, 0.11]], Lx=80., Ly=120., dt=0.03, v=1, eta=1.5, seed=1)\nIn this example, if particle i is from species two and particle j is from species three, and i and j are neighbors, there appears the following coupling term (depending on the model, potentially with a weighting prefactor depending on the number of neighbors of particle i):\nd/dt phi_i = 0.05*sin((phi_j-phi_i)*order) + ...\n"},
	{"aappp_free", (PyCFunction)VM_free, METH_VARARGS, "aappp_free(simulation)\n--\n\nfrees all memory that is allocated to the *simulation* data structure.\nAll data get lost.\n"},
	{"aappp_save", (PyCFunction)VM_save, METH_VARARGS, "aappp_save(simulation, filename)\n--\n\nsaves the *simulation* data structure containing all particles positions and orientations, all parameters, the state of the pseudo random number generator and the results of measurements into a file under path *filename*.\nAll data are saved in binary format from C-data types.\nNote that e.g. C doubles are not type safe.\nThus, performing simulation on unusual devices such as mirco processors might result in data files that are unreadable on a different device such as a desktop computer.\nHowever, such behavior is unexpected when using usual desktop computers or high performence computing clusters only.\nIn case that such problems occur it is recommended to analyze data on the same machine that performed the simulation.\nReturns 1 on success.\n"},
	{"aappp_load", (PyCFunction)VM_load, METH_VARARGS, "aappp_load(filename)\n--\n\nloads and returns a simulation data structure from a file under path *filename*, that was previously saved using aappp_save, see >>>help(aappp_save)<<<.\n"},
	{"aappp_set_state", (PyCFunction)VM_set_state, METH_VARARGS, "aappp_set_state(simulation, x, y, theta, lx, ly)\n--\n\nallows to change the state of the simulation.\n*x*, *y*, and *theta* are lists of x-positions, y-positions and orientations of N particles.\nthe length of *x*, *y*, *theta* must be exactly N.\nIt is not possible to modify the particle number of a given simulation.\nIf it is desired to change the particle number one needs to create a new simulation data structure using aappp_init, see >>>help(aappp_init)<<<.\n*lx* and *ly* are the new simulation box sizes in x- and y-direction.\nAll arguments must be given.\nIf particle positions are outside simulation box or orientations are outside [-pi, pi], they are put inside via periodic boundary conditions.\nReturns 1 on success, -1 when failing.\n"},
	{"aappp_reset_parameters", (PyCFunction)VM_reset_parameters, METH_VARARGS | METH_KEYWORDS, "aappp_reset_parameters(simulation, v, R, eta, omega, dt, kn, gamma, order, bx, by, N, weight_function, weight_vector_length)\n--\n\nIs used to change (some) simulation parameters for the following time steps.\nExcept for the simulation data structure (*simulation*) all arguments are otional and can be given by keyword.\nThe total particle number can not be changed.\nIn order to change the total particle number, a new simulation data structure has to be created using aappp_init, see >>>help(aappp_init)<<<.\nHowever, one can change the number of species and the number of particles in each species as long as the total particle number is not changed.\nThis can be done by giving a list for *N*=[N1, N2, ...].\nSimilarly, *v*, *omega* and *eta* can be lists specifying velocity, chirality and noise strength for each species.\n*gamma* can be a list of lists specifying the coupling matrix.\nReturns 1 on success.\n"},
	{"aappp_reset_measurement", (PyCFunction)VM_reset_measurement, METH_VARARGS | METH_KEYWORDS, "aappp_reset_measurement(simulation, binnum_theta, binnum_neighbors)\n--\n\nForgets all measured data such as moments of polar order, etc. and histograms.\nThe number of bins of orientation histograms (*binnum_theta*) or number of neighbor histograms (*binnum_neighbors*) for following measurements can be set optionally by keyword.\nReturns 1 on success.\n"},
	{"aappp_get_xythetalxly", (PyCFunction)VM_xythetalxly, METH_VARARGS, "aappp_get_xythetalxly(simulation)\n--\n\nreturns [x, y, theta, lx, ly],\nwhere x, y and theta are lists of length N containing x-position, y-position and orientation of all particles.\nlx and ly give the size of the simulation box in x- and y-direction.\n"},
	{"aappp_get_parameters", (PyCFunction)VM_print_parameters, METH_VARARGS, "aappp_get_parameters(simulation)\n--\n\nreturns: [v, eta, R, omega, lx, ly, gamma, dt, N, kn, order, binnum_theta, binnum_neighbors, boundary_x, boundary_y],\nwhere *v* is particle speed, *eta* is noise strength, *R* is interaction radius, *omega* is chirality, *lx* is simulation box size in x-direction, *ly* is simulation box size in y-direction, *gamma* is the interaction strength, *dt* is step size, *N* is particle number, *kn* is the number of distinct interacting neighbors in metric free/topological models, *order* is the interaction order for overdamped Langevin models, *binnum_theta* is the number of bins that are used for orientation histograms, *binnum_neighbors* is the number of bins used in number of neighbor histograms, *boundary_x* is a flag for boundary conditions in x-direction (0->periodic boundary conditions, 1->reflecting boundary conditions), *boundary_y* is a flag for boundary conditions in y-direction (0->periodic boundary conditions, 1->reflecting boundary conditions).\nNot all parameters are used in each model (dynamics).\nFor example, *dt*, *order*, *gamma* are not used in Vicsek typ models (VM, NVM, mfVM, mfNVM), where the time evolution is a discrete map.\n*kn* is not used in metric models (VM, NVM, additiveL, nonadditive_L)\nand *R* is not used in metric free/topological models (mfVM, mfNVM, mfL).\nFor periodic boundary conditions (here explained in x-direction) particles that leave x in [0, lx] are put to x->x+lx or x->x-lx such that they are back inside the simulation box.\nFor calculating the distance between particles i and j, the minimum of |x_i-x_j|, |x_i-x_j-lx| and |x_i-x_j+lx| is used.\nFor reflecting boundary conditions (here explained in x-direction) particles that leave x in [0, lx] are put to either -x or -x + 2*lx such that they are back inside the simulation box [0, lx],\n when this happens, the particle orientation is reflected to theta->-theta+pi.\nParticles interact with all real particles but also with image particles that are the mirror image under the above reflection (of position and orientation).\nIf there is more than one particle species, then *v*, *eta*, *omega*, *N* are lists containing the corresponding quantities of all particle species and *gamma* is a list of lists (coupling matrix).\n"},
	{"aappp_get_weights", (PyCFunction)aappp_weights, METH_VARARGS, "aappp_get_weights(simulation)\n--\n\nreturns list of weights that are used for nonadditive overdamped Langevin dynamics nonadditiveL_update_timesteps, see >>>help(nonadditiveL_update_timesteps)<<<.\nif an empty list is returned, the standard weight function weight(n)=1/(n+1) is used.\n"},
	{"aappp_get_results", (PyCFunction)VM_measurement_results, METH_VARARGS, "aappp_get_results(simulation)\n--\n\nif there is only one particle species it returns: [[ptheta, theta], [pn, n], polar, nematic, neighbors],\nwhere *ptheta* is the normalized histogram of particle orientations and the list *theta* is composed of the centers of the bins of the corresponding histogram.\nThe list *pn* is the normalized histogram of the number of neighbor distribution and the list *n* gives the corresponding histogram bins 0, 1, 2, ..., n_max.\nEvents of particles that have more than n_max neighbors are not recorded.\nThus *pn* is not normalized to one as soon as such events are missed.\nThe list *polar* contains the polar order parameter from the last time measured (after penultimate time step) in the first entry and the first four moments of the polar order parameter in the next four entries.\nThe list *nematics* contains the nematic order parameter in the first entry and the first four moments of the nematic order parameter in the next four entries.\nThe list *neighbors* contains the first four moments of the number of neighbor distribution.\nIf there is more than one particle species, the function returns:\[[[ptheta1, theta1], [ptheta2, theta2], ...], [[pn1, n1], [pn2, n2], ...], [polar1, polar2, ...], [nematic1, nematic2, ...], [neighbors1, neighbors2, ...]].\nThose quantities are analogous to the single species case, but with the same results for each species."},
	{"VM_update_timesteps", (PyCFunction)VM_update_timesteps, METH_VARARGS, "VM_update_timesteps(simulation, timesteps)\n--\n\nPerforms *timesteps* timesteps of the Vicsek model dynamics with parameters specified before.\nReturns polar order parameter after penultimate time step.\nBoth arguments must be given.\nThe first argument *simulation* must specify a simulation state that was created earlier by the aappp_init function, see >>>help(aappp_init)<<< or loaded via aappp_load, see >>>help(aappp_load)<<<.\nThe second argument *timesteps* specifies the number of timesteps that are performed. The final state after those steps is stored in the *simulation* data structure and can be accessed by the aappp_get_xythetalxly function, see >>>help(aappp_get_xythetalxly)<<<.\n\n The Vicsek dynamics that is performed consists of two parts: 1) alignment and 2) streaming. In each step 1) is performed first and 2) afterwards. In 1) for each particle i (with position x_i, y_i, orientation theta_i), all neighboring particles j are found such that (x_i-x_j)**2 + (y_i-y_j)**2<R**2. Note that this holdes true for particle i itself. Then a new angle is calculated as:\ntheta_new_i=arg(sum_{j is neighbor of i} exp(imaginary unit*theta_i)).\nThe new angle is disturbed by a random number:\n theta_new_i = theta_new_i + xi_i,\nwhere xi_i are independently drawn from a uniform distribution on [-eta*pi, eta*pi]. After updating all orientations 1) is finished. In 2) all positions are updated according to:\nx_i=x_i+v*cos(theta_new_i)\ny_i=y_i+v*sin(theta_new_i)\nand the old orientations can be forgotten:\ntheta_i=theta_new_i.\n"},
	{"VM_measurement_timesteps", (PyCFunction)VM_measurement_timesteps, METH_VARARGS, "VM_measurement_timesteps(simulation, timesteps)\n--\n\nfollows the Vicsek model dynamics for *timesteps* timesteps as VM_update_timesteps, see >>>help(VM_update_timesteps)<<<, however, certain quantities are measured in each step.\nThe measured quantities are:\npolar order parameter: p=|sum_{i}exp(imaginary unit*theta_i)|/N,\nnematic order paremeter:q=|sum_{i}exp(imaginary unit*2*theta_i)|/N,\nthe number of neighbors: for each particle the number of (different) neighbors is counted.\nfor polar and nematic order parameter the first to fourth moment is calculated from the time series, for the number of neighbors these four moments are also calculated averaging over both, all particles and the time series (useful to calculate Binder parameter).\nThe averaging is performed over *timesteps* timesteps counting the initial configuration (before the first update), but not the last configuration (after the last update) in order to avoid to count certain states twice when calling VM_measurement_timesteps multiple times.\nFor the particle orientations (angles) and the number of neighbors, additionally to the first four moments, also a histogram is recorded.\nThe angles are defined in such a way that they are in the interval [-pi, pi].\nThe number of neighbors can in theory be any number from 0 to N-1.\nThe number of bins for the histograms can be specified in aappp_init, see >>>help(aappp_init)<<< or aappp_reset_measurement, see >>>help(aappp_reset_measurement)<<<.\nIf there is more than one particle species, polar and nematic order, number of neighbors and corresponding histogram are measured for each particle species separately.\nThe results of the measurement are stored in the *simulation* data structure.\nThe measurement results can be accessed by calling aappp_get_results, see >>>help(aappp_get_results)<<<.\nCalling VM_measurement_timesteps mutliple times will average over all timesteps of the multiple calls (except for the very last step).\nNormalization of the moments and histograms is done only when calling aappp_get_results, see >>>help(aappp_get_results)<<<.\nThus calling VM_measurement_timesteps(simulation, timesteps) twice introduces no rounding errors compared to a single call of VM_measurement_timesteps(simulation, 2*timesteps).\nAfter a call of VM_measurement_timesteps one can delete the measurement results (in order to perform a new measurement later) if desired, by calling aappp_reset_measurement, see >>>help(aappp_reset_measurement)<<<.\n"},
	{"NVM_update_timesteps", (PyCFunction)NVM_update_timesteps, METH_VARARGS, "NVM_update_timesteps(simulation, timesteps)\n--\n\nas VM_update_timesteps, see >>>help(VM_update_timesteps)<<<, but with different alignment rule 1):\nfirst, the new angle is determined as before:\ntheta_new_i=arg(sum_{j is neighbor of i} exp(imaginary unit*theta_i))\nbut then the angle is rotated by pi, if the previous orientation is closer to this rotated angle than to theta_new_i:\ntheta_new_i=theta_new_i*sign(sin(theta_new_i)*sin(theta_i)+cos(theta_new_i)*cos(theta_i)).\n"},
	{"NVM_measurement_timesteps", (PyCFunction)NVM_measurement_timesteps, METH_VARARGS, "NVM_measurement_timesteps(simulation, timesteps)\n--\n\nperforms *timesteps* timesteps with the same dynamics as NVM_update_timesteps, see >>>help(NVM_update_timesteps)<<< and measures the same quantities as VM_measurement_timesteps, see >>>help(VM_measurement_timesteps)<<<.\n"},
	{"mfVM_update_timesteps", (PyCFunction)mfVM_update_timesteps, METH_VARARGS, "mfVM_update_timesteps(simulation, timesteps)\n--\n\nas VM_update_timesteps, see >>>help(VM_update_timesteps)<<<, but with a different definition of neighborhoods (called metric free or topological).\nHere particle i is considered as neighbor of particle i (as before).\nThe other neighbors are the *kn* particles that are closest to particle i.\nThus the parameter *kn* should have values 1, 2, 3, ..., see >>>help(aappp_init)<<< for initialization of the parameters.\n"},
	{"mfVM_measurement_timesteps", (PyCFunction)mfVM_measurement_timesteps, METH_VARARGS, "mfVM_measurement_timesteps(simulation, timesteps)\n--\n\nperforms *timesteps* timesteps with the same dynamics as mfVM_update_timesteps, see >>>help(mfVM_update_timesteps)<<< and measures the same quantities as VM_measurement_timesteps, see >>>help(VM_measurement_timesteps)<<<.\nNote that the number of neighbors measurements are pointless in this model because the number of neighbors is alsways the same.\nHowever, the number of neighbor measurements are nevertheless performed in the same way.\n"},
	{"mfNVM_update_timesteps", (PyCFunction)mfNVM_update_timesteps, METH_VARARGS, "mfNVM_update_timesteps(simulation, timesteps)\n--\n\nas VM_update_timesteps, see >>>help(VM_update_timesteps)<<<, but with metric free/topological definition of neighborhoods as in mfVM_update_timesteps, see >>>help(mfVM_update_timesteps)<<<, and with a nematic interaction rule as in NVM_update_timesteps, see >>>help(NVM_update_timesteps)<<<.\n"},
	{"mfNVM_measurement_timesteps", (PyCFunction)mfNVM_measurement_timesteps, METH_VARARGS, "mfNVM_measurement_timesteps(simulation, timesteps)\n--\n\nperforms *timesteps* timesteps with the same dynamics as mfNVM_update_timesteps, see >>>help(mfNVM_update_timesteps)<<< and measures the same quantities as VM_measurement_timesteps, see >>>help(VM_measurement_timesteps)<<<.\nNote that the number of neighbors measurements are pointless in this model because the number of neighbors is alsways the same.\nHowever, the number of neighbor measurements are nevertheless performed in the same way.\n"},
	{"additiveL_update_timesteps", (PyCFunction)additiveL_update_timesteps, METH_VARARGS, "additiveL_update_timesteps(simulation, timesteps)\n--\n\nperforms *timesteps* timesteps for the following overdamped Langevin dynamics:\nd/dt x_i=v*cos(phi_i)\nd/dt y_i=v*sin(phi_i)\nd/dt phi_i=gamma*sum_{j is neighbor of i}sin((phi_j-phi_i)*order) + omega + eta*xi_i\n\ntwo particles are considered to be neighbors if their distance is smaller than *R* as in the Vicsek model, see >>>help(VM_update_timesteps)<<<.\nIf there is more than one particle species, v, omega and eta can be different for each species and gamma can depend on the species of i and j (coupling matrix), see >>>help(aappp_init)<<< for parameter definitions.\nThe overdamped Langevin equations are integrated with Euler-Maruyama-scheme with step size *dt*.\n*timesteps* denotes the number of timesteps performed, such that the real time of time evolution is *timesteps*  *  *dt*.\nReturns polar order parameter before last step.\n"},
	{"additiveL_measurement_timesteps", (PyCFunction)additiveL_measurement_timesteps, METH_VARARGS, "additiveL_measurement_timesteps(simulation, timesteps)\n--\n\nperforms *timesteps* timesteps as additiveL_update_timesteps, see >>>help(additiveL_update_timesteps)<<< and measures observables as VM_measurement_timesteps, see >>>help(VM_measurement_timesteps)<<<.\n"},
	{"nonadditiveL_update_timesteps", (PyCFunction)nonadditiveL_update_timesteps, METH_VARARGS, "nonadditiveL_update_timesteps(simulation, timesteps)\n--\n\ndynamics is iterated as in additiveL_update_timesteps, see >>>help(additiveL_update_timesteps)<<<.\nhere, however, the equation of motion of the orientation is:\nd/dt phi_i=gamma*weight(n)*sum_{j is neighbor of i}sin((phi_j-phi_i)*order) + omega + eta*xi_i,\nwhere n is the number of distinct neighbors of particle i (not counting i itself).\nif no weight function weight(n) is specified, the function weight(n)=1/(n+1) is used for all occuring values of n.\nthe weight function weight(n) can be specified using aappp_init or aappp_reset_parameters, see >>>help(aappp_init)<<< or >>>help(aappp_reset_parameters)<<<.\nif a weight function is specified, also a value for *weight_vector_length* must be specified.\nin that case the weight function is used for all values of n in the range 0, 1, 2, ..., *weight_vector_length*-1\nfor n>=weight_vector_length, the weight weight(*weight_vector_length*-1) is used.\n"},
	{"nonadditiveL_measurement_timesteps", (PyCFunction)nonadditiveL_measurement_timesteps, METH_VARARGS, "additiveL_measurement_timesteps(simulation, timesteps)\n--\n\nperforms *timesteps* timesteps as nonadditiveL_update_timesteps, see >>>help(nonadditiveL_update_timesteps)<<< and measures observables as VM_measurement_timesteps, see >>>help(VM_measurement_timesteps)<<<.\n"},
	{"mfL_update_timesteps", (PyCFunction)mfL_update_timesteps, METH_VARARGS, "mfL_update_timesteps(simulation, timesteps)\n--\n\nevolves the system for *timesteps* timesteps as additiveL_update_timesteps, see >>>help(additiveL_update_timesteps)<<<, but with metric free/topological definition of neighborhoods, as in mfVM_update_timesteps, see >>>help(mfVM_update_timesteps)<<<.\n"},
	{"mfL_measurement_timesteps", (PyCFunction)mfL_measurement_timesteps, METH_VARARGS, "additiveL_measurement_timesteps(simulation, timesteps)\n--\n\nperforms *timesteps* timesteps as mfL_update_timesteps, see >>>help(mfL_update_timesteps)<<< and measures observables as VM_measurement_timesteps, see >>>help(VM_measurement_timesteps)<<<.\n"},
	{NULL, NULL, 0, NULL}
};

static struct PyModuleDef aappp =
{
    PyModuleDef_HEAD_INIT,
    "aappp", /* name of module */
    "aappp aligning active particles py package version 1.2\n\n---VERSION INFO---\nBefore version 1.2 the implementation of reflecting boundary conditions for models mfL, additiveL and nonadditiveL was not as intended.\nParticle orientations have not been reflected at the boundary for the aforementioned models (reflecting bc were correct for VM, NVM, mfVM, mfNVM).\nThis was corrected with this version.\nFrom version 1.0 to version 1.1 the measurement of moments of the number of neighbors was modified.\nIn version 1.0 the moments of the neighbor number have been obtained averaging over all particles and all measurement time steps.\nIn that way the moments could also be calculated from the number of neighbor histogram (when the histogram is large enough to not miss any events) yielding the exact same result.\nIn version 1.1 at each time step the average number of neighbors is calculated (averaging over all particles at a given time).\nIn this way, the average number of neighbors is an observable.\nFrom this ensemble average number, moments are calculated averaging over different time steps.\nThis change does not affect the first moment, because the two averages commute for the first moment.\nHowever, the second to fourth momemt are different.\nThe number of neighbor histogram was not changed.\nThus, the histogram is still taken over all particles and time steps (without taking the average before).\nIt is possible to use version 1.1 to analyze data produced with version 1.0.\nIn this case a warning is shown, however all data can be accesses in the same way as using version 1.0.\nIf moments of the number of neighbors are not considered, versions 1.1 and 1.0 are equivalent.\n------\n\nThis package provides functions for agent-based simulations of various models for aligning active particles in two dimensions.\nSupported models are: Vicsek model, nematic Vicsek model, metric free (topological) Vicsek model, metric free (topological) nematic Vicsek model and overdamped Langevin models of the type:\nd/dt x_i= v*cos(phi_i),\nd/dt y_i= v*sin(phi_i),\nd/dt phi_i=sum_{j is neighbor of i}*gamma*sin((phi_j-phi_i)*order)+omega+eta*xi_i,\n where the neighborhood is either defined by an interaction distance, or metric free (topological) interactions are considered, where each particle interacts with the kn nearest neighbors.\nIn each model, chirality can be added (omega) and one can simulate multiple particle species, where each species can have a different value for noise strength, velocity, chirality, coupling (=coupling matrix).\nUse function aappp_init to initialize a simulation, see >>>help(aappp_init)<<<.\nUse function ***_update_timesteps to evolve the system in time.\n*** can be VM, NVM, mfVM, mfNVM for Vicsek models and additiveL, nonadditiveL or mfL for overdamped Langevin models, see >>>help(***_update_timesteps)<<< for exact model definitions.\nUse function ***_measurement_timesteps to evolve the system in time and measure polar order (first 4 moments), nematic order (first 4 moments), number of neighbors (first four moments), histogram of the orientation and a histogram of the number of neighbors, see >>>help(***_measurement_timesteps)<<<.\nUse function aappp_get_results to obtain the results of the measurements, see >>>help(aappp_get_results)<<<.\nUse function aappp_get_xythetalxly to get positions and orientation of all particles and simulation box size, see >>>help(aappp_get_xythetalxly)<<<.\nUse function aappp_get_parameters to obtain the simulation parameters, see >>>help(aappp_get_parameters)<<<.\nUse function aappp_set_state to manipulate the positions and orientations of the particles and/or resize the simulation box, see >>>help(aappp_set_state)<<<.\nUse function aappp_reset_parameters to change some simulation parameters, see >>>help(aappp_reset_parameters)<<<.\nUse function aappp_reset_measurement to clear all measurements and possibly change the number of bins for histograms, see >>>help(aappp_reset_measurement)<<<.\nUse function aappp_save to save the state of the simulation (including measurement resuts and state of the pseudo random number generator) to a file, see >>>help(aappp_save)<<<.\nUse function aappp_load to load a previous simulation, see >>>help(aappp_load)<<<.\nUse function aappp_get_weights to obtain a list of weights used for nonadditive overdamped Langevin model, see >>>help(aappp_get_weights)<<<.\nParallel computation is not supported by now. However, performence and memory usage should be good.\n\nEXAMPLE of use\n\nmysim=aappp.aappp_init(N=10000, eta=0.2, v=1., R=1., Lx=40., Ly=80.)\naappp.VM_update_timesteps(mysim, 10000)\naappp.aappp_measurement_timesteps(mysim, 1000)\naappp.aappp_get_results(mysim)\nx, y, theta, lx, ly=aappp.aappp_get_xythetalxly(mysim)\naappp.aappp_save(mysim, 'mydata.dat')\naappp.aappp_free(mysim)\n\n",          /* module documentation, may be NULL */
    -1,          /* size of per-interpreter state of the module, or -1 if the module keeps state in global variables. */
    aappp_methods
};

PyMODINIT_FUNC PyInit_aappp(void){
	return PyModule_Create(&aappp);
}

