/*
 * =====================================================================================
 *
 *       Filename:  c_functions_aappp.c
 *
 *    Description:  
 *
 *        Version:  1.2
 *        Created:  11/29/2022
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

#include "c_functions_aappp.h"

//returning the minimum of two integers
gint32 min(gint32 i, gint32 j)
{
	if (i<j)
		return i;
	return j;
}

//returning a zero mean, unit variance Gaussian distributed random number
gdouble gauss(struct VM_simulation* simulation)
{
	if (simulation->state->box_muller==3)
	{
		gdouble u1=g_rand_double(simulation->state->prng);
		gdouble u2=g_rand_double(simulation->state->prng);
		simulation->state->gauss1=sqrt(-2.0*log(u1))*cos(2.0*M_PI*u2);
		simulation->state->gauss2=sqrt(-2.0*log(u1))*sin(2.0*M_PI*u2);
		simulation->state->box_muller=1;
	}
	if (simulation->state->box_muller==1)
	{
		simulation->state->box_muller=2;
		return simulation->state->gauss1;
	}
	if (simulation->state->box_muller==2)
	{
		simulation->state->box_muller=3;
	}
	return simulation->state->gauss2;
}

//returns the index of the box in which a particle with positions (xposition, yposition) is located for metric models (where box length is approximately R)
gint32 VM_get_box_index(gdouble xposition, gdouble yposition, struct VM_parameters* parameters)
{
	return ( (gint32)(xposition/parameters->box_size_x)) + ( (gint32)(yposition/parameters->box_size_y))*parameters->box_number_x;
}

//For system of multiple particle species: returns the species to which the particle with index index belongs to
gint32 get_species(guint64 index, struct VM_parameters * parameters)
{
	gint32 result=0;
	guint64 pnumber=0;
	while (result<parameters->particle_species)
	{
		pnumber+=*(parameters->species_particle_number+result);
		if (index<pnumber)
			break;
		result+=1;
	}
	return result;
}

//empties the data structure that contains lists of neighbors for metric free/topological models
void mfVM_clean_neighborhood(struct VM_simulation * simulation)
{
	gint32 i;
	gdouble lmax=simulation->state->length_x; //maximum of Lx and Ly
	if (simulation->state->length_y>lmax)
		lmax=simulation->state->length_y;
	for (i=0; i<simulation->parameters->mf_kn; i++)
	{
		*((simulation->state->nh).neighbor_indexes+i)=-1;
		*((simulation->state->nh).distances+i)=lmax;
		*((simulation->state->nh).mirror_flag+i)=0;
	}
	(simulation->state->nh).max_distance_checked=0.0;
	return;
}

//adds a particle to a neighbor list for metric free models
void mfVM_add_to_neighborhood(guint64 index, gdouble distance, gint32 flag, struct VM_simulation * simulation)
{
	gint32 i;
	guint64 upper_index;
	guint64 lower_index;
	gdouble upper_distance;
	gdouble lower_distance;
	gint32 upper_flag;
	gint32 lower_flag;
	//if distance is smaller than distance to the last particle in the list
	if (distance<*(simulation->state->nh.distances+simulation->parameters->mf_kn-1))
	{
		//append new neighbor at the end of the neighbor list
		*(simulation->state->nh.distances+simulation->parameters->mf_kn-1)=distance;
		*(simulation->state->nh.neighbor_indexes+simulation->parameters->mf_kn-1)=index;
		*(simulation->state->nh.mirror_flag+simulation->parameters->mf_kn-1)=flag;
		//sort neighbors: the one with largest distance last
		i=simulation->parameters->mf_kn-1;
		while (i>0)
		{
			upper_index=*(simulation->state->nh.neighbor_indexes+i);
			lower_index=*(simulation->state->nh.neighbor_indexes+i-1);
			upper_distance=*(simulation->state->nh.distances+i);
			lower_distance=*(simulation->state->nh.distances+i-1);
			upper_flag=*(simulation->state->nh.mirror_flag+i);
			lower_flag=*(simulation->state->nh.mirror_flag+i-1);
			if (lower_distance<upper_distance)
				break;
			else
			{
				*(simulation->state->nh.neighbor_indexes+i)=lower_index;
				*(simulation->state->nh.neighbor_indexes+i-1)=upper_index;
				*(simulation->state->nh.distances+i)=lower_distance;
				*(simulation->state->nh.distances+i-1)=upper_distance;
				*(simulation->state->nh.mirror_flag+i)=lower_flag;
				*(simulation->state->nh.mirror_flag+i-1)=upper_flag;
				i=i-1;
			}
		}
	}
	return;
}

//returns the index of the box in which a particle with positions (xposition, yposition) is located for metric free/topological models (where box length is approximately R_eff for some effective radius R_eff)
gint32 mfVM_get_box_index(gdouble xposition, gdouble yposition, struct VM_parameters* parameters)
{
	return ( (gint32)(xposition/parameters->mf_box_size_x)) + ( (gint32)(yposition/parameters->mf_box_size_y))*parameters->mf_box_number_x;
}

//add particle to a specific box for metric model
void VM_append_particle_to_box(guint64 particle_index, struct VM_simulation * simulation)
{
	//find index of the box
	gint32 box_index=VM_get_box_index((simulation->state->particles+particle_index)->x, (simulation->state->particles+particle_index)->y, simulation->parameters);
	//actually append the particle to the box
	(simulation->state->particles+particle_index)->next_particle=*(simulation->parameters->box+box_index);
	*(simulation->parameters->box+box_index)=simulation->state->particles+particle_index;
}

//add particle to a specific box for metric free/topological model
void mfVM_append_particle_to_box(guint64 particle_index, struct VM_simulation * simulation)
{
	//find index of the box
	gint32 box_index=mfVM_get_box_index((simulation->state->particles+particle_index)->x, (simulation->state->particles+particle_index)->y, simulation->parameters);
	//actually append the particle to the box
	(simulation->state->particles+particle_index)->next_particle=*(simulation->parameters->mf_box+box_index);
	*(simulation->parameters->mf_box+box_index)=simulation->state->particles+particle_index;
}

//sort all particles in boxes depending on their position for metric models
void VM_fill_boxes(struct VM_simulation * simulation)
{
	guint64 i;
	for(i=0; i<simulation->state->particle_number; i++)
		VM_append_particle_to_box(i, simulation);
}

//sort all particles in boxes depending on their position for metric free/topological models
void mfVM_fill_boxes(struct VM_simulation * simulation)
{
	guint64 i;
	for(i=0; i<simulation->state->particle_number; i++)
		mfVM_append_particle_to_box(i, simulation);
}

//delete information about particles contained in boxes for metric models
void VM_clean_boxes(struct VM_parameters* parameters)
{
	gint32 i;
	for (i=0; i<parameters->box_number_x*parameters->box_number_y; i++)
		*(parameters->box+i)=NULL;
}

//delete information about particles contained in boxes for metric free/topological models
void mfVM_clean_boxes(struct VM_parameters* parameters)
{
	gint32 i;
	for (i=0; i<parameters->mf_box_number_x*parameters->mf_box_number_y; i++)
		*(parameters->mf_box+i)=NULL;
}

//calculate the interactions of particle 'index' with all particles in a given box for Vicsek model
void VM_interaction_with_neighbors_in_box(gint32 box_index, struct VM_simulation * simulation, guint64 index, gdouble * velocity_x, gdouble * velocity_y, gint32 * neighbor_n)
{
	gdouble distance, distance_squared;
	gdouble xposition_focus=(simulation->state->particles+index)->x;
	gdouble yposition_focus=(simulation->state->particles+index)->y;
	struct VM_particle* box= *(simulation->parameters->box+box_index);
	while (box!=NULL)
	{
		distance_squared=0.0;
		distance=xposition_focus-(simulation->state->particles+box->index)->x;
		if (distance>0.5*simulation->state->length_x)
				distance -= simulation->state->length_x;
		if (distance<-0.5*simulation->state->length_x)
				distance += simulation->state->length_x;
		distance_squared += distance*distance;
		distance=yposition_focus-(simulation->state->particles+box->index)->y;
		if (distance>0.5*simulation->state->length_y)
				distance -= simulation->state->length_y;
		if (distance<-0.5*simulation->state->length_y)
				distance += simulation->state->length_y;
		distance_squared += distance*distance;
		if (distance_squared<((simulation->parameters->interaction_radius)*(simulation->parameters->interaction_radius)))
		{
			*velocity_x += cos((simulation->state->particles+box->index)->theta);
			*velocity_y += sin((simulation->state->particles+box->index)->theta);
			if (index!=box->index)
				*neighbor_n+=1;
		}
		box=box->next_particle;
	}
	return;
}

//calculate the interactions of particle 'index' with the mirror immage of all particles in a given box under reflection with x=0,Lx axis for Vicsek model
void VM_interaction_with_neighbors_in_box_mirrorx(gint32 box_index, struct VM_simulation * simulation, guint64 index, gdouble * velocity_x, gdouble * velocity_y, gint32 * neighbor_n)
{
	gdouble distance, distance_squared;
	gdouble xposition_focus=(simulation->state->particles+index)->x;
	gdouble yposition_focus=(simulation->state->particles+index)->y;
	struct VM_particle* box= *(simulation->parameters->box+box_index);
	while (box!=NULL)
	{
		distance_squared=0.0;
		distance=xposition_focus+(simulation->state->particles+box->index)->x;
		while (distance>0.5*simulation->state->length_x)
				distance -= simulation->state->length_x;
		while (distance<-0.5*simulation->state->length_x)
				distance += simulation->state->length_x;
		distance_squared += distance*distance;
		distance=yposition_focus-(simulation->state->particles+box->index)->y;
		if (distance>0.5*simulation->state->length_y)
				distance -= simulation->state->length_y;
		if (distance<-0.5*simulation->state->length_y)
				distance += simulation->state->length_y;
		distance_squared += distance*distance;
		if (distance_squared<((simulation->parameters->interaction_radius)*(simulation->parameters->interaction_radius)))
		{
			*velocity_x += cos(M_PI-(simulation->state->particles+box->index)->theta);
			*velocity_y += sin(M_PI-(simulation->state->particles+box->index)->theta);
			*neighbor_n+=1;
		}
		box=box->next_particle;
	}
	return;
}

//calculate the interactions of particle 'index' with the mirror immage of all particles in a given box under reflection with y=0,Ly axis for Vicsek model
void VM_interaction_with_neighbors_in_box_mirrory(gint32 box_index, struct VM_simulation * simulation, guint64 index, gdouble * velocity_x, gdouble * velocity_y, gint32 * neighbor_n)
{
	gdouble distance, distance_squared;
	gdouble xposition_focus=(simulation->state->particles+index)->x;
	gdouble yposition_focus=(simulation->state->particles+index)->y;
	struct VM_particle* box= *(simulation->parameters->box+box_index);
	while (box!=NULL)
	{
		distance_squared=0.0;
		distance=xposition_focus-(simulation->state->particles+box->index)->x;
		if (distance>0.5*simulation->state->length_x)
				distance -= simulation->state->length_x;
		if (distance<-0.5*simulation->state->length_x)
				distance += simulation->state->length_x;
		distance_squared += distance*distance;
		distance=yposition_focus+(simulation->state->particles+box->index)->y;
		while (distance>0.5*simulation->state->length_y)
				distance -= simulation->state->length_y;
		while (distance<-0.5*simulation->state->length_y)
				distance += simulation->state->length_y;
		distance_squared += distance*distance;
		if (distance_squared<((simulation->parameters->interaction_radius)*(simulation->parameters->interaction_radius)))
		{
			*velocity_x += cos(-(simulation->state->particles+box->index)->theta);
			*velocity_y += sin(-(simulation->state->particles+box->index)->theta);
			*neighbor_n+=1;
		}
		box=box->next_particle;
	}
	return;
}

//calculate the interactions of particle 'index' with the mirror immage of all particles in a given box under reflection with both axis for Vicsek model
void VM_interaction_with_neighbors_in_box_mirrorxy(gint32 box_index, struct VM_simulation * simulation, guint64 index, gdouble * velocity_x, gdouble * velocity_y, gint32 * neighbor_n)
{
	gdouble distance, distance_squared;
	gdouble xposition_focus=(simulation->state->particles+index)->x;
	gdouble yposition_focus=(simulation->state->particles+index)->y;
	struct VM_particle* box= *(simulation->parameters->box+box_index);
	while (box!=NULL)
	{
		distance_squared=0.0;
		distance=xposition_focus+(simulation->state->particles+box->index)->x;
		while (distance>0.5*simulation->state->length_x)
				distance -= simulation->state->length_x;
		while (distance<-0.5*simulation->state->length_x)
				distance += simulation->state->length_x;
		distance_squared += distance*distance;
		distance=yposition_focus+(simulation->state->particles+box->index)->y;
		while (distance>0.5*simulation->state->length_y)
				distance -= simulation->state->length_y;
		while (distance<-0.5*simulation->state->length_y)
				distance += simulation->state->length_y;
		distance_squared += distance*distance;
		if (distance_squared<((simulation->parameters->interaction_radius)*(simulation->parameters->interaction_radius)))
		{
			*velocity_x += cos(M_PI+(simulation->state->particles+box->index)->theta);
			*velocity_y += sin(M_PI+(simulation->state->particles+box->index)->theta);
			*neighbor_n+=1;
		}
		box=box->next_particle;
	}
	return;
}

//calculate the interactions of particle 'index' with all particles in a given box for nematic Vicsek model
void NVM_interaction_with_neighbors_in_box(gint32 box_index, struct VM_simulation * simulation, guint64 index, gdouble * velocity_x, gdouble * velocity_y, gint32 * neighbor_n)
{
	gdouble distance, distance_squared;
	gdouble xposition_focus=(simulation->state->particles+index)->x;
	gdouble yposition_focus=(simulation->state->particles+index)->y;
	struct VM_particle* box= *(simulation->parameters->box+box_index);
	while (box!=NULL)
	{
		distance_squared=0.0;
		distance=xposition_focus-(simulation->state->particles+box->index)->x;
		if (distance>0.5*simulation->state->length_x)
				distance -= simulation->state->length_x;
		if (distance<-0.5*simulation->state->length_x)
				distance += simulation->state->length_x;
		distance_squared += distance*distance;
		distance=yposition_focus-(simulation->state->particles+box->index)->y;
		if (distance>0.5*simulation->state->length_y)
				distance -= simulation->state->length_y;
		if (distance<-0.5*simulation->state->length_y)
				distance += simulation->state->length_y;
		distance_squared += distance*distance;
		if (distance_squared<((simulation->parameters->interaction_radius)*(simulation->parameters->interaction_radius)))
		{
			if (cos((simulation->state->particles+box->index)->theta-(simulation->state->particles+index)->theta)>0.0)
			{
				*velocity_x += cos((simulation->state->particles+box->index)->theta);
				*velocity_y += sin((simulation->state->particles+box->index)->theta);
			}
			else
			{
				*velocity_x -= cos((simulation->state->particles+box->index)->theta);
				*velocity_y -= sin((simulation->state->particles+box->index)->theta);
			}
			if (index!=box->index)
				*neighbor_n+=1;
		}
		box=box->next_particle;
	}
	return;
}

//calculate the interactions of particle 'index' with the mirror immage of all particles in a given box under reflection with x=0,Lx axis for nematic Vicsek model
void NVM_interaction_with_neighbors_in_box_mirrorx(gint32 box_index, struct VM_simulation * simulation, guint64 index, gdouble * velocity_x, gdouble * velocity_y, gint32 * neighbor_n)
{
	gdouble distance, distance_squared;
	gdouble xposition_focus=(simulation->state->particles+index)->x;
	gdouble yposition_focus=(simulation->state->particles+index)->y;
	struct VM_particle* box= *(simulation->parameters->box+box_index);
	while (box!=NULL)
	{
		distance_squared=0.0;
		distance=xposition_focus+(simulation->state->particles+box->index)->x;
		while (distance>0.5*simulation->state->length_x)
				distance -= simulation->state->length_x;
		while (distance<-0.5*simulation->state->length_x)
				distance += simulation->state->length_x;
		distance_squared += distance*distance;
		distance=yposition_focus-(simulation->state->particles+box->index)->y;
		if (distance>0.5*simulation->state->length_y)
				distance -= simulation->state->length_y;
		if (distance<-0.5*simulation->state->length_y)
				distance += simulation->state->length_y;
		distance_squared += distance*distance;
		if (distance_squared<((simulation->parameters->interaction_radius)*(simulation->parameters->interaction_radius)))
		{
			if (cos(M_PI-(simulation->state->particles+box->index)->theta-(simulation->state->particles+index)->theta)>0.0)
			{
				*velocity_x += cos(M_PI-(simulation->state->particles+box->index)->theta);
				*velocity_y += sin(M_PI-(simulation->state->particles+box->index)->theta);
			}
			else
			{
				*velocity_x -= cos(M_PI-(simulation->state->particles+box->index)->theta);
				*velocity_y -= sin(M_PI-(simulation->state->particles+box->index)->theta);
			}
			*neighbor_n+=1;
		}
		box=box->next_particle;
	}
	return;
}

//calculate the interactions of particle 'index' with the mirror immage of all particles in a given box under reflection with y=0,Ly axis for nematic Vicsek model
void NVM_interaction_with_neighbors_in_box_mirrory(gint32 box_index, struct VM_simulation * simulation, guint64 index, gdouble * velocity_x, gdouble * velocity_y, gint32 * neighbor_n)
{
	gdouble distance, distance_squared;
	gdouble xposition_focus=(simulation->state->particles+index)->x;
	gdouble yposition_focus=(simulation->state->particles+index)->y;
	struct VM_particle* box= *(simulation->parameters->box+box_index);
	while (box!=NULL)
	{
		distance_squared=0.0;
		distance=xposition_focus-(simulation->state->particles+box->index)->x;
		if (distance>0.5*simulation->state->length_x)
				distance -= simulation->state->length_x;
		if (distance<-0.5*simulation->state->length_x)
				distance += simulation->state->length_x;
		distance_squared += distance*distance;
		distance=yposition_focus+(simulation->state->particles+box->index)->y;
		while (distance>0.5*simulation->state->length_y)
				distance -= simulation->state->length_y;
		while (distance<-0.5*simulation->state->length_y)
				distance += simulation->state->length_y;
		distance_squared += distance*distance;
		if (distance_squared<((simulation->parameters->interaction_radius)*(simulation->parameters->interaction_radius)))
		{
			if (cos(-(simulation->state->particles+box->index)->theta-(simulation->state->particles+index)->theta)>0.0)
			{
				*velocity_x += cos(-(simulation->state->particles+box->index)->theta);
				*velocity_y += sin(-(simulation->state->particles+box->index)->theta);
			}
			else
			{
				*velocity_x -= cos(-(simulation->state->particles+box->index)->theta);
				*velocity_y -= sin(-(simulation->state->particles+box->index)->theta);
			}
			*neighbor_n+=1;
		}
		box=box->next_particle;
	}
	return;
}

//calculate the interactions of particle 'index' with the mirror immage of all particles in a given box under reflection with both axis for nematic Vicsek model
void NVM_interaction_with_neighbors_in_box_mirrorxy(gint32 box_index, struct VM_simulation * simulation, guint64 index, gdouble * velocity_x, gdouble * velocity_y, gint32 * neighbor_n)
{
	gdouble distance, distance_squared;
	gdouble xposition_focus=(simulation->state->particles+index)->x;
	gdouble yposition_focus=(simulation->state->particles+index)->y;
	struct VM_particle* box= *(simulation->parameters->box+box_index);
	while (box!=NULL)
	{
		distance_squared=0.0;
		distance=xposition_focus+(simulation->state->particles+box->index)->x;
		while (distance>0.5*simulation->state->length_x)
				distance -= simulation->state->length_x;
		while (distance<-0.5*simulation->state->length_x)
				distance += simulation->state->length_x;
		distance_squared += distance*distance;
		distance=yposition_focus+(simulation->state->particles+box->index)->y;
		while (distance>0.5*simulation->state->length_y)
				distance -= simulation->state->length_y;
		while (distance<-0.5*simulation->state->length_y)
				distance += simulation->state->length_y;
		distance_squared += distance*distance;
		if (distance_squared<((simulation->parameters->interaction_radius)*(simulation->parameters->interaction_radius)))
		{
			if (cos(M_PI+(simulation->state->particles+box->index)->theta-(simulation->state->particles+index)->theta)>0.0)
			{
				*velocity_x += cos(M_PI+(simulation->state->particles+box->index)->theta);
				*velocity_y += sin(M_PI+(simulation->state->particles+box->index)->theta);
			}
			else
			{
				*velocity_x -= cos(M_PI+(simulation->state->particles+box->index)->theta);
				*velocity_y -= sin(M_PI+(simulation->state->particles+box->index)->theta);
			}
			*neighbor_n+=1;
		}
		box=box->next_particle;
	}
	return;
}

//calculate the interactions of particle 'index' with all particles in a given box for Langevin models
void L_interaction_with_neighbors_in_box(gint32 box_index, struct VM_simulation * simulation, guint64 index, gdouble * delta_v, gint32 * neighbor_n)
{
	gdouble distance, distance_squared;
	gdouble xposition_focus=(simulation->state->particles+index)->x;
	gdouble yposition_focus=(simulation->state->particles+index)->y;
	struct VM_particle* box= *(simulation->parameters->box+box_index);
	while (box!=NULL)
	{
		if (index!=box->index)
		{
			distance_squared=0.0;
			distance=xposition_focus-(simulation->state->particles+box->index)->x;
			if (distance>0.5*simulation->state->length_x)
					distance -= simulation->state->length_x;
			if (distance<-0.5*simulation->state->length_x)
					distance += simulation->state->length_x;
			distance_squared += distance*distance;
			distance=yposition_focus-(simulation->state->particles+box->index)->y;
			if (distance>0.5*simulation->state->length_y)
					distance -= simulation->state->length_y;
			if (distance<-0.5*simulation->state->length_y)
					distance += simulation->state->length_y;
			distance_squared += distance*distance;
			if (distance_squared<((simulation->parameters->interaction_radius)*(simulation->parameters->interaction_radius)))
			{
				if (simulation->parameters->particle_species==1)
					*delta_v+=**simulation->parameters->coupling*sin( ((simulation->state->particles+box->index)->theta - (simulation->state->particles+index)->theta)*simulation->parameters->interaction_order  );
				else
					*delta_v+=*(*(simulation->parameters->coupling+get_species(index, simulation->parameters))+get_species(box->index, simulation->parameters))*sin( ((simulation->state->particles+box->index)->theta - (simulation->state->particles+index)->theta)*simulation->parameters->interaction_order  );
				*neighbor_n+=1;
			}
		}
		box=box->next_particle;
	}
	return;
}

//calculate the interactions of particle 'index' with the mirror immage of all particles in a given box under reflection with x=0,Lx axis for Langevin models
void L_interaction_with_neighbors_in_box_mirrorx(gint32 box_index, struct VM_simulation * simulation, guint64 index, gdouble * delta_v, gint32 * neighbor_n)
{
	gdouble distance, distance_squared;
	gdouble xposition_focus=(simulation->state->particles+index)->x;
	gdouble yposition_focus=(simulation->state->particles+index)->y;
	struct VM_particle* box= *(simulation->parameters->box+box_index);
	while (box!=NULL)
	{
		distance_squared=0.0;
		distance=xposition_focus+(simulation->state->particles+box->index)->x;
		while (distance>0.5*simulation->state->length_x)
				distance -= simulation->state->length_x;
		while (distance<-0.5*simulation->state->length_x)
				distance += simulation->state->length_x;
		distance_squared += distance*distance;
		distance=yposition_focus-(simulation->state->particles+box->index)->y;
		if (distance>0.5*simulation->state->length_y)
				distance -= simulation->state->length_y;
		if (distance<-0.5*simulation->state->length_y)
				distance += simulation->state->length_y;
		distance_squared += distance*distance;
		if (distance_squared<((simulation->parameters->interaction_radius)*(simulation->parameters->interaction_radius)))
		{
			if (simulation->parameters->particle_species==1)
				*delta_v+=**simulation->parameters->coupling*sin( (M_PI-(simulation->state->particles+box->index)->theta - (simulation->state->particles+index)->theta)*simulation->parameters->interaction_order  );
			else
				*delta_v+=*(*(simulation->parameters->coupling+get_species(index, simulation->parameters))+get_species(box->index, simulation->parameters))*sin( (M_PI-(simulation->state->particles+box->index)->theta - (simulation->state->particles+index)->theta)*simulation->parameters->interaction_order  );
			*neighbor_n+=1;
		}
		box=box->next_particle;
	}
	return;
}

//calculate the interactions of particle 'index' with the mirror immage of all particles in a given box under reflection with y=0,Ly axis for Langevin models
void L_interaction_with_neighbors_in_box_mirrory(gint32 box_index, struct VM_simulation * simulation, guint64 index, gdouble * delta_v, gint32 * neighbor_n)
{
	gdouble distance, distance_squared;
	gdouble xposition_focus=(simulation->state->particles+index)->x;
	gdouble yposition_focus=(simulation->state->particles+index)->y;
	struct VM_particle* box= *(simulation->parameters->box+box_index);
	while (box!=NULL)
	{
		distance_squared=0.0;
		distance=xposition_focus-(simulation->state->particles+box->index)->x;
		if (distance>0.5*simulation->state->length_x)
				distance -= simulation->state->length_x;
		if (distance<-0.5*simulation->state->length_x)
				distance += simulation->state->length_x;
		distance_squared += distance*distance;
		distance=yposition_focus+(simulation->state->particles+box->index)->y;
		while (distance>0.5*simulation->state->length_y)
				distance -= simulation->state->length_y;
		while (distance<-0.5*simulation->state->length_y)
				distance += simulation->state->length_y;
		distance_squared += distance*distance;
		if (distance_squared<((simulation->parameters->interaction_radius)*(simulation->parameters->interaction_radius)))
		{
			if (simulation->parameters->particle_species==1)
				*delta_v+=**simulation->parameters->coupling*sin( (-(simulation->state->particles+box->index)->theta - (simulation->state->particles+index)->theta)*simulation->parameters->interaction_order  );
			else
				*delta_v+=*(*(simulation->parameters->coupling+get_species(index, simulation->parameters))+get_species(box->index, simulation->parameters))*sin( (-(simulation->state->particles+box->index)->theta - (simulation->state->particles+index)->theta)*simulation->parameters->interaction_order  );
			*neighbor_n+=1;
		}
		box=box->next_particle;
	}
	return;
}

//calculate the interactions of particle 'index' with the mirror immage of all particles in a given box under reflection with both axis for Langevin models
void L_interaction_with_neighbors_in_box_mirrorxy(gint32 box_index, struct VM_simulation * simulation, guint64 index, gdouble * delta_v, gint32 * neighbor_n)
{
	gdouble distance, distance_squared;
	gdouble xposition_focus=(simulation->state->particles+index)->x;
	gdouble yposition_focus=(simulation->state->particles+index)->y;
	struct VM_particle* box= *(simulation->parameters->box+box_index);
	while (box!=NULL)
	{
		distance_squared=0.0;
		distance=xposition_focus+(simulation->state->particles+box->index)->x;
		while (distance>0.5*simulation->state->length_x)
				distance -= simulation->state->length_x;
		while (distance<-0.5*simulation->state->length_x)
				distance += simulation->state->length_x;
		distance_squared += distance*distance;
		distance=yposition_focus+(simulation->state->particles+box->index)->y;
		while (distance>0.5*simulation->state->length_y)
				distance -= simulation->state->length_y;
		while (distance<-0.5*simulation->state->length_y)
				distance += simulation->state->length_y;
		distance_squared += distance*distance;
		if (distance_squared<((simulation->parameters->interaction_radius)*(simulation->parameters->interaction_radius)))
		{
			if (simulation->parameters->particle_species==1)
				*delta_v+=**simulation->parameters->coupling*sin( (M_PI+(simulation->state->particles+box->index)->theta - (simulation->state->particles+index)->theta)*simulation->parameters->interaction_order  );
			else
				*delta_v+=*(*(simulation->parameters->coupling+get_species(index, simulation->parameters))+get_species(box->index, simulation->parameters))*sin( (M_PI+(simulation->state->particles+box->index)->theta - (simulation->state->particles+index)->theta)*simulation->parameters->interaction_order  );
			*neighbor_n+=1;
		}
		box=box->next_particle;
	}
	return;
}

//calculate interaction of particle 'index' with all its neighbors for Vicsek model
//returns number of neighbors of focus particle (not counting the focus particles itself)
gint32 VM_interaction_with_neighbors(guint64 index, struct VM_simulation * simulation)
{
	gdouble velocity_x=0.0;
	gdouble velocity_y=0.0;
	gint32 neighbor_n=0;
	gdouble xposition_focus=(simulation->state->particles+index)->x;
	gdouble yposition_focus=(simulation->state->particles+index)->y;
	gdouble virtual_position_x=xposition_focus;
	gdouble virtual_position_y=yposition_focus;
	//central box
	VM_interaction_with_neighbors_in_box(VM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index, &velocity_x, &velocity_y, &neighbor_n);
	//bottom box
	virtual_position_y=yposition_focus - simulation->parameters->box_size_y;
	if (virtual_position_y<0.0) //neighbor box is out of y-range -> apply boundary conditions
	{
		if (simulation->parameters->boundary_y==0) //use periodic boundary conditions in y-direction
		{
			virtual_position_y+=simulation->state->length_y;
			VM_interaction_with_neighbors_in_box(VM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index, &velocity_x, &velocity_y, &neighbor_n);
		}
		else //use reflecting boundary conditions in y-direction
			VM_interaction_with_neighbors_in_box_mirrory(VM_get_box_index(virtual_position_x, yposition_focus, simulation->parameters), simulation, index, &velocity_x, &velocity_y, &neighbor_n);
	}
	else
		VM_interaction_with_neighbors_in_box(VM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index, &velocity_x, &velocity_y, &neighbor_n);
	//top box
	virtual_position_y=yposition_focus + simulation->parameters->box_size_y;
	if (virtual_position_y>simulation->state->length_y) //neighbor box is out of y-range -> apply boundary conditions
	{
		if (simulation->parameters->boundary_y==0) //use periodic boundary conditions in y-direction
		{
			virtual_position_y-=simulation->state->length_y;
			VM_interaction_with_neighbors_in_box(VM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index, &velocity_x, &velocity_y, &neighbor_n);
		}
		else //use reflecting boundary conditions in y-direction
			VM_interaction_with_neighbors_in_box_mirrory(VM_get_box_index(virtual_position_x, yposition_focus, simulation->parameters), simulation, index, &velocity_x, &velocity_y, &neighbor_n);
	}
	else
		VM_interaction_with_neighbors_in_box(VM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index, &velocity_x, &velocity_y, &neighbor_n);
	//left box
	virtual_position_y=yposition_focus;
	virtual_position_x=xposition_focus - simulation->parameters->box_size_x;
	if (virtual_position_x<0.0) //neighbor box is out of x-range -> apply boundary conditions
	{
		if (simulation->parameters->boundary_x==0) //use periodic boundary conditions in x-direction
		{
			virtual_position_x+=simulation->state->length_x;
			VM_interaction_with_neighbors_in_box(VM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index, &velocity_x, &velocity_y, &neighbor_n);
		}
		else //use reflecting boundary conditions in x-direction
			VM_interaction_with_neighbors_in_box_mirrorx(VM_get_box_index(xposition_focus, virtual_position_y, simulation->parameters), simulation, index, &velocity_x, &velocity_y, &neighbor_n);
	}
	else
		VM_interaction_with_neighbors_in_box(VM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index, &velocity_x, &velocity_y, &neighbor_n);
	//left bottom box
	virtual_position_y=yposition_focus - simulation->parameters->box_size_y;
	if ((virtual_position_y<0.0) && !(virtual_position_x<0.0))
	{
		if (simulation->parameters->boundary_y==0)
		{
			virtual_position_y+=simulation->state->length_y;
			VM_interaction_with_neighbors_in_box(VM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index, &velocity_x, &velocity_y, &neighbor_n);
		}
		else
			VM_interaction_with_neighbors_in_box_mirrory(VM_get_box_index(virtual_position_x, yposition_focus, simulation->parameters), simulation, index, &velocity_x, &velocity_y, &neighbor_n);
	}
	else if (!(virtual_position_y<0.0) && (virtual_position_x<0.0))
	{
		if (simulation->parameters->boundary_x==0)
		{
			virtual_position_x+=simulation->state->length_x;
			VM_interaction_with_neighbors_in_box(VM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index, &velocity_x, &velocity_y, &neighbor_n);
		}
		else
			VM_interaction_with_neighbors_in_box_mirrorx(VM_get_box_index(xposition_focus, virtual_position_y, simulation->parameters), simulation, index, &velocity_x, &velocity_y, &neighbor_n);
	}
	else if ((virtual_position_y<0.0) && (virtual_position_x<0.0))
	{
		if ((simulation->parameters->boundary_x!=0) && (simulation->parameters->boundary_y!=0))
		{
			VM_interaction_with_neighbors_in_box_mirrorxy(VM_get_box_index(xposition_focus, yposition_focus, simulation->parameters), simulation, index, &velocity_x, &velocity_y, &neighbor_n);
		}
		else if ((simulation->parameters->boundary_x==0) && (simulation->parameters->boundary_y!=0))
		{
			virtual_position_x+=simulation->state->length_x;
			VM_interaction_with_neighbors_in_box_mirrory(VM_get_box_index(virtual_position_x, yposition_focus, simulation->parameters), simulation, index, &velocity_x, &velocity_y, &neighbor_n);
		}
		else if ((simulation->parameters->boundary_x!=0) && (simulation->parameters->boundary_y==0))
		{
			virtual_position_y+=simulation->state->length_y;
			VM_interaction_with_neighbors_in_box_mirrorx(VM_get_box_index(xposition_focus, virtual_position_y, simulation->parameters), simulation, index, &velocity_x, &velocity_y, &neighbor_n);
		}
		else if ((simulation->parameters->boundary_x==0) && (simulation->parameters->boundary_y==0))
		{
			virtual_position_y+=simulation->state->length_y;
			virtual_position_x+=simulation->state->length_x;
			VM_interaction_with_neighbors_in_box(VM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index, &velocity_x, &velocity_y, &neighbor_n);
		}
	}
	else
		VM_interaction_with_neighbors_in_box(VM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index, &velocity_x, &velocity_y, &neighbor_n);
	//left top box
	virtual_position_y=yposition_focus + simulation->parameters->box_size_y;
	if ((virtual_position_y>simulation->state->length_y) && !(virtual_position_x<0.0))
	{
		if (simulation->parameters->boundary_y==0)
		{
			virtual_position_y-=simulation->state->length_y;
			VM_interaction_with_neighbors_in_box(VM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index, &velocity_x, &velocity_y, &neighbor_n);
		}
		else
			VM_interaction_with_neighbors_in_box_mirrory(VM_get_box_index(virtual_position_x, yposition_focus, simulation->parameters), simulation, index, &velocity_x, &velocity_y, &neighbor_n);
	}
	else if (!(virtual_position_y>simulation->state->length_y) && (virtual_position_x<0.0))
	{
		if (simulation->parameters->boundary_x==0)
		{
			virtual_position_x+=simulation->state->length_x;
			VM_interaction_with_neighbors_in_box(VM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index, &velocity_x, &velocity_y, &neighbor_n);
		}
		else
			VM_interaction_with_neighbors_in_box_mirrorx(VM_get_box_index(xposition_focus, virtual_position_y, simulation->parameters), simulation, index, &velocity_x, &velocity_y, &neighbor_n);
	}
	else if ((virtual_position_y>simulation->state->length_y) && (virtual_position_x<0.0))
	{
		if ((simulation->parameters->boundary_x!=0) && (simulation->parameters->boundary_y!=0))
		{
			VM_interaction_with_neighbors_in_box_mirrorxy(VM_get_box_index(xposition_focus, yposition_focus, simulation->parameters), simulation, index, &velocity_x, &velocity_y, &neighbor_n);
		}
		else if ((simulation->parameters->boundary_x==0) && (simulation->parameters->boundary_y!=0))
		{
			virtual_position_x+=simulation->state->length_x;
			VM_interaction_with_neighbors_in_box_mirrory(VM_get_box_index(virtual_position_x, yposition_focus, simulation->parameters), simulation, index, &velocity_x, &velocity_y, &neighbor_n);
		}
		else if ((simulation->parameters->boundary_x!=0) && (simulation->parameters->boundary_y==0))
		{
			virtual_position_y-=simulation->state->length_y;
			VM_interaction_with_neighbors_in_box_mirrorx(VM_get_box_index(xposition_focus, virtual_position_y, simulation->parameters), simulation, index, &velocity_x, &velocity_y, &neighbor_n);
		}
		else if ((simulation->parameters->boundary_x==0) && (simulation->parameters->boundary_y==0))
		{
			virtual_position_y-=simulation->state->length_y;
			virtual_position_x+=simulation->state->length_x;
			VM_interaction_with_neighbors_in_box(VM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index, &velocity_x, &velocity_y, &neighbor_n);
		}
	}
	else
		VM_interaction_with_neighbors_in_box(VM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index, &velocity_x, &velocity_y, &neighbor_n);
	// right box
	virtual_position_y=yposition_focus;
	virtual_position_x=xposition_focus + simulation->parameters->box_size_x;
	if (virtual_position_x>simulation->state->length_x) //neighbor box is out of x-range -> apply boundary conditions
	{
		if (simulation->parameters->boundary_x==0) //use periodic boundary conditions in x-direction
		{
			virtual_position_x-=simulation->state->length_x;
			VM_interaction_with_neighbors_in_box(VM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index, &velocity_x, &velocity_y, &neighbor_n);
		}
		else  //use reflecting boundary conditions in x-direction
			VM_interaction_with_neighbors_in_box_mirrorx(VM_get_box_index(xposition_focus, virtual_position_y, simulation->parameters), simulation, index, &velocity_x, &velocity_y, &neighbor_n);
	}
	else
		VM_interaction_with_neighbors_in_box(VM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index, &velocity_x, &velocity_y, &neighbor_n);
	// right bottom box
	virtual_position_y=yposition_focus - simulation->parameters->box_size_y;
	if ((virtual_position_y<0.0) && !(virtual_position_x>simulation->state->length_x))
	{
		if (simulation->parameters->boundary_y==0)
		{
			virtual_position_y+=simulation->state->length_y;
			VM_interaction_with_neighbors_in_box(VM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index, &velocity_x, &velocity_y, &neighbor_n);
		}
		else
			VM_interaction_with_neighbors_in_box_mirrory(VM_get_box_index(virtual_position_x, yposition_focus, simulation->parameters), simulation, index, &velocity_x, &velocity_y, &neighbor_n);
	}
	else if (!(virtual_position_y<0.0) && (virtual_position_x>simulation->state->length_x))
	{
		if (simulation->parameters->boundary_x==0)
		{
			virtual_position_x-=simulation->state->length_x;
			VM_interaction_with_neighbors_in_box(VM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index, &velocity_x, &velocity_y, &neighbor_n);
		}
		else
			VM_interaction_with_neighbors_in_box_mirrorx(VM_get_box_index(xposition_focus, virtual_position_y, simulation->parameters), simulation, index, &velocity_x, &velocity_y, &neighbor_n);
	}
	else if ((virtual_position_y<0.0) && (virtual_position_x>simulation->state->length_x))
	{
		if ((simulation->parameters->boundary_x!=0) && (simulation->parameters->boundary_y!=0))
		{
			VM_interaction_with_neighbors_in_box_mirrorxy(VM_get_box_index(xposition_focus, yposition_focus, simulation->parameters), simulation, index, &velocity_x, &velocity_y, &neighbor_n);
		}
		else if ((simulation->parameters->boundary_x==0) && (simulation->parameters->boundary_y!=0))
		{
			virtual_position_x-=simulation->state->length_x;
			VM_interaction_with_neighbors_in_box_mirrory(VM_get_box_index(virtual_position_x, yposition_focus, simulation->parameters), simulation, index, &velocity_x, &velocity_y, &neighbor_n);
		}
		else if ((simulation->parameters->boundary_x!=0) && (simulation->parameters->boundary_y==0))
		{
			virtual_position_y+=simulation->state->length_y;
			VM_interaction_with_neighbors_in_box_mirrorx(VM_get_box_index(xposition_focus, virtual_position_y, simulation->parameters), simulation, index, &velocity_x, &velocity_y, &neighbor_n);
		}
		else if ((simulation->parameters->boundary_x==0) && (simulation->parameters->boundary_y==0))
		{
			virtual_position_y+=simulation->state->length_y;
			virtual_position_x-=simulation->state->length_x;
			VM_interaction_with_neighbors_in_box(VM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index, &velocity_x, &velocity_y, &neighbor_n);
		}
	}
	else
		VM_interaction_with_neighbors_in_box(VM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index, &velocity_x, &velocity_y, &neighbor_n);
	// right top box
	virtual_position_y=yposition_focus + simulation->parameters->box_size_y;
	if ((virtual_position_y>simulation->state->length_y) && !(virtual_position_x>simulation->state->length_x))
	{
		if (simulation->parameters->boundary_y==0)
		{
			virtual_position_y-=simulation->state->length_y;
			VM_interaction_with_neighbors_in_box(VM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index, &velocity_x, &velocity_y, &neighbor_n);
		}
		else
			VM_interaction_with_neighbors_in_box_mirrory(VM_get_box_index(virtual_position_x, yposition_focus, simulation->parameters), simulation, index, &velocity_x, &velocity_y, &neighbor_n);
	}
	else if (!(virtual_position_y>simulation->state->length_y) && (virtual_position_x>simulation->state->length_x))
	{
		if (simulation->parameters->boundary_x==0)
		{
			virtual_position_x-=simulation->state->length_x;
			VM_interaction_with_neighbors_in_box(VM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index, &velocity_x, &velocity_y, &neighbor_n);
		}
		else
			VM_interaction_with_neighbors_in_box_mirrorx(VM_get_box_index(xposition_focus, virtual_position_y, simulation->parameters), simulation, index, &velocity_x, &velocity_y, &neighbor_n);
	}
	else if ((virtual_position_y>simulation->state->length_y) && (virtual_position_x>simulation->state->length_x))
	{
		if ((simulation->parameters->boundary_x!=0) && (simulation->parameters->boundary_y!=0))
		{
			VM_interaction_with_neighbors_in_box_mirrorxy(VM_get_box_index(xposition_focus, yposition_focus, simulation->parameters), simulation, index, &velocity_x, &velocity_y, &neighbor_n);
		}
		else if ((simulation->parameters->boundary_x==0) && (simulation->parameters->boundary_y!=0))
		{
			virtual_position_x-=simulation->state->length_x;
			VM_interaction_with_neighbors_in_box_mirrory(VM_get_box_index(virtual_position_x, yposition_focus, simulation->parameters), simulation, index, &velocity_x, &velocity_y, &neighbor_n);
		}
		else if ((simulation->parameters->boundary_x!=0) && (simulation->parameters->boundary_y==0))
		{
			virtual_position_y-=simulation->state->length_y;
			VM_interaction_with_neighbors_in_box_mirrorx(VM_get_box_index(xposition_focus, virtual_position_y, simulation->parameters), simulation, index, &velocity_x, &velocity_y, &neighbor_n);
		}
		else if ((simulation->parameters->boundary_x==0) && (simulation->parameters->boundary_y==0))
		{
			virtual_position_y-=simulation->state->length_y;
			virtual_position_x-=simulation->state->length_x;
			VM_interaction_with_neighbors_in_box(VM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index, &velocity_x, &velocity_y, &neighbor_n);
		}
	}
	else
		VM_interaction_with_neighbors_in_box(VM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index, &velocity_x, &velocity_y, &neighbor_n);
	//update angle
	if (simulation->parameters->particle_species==1)
		(simulation->state->particles+index)->theta_new=atan2(velocity_y, velocity_x)+(g_rand_double(simulation->state->prng)*2*M_PI -M_PI)**simulation->parameters->noise_strength+*simulation->parameters->omega;
	else
		(simulation->state->particles+index)->theta_new=atan2(velocity_y, velocity_x)+(g_rand_double(simulation->state->prng)*2*M_PI -M_PI)**(simulation->parameters->noise_strength+get_species(index, simulation->parameters))+*(simulation->parameters->omega+get_species(index, simulation->parameters));
	//set -pi<angle<pi
	while ((simulation->state->particles+index)->theta_new<-M_PI)
		(simulation->state->particles+index)->theta_new+=2.0*M_PI;
	while ((simulation->state->particles+index)->theta_new>M_PI)
		(simulation->state->particles+index)->theta_new-=2.0*M_PI;
	return neighbor_n;
}

//calculate interaction of particle 'index' with all its neighbors for nematic Vicsek model
//returns number of neighbors of focus particle (not counting the focus particles itself)
gint32 NVM_interaction_with_neighbors(guint64 index, struct VM_simulation * simulation)
{
	gdouble velocity_x=0.0;
	gdouble velocity_y=0.0;
	gint32 neighbor_n=0;
	gdouble xposition_focus=(simulation->state->particles+index)->x;
	gdouble yposition_focus=(simulation->state->particles+index)->y;
	gdouble virtual_position_x=xposition_focus;
	gdouble virtual_position_y=yposition_focus;
	//central box
	NVM_interaction_with_neighbors_in_box(VM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index, &velocity_x, &velocity_y, &neighbor_n);
	//bottom box
	virtual_position_y=yposition_focus - simulation->parameters->box_size_y;
	if (virtual_position_y<0.0) //neighbor box is out of y-range -> apply boundary conditions
	{
		if (simulation->parameters->boundary_y==0) //use periodic boundary conditions in y-direction
		{
			virtual_position_y+=simulation->state->length_y;
			NVM_interaction_with_neighbors_in_box(VM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index, &velocity_x, &velocity_y, &neighbor_n);
		}
		else //use reflecting boundary conditions in y-direction
			NVM_interaction_with_neighbors_in_box_mirrory(VM_get_box_index(virtual_position_x, yposition_focus, simulation->parameters), simulation, index, &velocity_x, &velocity_y, &neighbor_n);
	}
	else
		NVM_interaction_with_neighbors_in_box(VM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index, &velocity_x, &velocity_y, &neighbor_n);
	//top box
	virtual_position_y=yposition_focus + simulation->parameters->box_size_y;
	if (virtual_position_y>simulation->state->length_y) //neighbor box is out of y-range -> apply boundary conditions
	{
		if (simulation->parameters->boundary_y==0) //use periodic boundary conditions in y-direction
		{
			virtual_position_y-=simulation->state->length_y;
			NVM_interaction_with_neighbors_in_box(VM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index, &velocity_x, &velocity_y, &neighbor_n);
		}
		else //use reflecting boundary conditions in y-direction
			NVM_interaction_with_neighbors_in_box_mirrory(VM_get_box_index(virtual_position_x, yposition_focus, simulation->parameters), simulation, index, &velocity_x, &velocity_y, &neighbor_n);
	}
	else
		NVM_interaction_with_neighbors_in_box(VM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index, &velocity_x, &velocity_y, &neighbor_n);
	//left box
	virtual_position_y=yposition_focus;
	virtual_position_x=xposition_focus - simulation->parameters->box_size_x;
	if (virtual_position_x<0.0) //neighbor box is out of x-range -> apply boundary conditions
	{
		if (simulation->parameters->boundary_x==0) //use periodic boundary conditions in x-direction
		{
			virtual_position_x+=simulation->state->length_x;
			NVM_interaction_with_neighbors_in_box(VM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index, &velocity_x, &velocity_y, &neighbor_n);
		}
		else //use reflecting boundary conditions in x-direction
			NVM_interaction_with_neighbors_in_box_mirrorx(VM_get_box_index(xposition_focus, virtual_position_y, simulation->parameters), simulation, index, &velocity_x, &velocity_y, &neighbor_n);
	}
	else
		NVM_interaction_with_neighbors_in_box(VM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index, &velocity_x, &velocity_y, &neighbor_n);
	//left bottom box
	virtual_position_y=yposition_focus - simulation->parameters->box_size_y;
	if ((virtual_position_y<0.0) && !(virtual_position_x<0.0))
	{
		if (simulation->parameters->boundary_y==0)
		{
			virtual_position_y+=simulation->state->length_y;
			NVM_interaction_with_neighbors_in_box(VM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index, &velocity_x, &velocity_y, &neighbor_n);
		}
		else
			NVM_interaction_with_neighbors_in_box_mirrory(VM_get_box_index(virtual_position_x, yposition_focus, simulation->parameters), simulation, index, &velocity_x, &velocity_y, &neighbor_n);
	}
	else if (!(virtual_position_y<0.0) && (virtual_position_x<0.0))
	{
		if (simulation->parameters->boundary_x==0)
		{
			virtual_position_x+=simulation->state->length_x;
			NVM_interaction_with_neighbors_in_box(VM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index, &velocity_x, &velocity_y, &neighbor_n);
		}
		else
			NVM_interaction_with_neighbors_in_box_mirrorx(VM_get_box_index(xposition_focus, virtual_position_y, simulation->parameters), simulation, index, &velocity_x, &velocity_y, &neighbor_n);
	}
	else if ((virtual_position_y<0.0) && (virtual_position_x<0.0))
	{
		if ((simulation->parameters->boundary_x!=0) && (simulation->parameters->boundary_y!=0))
		{
			NVM_interaction_with_neighbors_in_box_mirrorxy(VM_get_box_index(xposition_focus, yposition_focus, simulation->parameters), simulation, index, &velocity_x, &velocity_y, &neighbor_n);
		}
		else if ((simulation->parameters->boundary_x==0) && (simulation->parameters->boundary_y!=0))
		{
			virtual_position_x+=simulation->state->length_x;
			NVM_interaction_with_neighbors_in_box_mirrory(VM_get_box_index(virtual_position_x, yposition_focus, simulation->parameters), simulation, index, &velocity_x, &velocity_y, &neighbor_n);
		}
		else if ((simulation->parameters->boundary_x!=0) && (simulation->parameters->boundary_y==0))
		{
			virtual_position_y+=simulation->state->length_y;
			NVM_interaction_with_neighbors_in_box_mirrorx(VM_get_box_index(xposition_focus, virtual_position_y, simulation->parameters), simulation, index, &velocity_x, &velocity_y, &neighbor_n);
		}
		else if ((simulation->parameters->boundary_x==0) && (simulation->parameters->boundary_y==0))
		{
			virtual_position_y+=simulation->state->length_y;
			virtual_position_x+=simulation->state->length_x;
			NVM_interaction_with_neighbors_in_box(VM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index, &velocity_x, &velocity_y, &neighbor_n);
		}
	}
	else
		NVM_interaction_with_neighbors_in_box(VM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index, &velocity_x, &velocity_y, &neighbor_n);
	//left top box
	virtual_position_y=yposition_focus + simulation->parameters->box_size_y;
	if ((virtual_position_y>simulation->state->length_y) && !(virtual_position_x<0.0))
	{
		if (simulation->parameters->boundary_y==0)
		{
			virtual_position_y-=simulation->state->length_y;
			NVM_interaction_with_neighbors_in_box(VM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index, &velocity_x, &velocity_y, &neighbor_n);
		}
		else
			NVM_interaction_with_neighbors_in_box_mirrory(VM_get_box_index(virtual_position_x, yposition_focus, simulation->parameters), simulation, index, &velocity_x, &velocity_y, &neighbor_n);
	}
	else if (!(virtual_position_y>simulation->state->length_y) && (virtual_position_x<0.0))
	{
		if (simulation->parameters->boundary_x==0)
		{
			virtual_position_x+=simulation->state->length_x;
			NVM_interaction_with_neighbors_in_box(VM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index, &velocity_x, &velocity_y, &neighbor_n);
		}
		else
			NVM_interaction_with_neighbors_in_box_mirrorx(VM_get_box_index(xposition_focus, virtual_position_y, simulation->parameters), simulation, index, &velocity_x, &velocity_y, &neighbor_n);
	}
	else if ((virtual_position_y>simulation->state->length_y) && (virtual_position_x<0.0))
	{
		if ((simulation->parameters->boundary_x!=0) && (simulation->parameters->boundary_y!=0))
		{
			NVM_interaction_with_neighbors_in_box_mirrorxy(VM_get_box_index(xposition_focus, yposition_focus, simulation->parameters), simulation, index, &velocity_x, &velocity_y, &neighbor_n);
		}
		else if ((simulation->parameters->boundary_x==0) && (simulation->parameters->boundary_y!=0))
		{
			virtual_position_x+=simulation->state->length_x;
			NVM_interaction_with_neighbors_in_box_mirrory(VM_get_box_index(virtual_position_x, yposition_focus, simulation->parameters), simulation, index, &velocity_x, &velocity_y, &neighbor_n);
		}
		else if ((simulation->parameters->boundary_x!=0) && (simulation->parameters->boundary_y==0))
		{
			virtual_position_y-=simulation->state->length_y;
			NVM_interaction_with_neighbors_in_box_mirrorx(VM_get_box_index(xposition_focus, virtual_position_y, simulation->parameters), simulation, index, &velocity_x, &velocity_y, &neighbor_n);
		}
		else if ((simulation->parameters->boundary_x==0) && (simulation->parameters->boundary_y==0))
		{
			virtual_position_y-=simulation->state->length_y;
			virtual_position_x+=simulation->state->length_x;
			NVM_interaction_with_neighbors_in_box(VM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index, &velocity_x, &velocity_y, &neighbor_n);
		}
	}
	else
		NVM_interaction_with_neighbors_in_box(VM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index, &velocity_x, &velocity_y, &neighbor_n);
	// right box
	virtual_position_y=yposition_focus;
	virtual_position_x=xposition_focus + simulation->parameters->box_size_x;
	if (virtual_position_x>simulation->state->length_x) //neighbor box is out of x-range -> apply boundary conditions
	{
		if (simulation->parameters->boundary_x==0) //use periodic boundary conditions in x-direction
		{
			virtual_position_x-=simulation->state->length_x;
			NVM_interaction_with_neighbors_in_box(VM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index, &velocity_x, &velocity_y, &neighbor_n);
		}
		else  //use reflecting boundary conditions in x-direction
			NVM_interaction_with_neighbors_in_box_mirrorx(VM_get_box_index(xposition_focus, virtual_position_y, simulation->parameters), simulation, index, &velocity_x, &velocity_y, &neighbor_n);
	}
	else
		NVM_interaction_with_neighbors_in_box(VM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index, &velocity_x, &velocity_y, &neighbor_n);
	// right bottom box
	virtual_position_y=yposition_focus - simulation->parameters->box_size_y;
	if ((virtual_position_y<0.0) && !(virtual_position_x>simulation->state->length_x))
	{
		if (simulation->parameters->boundary_y==0)
		{
			virtual_position_y+=simulation->state->length_y;
			NVM_interaction_with_neighbors_in_box(VM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index, &velocity_x, &velocity_y, &neighbor_n);
		}
		else
			NVM_interaction_with_neighbors_in_box_mirrory(VM_get_box_index(virtual_position_x, yposition_focus, simulation->parameters), simulation, index, &velocity_x, &velocity_y, &neighbor_n);
	}
	else if (!(virtual_position_y<0.0) && (virtual_position_x>simulation->state->length_x))
	{
		if (simulation->parameters->boundary_x==0)
		{
			virtual_position_x-=simulation->state->length_x;
			NVM_interaction_with_neighbors_in_box(VM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index, &velocity_x, &velocity_y, &neighbor_n);
		}
		else
			NVM_interaction_with_neighbors_in_box_mirrorx(VM_get_box_index(xposition_focus, virtual_position_y, simulation->parameters), simulation, index, &velocity_x, &velocity_y, &neighbor_n);
	}
	else if ((virtual_position_y<0.0) && (virtual_position_x>simulation->state->length_x))
	{
		if ((simulation->parameters->boundary_x!=0) && (simulation->parameters->boundary_y!=0))
		{
			NVM_interaction_with_neighbors_in_box_mirrorxy(VM_get_box_index(xposition_focus, yposition_focus, simulation->parameters), simulation, index, &velocity_x, &velocity_y, &neighbor_n);
		}
		else if ((simulation->parameters->boundary_x==0) && (simulation->parameters->boundary_y!=0))
		{
			virtual_position_x-=simulation->state->length_x;
			NVM_interaction_with_neighbors_in_box_mirrory(VM_get_box_index(virtual_position_x, yposition_focus, simulation->parameters), simulation, index, &velocity_x, &velocity_y, &neighbor_n);
		}
		else if ((simulation->parameters->boundary_x!=0) && (simulation->parameters->boundary_y==0))
		{
			virtual_position_y+=simulation->state->length_y;
			NVM_interaction_with_neighbors_in_box_mirrorx(VM_get_box_index(xposition_focus, virtual_position_y, simulation->parameters), simulation, index, &velocity_x, &velocity_y, &neighbor_n);
		}
		else if ((simulation->parameters->boundary_x==0) && (simulation->parameters->boundary_y==0))
		{
			virtual_position_y+=simulation->state->length_y;
			virtual_position_x-=simulation->state->length_x;
			NVM_interaction_with_neighbors_in_box(VM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index, &velocity_x, &velocity_y, &neighbor_n);
		}
	}
	else
		NVM_interaction_with_neighbors_in_box(VM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index, &velocity_x, &velocity_y, &neighbor_n);
	// right top box
	virtual_position_y=yposition_focus + simulation->parameters->box_size_y;
	if ((virtual_position_y>simulation->state->length_y) && !(virtual_position_x>simulation->state->length_x))
	{
		if (simulation->parameters->boundary_y==0)
		{
			virtual_position_y-=simulation->state->length_y;
			NVM_interaction_with_neighbors_in_box(VM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index, &velocity_x, &velocity_y, &neighbor_n);
		}
		else
			NVM_interaction_with_neighbors_in_box_mirrory(VM_get_box_index(virtual_position_x, yposition_focus, simulation->parameters), simulation, index, &velocity_x, &velocity_y, &neighbor_n);
	}
	else if (!(virtual_position_y>simulation->state->length_y) && (virtual_position_x>simulation->state->length_x))
	{
		if (simulation->parameters->boundary_x==0)
		{
			virtual_position_x-=simulation->state->length_x;
			NVM_interaction_with_neighbors_in_box(VM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index, &velocity_x, &velocity_y, &neighbor_n);
		}
		else
			NVM_interaction_with_neighbors_in_box_mirrorx(VM_get_box_index(xposition_focus, virtual_position_y, simulation->parameters), simulation, index, &velocity_x, &velocity_y, &neighbor_n);
	}
	else if ((virtual_position_y>simulation->state->length_y) && (virtual_position_x>simulation->state->length_x))
	{
		if ((simulation->parameters->boundary_x!=0) && (simulation->parameters->boundary_y!=0))
		{
			NVM_interaction_with_neighbors_in_box_mirrorxy(VM_get_box_index(xposition_focus, yposition_focus, simulation->parameters), simulation, index, &velocity_x, &velocity_y, &neighbor_n);
		}
		else if ((simulation->parameters->boundary_x==0) && (simulation->parameters->boundary_y!=0))
		{
			virtual_position_x-=simulation->state->length_x;
			NVM_interaction_with_neighbors_in_box_mirrory(VM_get_box_index(virtual_position_x, yposition_focus, simulation->parameters), simulation, index, &velocity_x, &velocity_y, &neighbor_n);
		}
		else if ((simulation->parameters->boundary_x!=0) && (simulation->parameters->boundary_y==0))
		{
			virtual_position_y-=simulation->state->length_y;
			NVM_interaction_with_neighbors_in_box_mirrorx(VM_get_box_index(xposition_focus, virtual_position_y, simulation->parameters), simulation, index, &velocity_x, &velocity_y, &neighbor_n);
		}
		else if ((simulation->parameters->boundary_x==0) && (simulation->parameters->boundary_y==0))
		{
			virtual_position_y-=simulation->state->length_y;
			virtual_position_x-=simulation->state->length_x;
			NVM_interaction_with_neighbors_in_box(VM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index, &velocity_x, &velocity_y, &neighbor_n);
		}
	}
	else
		NVM_interaction_with_neighbors_in_box(VM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index, &velocity_x, &velocity_y, &neighbor_n);
	//update angle
	if (simulation->parameters->particle_species==1)
		(simulation->state->particles+index)->theta_new=atan2(velocity_y, velocity_x)+(g_rand_double(simulation->state->prng)*2*M_PI -M_PI)**simulation->parameters->noise_strength+*simulation->parameters->omega;
	else
		(simulation->state->particles+index)->theta_new=atan2(velocity_y, velocity_x)+(g_rand_double(simulation->state->prng)*2*M_PI -M_PI)**(simulation->parameters->noise_strength+get_species(index, simulation->parameters))+*(simulation->parameters->omega+get_species(index, simulation->parameters));
	//set -pi<angle<pi
	while ((simulation->state->particles+index)->theta_new<-M_PI)
		(simulation->state->particles+index)->theta_new+=2.0*M_PI;
	while ((simulation->state->particles+index)->theta_new>M_PI)
		(simulation->state->particles+index)->theta_new-=2.0*M_PI;
	return neighbor_n;
}

//calculate interaction of particle 'index' with all its neighbors for additive Langevin model
//returns number of neighbors of focus particle (not counting the focus particles itself)
gint32 additiveL_interaction_with_neighbors(guint64 index, struct VM_simulation * simulation)
{
	gdouble delta_v=0.0;
	gint32 neighbor_n=0;
	gdouble xposition_focus=(simulation->state->particles+index)->x;
	gdouble yposition_focus=(simulation->state->particles+index)->y;
	gdouble virtual_position_x=xposition_focus;
	gdouble virtual_position_y=yposition_focus;
	//central box
	L_interaction_with_neighbors_in_box(VM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index, &delta_v, &neighbor_n);
	//bottom box
	virtual_position_y=yposition_focus - simulation->parameters->box_size_y;
	if (virtual_position_y<0.0) //neighbor box is out of y-range -> apply boundary conditions
	{
		if (simulation->parameters->boundary_y==0) //use periodic boundary conditions in y-direction
		{
			virtual_position_y+=simulation->state->length_y;
			L_interaction_with_neighbors_in_box(VM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index, &delta_v, &neighbor_n);
		}
		else //use reflecting boundary conditions in y-direction
			L_interaction_with_neighbors_in_box_mirrory(VM_get_box_index(virtual_position_x, yposition_focus, simulation->parameters), simulation, index, &delta_v, &neighbor_n);
	}
	else
		L_interaction_with_neighbors_in_box(VM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index, &delta_v, &neighbor_n);
	//top box
	virtual_position_y=yposition_focus + simulation->parameters->box_size_y;
	if (virtual_position_y>simulation->state->length_y) //neighbor box is out of y-range -> apply boundary conditions
	{
		if (simulation->parameters->boundary_y==0) //use periodic boundary conditions in y-direction
		{
			virtual_position_y-=simulation->state->length_y;
			L_interaction_with_neighbors_in_box(VM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index, &delta_v, &neighbor_n);
		}
		else //use reflecting boundary conditions in y-direction
			L_interaction_with_neighbors_in_box_mirrory(VM_get_box_index(virtual_position_x, yposition_focus, simulation->parameters), simulation, index, &delta_v, &neighbor_n);
	}
	else
		L_interaction_with_neighbors_in_box(VM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index, &delta_v, &neighbor_n);
	//left box
	virtual_position_y=yposition_focus;
	virtual_position_x=xposition_focus - simulation->parameters->box_size_x;
	if (virtual_position_x<0.0) //neighbor box is out of x-range -> apply boundary conditions
	{
		if (simulation->parameters->boundary_x==0) //use periodic boundary conditions in x-direction
		{
			virtual_position_x+=simulation->state->length_x;
			L_interaction_with_neighbors_in_box(VM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index, &delta_v, &neighbor_n);
		}
		else //use reflecting boundary conditions in x-direction
			L_interaction_with_neighbors_in_box_mirrorx(VM_get_box_index(xposition_focus, virtual_position_y, simulation->parameters), simulation, index, &delta_v, &neighbor_n);
	}
	else
		L_interaction_with_neighbors_in_box(VM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index, &delta_v, &neighbor_n);
	//left bottom box
	virtual_position_y=yposition_focus - simulation->parameters->box_size_y;
	if ((virtual_position_y<0.0) && !(virtual_position_x<0.0))
	{
		if (simulation->parameters->boundary_y==0)
		{
			virtual_position_y+=simulation->state->length_y;
			L_interaction_with_neighbors_in_box(VM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index, &delta_v, &neighbor_n);
		}
		else
			L_interaction_with_neighbors_in_box_mirrory(VM_get_box_index(virtual_position_x, yposition_focus, simulation->parameters), simulation, index, &delta_v, &neighbor_n);
	}
	else if (!(virtual_position_y<0.0) && (virtual_position_x<0.0))
	{
		if (simulation->parameters->boundary_x==0)
		{
			virtual_position_x+=simulation->state->length_x;
			L_interaction_with_neighbors_in_box(VM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index, &delta_v, &neighbor_n);
		}
		else
			L_interaction_with_neighbors_in_box_mirrorx(VM_get_box_index(xposition_focus, virtual_position_y, simulation->parameters), simulation, index, &delta_v, &neighbor_n);
	}
	else if ((virtual_position_y<0.0) && (virtual_position_x<0.0))
	{
		if ((simulation->parameters->boundary_x!=0) && (simulation->parameters->boundary_y!=0))
		{
			L_interaction_with_neighbors_in_box_mirrorxy(VM_get_box_index(xposition_focus, yposition_focus, simulation->parameters), simulation, index, &delta_v, &neighbor_n);
		}
		else if ((simulation->parameters->boundary_x==0) && (simulation->parameters->boundary_y!=0))
		{
			virtual_position_x+=simulation->state->length_x;
			L_interaction_with_neighbors_in_box_mirrory(VM_get_box_index(virtual_position_x, yposition_focus, simulation->parameters), simulation, index, &delta_v, &neighbor_n);
		}
		else if ((simulation->parameters->boundary_x!=0) && (simulation->parameters->boundary_y==0))
		{
			virtual_position_y+=simulation->state->length_y;
			L_interaction_with_neighbors_in_box_mirrorx(VM_get_box_index(xposition_focus, virtual_position_y, simulation->parameters), simulation, index, &delta_v, &neighbor_n);
		}
		else if ((simulation->parameters->boundary_x==0) && (simulation->parameters->boundary_y==0))
		{
			virtual_position_y+=simulation->state->length_y;
			virtual_position_x+=simulation->state->length_x;
			L_interaction_with_neighbors_in_box(VM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index, &delta_v, &neighbor_n);
		}
	}
	else
		L_interaction_with_neighbors_in_box(VM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index, &delta_v, &neighbor_n);
	//left top box
	virtual_position_y=yposition_focus + simulation->parameters->box_size_y;
	if ((virtual_position_y>simulation->state->length_y) && !(virtual_position_x<0.0))
	{
		if (simulation->parameters->boundary_y==0)
		{
			virtual_position_y-=simulation->state->length_y;
			L_interaction_with_neighbors_in_box(VM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index, &delta_v, &neighbor_n);
		}
		else
			L_interaction_with_neighbors_in_box_mirrory(VM_get_box_index(virtual_position_x, yposition_focus, simulation->parameters), simulation, index, &delta_v, &neighbor_n);
	}
	else if (!(virtual_position_y>simulation->state->length_y) && (virtual_position_x<0.0))
	{
		if (simulation->parameters->boundary_x==0)
		{
			virtual_position_x+=simulation->state->length_x;
			L_interaction_with_neighbors_in_box(VM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index, &delta_v, &neighbor_n);
		}
		else
			L_interaction_with_neighbors_in_box_mirrorx(VM_get_box_index(xposition_focus, virtual_position_y, simulation->parameters), simulation, index, &delta_v, &neighbor_n);
	}
	else if ((virtual_position_y>simulation->state->length_y) && (virtual_position_x<0.0))
	{
		if ((simulation->parameters->boundary_x!=0) && (simulation->parameters->boundary_y!=0))
		{
			L_interaction_with_neighbors_in_box_mirrorxy(VM_get_box_index(xposition_focus, yposition_focus, simulation->parameters), simulation, index, &delta_v, &neighbor_n);
		}
		else if ((simulation->parameters->boundary_x==0) && (simulation->parameters->boundary_y!=0))
		{
			virtual_position_x+=simulation->state->length_x;
			L_interaction_with_neighbors_in_box_mirrory(VM_get_box_index(virtual_position_x, yposition_focus, simulation->parameters), simulation, index, &delta_v, &neighbor_n);
		}
		else if ((simulation->parameters->boundary_x!=0) && (simulation->parameters->boundary_y==0))
		{
			virtual_position_y-=simulation->state->length_y;
			L_interaction_with_neighbors_in_box_mirrorx(VM_get_box_index(xposition_focus, virtual_position_y, simulation->parameters), simulation, index, &delta_v, &neighbor_n);
		}
		else if ((simulation->parameters->boundary_x==0) && (simulation->parameters->boundary_y==0))
		{
			virtual_position_y-=simulation->state->length_y;
			virtual_position_x+=simulation->state->length_x;
			L_interaction_with_neighbors_in_box(VM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index, &delta_v, &neighbor_n);
		}
	}
	else
		L_interaction_with_neighbors_in_box(VM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index, &delta_v, &neighbor_n);
	// right box
	virtual_position_y=yposition_focus;
	virtual_position_x=xposition_focus + simulation->parameters->box_size_x;
	if (virtual_position_x>simulation->state->length_x) //neighbor box is out of x-range -> apply boundary conditions
	{
		if (simulation->parameters->boundary_x==0) //use periodic boundary conditions in x-direction
		{
			virtual_position_x-=simulation->state->length_x;
			L_interaction_with_neighbors_in_box(VM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index, &delta_v, &neighbor_n);
		}
		else  //use reflecting boundary conditions in x-direction
			L_interaction_with_neighbors_in_box_mirrorx(VM_get_box_index(xposition_focus, virtual_position_y, simulation->parameters), simulation, index, &delta_v, &neighbor_n);
	}
	else
		L_interaction_with_neighbors_in_box(VM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index, &delta_v, &neighbor_n);
	// right bottom box
	virtual_position_y=yposition_focus - simulation->parameters->box_size_y;
	if ((virtual_position_y<0.0) && !(virtual_position_x>simulation->state->length_x))
	{
		if (simulation->parameters->boundary_y==0)
		{
			virtual_position_y+=simulation->state->length_y;
			L_interaction_with_neighbors_in_box(VM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index, &delta_v, &neighbor_n);
		}
		else
			L_interaction_with_neighbors_in_box_mirrory(VM_get_box_index(virtual_position_x, yposition_focus, simulation->parameters), simulation, index, &delta_v, &neighbor_n);
	}
	else if (!(virtual_position_y<0.0) && (virtual_position_x>simulation->state->length_x))
	{
		if (simulation->parameters->boundary_x==0)
		{
			virtual_position_x-=simulation->state->length_x;
			L_interaction_with_neighbors_in_box(VM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index, &delta_v, &neighbor_n);
		}
		else
			L_interaction_with_neighbors_in_box_mirrorx(VM_get_box_index(xposition_focus, virtual_position_y, simulation->parameters), simulation, index, &delta_v, &neighbor_n);
	}
	else if ((virtual_position_y<0.0) && (virtual_position_x>simulation->state->length_x))
	{
		if ((simulation->parameters->boundary_x!=0) && (simulation->parameters->boundary_y!=0))
		{
			L_interaction_with_neighbors_in_box_mirrorxy(VM_get_box_index(xposition_focus, yposition_focus, simulation->parameters), simulation, index, &delta_v, &neighbor_n);
		}
		else if ((simulation->parameters->boundary_x==0) && (simulation->parameters->boundary_y!=0))
		{
			virtual_position_x-=simulation->state->length_x;
			L_interaction_with_neighbors_in_box_mirrory(VM_get_box_index(virtual_position_x, yposition_focus, simulation->parameters), simulation, index, &delta_v, &neighbor_n);
		}
		else if ((simulation->parameters->boundary_x!=0) && (simulation->parameters->boundary_y==0))
		{
			virtual_position_y+=simulation->state->length_y;
			L_interaction_with_neighbors_in_box_mirrorx(VM_get_box_index(xposition_focus, virtual_position_y, simulation->parameters), simulation, index, &delta_v, &neighbor_n);
		}
		else if ((simulation->parameters->boundary_x==0) && (simulation->parameters->boundary_y==0))
		{
			virtual_position_y+=simulation->state->length_y;
			virtual_position_x-=simulation->state->length_x;
			L_interaction_with_neighbors_in_box(VM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index, &delta_v, &neighbor_n);
		}
	}
	else
		L_interaction_with_neighbors_in_box(VM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index, &delta_v, &neighbor_n);
	// right top box
	virtual_position_y=yposition_focus + simulation->parameters->box_size_y;
	if ((virtual_position_y>simulation->state->length_y) && !(virtual_position_x>simulation->state->length_x))
	{
		if (simulation->parameters->boundary_y==0)
		{
			virtual_position_y-=simulation->state->length_y;
			L_interaction_with_neighbors_in_box(VM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index, &delta_v, &neighbor_n);
		}
		else
			L_interaction_with_neighbors_in_box_mirrory(VM_get_box_index(virtual_position_x, yposition_focus, simulation->parameters), simulation, index, &delta_v, &neighbor_n);
	}
	else if (!(virtual_position_y>simulation->state->length_y) && (virtual_position_x>simulation->state->length_x))
	{
		if (simulation->parameters->boundary_x==0)
		{
			virtual_position_x-=simulation->state->length_x;
			L_interaction_with_neighbors_in_box(VM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index, &delta_v, &neighbor_n);
		}
		else
			L_interaction_with_neighbors_in_box_mirrorx(VM_get_box_index(xposition_focus, virtual_position_y, simulation->parameters), simulation, index, &delta_v, &neighbor_n);
	}
	else if ((virtual_position_y>simulation->state->length_y) && (virtual_position_x>simulation->state->length_x))
	{
		if ((simulation->parameters->boundary_x!=0) && (simulation->parameters->boundary_y!=0))
		{
			L_interaction_with_neighbors_in_box_mirrorxy(VM_get_box_index(xposition_focus, yposition_focus, simulation->parameters), simulation, index, &delta_v, &neighbor_n);
		}
		else if ((simulation->parameters->boundary_x==0) && (simulation->parameters->boundary_y!=0))
		{
			virtual_position_x-=simulation->state->length_x;
			L_interaction_with_neighbors_in_box_mirrory(VM_get_box_index(virtual_position_x, yposition_focus, simulation->parameters), simulation, index, &delta_v, &neighbor_n);
		}
		else if ((simulation->parameters->boundary_x!=0) && (simulation->parameters->boundary_y==0))
		{
			virtual_position_y-=simulation->state->length_y;
			L_interaction_with_neighbors_in_box_mirrorx(VM_get_box_index(xposition_focus, virtual_position_y, simulation->parameters), simulation, index, &delta_v, &neighbor_n);
		}
		else if ((simulation->parameters->boundary_x==0) && (simulation->parameters->boundary_y==0))
		{
			virtual_position_y-=simulation->state->length_y;
			virtual_position_x-=simulation->state->length_x;
			L_interaction_with_neighbors_in_box(VM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index, &delta_v, &neighbor_n);
		}
	}
	else
		L_interaction_with_neighbors_in_box(VM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index, &delta_v, &neighbor_n);
	//update the angle
	if (simulation->parameters->particle_species==1)
		(simulation->state->particles+index)->theta_new=(simulation->state->particles+index)->theta + simulation->parameters->delta_t*(*simulation->parameters->omega+delta_v) +simulation->parameters->sqrt_delta_t**simulation->parameters->noise_strength*gauss(simulation);
	else
		(simulation->state->particles+index)->theta_new=(simulation->state->particles+index)->theta + simulation->parameters->delta_t*(*(simulation->parameters->omega+get_species(index, simulation->parameters))+delta_v) +simulation->parameters->sqrt_delta_t**(simulation->parameters->noise_strength+get_species(index, simulation->parameters))*gauss(simulation);
	while ((simulation->state->particles+index)->theta_new<-M_PI)
		(simulation->state->particles+index)->theta_new+=2.0*M_PI;
	while ((simulation->state->particles+index)->theta_new>M_PI)
		(simulation->state->particles+index)->theta_new-=2.0*M_PI;
	return neighbor_n;
}

//calculate interaction of particle 'index' with all its neighbors for nonadditive Langevin model
//returns number of neighbors of focus particle (not counting the focus particles itself)
gint32 nonadditiveL_interaction_with_neighbors(guint64 index, struct VM_simulation * simulation)
{
	gdouble delta_v=0.0;
	gint32 neighbor_n=0;
	gdouble xposition_focus=(simulation->state->particles+index)->x;
	gdouble yposition_focus=(simulation->state->particles+index)->y;
	gdouble virtual_position_x=xposition_focus;
	gdouble virtual_position_y=yposition_focus;
	//central box
	L_interaction_with_neighbors_in_box(VM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index, &delta_v, &neighbor_n);
	//bottom box
	virtual_position_y=yposition_focus - simulation->parameters->box_size_y;
	if (virtual_position_y<0.0) //neighbor box is out of y-range -> apply boundary conditions
	{
		if (simulation->parameters->boundary_y==0) //use periodic boundary conditions in y-direction
		{
			virtual_position_y+=simulation->state->length_y;
			L_interaction_with_neighbors_in_box(VM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index, &delta_v, &neighbor_n);
		}
		else //use reflecting boundary conditions in y-direction
			L_interaction_with_neighbors_in_box_mirrory(VM_get_box_index(virtual_position_x, yposition_focus, simulation->parameters), simulation, index, &delta_v, &neighbor_n);
	}
	else
		L_interaction_with_neighbors_in_box(VM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index, &delta_v, &neighbor_n);
	//top box
	virtual_position_y=yposition_focus + simulation->parameters->box_size_y;
	if (virtual_position_y>simulation->state->length_y) //neighbor box is out of y-range -> apply boundary conditions
	{
		if (simulation->parameters->boundary_y==0) //use periodic boundary conditions in y-direction
		{
			virtual_position_y-=simulation->state->length_y;
			L_interaction_with_neighbors_in_box(VM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index, &delta_v, &neighbor_n);
		}
		else //use reflecting boundary conditions in y-direction
			L_interaction_with_neighbors_in_box_mirrory(VM_get_box_index(virtual_position_x, yposition_focus, simulation->parameters), simulation, index, &delta_v, &neighbor_n);
	}
	else
		L_interaction_with_neighbors_in_box(VM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index, &delta_v, &neighbor_n);
	//left box
	virtual_position_y=yposition_focus;
	virtual_position_x=xposition_focus - simulation->parameters->box_size_x;
	if (virtual_position_x<0.0) //neighbor box is out of x-range -> apply boundary conditions
	{
		if (simulation->parameters->boundary_x==0) //use periodic boundary conditions in x-direction
		{
			virtual_position_x+=simulation->state->length_x;
			L_interaction_with_neighbors_in_box(VM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index, &delta_v, &neighbor_n);
		}
		else //use reflecting boundary conditions in x-direction
			L_interaction_with_neighbors_in_box_mirrorx(VM_get_box_index(xposition_focus, virtual_position_y, simulation->parameters), simulation, index, &delta_v, &neighbor_n);
	}
	else
		L_interaction_with_neighbors_in_box(VM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index, &delta_v, &neighbor_n);
	//left bottom box
	virtual_position_y=yposition_focus - simulation->parameters->box_size_y;
	if ((virtual_position_y<0.0) && !(virtual_position_x<0.0))
	{
		if (simulation->parameters->boundary_y==0)
		{
			virtual_position_y+=simulation->state->length_y;
			L_interaction_with_neighbors_in_box(VM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index, &delta_v, &neighbor_n);
		}
		else
			L_interaction_with_neighbors_in_box_mirrory(VM_get_box_index(virtual_position_x, yposition_focus, simulation->parameters), simulation, index, &delta_v, &neighbor_n);
	}
	else if (!(virtual_position_y<0.0) && (virtual_position_x<0.0))
	{
		if (simulation->parameters->boundary_x==0)
		{
			virtual_position_x+=simulation->state->length_x;
			L_interaction_with_neighbors_in_box(VM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index, &delta_v, &neighbor_n);
		}
		else
			L_interaction_with_neighbors_in_box_mirrorx(VM_get_box_index(xposition_focus, virtual_position_y, simulation->parameters), simulation, index, &delta_v, &neighbor_n);
	}
	else if ((virtual_position_y<0.0) && (virtual_position_x<0.0))
	{
		if ((simulation->parameters->boundary_x!=0) && (simulation->parameters->boundary_y!=0))
		{
			L_interaction_with_neighbors_in_box_mirrorxy(VM_get_box_index(xposition_focus, yposition_focus, simulation->parameters), simulation, index, &delta_v, &neighbor_n);
		}
		else if ((simulation->parameters->boundary_x==0) && (simulation->parameters->boundary_y!=0))
		{
			virtual_position_x+=simulation->state->length_x;
			L_interaction_with_neighbors_in_box_mirrory(VM_get_box_index(virtual_position_x, yposition_focus, simulation->parameters), simulation, index, &delta_v, &neighbor_n);
		}
		else if ((simulation->parameters->boundary_x!=0) && (simulation->parameters->boundary_y==0))
		{
			virtual_position_y+=simulation->state->length_y;
			L_interaction_with_neighbors_in_box_mirrorx(VM_get_box_index(xposition_focus, virtual_position_y, simulation->parameters), simulation, index, &delta_v, &neighbor_n);
		}
		else if ((simulation->parameters->boundary_x==0) && (simulation->parameters->boundary_y==0))
		{
			virtual_position_y+=simulation->state->length_y;
			virtual_position_x+=simulation->state->length_x;
			L_interaction_with_neighbors_in_box(VM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index, &delta_v, &neighbor_n);
		}
	}
	else
		L_interaction_with_neighbors_in_box(VM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index, &delta_v, &neighbor_n);
	//left top box
	virtual_position_y=yposition_focus + simulation->parameters->box_size_y;
	if ((virtual_position_y>simulation->state->length_y) && !(virtual_position_x<0.0))
	{
		if (simulation->parameters->boundary_y==0)
		{
			virtual_position_y-=simulation->state->length_y;
			L_interaction_with_neighbors_in_box(VM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index, &delta_v, &neighbor_n);
		}
		else
			L_interaction_with_neighbors_in_box_mirrory(VM_get_box_index(virtual_position_x, yposition_focus, simulation->parameters), simulation, index, &delta_v, &neighbor_n);
	}
	else if (!(virtual_position_y>simulation->state->length_y) && (virtual_position_x<0.0))
	{
		if (simulation->parameters->boundary_x==0)
		{
			virtual_position_x+=simulation->state->length_x;
			L_interaction_with_neighbors_in_box(VM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index, &delta_v, &neighbor_n);
		}
		else
			L_interaction_with_neighbors_in_box_mirrorx(VM_get_box_index(xposition_focus, virtual_position_y, simulation->parameters), simulation, index, &delta_v, &neighbor_n);
	}
	else if ((virtual_position_y>simulation->state->length_y) && (virtual_position_x<0.0))
	{
		if ((simulation->parameters->boundary_x!=0) && (simulation->parameters->boundary_y!=0))
		{
			L_interaction_with_neighbors_in_box_mirrorxy(VM_get_box_index(xposition_focus, yposition_focus, simulation->parameters), simulation, index, &delta_v, &neighbor_n);
		}
		else if ((simulation->parameters->boundary_x==0) && (simulation->parameters->boundary_y!=0))
		{
			virtual_position_x+=simulation->state->length_x;
			L_interaction_with_neighbors_in_box_mirrory(VM_get_box_index(virtual_position_x, yposition_focus, simulation->parameters), simulation, index, &delta_v, &neighbor_n);
		}
		else if ((simulation->parameters->boundary_x!=0) && (simulation->parameters->boundary_y==0))
		{
			virtual_position_y-=simulation->state->length_y;
			L_interaction_with_neighbors_in_box_mirrorx(VM_get_box_index(xposition_focus, virtual_position_y, simulation->parameters), simulation, index, &delta_v, &neighbor_n);
		}
		else if ((simulation->parameters->boundary_x==0) && (simulation->parameters->boundary_y==0))
		{
			virtual_position_y-=simulation->state->length_y;
			virtual_position_x+=simulation->state->length_x;
			L_interaction_with_neighbors_in_box(VM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index, &delta_v, &neighbor_n);
		}
	}
	else
		L_interaction_with_neighbors_in_box(VM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index, &delta_v, &neighbor_n);
	// right box
	virtual_position_y=yposition_focus;
	virtual_position_x=xposition_focus + simulation->parameters->box_size_x;
	if (virtual_position_x>simulation->state->length_x) //neighbor box is out of x-range -> apply boundary conditions
	{
		if (simulation->parameters->boundary_x==0) //use periodic boundary conditions in x-direction
		{
			virtual_position_x-=simulation->state->length_x;
			L_interaction_with_neighbors_in_box(VM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index, &delta_v, &neighbor_n);
		}
		else  //use reflecting boundary conditions in x-direction
			L_interaction_with_neighbors_in_box_mirrorx(VM_get_box_index(xposition_focus, virtual_position_y, simulation->parameters), simulation, index, &delta_v, &neighbor_n);
	}
	else
		L_interaction_with_neighbors_in_box(VM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index, &delta_v, &neighbor_n);
	// right bottom box
	virtual_position_y=yposition_focus - simulation->parameters->box_size_y;
	if ((virtual_position_y<0.0) && !(virtual_position_x>simulation->state->length_x))
	{
		if (simulation->parameters->boundary_y==0)
		{
			virtual_position_y+=simulation->state->length_y;
			L_interaction_with_neighbors_in_box(VM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index, &delta_v, &neighbor_n);
		}
		else
			L_interaction_with_neighbors_in_box_mirrory(VM_get_box_index(virtual_position_x, yposition_focus, simulation->parameters), simulation, index, &delta_v, &neighbor_n);
	}
	else if (!(virtual_position_y<0.0) && (virtual_position_x>simulation->state->length_x))
	{
		if (simulation->parameters->boundary_x==0)
		{
			virtual_position_x-=simulation->state->length_x;
			L_interaction_with_neighbors_in_box(VM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index, &delta_v, &neighbor_n);
		}
		else
			L_interaction_with_neighbors_in_box_mirrorx(VM_get_box_index(xposition_focus, virtual_position_y, simulation->parameters), simulation, index, &delta_v, &neighbor_n);
	}
	else if ((virtual_position_y<0.0) && (virtual_position_x>simulation->state->length_x))
	{
		if ((simulation->parameters->boundary_x!=0) && (simulation->parameters->boundary_y!=0))
		{
			L_interaction_with_neighbors_in_box_mirrorxy(VM_get_box_index(xposition_focus, yposition_focus, simulation->parameters), simulation, index, &delta_v, &neighbor_n);
		}
		else if ((simulation->parameters->boundary_x==0) && (simulation->parameters->boundary_y!=0))
		{
			virtual_position_x-=simulation->state->length_x;
			L_interaction_with_neighbors_in_box_mirrory(VM_get_box_index(virtual_position_x, yposition_focus, simulation->parameters), simulation, index, &delta_v, &neighbor_n);
		}
		else if ((simulation->parameters->boundary_x!=0) && (simulation->parameters->boundary_y==0))
		{
			virtual_position_y+=simulation->state->length_y;
			L_interaction_with_neighbors_in_box_mirrorx(VM_get_box_index(xposition_focus, virtual_position_y, simulation->parameters), simulation, index, &delta_v, &neighbor_n);
		}
		else if ((simulation->parameters->boundary_x==0) && (simulation->parameters->boundary_y==0))
		{
			virtual_position_y+=simulation->state->length_y;
			virtual_position_x-=simulation->state->length_x;
			L_interaction_with_neighbors_in_box(VM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index, &delta_v, &neighbor_n);
		}
	}
	else
		L_interaction_with_neighbors_in_box(VM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index, &delta_v, &neighbor_n);
	// right top box
	virtual_position_y=yposition_focus + simulation->parameters->box_size_y;
	if ((virtual_position_y>simulation->state->length_y) && !(virtual_position_x>simulation->state->length_x))
	{
		if (simulation->parameters->boundary_y==0)
		{
			virtual_position_y-=simulation->state->length_y;
			L_interaction_with_neighbors_in_box(VM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index, &delta_v, &neighbor_n);
		}
		else
			L_interaction_with_neighbors_in_box_mirrory(VM_get_box_index(virtual_position_x, yposition_focus, simulation->parameters), simulation, index, &delta_v, &neighbor_n);
	}
	else if (!(virtual_position_y>simulation->state->length_y) && (virtual_position_x>simulation->state->length_x))
	{
		if (simulation->parameters->boundary_x==0)
		{
			virtual_position_x-=simulation->state->length_x;
			L_interaction_with_neighbors_in_box(VM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index, &delta_v, &neighbor_n);
		}
		else
			L_interaction_with_neighbors_in_box_mirrorx(VM_get_box_index(xposition_focus, virtual_position_y, simulation->parameters), simulation, index, &delta_v, &neighbor_n);
	}
	else if ((virtual_position_y>simulation->state->length_y) && (virtual_position_x>simulation->state->length_x))
	{
		if ((simulation->parameters->boundary_x!=0) && (simulation->parameters->boundary_y!=0))
		{
			L_interaction_with_neighbors_in_box_mirrorxy(VM_get_box_index(xposition_focus, yposition_focus, simulation->parameters), simulation, index, &delta_v, &neighbor_n);
		}
		else if ((simulation->parameters->boundary_x==0) && (simulation->parameters->boundary_y!=0))
		{
			virtual_position_x-=simulation->state->length_x;
			L_interaction_with_neighbors_in_box_mirrory(VM_get_box_index(virtual_position_x, yposition_focus, simulation->parameters), simulation, index, &delta_v, &neighbor_n);
		}
		else if ((simulation->parameters->boundary_x!=0) && (simulation->parameters->boundary_y==0))
		{
			virtual_position_y-=simulation->state->length_y;
			L_interaction_with_neighbors_in_box_mirrorx(VM_get_box_index(xposition_focus, virtual_position_y, simulation->parameters), simulation, index, &delta_v, &neighbor_n);
		}
		else if ((simulation->parameters->boundary_x==0) && (simulation->parameters->boundary_y==0))
		{
			virtual_position_y-=simulation->state->length_y;
			virtual_position_x-=simulation->state->length_x;
			L_interaction_with_neighbors_in_box(VM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index, &delta_v, &neighbor_n);
		}
	}
	else
		L_interaction_with_neighbors_in_box(VM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index, &delta_v, &neighbor_n);
	//update the angle
	if (simulation->parameters->particle_species==1)
		if (simulation->state->weight_vector_length==0)
			(simulation->state->particles+index)->theta_new=(simulation->state->particles+index)->theta + simulation->parameters->delta_t*(*simulation->parameters->omega+delta_v/(neighbor_n+1)) +simulation->parameters->sqrt_delta_t**simulation->parameters->noise_strength*gauss(simulation);
		else
			(simulation->state->particles+index)->theta_new=(simulation->state->particles+index)->theta + simulation->parameters->delta_t*(*simulation->parameters->omega+delta_v**(simulation->state->weights+min(neighbor_n, simulation->state->weight_vector_length-1)) ) +simulation->parameters->sqrt_delta_t**simulation->parameters->noise_strength*gauss(simulation);
	else
		if (simulation->state->weight_vector_length==0)
			(simulation->state->particles+index)->theta_new=(simulation->state->particles+index)->theta + simulation->parameters->delta_t*(*(simulation->parameters->omega+get_species(index, simulation->parameters))+delta_v/(neighbor_n+1)) +simulation->parameters->sqrt_delta_t**(simulation->parameters->noise_strength+get_species(index, simulation->parameters))*gauss(simulation);
		else
			(simulation->state->particles+index)->theta_new=(simulation->state->particles+index)->theta + simulation->parameters->delta_t*(*(simulation->parameters->omega+get_species(index, simulation->parameters))+delta_v**(simulation->state->weights+min(neighbor_n, simulation->state->weight_vector_length-1)) ) +simulation->parameters->sqrt_delta_t**(simulation->parameters->noise_strength+get_species(index, simulation->parameters))*gauss(simulation);
	while ((simulation->state->particles+index)->theta_new<-M_PI)
		(simulation->state->particles+index)->theta_new+=2.0*M_PI;
	while ((simulation->state->particles+index)->theta_new>M_PI)
		(simulation->state->particles+index)->theta_new-=2.0*M_PI;
	return neighbor_n;
}

//checks all particles in a box to possibly include them in a neighbor list of particle 'index' for metric free/topological interactions
gint32 mfVM_sort_in_neighbors_in_box(gint32 box_index, struct VM_simulation * simulation, guint64 index)
{
	gdouble distance, distance_squared;
	gdouble xposition_focus=(simulation->state->particles+index)->x;
	gdouble yposition_focus=(simulation->state->particles+index)->y;
	struct VM_particle * box= *(simulation->parameters->mf_box+box_index);
	while (box!=NULL)
	{
		if (box->index != index)
		{
			distance_squared=0.0;
			distance=xposition_focus-(simulation->state->particles+box->index)->x;
			if (distance>0.5*simulation->state->length_x)
					distance -= simulation->state->length_x;
			if (distance<-0.5*simulation->state->length_x)
					distance += simulation->state->length_x;
			distance_squared += distance*distance;
			distance=yposition_focus-(simulation->state->particles+box->index)->y;
			if (distance>0.5*simulation->state->length_y)
					distance -= simulation->state->length_y;
			if (distance<-0.5*simulation->state->length_y)
					distance += simulation->state->length_y;
			distance_squared += distance*distance;
			mfVM_add_to_neighborhood(box->index, sqrt(distance_squared), 0, simulation);
		}
		box=box->next_particle;
	}
	return 0;
}

//checks mirror image under reflection at x=0,Lx of all particles in a box to possibly include them in a neighbor list of particle 'index' for metric free/topological interactions
gint32 mfVM_sort_in_neighbors_in_box_mirrorx(gint32 box_index, struct VM_simulation * simulation, guint64 index)
{
	gdouble distance, distance_squared;
	gdouble xposition_focus=(simulation->state->particles+index)->x;
	gdouble yposition_focus=(simulation->state->particles+index)->y;
	struct VM_particle * box= *(simulation->parameters->mf_box+box_index);
	while (box!=NULL)
	{
		if (box->index != index)
		{
			distance_squared=0.0;
			distance=xposition_focus+(simulation->state->particles+box->index)->x;
			while (distance>0.5*simulation->state->length_x)
					distance -= simulation->state->length_x;
			while (distance<-0.5*simulation->state->length_x)
					distance += simulation->state->length_x;
			distance_squared += distance*distance;
			distance=yposition_focus-(simulation->state->particles+box->index)->y;
			if (distance>0.5*simulation->state->length_y)
					distance -= simulation->state->length_y;
			if (distance<-0.5*simulation->state->length_y)
					distance += simulation->state->length_y;
			distance_squared += distance*distance;
			mfVM_add_to_neighborhood(box->index, sqrt(distance_squared), 1, simulation);
		}
		box=box->next_particle;
	}
	return 0;
}

//checks mirror image under reflection at y=0,Ly of all particles in a box to possibly include them in a neighbor list of particle 'index' for metric free/topological interactions
gint32 mfVM_sort_in_neighbors_in_box_mirrory(gint32 box_index, struct VM_simulation * simulation, guint64 index)
{
	gdouble distance, distance_squared;
	gdouble xposition_focus=(simulation->state->particles+index)->x;
	gdouble yposition_focus=(simulation->state->particles+index)->y;
	struct VM_particle * box= *(simulation->parameters->mf_box+box_index);
	while (box!=NULL)
	{
		if (box->index != index)
		{
			distance_squared=0.0;
			distance=xposition_focus-(simulation->state->particles+box->index)->x;
			if (distance>0.5*simulation->state->length_x)
					distance -= simulation->state->length_x;
			if (distance<-0.5*simulation->state->length_x)
					distance += simulation->state->length_x;
			distance_squared += distance*distance;
			distance=yposition_focus+(simulation->state->particles+box->index)->y;
			while (distance>0.5*simulation->state->length_y)
					distance -= simulation->state->length_y;
			while (distance<-0.5*simulation->state->length_y)
					distance += simulation->state->length_y;
			distance_squared += distance*distance;
			mfVM_add_to_neighborhood(box->index, sqrt(distance_squared), 2, simulation);
		}
		box=box->next_particle;
	}
	return 0;
}

//checks mirror image under reflection at both axis of all particles in a box to possibly include them in a neighbor list of particle 'index' for metric free/topological interactions
gint32 mfVM_sort_in_neighbors_in_box_mirrorxy(gint32 box_index, struct VM_simulation * simulation, guint64 index)
{
	gdouble distance, distance_squared;
	gdouble xposition_focus=(simulation->state->particles+index)->x;
	gdouble yposition_focus=(simulation->state->particles+index)->y;
	struct VM_particle * box= *(simulation->parameters->mf_box+box_index);
	while (box!=NULL)
	{
		if (box->index != index)
		{
			distance_squared=0.0;
			distance=xposition_focus+(simulation->state->particles+box->index)->x;
			while (distance>0.5*simulation->state->length_x)
					distance -= simulation->state->length_x;
			while (distance<-0.5*simulation->state->length_x)
					distance += simulation->state->length_x;
			distance_squared += distance*distance;
			distance=yposition_focus+(simulation->state->particles+box->index)->y;
			while (distance>0.5*simulation->state->length_y)
					distance -= simulation->state->length_y;
			while (distance<-0.5*simulation->state->length_y)
					distance += simulation->state->length_y;
			distance_squared += distance*distance;
			mfVM_add_to_neighborhood(box->index, sqrt(distance_squared), 3, simulation);
		}
		box=box->next_particle;
	}
	return 0;
}

//calculates interactions of particle 'index' for metric free/topological Vicsek model
gint32 mfVM_interaction_with_neighbors(guint64 index, struct VM_simulation * simulation)
{
	mfVM_clean_neighborhood(simulation);
	gint32 order=1; // this says which boxes are checked: if order=1 only the inner box is checked, if order=2 the surrounding 8 boxes are checked and so on
	gint32 i, j;
	gint32 kx, ky;
	gdouble velocity_x=cos((simulation->state->particles+index)->theta);
	gdouble velocity_y=sin((simulation->state->particles+index)->theta);
	gdouble virtual_position_x;
	gdouble virtual_position_y;
	gdouble position_x=(simulation->state->particles+index)->x;
	gdouble position_y=(simulation->state->particles+index)->y;
	gdouble dx=0.0;
	gdouble dy=0.0;
	//find mf_kn nearest neighbors of particle 'index'
	while (simulation->state->nh.max_distance_checked<*(simulation->state->nh.distances+simulation->parameters->mf_kn-1))
	{
		if (order==1)
		{
			i=0;
			j=0;
			virtual_position_x=(simulation->state->particles+index)->x+i*simulation->parameters->mf_box_size_x;
			if (virtual_position_x<0.0)
				virtual_position_x+=simulation->state->length_x;
			if (virtual_position_x>simulation->state->length_x)
				virtual_position_x-=simulation->state->length_x;
			virtual_position_y=(simulation->state->particles+index)->y+j*simulation->parameters->mf_box_size_y;
			if (virtual_position_y<0.0)
				virtual_position_y+=simulation->state->length_y;
			if (virtual_position_y>simulation->state->length_y)
				virtual_position_y-=simulation->state->length_y;
			mfVM_sort_in_neighbors_in_box(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
		}
		else
		{
			i=1-order;
			for (j=1-order; j<order; j++)
			{
				virtual_position_x=(simulation->state->particles+index)->x+i*simulation->parameters->mf_box_size_x;
				virtual_position_y=(simulation->state->particles+index)->y+j*simulation->parameters->mf_box_size_y;
				if ((virtual_position_x<0.0) && (virtual_position_y<0.0)) //out of simulation box to left bottom
				{
					if ((simulation->parameters->boundary_x==0)&& (simulation->parameters->boundary_y==0)) //periodic boundary conditions in both directions
					{
						virtual_position_x+=simulation->state->length_x;
						virtual_position_y+=simulation->state->length_y;
						mfVM_sort_in_neighbors_in_box(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else if((simulation->parameters->boundary_x==0)&& (simulation->parameters->boundary_y!=0))//periodic in x-direction and reflecting in y-direction
					{
						virtual_position_x+=simulation->state->length_x;
						ky=-1;
						while (virtual_position_y<0.0)
						{
							virtual_position_y+=simulation->parameters->mf_box_size_y;
							ky+=1;
						}
						while (ky>0)
						{
							virtual_position_y+=simulation->parameters->mf_box_size_y;
							ky-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrory(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else if((simulation->parameters->boundary_x!=0)&& (simulation->parameters->boundary_y==0)) //reflecting in x-direction, periodic in y-direction
					{
						virtual_position_y+=simulation->state->length_y;
						kx=-1;
						while (virtual_position_x<0.0)
						{
							virtual_position_x+=simulation->parameters->mf_box_size_x;
							kx+=1;
						}
						while (kx>0)
						{
							virtual_position_x+=simulation->parameters->mf_box_size_x;
							kx-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrorx(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else if((simulation->parameters->boundary_x!=0)&& (simulation->parameters->boundary_y!=0)) //reflecting boundary in both directions
					{
						kx=-1;
						while (virtual_position_x<0.0)
						{
							virtual_position_x+=simulation->parameters->mf_box_size_x;
							kx+=1;
						}
						while (kx>0)
						{
							virtual_position_x+=simulation->parameters->mf_box_size_x;
							kx-=1;
						}
						ky=-1;
						while (virtual_position_y<0.0)
						{
							virtual_position_y+=simulation->parameters->mf_box_size_y;
							ky+=1;
						}
						while (ky>0)
						{
							virtual_position_y+=simulation->parameters->mf_box_size_y;
							ky-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrorxy(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
				}
				else if ((virtual_position_x<0.0) && (virtual_position_y>simulation->state->length_y)) //out of simulation box to left top
				{
					if ((simulation->parameters->boundary_x==0)&& (simulation->parameters->boundary_y==0)) //periodic boundary conditions in both directions
					{
						virtual_position_x+=simulation->state->length_x;
						virtual_position_y-=simulation->state->length_y;
						mfVM_sort_in_neighbors_in_box(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else if((simulation->parameters->boundary_x==0)&& (simulation->parameters->boundary_y!=0))//periodic in x-direction and reflecting in y-direction
					{
						virtual_position_x+=simulation->state->length_x;
						ky=-1;
						while (virtual_position_y>simulation->state->length_y)
						{
							virtual_position_y-=simulation->parameters->mf_box_size_y;
							ky+=1;
						}
						while (ky>0)
						{
							virtual_position_y-=simulation->parameters->mf_box_size_y;
							ky-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrory(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else if((simulation->parameters->boundary_x!=0)&& (simulation->parameters->boundary_y==0)) //reflecting in x-direction, periodic in y-direction
					{
						virtual_position_y-=simulation->state->length_y;
						kx=-1;
						while (virtual_position_x<0.0)
						{
							virtual_position_x+=simulation->parameters->mf_box_size_x;
							kx+=1;
						}
						while (kx>0)
						{
							virtual_position_x+=simulation->parameters->mf_box_size_x;
							kx-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrorx(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else if((simulation->parameters->boundary_x!=0)&& (simulation->parameters->boundary_y!=0)) //reflecting boundary in both directions
					{
						kx=-1;
						while (virtual_position_x<0.0)
						{
							virtual_position_x+=simulation->parameters->mf_box_size_x;
							kx+=1;
						}
						while (kx>0)
						{
							virtual_position_x+=simulation->parameters->mf_box_size_x;
							kx-=1;
						}
						ky=-1;
						while (virtual_position_y>simulation->state->length_y)
						{
							virtual_position_y-=simulation->parameters->mf_box_size_y;
							ky+=1;
						}
						while (ky>0)
						{
							virtual_position_y-=simulation->parameters->mf_box_size_y;
							ky-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrorxy(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
				}
				else if ((virtual_position_x>simulation->state->length_x) && (virtual_position_y>simulation->state->length_y)) //out of simulation box to right top
				{
					if ((simulation->parameters->boundary_x==0)&& (simulation->parameters->boundary_y==0)) //periodic boundary conditions in both directions
					{
						virtual_position_x-=simulation->state->length_x;
						virtual_position_y-=simulation->state->length_y;
						mfVM_sort_in_neighbors_in_box(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else if((simulation->parameters->boundary_x==0)&& (simulation->parameters->boundary_y!=0))//periodic in x-direction and reflecting in y-direction
					{
						virtual_position_x-=simulation->state->length_x;
						ky=-1;
						while (virtual_position_y>simulation->state->length_y)
						{
							virtual_position_y-=simulation->parameters->mf_box_size_y;
							ky+=1;
						}
						while (ky>0)
						{
							virtual_position_y-=simulation->parameters->mf_box_size_y;
							ky-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrory(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else if((simulation->parameters->boundary_x!=0)&& (simulation->parameters->boundary_y==0)) //reflecting in x-direction, periodic in y-direction
					{
						virtual_position_y-=simulation->state->length_y;
						kx=-1;
						while (virtual_position_x>simulation->state->length_x)
						{
							virtual_position_x-=simulation->parameters->mf_box_size_x;
							kx+=1;
						}
						while (kx>0)
						{
							virtual_position_x-=simulation->parameters->mf_box_size_x;
							kx-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrorx(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else if((simulation->parameters->boundary_x!=0)&& (simulation->parameters->boundary_y!=0)) //reflecting boundary in both directions
					{
						kx=-1;
						while (virtual_position_x>simulation->state->length_x)
						{
							virtual_position_x-=simulation->parameters->mf_box_size_x;
							kx+=1;
						}
						while (kx>0)
						{
							virtual_position_x-=simulation->parameters->mf_box_size_x;
							kx-=1;
						}
						ky=-1;
						while (virtual_position_y>simulation->state->length_y)
						{
							virtual_position_y-=simulation->parameters->mf_box_size_y;
							ky+=1;
						}
						while (ky>0)
						{
							virtual_position_y-=simulation->parameters->mf_box_size_y;
							ky-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrorxy(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
				}
				else if ((virtual_position_x>simulation->state->length_x) && (virtual_position_y<0.0)) //out of simulation box to right bottom
				{
					if ((simulation->parameters->boundary_x==0)&& (simulation->parameters->boundary_y==0)) //periodic boundary conditions in both directions
					{
						virtual_position_x-=simulation->state->length_x;
						virtual_position_y+=simulation->state->length_y;
						mfVM_sort_in_neighbors_in_box(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else if((simulation->parameters->boundary_x==0)&& (simulation->parameters->boundary_y!=0))//periodic in x-direction and reflecting in y-direction
					{
						virtual_position_x-=simulation->state->length_x;
						ky=-1;
						while (virtual_position_y<0.0)
						{
							virtual_position_y+=simulation->parameters->mf_box_size_y;
							ky+=1;
						}
						while (ky>0)
						{
							virtual_position_y+=simulation->parameters->mf_box_size_y;
							ky-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrory(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else if((simulation->parameters->boundary_x!=0)&& (simulation->parameters->boundary_y==0)) //reflecting in x-direction, periodic in y-direction
					{
						virtual_position_y+=simulation->state->length_y;
						kx=-1;
						while (virtual_position_x>simulation->state->length_x)
						{
							virtual_position_x-=simulation->parameters->mf_box_size_x;
							kx+=1;
						}
						while (kx>0)
						{
							virtual_position_x-=simulation->parameters->mf_box_size_x;
							kx-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrorx(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else if((simulation->parameters->boundary_x!=0)&& (simulation->parameters->boundary_y!=0)) //reflecting boundary in both directions
					{
						kx=-1;
						while (virtual_position_x>simulation->state->length_x)
						{
							virtual_position_x-=simulation->parameters->mf_box_size_x;
							kx+=1;
						}
						while (kx>0)
						{
							virtual_position_x-=simulation->parameters->mf_box_size_x;
							kx-=1;
						}
						ky=-1;
						while (virtual_position_y<0.0)
						{
							virtual_position_y+=simulation->parameters->mf_box_size_y;
							ky+=1;
						}
						while (ky>0)
						{
							virtual_position_y+=simulation->parameters->mf_box_size_y;
							ky-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrorxy(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
				}
				else if (virtual_position_x<0.0) //out of simulation box to the left
				{
					if (simulation->parameters->boundary_x==0)
					{
						virtual_position_x+=simulation->state->length_x;
						mfVM_sort_in_neighbors_in_box(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else
					{
						kx=-1;
						while (virtual_position_x<0.0)
						{
							virtual_position_x+=simulation->parameters->mf_box_size_x;
							kx+=1;
						}
						while (kx>0)
						{
							virtual_position_x+=simulation->parameters->mf_box_size_x;
							kx-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrorx(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
				}
				else if (virtual_position_x>simulation->state->length_x) //out of simulation box to the right
				{
					if (simulation->parameters->boundary_x==0)
					{
						virtual_position_x-=simulation->state->length_x;
						mfVM_sort_in_neighbors_in_box(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else
					{
						kx=-1;
						while (virtual_position_x>simulation->state->length_x)
						{
							virtual_position_x-=simulation->parameters->mf_box_size_x;
							kx+=1;
						}
						while (kx>0)
						{
							virtual_position_x-=simulation->parameters->mf_box_size_x;
							kx-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrorx(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
				}
				else if (virtual_position_y<0.0) //out of simulation box to the bottom
				{
					if (simulation->parameters->boundary_y==0)
					{
						virtual_position_y+=simulation->state->length_y;
						mfVM_sort_in_neighbors_in_box(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else
					{
						ky=-1;
						while (virtual_position_y<0.0)
						{
							virtual_position_y+=simulation->parameters->mf_box_size_y;
							ky+=1;
						}
						while (ky>0)
						{
							virtual_position_y+=simulation->parameters->mf_box_size_y;
							ky-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrory(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
				}
				else if (virtual_position_y>simulation->state->length_y) //out of simulation box to the top
				{
					if (simulation->parameters->boundary_y==0)
					{
						virtual_position_y-=simulation->state->length_y;
						mfVM_sort_in_neighbors_in_box(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else
					{
						ky=-1;
						while (virtual_position_y>simulation->state->length_y)
						{
							virtual_position_y-=simulation->parameters->mf_box_size_y;
							ky+=1;
						}
						while (ky>0)
						{
							virtual_position_y-=simulation->parameters->mf_box_size_y;
							ky-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrorx(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
				}
				else  // virtual position still in simulation box
				{
					mfVM_sort_in_neighbors_in_box(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
				}
			}
			i=order-1;
			for (j=1-order; j<order; j++)
			{
				virtual_position_x=(simulation->state->particles+index)->x+i*simulation->parameters->mf_box_size_x;
				virtual_position_y=(simulation->state->particles+index)->y+j*simulation->parameters->mf_box_size_y;
				if ((virtual_position_x<0.0) && (virtual_position_y<0.0)) //out of simulation box to left bottom
				{
					if ((simulation->parameters->boundary_x==0)&& (simulation->parameters->boundary_y==0)) //periodic boundary conditions in both directions
					{
						virtual_position_x+=simulation->state->length_x;
						virtual_position_y+=simulation->state->length_y;
						mfVM_sort_in_neighbors_in_box(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else if((simulation->parameters->boundary_x==0)&& (simulation->parameters->boundary_y!=0))//periodic in x-direction and reflecting in y-direction
					{
						virtual_position_x+=simulation->state->length_x;
						ky=-1;
						while (virtual_position_y<0.0)
						{
							virtual_position_y+=simulation->parameters->mf_box_size_y;
							ky+=1;
						}
						while (ky>0)
						{
							virtual_position_y+=simulation->parameters->mf_box_size_y;
							ky-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrory(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else if((simulation->parameters->boundary_x!=0)&& (simulation->parameters->boundary_y==0)) //reflecting in x-direction, periodic in y-direction
					{
						virtual_position_y+=simulation->state->length_y;
						kx=-1;
						while (virtual_position_x<0.0)
						{
							virtual_position_x+=simulation->parameters->mf_box_size_x;
							kx+=1;
						}
						while (kx>0)
						{
							virtual_position_x+=simulation->parameters->mf_box_size_x;
							kx-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrorx(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else if((simulation->parameters->boundary_x!=0)&& (simulation->parameters->boundary_y!=0)) //reflecting boundary in both directions
					{
						kx=-1;
						while (virtual_position_x<0.0)
						{
							virtual_position_x+=simulation->parameters->mf_box_size_x;
							kx+=1;
						}
						while (kx>0)
						{
							virtual_position_x+=simulation->parameters->mf_box_size_x;
							kx-=1;
						}
						ky=-1;
						while (virtual_position_y<0.0)
						{
							virtual_position_y+=simulation->parameters->mf_box_size_y;
							ky+=1;
						}
						while (ky>0)
						{
							virtual_position_y+=simulation->parameters->mf_box_size_y;
							ky-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrorxy(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
				}
				else if ((virtual_position_x<0.0) && (virtual_position_y>simulation->state->length_y)) //out of simulation box to left top
				{
					if ((simulation->parameters->boundary_x==0)&& (simulation->parameters->boundary_y==0)) //periodic boundary conditions in both directions
					{
						virtual_position_x+=simulation->state->length_x;
						virtual_position_y-=simulation->state->length_y;
						mfVM_sort_in_neighbors_in_box(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else if((simulation->parameters->boundary_x==0)&& (simulation->parameters->boundary_y!=0))//periodic in x-direction and reflecting in y-direction
					{
						virtual_position_x+=simulation->state->length_x;
						ky=-1;
						while (virtual_position_y>simulation->state->length_y)
						{
							virtual_position_y-=simulation->parameters->mf_box_size_y;
							ky+=1;
						}
						while (ky>0)
						{
							virtual_position_y-=simulation->parameters->mf_box_size_y;
							ky-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrory(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else if((simulation->parameters->boundary_x!=0)&& (simulation->parameters->boundary_y==0)) //reflecting in x-direction, periodic in y-direction
					{
						virtual_position_y-=simulation->state->length_y;
						kx=-1;
						while (virtual_position_x<0.0)
						{
							virtual_position_x+=simulation->parameters->mf_box_size_x;
							kx+=1;
						}
						while (kx>0)
						{
							virtual_position_x+=simulation->parameters->mf_box_size_x;
							kx-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrorx(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else if((simulation->parameters->boundary_x!=0)&& (simulation->parameters->boundary_y!=0)) //reflecting boundary in both directions
					{
						kx=-1;
						while (virtual_position_x<0.0)
						{
							virtual_position_x+=simulation->parameters->mf_box_size_x;
							kx+=1;
						}
						while (kx>0)
						{
							virtual_position_x+=simulation->parameters->mf_box_size_x;
							kx-=1;
						}
						ky=-1;
						while (virtual_position_y>simulation->state->length_y)
						{
							virtual_position_y-=simulation->parameters->mf_box_size_y;
							ky+=1;
						}
						while (ky>0)
						{
							virtual_position_y-=simulation->parameters->mf_box_size_y;
							ky-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrorxy(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
				}
				else if ((virtual_position_x>simulation->state->length_x) && (virtual_position_y>simulation->state->length_y)) //out of simulation box to right top
				{
					if ((simulation->parameters->boundary_x==0)&& (simulation->parameters->boundary_y==0)) //periodic boundary conditions in both directions
					{
						virtual_position_x-=simulation->state->length_x;
						virtual_position_y-=simulation->state->length_y;
						mfVM_sort_in_neighbors_in_box(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else if((simulation->parameters->boundary_x==0)&& (simulation->parameters->boundary_y!=0))//periodic in x-direction and reflecting in y-direction
					{
						virtual_position_x-=simulation->state->length_x;
						ky=-1;
						while (virtual_position_y>simulation->state->length_y)
						{
							virtual_position_y-=simulation->parameters->mf_box_size_y;
							ky+=1;
						}
						while (ky>0)
						{
							virtual_position_y-=simulation->parameters->mf_box_size_y;
							ky-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrory(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else if((simulation->parameters->boundary_x!=0)&& (simulation->parameters->boundary_y==0)) //reflecting in x-direction, periodic in y-direction
					{
						virtual_position_y-=simulation->state->length_y;
						kx=-1;
						while (virtual_position_x>simulation->state->length_x)
						{
							virtual_position_x-=simulation->parameters->mf_box_size_x;
							kx+=1;
						}
						while (kx>0)
						{
							virtual_position_x-=simulation->parameters->mf_box_size_x;
							kx-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrorx(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else if((simulation->parameters->boundary_x!=0)&& (simulation->parameters->boundary_y!=0)) //reflecting boundary in both directions
					{
						kx=-1;
						while (virtual_position_x>simulation->state->length_x)
						{
							virtual_position_x-=simulation->parameters->mf_box_size_x;
							kx+=1;
						}
						while (kx>0)
						{
							virtual_position_x-=simulation->parameters->mf_box_size_x;
							kx-=1;
						}
						ky=-1;
						while (virtual_position_y>simulation->state->length_y)
						{
							virtual_position_y-=simulation->parameters->mf_box_size_y;
							ky+=1;
						}
						while (ky>0)
						{
							virtual_position_y-=simulation->parameters->mf_box_size_y;
							ky-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrorxy(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
				}
				else if ((virtual_position_x>simulation->state->length_x) && (virtual_position_y<0.0)) //out of simulation box to right bottom
				{
					if ((simulation->parameters->boundary_x==0)&& (simulation->parameters->boundary_y==0)) //periodic boundary conditions in both directions
					{
						virtual_position_x-=simulation->state->length_x;
						virtual_position_y+=simulation->state->length_y;
						mfVM_sort_in_neighbors_in_box(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else if((simulation->parameters->boundary_x==0)&& (simulation->parameters->boundary_y!=0))//periodic in x-direction and reflecting in y-direction
					{
						virtual_position_x-=simulation->state->length_x;
						ky=-1;
						while (virtual_position_y<0.0)
						{
							virtual_position_y+=simulation->parameters->mf_box_size_y;
							ky+=1;
						}
						while (ky>0)
						{
							virtual_position_y+=simulation->parameters->mf_box_size_y;
							ky-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrory(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else if((simulation->parameters->boundary_x!=0)&& (simulation->parameters->boundary_y==0)) //reflecting in x-direction, periodic in y-direction
					{
						virtual_position_y+=simulation->state->length_y;
						kx=-1;
						while (virtual_position_x>simulation->state->length_x)
						{
							virtual_position_x-=simulation->parameters->mf_box_size_x;
							kx+=1;
						}
						while (kx>0)
						{
							virtual_position_x-=simulation->parameters->mf_box_size_x;
							kx-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrorx(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else if((simulation->parameters->boundary_x!=0)&& (simulation->parameters->boundary_y!=0)) //reflecting boundary in both directions
					{
						kx=-1;
						while (virtual_position_x>simulation->state->length_x)
						{
							virtual_position_x-=simulation->parameters->mf_box_size_x;
							kx+=1;
						}
						while (kx>0)
						{
							virtual_position_x-=simulation->parameters->mf_box_size_x;
							kx-=1;
						}
						ky=-1;
						while (virtual_position_y<0.0)
						{
							virtual_position_y+=simulation->parameters->mf_box_size_y;
							ky+=1;
						}
						while (ky>0)
						{
							virtual_position_y+=simulation->parameters->mf_box_size_y;
							ky-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrorxy(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
				}
				else if (virtual_position_x<0.0) //out of simulation box to the left
				{
					if (simulation->parameters->boundary_x==0)
					{
						virtual_position_x+=simulation->state->length_x;
						mfVM_sort_in_neighbors_in_box(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else
					{
						kx=-1;
						while (virtual_position_x<0.0)
						{
							virtual_position_x+=simulation->parameters->mf_box_size_x;
							kx+=1;
						}
						while (kx>0)
						{
							virtual_position_x+=simulation->parameters->mf_box_size_x;
							kx-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrorx(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
				}
				else if (virtual_position_x>simulation->state->length_x) //out of simulation box to the right
				{
					if (simulation->parameters->boundary_x==0)
					{
						virtual_position_x-=simulation->state->length_x;
						mfVM_sort_in_neighbors_in_box(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else
					{
						kx=-1;
						while (virtual_position_x>simulation->state->length_x)
						{
							virtual_position_x-=simulation->parameters->mf_box_size_x;
							kx+=1;
						}
						while (kx>0)
						{
							virtual_position_x-=simulation->parameters->mf_box_size_x;
							kx-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrorx(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
				}
				else if (virtual_position_y<0.0) //out of simulation box to the bottom
				{
					if (simulation->parameters->boundary_y==0)
					{
						virtual_position_y+=simulation->state->length_y;
						mfVM_sort_in_neighbors_in_box(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else
					{
						ky=-1;
						while (virtual_position_y<0.0)
						{
							virtual_position_y+=simulation->parameters->mf_box_size_y;
							ky+=1;
						}
						while (ky>0)
						{
							virtual_position_y+=simulation->parameters->mf_box_size_y;
							ky-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrory(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
				}
				else if (virtual_position_y>simulation->state->length_y) //out of simulation box to the top
				{
					if (simulation->parameters->boundary_y==0)
					{
						virtual_position_y-=simulation->state->length_y;
						mfVM_sort_in_neighbors_in_box(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else
					{
						ky=-1;
						while (virtual_position_y>simulation->state->length_y)
						{
							virtual_position_y-=simulation->parameters->mf_box_size_y;
							ky+=1;
						}
						while (ky>0)
						{
							virtual_position_y-=simulation->parameters->mf_box_size_y;
							ky-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrorx(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
				}
				else  // virtual position still in simulation box
				{
					mfVM_sort_in_neighbors_in_box(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
				}
			}
			j=1-order;
			for (i=2-order; i<order-1; i++)
			{
				virtual_position_x=(simulation->state->particles+index)->x+i*simulation->parameters->mf_box_size_x;
				virtual_position_y=(simulation->state->particles+index)->y+j*simulation->parameters->mf_box_size_y;
				if ((virtual_position_x<0.0) && (virtual_position_y<0.0)) //out of simulation box to left bottom
				{
					if ((simulation->parameters->boundary_x==0)&& (simulation->parameters->boundary_y==0)) //periodic boundary conditions in both directions
					{
						virtual_position_x+=simulation->state->length_x;
						virtual_position_y+=simulation->state->length_y;
						mfVM_sort_in_neighbors_in_box(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else if((simulation->parameters->boundary_x==0)&& (simulation->parameters->boundary_y!=0))//periodic in x-direction and reflecting in y-direction
					{
						virtual_position_x+=simulation->state->length_x;
						ky=-1;
						while (virtual_position_y<0.0)
						{
							virtual_position_y+=simulation->parameters->mf_box_size_y;
							ky+=1;
						}
						while (ky>0)
						{
							virtual_position_y+=simulation->parameters->mf_box_size_y;
							ky-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrory(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else if((simulation->parameters->boundary_x!=0)&& (simulation->parameters->boundary_y==0)) //reflecting in x-direction, periodic in y-direction
					{
						virtual_position_y+=simulation->state->length_y;
						kx=-1;
						while (virtual_position_x<0.0)
						{
							virtual_position_x+=simulation->parameters->mf_box_size_x;
							kx+=1;
						}
						while (kx>0)
						{
							virtual_position_x+=simulation->parameters->mf_box_size_x;
							kx-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrorx(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else if((simulation->parameters->boundary_x!=0)&& (simulation->parameters->boundary_y!=0)) //reflecting boundary in both directions
					{
						kx=-1;
						while (virtual_position_x<0.0)
						{
							virtual_position_x+=simulation->parameters->mf_box_size_x;
							kx+=1;
						}
						while (kx>0)
						{
							virtual_position_x+=simulation->parameters->mf_box_size_x;
							kx-=1;
						}
						ky=-1;
						while (virtual_position_y<0.0)
						{
							virtual_position_y+=simulation->parameters->mf_box_size_y;
							ky+=1;
						}
						while (ky>0)
						{
							virtual_position_y+=simulation->parameters->mf_box_size_y;
							ky-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrorxy(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
				}
				else if ((virtual_position_x<0.0) && (virtual_position_y>simulation->state->length_y)) //out of simulation box to left top
				{
					if ((simulation->parameters->boundary_x==0)&& (simulation->parameters->boundary_y==0)) //periodic boundary conditions in both directions
					{
						virtual_position_x+=simulation->state->length_x;
						virtual_position_y-=simulation->state->length_y;
						mfVM_sort_in_neighbors_in_box(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else if((simulation->parameters->boundary_x==0)&& (simulation->parameters->boundary_y!=0))//periodic in x-direction and reflecting in y-direction
					{
						virtual_position_x+=simulation->state->length_x;
						ky=-1;
						while (virtual_position_y>simulation->state->length_y)
						{
							virtual_position_y-=simulation->parameters->mf_box_size_y;
							ky+=1;
						}
						while (ky>0)
						{
							virtual_position_y-=simulation->parameters->mf_box_size_y;
							ky-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrory(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else if((simulation->parameters->boundary_x!=0)&& (simulation->parameters->boundary_y==0)) //reflecting in x-direction, periodic in y-direction
					{
						virtual_position_y-=simulation->state->length_y;
						kx=-1;
						while (virtual_position_x<0.0)
						{
							virtual_position_x+=simulation->parameters->mf_box_size_x;
							kx+=1;
						}
						while (kx>0)
						{
							virtual_position_x+=simulation->parameters->mf_box_size_x;
							kx-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrorx(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else if((simulation->parameters->boundary_x!=0)&& (simulation->parameters->boundary_y!=0)) //reflecting boundary in both directions
					{
						kx=-1;
						while (virtual_position_x<0.0)
						{
							virtual_position_x+=simulation->parameters->mf_box_size_x;
							kx+=1;
						}
						while (kx>0)
						{
							virtual_position_x+=simulation->parameters->mf_box_size_x;
							kx-=1;
						}
						ky=-1;
						while (virtual_position_y>simulation->state->length_y)
						{
							virtual_position_y-=simulation->parameters->mf_box_size_y;
							ky+=1;
						}
						while (ky>0)
						{
							virtual_position_y-=simulation->parameters->mf_box_size_y;
							ky-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrorxy(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
				}
				else if ((virtual_position_x>simulation->state->length_x) && (virtual_position_y>simulation->state->length_y)) //out of simulation box to right top
				{
					if ((simulation->parameters->boundary_x==0)&& (simulation->parameters->boundary_y==0)) //periodic boundary conditions in both directions
					{
						virtual_position_x-=simulation->state->length_x;
						virtual_position_y-=simulation->state->length_y;
						mfVM_sort_in_neighbors_in_box(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else if((simulation->parameters->boundary_x==0)&& (simulation->parameters->boundary_y!=0))//periodic in x-direction and reflecting in y-direction
					{
						virtual_position_x-=simulation->state->length_x;
						ky=-1;
						while (virtual_position_y>simulation->state->length_y)
						{
							virtual_position_y-=simulation->parameters->mf_box_size_y;
							ky+=1;
						}
						while (ky>0)
						{
							virtual_position_y-=simulation->parameters->mf_box_size_y;
							ky-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrory(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else if((simulation->parameters->boundary_x!=0)&& (simulation->parameters->boundary_y==0)) //reflecting in x-direction, periodic in y-direction
					{
						virtual_position_y-=simulation->state->length_y;
						kx=-1;
						while (virtual_position_x>simulation->state->length_x)
						{
							virtual_position_x-=simulation->parameters->mf_box_size_x;
							kx+=1;
						}
						while (kx>0)
						{
							virtual_position_x-=simulation->parameters->mf_box_size_x;
							kx-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrorx(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else if((simulation->parameters->boundary_x!=0)&& (simulation->parameters->boundary_y!=0)) //reflecting boundary in both directions
					{
						kx=-1;
						while (virtual_position_x>simulation->state->length_x)
						{
							virtual_position_x-=simulation->parameters->mf_box_size_x;
							kx+=1;
						}
						while (kx>0)
						{
							virtual_position_x-=simulation->parameters->mf_box_size_x;
							kx-=1;
						}
						ky=-1;
						while (virtual_position_y>simulation->state->length_y)
						{
							virtual_position_y-=simulation->parameters->mf_box_size_y;
							ky+=1;
						}
						while (ky>0)
						{
							virtual_position_y-=simulation->parameters->mf_box_size_y;
							ky-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrorxy(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
				}
				else if ((virtual_position_x>simulation->state->length_x) && (virtual_position_y<0.0)) //out of simulation box to right bottom
				{
					if ((simulation->parameters->boundary_x==0)&& (simulation->parameters->boundary_y==0)) //periodic boundary conditions in both directions
					{
						virtual_position_x-=simulation->state->length_x;
						virtual_position_y+=simulation->state->length_y;
						mfVM_sort_in_neighbors_in_box(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else if((simulation->parameters->boundary_x==0)&& (simulation->parameters->boundary_y!=0))//periodic in x-direction and reflecting in y-direction
					{
						virtual_position_x-=simulation->state->length_x;
						ky=-1;
						while (virtual_position_y<0.0)
						{
							virtual_position_y+=simulation->parameters->mf_box_size_y;
							ky+=1;
						}
						while (ky>0)
						{
							virtual_position_y+=simulation->parameters->mf_box_size_y;
							ky-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrory(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else if((simulation->parameters->boundary_x!=0)&& (simulation->parameters->boundary_y==0)) //reflecting in x-direction, periodic in y-direction
					{
						virtual_position_y+=simulation->state->length_y;
						kx=-1;
						while (virtual_position_x>simulation->state->length_x)
						{
							virtual_position_x-=simulation->parameters->mf_box_size_x;
							kx+=1;
						}
						while (kx>0)
						{
							virtual_position_x-=simulation->parameters->mf_box_size_x;
							kx-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrorx(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else if((simulation->parameters->boundary_x!=0)&& (simulation->parameters->boundary_y!=0)) //reflecting boundary in both directions
					{
						kx=-1;
						while (virtual_position_x>simulation->state->length_x)
						{
							virtual_position_x-=simulation->parameters->mf_box_size_x;
							kx+=1;
						}
						while (kx>0)
						{
							virtual_position_x-=simulation->parameters->mf_box_size_x;
							kx-=1;
						}
						ky=-1;
						while (virtual_position_y<0.0)
						{
							virtual_position_y+=simulation->parameters->mf_box_size_y;
							ky+=1;
						}
						while (ky>0)
						{
							virtual_position_y+=simulation->parameters->mf_box_size_y;
							ky-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrorxy(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
				}
				else if (virtual_position_x<0.0) //out of simulation box to the left
				{
					if (simulation->parameters->boundary_x==0)
					{
						virtual_position_x+=simulation->state->length_x;
						mfVM_sort_in_neighbors_in_box(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else
					{
						kx=-1;
						while (virtual_position_x<0.0)
						{
							virtual_position_x+=simulation->parameters->mf_box_size_x;
							kx+=1;
						}
						while (kx>0)
						{
							virtual_position_x+=simulation->parameters->mf_box_size_x;
							kx-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrorx(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
				}
				else if (virtual_position_x>simulation->state->length_x) //out of simulation box to the right
				{
					if (simulation->parameters->boundary_x==0)
					{
						virtual_position_x-=simulation->state->length_x;
						mfVM_sort_in_neighbors_in_box(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else
					{
						kx=-1;
						while (virtual_position_x>simulation->state->length_x)
						{
							virtual_position_x-=simulation->parameters->mf_box_size_x;
							kx+=1;
						}
						while (kx>0)
						{
							virtual_position_x-=simulation->parameters->mf_box_size_x;
							kx-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrorx(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
				}
				else if (virtual_position_y<0.0) //out of simulation box to the bottom
				{
					if (simulation->parameters->boundary_y==0)
					{
						virtual_position_y+=simulation->state->length_y;
						mfVM_sort_in_neighbors_in_box(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else
					{
						ky=-1;
						while (virtual_position_y<0.0)
						{
							virtual_position_y+=simulation->parameters->mf_box_size_y;
							ky+=1;
						}
						while (ky>0)
						{
							virtual_position_y+=simulation->parameters->mf_box_size_y;
							ky-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrory(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
				}
				else if (virtual_position_y>simulation->state->length_y) //out of simulation box to the top
				{
					if (simulation->parameters->boundary_y==0)
					{
						virtual_position_y-=simulation->state->length_y;
						mfVM_sort_in_neighbors_in_box(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else
					{
						ky=-1;
						while (virtual_position_y>simulation->state->length_y)
						{
							virtual_position_y-=simulation->parameters->mf_box_size_y;
							ky+=1;
						}
						while (ky>0)
						{
							virtual_position_y-=simulation->parameters->mf_box_size_y;
							ky-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrorx(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
				}
				else  // virtual position still in simulation box
				{
					mfVM_sort_in_neighbors_in_box(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
				}
			}
			j=order-1;
			for (i=2-order; i<order-1; i++)
			{
				virtual_position_x=(simulation->state->particles+index)->x+i*simulation->parameters->mf_box_size_x;
				virtual_position_y=(simulation->state->particles+index)->y+j*simulation->parameters->mf_box_size_y;
				if ((virtual_position_x<0.0) && (virtual_position_y<0.0)) //out of simulation box to left bottom
				{
					if ((simulation->parameters->boundary_x==0)&& (simulation->parameters->boundary_y==0)) //periodic boundary conditions in both directions
					{
						virtual_position_x+=simulation->state->length_x;
						virtual_position_y+=simulation->state->length_y;
						mfVM_sort_in_neighbors_in_box(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else if((simulation->parameters->boundary_x==0)&& (simulation->parameters->boundary_y!=0))//periodic in x-direction and reflecting in y-direction
					{
						virtual_position_x+=simulation->state->length_x;
						ky=-1;
						while (virtual_position_y<0.0)
						{
							virtual_position_y+=simulation->parameters->mf_box_size_y;
							ky+=1;
						}
						while (ky>0)
						{
							virtual_position_y+=simulation->parameters->mf_box_size_y;
							ky-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrory(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else if((simulation->parameters->boundary_x!=0)&& (simulation->parameters->boundary_y==0)) //reflecting in x-direction, periodic in y-direction
					{
						virtual_position_y+=simulation->state->length_y;
						kx=-1;
						while (virtual_position_x<0.0)
						{
							virtual_position_x+=simulation->parameters->mf_box_size_x;
							kx+=1;
						}
						while (kx>0)
						{
							virtual_position_x+=simulation->parameters->mf_box_size_x;
							kx-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrorx(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else if((simulation->parameters->boundary_x!=0)&& (simulation->parameters->boundary_y!=0)) //reflecting boundary in both directions
					{
						kx=-1;
						while (virtual_position_x<0.0)
						{
							virtual_position_x+=simulation->parameters->mf_box_size_x;
							kx+=1;
						}
						while (kx>0)
						{
							virtual_position_x+=simulation->parameters->mf_box_size_x;
							kx-=1;
						}
						ky=-1;
						while (virtual_position_y<0.0)
						{
							virtual_position_y+=simulation->parameters->mf_box_size_y;
							ky+=1;
						}
						while (ky>0)
						{
							virtual_position_y+=simulation->parameters->mf_box_size_y;
							ky-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrorxy(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
				}
				else if ((virtual_position_x<0.0) && (virtual_position_y>simulation->state->length_y)) //out of simulation box to left top
				{
					if ((simulation->parameters->boundary_x==0)&& (simulation->parameters->boundary_y==0)) //periodic boundary conditions in both directions
					{
						virtual_position_x+=simulation->state->length_x;
						virtual_position_y-=simulation->state->length_y;
						mfVM_sort_in_neighbors_in_box(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else if((simulation->parameters->boundary_x==0)&& (simulation->parameters->boundary_y!=0))//periodic in x-direction and reflecting in y-direction
					{
						virtual_position_x+=simulation->state->length_x;
						ky=-1;
						while (virtual_position_y>simulation->state->length_y)
						{
							virtual_position_y-=simulation->parameters->mf_box_size_y;
							ky+=1;
						}
						while (ky>0)
						{
							virtual_position_y-=simulation->parameters->mf_box_size_y;
							ky-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrory(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else if((simulation->parameters->boundary_x!=0)&& (simulation->parameters->boundary_y==0)) //reflecting in x-direction, periodic in y-direction
					{
						virtual_position_y-=simulation->state->length_y;
						kx=-1;
						while (virtual_position_x<0.0)
						{
							virtual_position_x+=simulation->parameters->mf_box_size_x;
							kx+=1;
						}
						while (kx>0)
						{
							virtual_position_x+=simulation->parameters->mf_box_size_x;
							kx-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrorx(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else if((simulation->parameters->boundary_x!=0)&& (simulation->parameters->boundary_y!=0)) //reflecting boundary in both directions
					{
						kx=-1;
						while (virtual_position_x<0.0)
						{
							virtual_position_x+=simulation->parameters->mf_box_size_x;
							kx+=1;
						}
						while (kx>0)
						{
							virtual_position_x+=simulation->parameters->mf_box_size_x;
							kx-=1;
						}
						ky=-1;
						while (virtual_position_y>simulation->state->length_y)
						{
							virtual_position_y-=simulation->parameters->mf_box_size_y;
							ky+=1;
						}
						while (ky>0)
						{
							virtual_position_y-=simulation->parameters->mf_box_size_y;
							ky-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrorxy(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
				}
				else if ((virtual_position_x>simulation->state->length_x) && (virtual_position_y>simulation->state->length_y)) //out of simulation box to right top
				{
					if ((simulation->parameters->boundary_x==0)&& (simulation->parameters->boundary_y==0)) //periodic boundary conditions in both directions
					{
						virtual_position_x-=simulation->state->length_x;
						virtual_position_y-=simulation->state->length_y;
						mfVM_sort_in_neighbors_in_box(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else if((simulation->parameters->boundary_x==0)&& (simulation->parameters->boundary_y!=0))//periodic in x-direction and reflecting in y-direction
					{
						virtual_position_x-=simulation->state->length_x;
						ky=-1;
						while (virtual_position_y>simulation->state->length_y)
						{
							virtual_position_y-=simulation->parameters->mf_box_size_y;
							ky+=1;
						}
						while (ky>0)
						{
							virtual_position_y-=simulation->parameters->mf_box_size_y;
							ky-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrory(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else if((simulation->parameters->boundary_x!=0)&& (simulation->parameters->boundary_y==0)) //reflecting in x-direction, periodic in y-direction
					{
						virtual_position_y-=simulation->state->length_y;
						kx=-1;
						while (virtual_position_x>simulation->state->length_x)
						{
							virtual_position_x-=simulation->parameters->mf_box_size_x;
							kx+=1;
						}
						while (kx>0)
						{
							virtual_position_x-=simulation->parameters->mf_box_size_x;
							kx-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrorx(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else if((simulation->parameters->boundary_x!=0)&& (simulation->parameters->boundary_y!=0)) //reflecting boundary in both directions
					{
						kx=-1;
						while (virtual_position_x>simulation->state->length_x)
						{
							virtual_position_x-=simulation->parameters->mf_box_size_x;
							kx+=1;
						}
						while (kx>0)
						{
							virtual_position_x-=simulation->parameters->mf_box_size_x;
							kx-=1;
						}
						ky=-1;
						while (virtual_position_y>simulation->state->length_y)
						{
							virtual_position_y-=simulation->parameters->mf_box_size_y;
							ky+=1;
						}
						while (ky>0)
						{
							virtual_position_y-=simulation->parameters->mf_box_size_y;
							ky-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrorxy(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
				}
				else if ((virtual_position_x>simulation->state->length_x) && (virtual_position_y<0.0)) //out of simulation box to right bottom
				{
					if ((simulation->parameters->boundary_x==0)&& (simulation->parameters->boundary_y==0)) //periodic boundary conditions in both directions
					{
						virtual_position_x-=simulation->state->length_x;
						virtual_position_y+=simulation->state->length_y;
						mfVM_sort_in_neighbors_in_box(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else if((simulation->parameters->boundary_x==0)&& (simulation->parameters->boundary_y!=0))//periodic in x-direction and reflecting in y-direction
					{
						virtual_position_x-=simulation->state->length_x;
						ky=-1;
						while (virtual_position_y<0.0)
						{
							virtual_position_y+=simulation->parameters->mf_box_size_y;
							ky+=1;
						}
						while (ky>0)
						{
							virtual_position_y+=simulation->parameters->mf_box_size_y;
							ky-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrory(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else if((simulation->parameters->boundary_x!=0)&& (simulation->parameters->boundary_y==0)) //reflecting in x-direction, periodic in y-direction
					{
						virtual_position_y+=simulation->state->length_y;
						kx=-1;
						while (virtual_position_x>simulation->state->length_x)
						{
							virtual_position_x-=simulation->parameters->mf_box_size_x;
							kx+=1;
						}
						while (kx>0)
						{
							virtual_position_x-=simulation->parameters->mf_box_size_x;
							kx-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrorx(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else if((simulation->parameters->boundary_x!=0)&& (simulation->parameters->boundary_y!=0)) //reflecting boundary in both directions
					{
						kx=-1;
						while (virtual_position_x>simulation->state->length_x)
						{
							virtual_position_x-=simulation->parameters->mf_box_size_x;
							kx+=1;
						}
						while (kx>0)
						{
							virtual_position_x-=simulation->parameters->mf_box_size_x;
							kx-=1;
						}
						ky=-1;
						while (virtual_position_y<0.0)
						{
							virtual_position_y+=simulation->parameters->mf_box_size_y;
							ky+=1;
						}
						while (ky>0)
						{
							virtual_position_y+=simulation->parameters->mf_box_size_y;
							ky-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrorxy(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
				}
				else if (virtual_position_x<0.0) //out of simulation box to the left
				{
					if (simulation->parameters->boundary_x==0)
					{
						virtual_position_x+=simulation->state->length_x;
						mfVM_sort_in_neighbors_in_box(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else
					{
						kx=-1;
						while (virtual_position_x<0.0)
						{
							virtual_position_x+=simulation->parameters->mf_box_size_x;
							kx+=1;
						}
						while (kx>0)
						{
							virtual_position_x+=simulation->parameters->mf_box_size_x;
							kx-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrorx(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
				}
				else if (virtual_position_x>simulation->state->length_x) //out of simulation box to the right
				{
					if (simulation->parameters->boundary_x==0)
					{
						virtual_position_x-=simulation->state->length_x;
						mfVM_sort_in_neighbors_in_box(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else
					{
						kx=-1;
						while (virtual_position_x>simulation->state->length_x)
						{
							virtual_position_x-=simulation->parameters->mf_box_size_x;
							kx+=1;
						}
						while (kx>0)
						{
							virtual_position_x-=simulation->parameters->mf_box_size_x;
							kx-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrorx(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
				}
				else if (virtual_position_y<0.0) //out of simulation box to the bottom
				{
					if (simulation->parameters->boundary_y==0)
					{
						virtual_position_y+=simulation->state->length_y;
						mfVM_sort_in_neighbors_in_box(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else
					{
						ky=-1;
						while (virtual_position_y<0.0)
						{
							virtual_position_y+=simulation->parameters->mf_box_size_y;
							ky+=1;
						}
						while (ky>0)
						{
							virtual_position_y+=simulation->parameters->mf_box_size_y;
							ky-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrory(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
				}
				else if (virtual_position_y>simulation->state->length_y) //out of simulation box to the top
				{
					if (simulation->parameters->boundary_y==0)
					{
						virtual_position_y-=simulation->state->length_y;
						mfVM_sort_in_neighbors_in_box(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else
					{
						ky=-1;
						while (virtual_position_y>simulation->state->length_y)
						{
							virtual_position_y-=simulation->parameters->mf_box_size_y;
							ky+=1;
						}
						while (ky>0)
						{
							virtual_position_y-=simulation->parameters->mf_box_size_y;
							ky-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrorx(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
				}
				else  // virtual position still in simulation box
				{
					mfVM_sort_in_neighbors_in_box(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
				}
			}
		}
		dx=position_x-((gint32)(position_x/simulation->parameters->mf_box_size_x))*simulation->parameters->mf_box_size_x;
		if (dx>simulation->parameters->mf_box_size_x*0.5)
			dx=simulation->parameters->mf_box_size_x-dx;
		dx=dx+(order-1)*simulation->parameters->mf_box_size_x;
		dy=position_y-((gint32)(position_y/simulation->parameters->mf_box_size_y))*simulation->parameters->mf_box_size_y;
		if (dy>simulation->parameters->mf_box_size_y*0.5)
			dy=simulation->parameters->mf_box_size_y-dy;
		dy=dy+(order-1)*simulation->parameters->mf_box_size_y;
		if (dx<dy)
			simulation->state->nh.max_distance_checked=dx;
		else
			simulation->state->nh.max_distance_checked=dy;
		order++;
	}
	//calculate orientational interaction with nearest neighbors
	for (i=0; i<simulation->parameters->mf_kn; i++)
	{
		//the neighboring particle is not 'behind' a reflecting boundary
		if (*(simulation->state->nh.mirror_flag+i)==0)
		{
			velocity_x+=cos((simulation->state->particles+*(simulation->state->nh.neighbor_indexes+i))->theta);
			velocity_y+=sin((simulation->state->particles+*(simulation->state->nh.neighbor_indexes+i))->theta);
		}
		//the neighboring particle is 'behind' a reflecting boundary in x direction -> the angle has to be reflected accordingly: theta->pi-theta
		else if (*(simulation->state->nh.mirror_flag+i)==1)
		{
			velocity_x+=cos(M_PI-(simulation->state->particles+*(simulation->state->nh.neighbor_indexes+i))->theta);
			velocity_y+=sin(M_PI-(simulation->state->particles+*(simulation->state->nh.neighbor_indexes+i))->theta);
		}
		//the neighboring particle is 'behind' a reflecting boundary in y direction -> the angle has to be reflected accordingly: theta->-theta
		else if (*(simulation->state->nh.mirror_flag+i)==2)
		{
			velocity_x+=cos(-(simulation->state->particles+*(simulation->state->nh.neighbor_indexes+i))->theta);
			velocity_y+=sin(-(simulation->state->particles+*(simulation->state->nh.neighbor_indexes+i))->theta);
		}
		//the neighboring particle is 'behind' a reflecting boundary in both, x- and y-direction -> the angle has to be reflected accordingly: theta->pi+theta
		else if (*(simulation->state->nh.mirror_flag+i)==3)
		{
			velocity_x+=cos(M_PI+(simulation->state->particles+*(simulation->state->nh.neighbor_indexes+i))->theta);
			velocity_y+=sin(M_PI+(simulation->state->particles+*(simulation->state->nh.neighbor_indexes+i))->theta);
		}
	}
	//calculate the new angle, considering the effect of interactions (calculated before), noise and chirality
	if (simulation->parameters->particle_species==1)
		(simulation->state->particles+index)->theta_new=atan2(velocity_y, velocity_x)+(g_rand_double(simulation->state->prng)*2*M_PI -M_PI)**simulation->parameters->noise_strength+*simulation->parameters->omega;
	else
		(simulation->state->particles+index)->theta_new=atan2(velocity_y, velocity_x)+(g_rand_double(simulation->state->prng)*2*M_PI -M_PI)**(simulation->parameters->noise_strength+get_species(index, simulation->parameters))+*(simulation->parameters->omega+get_species(index, simulation->parameters));
	while ((simulation->state->particles+index)->theta_new<-M_PI)
		(simulation->state->particles+index)->theta_new+=2.0*M_PI;
	while ((simulation->state->particles+index)->theta_new>M_PI)
		(simulation->state->particles+index)->theta_new-=2.0*M_PI;
	return simulation->parameters->mf_kn;
}

//calculates interactions of particle 'index' for metric free/topological nematic Vicsek model
gint32 mfNVM_interaction_with_neighbors(guint64 index, struct VM_simulation * simulation)
{
	mfVM_clean_neighborhood(simulation);
	gint32 order=1; // this says which boxes are checked: if order=1 only the inner box is checked, if order=2 the surrounding 8 boxes are checked and so on
	gint32 i, j;
	gint32 kx, ky;
	gdouble velocity_x=cos((simulation->state->particles+index)->theta);
	gdouble velocity_y=sin((simulation->state->particles+index)->theta);
	gdouble virtual_position_x;
	gdouble virtual_position_y;
	gdouble position_x=(simulation->state->particles+index)->x;
	gdouble position_y=(simulation->state->particles+index)->y;
	gdouble dx=0.0;
	gdouble dy=0.0;
	//find mf_kn nearest neighbors of particle 'index'
	while (simulation->state->nh.max_distance_checked<*(simulation->state->nh.distances+simulation->parameters->mf_kn-1))
	{
		if (order==1)
		{
			i=0;
			j=0;
			virtual_position_x=(simulation->state->particles+index)->x+i*simulation->parameters->mf_box_size_x;
			if (virtual_position_x<0.0)
				virtual_position_x+=simulation->state->length_x;
			if (virtual_position_x>simulation->state->length_x)
				virtual_position_x-=simulation->state->length_x;
			virtual_position_y=(simulation->state->particles+index)->y+j*simulation->parameters->mf_box_size_y;
			if (virtual_position_y<0.0)
				virtual_position_y+=simulation->state->length_y;
			if (virtual_position_y>simulation->state->length_y)
				virtual_position_y-=simulation->state->length_y;
			mfVM_sort_in_neighbors_in_box(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
		}
		else
		{
			i=1-order;
			for (j=1-order; j<order; j++)
			{
				virtual_position_x=(simulation->state->particles+index)->x+i*simulation->parameters->mf_box_size_x;
				virtual_position_y=(simulation->state->particles+index)->y+j*simulation->parameters->mf_box_size_y;
				if ((virtual_position_x<0.0) && (virtual_position_y<0.0)) //out of simulation box to left bottom
				{
					if ((simulation->parameters->boundary_x==0)&& (simulation->parameters->boundary_y==0)) //periodic boundary conditions in both directions
					{
						virtual_position_x+=simulation->state->length_x;
						virtual_position_y+=simulation->state->length_y;
						mfVM_sort_in_neighbors_in_box(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else if((simulation->parameters->boundary_x==0)&& (simulation->parameters->boundary_y!=0))//periodic in x-direction and reflecting in y-direction
					{
						virtual_position_x+=simulation->state->length_x;
						ky=-1;
						while (virtual_position_y<0.0)
						{
							virtual_position_y+=simulation->parameters->mf_box_size_y;
							ky+=1;
						}
						while (ky>0)
						{
							virtual_position_y+=simulation->parameters->mf_box_size_y;
							ky-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrory(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else if((simulation->parameters->boundary_x!=0)&& (simulation->parameters->boundary_y==0)) //reflecting in x-direction, periodic in y-direction
					{
						virtual_position_y+=simulation->state->length_y;
						kx=-1;
						while (virtual_position_x<0.0)
						{
							virtual_position_x+=simulation->parameters->mf_box_size_x;
							kx+=1;
						}
						while (kx>0)
						{
							virtual_position_x+=simulation->parameters->mf_box_size_x;
							kx-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrorx(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else if((simulation->parameters->boundary_x!=0)&& (simulation->parameters->boundary_y!=0)) //reflecting boundary in both directions
					{
						kx=-1;
						while (virtual_position_x<0.0)
						{
							virtual_position_x+=simulation->parameters->mf_box_size_x;
							kx+=1;
						}
						while (kx>0)
						{
							virtual_position_x+=simulation->parameters->mf_box_size_x;
							kx-=1;
						}
						ky=-1;
						while (virtual_position_y<0.0)
						{
							virtual_position_y+=simulation->parameters->mf_box_size_y;
							ky+=1;
						}
						while (ky>0)
						{
							virtual_position_y+=simulation->parameters->mf_box_size_y;
							ky-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrorxy(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
				}
				else if ((virtual_position_x<0.0) && (virtual_position_y>simulation->state->length_y)) //out of simulation box to left top
				{
					if ((simulation->parameters->boundary_x==0)&& (simulation->parameters->boundary_y==0)) //periodic boundary conditions in both directions
					{
						virtual_position_x+=simulation->state->length_x;
						virtual_position_y-=simulation->state->length_y;
						mfVM_sort_in_neighbors_in_box(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else if((simulation->parameters->boundary_x==0)&& (simulation->parameters->boundary_y!=0))//periodic in x-direction and reflecting in y-direction
					{
						virtual_position_x+=simulation->state->length_x;
						ky=-1;
						while (virtual_position_y>simulation->state->length_y)
						{
							virtual_position_y-=simulation->parameters->mf_box_size_y;
							ky+=1;
						}
						while (ky>0)
						{
							virtual_position_y-=simulation->parameters->mf_box_size_y;
							ky-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrory(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else if((simulation->parameters->boundary_x!=0)&& (simulation->parameters->boundary_y==0)) //reflecting in x-direction, periodic in y-direction
					{
						virtual_position_y-=simulation->state->length_y;
						kx=-1;
						while (virtual_position_x<0.0)
						{
							virtual_position_x+=simulation->parameters->mf_box_size_x;
							kx+=1;
						}
						while (kx>0)
						{
							virtual_position_x+=simulation->parameters->mf_box_size_x;
							kx-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrorx(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else if((simulation->parameters->boundary_x!=0)&& (simulation->parameters->boundary_y!=0)) //reflecting boundary in both directions
					{
						kx=-1;
						while (virtual_position_x<0.0)
						{
							virtual_position_x+=simulation->parameters->mf_box_size_x;
							kx+=1;
						}
						while (kx>0)
						{
							virtual_position_x+=simulation->parameters->mf_box_size_x;
							kx-=1;
						}
						ky=-1;
						while (virtual_position_y>simulation->state->length_y)
						{
							virtual_position_y-=simulation->parameters->mf_box_size_y;
							ky+=1;
						}
						while (ky>0)
						{
							virtual_position_y-=simulation->parameters->mf_box_size_y;
							ky-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrorxy(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
				}
				else if ((virtual_position_x>simulation->state->length_x) && (virtual_position_y>simulation->state->length_y)) //out of simulation box to right top
				{
					if ((simulation->parameters->boundary_x==0)&& (simulation->parameters->boundary_y==0)) //periodic boundary conditions in both directions
					{
						virtual_position_x-=simulation->state->length_x;
						virtual_position_y-=simulation->state->length_y;
						mfVM_sort_in_neighbors_in_box(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else if((simulation->parameters->boundary_x==0)&& (simulation->parameters->boundary_y!=0))//periodic in x-direction and reflecting in y-direction
					{
						virtual_position_x-=simulation->state->length_x;
						ky=-1;
						while (virtual_position_y>simulation->state->length_y)
						{
							virtual_position_y-=simulation->parameters->mf_box_size_y;
							ky+=1;
						}
						while (ky>0)
						{
							virtual_position_y-=simulation->parameters->mf_box_size_y;
							ky-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrory(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else if((simulation->parameters->boundary_x!=0)&& (simulation->parameters->boundary_y==0)) //reflecting in x-direction, periodic in y-direction
					{
						virtual_position_y-=simulation->state->length_y;
						kx=-1;
						while (virtual_position_x>simulation->state->length_x)
						{
							virtual_position_x-=simulation->parameters->mf_box_size_x;
							kx+=1;
						}
						while (kx>0)
						{
							virtual_position_x-=simulation->parameters->mf_box_size_x;
							kx-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrorx(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else if((simulation->parameters->boundary_x!=0)&& (simulation->parameters->boundary_y!=0)) //reflecting boundary in both directions
					{
						kx=-1;
						while (virtual_position_x>simulation->state->length_x)
						{
							virtual_position_x-=simulation->parameters->mf_box_size_x;
							kx+=1;
						}
						while (kx>0)
						{
							virtual_position_x-=simulation->parameters->mf_box_size_x;
							kx-=1;
						}
						ky=-1;
						while (virtual_position_y>simulation->state->length_y)
						{
							virtual_position_y-=simulation->parameters->mf_box_size_y;
							ky+=1;
						}
						while (ky>0)
						{
							virtual_position_y-=simulation->parameters->mf_box_size_y;
							ky-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrorxy(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
				}
				else if ((virtual_position_x>simulation->state->length_x) && (virtual_position_y<0.0)) //out of simulation box to right bottom
				{
					if ((simulation->parameters->boundary_x==0)&& (simulation->parameters->boundary_y==0)) //periodic boundary conditions in both directions
					{
						virtual_position_x-=simulation->state->length_x;
						virtual_position_y+=simulation->state->length_y;
						mfVM_sort_in_neighbors_in_box(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else if((simulation->parameters->boundary_x==0)&& (simulation->parameters->boundary_y!=0))//periodic in x-direction and reflecting in y-direction
					{
						virtual_position_x-=simulation->state->length_x;
						ky=-1;
						while (virtual_position_y<0.0)
						{
							virtual_position_y+=simulation->parameters->mf_box_size_y;
							ky+=1;
						}
						while (ky>0)
						{
							virtual_position_y+=simulation->parameters->mf_box_size_y;
							ky-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrory(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else if((simulation->parameters->boundary_x!=0)&& (simulation->parameters->boundary_y==0)) //reflecting in x-direction, periodic in y-direction
					{
						virtual_position_y+=simulation->state->length_y;
						kx=-1;
						while (virtual_position_x>simulation->state->length_x)
						{
							virtual_position_x-=simulation->parameters->mf_box_size_x;
							kx+=1;
						}
						while (kx>0)
						{
							virtual_position_x-=simulation->parameters->mf_box_size_x;
							kx-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrorx(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else if((simulation->parameters->boundary_x!=0)&& (simulation->parameters->boundary_y!=0)) //reflecting boundary in both directions
					{
						kx=-1;
						while (virtual_position_x>simulation->state->length_x)
						{
							virtual_position_x-=simulation->parameters->mf_box_size_x;
							kx+=1;
						}
						while (kx>0)
						{
							virtual_position_x-=simulation->parameters->mf_box_size_x;
							kx-=1;
						}
						ky=-1;
						while (virtual_position_y<0.0)
						{
							virtual_position_y+=simulation->parameters->mf_box_size_y;
							ky+=1;
						}
						while (ky>0)
						{
							virtual_position_y+=simulation->parameters->mf_box_size_y;
							ky-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrorxy(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
				}
				else if (virtual_position_x<0.0) //out of simulation box to the left
				{
					if (simulation->parameters->boundary_x==0)
					{
						virtual_position_x+=simulation->state->length_x;
						mfVM_sort_in_neighbors_in_box(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else
					{
						kx=-1;
						while (virtual_position_x<0.0)
						{
							virtual_position_x+=simulation->parameters->mf_box_size_x;
							kx+=1;
						}
						while (kx>0)
						{
							virtual_position_x+=simulation->parameters->mf_box_size_x;
							kx-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrorx(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
				}
				else if (virtual_position_x>simulation->state->length_x) //out of simulation box to the right
				{
					if (simulation->parameters->boundary_x==0)
					{
						virtual_position_x-=simulation->state->length_x;
						mfVM_sort_in_neighbors_in_box(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else
					{
						kx=-1;
						while (virtual_position_x>simulation->state->length_x)
						{
							virtual_position_x-=simulation->parameters->mf_box_size_x;
							kx+=1;
						}
						while (kx>0)
						{
							virtual_position_x-=simulation->parameters->mf_box_size_x;
							kx-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrorx(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
				}
				else if (virtual_position_y<0.0) //out of simulation box to the bottom
				{
					if (simulation->parameters->boundary_y==0)
					{
						virtual_position_y+=simulation->state->length_y;
						mfVM_sort_in_neighbors_in_box(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else
					{
						ky=-1;
						while (virtual_position_y<0.0)
						{
							virtual_position_y+=simulation->parameters->mf_box_size_y;
							ky+=1;
						}
						while (ky>0)
						{
							virtual_position_y+=simulation->parameters->mf_box_size_y;
							ky-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrory(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
				}
				else if (virtual_position_y>simulation->state->length_y) //out of simulation box to the top
				{
					if (simulation->parameters->boundary_y==0)
					{
						virtual_position_y-=simulation->state->length_y;
						mfVM_sort_in_neighbors_in_box(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else
					{
						ky=-1;
						while (virtual_position_y>simulation->state->length_y)
						{
							virtual_position_y-=simulation->parameters->mf_box_size_y;
							ky+=1;
						}
						while (ky>0)
						{
							virtual_position_y-=simulation->parameters->mf_box_size_y;
							ky-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrorx(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
				}
				else  // virtual position still in simulation box
				{
					mfVM_sort_in_neighbors_in_box(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
				}
			}
			i=order-1;
			for (j=1-order; j<order; j++)
			{
				virtual_position_x=(simulation->state->particles+index)->x+i*simulation->parameters->mf_box_size_x;
				virtual_position_y=(simulation->state->particles+index)->y+j*simulation->parameters->mf_box_size_y;
				if ((virtual_position_x<0.0) && (virtual_position_y<0.0)) //out of simulation box to left bottom
				{
					if ((simulation->parameters->boundary_x==0)&& (simulation->parameters->boundary_y==0)) //periodic boundary conditions in both directions
					{
						virtual_position_x+=simulation->state->length_x;
						virtual_position_y+=simulation->state->length_y;
						mfVM_sort_in_neighbors_in_box(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else if((simulation->parameters->boundary_x==0)&& (simulation->parameters->boundary_y!=0))//periodic in x-direction and reflecting in y-direction
					{
						virtual_position_x+=simulation->state->length_x;
						ky=-1;
						while (virtual_position_y<0.0)
						{
							virtual_position_y+=simulation->parameters->mf_box_size_y;
							ky+=1;
						}
						while (ky>0)
						{
							virtual_position_y+=simulation->parameters->mf_box_size_y;
							ky-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrory(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else if((simulation->parameters->boundary_x!=0)&& (simulation->parameters->boundary_y==0)) //reflecting in x-direction, periodic in y-direction
					{
						virtual_position_y+=simulation->state->length_y;
						kx=-1;
						while (virtual_position_x<0.0)
						{
							virtual_position_x+=simulation->parameters->mf_box_size_x;
							kx+=1;
						}
						while (kx>0)
						{
							virtual_position_x+=simulation->parameters->mf_box_size_x;
							kx-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrorx(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else if((simulation->parameters->boundary_x!=0)&& (simulation->parameters->boundary_y!=0)) //reflecting boundary in both directions
					{
						kx=-1;
						while (virtual_position_x<0.0)
						{
							virtual_position_x+=simulation->parameters->mf_box_size_x;
							kx+=1;
						}
						while (kx>0)
						{
							virtual_position_x+=simulation->parameters->mf_box_size_x;
							kx-=1;
						}
						ky=-1;
						while (virtual_position_y<0.0)
						{
							virtual_position_y+=simulation->parameters->mf_box_size_y;
							ky+=1;
						}
						while (ky>0)
						{
							virtual_position_y+=simulation->parameters->mf_box_size_y;
							ky-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrorxy(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
				}
				else if ((virtual_position_x<0.0) && (virtual_position_y>simulation->state->length_y)) //out of simulation box to left top
				{
					if ((simulation->parameters->boundary_x==0)&& (simulation->parameters->boundary_y==0)) //periodic boundary conditions in both directions
					{
						virtual_position_x+=simulation->state->length_x;
						virtual_position_y-=simulation->state->length_y;
						mfVM_sort_in_neighbors_in_box(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else if((simulation->parameters->boundary_x==0)&& (simulation->parameters->boundary_y!=0))//periodic in x-direction and reflecting in y-direction
					{
						virtual_position_x+=simulation->state->length_x;
						ky=-1;
						while (virtual_position_y>simulation->state->length_y)
						{
							virtual_position_y-=simulation->parameters->mf_box_size_y;
							ky+=1;
						}
						while (ky>0)
						{
							virtual_position_y-=simulation->parameters->mf_box_size_y;
							ky-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrory(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else if((simulation->parameters->boundary_x!=0)&& (simulation->parameters->boundary_y==0)) //reflecting in x-direction, periodic in y-direction
					{
						virtual_position_y-=simulation->state->length_y;
						kx=-1;
						while (virtual_position_x<0.0)
						{
							virtual_position_x+=simulation->parameters->mf_box_size_x;
							kx+=1;
						}
						while (kx>0)
						{
							virtual_position_x+=simulation->parameters->mf_box_size_x;
							kx-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrorx(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else if((simulation->parameters->boundary_x!=0)&& (simulation->parameters->boundary_y!=0)) //reflecting boundary in both directions
					{
						kx=-1;
						while (virtual_position_x<0.0)
						{
							virtual_position_x+=simulation->parameters->mf_box_size_x;
							kx+=1;
						}
						while (kx>0)
						{
							virtual_position_x+=simulation->parameters->mf_box_size_x;
							kx-=1;
						}
						ky=-1;
						while (virtual_position_y>simulation->state->length_y)
						{
							virtual_position_y-=simulation->parameters->mf_box_size_y;
							ky+=1;
						}
						while (ky>0)
						{
							virtual_position_y-=simulation->parameters->mf_box_size_y;
							ky-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrorxy(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
				}
				else if ((virtual_position_x>simulation->state->length_x) && (virtual_position_y>simulation->state->length_y)) //out of simulation box to right top
				{
					if ((simulation->parameters->boundary_x==0)&& (simulation->parameters->boundary_y==0)) //periodic boundary conditions in both directions
					{
						virtual_position_x-=simulation->state->length_x;
						virtual_position_y-=simulation->state->length_y;
						mfVM_sort_in_neighbors_in_box(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else if((simulation->parameters->boundary_x==0)&& (simulation->parameters->boundary_y!=0))//periodic in x-direction and reflecting in y-direction
					{
						virtual_position_x-=simulation->state->length_x;
						ky=-1;
						while (virtual_position_y>simulation->state->length_y)
						{
							virtual_position_y-=simulation->parameters->mf_box_size_y;
							ky+=1;
						}
						while (ky>0)
						{
							virtual_position_y-=simulation->parameters->mf_box_size_y;
							ky-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrory(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else if((simulation->parameters->boundary_x!=0)&& (simulation->parameters->boundary_y==0)) //reflecting in x-direction, periodic in y-direction
					{
						virtual_position_y-=simulation->state->length_y;
						kx=-1;
						while (virtual_position_x>simulation->state->length_x)
						{
							virtual_position_x-=simulation->parameters->mf_box_size_x;
							kx+=1;
						}
						while (kx>0)
						{
							virtual_position_x-=simulation->parameters->mf_box_size_x;
							kx-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrorx(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else if((simulation->parameters->boundary_x!=0)&& (simulation->parameters->boundary_y!=0)) //reflecting boundary in both directions
					{
						kx=-1;
						while (virtual_position_x>simulation->state->length_x)
						{
							virtual_position_x-=simulation->parameters->mf_box_size_x;
							kx+=1;
						}
						while (kx>0)
						{
							virtual_position_x-=simulation->parameters->mf_box_size_x;
							kx-=1;
						}
						ky=-1;
						while (virtual_position_y>simulation->state->length_y)
						{
							virtual_position_y-=simulation->parameters->mf_box_size_y;
							ky+=1;
						}
						while (ky>0)
						{
							virtual_position_y-=simulation->parameters->mf_box_size_y;
							ky-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrorxy(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
				}
				else if ((virtual_position_x>simulation->state->length_x) && (virtual_position_y<0.0)) //out of simulation box to right bottom
				{
					if ((simulation->parameters->boundary_x==0)&& (simulation->parameters->boundary_y==0)) //periodic boundary conditions in both directions
					{
						virtual_position_x-=simulation->state->length_x;
						virtual_position_y+=simulation->state->length_y;
						mfVM_sort_in_neighbors_in_box(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else if((simulation->parameters->boundary_x==0)&& (simulation->parameters->boundary_y!=0))//periodic in x-direction and reflecting in y-direction
					{
						virtual_position_x-=simulation->state->length_x;
						ky=-1;
						while (virtual_position_y<0.0)
						{
							virtual_position_y+=simulation->parameters->mf_box_size_y;
							ky+=1;
						}
						while (ky>0)
						{
							virtual_position_y+=simulation->parameters->mf_box_size_y;
							ky-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrory(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else if((simulation->parameters->boundary_x!=0)&& (simulation->parameters->boundary_y==0)) //reflecting in x-direction, periodic in y-direction
					{
						virtual_position_y+=simulation->state->length_y;
						kx=-1;
						while (virtual_position_x>simulation->state->length_x)
						{
							virtual_position_x-=simulation->parameters->mf_box_size_x;
							kx+=1;
						}
						while (kx>0)
						{
							virtual_position_x-=simulation->parameters->mf_box_size_x;
							kx-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrorx(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else if((simulation->parameters->boundary_x!=0)&& (simulation->parameters->boundary_y!=0)) //reflecting boundary in both directions
					{
						kx=-1;
						while (virtual_position_x>simulation->state->length_x)
						{
							virtual_position_x-=simulation->parameters->mf_box_size_x;
							kx+=1;
						}
						while (kx>0)
						{
							virtual_position_x-=simulation->parameters->mf_box_size_x;
							kx-=1;
						}
						ky=-1;
						while (virtual_position_y<0.0)
						{
							virtual_position_y+=simulation->parameters->mf_box_size_y;
							ky+=1;
						}
						while (ky>0)
						{
							virtual_position_y+=simulation->parameters->mf_box_size_y;
							ky-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrorxy(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
				}
				else if (virtual_position_x<0.0) //out of simulation box to the left
				{
					if (simulation->parameters->boundary_x==0)
					{
						virtual_position_x+=simulation->state->length_x;
						mfVM_sort_in_neighbors_in_box(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else
					{
						kx=-1;
						while (virtual_position_x<0.0)
						{
							virtual_position_x+=simulation->parameters->mf_box_size_x;
							kx+=1;
						}
						while (kx>0)
						{
							virtual_position_x+=simulation->parameters->mf_box_size_x;
							kx-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrorx(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
				}
				else if (virtual_position_x>simulation->state->length_x) //out of simulation box to the right
				{
					if (simulation->parameters->boundary_x==0)
					{
						virtual_position_x-=simulation->state->length_x;
						mfVM_sort_in_neighbors_in_box(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else
					{
						kx=-1;
						while (virtual_position_x>simulation->state->length_x)
						{
							virtual_position_x-=simulation->parameters->mf_box_size_x;
							kx+=1;
						}
						while (kx>0)
						{
							virtual_position_x-=simulation->parameters->mf_box_size_x;
							kx-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrorx(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
				}
				else if (virtual_position_y<0.0) //out of simulation box to the bottom
				{
					if (simulation->parameters->boundary_y==0)
					{
						virtual_position_y+=simulation->state->length_y;
						mfVM_sort_in_neighbors_in_box(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else
					{
						ky=-1;
						while (virtual_position_y<0.0)
						{
							virtual_position_y+=simulation->parameters->mf_box_size_y;
							ky+=1;
						}
						while (ky>0)
						{
							virtual_position_y+=simulation->parameters->mf_box_size_y;
							ky-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrory(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
				}
				else if (virtual_position_y>simulation->state->length_y) //out of simulation box to the top
				{
					if (simulation->parameters->boundary_y==0)
					{
						virtual_position_y-=simulation->state->length_y;
						mfVM_sort_in_neighbors_in_box(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else
					{
						ky=-1;
						while (virtual_position_y>simulation->state->length_y)
						{
							virtual_position_y-=simulation->parameters->mf_box_size_y;
							ky+=1;
						}
						while (ky>0)
						{
							virtual_position_y-=simulation->parameters->mf_box_size_y;
							ky-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrorx(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
				}
				else  // virtual position still in simulation box
				{
					mfVM_sort_in_neighbors_in_box(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
				}
			}
			j=1-order;
			for (i=2-order; i<order-1; i++)
			{
				virtual_position_x=(simulation->state->particles+index)->x+i*simulation->parameters->mf_box_size_x;
				virtual_position_y=(simulation->state->particles+index)->y+j*simulation->parameters->mf_box_size_y;
				if ((virtual_position_x<0.0) && (virtual_position_y<0.0)) //out of simulation box to left bottom
				{
					if ((simulation->parameters->boundary_x==0)&& (simulation->parameters->boundary_y==0)) //periodic boundary conditions in both directions
					{
						virtual_position_x+=simulation->state->length_x;
						virtual_position_y+=simulation->state->length_y;
						mfVM_sort_in_neighbors_in_box(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else if((simulation->parameters->boundary_x==0)&& (simulation->parameters->boundary_y!=0))//periodic in x-direction and reflecting in y-direction
					{
						virtual_position_x+=simulation->state->length_x;
						ky=-1;
						while (virtual_position_y<0.0)
						{
							virtual_position_y+=simulation->parameters->mf_box_size_y;
							ky+=1;
						}
						while (ky>0)
						{
							virtual_position_y+=simulation->parameters->mf_box_size_y;
							ky-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrory(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else if((simulation->parameters->boundary_x!=0)&& (simulation->parameters->boundary_y==0)) //reflecting in x-direction, periodic in y-direction
					{
						virtual_position_y+=simulation->state->length_y;
						kx=-1;
						while (virtual_position_x<0.0)
						{
							virtual_position_x+=simulation->parameters->mf_box_size_x;
							kx+=1;
						}
						while (kx>0)
						{
							virtual_position_x+=simulation->parameters->mf_box_size_x;
							kx-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrorx(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else if((simulation->parameters->boundary_x!=0)&& (simulation->parameters->boundary_y!=0)) //reflecting boundary in both directions
					{
						kx=-1;
						while (virtual_position_x<0.0)
						{
							virtual_position_x+=simulation->parameters->mf_box_size_x;
							kx+=1;
						}
						while (kx>0)
						{
							virtual_position_x+=simulation->parameters->mf_box_size_x;
							kx-=1;
						}
						ky=-1;
						while (virtual_position_y<0.0)
						{
							virtual_position_y+=simulation->parameters->mf_box_size_y;
							ky+=1;
						}
						while (ky>0)
						{
							virtual_position_y+=simulation->parameters->mf_box_size_y;
							ky-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrorxy(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
				}
				else if ((virtual_position_x<0.0) && (virtual_position_y>simulation->state->length_y)) //out of simulation box to left top
				{
					if ((simulation->parameters->boundary_x==0)&& (simulation->parameters->boundary_y==0)) //periodic boundary conditions in both directions
					{
						virtual_position_x+=simulation->state->length_x;
						virtual_position_y-=simulation->state->length_y;
						mfVM_sort_in_neighbors_in_box(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else if((simulation->parameters->boundary_x==0)&& (simulation->parameters->boundary_y!=0))//periodic in x-direction and reflecting in y-direction
					{
						virtual_position_x+=simulation->state->length_x;
						ky=-1;
						while (virtual_position_y>simulation->state->length_y)
						{
							virtual_position_y-=simulation->parameters->mf_box_size_y;
							ky+=1;
						}
						while (ky>0)
						{
							virtual_position_y-=simulation->parameters->mf_box_size_y;
							ky-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrory(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else if((simulation->parameters->boundary_x!=0)&& (simulation->parameters->boundary_y==0)) //reflecting in x-direction, periodic in y-direction
					{
						virtual_position_y-=simulation->state->length_y;
						kx=-1;
						while (virtual_position_x<0.0)
						{
							virtual_position_x+=simulation->parameters->mf_box_size_x;
							kx+=1;
						}
						while (kx>0)
						{
							virtual_position_x+=simulation->parameters->mf_box_size_x;
							kx-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrorx(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else if((simulation->parameters->boundary_x!=0)&& (simulation->parameters->boundary_y!=0)) //reflecting boundary in both directions
					{
						kx=-1;
						while (virtual_position_x<0.0)
						{
							virtual_position_x+=simulation->parameters->mf_box_size_x;
							kx+=1;
						}
						while (kx>0)
						{
							virtual_position_x+=simulation->parameters->mf_box_size_x;
							kx-=1;
						}
						ky=-1;
						while (virtual_position_y>simulation->state->length_y)
						{
							virtual_position_y-=simulation->parameters->mf_box_size_y;
							ky+=1;
						}
						while (ky>0)
						{
							virtual_position_y-=simulation->parameters->mf_box_size_y;
							ky-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrorxy(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
				}
				else if ((virtual_position_x>simulation->state->length_x) && (virtual_position_y>simulation->state->length_y)) //out of simulation box to right top
				{
					if ((simulation->parameters->boundary_x==0)&& (simulation->parameters->boundary_y==0)) //periodic boundary conditions in both directions
					{
						virtual_position_x-=simulation->state->length_x;
						virtual_position_y-=simulation->state->length_y;
						mfVM_sort_in_neighbors_in_box(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else if((simulation->parameters->boundary_x==0)&& (simulation->parameters->boundary_y!=0))//periodic in x-direction and reflecting in y-direction
					{
						virtual_position_x-=simulation->state->length_x;
						ky=-1;
						while (virtual_position_y>simulation->state->length_y)
						{
							virtual_position_y-=simulation->parameters->mf_box_size_y;
							ky+=1;
						}
						while (ky>0)
						{
							virtual_position_y-=simulation->parameters->mf_box_size_y;
							ky-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrory(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else if((simulation->parameters->boundary_x!=0)&& (simulation->parameters->boundary_y==0)) //reflecting in x-direction, periodic in y-direction
					{
						virtual_position_y-=simulation->state->length_y;
						kx=-1;
						while (virtual_position_x>simulation->state->length_x)
						{
							virtual_position_x-=simulation->parameters->mf_box_size_x;
							kx+=1;
						}
						while (kx>0)
						{
							virtual_position_x-=simulation->parameters->mf_box_size_x;
							kx-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrorx(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else if((simulation->parameters->boundary_x!=0)&& (simulation->parameters->boundary_y!=0)) //reflecting boundary in both directions
					{
						kx=-1;
						while (virtual_position_x>simulation->state->length_x)
						{
							virtual_position_x-=simulation->parameters->mf_box_size_x;
							kx+=1;
						}
						while (kx>0)
						{
							virtual_position_x-=simulation->parameters->mf_box_size_x;
							kx-=1;
						}
						ky=-1;
						while (virtual_position_y>simulation->state->length_y)
						{
							virtual_position_y-=simulation->parameters->mf_box_size_y;
							ky+=1;
						}
						while (ky>0)
						{
							virtual_position_y-=simulation->parameters->mf_box_size_y;
							ky-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrorxy(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
				}
				else if ((virtual_position_x>simulation->state->length_x) && (virtual_position_y<0.0)) //out of simulation box to right bottom
				{
					if ((simulation->parameters->boundary_x==0)&& (simulation->parameters->boundary_y==0)) //periodic boundary conditions in both directions
					{
						virtual_position_x-=simulation->state->length_x;
						virtual_position_y+=simulation->state->length_y;
						mfVM_sort_in_neighbors_in_box(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else if((simulation->parameters->boundary_x==0)&& (simulation->parameters->boundary_y!=0))//periodic in x-direction and reflecting in y-direction
					{
						virtual_position_x-=simulation->state->length_x;
						ky=-1;
						while (virtual_position_y<0.0)
						{
							virtual_position_y+=simulation->parameters->mf_box_size_y;
							ky+=1;
						}
						while (ky>0)
						{
							virtual_position_y+=simulation->parameters->mf_box_size_y;
							ky-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrory(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else if((simulation->parameters->boundary_x!=0)&& (simulation->parameters->boundary_y==0)) //reflecting in x-direction, periodic in y-direction
					{
						virtual_position_y+=simulation->state->length_y;
						kx=-1;
						while (virtual_position_x>simulation->state->length_x)
						{
							virtual_position_x-=simulation->parameters->mf_box_size_x;
							kx+=1;
						}
						while (kx>0)
						{
							virtual_position_x-=simulation->parameters->mf_box_size_x;
							kx-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrorx(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else if((simulation->parameters->boundary_x!=0)&& (simulation->parameters->boundary_y!=0)) //reflecting boundary in both directions
					{
						kx=-1;
						while (virtual_position_x>simulation->state->length_x)
						{
							virtual_position_x-=simulation->parameters->mf_box_size_x;
							kx+=1;
						}
						while (kx>0)
						{
							virtual_position_x-=simulation->parameters->mf_box_size_x;
							kx-=1;
						}
						ky=-1;
						while (virtual_position_y<0.0)
						{
							virtual_position_y+=simulation->parameters->mf_box_size_y;
							ky+=1;
						}
						while (ky>0)
						{
							virtual_position_y+=simulation->parameters->mf_box_size_y;
							ky-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrorxy(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
				}
				else if (virtual_position_x<0.0) //out of simulation box to the left
				{
					if (simulation->parameters->boundary_x==0)
					{
						virtual_position_x+=simulation->state->length_x;
						mfVM_sort_in_neighbors_in_box(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else
					{
						kx=-1;
						while (virtual_position_x<0.0)
						{
							virtual_position_x+=simulation->parameters->mf_box_size_x;
							kx+=1;
						}
						while (kx>0)
						{
							virtual_position_x+=simulation->parameters->mf_box_size_x;
							kx-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrorx(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
				}
				else if (virtual_position_x>simulation->state->length_x) //out of simulation box to the right
				{
					if (simulation->parameters->boundary_x==0)
					{
						virtual_position_x-=simulation->state->length_x;
						mfVM_sort_in_neighbors_in_box(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else
					{
						kx=-1;
						while (virtual_position_x>simulation->state->length_x)
						{
							virtual_position_x-=simulation->parameters->mf_box_size_x;
							kx+=1;
						}
						while (kx>0)
						{
							virtual_position_x-=simulation->parameters->mf_box_size_x;
							kx-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrorx(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
				}
				else if (virtual_position_y<0.0) //out of simulation box to the bottom
				{
					if (simulation->parameters->boundary_y==0)
					{
						virtual_position_y+=simulation->state->length_y;
						mfVM_sort_in_neighbors_in_box(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else
					{
						ky=-1;
						while (virtual_position_y<0.0)
						{
							virtual_position_y+=simulation->parameters->mf_box_size_y;
							ky+=1;
						}
						while (ky>0)
						{
							virtual_position_y+=simulation->parameters->mf_box_size_y;
							ky-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrory(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
				}
				else if (virtual_position_y>simulation->state->length_y) //out of simulation box to the top
				{
					if (simulation->parameters->boundary_y==0)
					{
						virtual_position_y-=simulation->state->length_y;
						mfVM_sort_in_neighbors_in_box(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else
					{
						ky=-1;
						while (virtual_position_y>simulation->state->length_y)
						{
							virtual_position_y-=simulation->parameters->mf_box_size_y;
							ky+=1;
						}
						while (ky>0)
						{
							virtual_position_y-=simulation->parameters->mf_box_size_y;
							ky-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrorx(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
				}
				else  // virtual position still in simulation box
				{
					mfVM_sort_in_neighbors_in_box(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
				}
			}
			j=order-1;
			for (i=2-order; i<order-1; i++)
			{
				virtual_position_x=(simulation->state->particles+index)->x+i*simulation->parameters->mf_box_size_x;
				virtual_position_y=(simulation->state->particles+index)->y+j*simulation->parameters->mf_box_size_y;
				if ((virtual_position_x<0.0) && (virtual_position_y<0.0)) //out of simulation box to left bottom
				{
					if ((simulation->parameters->boundary_x==0)&& (simulation->parameters->boundary_y==0)) //periodic boundary conditions in both directions
					{
						virtual_position_x+=simulation->state->length_x;
						virtual_position_y+=simulation->state->length_y;
						mfVM_sort_in_neighbors_in_box(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else if((simulation->parameters->boundary_x==0)&& (simulation->parameters->boundary_y!=0))//periodic in x-direction and reflecting in y-direction
					{
						virtual_position_x+=simulation->state->length_x;
						ky=-1;
						while (virtual_position_y<0.0)
						{
							virtual_position_y+=simulation->parameters->mf_box_size_y;
							ky+=1;
						}
						while (ky>0)
						{
							virtual_position_y+=simulation->parameters->mf_box_size_y;
							ky-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrory(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else if((simulation->parameters->boundary_x!=0)&& (simulation->parameters->boundary_y==0)) //reflecting in x-direction, periodic in y-direction
					{
						virtual_position_y+=simulation->state->length_y;
						kx=-1;
						while (virtual_position_x<0.0)
						{
							virtual_position_x+=simulation->parameters->mf_box_size_x;
							kx+=1;
						}
						while (kx>0)
						{
							virtual_position_x+=simulation->parameters->mf_box_size_x;
							kx-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrorx(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else if((simulation->parameters->boundary_x!=0)&& (simulation->parameters->boundary_y!=0)) //reflecting boundary in both directions
					{
						kx=-1;
						while (virtual_position_x<0.0)
						{
							virtual_position_x+=simulation->parameters->mf_box_size_x;
							kx+=1;
						}
						while (kx>0)
						{
							virtual_position_x+=simulation->parameters->mf_box_size_x;
							kx-=1;
						}
						ky=-1;
						while (virtual_position_y<0.0)
						{
							virtual_position_y+=simulation->parameters->mf_box_size_y;
							ky+=1;
						}
						while (ky>0)
						{
							virtual_position_y+=simulation->parameters->mf_box_size_y;
							ky-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrorxy(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
				}
				else if ((virtual_position_x<0.0) && (virtual_position_y>simulation->state->length_y)) //out of simulation box to left top
				{
					if ((simulation->parameters->boundary_x==0)&& (simulation->parameters->boundary_y==0)) //periodic boundary conditions in both directions
					{
						virtual_position_x+=simulation->state->length_x;
						virtual_position_y-=simulation->state->length_y;
						mfVM_sort_in_neighbors_in_box(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else if((simulation->parameters->boundary_x==0)&& (simulation->parameters->boundary_y!=0))//periodic in x-direction and reflecting in y-direction
					{
						virtual_position_x+=simulation->state->length_x;
						ky=-1;
						while (virtual_position_y>simulation->state->length_y)
						{
							virtual_position_y-=simulation->parameters->mf_box_size_y;
							ky+=1;
						}
						while (ky>0)
						{
							virtual_position_y-=simulation->parameters->mf_box_size_y;
							ky-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrory(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else if((simulation->parameters->boundary_x!=0)&& (simulation->parameters->boundary_y==0)) //reflecting in x-direction, periodic in y-direction
					{
						virtual_position_y-=simulation->state->length_y;
						kx=-1;
						while (virtual_position_x<0.0)
						{
							virtual_position_x+=simulation->parameters->mf_box_size_x;
							kx+=1;
						}
						while (kx>0)
						{
							virtual_position_x+=simulation->parameters->mf_box_size_x;
							kx-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrorx(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else if((simulation->parameters->boundary_x!=0)&& (simulation->parameters->boundary_y!=0)) //reflecting boundary in both directions
					{
						kx=-1;
						while (virtual_position_x<0.0)
						{
							virtual_position_x+=simulation->parameters->mf_box_size_x;
							kx+=1;
						}
						while (kx>0)
						{
							virtual_position_x+=simulation->parameters->mf_box_size_x;
							kx-=1;
						}
						ky=-1;
						while (virtual_position_y>simulation->state->length_y)
						{
							virtual_position_y-=simulation->parameters->mf_box_size_y;
							ky+=1;
						}
						while (ky>0)
						{
							virtual_position_y-=simulation->parameters->mf_box_size_y;
							ky-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrorxy(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
				}
				else if ((virtual_position_x>simulation->state->length_x) && (virtual_position_y>simulation->state->length_y)) //out of simulation box to right top
				{
					if ((simulation->parameters->boundary_x==0)&& (simulation->parameters->boundary_y==0)) //periodic boundary conditions in both directions
					{
						virtual_position_x-=simulation->state->length_x;
						virtual_position_y-=simulation->state->length_y;
						mfVM_sort_in_neighbors_in_box(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else if((simulation->parameters->boundary_x==0)&& (simulation->parameters->boundary_y!=0))//periodic in x-direction and reflecting in y-direction
					{
						virtual_position_x-=simulation->state->length_x;
						ky=-1;
						while (virtual_position_y>simulation->state->length_y)
						{
							virtual_position_y-=simulation->parameters->mf_box_size_y;
							ky+=1;
						}
						while (ky>0)
						{
							virtual_position_y-=simulation->parameters->mf_box_size_y;
							ky-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrory(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else if((simulation->parameters->boundary_x!=0)&& (simulation->parameters->boundary_y==0)) //reflecting in x-direction, periodic in y-direction
					{
						virtual_position_y-=simulation->state->length_y;
						kx=-1;
						while (virtual_position_x>simulation->state->length_x)
						{
							virtual_position_x-=simulation->parameters->mf_box_size_x;
							kx+=1;
						}
						while (kx>0)
						{
							virtual_position_x-=simulation->parameters->mf_box_size_x;
							kx-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrorx(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else if((simulation->parameters->boundary_x!=0)&& (simulation->parameters->boundary_y!=0)) //reflecting boundary in both directions
					{
						kx=-1;
						while (virtual_position_x>simulation->state->length_x)
						{
							virtual_position_x-=simulation->parameters->mf_box_size_x;
							kx+=1;
						}
						while (kx>0)
						{
							virtual_position_x-=simulation->parameters->mf_box_size_x;
							kx-=1;
						}
						ky=-1;
						while (virtual_position_y>simulation->state->length_y)
						{
							virtual_position_y-=simulation->parameters->mf_box_size_y;
							ky+=1;
						}
						while (ky>0)
						{
							virtual_position_y-=simulation->parameters->mf_box_size_y;
							ky-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrorxy(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
				}
				else if ((virtual_position_x>simulation->state->length_x) && (virtual_position_y<0.0)) //out of simulation box to right bottom
				{
					if ((simulation->parameters->boundary_x==0)&& (simulation->parameters->boundary_y==0)) //periodic boundary conditions in both directions
					{
						virtual_position_x-=simulation->state->length_x;
						virtual_position_y+=simulation->state->length_y;
						mfVM_sort_in_neighbors_in_box(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else if((simulation->parameters->boundary_x==0)&& (simulation->parameters->boundary_y!=0))//periodic in x-direction and reflecting in y-direction
					{
						virtual_position_x-=simulation->state->length_x;
						ky=-1;
						while (virtual_position_y<0.0)
						{
							virtual_position_y+=simulation->parameters->mf_box_size_y;
							ky+=1;
						}
						while (ky>0)
						{
							virtual_position_y+=simulation->parameters->mf_box_size_y;
							ky-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrory(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else if((simulation->parameters->boundary_x!=0)&& (simulation->parameters->boundary_y==0)) //reflecting in x-direction, periodic in y-direction
					{
						virtual_position_y+=simulation->state->length_y;
						kx=-1;
						while (virtual_position_x>simulation->state->length_x)
						{
							virtual_position_x-=simulation->parameters->mf_box_size_x;
							kx+=1;
						}
						while (kx>0)
						{
							virtual_position_x-=simulation->parameters->mf_box_size_x;
							kx-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrorx(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else if((simulation->parameters->boundary_x!=0)&& (simulation->parameters->boundary_y!=0)) //reflecting boundary in both directions
					{
						kx=-1;
						while (virtual_position_x>simulation->state->length_x)
						{
							virtual_position_x-=simulation->parameters->mf_box_size_x;
							kx+=1;
						}
						while (kx>0)
						{
							virtual_position_x-=simulation->parameters->mf_box_size_x;
							kx-=1;
						}
						ky=-1;
						while (virtual_position_y<0.0)
						{
							virtual_position_y+=simulation->parameters->mf_box_size_y;
							ky+=1;
						}
						while (ky>0)
						{
							virtual_position_y+=simulation->parameters->mf_box_size_y;
							ky-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrorxy(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
				}
				else if (virtual_position_x<0.0) //out of simulation box to the left
				{
					if (simulation->parameters->boundary_x==0)
					{
						virtual_position_x+=simulation->state->length_x;
						mfVM_sort_in_neighbors_in_box(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else
					{
						kx=-1;
						while (virtual_position_x<0.0)
						{
							virtual_position_x+=simulation->parameters->mf_box_size_x;
							kx+=1;
						}
						while (kx>0)
						{
							virtual_position_x+=simulation->parameters->mf_box_size_x;
							kx-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrorx(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
				}
				else if (virtual_position_x>simulation->state->length_x) //out of simulation box to the right
				{
					if (simulation->parameters->boundary_x==0)
					{
						virtual_position_x-=simulation->state->length_x;
						mfVM_sort_in_neighbors_in_box(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else
					{
						kx=-1;
						while (virtual_position_x>simulation->state->length_x)
						{
							virtual_position_x-=simulation->parameters->mf_box_size_x;
							kx+=1;
						}
						while (kx>0)
						{
							virtual_position_x-=simulation->parameters->mf_box_size_x;
							kx-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrorx(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
				}
				else if (virtual_position_y<0.0) //out of simulation box to the bottom
				{
					if (simulation->parameters->boundary_y==0)
					{
						virtual_position_y+=simulation->state->length_y;
						mfVM_sort_in_neighbors_in_box(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else
					{
						ky=-1;
						while (virtual_position_y<0.0)
						{
							virtual_position_y+=simulation->parameters->mf_box_size_y;
							ky+=1;
						}
						while (ky>0)
						{
							virtual_position_y+=simulation->parameters->mf_box_size_y;
							ky-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrory(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
				}
				else if (virtual_position_y>simulation->state->length_y) //out of simulation box to the top
				{
					if (simulation->parameters->boundary_y==0)
					{
						virtual_position_y-=simulation->state->length_y;
						mfVM_sort_in_neighbors_in_box(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else
					{
						ky=-1;
						while (virtual_position_y>simulation->state->length_y)
						{
							virtual_position_y-=simulation->parameters->mf_box_size_y;
							ky+=1;
						}
						while (ky>0)
						{
							virtual_position_y-=simulation->parameters->mf_box_size_y;
							ky-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrorx(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
				}
				else  // virtual position still in simulation box
				{
					mfVM_sort_in_neighbors_in_box(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
				}
			}
		}
		dx=position_x-((gint32)(position_x/simulation->parameters->mf_box_size_x))*simulation->parameters->mf_box_size_x;
		if (dx>simulation->parameters->mf_box_size_x*0.5)
			dx=simulation->parameters->mf_box_size_x-dx;
		dx=dx+(order-1)*simulation->parameters->mf_box_size_x;
		dy=position_y-((gint32)(position_y/simulation->parameters->mf_box_size_y))*simulation->parameters->mf_box_size_y;
		if (dy>simulation->parameters->mf_box_size_y*0.5)
			dy=simulation->parameters->mf_box_size_y-dy;
		dy=dy+(order-1)*simulation->parameters->mf_box_size_y;
		if (dx<dy)
			simulation->state->nh.max_distance_checked=dx;
		else
			simulation->state->nh.max_distance_checked=dy;
		order++;
	}
	//calculate orientational interaction with nearest neighbors
	for (i=0; i<simulation->parameters->mf_kn; i++)
	{
		//the neighboring particle is not 'behind' a reflecting boundary
		if (*(simulation->state->nh.mirror_flag+i)==0)
		{
			if (cos( (simulation->state->particles+index)->theta - (simulation->state->particles+*(simulation->state->nh.neighbor_indexes+i))->theta )>0.0 )
			{
				velocity_x+=cos((simulation->state->particles+*(simulation->state->nh.neighbor_indexes+i))->theta);
				velocity_y+=sin((simulation->state->particles+*(simulation->state->nh.neighbor_indexes+i))->theta);
			}
			else
			{
				velocity_x-=cos((simulation->state->particles+*(simulation->state->nh.neighbor_indexes+i))->theta);
				velocity_y-=sin((simulation->state->particles+*(simulation->state->nh.neighbor_indexes+i))->theta);
			}
		}
		//the neighboring particle is 'behind' a reflecting boundary in x direction -> the angle has to be reflected accordingly: theta->pi-theta
		else if (*(simulation->state->nh.mirror_flag+i)==1)
		{
			if (cos( (simulation->state->particles+index)->theta -M_PI+ (simulation->state->particles+*(simulation->state->nh.neighbor_indexes+i))->theta )>0.0 )
			{
				velocity_x+=cos(M_PI-(simulation->state->particles+*(simulation->state->nh.neighbor_indexes+i))->theta);
				velocity_y+=sin(M_PI-(simulation->state->particles+*(simulation->state->nh.neighbor_indexes+i))->theta);
			}
			else
			{
				velocity_x-=cos(M_PI-(simulation->state->particles+*(simulation->state->nh.neighbor_indexes+i))->theta);
				velocity_y-=sin(M_PI-(simulation->state->particles+*(simulation->state->nh.neighbor_indexes+i))->theta);
			}
		}
		//the neighboring particle is 'behind' a reflecting boundary in y direction -> the angle has to be reflected accordingly: theta->-theta
		else if (*(simulation->state->nh.mirror_flag+i)==2)
		{
			if (cos( (simulation->state->particles+index)->theta + (simulation->state->particles+*(simulation->state->nh.neighbor_indexes+i))->theta )>0.0 )
			{
				velocity_x+=cos(-(simulation->state->particles+*(simulation->state->nh.neighbor_indexes+i))->theta);
				velocity_y+=sin(-(simulation->state->particles+*(simulation->state->nh.neighbor_indexes+i))->theta);
			}
			else
			{
				velocity_x-=cos(-(simulation->state->particles+*(simulation->state->nh.neighbor_indexes+i))->theta);
				velocity_y-=sin(-(simulation->state->particles+*(simulation->state->nh.neighbor_indexes+i))->theta);
			}
		}
		//the neighboring particle is 'behind' a reflecting boundary in both, x- and y-direction -> the angle has to be reflected accordingly: theta->pi+theta
		else if (*(simulation->state->nh.mirror_flag+i)==3)
		{
			if (cos( (simulation->state->particles+index)->theta -M_PI - (simulation->state->particles+*(simulation->state->nh.neighbor_indexes+i))->theta )>0.0 )
			{
				velocity_x+=cos(M_PI+(simulation->state->particles+*(simulation->state->nh.neighbor_indexes+i))->theta);
				velocity_y+=sin(M_PI+(simulation->state->particles+*(simulation->state->nh.neighbor_indexes+i))->theta);
			}
			else
			{
				velocity_x-=cos(M_PI+(simulation->state->particles+*(simulation->state->nh.neighbor_indexes+i))->theta);
				velocity_y-=sin(M_PI+(simulation->state->particles+*(simulation->state->nh.neighbor_indexes+i))->theta);
			}
		}
	}
	//calculate the new angle, considering the effect of interactions (calculated before), noise and chirality
	if (simulation->parameters->particle_species==1)
		(simulation->state->particles+index)->theta_new=atan2(velocity_y, velocity_x)+(g_rand_double(simulation->state->prng)*2*M_PI -M_PI)**simulation->parameters->noise_strength+*simulation->parameters->omega;
	else
		(simulation->state->particles+index)->theta_new=atan2(velocity_y, velocity_x)+(g_rand_double(simulation->state->prng)*2*M_PI -M_PI)**(simulation->parameters->noise_strength+get_species(index, simulation->parameters))+*(simulation->parameters->omega+get_species(index, simulation->parameters));
	while ((simulation->state->particles+index)->theta_new<-M_PI)
		(simulation->state->particles+index)->theta_new+=2.0*M_PI;
	while ((simulation->state->particles+index)->theta_new>M_PI)
		(simulation->state->particles+index)->theta_new-=2.0*M_PI;
	return simulation->parameters->mf_kn;
}

//calculates interactions of particle 'index' for metric free/topological Langevin model
gint32 mfL_interaction_with_neighbors(guint64 index, struct VM_simulation * simulation)
{
	mfVM_clean_neighborhood(simulation);
	gint32 order=1; // this says which boxes are checked: if order=1 only the inner box is checked, if order=2 the surrounding 8 boxes are checked and so on
	gint32 i, j;
	gint32 kx, ky;
	gdouble delta_v=0.0;
	gdouble virtual_position_x;
	gdouble virtual_position_y;
	gdouble position_x=(simulation->state->particles+index)->x;
	gdouble position_y=(simulation->state->particles+index)->y;
	gdouble dx=0.0;
	gdouble dy=0.0;
	//find mf_kn nearest neighbors of particle 'index'
	while (simulation->state->nh.max_distance_checked<*(simulation->state->nh.distances+simulation->parameters->mf_kn-1))
	{
		if (order==1)
		{
			i=0;
			j=0;
			virtual_position_x=(simulation->state->particles+index)->x+i*simulation->parameters->mf_box_size_x;
			if (virtual_position_x<0.0)
				virtual_position_x+=simulation->state->length_x;
			if (virtual_position_x>simulation->state->length_x)
				virtual_position_x-=simulation->state->length_x;
			virtual_position_y=(simulation->state->particles+index)->y+j*simulation->parameters->mf_box_size_y;
			if (virtual_position_y<0.0)
				virtual_position_y+=simulation->state->length_y;
			if (virtual_position_y>simulation->state->length_y)
				virtual_position_y-=simulation->state->length_y;
			mfVM_sort_in_neighbors_in_box(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
		}
		else
		{
			i=1-order;
			for (j=1-order; j<order; j++)
			{
				virtual_position_x=(simulation->state->particles+index)->x+i*simulation->parameters->mf_box_size_x;
				virtual_position_y=(simulation->state->particles+index)->y+j*simulation->parameters->mf_box_size_y;
				if ((virtual_position_x<0.0) && (virtual_position_y<0.0)) //out of simulation box to left bottom
				{
					if ((simulation->parameters->boundary_x==0)&& (simulation->parameters->boundary_y==0)) //periodic boundary conditions in both directions
					{
						virtual_position_x+=simulation->state->length_x;
						virtual_position_y+=simulation->state->length_y;
						mfVM_sort_in_neighbors_in_box(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else if((simulation->parameters->boundary_x==0)&& (simulation->parameters->boundary_y!=0))//periodic in x-direction and reflecting in y-direction
					{
						virtual_position_x+=simulation->state->length_x;
						ky=-1;
						while (virtual_position_y<0.0)
						{
							virtual_position_y+=simulation->parameters->mf_box_size_y;
							ky+=1;
						}
						while (ky>0)
						{
							virtual_position_y+=simulation->parameters->mf_box_size_y;
							ky-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrory(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else if((simulation->parameters->boundary_x!=0)&& (simulation->parameters->boundary_y==0)) //reflecting in x-direction, periodic in y-direction
					{
						virtual_position_y+=simulation->state->length_y;
						kx=-1;
						while (virtual_position_x<0.0)
						{
							virtual_position_x+=simulation->parameters->mf_box_size_x;
							kx+=1;
						}
						while (kx>0)
						{
							virtual_position_x+=simulation->parameters->mf_box_size_x;
							kx-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrorx(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else if((simulation->parameters->boundary_x!=0)&& (simulation->parameters->boundary_y!=0)) //reflecting boundary in both directions
					{
						kx=-1;
						while (virtual_position_x<0.0)
						{
							virtual_position_x+=simulation->parameters->mf_box_size_x;
							kx+=1;
						}
						while (kx>0)
						{
							virtual_position_x+=simulation->parameters->mf_box_size_x;
							kx-=1;
						}
						ky=-1;
						while (virtual_position_y<0.0)
						{
							virtual_position_y+=simulation->parameters->mf_box_size_y;
							ky+=1;
						}
						while (ky>0)
						{
							virtual_position_y+=simulation->parameters->mf_box_size_y;
							ky-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrorxy(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
				}
				else if ((virtual_position_x<0.0) && (virtual_position_y>simulation->state->length_y)) //out of simulation box to left top
				{
					if ((simulation->parameters->boundary_x==0)&& (simulation->parameters->boundary_y==0)) //periodic boundary conditions in both directions
					{
						virtual_position_x+=simulation->state->length_x;
						virtual_position_y-=simulation->state->length_y;
						mfVM_sort_in_neighbors_in_box(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else if((simulation->parameters->boundary_x==0)&& (simulation->parameters->boundary_y!=0))//periodic in x-direction and reflecting in y-direction
					{
						virtual_position_x+=simulation->state->length_x;
						ky=-1;
						while (virtual_position_y>simulation->state->length_y)
						{
							virtual_position_y-=simulation->parameters->mf_box_size_y;
							ky+=1;
						}
						while (ky>0)
						{
							virtual_position_y-=simulation->parameters->mf_box_size_y;
							ky-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrory(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else if((simulation->parameters->boundary_x!=0)&& (simulation->parameters->boundary_y==0)) //reflecting in x-direction, periodic in y-direction
					{
						virtual_position_y-=simulation->state->length_y;
						kx=-1;
						while (virtual_position_x<0.0)
						{
							virtual_position_x+=simulation->parameters->mf_box_size_x;
							kx+=1;
						}
						while (kx>0)
						{
							virtual_position_x+=simulation->parameters->mf_box_size_x;
							kx-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrorx(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else if((simulation->parameters->boundary_x!=0)&& (simulation->parameters->boundary_y!=0)) //reflecting boundary in both directions
					{
						kx=-1;
						while (virtual_position_x<0.0)
						{
							virtual_position_x+=simulation->parameters->mf_box_size_x;
							kx+=1;
						}
						while (kx>0)
						{
							virtual_position_x+=simulation->parameters->mf_box_size_x;
							kx-=1;
						}
						ky=-1;
						while (virtual_position_y>simulation->state->length_y)
						{
							virtual_position_y-=simulation->parameters->mf_box_size_y;
							ky+=1;
						}
						while (ky>0)
						{
							virtual_position_y-=simulation->parameters->mf_box_size_y;
							ky-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrorxy(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
				}
				else if ((virtual_position_x>simulation->state->length_x) && (virtual_position_y>simulation->state->length_y)) //out of simulation box to right top
				{
					if ((simulation->parameters->boundary_x==0)&& (simulation->parameters->boundary_y==0)) //periodic boundary conditions in both directions
					{
						virtual_position_x-=simulation->state->length_x;
						virtual_position_y-=simulation->state->length_y;
						mfVM_sort_in_neighbors_in_box(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else if((simulation->parameters->boundary_x==0)&& (simulation->parameters->boundary_y!=0))//periodic in x-direction and reflecting in y-direction
					{
						virtual_position_x-=simulation->state->length_x;
						ky=-1;
						while (virtual_position_y>simulation->state->length_y)
						{
							virtual_position_y-=simulation->parameters->mf_box_size_y;
							ky+=1;
						}
						while (ky>0)
						{
							virtual_position_y-=simulation->parameters->mf_box_size_y;
							ky-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrory(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else if((simulation->parameters->boundary_x!=0)&& (simulation->parameters->boundary_y==0)) //reflecting in x-direction, periodic in y-direction
					{
						virtual_position_y-=simulation->state->length_y;
						kx=-1;
						while (virtual_position_x>simulation->state->length_x)
						{
							virtual_position_x-=simulation->parameters->mf_box_size_x;
							kx+=1;
						}
						while (kx>0)
						{
							virtual_position_x-=simulation->parameters->mf_box_size_x;
							kx-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrorx(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else if((simulation->parameters->boundary_x!=0)&& (simulation->parameters->boundary_y!=0)) //reflecting boundary in both directions
					{
						kx=-1;
						while (virtual_position_x>simulation->state->length_x)
						{
							virtual_position_x-=simulation->parameters->mf_box_size_x;
							kx+=1;
						}
						while (kx>0)
						{
							virtual_position_x-=simulation->parameters->mf_box_size_x;
							kx-=1;
						}
						ky=-1;
						while (virtual_position_y>simulation->state->length_y)
						{
							virtual_position_y-=simulation->parameters->mf_box_size_y;
							ky+=1;
						}
						while (ky>0)
						{
							virtual_position_y-=simulation->parameters->mf_box_size_y;
							ky-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrorxy(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
				}
				else if ((virtual_position_x>simulation->state->length_x) && (virtual_position_y<0.0)) //out of simulation box to right bottom
				{
					if ((simulation->parameters->boundary_x==0)&& (simulation->parameters->boundary_y==0)) //periodic boundary conditions in both directions
					{
						virtual_position_x-=simulation->state->length_x;
						virtual_position_y+=simulation->state->length_y;
						mfVM_sort_in_neighbors_in_box(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else if((simulation->parameters->boundary_x==0)&& (simulation->parameters->boundary_y!=0))//periodic in x-direction and reflecting in y-direction
					{
						virtual_position_x-=simulation->state->length_x;
						ky=-1;
						while (virtual_position_y<0.0)
						{
							virtual_position_y+=simulation->parameters->mf_box_size_y;
							ky+=1;
						}
						while (ky>0)
						{
							virtual_position_y+=simulation->parameters->mf_box_size_y;
							ky-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrory(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else if((simulation->parameters->boundary_x!=0)&& (simulation->parameters->boundary_y==0)) //reflecting in x-direction, periodic in y-direction
					{
						virtual_position_y+=simulation->state->length_y;
						kx=-1;
						while (virtual_position_x>simulation->state->length_x)
						{
							virtual_position_x-=simulation->parameters->mf_box_size_x;
							kx+=1;
						}
						while (kx>0)
						{
							virtual_position_x-=simulation->parameters->mf_box_size_x;
							kx-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrorx(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else if((simulation->parameters->boundary_x!=0)&& (simulation->parameters->boundary_y!=0)) //reflecting boundary in both directions
					{
						kx=-1;
						while (virtual_position_x>simulation->state->length_x)
						{
							virtual_position_x-=simulation->parameters->mf_box_size_x;
							kx+=1;
						}
						while (kx>0)
						{
							virtual_position_x-=simulation->parameters->mf_box_size_x;
							kx-=1;
						}
						ky=-1;
						while (virtual_position_y<0.0)
						{
							virtual_position_y+=simulation->parameters->mf_box_size_y;
							ky+=1;
						}
						while (ky>0)
						{
							virtual_position_y+=simulation->parameters->mf_box_size_y;
							ky-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrorxy(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
				}
				else if (virtual_position_x<0.0) //out of simulation box to the left
				{
					if (simulation->parameters->boundary_x==0)
					{
						virtual_position_x+=simulation->state->length_x;
						mfVM_sort_in_neighbors_in_box(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else
					{
						kx=-1;
						while (virtual_position_x<0.0)
						{
							virtual_position_x+=simulation->parameters->mf_box_size_x;
							kx+=1;
						}
						while (kx>0)
						{
							virtual_position_x+=simulation->parameters->mf_box_size_x;
							kx-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrorx(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
				}
				else if (virtual_position_x>simulation->state->length_x) //out of simulation box to the right
				{
					if (simulation->parameters->boundary_x==0)
					{
						virtual_position_x-=simulation->state->length_x;
						mfVM_sort_in_neighbors_in_box(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else
					{
						kx=-1;
						while (virtual_position_x>simulation->state->length_x)
						{
							virtual_position_x-=simulation->parameters->mf_box_size_x;
							kx+=1;
						}
						while (kx>0)
						{
							virtual_position_x-=simulation->parameters->mf_box_size_x;
							kx-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrorx(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
				}
				else if (virtual_position_y<0.0) //out of simulation box to the bottom
				{
					if (simulation->parameters->boundary_y==0)
					{
						virtual_position_y+=simulation->state->length_y;
						mfVM_sort_in_neighbors_in_box(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else
					{
						ky=-1;
						while (virtual_position_y<0.0)
						{
							virtual_position_y+=simulation->parameters->mf_box_size_y;
							ky+=1;
						}
						while (ky>0)
						{
							virtual_position_y+=simulation->parameters->mf_box_size_y;
							ky-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrory(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
				}
				else if (virtual_position_y>simulation->state->length_y) //out of simulation box to the top
				{
					if (simulation->parameters->boundary_y==0)
					{
						virtual_position_y-=simulation->state->length_y;
						mfVM_sort_in_neighbors_in_box(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else
					{
						ky=-1;
						while (virtual_position_y>simulation->state->length_y)
						{
							virtual_position_y-=simulation->parameters->mf_box_size_y;
							ky+=1;
						}
						while (ky>0)
						{
							virtual_position_y-=simulation->parameters->mf_box_size_y;
							ky-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrorx(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
				}
				else  // virtual position still in simulation box
				{
					mfVM_sort_in_neighbors_in_box(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
				}
			}
			i=order-1;
			for (j=1-order; j<order; j++)
			{
				virtual_position_x=(simulation->state->particles+index)->x+i*simulation->parameters->mf_box_size_x;
				virtual_position_y=(simulation->state->particles+index)->y+j*simulation->parameters->mf_box_size_y;
				if ((virtual_position_x<0.0) && (virtual_position_y<0.0)) //out of simulation box to left bottom
				{
					if ((simulation->parameters->boundary_x==0)&& (simulation->parameters->boundary_y==0)) //periodic boundary conditions in both directions
					{
						virtual_position_x+=simulation->state->length_x;
						virtual_position_y+=simulation->state->length_y;
						mfVM_sort_in_neighbors_in_box(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else if((simulation->parameters->boundary_x==0)&& (simulation->parameters->boundary_y!=0))//periodic in x-direction and reflecting in y-direction
					{
						virtual_position_x+=simulation->state->length_x;
						ky=-1;
						while (virtual_position_y<0.0)
						{
							virtual_position_y+=simulation->parameters->mf_box_size_y;
							ky+=1;
						}
						while (ky>0)
						{
							virtual_position_y+=simulation->parameters->mf_box_size_y;
							ky-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrory(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else if((simulation->parameters->boundary_x!=0)&& (simulation->parameters->boundary_y==0)) //reflecting in x-direction, periodic in y-direction
					{
						virtual_position_y+=simulation->state->length_y;
						kx=-1;
						while (virtual_position_x<0.0)
						{
							virtual_position_x+=simulation->parameters->mf_box_size_x;
							kx+=1;
						}
						while (kx>0)
						{
							virtual_position_x+=simulation->parameters->mf_box_size_x;
							kx-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrorx(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else if((simulation->parameters->boundary_x!=0)&& (simulation->parameters->boundary_y!=0)) //reflecting boundary in both directions
					{
						kx=-1;
						while (virtual_position_x<0.0)
						{
							virtual_position_x+=simulation->parameters->mf_box_size_x;
							kx+=1;
						}
						while (kx>0)
						{
							virtual_position_x+=simulation->parameters->mf_box_size_x;
							kx-=1;
						}
						ky=-1;
						while (virtual_position_y<0.0)
						{
							virtual_position_y+=simulation->parameters->mf_box_size_y;
							ky+=1;
						}
						while (ky>0)
						{
							virtual_position_y+=simulation->parameters->mf_box_size_y;
							ky-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrorxy(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
				}
				else if ((virtual_position_x<0.0) && (virtual_position_y>simulation->state->length_y)) //out of simulation box to left top
				{
					if ((simulation->parameters->boundary_x==0)&& (simulation->parameters->boundary_y==0)) //periodic boundary conditions in both directions
					{
						virtual_position_x+=simulation->state->length_x;
						virtual_position_y-=simulation->state->length_y;
						mfVM_sort_in_neighbors_in_box(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else if((simulation->parameters->boundary_x==0)&& (simulation->parameters->boundary_y!=0))//periodic in x-direction and reflecting in y-direction
					{
						virtual_position_x+=simulation->state->length_x;
						ky=-1;
						while (virtual_position_y>simulation->state->length_y)
						{
							virtual_position_y-=simulation->parameters->mf_box_size_y;
							ky+=1;
						}
						while (ky>0)
						{
							virtual_position_y-=simulation->parameters->mf_box_size_y;
							ky-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrory(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else if((simulation->parameters->boundary_x!=0)&& (simulation->parameters->boundary_y==0)) //reflecting in x-direction, periodic in y-direction
					{
						virtual_position_y-=simulation->state->length_y;
						kx=-1;
						while (virtual_position_x<0.0)
						{
							virtual_position_x+=simulation->parameters->mf_box_size_x;
							kx+=1;
						}
						while (kx>0)
						{
							virtual_position_x+=simulation->parameters->mf_box_size_x;
							kx-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrorx(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else if((simulation->parameters->boundary_x!=0)&& (simulation->parameters->boundary_y!=0)) //reflecting boundary in both directions
					{
						kx=-1;
						while (virtual_position_x<0.0)
						{
							virtual_position_x+=simulation->parameters->mf_box_size_x;
							kx+=1;
						}
						while (kx>0)
						{
							virtual_position_x+=simulation->parameters->mf_box_size_x;
							kx-=1;
						}
						ky=-1;
						while (virtual_position_y>simulation->state->length_y)
						{
							virtual_position_y-=simulation->parameters->mf_box_size_y;
							ky+=1;
						}
						while (ky>0)
						{
							virtual_position_y-=simulation->parameters->mf_box_size_y;
							ky-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrorxy(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
				}
				else if ((virtual_position_x>simulation->state->length_x) && (virtual_position_y>simulation->state->length_y)) //out of simulation box to right top
				{
					if ((simulation->parameters->boundary_x==0)&& (simulation->parameters->boundary_y==0)) //periodic boundary conditions in both directions
					{
						virtual_position_x-=simulation->state->length_x;
						virtual_position_y-=simulation->state->length_y;
						mfVM_sort_in_neighbors_in_box(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else if((simulation->parameters->boundary_x==0)&& (simulation->parameters->boundary_y!=0))//periodic in x-direction and reflecting in y-direction
					{
						virtual_position_x-=simulation->state->length_x;
						ky=-1;
						while (virtual_position_y>simulation->state->length_y)
						{
							virtual_position_y-=simulation->parameters->mf_box_size_y;
							ky+=1;
						}
						while (ky>0)
						{
							virtual_position_y-=simulation->parameters->mf_box_size_y;
							ky-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrory(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else if((simulation->parameters->boundary_x!=0)&& (simulation->parameters->boundary_y==0)) //reflecting in x-direction, periodic in y-direction
					{
						virtual_position_y-=simulation->state->length_y;
						kx=-1;
						while (virtual_position_x>simulation->state->length_x)
						{
							virtual_position_x-=simulation->parameters->mf_box_size_x;
							kx+=1;
						}
						while (kx>0)
						{
							virtual_position_x-=simulation->parameters->mf_box_size_x;
							kx-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrorx(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else if((simulation->parameters->boundary_x!=0)&& (simulation->parameters->boundary_y!=0)) //reflecting boundary in both directions
					{
						kx=-1;
						while (virtual_position_x>simulation->state->length_x)
						{
							virtual_position_x-=simulation->parameters->mf_box_size_x;
							kx+=1;
						}
						while (kx>0)
						{
							virtual_position_x-=simulation->parameters->mf_box_size_x;
							kx-=1;
						}
						ky=-1;
						while (virtual_position_y>simulation->state->length_y)
						{
							virtual_position_y-=simulation->parameters->mf_box_size_y;
							ky+=1;
						}
						while (ky>0)
						{
							virtual_position_y-=simulation->parameters->mf_box_size_y;
							ky-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrorxy(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
				}
				else if ((virtual_position_x>simulation->state->length_x) && (virtual_position_y<0.0)) //out of simulation box to right bottom
				{
					if ((simulation->parameters->boundary_x==0)&& (simulation->parameters->boundary_y==0)) //periodic boundary conditions in both directions
					{
						virtual_position_x-=simulation->state->length_x;
						virtual_position_y+=simulation->state->length_y;
						mfVM_sort_in_neighbors_in_box(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else if((simulation->parameters->boundary_x==0)&& (simulation->parameters->boundary_y!=0))//periodic in x-direction and reflecting in y-direction
					{
						virtual_position_x-=simulation->state->length_x;
						ky=-1;
						while (virtual_position_y<0.0)
						{
							virtual_position_y+=simulation->parameters->mf_box_size_y;
							ky+=1;
						}
						while (ky>0)
						{
							virtual_position_y+=simulation->parameters->mf_box_size_y;
							ky-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrory(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else if((simulation->parameters->boundary_x!=0)&& (simulation->parameters->boundary_y==0)) //reflecting in x-direction, periodic in y-direction
					{
						virtual_position_y+=simulation->state->length_y;
						kx=-1;
						while (virtual_position_x>simulation->state->length_x)
						{
							virtual_position_x-=simulation->parameters->mf_box_size_x;
							kx+=1;
						}
						while (kx>0)
						{
							virtual_position_x-=simulation->parameters->mf_box_size_x;
							kx-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrorx(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else if((simulation->parameters->boundary_x!=0)&& (simulation->parameters->boundary_y!=0)) //reflecting boundary in both directions
					{
						kx=-1;
						while (virtual_position_x>simulation->state->length_x)
						{
							virtual_position_x-=simulation->parameters->mf_box_size_x;
							kx+=1;
						}
						while (kx>0)
						{
							virtual_position_x-=simulation->parameters->mf_box_size_x;
							kx-=1;
						}
						ky=-1;
						while (virtual_position_y<0.0)
						{
							virtual_position_y+=simulation->parameters->mf_box_size_y;
							ky+=1;
						}
						while (ky>0)
						{
							virtual_position_y+=simulation->parameters->mf_box_size_y;
							ky-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrorxy(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
				}
				else if (virtual_position_x<0.0) //out of simulation box to the left
				{
					if (simulation->parameters->boundary_x==0)
					{
						virtual_position_x+=simulation->state->length_x;
						mfVM_sort_in_neighbors_in_box(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else
					{
						kx=-1;
						while (virtual_position_x<0.0)
						{
							virtual_position_x+=simulation->parameters->mf_box_size_x;
							kx+=1;
						}
						while (kx>0)
						{
							virtual_position_x+=simulation->parameters->mf_box_size_x;
							kx-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrorx(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
				}
				else if (virtual_position_x>simulation->state->length_x) //out of simulation box to the right
				{
					if (simulation->parameters->boundary_x==0)
					{
						virtual_position_x-=simulation->state->length_x;
						mfVM_sort_in_neighbors_in_box(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else
					{
						kx=-1;
						while (virtual_position_x>simulation->state->length_x)
						{
							virtual_position_x-=simulation->parameters->mf_box_size_x;
							kx+=1;
						}
						while (kx>0)
						{
							virtual_position_x-=simulation->parameters->mf_box_size_x;
							kx-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrorx(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
				}
				else if (virtual_position_y<0.0) //out of simulation box to the bottom
				{
					if (simulation->parameters->boundary_y==0)
					{
						virtual_position_y+=simulation->state->length_y;
						mfVM_sort_in_neighbors_in_box(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else
					{
						ky=-1;
						while (virtual_position_y<0.0)
						{
							virtual_position_y+=simulation->parameters->mf_box_size_y;
							ky+=1;
						}
						while (ky>0)
						{
							virtual_position_y+=simulation->parameters->mf_box_size_y;
							ky-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrory(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
				}
				else if (virtual_position_y>simulation->state->length_y) //out of simulation box to the top
				{
					if (simulation->parameters->boundary_y==0)
					{
						virtual_position_y-=simulation->state->length_y;
						mfVM_sort_in_neighbors_in_box(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else
					{
						ky=-1;
						while (virtual_position_y>simulation->state->length_y)
						{
							virtual_position_y-=simulation->parameters->mf_box_size_y;
							ky+=1;
						}
						while (ky>0)
						{
							virtual_position_y-=simulation->parameters->mf_box_size_y;
							ky-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrorx(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
				}
				else  // virtual position still in simulation box
				{
					mfVM_sort_in_neighbors_in_box(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
				}
			}
			j=1-order;
			for (i=2-order; i<order-1; i++)
			{
				virtual_position_x=(simulation->state->particles+index)->x+i*simulation->parameters->mf_box_size_x;
				virtual_position_y=(simulation->state->particles+index)->y+j*simulation->parameters->mf_box_size_y;
				if ((virtual_position_x<0.0) && (virtual_position_y<0.0)) //out of simulation box to left bottom
				{
					if ((simulation->parameters->boundary_x==0)&& (simulation->parameters->boundary_y==0)) //periodic boundary conditions in both directions
					{
						virtual_position_x+=simulation->state->length_x;
						virtual_position_y+=simulation->state->length_y;
						mfVM_sort_in_neighbors_in_box(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else if((simulation->parameters->boundary_x==0)&& (simulation->parameters->boundary_y!=0))//periodic in x-direction and reflecting in y-direction
					{
						virtual_position_x+=simulation->state->length_x;
						ky=-1;
						while (virtual_position_y<0.0)
						{
							virtual_position_y+=simulation->parameters->mf_box_size_y;
							ky+=1;
						}
						while (ky>0)
						{
							virtual_position_y+=simulation->parameters->mf_box_size_y;
							ky-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrory(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else if((simulation->parameters->boundary_x!=0)&& (simulation->parameters->boundary_y==0)) //reflecting in x-direction, periodic in y-direction
					{
						virtual_position_y+=simulation->state->length_y;
						kx=-1;
						while (virtual_position_x<0.0)
						{
							virtual_position_x+=simulation->parameters->mf_box_size_x;
							kx+=1;
						}
						while (kx>0)
						{
							virtual_position_x+=simulation->parameters->mf_box_size_x;
							kx-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrorx(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else if((simulation->parameters->boundary_x!=0)&& (simulation->parameters->boundary_y!=0)) //reflecting boundary in both directions
					{
						kx=-1;
						while (virtual_position_x<0.0)
						{
							virtual_position_x+=simulation->parameters->mf_box_size_x;
							kx+=1;
						}
						while (kx>0)
						{
							virtual_position_x+=simulation->parameters->mf_box_size_x;
							kx-=1;
						}
						ky=-1;
						while (virtual_position_y<0.0)
						{
							virtual_position_y+=simulation->parameters->mf_box_size_y;
							ky+=1;
						}
						while (ky>0)
						{
							virtual_position_y+=simulation->parameters->mf_box_size_y;
							ky-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrorxy(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
				}
				else if ((virtual_position_x<0.0) && (virtual_position_y>simulation->state->length_y)) //out of simulation box to left top
				{
					if ((simulation->parameters->boundary_x==0)&& (simulation->parameters->boundary_y==0)) //periodic boundary conditions in both directions
					{
						virtual_position_x+=simulation->state->length_x;
						virtual_position_y-=simulation->state->length_y;
						mfVM_sort_in_neighbors_in_box(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else if((simulation->parameters->boundary_x==0)&& (simulation->parameters->boundary_y!=0))//periodic in x-direction and reflecting in y-direction
					{
						virtual_position_x+=simulation->state->length_x;
						ky=-1;
						while (virtual_position_y>simulation->state->length_y)
						{
							virtual_position_y-=simulation->parameters->mf_box_size_y;
							ky+=1;
						}
						while (ky>0)
						{
							virtual_position_y-=simulation->parameters->mf_box_size_y;
							ky-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrory(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else if((simulation->parameters->boundary_x!=0)&& (simulation->parameters->boundary_y==0)) //reflecting in x-direction, periodic in y-direction
					{
						virtual_position_y-=simulation->state->length_y;
						kx=-1;
						while (virtual_position_x<0.0)
						{
							virtual_position_x+=simulation->parameters->mf_box_size_x;
							kx+=1;
						}
						while (kx>0)
						{
							virtual_position_x+=simulation->parameters->mf_box_size_x;
							kx-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrorx(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else if((simulation->parameters->boundary_x!=0)&& (simulation->parameters->boundary_y!=0)) //reflecting boundary in both directions
					{
						kx=-1;
						while (virtual_position_x<0.0)
						{
							virtual_position_x+=simulation->parameters->mf_box_size_x;
							kx+=1;
						}
						while (kx>0)
						{
							virtual_position_x+=simulation->parameters->mf_box_size_x;
							kx-=1;
						}
						ky=-1;
						while (virtual_position_y>simulation->state->length_y)
						{
							virtual_position_y-=simulation->parameters->mf_box_size_y;
							ky+=1;
						}
						while (ky>0)
						{
							virtual_position_y-=simulation->parameters->mf_box_size_y;
							ky-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrorxy(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
				}
				else if ((virtual_position_x>simulation->state->length_x) && (virtual_position_y>simulation->state->length_y)) //out of simulation box to right top
				{
					if ((simulation->parameters->boundary_x==0)&& (simulation->parameters->boundary_y==0)) //periodic boundary conditions in both directions
					{
						virtual_position_x-=simulation->state->length_x;
						virtual_position_y-=simulation->state->length_y;
						mfVM_sort_in_neighbors_in_box(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else if((simulation->parameters->boundary_x==0)&& (simulation->parameters->boundary_y!=0))//periodic in x-direction and reflecting in y-direction
					{
						virtual_position_x-=simulation->state->length_x;
						ky=-1;
						while (virtual_position_y>simulation->state->length_y)
						{
							virtual_position_y-=simulation->parameters->mf_box_size_y;
							ky+=1;
						}
						while (ky>0)
						{
							virtual_position_y-=simulation->parameters->mf_box_size_y;
							ky-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrory(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else if((simulation->parameters->boundary_x!=0)&& (simulation->parameters->boundary_y==0)) //reflecting in x-direction, periodic in y-direction
					{
						virtual_position_y-=simulation->state->length_y;
						kx=-1;
						while (virtual_position_x>simulation->state->length_x)
						{
							virtual_position_x-=simulation->parameters->mf_box_size_x;
							kx+=1;
						}
						while (kx>0)
						{
							virtual_position_x-=simulation->parameters->mf_box_size_x;
							kx-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrorx(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else if((simulation->parameters->boundary_x!=0)&& (simulation->parameters->boundary_y!=0)) //reflecting boundary in both directions
					{
						kx=-1;
						while (virtual_position_x>simulation->state->length_x)
						{
							virtual_position_x-=simulation->parameters->mf_box_size_x;
							kx+=1;
						}
						while (kx>0)
						{
							virtual_position_x-=simulation->parameters->mf_box_size_x;
							kx-=1;
						}
						ky=-1;
						while (virtual_position_y>simulation->state->length_y)
						{
							virtual_position_y-=simulation->parameters->mf_box_size_y;
							ky+=1;
						}
						while (ky>0)
						{
							virtual_position_y-=simulation->parameters->mf_box_size_y;
							ky-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrorxy(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
				}
				else if ((virtual_position_x>simulation->state->length_x) && (virtual_position_y<0.0)) //out of simulation box to right bottom
				{
					if ((simulation->parameters->boundary_x==0)&& (simulation->parameters->boundary_y==0)) //periodic boundary conditions in both directions
					{
						virtual_position_x-=simulation->state->length_x;
						virtual_position_y+=simulation->state->length_y;
						mfVM_sort_in_neighbors_in_box(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else if((simulation->parameters->boundary_x==0)&& (simulation->parameters->boundary_y!=0))//periodic in x-direction and reflecting in y-direction
					{
						virtual_position_x-=simulation->state->length_x;
						ky=-1;
						while (virtual_position_y<0.0)
						{
							virtual_position_y+=simulation->parameters->mf_box_size_y;
							ky+=1;
						}
						while (ky>0)
						{
							virtual_position_y+=simulation->parameters->mf_box_size_y;
							ky-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrory(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else if((simulation->parameters->boundary_x!=0)&& (simulation->parameters->boundary_y==0)) //reflecting in x-direction, periodic in y-direction
					{
						virtual_position_y+=simulation->state->length_y;
						kx=-1;
						while (virtual_position_x>simulation->state->length_x)
						{
							virtual_position_x-=simulation->parameters->mf_box_size_x;
							kx+=1;
						}
						while (kx>0)
						{
							virtual_position_x-=simulation->parameters->mf_box_size_x;
							kx-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrorx(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else if((simulation->parameters->boundary_x!=0)&& (simulation->parameters->boundary_y!=0)) //reflecting boundary in both directions
					{
						kx=-1;
						while (virtual_position_x>simulation->state->length_x)
						{
							virtual_position_x-=simulation->parameters->mf_box_size_x;
							kx+=1;
						}
						while (kx>0)
						{
							virtual_position_x-=simulation->parameters->mf_box_size_x;
							kx-=1;
						}
						ky=-1;
						while (virtual_position_y<0.0)
						{
							virtual_position_y+=simulation->parameters->mf_box_size_y;
							ky+=1;
						}
						while (ky>0)
						{
							virtual_position_y+=simulation->parameters->mf_box_size_y;
							ky-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrorxy(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
				}
				else if (virtual_position_x<0.0) //out of simulation box to the left
				{
					if (simulation->parameters->boundary_x==0)
					{
						virtual_position_x+=simulation->state->length_x;
						mfVM_sort_in_neighbors_in_box(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else
					{
						kx=-1;
						while (virtual_position_x<0.0)
						{
							virtual_position_x+=simulation->parameters->mf_box_size_x;
							kx+=1;
						}
						while (kx>0)
						{
							virtual_position_x+=simulation->parameters->mf_box_size_x;
							kx-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrorx(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
				}
				else if (virtual_position_x>simulation->state->length_x) //out of simulation box to the right
				{
					if (simulation->parameters->boundary_x==0)
					{
						virtual_position_x-=simulation->state->length_x;
						mfVM_sort_in_neighbors_in_box(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else
					{
						kx=-1;
						while (virtual_position_x>simulation->state->length_x)
						{
							virtual_position_x-=simulation->parameters->mf_box_size_x;
							kx+=1;
						}
						while (kx>0)
						{
							virtual_position_x-=simulation->parameters->mf_box_size_x;
							kx-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrorx(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
				}
				else if (virtual_position_y<0.0) //out of simulation box to the bottom
				{
					if (simulation->parameters->boundary_y==0)
					{
						virtual_position_y+=simulation->state->length_y;
						mfVM_sort_in_neighbors_in_box(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else
					{
						ky=-1;
						while (virtual_position_y<0.0)
						{
							virtual_position_y+=simulation->parameters->mf_box_size_y;
							ky+=1;
						}
						while (ky>0)
						{
							virtual_position_y+=simulation->parameters->mf_box_size_y;
							ky-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrory(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
				}
				else if (virtual_position_y>simulation->state->length_y) //out of simulation box to the top
				{
					if (simulation->parameters->boundary_y==0)
					{
						virtual_position_y-=simulation->state->length_y;
						mfVM_sort_in_neighbors_in_box(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else
					{
						ky=-1;
						while (virtual_position_y>simulation->state->length_y)
						{
							virtual_position_y-=simulation->parameters->mf_box_size_y;
							ky+=1;
						}
						while (ky>0)
						{
							virtual_position_y-=simulation->parameters->mf_box_size_y;
							ky-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrorx(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
				}
				else  // virtual position still in simulation box
				{
					mfVM_sort_in_neighbors_in_box(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
				}
			}
			j=order-1;
			for (i=2-order; i<order-1; i++)
			{
				virtual_position_x=(simulation->state->particles+index)->x+i*simulation->parameters->mf_box_size_x;
				virtual_position_y=(simulation->state->particles+index)->y+j*simulation->parameters->mf_box_size_y;
				if ((virtual_position_x<0.0) && (virtual_position_y<0.0)) //out of simulation box to left bottom
				{
					if ((simulation->parameters->boundary_x==0)&& (simulation->parameters->boundary_y==0)) //periodic boundary conditions in both directions
					{
						virtual_position_x+=simulation->state->length_x;
						virtual_position_y+=simulation->state->length_y;
						mfVM_sort_in_neighbors_in_box(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else if((simulation->parameters->boundary_x==0)&& (simulation->parameters->boundary_y!=0))//periodic in x-direction and reflecting in y-direction
					{
						virtual_position_x+=simulation->state->length_x;
						ky=-1;
						while (virtual_position_y<0.0)
						{
							virtual_position_y+=simulation->parameters->mf_box_size_y;
							ky+=1;
						}
						while (ky>0)
						{
							virtual_position_y+=simulation->parameters->mf_box_size_y;
							ky-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrory(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else if((simulation->parameters->boundary_x!=0)&& (simulation->parameters->boundary_y==0)) //reflecting in x-direction, periodic in y-direction
					{
						virtual_position_y+=simulation->state->length_y;
						kx=-1;
						while (virtual_position_x<0.0)
						{
							virtual_position_x+=simulation->parameters->mf_box_size_x;
							kx+=1;
						}
						while (kx>0)
						{
							virtual_position_x+=simulation->parameters->mf_box_size_x;
							kx-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrorx(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else if((simulation->parameters->boundary_x!=0)&& (simulation->parameters->boundary_y!=0)) //reflecting boundary in both directions
					{
						kx=-1;
						while (virtual_position_x<0.0)
						{
							virtual_position_x+=simulation->parameters->mf_box_size_x;
							kx+=1;
						}
						while (kx>0)
						{
							virtual_position_x+=simulation->parameters->mf_box_size_x;
							kx-=1;
						}
						ky=-1;
						while (virtual_position_y<0.0)
						{
							virtual_position_y+=simulation->parameters->mf_box_size_y;
							ky+=1;
						}
						while (ky>0)
						{
							virtual_position_y+=simulation->parameters->mf_box_size_y;
							ky-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrorxy(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
				}
				else if ((virtual_position_x<0.0) && (virtual_position_y>simulation->state->length_y)) //out of simulation box to left top
				{
					if ((simulation->parameters->boundary_x==0)&& (simulation->parameters->boundary_y==0)) //periodic boundary conditions in both directions
					{
						virtual_position_x+=simulation->state->length_x;
						virtual_position_y-=simulation->state->length_y;
						mfVM_sort_in_neighbors_in_box(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else if((simulation->parameters->boundary_x==0)&& (simulation->parameters->boundary_y!=0))//periodic in x-direction and reflecting in y-direction
					{
						virtual_position_x+=simulation->state->length_x;
						ky=-1;
						while (virtual_position_y>simulation->state->length_y)
						{
							virtual_position_y-=simulation->parameters->mf_box_size_y;
							ky+=1;
						}
						while (ky>0)
						{
							virtual_position_y-=simulation->parameters->mf_box_size_y;
							ky-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrory(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else if((simulation->parameters->boundary_x!=0)&& (simulation->parameters->boundary_y==0)) //reflecting in x-direction, periodic in y-direction
					{
						virtual_position_y-=simulation->state->length_y;
						kx=-1;
						while (virtual_position_x<0.0)
						{
							virtual_position_x+=simulation->parameters->mf_box_size_x;
							kx+=1;
						}
						while (kx>0)
						{
							virtual_position_x+=simulation->parameters->mf_box_size_x;
							kx-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrorx(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else if((simulation->parameters->boundary_x!=0)&& (simulation->parameters->boundary_y!=0)) //reflecting boundary in both directions
					{
						kx=-1;
						while (virtual_position_x<0.0)
						{
							virtual_position_x+=simulation->parameters->mf_box_size_x;
							kx+=1;
						}
						while (kx>0)
						{
							virtual_position_x+=simulation->parameters->mf_box_size_x;
							kx-=1;
						}
						ky=-1;
						while (virtual_position_y>simulation->state->length_y)
						{
							virtual_position_y-=simulation->parameters->mf_box_size_y;
							ky+=1;
						}
						while (ky>0)
						{
							virtual_position_y-=simulation->parameters->mf_box_size_y;
							ky-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrorxy(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
				}
				else if ((virtual_position_x>simulation->state->length_x) && (virtual_position_y>simulation->state->length_y)) //out of simulation box to right top
				{
					if ((simulation->parameters->boundary_x==0)&& (simulation->parameters->boundary_y==0)) //periodic boundary conditions in both directions
					{
						virtual_position_x-=simulation->state->length_x;
						virtual_position_y-=simulation->state->length_y;
						mfVM_sort_in_neighbors_in_box(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else if((simulation->parameters->boundary_x==0)&& (simulation->parameters->boundary_y!=0))//periodic in x-direction and reflecting in y-direction
					{
						virtual_position_x-=simulation->state->length_x;
						ky=-1;
						while (virtual_position_y>simulation->state->length_y)
						{
							virtual_position_y-=simulation->parameters->mf_box_size_y;
							ky+=1;
						}
						while (ky>0)
						{
							virtual_position_y-=simulation->parameters->mf_box_size_y;
							ky-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrory(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else if((simulation->parameters->boundary_x!=0)&& (simulation->parameters->boundary_y==0)) //reflecting in x-direction, periodic in y-direction
					{
						virtual_position_y-=simulation->state->length_y;
						kx=-1;
						while (virtual_position_x>simulation->state->length_x)
						{
							virtual_position_x-=simulation->parameters->mf_box_size_x;
							kx+=1;
						}
						while (kx>0)
						{
							virtual_position_x-=simulation->parameters->mf_box_size_x;
							kx-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrorx(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else if((simulation->parameters->boundary_x!=0)&& (simulation->parameters->boundary_y!=0)) //reflecting boundary in both directions
					{
						kx=-1;
						while (virtual_position_x>simulation->state->length_x)
						{
							virtual_position_x-=simulation->parameters->mf_box_size_x;
							kx+=1;
						}
						while (kx>0)
						{
							virtual_position_x-=simulation->parameters->mf_box_size_x;
							kx-=1;
						}
						ky=-1;
						while (virtual_position_y>simulation->state->length_y)
						{
							virtual_position_y-=simulation->parameters->mf_box_size_y;
							ky+=1;
						}
						while (ky>0)
						{
							virtual_position_y-=simulation->parameters->mf_box_size_y;
							ky-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrorxy(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
				}
				else if ((virtual_position_x>simulation->state->length_x) && (virtual_position_y<0.0)) //out of simulation box to right bottom
				{
					if ((simulation->parameters->boundary_x==0)&& (simulation->parameters->boundary_y==0)) //periodic boundary conditions in both directions
					{
						virtual_position_x-=simulation->state->length_x;
						virtual_position_y+=simulation->state->length_y;
						mfVM_sort_in_neighbors_in_box(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else if((simulation->parameters->boundary_x==0)&& (simulation->parameters->boundary_y!=0))//periodic in x-direction and reflecting in y-direction
					{
						virtual_position_x-=simulation->state->length_x;
						ky=-1;
						while (virtual_position_y<0.0)
						{
							virtual_position_y+=simulation->parameters->mf_box_size_y;
							ky+=1;
						}
						while (ky>0)
						{
							virtual_position_y+=simulation->parameters->mf_box_size_y;
							ky-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrory(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else if((simulation->parameters->boundary_x!=0)&& (simulation->parameters->boundary_y==0)) //reflecting in x-direction, periodic in y-direction
					{
						virtual_position_y+=simulation->state->length_y;
						kx=-1;
						while (virtual_position_x>simulation->state->length_x)
						{
							virtual_position_x-=simulation->parameters->mf_box_size_x;
							kx+=1;
						}
						while (kx>0)
						{
							virtual_position_x-=simulation->parameters->mf_box_size_x;
							kx-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrorx(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else if((simulation->parameters->boundary_x!=0)&& (simulation->parameters->boundary_y!=0)) //reflecting boundary in both directions
					{
						kx=-1;
						while (virtual_position_x>simulation->state->length_x)
						{
							virtual_position_x-=simulation->parameters->mf_box_size_x;
							kx+=1;
						}
						while (kx>0)
						{
							virtual_position_x-=simulation->parameters->mf_box_size_x;
							kx-=1;
						}
						ky=-1;
						while (virtual_position_y<0.0)
						{
							virtual_position_y+=simulation->parameters->mf_box_size_y;
							ky+=1;
						}
						while (ky>0)
						{
							virtual_position_y+=simulation->parameters->mf_box_size_y;
							ky-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrorxy(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
				}
				else if (virtual_position_x<0.0) //out of simulation box to the left
				{
					if (simulation->parameters->boundary_x==0)
					{
						virtual_position_x+=simulation->state->length_x;
						mfVM_sort_in_neighbors_in_box(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else
					{
						kx=-1;
						while (virtual_position_x<0.0)
						{
							virtual_position_x+=simulation->parameters->mf_box_size_x;
							kx+=1;
						}
						while (kx>0)
						{
							virtual_position_x+=simulation->parameters->mf_box_size_x;
							kx-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrorx(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
				}
				else if (virtual_position_x>simulation->state->length_x) //out of simulation box to the right
				{
					if (simulation->parameters->boundary_x==0)
					{
						virtual_position_x-=simulation->state->length_x;
						mfVM_sort_in_neighbors_in_box(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else
					{
						kx=-1;
						while (virtual_position_x>simulation->state->length_x)
						{
							virtual_position_x-=simulation->parameters->mf_box_size_x;
							kx+=1;
						}
						while (kx>0)
						{
							virtual_position_x-=simulation->parameters->mf_box_size_x;
							kx-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrorx(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
				}
				else if (virtual_position_y<0.0) //out of simulation box to the bottom
				{
					if (simulation->parameters->boundary_y==0)
					{
						virtual_position_y+=simulation->state->length_y;
						mfVM_sort_in_neighbors_in_box(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else
					{
						ky=-1;
						while (virtual_position_y<0.0)
						{
							virtual_position_y+=simulation->parameters->mf_box_size_y;
							ky+=1;
						}
						while (ky>0)
						{
							virtual_position_y+=simulation->parameters->mf_box_size_y;
							ky-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrory(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
				}
				else if (virtual_position_y>simulation->state->length_y) //out of simulation box to the top
				{
					if (simulation->parameters->boundary_y==0)
					{
						virtual_position_y-=simulation->state->length_y;
						mfVM_sort_in_neighbors_in_box(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
					else
					{
						ky=-1;
						while (virtual_position_y>simulation->state->length_y)
						{
							virtual_position_y-=simulation->parameters->mf_box_size_y;
							ky+=1;
						}
						while (ky>0)
						{
							virtual_position_y-=simulation->parameters->mf_box_size_y;
							ky-=1;
						}
						mfVM_sort_in_neighbors_in_box_mirrorx(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
					}
				}
				else  // virtual position still in simulation box
				{
					mfVM_sort_in_neighbors_in_box(mfVM_get_box_index(virtual_position_x, virtual_position_y, simulation->parameters), simulation, index);
				}
			}
		}
		dx=position_x-((gint32)(position_x/simulation->parameters->mf_box_size_x))*simulation->parameters->mf_box_size_x;
		if (dx>simulation->parameters->mf_box_size_x*0.5)
			dx=simulation->parameters->mf_box_size_x-dx;
		dx=dx+(order-1)*simulation->parameters->mf_box_size_x;
		dy=position_y-((gint32)(position_y/simulation->parameters->mf_box_size_y))*simulation->parameters->mf_box_size_y;
		if (dy>simulation->parameters->mf_box_size_y*0.5)
			dy=simulation->parameters->mf_box_size_y-dy;
		dy=dy+(order-1)*simulation->parameters->mf_box_size_y;
		if (dx<dy)
			simulation->state->nh.max_distance_checked=dx;
		else
			simulation->state->nh.max_distance_checked=dy;
		order++;
	}
	//calculate orientational interaction with nearest neighbors
	for (i=0; i<simulation->parameters->mf_kn; i++)
	{
		//the neighboring particle is not 'behind' a reflecting boundary
		if (*(simulation->state->nh.mirror_flag+i)==0)
		{
			if (simulation->parameters->particle_species==1)
				delta_v+=**simulation->parameters->coupling*sin( ((simulation->state->particles+*(simulation->state->nh.neighbor_indexes+i))->theta - (simulation->state->particles+index)->theta)*simulation->parameters->interaction_order  );
			else
				delta_v+=*(*(simulation->parameters->coupling+get_species(index,simulation->parameters))+get_species(*(simulation->state->nh.neighbor_indexes+i), simulation->parameters))*sin( ((simulation->state->particles+*(simulation->state->nh.neighbor_indexes+i))->theta - (simulation->state->particles+index)->theta)*simulation->parameters->interaction_order  );
		}
		//the neighboring particle is 'behind' a reflecting boundary in x direction -> the angle has to be reflected accordingly: theta->pi-theta
		else if (*(simulation->state->nh.mirror_flag+i)==1)
		{
			if (simulation->parameters->particle_species==1)
				delta_v+=**simulation->parameters->coupling*sin( (M_PI-(simulation->state->particles+*(simulation->state->nh.neighbor_indexes+i))->theta - (simulation->state->particles+index)->theta)*simulation->parameters->interaction_order  );
			else
				delta_v+=*(*(simulation->parameters->coupling+get_species(index,simulation->parameters))+get_species(*(simulation->state->nh.neighbor_indexes+i), simulation->parameters))*sin( (M_PI-(simulation->state->particles+*(simulation->state->nh.neighbor_indexes+i))->theta - (simulation->state->particles+index)->theta)*simulation->parameters->interaction_order  );
		}
		//the neighboring particle is 'behind' a reflecting boundary in y direction -> the angle has to be reflected accordingly: theta->-theta
		else if (*(simulation->state->nh.mirror_flag+i)==2)
		{
			if (simulation->parameters->particle_species==1)
				delta_v+=**simulation->parameters->coupling*sin( (-(simulation->state->particles+*(simulation->state->nh.neighbor_indexes+i))->theta - (simulation->state->particles+index)->theta)*simulation->parameters->interaction_order  );
			else
				delta_v+=*(*(simulation->parameters->coupling+get_species(index,simulation->parameters))+get_species(*(simulation->state->nh.neighbor_indexes+i), simulation->parameters))*sin( (-(simulation->state->particles+*(simulation->state->nh.neighbor_indexes+i))->theta - (simulation->state->particles+index)->theta)*simulation->parameters->interaction_order  );
		}
		//the neighboring particle is 'behind' a reflecting boundary in both, x- and y-direction -> the angle has to be reflected accordingly: theta->pi+theta
		else if (*(simulation->state->nh.mirror_flag+i)==3)
		{
			if (simulation->parameters->particle_species==1)
				delta_v+=**simulation->parameters->coupling*sin( (M_PI+(simulation->state->particles+*(simulation->state->nh.neighbor_indexes+i))->theta - (simulation->state->particles+index)->theta)*simulation->parameters->interaction_order  );
			else
				delta_v+=*(*(simulation->parameters->coupling+get_species(index,simulation->parameters))+get_species(*(simulation->state->nh.neighbor_indexes+i), simulation->parameters))*sin( (M_PI+(simulation->state->particles+*(simulation->state->nh.neighbor_indexes+i))->theta - (simulation->state->particles+index)->theta)*simulation->parameters->interaction_order  );
		}
	}
	//update the orientations according to interaction term, noise and chirality
	if (simulation->parameters->particle_species==1)
		(simulation->state->particles+index)->theta_new=(simulation->state->particles+index)->theta + simulation->parameters->delta_t*(*simulation->parameters->omega+delta_v) +simulation->parameters->sqrt_delta_t**simulation->parameters->noise_strength*gauss(simulation);
	else
		(simulation->state->particles+index)->theta_new=(simulation->state->particles+index)->theta + simulation->parameters->delta_t*(*(simulation->parameters->omega+get_species(index, simulation->parameters))+delta_v) +simulation->parameters->sqrt_delta_t**(simulation->parameters->noise_strength+get_species(index, simulation->parameters))*gauss(simulation);
	while ((simulation->state->particles+index)->theta_new<-M_PI)
		(simulation->state->particles+index)->theta_new+=2.0*M_PI;
	while ((simulation->state->particles+index)->theta_new>M_PI)
		(simulation->state->particles+index)->theta_new-=2.0*M_PI;
	return simulation->parameters->mf_kn;
}

//performs time evolution of one step for Vicsek model
void VM_one_step(struct VM_simulation * simulation)
{
	gint32 i;
	struct VM_particle * box;
	gint32 k;
	gdouble * polar_x=malloc(simulation->parameters->particle_species*sizeof(gdouble));
	gdouble * polar_y=malloc(simulation->parameters->particle_species*sizeof(gdouble));
	gdouble * nematic_x=malloc(simulation->parameters->particle_species*sizeof(gdouble));
	gdouble * nematic_y=malloc(simulation->parameters->particle_species*sizeof(gdouble));
	for (k=0; k<simulation->parameters->particle_species; k++)
	{
		*(polar_x+k)=0.0;
		*(polar_y+k)=0.0;
		*(nematic_x+k)=0.0;
		*(nematic_y+k)=0.0;
	}
	VM_fill_boxes(simulation);
	for (i=0; i<simulation->parameters->box_number_x*simulation->parameters->box_number_y; i++)
	{
		box=*(simulation->parameters->box+i);
		while (box !=NULL)
		{
			//alignment interaction
			VM_interaction_with_neighbors(box->index, simulation);
			box=box->next_particle;
		}
	}
	for (i=0; i<simulation->parameters->box_number_x*simulation->parameters->box_number_y; i++)
	{
		box=*(simulation->parameters->box+i);
		while (box !=NULL)
		{
			//calculate contribution to polar/nematic order parameter
			if (simulation->parameters->particle_species==1)
			{
				*polar_x+=cos(box->theta);
				*polar_y+=sin(box->theta);
				*nematic_x+=cos(2.0*box->theta);
				*nematic_y+=sin(2.0*box->theta);
			}
			else
			{
				*(polar_x+get_species(box->index, simulation->parameters))+=cos(box->theta);
				*(polar_y+get_species(box->index, simulation->parameters))+=sin(box->theta);
				*(nematic_x+get_species(box->index, simulation->parameters))+=cos(2.0*box->theta);
				*(nematic_y+get_species(box->index, simulation->parameters))+=sin(2.0*box->theta);
			}
			//update angle
			box->theta=box->theta_new;
			//streaming -> movement of the particles
			if (simulation->parameters->particle_species==1)
				box->x+=*simulation->parameters->speed*cos(box->theta);
			else
				box->x+=*(simulation->parameters->speed+get_species(box->index, simulation->parameters))*cos(box->theta);
			if (simulation->parameters->boundary_x!=0) //reflecting boundary in x-direction
				if ((box->x<0.0) || (box->x>simulation->state->length_x))
				{
					box->x= -box->x;
					box->theta= M_PI-box->theta;
					if (box->theta>M_PI)
						box->theta-= 2.0*M_PI;
				}
			while (box->x <0.0)
				box->x+=simulation->state->length_x;
			while (box->x >simulation->state->length_x)
				box->x-=simulation->state->length_x;
			if (simulation->parameters->particle_species==1)
				box->y+=*simulation->parameters->speed*sin(box->theta);
			else
				box->y+=*(simulation->parameters->speed+get_species(box->index, simulation->parameters))*sin(box->theta);
			if (simulation->parameters->boundary_y!=0) //reflecting boundary in y-direction
				if ((box->y<0.0) || (box->y>simulation->state->length_y))
				{
					box->y= -box->y;
					box->theta= -box->theta;
				}
			while (box->y <0.0)
				box->y+=simulation->state->length_y;
			while (box->y >simulation->state->length_y)
				box->y-=simulation->state->length_y;
			box=box->next_particle;
		}
	}
	VM_clean_boxes(simulation->parameters);
	if (simulation->parameters->particle_species==1)
	{
		*polar_x=*polar_x/simulation->state->particle_number;
		*polar_y=*polar_y/simulation->state->particle_number;
		*nematic_x=*nematic_x/simulation->state->particle_number;
		*nematic_y=*nematic_y/simulation->state->particle_number;
		*simulation->observables->polar_order=sqrt(*polar_x**polar_x+*polar_y**polar_y);
		*simulation->observables->nematic_order=sqrt(*nematic_x**nematic_x+*nematic_y**nematic_y);
	}
	else
		for (k=0; k<simulation->parameters->particle_species; k++)
		{
			*(polar_x+k)=*(polar_x+k)/ *(simulation->parameters->species_particle_number+k);
			*(polar_y+k)=*(polar_y+k)/ *(simulation->parameters->species_particle_number+k);
			*(nematic_x+k)=*(nematic_x+k)/ *(simulation->parameters->species_particle_number+k);
			*(nematic_y+k)=*(nematic_y+k)/ *(simulation->parameters->species_particle_number+k);
			*(simulation->observables->polar_order+k)=sqrt(*(polar_x+k)**(polar_x+k)+*(polar_y+k)**(polar_y+k));
			*(simulation->observables->nematic_order+k)=sqrt(*(nematic_x+k)**(nematic_x+k)+*(nematic_y+k)**(nematic_y+k));
		}
	free(polar_x);
	free(polar_y);
	free(nematic_x);
	free(nematic_y);
	return;
}

//performs time evolution of one step for nematic Vicsek model
void NVM_one_step(struct VM_simulation * simulation)
{
	gint32 i;
	struct VM_particle * box;
	gint32 k;
	gdouble * polar_x=malloc(simulation->parameters->particle_species*sizeof(gdouble));
	gdouble * polar_y=malloc(simulation->parameters->particle_species*sizeof(gdouble));
	gdouble * nematic_x=malloc(simulation->parameters->particle_species*sizeof(gdouble));
	gdouble * nematic_y=malloc(simulation->parameters->particle_species*sizeof(gdouble));
	for (k=0; k<simulation->parameters->particle_species; k++)
	{
		*(polar_x+k)=0.0;
		*(polar_y+k)=0.0;
		*(nematic_x+k)=0.0;
		*(nematic_y+k)=0.0;
	}
	VM_fill_boxes(simulation);
	for (i=0; i<simulation->parameters->box_number_x*simulation->parameters->box_number_y; i++)
	{
		box=*(simulation->parameters->box+i);
		while (box !=NULL)
		{
			//alignment interaction
			NVM_interaction_with_neighbors(box->index, simulation);
			box=box->next_particle;
		}
	}
	for (i=0; i<simulation->parameters->box_number_x*simulation->parameters->box_number_y; i++)
	{
		box=*(simulation->parameters->box+i);
		while (box !=NULL)
		{
			//calculate contributions to polar and nematic order parameters
			if (simulation->parameters->particle_species==1)
			{
				*polar_x+=cos(box->theta);
				*polar_y+=sin(box->theta);
				*nematic_x+=cos(2.0*box->theta);
				*nematic_y+=sin(2.0*box->theta);
			}
			else
			{
				*(polar_x+get_species(box->index, simulation->parameters))+=cos(box->theta);
				*(polar_y+get_species(box->index, simulation->parameters))+=sin(box->theta);
				*(nematic_x+get_species(box->index, simulation->parameters))+=cos(2.0*box->theta);
				*(nematic_y+get_species(box->index, simulation->parameters))+=sin(2.0*box->theta);
			}
			//set angle to value after interaction
			box->theta=box->theta_new;
			//streaming -> movement of the particles
			if (simulation->parameters->particle_species==1)
				box->x+=*simulation->parameters->speed*cos(box->theta);
			else
				box->x+=*(simulation->parameters->speed+get_species(box->index, simulation->parameters))*cos(box->theta);
			if (simulation->parameters->boundary_x!=0) //reflecting boundary in x-direction
				if ((box->x<0.0) || (box->x>simulation->state->length_x))
				{
					box->x= -box->x;
					box->theta= M_PI-box->theta;
					if (box->theta>M_PI)
						box->theta-= 2.0*M_PI;
				}
			while (box->x <0.0)
				box->x+=simulation->state->length_x;
			while (box->x >simulation->state->length_x)
				box->x-=simulation->state->length_x;
			if (simulation->parameters->particle_species==1)
				box->y+=*simulation->parameters->speed*sin(box->theta);
			else
				box->y+=*(simulation->parameters->speed+get_species(box->index, simulation->parameters))*sin(box->theta);
			if (simulation->parameters->boundary_y!=0) //reflecting boundary in y-direction
				if ((box->y<0.0) || (box->y>simulation->state->length_y))
				{
					box->y= -box->y;
					box->theta= -box->theta;
				}
			while (box->y <0.0)
				box->y+=simulation->state->length_y;
			while (box->y >simulation->state->length_y)
				box->y-=simulation->state->length_y;
			box=box->next_particle;
		}
	}
	VM_clean_boxes(simulation->parameters);
	if (simulation->parameters->particle_species==1)
	{
		*polar_x=*polar_x/simulation->state->particle_number;
		*polar_y=*polar_y/simulation->state->particle_number;
		*nematic_x=*nematic_x/simulation->state->particle_number;
		*nematic_y=*nematic_y/simulation->state->particle_number;
		*simulation->observables->polar_order=sqrt(*polar_x**polar_x+*polar_y**polar_y);
		*simulation->observables->nematic_order=sqrt(*nematic_x**nematic_x+*nematic_y**nematic_y);
	}
	else
		for (k=0; k<simulation->parameters->particle_species; k++)
		{
			*(polar_x+k)=*(polar_x+k)/ *(simulation->parameters->species_particle_number+k);
			*(polar_y+k)=*(polar_y+k)/ *(simulation->parameters->species_particle_number+k);
			*(nematic_x+k)=*(nematic_x+k)/ *(simulation->parameters->species_particle_number+k);
			*(nematic_y+k)=*(nematic_y+k)/ *(simulation->parameters->species_particle_number+k);
			*(simulation->observables->polar_order+k)=sqrt(*(polar_x+k)**(polar_x+k)+*(polar_y+k)**(polar_y+k));
			*(simulation->observables->nematic_order+k)=sqrt(*(nematic_x+k)**(nematic_x+k)+*(nematic_y+k)**(nematic_y+k));
		}
	free(polar_x);
	free(polar_y);
	free(nematic_x);
	free(nematic_y);
	return;
}

//performs time evolution of one step for metric free/topological Vicsek model
void mfVM_one_step(struct VM_simulation * simulation)
{
	gint32 i;
	struct VM_particle * box;
	gint32 k;
	gdouble * polar_x=malloc(simulation->parameters->particle_species*sizeof(gdouble));
	gdouble * polar_y=malloc(simulation->parameters->particle_species*sizeof(gdouble));
	gdouble * nematic_x=malloc(simulation->parameters->particle_species*sizeof(gdouble));
	gdouble * nematic_y=malloc(simulation->parameters->particle_species*sizeof(gdouble));
	for (k=0; k<simulation->parameters->particle_species; k++)
	{
		*(polar_x+k)=0.0;
		*(polar_y+k)=0.0;
		*(nematic_x+k)=0.0;
		*(nematic_y+k)=0.0;
	}
	mfVM_fill_boxes(simulation);
	for (i=0; i<simulation->parameters->mf_box_number_x*simulation->parameters->mf_box_number_y; i++)
	{
		box=*(simulation->parameters->mf_box+i);
		while (box !=NULL)
		{
			//alignment interaction
			mfVM_interaction_with_neighbors(box->index, simulation);
			box=box->next_particle;
		}
	}
	for (i=0; i<simulation->parameters->mf_box_number_x*simulation->parameters->mf_box_number_y; i++)
	{
		box=*(simulation->parameters->mf_box+i);
		while (box !=NULL)
		{
			//calculate contributions to polar and nematic order parameters
			if (simulation->parameters->particle_species==1)
			{
				*polar_x+=cos(box->theta);
				*polar_y+=sin(box->theta);
				*nematic_x+=cos(2.0*box->theta);
				*nematic_y+=sin(2.0*box->theta);
			}
			else
			{
				*(polar_x+get_species(box->index, simulation->parameters))+=cos(box->theta);
				*(polar_y+get_species(box->index, simulation->parameters))+=sin(box->theta);
				*(nematic_x+get_species(box->index, simulation->parameters))+=cos(2.0*box->theta);
				*(nematic_y+get_species(box->index, simulation->parameters))+=sin(2.0*box->theta);
			}
			//set angle to value after interaction
			box->theta=box->theta_new;
			//streaming -> movement of the particles
			if (simulation->parameters->particle_species==1)
				box->x+=*simulation->parameters->speed*cos(box->theta);
			else
				box->x+=*(simulation->parameters->speed+get_species(box->index, simulation->parameters))*cos(box->theta);
			if (simulation->parameters->boundary_x!=0) //reflecting boundary in x-direction
				if ((box->x<0.0) || (box->x>simulation->state->length_x))
				{
					box->x= -box->x;
					box->theta= M_PI-box->theta;
					if (box->theta>M_PI)
						box->theta-= 2.0*M_PI;
				}
			while (box->x <0.0)
				box->x+=simulation->state->length_x;
			while (box->x >simulation->state->length_x)
				box->x-=simulation->state->length_x;
			if (simulation->parameters->particle_species==1)
				box->y+=*simulation->parameters->speed*sin(box->theta);
			else
				box->y+=*(simulation->parameters->speed+get_species(box->index, simulation->parameters))*sin(box->theta);
			if (simulation->parameters->boundary_y!=0) //reflecting boundary in y-direction
				if ((box->y<0.0) || (box->y>simulation->state->length_y))
				{
					box->y= -box->y;
					box->theta= -box->theta;
				}
			while (box->y <0.0)
				box->y+=simulation->state->length_y;
			while (box->y >simulation->state->length_y)
				box->y-=simulation->state->length_y;
			box=box->next_particle;
		}
	}
	mfVM_clean_boxes(simulation->parameters);
	if (simulation->parameters->particle_species==1)
	{
		*polar_x=*polar_x/simulation->state->particle_number;
		*polar_y=*polar_y/simulation->state->particle_number;
		*nematic_x=*nematic_x/simulation->state->particle_number;
		*nematic_y=*nematic_y/simulation->state->particle_number;
		*simulation->observables->polar_order=sqrt(*polar_x**polar_x+*polar_y**polar_y);
		*simulation->observables->nematic_order=sqrt(*nematic_x**nematic_x+*nematic_y**nematic_y);
	}
	else
		for (k=0; k<simulation->parameters->particle_species; k++)
		{
			*(polar_x+k)=*(polar_x+k)/ *(simulation->parameters->species_particle_number+k);
			*(polar_y+k)=*(polar_y+k)/ *(simulation->parameters->species_particle_number+k);
			*(nematic_x+k)=*(nematic_x+k)/ *(simulation->parameters->species_particle_number+k);
			*(nematic_y+k)=*(nematic_y+k)/ *(simulation->parameters->species_particle_number+k);
			*(simulation->observables->polar_order+k)=sqrt(*(polar_x+k)**(polar_x+k)+*(polar_y+k)**(polar_y+k));
			*(simulation->observables->nematic_order+k)=sqrt(*(nematic_x+k)**(nematic_x+k)+*(nematic_y+k)**(nematic_y+k));
		}
	free(polar_x);
	free(polar_y);
	free(nematic_x);
	free(nematic_y);
	return;
}

//performs time evolution of one step for metric free/topological nematic Vicsek model
void mfNVM_one_step(struct VM_simulation * simulation)
{
	gint32 i;
	struct VM_particle * box;
	gint32 k;
	gdouble * polar_x=malloc(simulation->parameters->particle_species*sizeof(gdouble));
	gdouble * polar_y=malloc(simulation->parameters->particle_species*sizeof(gdouble));
	gdouble * nematic_x=malloc(simulation->parameters->particle_species*sizeof(gdouble));
	gdouble * nematic_y=malloc(simulation->parameters->particle_species*sizeof(gdouble));
	for (k=0; k<simulation->parameters->particle_species; k++)
	{
		*(polar_x+k)=0.0;
		*(polar_y+k)=0.0;
		*(nematic_x+k)=0.0;
		*(nematic_y+k)=0.0;
	}
	mfVM_fill_boxes(simulation);
	for (i=0; i<simulation->parameters->mf_box_number_x*simulation->parameters->mf_box_number_y; i++)
	{
		box=*(simulation->parameters->mf_box+i);
		while (box !=NULL)
		{
			//alignment interaction
			mfNVM_interaction_with_neighbors(box->index, simulation);
			box=box->next_particle;
		}
	}
	for (i=0; i<simulation->parameters->mf_box_number_x*simulation->parameters->mf_box_number_y; i++)
	{
		box=*(simulation->parameters->mf_box+i);
		while (box !=NULL)
		{
			//calculate contributions to polar and nematic order parameters
			if (simulation->parameters->particle_species==1)
			{
				*polar_x+=cos(box->theta);
				*polar_y+=sin(box->theta);
				*nematic_x+=cos(2.0*box->theta);
				*nematic_y+=sin(2.0*box->theta);
			}
			else
			{
				*(polar_x+get_species(box->index, simulation->parameters))+=cos(box->theta);
				*(polar_y+get_species(box->index, simulation->parameters))+=sin(box->theta);
				*(nematic_x+get_species(box->index, simulation->parameters))+=cos(2.0*box->theta);
				*(nematic_y+get_species(box->index, simulation->parameters))+=sin(2.0*box->theta);
			}
			//set angle to value after interaction
			box->theta=box->theta_new;
			//streaming -> movement of the particles
			if (simulation->parameters->particle_species==1)
				box->x+=*simulation->parameters->speed*cos(box->theta);
			else
				box->x+=*(simulation->parameters->speed+get_species(box->index, simulation->parameters))*cos(box->theta);
			if (simulation->parameters->boundary_x!=0) //reflecting boundary in x-direction
				if ((box->x<0.0) || (box->x>simulation->state->length_x))
				{
					box->x= -box->x;
					box->theta= M_PI-box->theta;
					if (box->theta>M_PI)
						box->theta-= 2.0*M_PI;
				}
			while (box->x <0.0)
				box->x+=simulation->state->length_x;
			while (box->x >simulation->state->length_x)
				box->x-=simulation->state->length_x;
			if (simulation->parameters->particle_species==1)
				box->y+=*simulation->parameters->speed*sin(box->theta);
			else
				box->y+=*(simulation->parameters->speed+get_species(box->index, simulation->parameters))*sin(box->theta);
			if (simulation->parameters->boundary_y!=0) //reflecting boundary in y-direction
				if ((box->y<0.0) || (box->y>simulation->state->length_y))
				{
					box->y= -box->y;
					box->theta= -box->theta;
				}
			while (box->y <0.0)
				box->y+=simulation->state->length_y;
			while (box->y >simulation->state->length_y)
				box->y-=simulation->state->length_y;
			box=box->next_particle;
		}
	}
	mfVM_clean_boxes(simulation->parameters);
	if (simulation->parameters->particle_species==1)
	{
		*polar_x=*polar_x/simulation->state->particle_number;
		*polar_y=*polar_y/simulation->state->particle_number;
		*nematic_x=*nematic_x/simulation->state->particle_number;
		*nematic_y=*nematic_y/simulation->state->particle_number;
		*simulation->observables->polar_order=sqrt(*polar_x**polar_x+*polar_y**polar_y);
		*simulation->observables->nematic_order=sqrt(*nematic_x**nematic_x+*nematic_y**nematic_y);
	}
	else
		for (k=0; k<simulation->parameters->particle_species; k++)
		{
			*(polar_x+k)=*(polar_x+k)/ *(simulation->parameters->species_particle_number+k);
			*(polar_y+k)=*(polar_y+k)/ *(simulation->parameters->species_particle_number+k);
			*(nematic_x+k)=*(nematic_x+k)/ *(simulation->parameters->species_particle_number+k);
			*(nematic_y+k)=*(nematic_y+k)/ *(simulation->parameters->species_particle_number+k);
			*(simulation->observables->polar_order+k)=sqrt(*(polar_x+k)**(polar_x+k)+*(polar_y+k)**(polar_y+k));
			*(simulation->observables->nematic_order+k)=sqrt(*(nematic_x+k)**(nematic_x+k)+*(nematic_y+k)**(nematic_y+k));
		}
	free(polar_x);
	free(polar_y);
	free(nematic_x);
	free(nematic_y);
	return;
}

//performs time evolution of one step for additive Langevin model
void additiveL_one_step(struct VM_simulation * simulation)
{
	gint32 i;
	struct VM_particle * box;
	gint32 k;
	gdouble * polar_x=malloc(simulation->parameters->particle_species*sizeof(gdouble));
	gdouble * polar_y=malloc(simulation->parameters->particle_species*sizeof(gdouble));
	gdouble * nematic_x=malloc(simulation->parameters->particle_species*sizeof(gdouble));
	gdouble * nematic_y=malloc(simulation->parameters->particle_species*sizeof(gdouble));
	for (k=0; k<simulation->parameters->particle_species; k++)
	{
		*(polar_x+k)=0.0;
		*(polar_y+k)=0.0;
		*(nematic_x+k)=0.0;
		*(nematic_y+k)=0.0;
	}
	VM_fill_boxes(simulation);
	for (i=0; i<simulation->parameters->box_number_x*simulation->parameters->box_number_y; i++)
	{
		box=*(simulation->parameters->box+i);
		while (box !=NULL)
		{
			//alignment interaction
			additiveL_interaction_with_neighbors(box->index, simulation);
			box=box->next_particle;
		}
	}
	for (i=0; i<simulation->parameters->box_number_x*simulation->parameters->box_number_y; i++)
	{
		box=*(simulation->parameters->box+i);
		while (box !=NULL)
		{
			//streaming -> movement of the particles
			if (simulation->parameters->particle_species==1)
				box->x+=simulation->parameters->delta_t**simulation->parameters->speed*cos(box->theta);
			else
				box->x+=simulation->parameters->delta_t**(simulation->parameters->speed+get_species(box->index, simulation->parameters))*cos(box->theta);
			if (simulation->parameters->particle_species==1)
				box->y+=simulation->parameters->delta_t**simulation->parameters->speed*sin(box->theta);
			else
				box->y+=simulation->parameters->delta_t**(simulation->parameters->speed+get_species(box->index, simulation->parameters))*sin(box->theta);
			//apply boundary conditions
			if (simulation->parameters->boundary_x!=0) //reflecting boundary in x-direction
				if ((box->x<0.0) || (box->x>simulation->state->length_x))
				{
					box->x= -box->x;
					//mistake in implementation for reflecting bc corrected in version 1.2
					box->theta_new= M_PI-box->theta_new;
					if (box->theta_new>M_PI)
						box->theta_new-= 2.0*M_PI;
				}
			while (box->x <0.0)
				box->x+=simulation->state->length_x;
			while (box->x >simulation->state->length_x)
				box->x-=simulation->state->length_x;
			if (simulation->parameters->boundary_y!=0) //reflecting boundary in y-direction
				if ((box->y<0.0) || (box->y>simulation->state->length_y))
				{
					box->y= -box->y;
					//mistake in implementation for reflecting bc corrected in version 1.2
					box->theta_new= -box->theta_new;
				}
			while (box->y <0.0)
				box->y+=simulation->state->length_y;
			while (box->y >simulation->state->length_y)
				box->y-=simulation->state->length_y;
			//calculate contribution to polar/nematic order parameter
			if (simulation->parameters->particle_species==1)
			{
				*polar_x+=cos(box->theta);
				*polar_y+=sin(box->theta);
				*nematic_x+=cos(2.0*box->theta);
				*nematic_y+=sin(2.0*box->theta);
			}
			else
			{
				*(polar_x+get_species(box->index, simulation->parameters))+=cos(box->theta);
				*(polar_y+get_species(box->index, simulation->parameters))+=sin(box->theta);
				*(nematic_x+get_species(box->index, simulation->parameters))+=cos(2.0*box->theta);
				*(nematic_y+get_species(box->index, simulation->parameters))+=sin(2.0*box->theta);
			}
			//update angles
			box->theta=box->theta_new;
			box=box->next_particle;
		}
	}
	VM_clean_boxes(simulation->parameters);
	if (simulation->parameters->particle_species==1)
	{
		*polar_x=*polar_x/simulation->state->particle_number;
		*polar_y=*polar_y/simulation->state->particle_number;
		*nematic_x=*nematic_x/simulation->state->particle_number;
		*nematic_y=*nematic_y/simulation->state->particle_number;
		*simulation->observables->polar_order=sqrt(*polar_x**polar_x+*polar_y**polar_y);
		*simulation->observables->nematic_order=sqrt(*nematic_x**nematic_x+*nematic_y**nematic_y);
	}
	else
		for (k=0; k<simulation->parameters->particle_species; k++)
		{
			*(polar_x+k)=*(polar_x+k)/ *(simulation->parameters->species_particle_number+k);
			*(polar_y+k)=*(polar_y+k)/ *(simulation->parameters->species_particle_number+k);
			*(nematic_x+k)=*(nematic_x+k)/ *(simulation->parameters->species_particle_number+k);
			*(nematic_y+k)=*(nematic_y+k)/ *(simulation->parameters->species_particle_number+k);
			*(simulation->observables->polar_order+k)=sqrt(*(polar_x+k)**(polar_x+k)+*(polar_y+k)**(polar_y+k));
			*(simulation->observables->nematic_order+k)=sqrt(*(nematic_x+k)**(nematic_x+k)+*(nematic_y+k)**(nematic_y+k));
		}
	free(polar_x);
	free(polar_y);
	free(nematic_x);
	free(nematic_y);
	return;
}

//performs time evolution of one step for nonadditive Langevin model
void nonadditiveL_one_step(struct VM_simulation * simulation)
{
	gint32 i;
	struct VM_particle * box;
	gint32 k;
	gdouble * polar_x=malloc(simulation->parameters->particle_species*sizeof(gdouble));
	gdouble * polar_y=malloc(simulation->parameters->particle_species*sizeof(gdouble));
	gdouble * nematic_x=malloc(simulation->parameters->particle_species*sizeof(gdouble));
	gdouble * nematic_y=malloc(simulation->parameters->particle_species*sizeof(gdouble));
	for (k=0; k<simulation->parameters->particle_species; k++)
	{
		*(polar_x+k)=0.0;
		*(polar_y+k)=0.0;
		*(nematic_x+k)=0.0;
		*(nematic_y+k)=0.0;
	}
	VM_fill_boxes(simulation);
	for (i=0; i<simulation->parameters->box_number_x*simulation->parameters->box_number_y; i++)
	{
		box=*(simulation->parameters->box+i);
		while (box !=NULL)
		{
			//alignment interaction
			nonadditiveL_interaction_with_neighbors(box->index, simulation);
			box=box->next_particle;
		}
	}
	for (i=0; i<simulation->parameters->box_number_x*simulation->parameters->box_number_y; i++)
	{
		box=*(simulation->parameters->box+i);
		while (box !=NULL)
		{
			//streaming -> movement of the particles
			if (simulation->parameters->particle_species==1)
				box->x+=simulation->parameters->delta_t**simulation->parameters->speed*cos(box->theta);
			else
				box->x+=simulation->parameters->delta_t**(simulation->parameters->speed+get_species(box->index, simulation->parameters))*cos(box->theta);
			if (simulation->parameters->particle_species==1)
				box->y+=simulation->parameters->delta_t**simulation->parameters->speed*sin(box->theta);
			else
				box->y+=simulation->parameters->delta_t**(simulation->parameters->speed+get_species(box->index, simulation->parameters))*sin(box->theta);
			//apply boundary conditions
			if (simulation->parameters->boundary_x!=0) //reflecting boundary in x-direction
				if ((box->x<0.0) || (box->x>simulation->state->length_x))
				{
					box->x= -box->x;
					//mistake in implementation for reflecting bc corrected in version 1.2
					box->theta_new= M_PI-box->theta_new;
					if (box->theta_new>M_PI)
						box->theta_new-= 2.0*M_PI;
				}
			while (box->x <0.0)
				box->x+=simulation->state->length_x;
			while (box->x >simulation->state->length_x)
				box->x-=simulation->state->length_x;
			if (simulation->parameters->boundary_y!=0) //reflecting boundary in y-direction
				if ((box->y<0.0) || (box->y>simulation->state->length_y))
				{
					box->y= -box->y;
					//mistake in implementation for reflecting bc corrected in version 1.2
					box->theta_new= -box->theta_new;
				}
			while (box->y <0.0)
				box->y+=simulation->state->length_y;
			while (box->y >simulation->state->length_y)
				box->y-=simulation->state->length_y;
			//calculate contribution to polar/nematic order parameter
			if (simulation->parameters->particle_species==1)
			{
				*polar_x+=cos(box->theta);
				*polar_y+=sin(box->theta);
				*nematic_x+=cos(2.0*box->theta);
				*nematic_y+=sin(2.0*box->theta);
			}
			else
			{
				*(polar_x+get_species(box->index, simulation->parameters))+=cos(box->theta);
				*(polar_y+get_species(box->index, simulation->parameters))+=sin(box->theta);
				*(nematic_x+get_species(box->index, simulation->parameters))+=cos(2.0*box->theta);
				*(nematic_y+get_species(box->index, simulation->parameters))+=sin(2.0*box->theta);
			}
			//update angles
			box->theta=box->theta_new;
			box=box->next_particle;
		}
	}
	VM_clean_boxes(simulation->parameters);
	if (simulation->parameters->particle_species==1)
	{
		*polar_x=*polar_x/simulation->state->particle_number;
		*polar_y=*polar_y/simulation->state->particle_number;
		*nematic_x=*nematic_x/simulation->state->particle_number;
		*nematic_y=*nematic_y/simulation->state->particle_number;
		*simulation->observables->polar_order=sqrt(*polar_x**polar_x+*polar_y**polar_y);
		*simulation->observables->nematic_order=sqrt(*nematic_x**nematic_x+*nematic_y**nematic_y);
	}
	else
		for (k=0; k<simulation->parameters->particle_species; k++)
		{
			*(polar_x+k)=*(polar_x+k)/ *(simulation->parameters->species_particle_number+k);
			*(polar_y+k)=*(polar_y+k)/ *(simulation->parameters->species_particle_number+k);
			*(nematic_x+k)=*(nematic_x+k)/ *(simulation->parameters->species_particle_number+k);
			*(nematic_y+k)=*(nematic_y+k)/ *(simulation->parameters->species_particle_number+k);
			*(simulation->observables->polar_order+k)=sqrt(*(polar_x+k)**(polar_x+k)+*(polar_y+k)**(polar_y+k));
			*(simulation->observables->nematic_order+k)=sqrt(*(nematic_x+k)**(nematic_x+k)+*(nematic_y+k)**(nematic_y+k));
		}
	free(polar_x);
	free(polar_y);
	free(nematic_x);
	free(nematic_y);
	return;
}

//performs time evolution of one step for metric free/topological Langevin model
void mfL_one_step(struct VM_simulation * simulation)
{
	gint32 i;
	struct VM_particle * box;
	gint32 k;
	gdouble * polar_x=malloc(simulation->parameters->particle_species*sizeof(gdouble));
	gdouble * polar_y=malloc(simulation->parameters->particle_species*sizeof(gdouble));
	gdouble * nematic_x=malloc(simulation->parameters->particle_species*sizeof(gdouble));
	gdouble * nematic_y=malloc(simulation->parameters->particle_species*sizeof(gdouble));
	for (k=0; k<simulation->parameters->particle_species; k++)
	{
		*(polar_x+k)=0.0;
		*(polar_y+k)=0.0;
		*(nematic_x+k)=0.0;
		*(nematic_y+k)=0.0;
	}
	mfVM_fill_boxes(simulation);
	for (i=0; i<simulation->parameters->mf_box_number_x*simulation->parameters->mf_box_number_y; i++)
	{
		box=*(simulation->parameters->mf_box+i);
		while (box !=NULL)
		{
			//alignment interaction
			mfL_interaction_with_neighbors(box->index, simulation);
			box=box->next_particle;
		}
	}
	for (i=0; i<simulation->parameters->mf_box_number_x*simulation->parameters->mf_box_number_y; i++)
	{
		box=*(simulation->parameters->mf_box+i);
		while (box !=NULL)
		{
			//streaming -> movement of the particles
			if (simulation->parameters->particle_species==1)
				box->x+=simulation->parameters->delta_t**simulation->parameters->speed*cos(box->theta);
			else
				box->x+=simulation->parameters->delta_t**(simulation->parameters->speed+get_species(box->index, simulation->parameters))*cos(box->theta);
			if (simulation->parameters->particle_species==1)
				box->y+=simulation->parameters->delta_t**simulation->parameters->speed*sin(box->theta);
			else
				box->y+=simulation->parameters->delta_t**(simulation->parameters->speed+get_species(box->index, simulation->parameters))*sin(box->theta);
			//apply boundary conditions
			if (simulation->parameters->boundary_x!=0) //reflecting boundary in x-direction
				if ((box->x<0.0) || (box->x>simulation->state->length_x))
				{
					box->x= -box->x;
					//mistake in implementation for reflecting bc corrected in version 1.2
					box->theta_new= M_PI-box->theta_new;
					if (box->theta_new>M_PI)
						box->theta_new-= 2.0*M_PI;
				}
			while (box->x <0.0)
				box->x+=simulation->state->length_x;
			while (box->x >simulation->state->length_x)
				box->x-=simulation->state->length_x;
			if (simulation->parameters->boundary_y!=0) //reflecting boundary in y-direction
				if ((box->y<0.0) || (box->y>simulation->state->length_y))
				{
					box->y= -box->y;
					//mistake in implementation for reflecting bc corrected in version 1.2
					box->theta_new= -box->theta_new;
				}
			while (box->y <0.0)
				box->y+=simulation->state->length_y;
			while (box->y >simulation->state->length_y)
				box->y-=simulation->state->length_y;
			//calculate contribution to polar/nematic order parameter
			if (simulation->parameters->particle_species==1)
			{
				*polar_x+=cos(box->theta);
				*polar_y+=sin(box->theta);
				*nematic_x+=cos(2.0*box->theta);
				*nematic_y+=sin(2.0*box->theta);
			}
			else
			{
				*(polar_x+get_species(box->index, simulation->parameters))+=cos(box->theta);
				*(polar_y+get_species(box->index, simulation->parameters))+=sin(box->theta);
				*(nematic_x+get_species(box->index, simulation->parameters))+=cos(2.0*box->theta);
				*(nematic_y+get_species(box->index, simulation->parameters))+=sin(2.0*box->theta);
			}
			//update angles
			box->theta=box->theta_new;
			box=box->next_particle;
		}
	}
	mfVM_clean_boxes(simulation->parameters);
	if (simulation->parameters->particle_species==1)
	{
		*polar_x=*polar_x/simulation->state->particle_number;
		*polar_y=*polar_y/simulation->state->particle_number;
		*nematic_x=*nematic_x/simulation->state->particle_number;
		*nematic_y=*nematic_y/simulation->state->particle_number;
		*simulation->observables->polar_order=sqrt(*polar_x**polar_x+*polar_y**polar_y);
		*simulation->observables->nematic_order=sqrt(*nematic_x**nematic_x+*nematic_y**nematic_y);
	}
	else
		for (k=0; k<simulation->parameters->particle_species; k++)
		{
			*(polar_x+k)=*(polar_x+k)/ *(simulation->parameters->species_particle_number+k);
			*(polar_y+k)=*(polar_y+k)/ *(simulation->parameters->species_particle_number+k);
			*(nematic_x+k)=*(nematic_x+k)/ *(simulation->parameters->species_particle_number+k);
			*(nematic_y+k)=*(nematic_y+k)/ *(simulation->parameters->species_particle_number+k);
			*(simulation->observables->polar_order+k)=sqrt(*(polar_x+k)**(polar_x+k)+*(polar_y+k)**(polar_y+k));
			*(simulation->observables->nematic_order+k)=sqrt(*(nematic_x+k)**(nematic_x+k)+*(nematic_y+k)**(nematic_y+k));
		}
	free(polar_x);
	free(polar_y);
	free(nematic_x);
	free(nematic_y);
	return;
}

//performs time evolution of one step for Vicsek model and measures observables
void VM_one_step_measurement(struct VM_simulation * simulation)
{
	gint32 i;
	gint32 neighbor_n;
	struct VM_particle * box;
	gint32 k;
	gdouble * polar_x=malloc(simulation->parameters->particle_species*sizeof(gdouble));
	gdouble * polar_y=malloc(simulation->parameters->particle_species*sizeof(gdouble));
	gdouble * nematic_x=malloc(simulation->parameters->particle_species*sizeof(gdouble));
	gdouble * nematic_y=malloc(simulation->parameters->particle_species*sizeof(gdouble));
	//new in version 1.1 calculate neighbor moments from ensemble average
	gdouble * neighbor_number=malloc(simulation->parameters->particle_species*sizeof(gdouble));
	for (i=0; i<simulation->parameters->particle_species; i++)
		*(neighbor_number+i)=0.1;
	for (k=0; k<simulation->parameters->particle_species; k++)
	{
		*(polar_x+k)=0.0;
		*(polar_y+k)=0.0;
		*(nematic_x+k)=0.0;
		*(nematic_y+k)=0.0;
	}
	VM_fill_boxes(simulation);
	for (i=0; i<simulation->parameters->box_number_x*simulation->parameters->box_number_y; i++)
	{
		box=*(simulation->parameters->box+i);
		while (box !=NULL)
		{
			//alignment interaction
			neighbor_n=VM_interaction_with_neighbors(box->index, simulation);
			//analyze number of neighbors statistics
			if (simulation->parameters->particle_species==1)
			{
				//changed in version 1.1
				*neighbor_number+=1.0*neighbor_n;
				if (neighbor_n<simulation->observables->binnum_neighbors)
					*(*simulation->observables->hist_neighbors+neighbor_n)+=1;
			}
			else
			{
				//changed in version 1.1
				*(neighbor_number+get_species(box->index, simulation->parameters))+=1.0*neighbor_n;
				if (neighbor_n<simulation->observables->binnum_neighbors)
					*(*(simulation->observables->hist_neighbors+get_species(box->index, simulation->parameters))+neighbor_n)+=1;
			}
			box=box->next_particle;
		}
	}
	//normalize ensemble averaged neighbor number
	//changed in version 1.1
	if (simulation->parameters->particle_species==1)
		*neighbor_number=*neighbor_number/simulation->state->particle_number;
	else
		for (i=0; i<simulation->parameters->particle_species; i++)
			*(neighbor_number+i)=*(neighbor_number+i)/ *(simulation->parameters->species_particle_number+i);
	// update moments of ensemble averaged neighbor number
	//changed in version 1.1
	if (simulation->parameters->particle_species==1)
	{
		*simulation->observables->neighbors_moment1+=*(neighbor_number);
		*simulation->observables->neighbors_moment2+=*(neighbor_number)**(neighbor_number);
		*simulation->observables->neighbors_moment3+=*(neighbor_number)**(neighbor_number)**(neighbor_number);
		*simulation->observables->neighbors_moment4+=*(neighbor_number)**(neighbor_number)**(neighbor_number)**(neighbor_number);
	}
	else
		for (i=0; i<simulation->parameters->particle_species; i++)
		{
			*(simulation->observables->neighbors_moment1+i)+=*(neighbor_number+i);
			*(simulation->observables->neighbors_moment2+i)+=*(neighbor_number+i)**(neighbor_number+i);
			*(simulation->observables->neighbors_moment3+i)+=*(neighbor_number+i)**(neighbor_number+i)**(neighbor_number+i);
			*(simulation->observables->neighbors_moment4+i)+=*(neighbor_number+i)**(neighbor_number+i)**(neighbor_number+i)**(neighbor_number+i);
		}
	for (i=0; i<simulation->parameters->box_number_x*simulation->parameters->box_number_y; i++)
	{
		box=*(simulation->parameters->box+i);
		while (box !=NULL)
		{
			//calculate contributions to polar and nematic order parameters
			if (simulation->parameters->particle_species==1)
			{
				*polar_x+=cos(box->theta);
				*polar_y+=sin(box->theta);
				*nematic_x+=cos(2.0*box->theta);
				*nematic_y+=sin(2.0*box->theta);
				*(*simulation->observables->hist_theta+((gint32)((M_PI+box->theta)/(2.0*M_PI)*simulation->observables->binnum_theta)))+=1;
			}
			else
			{
				*(polar_x+get_species(box->index, simulation->parameters))+=cos(box->theta);
				*(polar_y+get_species(box->index, simulation->parameters))+=sin(box->theta);
				*(nematic_x+get_species(box->index, simulation->parameters))+=cos(2.0*box->theta);
				*(nematic_y+get_species(box->index, simulation->parameters))+=sin(2.0*box->theta);
				*(*(simulation->observables->hist_theta+get_species(box->index, simulation->parameters))+((gint32)((M_PI+box->theta)/(2.0*M_PI)*simulation->observables->binnum_theta)))+=1;
			}
			//set angle to value after interaction
			box->theta=box->theta_new;
			//streaming -> movement of the particles
			if (simulation->parameters->particle_species==1)
				box->x+=*simulation->parameters->speed*cos(box->theta);
			else
				box->x+=*(simulation->parameters->speed+get_species(box->index, simulation->parameters))*cos(box->theta);
			if (simulation->parameters->boundary_x!=0) //reflecting boundary in x-direction
				if ((box->x<0.0) || (box->x>simulation->state->length_x))
				{
					box->x= -box->x;
					box->theta=M_PI-box->theta;
					if (box->theta>M_PI)
						box->theta-= 2.0*M_PI;
				}
			while (box->x <0.0)
				box->x+=simulation->state->length_x;
			while (box->x >simulation->state->length_x)
				box->x-=simulation->state->length_x;
			if (simulation->parameters->particle_species==1)
				box->y+=*simulation->parameters->speed*sin(box->theta);
			else
				box->y+=*(simulation->parameters->speed+get_species(box->index, simulation->parameters))*sin(box->theta);
			if (simulation->parameters->boundary_y!=0) //reflecting boundary in y-direction
				if ((box->y<0.0) || (box->y>simulation->state->length_y))
				{
					box->y= -box->y;
					box->theta=-box->theta;
				}
			while (box->y <0.0)
				box->y+=simulation->state->length_y;
			while (box->y >simulation->state->length_y)
				box->y-=simulation->state->length_y;
			//go to next particle in box
			box=box->next_particle;
		}
	}
	VM_clean_boxes(simulation->parameters);
	if (simulation->parameters->particle_species==1)
	{
		*polar_x=*polar_x/simulation->state->particle_number;
		*polar_y=*polar_y/simulation->state->particle_number;
		*nematic_x=*nematic_x/simulation->state->particle_number;
		*nematic_y=*nematic_y/simulation->state->particle_number;
		*simulation->observables->polar_order=sqrt(*polar_x**polar_x+*polar_y**polar_y);
		*simulation->observables->nematic_order=sqrt(*nematic_x**nematic_x+*nematic_y**nematic_y);
		*simulation->observables->polar_order_moment1+=*simulation->observables->polar_order;
		*simulation->observables->polar_order_moment2+=*simulation->observables->polar_order**simulation->observables->polar_order;
		*simulation->observables->polar_order_moment3+=*simulation->observables->polar_order**simulation->observables->polar_order**simulation->observables->polar_order;
		*simulation->observables->polar_order_moment4+=*simulation->observables->polar_order**simulation->observables->polar_order**simulation->observables->polar_order**simulation->observables->polar_order;
		*simulation->observables->nematic_order_moment1+=*simulation->observables->nematic_order;
		*simulation->observables->nematic_order_moment2+=*simulation->observables->nematic_order**simulation->observables->nematic_order;
		*simulation->observables->nematic_order_moment3+=*simulation->observables->nematic_order**simulation->observables->nematic_order**simulation->observables->nematic_order;
		*simulation->observables->nematic_order_moment4+=*simulation->observables->nematic_order**simulation->observables->nematic_order**simulation->observables->nematic_order**simulation->observables->nematic_order;
	}
	else
		for (k=0; k<simulation->parameters->particle_species; k++)
		{
			*(polar_x+k)=*(polar_x+k)/ *(simulation->parameters->species_particle_number+k);
			*(polar_y+k)=*(polar_y+k)/ *(simulation->parameters->species_particle_number+k);
			*(nematic_x+k)=*(nematic_x+k)/ *(simulation->parameters->species_particle_number+k);
			*(nematic_y+k)=*(nematic_y+k)/ *(simulation->parameters->species_particle_number+k);
			*(simulation->observables->polar_order+k)=sqrt(*(polar_x+k)**(polar_x+k)+*(polar_y+k)**(polar_y+k));
			*(simulation->observables->nematic_order+k)=sqrt(*(nematic_x+k)**(nematic_x+k)+*(nematic_y+k)**(nematic_y+k));
			*(simulation->observables->polar_order_moment1+k)+=*(simulation->observables->polar_order+k);
			*(simulation->observables->polar_order_moment2+k)+=*(simulation->observables->polar_order+k)**(simulation->observables->polar_order+k);
			*(simulation->observables->polar_order_moment3+k)+=*(simulation->observables->polar_order+k)**(simulation->observables->polar_order+k)**(simulation->observables->polar_order+k);
			*(simulation->observables->polar_order_moment4+k)+=*(simulation->observables->polar_order+k)**(simulation->observables->polar_order+k)**(simulation->observables->polar_order+k)**(simulation->observables->polar_order+k);
			*(simulation->observables->nematic_order_moment1+k)+=*(simulation->observables->nematic_order+k);
			*(simulation->observables->nematic_order_moment2+k)+=*(simulation->observables->nematic_order+k)**(simulation->observables->nematic_order+k);
			*(simulation->observables->nematic_order_moment3+k)+=*(simulation->observables->nematic_order+k)**(simulation->observables->nematic_order+k)**(simulation->observables->nematic_order+k);
			*(simulation->observables->nematic_order_moment4+k)+=*(simulation->observables->nematic_order+k)**(simulation->observables->nematic_order+k)**(simulation->observables->nematic_order+k)**(simulation->observables->nematic_order+k);
		}
	simulation->observables->measurement_steps+=1;
	free(polar_x);
	free(polar_y);
	free(nematic_x);
	free(nematic_y);
	free(neighbor_number);
	return;
}

//performs time evolution of one step for nematic Vicsek model and measures observables
void NVM_one_step_measurement(struct VM_simulation * simulation)
{
	gint32 i;
	gint32 neighbor_n;
	struct VM_particle * box;
	gint32 k;
	gdouble * polar_x=malloc(simulation->parameters->particle_species*sizeof(gdouble));
	gdouble * polar_y=malloc(simulation->parameters->particle_species*sizeof(gdouble));
	gdouble * nematic_x=malloc(simulation->parameters->particle_species*sizeof(gdouble));
	gdouble * nematic_y=malloc(simulation->parameters->particle_species*sizeof(gdouble));
	//new in version 1.1 calculate neighbor moments from ensemble average
	gdouble * neighbor_number=malloc(simulation->parameters->particle_species*sizeof(gdouble));
	for (i=0; i<simulation->parameters->particle_species; i++)
		*(neighbor_number+i)=0.1;
	for (k=0; k<simulation->parameters->particle_species; k++)
	{
		*(polar_x+k)=0.0;
		*(polar_y+k)=0.0;
		*(nematic_x+k)=0.0;
		*(nematic_y+k)=0.0;
	}
	VM_fill_boxes(simulation);
	for (i=0; i<simulation->parameters->box_number_x*simulation->parameters->box_number_y; i++)
	{
		box=*(simulation->parameters->box+i);
		while (box !=NULL)
		{
			//alignment interaction
			neighbor_n=NVM_interaction_with_neighbors(box->index, simulation);
			//analyze number of neighbors statistics
			if (simulation->parameters->particle_species==1)
			{
				//changed in version 1.1
				*neighbor_number+=1.0*neighbor_n;
				if (neighbor_n<simulation->observables->binnum_neighbors)
					*(*simulation->observables->hist_neighbors+neighbor_n)+=1;
			}
			else
			{
				//changed in version 1.1
				*(neighbor_number+get_species(box->index, simulation->parameters))+=1.0*neighbor_n;
				if (neighbor_n<simulation->observables->binnum_neighbors)
					*(*(simulation->observables->hist_neighbors+get_species(box->index, simulation->parameters))+neighbor_n)+=1;
			}
			box=box->next_particle;
		}
	}
	//normalize ensemble averaged neighbor number
	//changed in version 1.1
	if (simulation->parameters->particle_species==1)
		*neighbor_number=*neighbor_number/simulation->state->particle_number;
	else
		for (i=0; i<simulation->parameters->particle_species; i++)
			*(neighbor_number+i)=*(neighbor_number+i)/ *(simulation->parameters->species_particle_number+i);
	// update moments of ensemble averaged neighbor number
	//changed in version 1.1
	if (simulation->parameters->particle_species==1)
	{
		*simulation->observables->neighbors_moment1+=*(neighbor_number);
		*simulation->observables->neighbors_moment2+=*(neighbor_number)**(neighbor_number);
		*simulation->observables->neighbors_moment3+=*(neighbor_number)**(neighbor_number)**(neighbor_number);
		*simulation->observables->neighbors_moment4+=*(neighbor_number)**(neighbor_number)**(neighbor_number)**(neighbor_number);
	}
	else
		for (i=0; i<simulation->parameters->particle_species; i++)
		{
			*(simulation->observables->neighbors_moment1+i)+=*(neighbor_number+i);
			*(simulation->observables->neighbors_moment2+i)+=*(neighbor_number+i)**(neighbor_number+i);
			*(simulation->observables->neighbors_moment3+i)+=*(neighbor_number+i)**(neighbor_number+i)**(neighbor_number+i);
			*(simulation->observables->neighbors_moment4+i)+=*(neighbor_number+i)**(neighbor_number+i)**(neighbor_number+i)**(neighbor_number+i);
		}
	for (i=0; i<simulation->parameters->box_number_x*simulation->parameters->box_number_y; i++)
	{
		box=*(simulation->parameters->box+i);
		while (box !=NULL)
		{
			//calculate contributions to polar and nematic order parameters
			if (simulation->parameters->particle_species==1)
			{
				*polar_x+=cos(box->theta);
				*polar_y+=sin(box->theta);
				*nematic_x+=cos(2.0*box->theta);
				*nematic_y+=sin(2.0*box->theta);
				*(*simulation->observables->hist_theta+((gint32)((M_PI+box->theta)/(2.0*M_PI)*simulation->observables->binnum_theta)))+=1;
			}
			else
			{
				*(polar_x+get_species(box->index, simulation->parameters))+=cos(box->theta);
				*(polar_y+get_species(box->index, simulation->parameters))+=sin(box->theta);
				*(nematic_x+get_species(box->index, simulation->parameters))+=cos(2.0*box->theta);
				*(nematic_y+get_species(box->index, simulation->parameters))+=sin(2.0*box->theta);
				*(*(simulation->observables->hist_theta+get_species(box->index, simulation->parameters))+((gint32)((M_PI+box->theta)/(2.0*M_PI)*simulation->observables->binnum_theta)))+=1;
			}
			//set angle to value after interaction
			box->theta=box->theta_new;
			//streaming -> movement of the particles
			if (simulation->parameters->particle_species==1)
				box->x+=*simulation->parameters->speed*cos(box->theta);
			else
				box->x+=*(simulation->parameters->speed+get_species(box->index, simulation->parameters))*cos(box->theta);
			if (simulation->parameters->boundary_x!=0) //reflecting boundary in x-direction
				if ((box->x<0.0) || (box->x>simulation->state->length_x))
				{
					box->x= -box->x;
					box->theta= M_PI-box->theta;
					if (box->theta>M_PI)
						box->theta-= 2.0*M_PI;
				}
			while (box->x <0.0)
				box->x+=simulation->state->length_x;
			while (box->x >simulation->state->length_x)
				box->x-=simulation->state->length_x;
			if (simulation->parameters->particle_species==1)
				box->y+=*simulation->parameters->speed*sin(box->theta);
			else
				box->y+=*(simulation->parameters->speed+get_species(box->index, simulation->parameters))*sin(box->theta);
			if (simulation->parameters->boundary_y!=0) //reflecting boundary in y-direction
				if ((box->y<0.0) || (box->y>simulation->state->length_y))
				{
					box->y= -box->y;
					box->theta= -box->theta;
				}
			while (box->y <0.0)
				box->y+=simulation->state->length_y;
			while (box->y >simulation->state->length_y)
				box->y-=simulation->state->length_y;
			box=box->next_particle;
		}
	}
	VM_clean_boxes(simulation->parameters);
	if (simulation->parameters->particle_species==1)
	{
		*polar_x=*polar_x/simulation->state->particle_number;
		*polar_y=*polar_y/simulation->state->particle_number;
		*nematic_x=*nematic_x/simulation->state->particle_number;
		*nematic_y=*nematic_y/simulation->state->particle_number;
		*simulation->observables->polar_order=sqrt(*polar_x**polar_x+*polar_y**polar_y);
		*simulation->observables->nematic_order=sqrt(*nematic_x**nematic_x+*nematic_y**nematic_y);
		*simulation->observables->polar_order_moment1+=*simulation->observables->polar_order;
		*simulation->observables->polar_order_moment2+=*simulation->observables->polar_order**simulation->observables->polar_order;
		*simulation->observables->polar_order_moment3+=*simulation->observables->polar_order**simulation->observables->polar_order**simulation->observables->polar_order;
		*simulation->observables->polar_order_moment4+=*simulation->observables->polar_order**simulation->observables->polar_order**simulation->observables->polar_order**simulation->observables->polar_order;
		*simulation->observables->nematic_order_moment1+=*simulation->observables->nematic_order;
		*simulation->observables->nematic_order_moment2+=*simulation->observables->nematic_order**simulation->observables->nematic_order;
		*simulation->observables->nematic_order_moment3+=*simulation->observables->nematic_order**simulation->observables->nematic_order**simulation->observables->nematic_order;
		*simulation->observables->nematic_order_moment4+=*simulation->observables->nematic_order**simulation->observables->nematic_order**simulation->observables->nematic_order**simulation->observables->nematic_order;
	}
	else
		for (k=0; k<simulation->parameters->particle_species; k++)
		{
			*(polar_x+k)=*(polar_x+k)/ *(simulation->parameters->species_particle_number+k);
			*(polar_y+k)=*(polar_y+k)/ *(simulation->parameters->species_particle_number+k);
			*(nematic_x+k)=*(nematic_x+k)/ *(simulation->parameters->species_particle_number+k);
			*(nematic_y+k)=*(nematic_y+k)/ *(simulation->parameters->species_particle_number+k);
			*(simulation->observables->polar_order+k)=sqrt(*(polar_x+k)**(polar_x+k)+*(polar_y+k)**(polar_y+k));
			*(simulation->observables->nematic_order+k)=sqrt(*(nematic_x+k)**(nematic_x+k)+*(nematic_y+k)**(nematic_y+k));
			*(simulation->observables->polar_order_moment1+k)+=*(simulation->observables->polar_order+k);
			*(simulation->observables->polar_order_moment2+k)+=*(simulation->observables->polar_order+k)**(simulation->observables->polar_order+k);
			*(simulation->observables->polar_order_moment3+k)+=*(simulation->observables->polar_order+k)**(simulation->observables->polar_order+k)**(simulation->observables->polar_order+k);
			*(simulation->observables->polar_order_moment4+k)+=*(simulation->observables->polar_order+k)**(simulation->observables->polar_order+k)**(simulation->observables->polar_order+k)**(simulation->observables->polar_order+k);
			*(simulation->observables->nematic_order_moment1+k)+=*(simulation->observables->nematic_order+k);
			*(simulation->observables->nematic_order_moment2+k)+=*(simulation->observables->nematic_order+k)**(simulation->observables->nematic_order+k);
			*(simulation->observables->nematic_order_moment3+k)+=*(simulation->observables->nematic_order+k)**(simulation->observables->nematic_order+k)**(simulation->observables->nematic_order+k);
			*(simulation->observables->nematic_order_moment4+k)+=*(simulation->observables->nematic_order+k)**(simulation->observables->nematic_order+k)**(simulation->observables->nematic_order+k)**(simulation->observables->nematic_order+k);
		}
	simulation->observables->measurement_steps+=1;
	free(polar_x);
	free(polar_y);
	free(nematic_x);
	free(nematic_y);
	free(neighbor_number);
	return;
}

//performs time evolution of one step for metric free/topological Vicsek model and measures observables
void mfVM_one_step_measurement(struct VM_simulation * simulation)
{
	gint32 i;
	gint32 neighbor_n;
	struct VM_particle * box;
	gint32 k;
	gdouble * polar_x=malloc(simulation->parameters->particle_species*sizeof(gdouble));
	gdouble * polar_y=malloc(simulation->parameters->particle_species*sizeof(gdouble));
	gdouble * nematic_x=malloc(simulation->parameters->particle_species*sizeof(gdouble));
	gdouble * nematic_y=malloc(simulation->parameters->particle_species*sizeof(gdouble));
	//new in version 1.1 calculate neighbor moments from ensemble average
	gdouble * neighbor_number=malloc(simulation->parameters->particle_species*sizeof(gdouble));
	for (i=0; i<simulation->parameters->particle_species; i++)
		*(neighbor_number+i)=0.1;
	for (k=0; k<simulation->parameters->particle_species; k++)
	{
		*(polar_x+k)=0.0;
		*(polar_y+k)=0.0;
		*(nematic_x+k)=0.0;
		*(nematic_y+k)=0.0;
	}
	mfVM_fill_boxes(simulation);
	for (i=0; i<simulation->parameters->mf_box_number_x*simulation->parameters->mf_box_number_y; i++)
	{
		box=*(simulation->parameters->mf_box+i);
		while (box !=NULL)
		{
			//alignment interaction
			neighbor_n=mfVM_interaction_with_neighbors(box->index, simulation);
			//analyze number of neighbors statistics
			if (simulation->parameters->particle_species==1)
			{
				//changed in version 1.1
				*neighbor_number+=1.0*neighbor_n;
				if (neighbor_n<simulation->observables->binnum_neighbors)
					*(*simulation->observables->hist_neighbors+neighbor_n)+=1;
			}
			else
			{
				//changed in version 1.1
				*(neighbor_number+get_species(box->index, simulation->parameters))+=1.0*neighbor_n;
				if (neighbor_n<simulation->observables->binnum_neighbors)
					*(*(simulation->observables->hist_neighbors+get_species(box->index, simulation->parameters))+neighbor_n)+=1;
			}
			box=box->next_particle;
		}
	}
	//normalize ensemble averaged neighbor number
	//changed in version 1.1
	if (simulation->parameters->particle_species==1)
		*neighbor_number=*neighbor_number/simulation->state->particle_number;
	else
		for (i=0; i<simulation->parameters->particle_species; i++)
			*(neighbor_number+i)=*(neighbor_number+i)/ *(simulation->parameters->species_particle_number+i);
	// update moments of ensemble averaged neighbor number
	//changed in version 1.1
	if (simulation->parameters->particle_species==1)
	{
		*simulation->observables->neighbors_moment1+=*(neighbor_number);
		*simulation->observables->neighbors_moment2+=*(neighbor_number)**(neighbor_number);
		*simulation->observables->neighbors_moment3+=*(neighbor_number)**(neighbor_number)**(neighbor_number);
		*simulation->observables->neighbors_moment4+=*(neighbor_number)**(neighbor_number)**(neighbor_number)**(neighbor_number);
	}
	else
		for (i=0; i<simulation->parameters->particle_species; i++)
		{
			*(simulation->observables->neighbors_moment1+i)+=*(neighbor_number+i);
			*(simulation->observables->neighbors_moment2+i)+=*(neighbor_number+i)**(neighbor_number+i);
			*(simulation->observables->neighbors_moment3+i)+=*(neighbor_number+i)**(neighbor_number+i)**(neighbor_number+i);
			*(simulation->observables->neighbors_moment4+i)+=*(neighbor_number+i)**(neighbor_number+i)**(neighbor_number+i)**(neighbor_number+i);
		}
	for (i=0; i<simulation->parameters->mf_box_number_x*simulation->parameters->mf_box_number_y; i++)
	{
		box=*(simulation->parameters->mf_box+i);
		while (box !=NULL)
		{
			//calculate contributions to polar and nematic order parameters
			if (simulation->parameters->particle_species==1)
			{
				*polar_x+=cos(box->theta);
				*polar_y+=sin(box->theta);
				*nematic_x+=cos(2.0*box->theta);
				*nematic_y+=sin(2.0*box->theta);
				*(*simulation->observables->hist_theta+((gint32)((M_PI+box->theta)/(2.0*M_PI)*simulation->observables->binnum_theta)))+=1;
			}
			else
			{
				*(polar_x+get_species(box->index, simulation->parameters))+=cos(box->theta);
				*(polar_y+get_species(box->index, simulation->parameters))+=sin(box->theta);
				*(nematic_x+get_species(box->index, simulation->parameters))+=cos(2.0*box->theta);
				*(nematic_y+get_species(box->index, simulation->parameters))+=sin(2.0*box->theta);
				*(*(simulation->observables->hist_theta+get_species(box->index, simulation->parameters))+((gint32)((M_PI+box->theta)/(2.0*M_PI)*simulation->observables->binnum_theta)))+=1;
			}
			//set angle to value after interaction
			box->theta=box->theta_new;
			//streaming -> movement of the particles
			if (simulation->parameters->particle_species==1)
				box->x+=*simulation->parameters->speed*cos(box->theta);
			else
				box->x+=*(simulation->parameters->speed+get_species(box->index, simulation->parameters))*cos(box->theta);
			if (simulation->parameters->boundary_x!=0) //reflecting boundary in x-direction
				if ((box->x<0.0) || (box->x>simulation->state->length_x))
				{
					box->x= -box->x;
					box->theta=M_PI-box->theta;
					if (box->theta>M_PI)
						box->theta-= 2.0*M_PI;
				}
			while (box->x <0.0)
				box->x+=simulation->state->length_x;
			while (box->x >simulation->state->length_x)
				box->x-=simulation->state->length_x;
			if (simulation->parameters->particle_species==1)
				box->y+=*simulation->parameters->speed*sin(box->theta);
			else
				box->y+=*(simulation->parameters->speed+get_species(box->index, simulation->parameters))*sin(box->theta);
			if (simulation->parameters->boundary_y!=0) //reflecting boundary in y-direction
				if ((box->y<0.0) || (box->y>simulation->state->length_y))
				{
					box->y= -box->y;
					box->theta=-box->theta;
				}
			while (box->y <0.0)
				box->y+=simulation->state->length_y;
			while (box->y >simulation->state->length_y)
				box->y-=simulation->state->length_y;
			//go to next particle in box
			box=box->next_particle;
		}
	}
	mfVM_clean_boxes(simulation->parameters);
	if (simulation->parameters->particle_species==1)
	{
		*polar_x=*polar_x/simulation->state->particle_number;
		*polar_y=*polar_y/simulation->state->particle_number;
		*nematic_x=*nematic_x/simulation->state->particle_number;
		*nematic_y=*nematic_y/simulation->state->particle_number;
		*simulation->observables->polar_order=sqrt(*polar_x**polar_x+*polar_y**polar_y);
		*simulation->observables->nematic_order=sqrt(*nematic_x**nematic_x+*nematic_y**nematic_y);
		*simulation->observables->polar_order_moment1+=*simulation->observables->polar_order;
		*simulation->observables->polar_order_moment2+=*simulation->observables->polar_order**simulation->observables->polar_order;
		*simulation->observables->polar_order_moment3+=*simulation->observables->polar_order**simulation->observables->polar_order**simulation->observables->polar_order;
		*simulation->observables->polar_order_moment4+=*simulation->observables->polar_order**simulation->observables->polar_order**simulation->observables->polar_order**simulation->observables->polar_order;
		*simulation->observables->nematic_order_moment1+=*simulation->observables->nematic_order;
		*simulation->observables->nematic_order_moment2+=*simulation->observables->nematic_order**simulation->observables->nematic_order;
		*simulation->observables->nematic_order_moment3+=*simulation->observables->nematic_order**simulation->observables->nematic_order**simulation->observables->nematic_order;
		*simulation->observables->nematic_order_moment4+=*simulation->observables->nematic_order**simulation->observables->nematic_order**simulation->observables->nematic_order**simulation->observables->nematic_order;
	}
	else
		for (k=0; k<simulation->parameters->particle_species; k++)
		{
			*(polar_x+k)=*(polar_x+k)/ *(simulation->parameters->species_particle_number+k);
			*(polar_y+k)=*(polar_y+k)/ *(simulation->parameters->species_particle_number+k);
			*(nematic_x+k)=*(nematic_x+k)/ *(simulation->parameters->species_particle_number+k);
			*(nematic_y+k)=*(nematic_y+k)/ *(simulation->parameters->species_particle_number+k);
			*(simulation->observables->polar_order+k)=sqrt(*(polar_x+k)**(polar_x+k)+*(polar_y+k)**(polar_y+k));
			*(simulation->observables->nematic_order+k)=sqrt(*(nematic_x+k)**(nematic_x+k)+*(nematic_y+k)**(nematic_y+k));
			*(simulation->observables->polar_order_moment1+k)+=*(simulation->observables->polar_order+k);
			*(simulation->observables->polar_order_moment2+k)+=*(simulation->observables->polar_order+k)**(simulation->observables->polar_order+k);
			*(simulation->observables->polar_order_moment3+k)+=*(simulation->observables->polar_order+k)**(simulation->observables->polar_order+k)**(simulation->observables->polar_order+k);
			*(simulation->observables->polar_order_moment4+k)+=*(simulation->observables->polar_order+k)**(simulation->observables->polar_order+k)**(simulation->observables->polar_order+k)**(simulation->observables->polar_order+k);
			*(simulation->observables->nematic_order_moment1+k)+=*(simulation->observables->nematic_order+k);
			*(simulation->observables->nematic_order_moment2+k)+=*(simulation->observables->nematic_order+k)**(simulation->observables->nematic_order+k);
			*(simulation->observables->nematic_order_moment3+k)+=*(simulation->observables->nematic_order+k)**(simulation->observables->nematic_order+k)**(simulation->observables->nematic_order+k);
			*(simulation->observables->nematic_order_moment4+k)+=*(simulation->observables->nematic_order+k)**(simulation->observables->nematic_order+k)**(simulation->observables->nematic_order+k)**(simulation->observables->nematic_order+k);
		}
	simulation->observables->measurement_steps+=1;
	free(polar_x);
	free(polar_y);
	free(nematic_x);
	free(nematic_y);
	free(neighbor_number);
	return;
}

//performs time evolution of one step for metric free/topological nematic Vicsek model and measures observables
void mfNVM_one_step_measurement(struct VM_simulation * simulation)
{
	gint32 i;
	gint32 neighbor_n;
	struct VM_particle * box;
	gint32 k;
	gdouble * polar_x=malloc(simulation->parameters->particle_species*sizeof(gdouble));
	gdouble * polar_y=malloc(simulation->parameters->particle_species*sizeof(gdouble));
	gdouble * nematic_x=malloc(simulation->parameters->particle_species*sizeof(gdouble));
	gdouble * nematic_y=malloc(simulation->parameters->particle_species*sizeof(gdouble));
	//new in version 1.1 calculate neighbor moments from ensemble average
	gdouble * neighbor_number=malloc(simulation->parameters->particle_species*sizeof(gdouble));
	for (i=0; i<simulation->parameters->particle_species; i++)
		*(neighbor_number+i)=0.1;
	for (k=0; k<simulation->parameters->particle_species; k++)
	{
		*(polar_x+k)=0.0;
		*(polar_y+k)=0.0;
		*(nematic_x+k)=0.0;
		*(nematic_y+k)=0.0;
	}
	mfVM_fill_boxes(simulation);
	for (i=0; i<simulation->parameters->mf_box_number_x*simulation->parameters->mf_box_number_y; i++)
	{
		box=*(simulation->parameters->mf_box+i);
		while (box !=NULL)
		{
			//alignment interaction
			neighbor_n=mfNVM_interaction_with_neighbors(box->index, simulation);
			//analyze number of neighbors statistics
			if (simulation->parameters->particle_species==1)
			{
				//changed in version 1.1
				*neighbor_number+=1.0*neighbor_n;
				if (neighbor_n<simulation->observables->binnum_neighbors)
					*(*simulation->observables->hist_neighbors+neighbor_n)+=1;
			}
			else
			{
				//changed in version 1.1
				*(neighbor_number+get_species(box->index, simulation->parameters))+=1.0*neighbor_n;
				if (neighbor_n<simulation->observables->binnum_neighbors)
					*(*(simulation->observables->hist_neighbors+get_species(box->index, simulation->parameters))+neighbor_n)+=1;
			}
			box=box->next_particle;
		}
	}
	//normalize ensemble averaged neighbor number
	//changed in version 1.1
	if (simulation->parameters->particle_species==1)
		*neighbor_number=*neighbor_number/simulation->state->particle_number;
	else
		for (i=0; i<simulation->parameters->particle_species; i++)
			*(neighbor_number+i)=*(neighbor_number+i)/ *(simulation->parameters->species_particle_number+i);
	// update moments of ensemble averaged neighbor number
	//changed in version 1.1
	if (simulation->parameters->particle_species==1)
	{
		*simulation->observables->neighbors_moment1+=*(neighbor_number);
		*simulation->observables->neighbors_moment2+=*(neighbor_number)**(neighbor_number);
		*simulation->observables->neighbors_moment3+=*(neighbor_number)**(neighbor_number)**(neighbor_number);
		*simulation->observables->neighbors_moment4+=*(neighbor_number)**(neighbor_number)**(neighbor_number)**(neighbor_number);
	}
	else
		for (i=0; i<simulation->parameters->particle_species; i++)
		{
			*(simulation->observables->neighbors_moment1+i)+=*(neighbor_number+i);
			*(simulation->observables->neighbors_moment2+i)+=*(neighbor_number+i)**(neighbor_number+i);
			*(simulation->observables->neighbors_moment3+i)+=*(neighbor_number+i)**(neighbor_number+i)**(neighbor_number+i);
			*(simulation->observables->neighbors_moment4+i)+=*(neighbor_number+i)**(neighbor_number+i)**(neighbor_number+i)**(neighbor_number+i);
		}
	for (i=0; i<simulation->parameters->mf_box_number_x*simulation->parameters->mf_box_number_y; i++)
	{
		box=*(simulation->parameters->mf_box+i);
		while (box !=NULL)
		{
			//calculate contributions to polar and nematic order parameters
			if (simulation->parameters->particle_species==1)
			{
				*polar_x+=cos(box->theta);
				*polar_y+=sin(box->theta);
				*nematic_x+=cos(2.0*box->theta);
				*nematic_y+=sin(2.0*box->theta);
				*(*simulation->observables->hist_theta+((gint32)((M_PI+box->theta)/(2.0*M_PI)*simulation->observables->binnum_theta)))+=1;
			}
			else
			{
				*(polar_x+get_species(box->index, simulation->parameters))+=cos(box->theta);
				*(polar_y+get_species(box->index, simulation->parameters))+=sin(box->theta);
				*(nematic_x+get_species(box->index, simulation->parameters))+=cos(2.0*box->theta);
				*(nematic_y+get_species(box->index, simulation->parameters))+=sin(2.0*box->theta);
				*(*(simulation->observables->hist_theta+get_species(box->index, simulation->parameters))+((gint32)((M_PI+box->theta)/(2.0*M_PI)*simulation->observables->binnum_theta)))+=1;
			}
			//set angle to value after interaction
			box->theta=box->theta_new;
			//streaming -> movement of the particles
			if (simulation->parameters->particle_species==1)
				box->x+=*simulation->parameters->speed*cos(box->theta);
			else
				box->x+=*(simulation->parameters->speed+get_species(box->index, simulation->parameters))*cos(box->theta);
			if (simulation->parameters->boundary_x!=0) //reflecting boundary in x-direction
				if ((box->x<0.0) || (box->x>simulation->state->length_x))
				{
					box->x= -box->x;
					box->theta=M_PI-box->theta;
					if (box->theta>M_PI)
						box->theta-= 2.0*M_PI;
				}
			while (box->x <0.0)
				box->x+=simulation->state->length_x;
			while (box->x >simulation->state->length_x)
				box->x-=simulation->state->length_x;
			if (simulation->parameters->particle_species==1)
				box->y+=*simulation->parameters->speed*sin(box->theta);
			else
				box->y+=*(simulation->parameters->speed+get_species(box->index, simulation->parameters))*sin(box->theta);
			if (simulation->parameters->boundary_y!=0) //reflecting boundary in y-direction
				if ((box->y<0.0) || (box->y>simulation->state->length_y))
				{
					box->y= -box->y;
					box->theta=-box->theta;
				}
			while (box->y <0.0)
				box->y+=simulation->state->length_y;
			while (box->y >simulation->state->length_y)
				box->y-=simulation->state->length_y;
			//go to next particle in box
			box=box->next_particle;
		}
	}
	mfVM_clean_boxes(simulation->parameters);
	if (simulation->parameters->particle_species==1)
	{
		*polar_x=*polar_x/simulation->state->particle_number;
		*polar_y=*polar_y/simulation->state->particle_number;
		*nematic_x=*nematic_x/simulation->state->particle_number;
		*nematic_y=*nematic_y/simulation->state->particle_number;
		*simulation->observables->polar_order=sqrt(*polar_x**polar_x+*polar_y**polar_y);
		*simulation->observables->nematic_order=sqrt(*nematic_x**nematic_x+*nematic_y**nematic_y);
		*simulation->observables->polar_order_moment1+=*simulation->observables->polar_order;
		*simulation->observables->polar_order_moment2+=*simulation->observables->polar_order**simulation->observables->polar_order;
		*simulation->observables->polar_order_moment3+=*simulation->observables->polar_order**simulation->observables->polar_order**simulation->observables->polar_order;
		*simulation->observables->polar_order_moment4+=*simulation->observables->polar_order**simulation->observables->polar_order**simulation->observables->polar_order**simulation->observables->polar_order;
		*simulation->observables->nematic_order_moment1+=*simulation->observables->nematic_order;
		*simulation->observables->nematic_order_moment2+=*simulation->observables->nematic_order**simulation->observables->nematic_order;
		*simulation->observables->nematic_order_moment3+=*simulation->observables->nematic_order**simulation->observables->nematic_order**simulation->observables->nematic_order;
		*simulation->observables->nematic_order_moment4+=*simulation->observables->nematic_order**simulation->observables->nematic_order**simulation->observables->nematic_order**simulation->observables->nematic_order;
	}
	else
		for (k=0; k<simulation->parameters->particle_species; k++)
		{
			*(polar_x+k)=*(polar_x+k)/ *(simulation->parameters->species_particle_number+k);
			*(polar_y+k)=*(polar_y+k)/ *(simulation->parameters->species_particle_number+k);
			*(nematic_x+k)=*(nematic_x+k)/ *(simulation->parameters->species_particle_number+k);
			*(nematic_y+k)=*(nematic_y+k)/ *(simulation->parameters->species_particle_number+k);
			*(simulation->observables->polar_order+k)=sqrt(*(polar_x+k)**(polar_x+k)+*(polar_y+k)**(polar_y+k));
			*(simulation->observables->nematic_order+k)=sqrt(*(nematic_x+k)**(nematic_x+k)+*(nematic_y+k)**(nematic_y+k));
			*(simulation->observables->polar_order_moment1+k)+=*(simulation->observables->polar_order+k);
			*(simulation->observables->polar_order_moment2+k)+=*(simulation->observables->polar_order+k)**(simulation->observables->polar_order+k);
			*(simulation->observables->polar_order_moment3+k)+=*(simulation->observables->polar_order+k)**(simulation->observables->polar_order+k)**(simulation->observables->polar_order+k);
			*(simulation->observables->polar_order_moment4+k)+=*(simulation->observables->polar_order+k)**(simulation->observables->polar_order+k)**(simulation->observables->polar_order+k)**(simulation->observables->polar_order+k);
			*(simulation->observables->nematic_order_moment1+k)+=*(simulation->observables->nematic_order+k);
			*(simulation->observables->nematic_order_moment2+k)+=*(simulation->observables->nematic_order+k)**(simulation->observables->nematic_order+k);
			*(simulation->observables->nematic_order_moment3+k)+=*(simulation->observables->nematic_order+k)**(simulation->observables->nematic_order+k)**(simulation->observables->nematic_order+k);
			*(simulation->observables->nematic_order_moment4+k)+=*(simulation->observables->nematic_order+k)**(simulation->observables->nematic_order+k)**(simulation->observables->nematic_order+k)**(simulation->observables->nematic_order+k);
		}
	simulation->observables->measurement_steps+=1;
	free(polar_x);
	free(polar_y);
	free(nematic_x);
	free(nematic_y);
	free(neighbor_number);
	return;
}

//performs time evolution of one step for additive Langevin model and measures observables
void additiveL_one_step_measurement(struct VM_simulation * simulation)
{
	gint32 i;
	gint32 neighbor_n;
	struct VM_particle * box;
	gint32 k;
	gdouble * polar_x=malloc(simulation->parameters->particle_species*sizeof(gdouble));
	gdouble * polar_y=malloc(simulation->parameters->particle_species*sizeof(gdouble));
	gdouble * nematic_x=malloc(simulation->parameters->particle_species*sizeof(gdouble));
	gdouble * nematic_y=malloc(simulation->parameters->particle_species*sizeof(gdouble));
	//new in version 1.1 calculate neighbor moments from ensemble average
	gdouble * neighbor_number=malloc(simulation->parameters->particle_species*sizeof(gdouble));
	for (i=0; i<simulation->parameters->particle_species; i++)
		*(neighbor_number+i)=0.1;
	for (k=0; k<simulation->parameters->particle_species; k++)
	{
		*(polar_x+k)=0.0;
		*(polar_y+k)=0.0;
		*(nematic_x+k)=0.0;
		*(nematic_y+k)=0.0;
	}
	VM_fill_boxes(simulation);
	for (i=0; i<simulation->parameters->box_number_x*simulation->parameters->box_number_y; i++)
	{
		box=*(simulation->parameters->box+i);
		while (box !=NULL)
		{
			//alignment interaction
			neighbor_n=additiveL_interaction_with_neighbors(box->index, simulation);
			//analyze number of neighbors statistics
			if (simulation->parameters->particle_species==1)
			{
				//changed in version 1.1
				*neighbor_number+=1.0*neighbor_n;
				if (neighbor_n<simulation->observables->binnum_neighbors)
					*(*simulation->observables->hist_neighbors+neighbor_n)+=1;
			}
			else
			{
				//changed in version 1.1
				*(neighbor_number+get_species(box->index, simulation->parameters))+=1.0*neighbor_n;
				if (neighbor_n<simulation->observables->binnum_neighbors)
					*(*(simulation->observables->hist_neighbors+get_species(box->index, simulation->parameters))+neighbor_n)+=1;
			}
			box=box->next_particle;
		}
	}
	//normalize ensemble averaged neighbor number
	//changed in version 1.1
	if (simulation->parameters->particle_species==1)
		*neighbor_number=*neighbor_number/simulation->state->particle_number;
	else
		for (i=0; i<simulation->parameters->particle_species; i++)
			*(neighbor_number+i)=*(neighbor_number+i)/ *(simulation->parameters->species_particle_number+i);
	// update moments of ensemble averaged neighbor number
	//changed in version 1.1
	if (simulation->parameters->particle_species==1)
	{
		*simulation->observables->neighbors_moment1+=*(neighbor_number);
		*simulation->observables->neighbors_moment2+=*(neighbor_number)**(neighbor_number);
		*simulation->observables->neighbors_moment3+=*(neighbor_number)**(neighbor_number)**(neighbor_number);
		*simulation->observables->neighbors_moment4+=*(neighbor_number)**(neighbor_number)**(neighbor_number)**(neighbor_number);
	}
	else
		for (i=0; i<simulation->parameters->particle_species; i++)
		{
			*(simulation->observables->neighbors_moment1+i)+=*(neighbor_number+i);
			*(simulation->observables->neighbors_moment2+i)+=*(neighbor_number+i)**(neighbor_number+i);
			*(simulation->observables->neighbors_moment3+i)+=*(neighbor_number+i)**(neighbor_number+i)**(neighbor_number+i);
			*(simulation->observables->neighbors_moment4+i)+=*(neighbor_number+i)**(neighbor_number+i)**(neighbor_number+i)**(neighbor_number+i);
		}
	for (i=0; i<simulation->parameters->box_number_x*simulation->parameters->box_number_y; i++)
	{
		box=*(simulation->parameters->box+i);
		while (box !=NULL)
		{
			//streaming -> movement of the particles
			if (simulation->parameters->particle_species==1)
				box->x+=simulation->parameters->delta_t**simulation->parameters->speed*cos(box->theta);
			else
				box->x+=simulation->parameters->delta_t**(simulation->parameters->speed+get_species(box->index, simulation->parameters))*cos(box->theta);
			if (simulation->parameters->particle_species==1)
				box->y+=simulation->parameters->delta_t**simulation->parameters->speed*sin(box->theta);
			else
				box->y+=simulation->parameters->delta_t**(simulation->parameters->speed+get_species(box->index, simulation->parameters))*sin(box->theta);
			//apply boundary conditions
			if (simulation->parameters->boundary_x!=0) //reflecting boundary in x-direction
				if ((box->x<0.0) || (box->x>simulation->state->length_x))
				{
					box->x= -box->x;
					//mistake in implementation for reflecting bc corrected in version 1.2
					box->theta_new= M_PI-box->theta_new;
					if (box->theta_new>M_PI)
						box->theta_new-= 2.0*M_PI;
				}
			while (box->x <0.0)
				box->x+=simulation->state->length_x;
			while (box->x >simulation->state->length_x)
				box->x-=simulation->state->length_x;
			if (simulation->parameters->boundary_y!=0) //reflecting boundary in y-direction
				if ((box->y<0.0) || (box->y>simulation->state->length_y))
				{
					box->y= -box->y;
					//mistake in implementation for reflecting bc corrected in version 1.2
					box->theta_new= -box->theta_new;
				}
			while (box->y <0.0)
				box->y+=simulation->state->length_y;
			while (box->y >simulation->state->length_y)
				box->y-=simulation->state->length_y;
			//calculate contribution to polar/nematic order parameter
			if (simulation->parameters->particle_species==1)
			{
				*polar_x+=cos(box->theta);
				*polar_y+=sin(box->theta);
				*nematic_x+=cos(2.0*box->theta);
				*nematic_y+=sin(2.0*box->theta);
				*(*simulation->observables->hist_theta+((gint32)((M_PI+box->theta)/(2.0*M_PI)*simulation->observables->binnum_theta)))+=1;
			}
			else
			{
				*(polar_x+get_species(box->index, simulation->parameters))+=cos(box->theta);
				*(polar_y+get_species(box->index, simulation->parameters))+=sin(box->theta);
				*(nematic_x+get_species(box->index, simulation->parameters))+=cos(2.0*box->theta);
				*(nematic_y+get_species(box->index, simulation->parameters))+=sin(2.0*box->theta);
				*(*(simulation->observables->hist_theta+get_species(box->index, simulation->parameters))+((gint32)((M_PI+box->theta)/(2.0*M_PI)*simulation->observables->binnum_theta)))+=1;
			}
			//update angles
			box->theta=box->theta_new;
			box=box->next_particle;
		}
	}
	VM_clean_boxes(simulation->parameters);
	if (simulation->parameters->particle_species==1)
	{
		*polar_x=*polar_x/simulation->state->particle_number;
		*polar_y=*polar_y/simulation->state->particle_number;
		*nematic_x=*nematic_x/simulation->state->particle_number;
		*nematic_y=*nematic_y/simulation->state->particle_number;
		*simulation->observables->polar_order=sqrt(*polar_x**polar_x+*polar_y**polar_y);
		*simulation->observables->nematic_order=sqrt(*nematic_x**nematic_x+*nematic_y**nematic_y);
		*simulation->observables->polar_order_moment1+=*simulation->observables->polar_order;
		*simulation->observables->polar_order_moment2+=*simulation->observables->polar_order**simulation->observables->polar_order;
		*simulation->observables->polar_order_moment3+=*simulation->observables->polar_order**simulation->observables->polar_order**simulation->observables->polar_order;
		*simulation->observables->polar_order_moment4+=*simulation->observables->polar_order**simulation->observables->polar_order**simulation->observables->polar_order**simulation->observables->polar_order;
		*simulation->observables->nematic_order_moment1+=*simulation->observables->nematic_order;
		*simulation->observables->nematic_order_moment2+=*simulation->observables->nematic_order**simulation->observables->nematic_order;
		*simulation->observables->nematic_order_moment3+=*simulation->observables->nematic_order**simulation->observables->nematic_order**simulation->observables->nematic_order;
		*simulation->observables->nematic_order_moment4+=*simulation->observables->nematic_order**simulation->observables->nematic_order**simulation->observables->nematic_order**simulation->observables->nematic_order;
	}
	else
		for (k=0; k<simulation->parameters->particle_species; k++)
		{
			*(polar_x+k)=*(polar_x+k)/ *(simulation->parameters->species_particle_number+k);
			*(polar_y+k)=*(polar_y+k)/ *(simulation->parameters->species_particle_number+k);
			*(nematic_x+k)=*(nematic_x+k)/ *(simulation->parameters->species_particle_number+k);
			*(nematic_y+k)=*(nematic_y+k)/ *(simulation->parameters->species_particle_number+k);
			*(simulation->observables->polar_order+k)=sqrt(*(polar_x+k)**(polar_x+k)+*(polar_y+k)**(polar_y+k));
			*(simulation->observables->nematic_order+k)=sqrt(*(nematic_x+k)**(nematic_x+k)+*(nematic_y+k)**(nematic_y+k));
			*(simulation->observables->polar_order_moment1+k)+=*(simulation->observables->polar_order+k);
			*(simulation->observables->polar_order_moment2+k)+=*(simulation->observables->polar_order+k)**(simulation->observables->polar_order+k);
			*(simulation->observables->polar_order_moment3+k)+=*(simulation->observables->polar_order+k)**(simulation->observables->polar_order+k)**(simulation->observables->polar_order+k);
			*(simulation->observables->polar_order_moment4+k)+=*(simulation->observables->polar_order+k)**(simulation->observables->polar_order+k)**(simulation->observables->polar_order+k)**(simulation->observables->polar_order+k);
			*(simulation->observables->nematic_order_moment1+k)+=*(simulation->observables->nematic_order+k);
			*(simulation->observables->nematic_order_moment2+k)+=*(simulation->observables->nematic_order+k)**(simulation->observables->nematic_order+k);
			*(simulation->observables->nematic_order_moment3+k)+=*(simulation->observables->nematic_order+k)**(simulation->observables->nematic_order+k)**(simulation->observables->nematic_order+k);
			*(simulation->observables->nematic_order_moment4+k)+=*(simulation->observables->nematic_order+k)**(simulation->observables->nematic_order+k)**(simulation->observables->nematic_order+k)**(simulation->observables->nematic_order+k);
		}
	simulation->observables->measurement_steps+=1;
	free(polar_x);
	free(polar_y);
	free(nematic_x);
	free(nematic_y);
	free(neighbor_number);
	return;
}

//performs time evolution of one step for nonadditive Langevin model and measures observables
void nonadditiveL_one_step_measurement(struct VM_simulation * simulation)
{
	gint32 i;
	gint32 neighbor_n;
	struct VM_particle * box;
	gint32 k;
	gdouble * polar_x=malloc(simulation->parameters->particle_species*sizeof(gdouble));
	gdouble * polar_y=malloc(simulation->parameters->particle_species*sizeof(gdouble));
	gdouble * nematic_x=malloc(simulation->parameters->particle_species*sizeof(gdouble));
	gdouble * nematic_y=malloc(simulation->parameters->particle_species*sizeof(gdouble));
	//new in version 1.1 calculate neighbor moments from ensemble average
	gdouble * neighbor_number=malloc(simulation->parameters->particle_species*sizeof(gdouble));
	for (i=0; i<simulation->parameters->particle_species; i++)
		*(neighbor_number+i)=0.1;
	for (k=0; k<simulation->parameters->particle_species; k++)
	{
		*(polar_x+k)=0.0;
		*(polar_y+k)=0.0;
		*(nematic_x+k)=0.0;
		*(nematic_y+k)=0.0;
	}
	VM_fill_boxes(simulation);
	for (i=0; i<simulation->parameters->box_number_x*simulation->parameters->box_number_y; i++)
	{
		box=*(simulation->parameters->box+i);
		while (box !=NULL)
		{
			//alignment interaction
			neighbor_n=nonadditiveL_interaction_with_neighbors(box->index, simulation);
			//analyze number of neighbors statistics
			if (simulation->parameters->particle_species==1)
			{
				//changed in version 1.1
				*neighbor_number+=1.0*neighbor_n;
				if (neighbor_n<simulation->observables->binnum_neighbors)
					*(*simulation->observables->hist_neighbors+neighbor_n)+=1;
			}
			else
			{
				//changed in version 1.1
				*(neighbor_number+get_species(box->index, simulation->parameters))+=1.0*neighbor_n;
				if (neighbor_n<simulation->observables->binnum_neighbors)
					*(*(simulation->observables->hist_neighbors+get_species(box->index, simulation->parameters))+neighbor_n)+=1;
			}
			box=box->next_particle;
		}
	}
	//normalize ensemble averaged neighbor number
	//changed in version 1.1
	if (simulation->parameters->particle_species==1)
		*neighbor_number=*neighbor_number/simulation->state->particle_number;
	else
		for (i=0; i<simulation->parameters->particle_species; i++)
			*(neighbor_number+i)=*(neighbor_number+i)/ *(simulation->parameters->species_particle_number+i);
	// update moments of ensemble averaged neighbor number
	//changed in version 1.1
	if (simulation->parameters->particle_species==1)
	{
		*simulation->observables->neighbors_moment1+=*(neighbor_number);
		*simulation->observables->neighbors_moment2+=*(neighbor_number)**(neighbor_number);
		*simulation->observables->neighbors_moment3+=*(neighbor_number)**(neighbor_number)**(neighbor_number);
		*simulation->observables->neighbors_moment4+=*(neighbor_number)**(neighbor_number)**(neighbor_number)**(neighbor_number);
	}
	else
		for (i=0; i<simulation->parameters->particle_species; i++)
		{
			*(simulation->observables->neighbors_moment1+i)+=*(neighbor_number+i);
			*(simulation->observables->neighbors_moment2+i)+=*(neighbor_number+i)**(neighbor_number+i);
			*(simulation->observables->neighbors_moment3+i)+=*(neighbor_number+i)**(neighbor_number+i)**(neighbor_number+i);
			*(simulation->observables->neighbors_moment4+i)+=*(neighbor_number+i)**(neighbor_number+i)**(neighbor_number+i)**(neighbor_number+i);
		}
	for (i=0; i<simulation->parameters->box_number_x*simulation->parameters->box_number_y; i++)
	{
		box=*(simulation->parameters->box+i);
		while (box !=NULL)
		{
			//streaming -> movement of the particles
			if (simulation->parameters->particle_species==1)
				box->x+=simulation->parameters->delta_t**simulation->parameters->speed*cos(box->theta);
			else
				box->x+=simulation->parameters->delta_t**(simulation->parameters->speed+get_species(box->index, simulation->parameters))*cos(box->theta);
			if (simulation->parameters->particle_species==1)
				box->y+=simulation->parameters->delta_t**simulation->parameters->speed*sin(box->theta);
			else
				box->y+=simulation->parameters->delta_t**(simulation->parameters->speed+get_species(box->index, simulation->parameters))*sin(box->theta);
			//apply boundary conditions
			if (simulation->parameters->boundary_x!=0) //reflecting boundary in x-direction
				if ((box->x<0.0) || (box->x>simulation->state->length_x))
				{
					box->x= -box->x;
					//mistake in implementation for reflecting bc corrected in version 1.2
					box->theta_new= M_PI-box->theta_new;
					if (box->theta_new>M_PI)
						box->theta_new-= 2.0*M_PI;
				}
			while (box->x <0.0)
				box->x+=simulation->state->length_x;
			while (box->x >simulation->state->length_x)
				box->x-=simulation->state->length_x;
			if (simulation->parameters->boundary_y!=0) //reflecting boundary in y-direction
				if ((box->y<0.0) || (box->y>simulation->state->length_y))
				{
					box->y= -box->y;
					//mistake in implementation for reflecting bc corrected in version 1.2
					box->theta_new= -box->theta_new;
				}
			while (box->y <0.0)
				box->y+=simulation->state->length_y;
			while (box->y >simulation->state->length_y)
				box->y-=simulation->state->length_y;
			//calculate contribution to polar/nematic order parameter
			if (simulation->parameters->particle_species==1)
			{
				*polar_x+=cos(box->theta);
				*polar_y+=sin(box->theta);
				*nematic_x+=cos(2.0*box->theta);
				*nematic_y+=sin(2.0*box->theta);
				*(*simulation->observables->hist_theta+((gint32)((M_PI+box->theta)/(2.0*M_PI)*simulation->observables->binnum_theta)))+=1;
			}
			else
			{
				*(polar_x+get_species(box->index, simulation->parameters))+=cos(box->theta);
				*(polar_y+get_species(box->index, simulation->parameters))+=sin(box->theta);
				*(nematic_x+get_species(box->index, simulation->parameters))+=cos(2.0*box->theta);
				*(nematic_y+get_species(box->index, simulation->parameters))+=sin(2.0*box->theta);
				*(*(simulation->observables->hist_theta+get_species(box->index, simulation->parameters))+((gint32)((M_PI+box->theta)/(2.0*M_PI)*simulation->observables->binnum_theta)))+=1;
			}
			//update angles
			box->theta=box->theta_new;
			box=box->next_particle;
		}
	}
	VM_clean_boxes(simulation->parameters);
	if (simulation->parameters->particle_species==1)
	{
		*polar_x=*polar_x/simulation->state->particle_number;
		*polar_y=*polar_y/simulation->state->particle_number;
		*nematic_x=*nematic_x/simulation->state->particle_number;
		*nematic_y=*nematic_y/simulation->state->particle_number;
		*simulation->observables->polar_order=sqrt(*polar_x**polar_x+*polar_y**polar_y);
		*simulation->observables->nematic_order=sqrt(*nematic_x**nematic_x+*nematic_y**nematic_y);
		*simulation->observables->polar_order_moment1+=*simulation->observables->polar_order;
		*simulation->observables->polar_order_moment2+=*simulation->observables->polar_order**simulation->observables->polar_order;
		*simulation->observables->polar_order_moment3+=*simulation->observables->polar_order**simulation->observables->polar_order**simulation->observables->polar_order;
		*simulation->observables->polar_order_moment4+=*simulation->observables->polar_order**simulation->observables->polar_order**simulation->observables->polar_order**simulation->observables->polar_order;
		*simulation->observables->nematic_order_moment1+=*simulation->observables->nematic_order;
		*simulation->observables->nematic_order_moment2+=*simulation->observables->nematic_order**simulation->observables->nematic_order;
		*simulation->observables->nematic_order_moment3+=*simulation->observables->nematic_order**simulation->observables->nematic_order**simulation->observables->nematic_order;
		*simulation->observables->nematic_order_moment4+=*simulation->observables->nematic_order**simulation->observables->nematic_order**simulation->observables->nematic_order**simulation->observables->nematic_order;
	}
	else
		for (k=0; k<simulation->parameters->particle_species; k++)
		{
			*(polar_x+k)=*(polar_x+k)/ *(simulation->parameters->species_particle_number+k);
			*(polar_y+k)=*(polar_y+k)/ *(simulation->parameters->species_particle_number+k);
			*(nematic_x+k)=*(nematic_x+k)/ *(simulation->parameters->species_particle_number+k);
			*(nematic_y+k)=*(nematic_y+k)/ *(simulation->parameters->species_particle_number+k);
			*(simulation->observables->polar_order+k)=sqrt(*(polar_x+k)**(polar_x+k)+*(polar_y+k)**(polar_y+k));
			*(simulation->observables->nematic_order+k)=sqrt(*(nematic_x+k)**(nematic_x+k)+*(nematic_y+k)**(nematic_y+k));
			*(simulation->observables->polar_order_moment1+k)+=*(simulation->observables->polar_order+k);
			*(simulation->observables->polar_order_moment2+k)+=*(simulation->observables->polar_order+k)**(simulation->observables->polar_order+k);
			*(simulation->observables->polar_order_moment3+k)+=*(simulation->observables->polar_order+k)**(simulation->observables->polar_order+k)**(simulation->observables->polar_order+k);
			*(simulation->observables->polar_order_moment4+k)+=*(simulation->observables->polar_order+k)**(simulation->observables->polar_order+k)**(simulation->observables->polar_order+k)**(simulation->observables->polar_order+k);
			*(simulation->observables->nematic_order_moment1+k)+=*(simulation->observables->nematic_order+k);
			*(simulation->observables->nematic_order_moment2+k)+=*(simulation->observables->nematic_order+k)**(simulation->observables->nematic_order+k);
			*(simulation->observables->nematic_order_moment3+k)+=*(simulation->observables->nematic_order+k)**(simulation->observables->nematic_order+k)**(simulation->observables->nematic_order+k);
			*(simulation->observables->nematic_order_moment4+k)+=*(simulation->observables->nematic_order+k)**(simulation->observables->nematic_order+k)**(simulation->observables->nematic_order+k)**(simulation->observables->nematic_order+k);
		}
	simulation->observables->measurement_steps+=1;
	free(polar_x);
	free(polar_y);
	free(nematic_x);
	free(nematic_y);
	free(neighbor_number);
	return;
}

//performs time evolution of one step for metric free/topological Langevin model and measures observables
void mfL_one_step_measurement(struct VM_simulation * simulation)
{
	gint32 i;
	gint32 neighbor_n;
	struct VM_particle * box;
	gint32 k;
	gdouble * polar_x=malloc(simulation->parameters->particle_species*sizeof(gdouble));
	gdouble * polar_y=malloc(simulation->parameters->particle_species*sizeof(gdouble));
	gdouble * nematic_x=malloc(simulation->parameters->particle_species*sizeof(gdouble));
	gdouble * nematic_y=malloc(simulation->parameters->particle_species*sizeof(gdouble));
	//new in version 1.1 calculate neighbor moments from ensemble average
	gdouble * neighbor_number=malloc(simulation->parameters->particle_species*sizeof(gdouble));
	for (i=0; i<simulation->parameters->particle_species; i++)
		*(neighbor_number+i)=0.1;
	for (k=0; k<simulation->parameters->particle_species; k++)
	{
		*(polar_x+k)=0.0;
		*(polar_y+k)=0.0;
		*(nematic_x+k)=0.0;
		*(nematic_y+k)=0.0;
	}
	mfVM_fill_boxes(simulation);
	for (i=0; i<simulation->parameters->mf_box_number_x*simulation->parameters->mf_box_number_y; i++)
	{
		box=*(simulation->parameters->mf_box+i);
		while (box !=NULL)
		{
			//alignment interaction
			neighbor_n=mfL_interaction_with_neighbors(box->index, simulation);
			//analyze number of neighbors statistics
			if (simulation->parameters->particle_species==1)
			{
				//changed in version 1.1
				*neighbor_number+=1.0*neighbor_n;
				if (neighbor_n<simulation->observables->binnum_neighbors)
					*(*simulation->observables->hist_neighbors+neighbor_n)+=1;
			}
			else
			{
				//changed in version 1.1
				*(neighbor_number+get_species(box->index, simulation->parameters))+=1.0*neighbor_n;
				if (neighbor_n<simulation->observables->binnum_neighbors)
					*(*(simulation->observables->hist_neighbors+get_species(box->index, simulation->parameters))+neighbor_n)+=1;
			}
			box=box->next_particle;
		}
	}
	//normalize ensemble averaged neighbor number
	//changed in version 1.1
	if (simulation->parameters->particle_species==1)
		*neighbor_number=*neighbor_number/simulation->state->particle_number;
	else
		for (i=0; i<simulation->parameters->particle_species; i++)
			*(neighbor_number+i)=*(neighbor_number+i)/ *(simulation->parameters->species_particle_number+i);
	// update moments of ensemble averaged neighbor number
	//changed in version 1.1
	if (simulation->parameters->particle_species==1)
	{
		*simulation->observables->neighbors_moment1+=*(neighbor_number);
		*simulation->observables->neighbors_moment2+=*(neighbor_number)**(neighbor_number);
		*simulation->observables->neighbors_moment3+=*(neighbor_number)**(neighbor_number)**(neighbor_number);
		*simulation->observables->neighbors_moment4+=*(neighbor_number)**(neighbor_number)**(neighbor_number)**(neighbor_number);
	}
	else
		for (i=0; i<simulation->parameters->particle_species; i++)
		{
			*(simulation->observables->neighbors_moment1+i)+=*(neighbor_number+i);
			*(simulation->observables->neighbors_moment2+i)+=*(neighbor_number+i)**(neighbor_number+i);
			*(simulation->observables->neighbors_moment3+i)+=*(neighbor_number+i)**(neighbor_number+i)**(neighbor_number+i);
			*(simulation->observables->neighbors_moment4+i)+=*(neighbor_number+i)**(neighbor_number+i)**(neighbor_number+i)**(neighbor_number+i);
		}
	for (i=0; i<simulation->parameters->mf_box_number_x*simulation->parameters->mf_box_number_y; i++)
	{
		box=*(simulation->parameters->mf_box+i);
		while (box !=NULL)
		{
			//streaming -> movement of the particles
			if (simulation->parameters->particle_species==1)
				box->x+=simulation->parameters->delta_t**simulation->parameters->speed*cos(box->theta);
			else
				box->x+=simulation->parameters->delta_t**(simulation->parameters->speed+get_species(box->index, simulation->parameters))*cos(box->theta);
			if (simulation->parameters->particle_species==1)
				box->y+=simulation->parameters->delta_t**simulation->parameters->speed*sin(box->theta);
			else
				box->y+=simulation->parameters->delta_t**(simulation->parameters->speed+get_species(box->index, simulation->parameters))*sin(box->theta);
			//apply boundary conditions
			if (simulation->parameters->boundary_x!=0) //reflecting boundary in x-direction
				if ((box->x<0.0) || (box->x>simulation->state->length_x))
				{
					box->x= -box->x;
					//mistake in implementation for reflecting bc corrected in version 1.2
					box->theta_new= M_PI-box->theta_new;
					if (box->theta_new>M_PI)
						box->theta_new-= 2.0*M_PI;
				}
			while (box->x <0.0)
				box->x+=simulation->state->length_x;
			while (box->x >simulation->state->length_x)
				box->x-=simulation->state->length_x;
			if (simulation->parameters->boundary_y!=0) //reflecting boundary in y-direction
				if ((box->y<0.0) || (box->y>simulation->state->length_y))
				{
					box->y= -box->y;
					//mistake in implementation for reflecting bc corrected in version 1.2
					box->theta_new= -box->theta_new;
				}
			while (box->y <0.0)
				box->y+=simulation->state->length_y;
			while (box->y >simulation->state->length_y)
				box->y-=simulation->state->length_y;
			//calculate contribution to polar/nematic order parameter
			if (simulation->parameters->particle_species==1)
			{
				*polar_x+=cos(box->theta);
				*polar_y+=sin(box->theta);
				*nematic_x+=cos(2.0*box->theta);
				*nematic_y+=sin(2.0*box->theta);
				*(*simulation->observables->hist_theta+((gint32)((M_PI+box->theta)/(2.0*M_PI)*simulation->observables->binnum_theta)))+=1;
			}
			else
			{
				*(polar_x+get_species(box->index, simulation->parameters))+=cos(box->theta);
				*(polar_y+get_species(box->index, simulation->parameters))+=sin(box->theta);
				*(nematic_x+get_species(box->index, simulation->parameters))+=cos(2.0*box->theta);
				*(nematic_y+get_species(box->index, simulation->parameters))+=sin(2.0*box->theta);
				*(*(simulation->observables->hist_theta+get_species(box->index, simulation->parameters))+((gint32)((M_PI+box->theta)/(2.0*M_PI)*simulation->observables->binnum_theta)))+=1;
			}
			//update angles
			box->theta=box->theta_new;
			box=box->next_particle;
		}
	}
	mfVM_clean_boxes(simulation->parameters);
	if (simulation->parameters->particle_species==1)
	{
		*polar_x=*polar_x/simulation->state->particle_number;
		*polar_y=*polar_y/simulation->state->particle_number;
		*nematic_x=*nematic_x/simulation->state->particle_number;
		*nematic_y=*nematic_y/simulation->state->particle_number;
		*simulation->observables->polar_order=sqrt(*polar_x**polar_x+*polar_y**polar_y);
		*simulation->observables->nematic_order=sqrt(*nematic_x**nematic_x+*nematic_y**nematic_y);
		*simulation->observables->polar_order_moment1+=*simulation->observables->polar_order;
		*simulation->observables->polar_order_moment2+=*simulation->observables->polar_order**simulation->observables->polar_order;
		*simulation->observables->polar_order_moment3+=*simulation->observables->polar_order**simulation->observables->polar_order**simulation->observables->polar_order;
		*simulation->observables->polar_order_moment4+=*simulation->observables->polar_order**simulation->observables->polar_order**simulation->observables->polar_order**simulation->observables->polar_order;
		*simulation->observables->nematic_order_moment1+=*simulation->observables->nematic_order;
		*simulation->observables->nematic_order_moment2+=*simulation->observables->nematic_order**simulation->observables->nematic_order;
		*simulation->observables->nematic_order_moment3+=*simulation->observables->nematic_order**simulation->observables->nematic_order**simulation->observables->nematic_order;
		*simulation->observables->nematic_order_moment4+=*simulation->observables->nematic_order**simulation->observables->nematic_order**simulation->observables->nematic_order**simulation->observables->nematic_order;
	}
	else
		for (k=0; k<simulation->parameters->particle_species; k++)
		{
			*(polar_x+k)=*(polar_x+k)/ *(simulation->parameters->species_particle_number+k);
			*(polar_y+k)=*(polar_y+k)/ *(simulation->parameters->species_particle_number+k);
			*(nematic_x+k)=*(nematic_x+k)/ *(simulation->parameters->species_particle_number+k);
			*(nematic_y+k)=*(nematic_y+k)/ *(simulation->parameters->species_particle_number+k);
			*(simulation->observables->polar_order+k)=sqrt(*(polar_x+k)**(polar_x+k)+*(polar_y+k)**(polar_y+k));
			*(simulation->observables->nematic_order+k)=sqrt(*(nematic_x+k)**(nematic_x+k)+*(nematic_y+k)**(nematic_y+k));
			*(simulation->observables->polar_order_moment1+k)+=*(simulation->observables->polar_order+k);
			*(simulation->observables->polar_order_moment2+k)+=*(simulation->observables->polar_order+k)**(simulation->observables->polar_order+k);
			*(simulation->observables->polar_order_moment3+k)+=*(simulation->observables->polar_order+k)**(simulation->observables->polar_order+k)**(simulation->observables->polar_order+k);
			*(simulation->observables->polar_order_moment4+k)+=*(simulation->observables->polar_order+k)**(simulation->observables->polar_order+k)**(simulation->observables->polar_order+k)**(simulation->observables->polar_order+k);
			*(simulation->observables->nematic_order_moment1+k)+=*(simulation->observables->nematic_order+k);
			*(simulation->observables->nematic_order_moment2+k)+=*(simulation->observables->nematic_order+k)**(simulation->observables->nematic_order+k);
			*(simulation->observables->nematic_order_moment3+k)+=*(simulation->observables->nematic_order+k)**(simulation->observables->nematic_order+k)**(simulation->observables->nematic_order+k);
			*(simulation->observables->nematic_order_moment4+k)+=*(simulation->observables->nematic_order+k)**(simulation->observables->nematic_order+k)**(simulation->observables->nematic_order+k)**(simulation->observables->nematic_order+k);
		}
	simulation->observables->measurement_steps+=1;
	free(polar_x);
	free(polar_y);
	free(nematic_x);
	free(nematic_y);
	free(neighbor_number);
	return;
}

//saves the state of the simulation in binary format to file 'filename'
gint32 VM_save_state(char * filename, struct VM_simulation * simulation)
{
	gint32 myfile;
	guint64 i;
	gint32 j, k;
	//open file for writing, create if it does not exist, give user reading and writing rights
	if ((myfile=open(filename, O_WRONLY|O_CREAT|O_TRUNC, S_IRUSR|S_IWUSR))==-1)
	{
		PyErr_SetString(PyExc_TypeError, "Error when opening file for saving simulation data\n");
		return -1;
	}
	//write a version flag, this is version 1.1, it writes 1, it can only read/write files that have this flag=1, or the previous flag=0, later versions might write different data files and indicate this by using another value for this flag
	k=1;
	if ((write(myfile, &k, sizeof(gint32)))==-1)
	{
		PyErr_SetString(PyExc_TypeError, "Error when writing version flag\n");
		return -1;
	}
	//write some parameters
		//write particle number
	if ((write(myfile, &(simulation->state->particle_number), sizeof(guint64)))==-1)
	{
		PyErr_SetString(PyExc_TypeError, "Error when writing total particle number\n");
		return -1;
	}
		//write particle species number
	if ((write(myfile, &(simulation->parameters->particle_species), sizeof(gint32)))==-1)
	{
		PyErr_SetString(PyExc_TypeError, "Error when writing number of particle species\n");
		return -1;
	}
		//if there is more than one species write number of particles for each species
	if (simulation->parameters->particle_species > 1)
		for(k=0; k<simulation->parameters->particle_species; k++)
		{
			if ((write(myfile, simulation->parameters->species_particle_number+k, sizeof(guint64)))==-1)
			{
				PyErr_SetString(PyExc_TypeError, "Error when writing particle number for each species\n");
				return -1;
			}
		}
		//write interaction radius
	if ((write(myfile, &(simulation->parameters->interaction_radius), sizeof(gdouble)))==-1)
	{
		PyErr_SetString(PyExc_TypeError, "Error when writing interaction radius\n");
		return -1;
	}
		//write step size
	if ((write(myfile, &(simulation->parameters->delta_t), sizeof(gdouble)))==-1)
	{
		PyErr_SetString(PyExc_TypeError, "Error when writing step size\n");
		return -1;
	}
		//write metric free neighbor number
	if ((write(myfile, &(simulation->parameters->mf_kn), sizeof(gint32)))==-1)
	{
		PyErr_SetString(PyExc_TypeError, "Error when writing interaction neighbor number for metric free models\n");
		return -1;
	}
		//write interaction order
	if ((write(myfile, &(simulation->parameters->interaction_order), sizeof(gint32)))==-1)
	{
		PyErr_SetString(PyExc_TypeError, "Error when writing interaction order (for overdamped Langevin models)\n");
		return -1;
	}
		//write boundary flag for x-direction
	if ((write(myfile, &(simulation->parameters->boundary_x), sizeof(gint32)))==-1)
	{
		PyErr_SetString(PyExc_TypeError, "Error when writing boundary condition flag for x-direction\n");
		return -1;
	}
		//write boundary flag for y-direction
	if ((write(myfile, &(simulation->parameters->boundary_y), sizeof(gint32)))==-1)
	{
		PyErr_SetString(PyExc_TypeError, "Error when writing boundary condition flag for y-direction\n");
		return -1;
	}
		//write natural rotation frequency(ies) omega
	for(k=0; k<simulation->parameters->particle_species; k++)
	{
		if ((write(myfile, simulation->parameters->omega+k, sizeof(gdouble)))==-1)
		{
			PyErr_SetString(PyExc_TypeError, "Error when writing omega(s)\n");
			return -1;
		}
	}
		//write noise strength(s) eta
	for(k=0; k<simulation->parameters->particle_species; k++)
	{
		if ((write(myfile, simulation->parameters->noise_strength+k, sizeof(gdouble)))==-1)
		{
			PyErr_SetString(PyExc_TypeError, "Error when writing eta(s)\n");
			return -1;
		}
	}
		//write velocity(ies) v
	for(k=0; k<simulation->parameters->particle_species; k++)
	{
		if ((write(myfile, simulation->parameters->speed+k, sizeof(gdouble)))==-1)
		{
			PyErr_SetString(PyExc_TypeError, "Error when writing velocity(ies)\n");
			return -1;
		}
	}
		//write coupling(s) gamma
	for (j=0; j<simulation->parameters->particle_species; j++)
	{
		for (k=0; k<simulation->parameters->particle_species; k++)
		{
			if ((write(myfile, *(simulation->parameters->coupling+j)+k, sizeof(gdouble)))==-1)
			{
				PyErr_SetString(PyExc_TypeError, "Error when writing coupling matrix to file\n");
				return -1;
			}
		}
	}
		//write simulation box size in x-direction Lx
	if ((write(myfile, &(simulation->state->length_x), sizeof(gdouble)))==-1)
	{
		PyErr_SetString(PyExc_TypeError, "Error when writing Lx\n");
		return -1;
	}
		//write simulation box size in y-direction Ly
	if ((write(myfile, &(simulation->state->length_y), sizeof(gdouble)))==-1)
	{
		PyErr_SetString(PyExc_TypeError, "Error when writing Ly\n");
		return -1;
	}
	//write internal states related to random number generation and weight function
		//write gauss1
	if ((write(myfile, &(simulation->state->gauss1), sizeof(gdouble)))==-1)
	{
		PyErr_SetString(PyExc_TypeError, "Error when writing gauss1\n");
		return -1;
	}
		//write gauss2
	if ((write(myfile, &(simulation->state->gauss2), sizeof(gdouble)))==-1)
	{
		PyErr_SetString(PyExc_TypeError, "Error when writing gauss2\n");
		return -1;
	}
		//write box_muller flag
	if ((write(myfile, &(simulation->state->box_muller), sizeof(gint32)))==-1)
	{
		PyErr_SetString(PyExc_TypeError, "Error when writing box muller flag\n");
		return -1;
	}
		//write pseuda random number state variable
	if ((write(myfile, simulation->state->prng, 624*sizeof(guint32)+sizeof(guint)))==-1)
	{
		PyErr_SetString(PyExc_TypeError, "Error when writing pseudo random number generator state\n");
		return -1;
	}
		//write weight vector length
	if ((write(myfile, &(simulation->state->weight_vector_length), sizeof(gint32)))==-1)
	{
		PyErr_SetString(PyExc_TypeError, "Error when writing weight vector length\n");
		return -1;
	}
		//write weights (if there are some)
	if (simulation->state->weight_vector_length>0)
		if ((write(myfile, simulation->state->weights, simulation->state->weight_vector_length*sizeof(gdouble)))==-1)
		{
			PyErr_SetString(PyExc_TypeError, "Error when writing weight vector\n");
			return -1;
		}
	//write the states of all particles
		//write x-positions
	for(i=0; i<simulation->state->particle_number; i++)
	{
		if ((write(myfile, &((simulation->state->particles+i)->x), sizeof(gdouble)))==-1)
		{
			PyErr_SetString(PyExc_TypeError, "Error when writing x-positions\n");
			return -1;
		}
	}
		//write y-positions
	for(i=0; i<simulation->state->particle_number; i++)
	{
		if ((write(myfile, &((simulation->state->particles+i)->y), sizeof(gdouble)))==-1)
		{
			PyErr_SetString(PyExc_TypeError, "Error when writing y-positions\n");
			return -1;
		}
	}
		//write orientations
	for(i=0; i<simulation->state->particle_number; i++)
	{
		if ((write(myfile, &((simulation->state->particles+i)->theta), sizeof(gdouble)))==-1)
		{
			PyErr_SetString(PyExc_TypeError, "Error when writing orientations\n");
			return -1;
		}
	}
	//write the results of measurements
		//write polar order(s)
	for (k=0; k<simulation->parameters->particle_species; k++)
	{
		if ((write(myfile, (simulation->observables->polar_order+k), sizeof(gdouble)))==-1)
		{
			PyErr_SetString(PyExc_TypeError, "Error when writing polar order(s)\n");
			return -1;
		}
	}
		//write first moment of polar order(s)
	for (k=0; k<simulation->parameters->particle_species; k++)
	{
		if ((write(myfile, (simulation->observables->polar_order_moment1+k), sizeof(gdouble)))==-1)
		{
			PyErr_SetString(PyExc_TypeError, "Error when writing polar order first moment(s)\n");
			return -1;
		}
	}
		//write second moment of polar order(s)
	for (k=0; k<simulation->parameters->particle_species; k++)
	{
		if ((write(myfile, (simulation->observables->polar_order_moment2+k), sizeof(gdouble)))==-1)
		{
			PyErr_SetString(PyExc_TypeError, "Error when writing polar order second moment(s)\n");
			return -1;
		}
	}
		//write third moment of polar order(s)
	for (k=0; k<simulation->parameters->particle_species; k++)
	{
		if ((write(myfile, (simulation->observables->polar_order_moment3+k), sizeof(gdouble)))==-1)
		{
			PyErr_SetString(PyExc_TypeError, "Error when writing polar order third moment(s)\n");
			return -1;
		}
	}
		//write fourth moment of polar order(s)
	for (k=0; k<simulation->parameters->particle_species; k++)
	{
		if ((write(myfile, (simulation->observables->polar_order_moment4+k), sizeof(gdouble)))==-1)
		{
			PyErr_SetString(PyExc_TypeError, "Error when writing polar oder fourth moment(s)\n");
			return -1;
		}
	}
		//write nematic order(s)
	for (k=0; k<simulation->parameters->particle_species; k++)
	{
		if ((write(myfile, (simulation->observables->nematic_order+k), sizeof(gdouble)))==-1)
		{
			PyErr_SetString(PyExc_TypeError, "Error when writing nematic order(s)\n");
			return -1;
		}
	}
		//write first moment of nematic order(s)
	for (k=0; k<simulation->parameters->particle_species; k++)
	{
		if ((write(myfile, (simulation->observables->nematic_order_moment1+k), sizeof(gdouble)))==-1)
		{
			PyErr_SetString(PyExc_TypeError, "Error when writing nematic order first moment(s)\n");
			return -1;
		}
	}
		//write secondmoment of nematic order(s)
	for (k=0; k<simulation->parameters->particle_species; k++)
	{
		if ((write(myfile, (simulation->observables->nematic_order_moment2+k), sizeof(gdouble)))==-1)
		{
			PyErr_SetString(PyExc_TypeError, "Error when writing nematic order second moment(s)\n");
			return -1;
		}
	}
		//write third moment of nematic order(s)
	for (k=0; k<simulation->parameters->particle_species; k++)
	{
		if ((write(myfile, (simulation->observables->nematic_order_moment3+k), sizeof(gdouble)))==-1)
		{
			PyErr_SetString(PyExc_TypeError, "Error when writing nematic order third moment(s)\n");
			return -1;
		}
	}
		//write fourth moment of nematic order(s)
	for (k=0; k<simulation->parameters->particle_species; k++)
	{
		if ((write(myfile, (simulation->observables->nematic_order_moment4+k), sizeof(gdouble)))==-1)
		{
			PyErr_SetString(PyExc_TypeError, "Error when writing nematic order fourth moment(s)\n");
			return -1;
		}
	}
		//write first moment of neighbor number(s)
	for (k=0; k<simulation->parameters->particle_species; k++)
	{
		if ((write(myfile, (simulation->observables->neighbors_moment1+k), sizeof(gdouble)))==-1)
		{
			PyErr_SetString(PyExc_TypeError, "Error when writing neighbors first moment(s)\n");
			return -1;
		}
	}
		//write second moment of neighbor number(s)
	for (k=0; k<simulation->parameters->particle_species; k++)
	{
		if ((write(myfile, (simulation->observables->neighbors_moment2+k), sizeof(gdouble)))==-1)
		{
			PyErr_SetString(PyExc_TypeError, "Error when writing neighbors second moment(s)\n");
			return -1;
		}
	}
		//write third moment of neighbor number(s)
	for (k=0; k<simulation->parameters->particle_species; k++)
	{
		if ((write(myfile, (simulation->observables->neighbors_moment3+k), sizeof(gdouble)))==-1)
		{
			PyErr_SetString(PyExc_TypeError, "Error when writing neighbors third moment(s)\n");
			return -1;
		}
	}
		//write fourth moment of neighbor number(s)
	for (k=0; k<simulation->parameters->particle_species; k++)
	{
		if ((write(myfile, (simulation->observables->neighbors_moment4+k), sizeof(gdouble)))==-1)
		{
			PyErr_SetString(PyExc_TypeError, "Error when writing neighbors fourth moment(s)\n");
			return -1;
		}
	}
		//write number of measurement time steps
	if ((write(myfile, &(simulation->observables->measurement_steps), sizeof(guint64)))==-1)
	{
		PyErr_SetString(PyExc_TypeError, "Error when writing number of measurement steps\n");
		return -1;
	}
		//write number of histogram bins for orientation
	if ((write(myfile, &(simulation->observables->binnum_theta), sizeof(gint32)))==-1)
	{
		PyErr_SetString(PyExc_TypeError, "Error when writing number histogram bins for orientation\n");
		return -1;
	}
		//write number of histogram bins for neighbor number
	if ((write(myfile, &(simulation->observables->binnum_neighbors), sizeof(gint32)))==-1)
	{
		PyErr_SetString(PyExc_TypeError, "Error when writing number histogram bins for neighbor number\n");
		return -1;
	}
		//write histogram(s) for orientations
	for (k=0; k<simulation->parameters->particle_species; k++)
		for (j=0; j<simulation->observables->binnum_theta; j++)
			if ((write(myfile, (*(simulation->observables->hist_theta+k)+j), sizeof(guint64)))==-1)
			{
				PyErr_SetString(PyExc_TypeError, "Error when writing orientation histogram\n");
				return -1;
			}
		//write histogram(s) for neighbor numbers
	for (k=0; k<simulation->parameters->particle_species; k++)
		for (j=0; j<simulation->observables->binnum_neighbors; j++)
			if ((write(myfile, (*(simulation->observables->hist_neighbors+k)+j), sizeof(guint64)))==-1)
			{
				PyErr_SetString(PyExc_TypeError, "Error when writing neighbor histogram\n");
				return -1;
			}
	close(myfile);
	return 1;
}

