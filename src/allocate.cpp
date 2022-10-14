#include "global.h"


/* ---------------------------------------------------------------------------------  allocate */
int allocate(struct System *sys)
/*
allocate :
    allocates memory for arrays that are unrelated to the different
	distribution (ie. f and g).

 INPUT  : void

 OUTPUT : (boolean) GOOD : no error allocating memory
                    BAD  : error allocating memory
*/
{
        struct Node *node = sys->node;
        struct Step *step = sys->step;

	/*** MICROSCOPIC ***/
	/* ev, velocity basis allocation */
    //if ( !( sys->ev= (real *) calloc(sys->n_Q*sys->n_D, sizeof(real)) ) ) return BAD;
	
	if ( !( sys->ev_ind = (int *) calloc(sys->n_Q, sizeof(int)) ) ) return BAD;
	
	/* ei, direction basis allocation */
	// moved if ( !( sys->ei= (int *) calloc(sys->n_Q*sys->n_D, sizeof(int)) ) ) return BAD;

	/* directional weights */
	// moved to read_basis if ( !( sys->w= (real *) calloc(sys->n_Q, sizeof(real)) ) ) return BAD;

	/* node->mode, nodal mode */
	if ( !( node->mode= (char *) calloc(sys->max_n_tot, sizeof(char)) ) ) return BAD;

    /* node->Rmask, soruce on/off AMIN */
    if ( !( node->Rmask= (char *) calloc(sys->max_n_tot, sizeof(char)) ) ) return BAD;

	/* node->Bmask, soruce on/off AMIN */
    if ( !( node->Bmask= (char *) calloc(sys->max_n_tot, sizeof(char)) ) ) return BAD;
	
    /* biiozement inert wall*/
    if ( !( node->inert= (char *) calloc(sys->max_n_tot, sizeof(char)) ) ) return BAD;

#ifdef _MPI_ 
	/* sys->mpi->is_ghost, if node is ghost (1) or not (0) */
	if ( !( node->is_ghost= (char *) calloc(sys->max_n_tot, sizeof(char)) ) ) return BAD;
	if ( !( node->is_ghost_or_periodic= (char *) calloc(sys->max_n_tot, sizeof(char)) ) ) return BAD;
#endif
	/*** NUMERICS ***/
	if ( !( step->f_rel_neig = (int *) calloc(sys->n_Q, sizeof(int)) )  ) return BAD;
	if ( !( step->n_rel_neig  = (int *) calloc(sys->n_Q, sizeof(int)) )  ) return BAD;
	
	/* Bounce back directions */
	// moved if ( !( sys->k_bb = (int *) calloc(sys->n_Q, sizeof(int)) )  ) return BAD;


	/* TEMP VARIABLES */
	/* MOVED TO allocate_fg() below */
	/* /\* f_nc_tmp, local distirubtion function for one phase allocation *\/ */
	/* if ( !( field->nc_tmp = (real *) calloc(sys->n_Q, sizeof(real)) ) )  return BAD; */

	/* /\* field->u_eq, local equilibrium velocity *\/ */
	/* if ( !( field->u_eq = (real *) calloc(sys->n_D, sizeof(real)) ) )  return BAD; */

	return GOOD;
}


/* ---------------------------------------------------------------------------------  allocate_fg */
int allocate_fg(struct Field *field, struct System *sys)
/*
allocate_fg :
	allocates memory for arrays that are related to velocity distributions

  INPUT : 1) fg     : pointer to distribution
          2) fg_tmp : pointer to temprorary distribution
		  3) rho    : pointer to density array
		  4) u      : pointer to velocity array
		  5) n_c    : number of phases

  OUTPUT : (boolean) GOOD : no error allocating memory
                     BAD  : error allocating memory

*/
{
  /*** TEMP VARIABLES ***/
  /* f_nc_tmp, local distirubtion function for one phase allocation */
  if ( !( field->nc_tmp = (real *) calloc(sys->n_Q, sizeof(real)) ) )  return BAD;
  
  /* field->u_eq, local equilibrium velocity */
  if ( !( field->u_eq = (real *) calloc(sys->n_D, sizeof(real)) ) )  return BAD;
  
  /*** MICROSCOPIC ***/
  /* g, distribution function allocation */
  if ( !( field->f = (real *) calloc(sys->max_n_tot * field->n_c * sys->n_Q, sizeof(real)) ) )  return BAD;
  
  
  /*** MACROSCOPIC ***/
  /* rho, density allocation */
  if ( !( field->rho =     (real *) calloc(sys->max_n_tot * field->n_c, sizeof(real)) )  ) return BAD;
  // used in convergence test
  if ( !( field->rho_old = (real *) calloc(sys->max_n_tot * field->n_c, sizeof(real)) ) )  return BAD;

  /* u, velocity allocation */
  if ( !( field->u =     (real *) calloc(sys->max_n_tot * field->n_c * sys->n_D, sizeof(real)) )  ) return BAD;
  /* u_old, velocity convergence, ONLY X-COMPONENT */
  if ( !( field->u_old = (real *) calloc(sys->max_n_tot * field->n_c           , sizeof(real)) )  ) return BAD;
  

  /*** TEMP VARAIBLES ***/
  /* f_tmp, distribution function allocation */
  if ( !( field->f_tmp = (real *) calloc(sys->max_n_tot * field->n_c * sys->n_Q, sizeof(real)) ) )  return BAD;
  
  // effluent concentrations for each species
  if ( !( field->effluent_conc = (real *) calloc(field->n_c, sizeof(real)) )  ) return BAD;

  return GOOD;
}


/* ---------------------------------------------------------------------------------  allocate_wall_attr_real */
int allocate_wall_attr_real(real **attr, int n_c, int nr_wall)
/*
allocate_wall_attr_real :
	allocates memory for the wall attributes of type 'real'

  INPUT  : 1) attr : pointer to wall attribute
		   2) n_c  : number of attributes per wall site

  OUTPUT : (boolean) GOOD : no error allocating memory
                     BAD  : error allocating memory
*/
{
	if ( !( *attr = (real *) calloc(nr_wall*n_c, sizeof(real)) ) )  return BAD;

	return GOOD;
}


/* ---------------------------------------------------------------------------------  allocate_wall_attr_int */
int allocate_wall_attr_int(int **attr, int n_c, int nr_wall)
/*
allocate_wall_attr_int :
	allocates memory for the wall attributes of type 'int'

  INPUT  : 1) attr : pointer to wall attribute
		   2) n_c  : number of attributes per wall site

  OUTPUT : (boolean) GOOD : no error allocating memory
                     BAD  : error allocating memory
*/
{
	if ( !( *attr = (int *) calloc(nr_wall*n_c, sizeof(int)) ) )  return BAD;

	return GOOD;
}


/* ---------------------------------------------------------------------------------  allocate_wall_tmp_real */
int allocate_wall_tmp_real(real **attr, int n_c)
/*
allocate_wall_tmp_real :
	allocates memory for the 'real' tempororay wall attributes

  INPUT  : 1) attr : pointer to wall attribute
		   2) n_c  : number of attributes per wall site

  OUTPUT : (boolean) GOOD : no error allocating memory
                     BAD  : error allocating memory
*/
{
	if ( !( *attr = (real *) calloc(n_c, sizeof(real)) ) )  return BAD;

	return GOOD;
}


/* ---------------------------------------------------------------------------------  allocate_wall_tmp_int */
int allocate_wall_tmp_int(int **attr, int n_c)
/*
allocate_wall_tmp_int :
	allocates memory for the 'int' tempororay wall attributes

  INPUT  : 1) attr : pointer to wall attribute
		   2) n_c  : number of attributes per wall site

  OUTPUT : (boolean) GOOD : no error allocating memory
                     BAD  : error allocating memory
*/
{
	if ( !( *attr = (int *) calloc(n_c, sizeof(int)) ) )  return BAD;

	return GOOD;
}


/* ---------------------------------------------------------------------------------  allocate_vel */
int allocate_vel(real **u, int max_n_tot, int n_D)
/*
allocate_vel :
	allocates a velocity vector.

  INPUT  : 1) u : pointer to a distribtuion 

  OUTPUT : (boolean) GOOD : no error allocating memory
                     BAD  : error allocating memory

*/
{
	if ( !( *u = (real *) calloc(max_n_tot*n_D, sizeof(real)) )  ) return BAD;

	return GOOD;
}


/* ---------------------------------------------------------------------------------  allocate_rho */
int allocate_rho(real **rho, int max_n_tot)
/*
allocate_rho :
	allocates a density vector.

  INPUT  : 1) rho : pointer to a distribtuion 

  OUTPUT : (boolean) GOOD : no error allocating memory
                     BAD  : error allocating memory

*/
{
	if ( !( *rho = (real *) calloc(max_n_tot, sizeof(real)) )  ) return BAD;

	return GOOD;
}


/* ---------------------------------------------------------------------------------  allocate_phase_forces */
//int allocate_phase_forces(struct Field *fluid, struct System *sys)
///*
//allocate_phase_forces :
//	allocates memory for phase forces, including wetability (i.e. rock-fluid forces )
//
//  INPUT  : 1) flforces : pointer to fluid-fluid forces
//		   2) glforces : pointer to fluid-rock forces
//		   3) n_c      : number of phases
//
//  OUTPUT : (boolean) GOOD : no error allocating memory
//                     BAD  : error allocating memory
//*/
//{
//	if ( !( fluid->flforce = (real *) calloc(sys->max_n_tot * sys->n_D * fluid->n_c, sizeof(real)) )  ) return BAD;
//
//	if ( !( fluid->glforce = (real *) calloc(sys->max_n_tot * sys->n_D * fluid->n_c, sizeof(real)) )  ) return BAD;
//
//	return GOOD;
//}

/* ---------------------------------------------------------------------------------  allocate_outlet */ 
int allocate_outlet(struct Boundary *bndry, struct System *sys)
/*
allocate_outlet :
	allocates memory for directions in to and out of the outlet nodes and 
	initiates the directions. 
	Assumes that the outlet plane has an normal in the positive x-direction

  INPUT  : void

  OUTPUT : (boolean) GOOD : no error allocating memory
                     BAD  : error allocating memory
*/
{
	int ck;  /* Counter direction */
	int cc;  /* Counter */

	/* Number of direction */
	bndry->n_dir_inout = 0;
	for( ck = 0; ck < sys->n_Q; ++ck )
		if( sys->ei[sys->n_D*ck] >= 0 )
			bndry->n_dir_inout += 1;

	/* Allocate memory */ 	
	if( !( bndry->dir_in  = (int *) calloc(bndry->n_dir_inout, sizeof(int)) )  ) return BAD;
	if( !( bndry->dir_out = (int *) calloc(bndry->n_dir_inout, sizeof(int)) )  ) return BAD;

	/* Set in - directions */
	cc = 0;
	for( ck = 0; ck < sys->n_Q; ++ck )
		if( sys->ei[sys->n_D*ck] >= 0 )
			bndry->dir_in[cc++] = ck;

	/* Set out - directions */
	cc = 0;
	for( ck = 0; ck < sys->n_Q; ++ck )
		if( sys->ei[sys->n_D*ck] <= 0 )
			bndry->dir_out[cc++] = ck;



	return GOOD;
}


//int allocate_boundary_fluxes(struct Boundary *bndry, struct Field *species, struct Node *node)
///*
//allocate_boundary_fluxes :
//	allocates memory for the flux counters at the boundaries of the system:
//
//  INPUT  : void
//
//  OUTPUT : (boolean) GOOD : no error allocating memory
//                     BAD  : error allocating memory
//*/
//{
//	/* Allocate memory */
//	if( !( bndry->outlet_flux  = (real *) calloc(species->n_c*node->nr_out, sizeof(real)) )  ) return BAD;
//	if( !( bndry->inlet_flux   = (real *) calloc(species->n_c*node->nr_in, sizeof(real)) )  )  return BAD;
//
//	return GOOD;
//}
