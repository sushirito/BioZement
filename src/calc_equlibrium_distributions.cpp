#include "global.h"


/* ---------------------------------------------------------------------------------  calc_feq */
double calc_feq(int ck, real rho, real *u, struct System *sys)
/* 
calc_feq: 
	Calculates the velocity equilibrium distributions for a fluid node

  INPUT  : 1) ck  : direction
		   2) rho : density
		   3) u   : velocity

  OUTPUT : (real) distribution in the 'ck' direction
*/
{
  int cd; /* Counter: node, phse, direction, dimension */
  int ev_ind; /* velocity basis index */
  real u_ev, u_u; /* different velocity products */
  
  ev_ind = sys->step->stp_e_k*ck;
	
  /* calculate square of velocity */
  u_u = 0; 
  for( cd = 0; cd < sys->n_D; ++cd) {
    u_u += u[cd]*u[cd]; 
  }
  
  /* calculate vector product */
  u_ev = 0;
  for( cd = 0; cd < sys->n_D; ++cd ) 					
    u_ev += u[cd]*sys->ev[ev_ind + cd];
  
  /* Calculate equilibrium distribution */
  return rho*sys->w[ck]*(1. + sys->eq_cst[0]*u_ev + sys->eq_cst[1]*u_ev*u_ev + sys->eq_cst[2]*u_u);
}



/* ---------------------------------------------------------------------------------  calc_geq */
double calc_geq(int ck, real rho, real *u, struct System *sys)
/* 
calc_geq:
	calculates the equlibrium distribution for the diffusive components

  INPUT : 1) ck  : direction
          2) rho : node density
		  3) u   : velocity

  OUTPUT : (real) equlibrium for the given direction 'ck'
*/
{
	int cd;     /* Counter: node, phse, direction, dimension */
	int ev_ind; /* velocity basis index */
	real u_ev;  /* different velocity products */

#ifdef FLUID_OFF
	return rho*sys->w[ck];
#endif
	ev_ind = sys->step->stp_e_k*ck;

	/* calculate vector product */
	u_ev = 0;
	for( cd = 0; cd < sys->n_D; ++cd ) 					
	  u_ev += u[cd]*sys->ev[ev_ind + cd];
	
	
	/* Calculate equilibrium distribution */
	return rho*sys->w[ck]*(1. + sys->eq_cst[0]*u_ev);
}

