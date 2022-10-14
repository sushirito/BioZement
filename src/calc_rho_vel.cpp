#include "global.h"

#ifdef _USE_COLOR_GRAD_
void correction_rho_vel_color_grad_fluid(Node *node, Field *fluid, System *sys)
/* Add corrections to macroscopic varaibles.
   NB: Should be called after the calc_rho_vel procedures for the fluid and speceis.
       But before the collision.
 */
{
  int u_ind, f_ind; /* index in for the velocity and velocity distribution */
  int cn_neig;
  struct Step *step = sys->step;
  double *u_tot, rho_tot;
  double *F;
  double max_vel = 0.0;
  double vel_norm;
  double delta_color_wall = 0.0;
  double delta_w_wall = 0.0;
  double *w = sys->w;

  double sum_rho;
  double sum_rho_cc;

  sum_rho = 0.0;
  sum_rho_cc = 0;

  F = fluid->gravity;

  u_ind = 0;
  f_ind = 0;
  for (int cn = 0; cn < sys->max_n_tot; ++cn) { /* nodes */
    if ( (node->mode[cn] == FLUID) && !(node->is_ghost_or_periodic[cn])) {

      u_tot = &fluid->u_tot[u_ind];
      rho_tot = fluid->rho_tot[cn];
      /* Add external force correction */
      for (int nd = 0; nd < sys->n_D; ++nd) {
        u_tot[nd] += 0.5*F[nd]/rho_tot;
      }
    }
    /* Update indices */
    u_ind += step->u_node;
    f_ind += step->f_node;
  }
}


void correction_rho_vel_color_grad(struct Node *node, struct Field *fluid, struct Field *species, struct System *sys)
/* Add corrections to macroscopic varaibles.
   NB: Should be called after the calc_rho_vel procedures for the fluid and speceis.
       But before the collision.
 */
{
  int u_ind, f_ind; /* index in for the velocity and velocity distribution */
  int cn_neig;
  struct Step *step = sys->step;
  double *u_tot, rho_tot;
  double *F;
  double max_vel = 0.0;
  double vel_norm;
  double delta_color_wall = 0.0;
  double delta_w_wall = 0.0;
  double *w = sys->w;

  double sum_rho;
  double sum_rho_cc;

  sum_rho = 0.0;
  sum_rho_cc = 0;

  F = fluid->gravity;

  u_ind = 0;
  f_ind = 0;
  for (int cn = 0; cn < sys->max_n_tot; ++cn) { /* nodes */
    if ( (node->mode[cn] == FLUID) && !(node->is_ghost_or_periodic[cn])) {

      u_tot = &fluid->u_tot[u_ind];
      rho_tot = fluid->rho_tot[cn];
      /* Add external force correction */
      for (int nd = 0; nd < sys->n_D; ++nd) {
        u_tot[nd] += 0.5*F[nd]/rho_tot;
      }
      vel_norm = u_tot[0]*u_tot[0] + u_tot[1]*u_tot[1] + u_tot[2]*u_tot[2];
      if (max_vel < vel_norm) {
	max_vel = vel_norm;
      }
      /* Calculate the psi field psi = rho_species/rho_fluid_tot */
      for (int cc=0; cc < species->n_c; cc++) {
	sum_rho += species->rho[step->rho_phase*cc + cn];
	sum_rho_cc += fluid->rho[step->rho_phase * 0 + cn];
	species->psi[step->rho_phase*cc + cn] = species->rho[step->rho_phase*cc + cn]/rho_tot;
      }
    } else if ((node->mode[cn]) == WALL && !(node->is_ghost_or_periodic[cn])) {  /* Add wetting corrections due to species value */
      delta_color_wall = 0.0;
      delta_w_wall = 0.0;

      for (int ck = 1; ck < sys->n_Q; ++ck) { /* For each lattice direction */
	cn_neig = cn - step->n_rel_neig[ck];
	if (node->mode[cn_neig] == FLUID) {
	  delta_color_wall += w[ck]*species->rho[cn_neig]/fluid->rho_tot[cn_neig];
	  delta_w_wall += w[ck];
	}
      } /* end Lattice Direction*/
      /* change wall wettability */
      delta_color_wall /= delta_w_wall;
      fluid->rho[step->rho_phase * 0 + cn] += 0.5*delta_color_wall;
      fluid->rho[step->rho_phase * 1 + cn] -= 0.5*delta_color_wall;
    } /* end WALL */

    /* Update indices */
    u_ind += step->u_node;
    f_ind += step->f_node;
  } /* END for each node */

  printf("max_vel = %g\n", sqrt(max_vel));
  printf("sum_tot/sum_rho_cc = %g / %g\n", sum_rho, sum_rho_cc);
}

void calc_rho_vel_f_color_grad(struct Node *node, struct Field *fluid, struct System *sys)
/* Calculates the macroscopic variabels for the mulitphase color gradient method  
   - calculates the density of each phase
   - calculates the total density and savs it in field->rho_tot
   - calculates the total velocity and savs it in field->u_tot
*/
{
  int u_ind, f_ind; /* index in for the velocity and velocity distribution */
  struct Step *step = sys->step;
  int *ev_ind = sys->ev_ind;
  double *ev = sys->ev;

  double * f_tmp_cc;
  double *rho_cc, * u_cc;

  /* step counters */
  u_ind = 0;
  f_ind = 0;

  for (int cn = 0; cn < sys->max_n_tot; ++cn) { /* nodes */
    /* Init the total fluid velcocity and density */
    fluid->rho_tot[cn] = 0.0;
    for (int cd = 0; cd < sys->n_D; ++cd) 
      fluid->u_tot[u_ind + cd] = 0.0;

    for (int cc = 0; cc < fluid->n_c; ++cc) { /* Phases */
      rho_cc = &fluid->rho[step->rho_phase*cc];
      u_cc = &fluid->u[step->u_phase*cc];
      f_tmp_cc = &fluid->f_tmp[step->f_phase*cc];
      /* Calculate rho and velocity for each phase */
      rho_cc[cn] = 0.0;
      for (int cd = 0; cd < sys->n_D; ++cd) 
        u_cc[u_ind + cd] = 0.0;
      for (int ck = 0; ck < sys->n_Q; ++ck) { /* Lattice direction */
        /* increase density */
        rho_cc[cn] += f_tmp_cc[f_ind + ck];
        /* increase velocity */
        for (int cd = 0; cd < sys->n_D; ++cd) {
          u_cc[u_ind + cd] += ev[ev_ind[ck]+cd]*f_tmp_cc[f_ind + ck];
        }	
      } /* END for each lattice direction*/

      fluid->rho_tot[cn] += rho_cc[cn];
      for (int cd = 0; cd < sys->n_D; ++cd) {
	fluid->u_tot[u_ind + cd] += u_cc[u_ind + cd];
      }	

      if ( (node->mode[cn] % PERIODIC) == WALL ) {
        rho_cc[cn] = fluid->rho_wall[cc];
      }

    } /* END for each phase */

    /* Set velocity:  */
    /* NB: assume that force corrections are added in the collision function */
    if ( (node->mode[cn] == FLUID) && !(node->is_ghost_or_periodic[cn])) {
      for (int cd = 0; cd < sys->n_D; ++cd) {
	fluid->u_tot[u_ind + cd] /= fluid->rho_tot[cn];
      }
    }

    /* Update indices */
    u_ind += step->u_node;
    f_ind += step->f_node;
  } /* END for each node */
}
#endif


void calc_rho_vel_newconc(real *fg, real *rho, real *u, real *q, real *F, struct System *sys)
{
  int ck, cd; /* counters: directions, dimensions */
  int ev_ind; 
  
  /* Calculate the density and the velocity */ 
  *rho = 0; 
  for( cd = 0; cd < sys->n_D; ++cd) 
    u[cd] = 0;
  
  ev_ind = 0;
  for( ck = 0; ck < sys->n_Q; ++ck )
    {
      *rho += fg[ck];
      for (cd = 0; cd < sys->n_D; ++cd)
	u[cd] += sys->ev[ev_ind + cd]*fg[ck];
      
      ev_ind += sys->step->stp_e_k;
    }

  /* Correct the veloicty. 
   * This must be done before corrections are added
   * to the density 
   */
  for (cd = 0; cd < sys->n_D; ++cd)
    u[cd] += 0.5*F[cd];
  
  for( cd = 0; cd < sys->n_D; ++cd)
    u[cd] /= *rho;

  /* correct the density */
  *rho += 0.5*(*q);

}


/* ---------------------------------------------------------------------------------  calc_rho_vel */
void calc_rho_vel(real *fg, real *rho, real *u, struct System *sys)
/*
calc_rho_vel :
	calculates the density and velocity

  INPUT : 1) fg  : distribution at the current node
          2) rho : density to be calculated
		  3) u   : velocity to be calculated

  OUTPUT : void
	 
*/
{
	int ck, cd; /* counters: directions, dimensions */
	int ev_ind; 

	/* Calculate the density and the velocity */ 
	*rho = 0; 
	for( cd = 0; cd < sys->n_D; ++cd) 
		u[cd] = 0;
	
	ev_ind = 0;
	for( ck = 0; ck < sys->n_Q; ++ck )
	{
	  *rho += fg[ck];
	  for (cd = 0; cd < sys->n_D; ++cd)
	    u[cd] += sys->ev[ev_ind + cd]*fg[ck];

	  ev_ind += sys->step->stp_e_k;
	}

	for( cd = 0; cd < sys->n_D; ++cd)
	  u[cd] /= *rho;
}




void calc_rho_vel_loop(double *f_dist, struct Field *fluid, struct System *sys, double *q)
{
  // calculate rho and u for 0'th output-file
  for( int cc = 0; cc < fluid->n_c; ++cc ) { /* phases */
    double *f     = &f_dist    [sys->step->f_phase*cc];
    double *rho   = &fluid->rho[sys->step->rho_phase*cc];
    double *u     = &fluid->u  [sys->step->u_phase*cc];
    int f_ind = 0, u_ind = 0;
    double q = 0;
    for( int cn = 0; cn < sys->max_n_tot; ++cn) {
      if( sys->node->mode[cn]==FLUID || sys->node->mode[cn]==INLET || sys->node->mode[cn]==OUTLET ) {
        calc_rho_vel_newconc(&f[f_ind], &rho[cn], &u[u_ind], &q, fluid->gravity, sys);
      }
      u_ind += sys->step->u_node;
      f_ind += sys->step->f_node;
    }
  }

}


/* ---------------------------------------------------------------------------------  calc_rho */
void calc_rho(real *fg, real *rho, struct System *sys)
/*
calc_rho :
	calculates the density

  INPUT : 1) fg  : distribution at the current node
          2) rho : density to be calculated

  OUTPUT : void
*/
{
	int ck; /* counters: directions */

	/* Calculate the density and the velocity */ 
	*rho = 0; 
	for( ck = 0; ck < sys->n_Q; ++ck )
		*rho += fg[ck];
}


/* ---------------------------------------------------------------------------------  calc_mom */
void calc_mom(real *fg, real *p, struct System *sys)
/*
calc_vel :
	calculates the momentum at a node

  INPUT  : 1) fg : distribution at the current node
		   2) p  : momentum at the given node

  OUTPUT :	void
*/
{
	int ck, cd; /* counters: directions, dimensions */
	int ev_ind; /* index in the base velocity array */

	/* Initiate momentum */
	for( cd = 0; cd < sys->n_D; ++cd) 
		p[cd] = 0.0;

	ev_ind = 0;
	for( ck = 0; ck < sys->n_Q; ++ck )
	{
		for (cd = 0; cd < sys->n_D; ++cd)
			p[cd] += sys->ev[ev_ind + cd]*fg[ck];

		ev_ind += sys->step->stp_e_k;
	}
}



/* ---------------------------------------------------------------------------------  calc_rho_mom_all */
/* void calc_rho_mom_all(real *fg_c, real *rho_c, real *p_c, int n_c) */
void calc_rho_mom_all(struct Field *fluid, struct System *sys)
/*
calc_rho_mom_all ;
	calculates the density and momentum at all nodes

  INPUT  : 1) fg_c  : velocty density
		   2) rho_c : density (to be calculated)
		   3) p_c   : momentum (to be calculated)
		   4) n_c   : number of phases

  OUTPUT : void
*/
{
	int cc, cn, cd, ck; /* Conter phase, node, dimension, direction  */
	int ind_fg, ind_rho, ind_p, ind_ev; /* Index for a distribution, density and velocity, basis velocities */
	real *fg, *rho, *p; /* distribution, density and momentum for one phase */
	real fg_k; /* temporal distirbution for one nodal direction */

	struct Step *step = sys->step;
		
	for( cc = 0; cc < fluid->n_c; ++cc )
	{
		fg  = &fluid->f_tmp[step->f_phase*cc];
		rho = &fluid->rho[step->rho_phase*cc];
		p   = &fluid->u[step->u_phase*cc];

		ind_fg  = 0;
		ind_rho = 0;
		ind_p   = 0;

		for( cn = 0; cn < sys->max_n_tot; ++cn )
		{
			/* Initialize density and momentum */
			rho[ind_rho] = 0.0;
			for( cd = 0; cd < sys->n_D; ++cd)
				p[ind_p + cd] = 0.0;

			/* Calculate the hydrovar */
			ind_ev = 0;
			for( ck = 0; ck < sys->n_Q; ++ck )
			{
				fg_k = fg[ind_fg + ck];

				rho[ind_rho] += fg_k;
				for( cd = 0; cd < sys->n_D; ++cd)
					p[ind_p + cd] += sys->ev[ind_ev + cd]*fg_k;

				ind_ev += step->stp_e_k;
			}

			ind_fg  += step->f_node;
			ind_rho += step->rho_node;
			ind_p   += step->u_node;
		}
	}
}




/* ---------------------------------------------------------------------------------  shanchen_f */
//void shanchen_f(struct Field *fluid, struct System *sys)
///*
//shanchen_f :
//calculates the forces between phases, and phases and solid. For the fluid
//
//  INPUT  : void
//  OUTPUT : void
//  */
//{
//	int cc, cn, cd, ck; /* Coutners for phase, node, dimension and direction */
//	int ind_rho, ind_lf, ind_ei, ind_neig; /* index: density, force, position, neighbor  */
//	real *rho, *flf, *glf; /* pointer to density, force fluid-fluid, force fluid wall */
//	real norm[15] = {0., 1., 1., 1., 1., 1., 1., ISQRT3, ISQRT3, ISQRT3, ISQRT3, ISQRT3, ISQRT3, ISQRT3, ISQRT3};  /* Inverse(?) square(?) norm for each direction */
//
//	struct Step *step = sys->step;
//	struct Node *node = sys->node;
//
//	for( cc = 0; cc < fluid->n_c; ++cc )
//	{
//		flf = &fluid->flforce[step->u_phase*cc];
//		glf = &fluid->glforce[step->u_phase*cc];
//		rho = &fluid->rho[step->rho_phase*cc];
//
//		ind_lf  = 0;
//		ind_rho = 0;
//
//		for( cn = 0; cn < sys->max_n_tot; ++cn )
//		{
//			/* initiate the phase forces */
//			for( cd = 0; cd < sys->n_D; ++cd )
//				flf[ind_lf + cd] = glf[ind_lf + cd] = 0.;
//
//			if( node->mode[cn] == FLUID )
//			{
//				/* Initiation due to ck = 1 */
//				ind_ei = step->stp_e_k;
//				for( ck = 1; ck < sys->n_Q; ++ck )
//				{
//					/* NB: use '- rel_neig_ind_' as the relative position of the neighbor, this is the
//					       "inverse" of the definition in 'propagate' */
//					ind_neig = cn - step->n_rel_neig[ck];
//
//					if( node->mode[ind_neig] == FLUID || node->mode[ind_neig] == (PERIODIC + FLUID) )
//					{
//						for( cd = 0; cd < sys->n_D; ++cd )
//							flf[ind_lf + cd] += sys->ei[ind_ei + cd]*norm[ck]*rho[ind_neig];
//
//					} else if( node->mode[ind_neig] == WALL || node->mode[ind_neig] == (PERIODIC + WALL) ) {
//						for( cd = 0; cd < sys->n_D; ++cd )
//							glf[ind_lf + cd] += sys->ei[ind_ei + cd]*norm[ck];
//					}
//
//					ind_ei += step->stp_e_k;
//				}
//
//			}
//
//
//			ind_lf  += step->u_node;
//			ind_rho += step->rho_node;
//		}
//	}
//}


/* ---------------------------------------------------------------------------------  calc_rho_vel_cn_f */
/* void calc_rho_vel_cn_f(int cn, real *u_norm, real *u_eq) */
void calc_rho_vel_cn_f(int cn, real *u_norm, struct Field *fluid, struct System *sys)
/*
calc_rho_vel_cn_f :
	calculates a normalized baricentric velocity and a noralization factor

  INPUT  : 1) cn     : node number
		   2) u_norm : velocity normalization factor
		   3) u_eq   : equlibrium velocity

  OUTPUT : void
*/
{
	int cc, cd;        /* Counter phase direction */
	real tau, tau_inv; /* relaxation time, inverse of */
	real *u, rho;     /* pointer to velcity and density */

	struct Step *step = sys->step;

	/* Initiation */
	(*u_norm) = 0.0;
	for( cd = 0; cd < sys->n_D; ++cd )
		fluid->u_eq[cd] = 0.0;

	for( cc = 0; cc < fluid->n_c; ++cc )
	{
		/* Relaxation time */
		tau     = fluid->tau[cc];
		tau_inv = 1./tau;

		/* Pointer update */
		u   = &fluid->u[step->u_phase*cc + step->u_node*cn];
		rho = fluid->rho[step->rho_phase*cc + step->rho_node*cn];

		/* update the normalization factor, and equlibrium velcotiy */
		(*u_norm) += rho*tau_inv;

		for( cd = 0; cd < sys->n_D; ++cd )
			fluid->u_eq[cd] += u[cd]*tau_inv;
	}
}
