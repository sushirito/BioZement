#include "global.h"


///* ---------------------------------------------------------------------------------  inlet_bc_f */
//
//void inlet_bc_f(struct Field *fluid, struct System *sys, struct Boundary *bndry)
///*
//  inlet_bc_f:
//  sets the inlet boundary condition to a predefined density and fluid velocity
//
//  INPUT  : void
//  OUTPUT : void
//
//*/
//{
//  int cc, ci, ck;             /* Counter phase, inlet node, direction */
//  int in_ind;                 /* inlet node index (in f) */
//  real *f;                    /* pointer to the first entry of a phase in fluid->f*/
//  real rho;                   /* Initial density */
//  real *u;                    /* Initial velocity */
//
//  struct Step *step = sys->step;
//  struct Node *node = sys->node;
//
//  /* Initiate the influx */
//  bndry->J_influx = 0.0;
//
//  for( cc = 0; cc < fluid->n_c; ++cc) /* Phases */
//    {
//      f = &fluid->f[step->f_phase*cc];
//      rho = fluid->rho_inlet[cc];
//      u   = &fluid->u_inlet[sys->n_D*cc];
//
//      for( ci = 0; ci < node->nr_in; ++ci ) /* Inlet nodes */
//	{
//	  /* Update influx */
//	  bndry->J_influx += rho*u[0];
//	  /* Set the index in fluid->f for the inlet node*/
//	  in_ind = step->f_node*node->list_in[ci];
//
//	  for( ck = 0; ck < sys->n_Q; ++ck ) /* Inlet directions */
//	    f[in_ind + ck] = calc_feq(ck, rho, &u[0], sys);
//	}
//    }
//
//}
//


/* ---------------------------------------------------------------------------------  inlet_bc_g */
void inlet_bc_g(struct Field *species, struct Field *fluid, struct System *sys, struct Boundary *bndry)
/*
  inlet_bc_g:
  sets the inlet boundary condition to a predefined density
  Accounts for the diffusive flux out of the inlet side of the core

  INPUT  : void
  OUTPUT : void

*/
{
  int cc, ci, ck;         /* Counter phase, inlet node, direction */
  int nn_ind;             /* node index */
  int fn_ind, fn_nb;             /* inlet node index (in g) */
  int un_ind, un_nb;             /* veloctiy index */
  int rn_ind, rn_nb;             // rho index
  real *g;//, *f, *f_tmp;                /* pointer to the first entry of a phase in species->f*/
  real *u, *rho;                /* */
  real conc = 0. ;              /* Initial density */
#ifndef NEWCONC
  real u_mean;            /* Half way velociyt at the inlet */
#endif
  //  real rho_u_mean;            // half-way f-distribution along a link at the inlet
  
  //static int nneg = 0;
  //static int ncall = 0;
  struct Step *step = sys->step;
  struct Node *node = sys->node;
  
  //ncall++;
  
  u = &fluid->u[0];
  rho = &fluid->rho[0];
  //f     = &fluid->f    [0];
  //f_tmp = &fluid->f_tmp[0];
  
  for( cc = 0; cc < species->n_c; ++cc) { /* Phases */
    g = &species->f[step->f_phase*cc];
#ifdef BB_INLET
    conc = 0.;
#else
    conc = species->rho_inlet[cc];
#endif
    //conc = 0;
    for( ci = 0; ci < node->nr_in; ++ci ) { /* Inlet nodes */
      /* Set the index in fluid->f for the inlet node*/
      nn_ind = node->list_in[ci];
      fn_ind = step->f_node*nn_ind;
      un_ind = step->u_node*nn_ind;
      rn_ind = step->rho_node*nn_ind;
      
      for( ck = 0; ck < sys->n_Q; ++ck ) { /* Inlet directions */  
        if (node->mode[nn_ind - step->n_rel_neig[ck]] == FLUID) {
#ifdef NEWCONC
	  rn_nb = rn_ind - step->rho_node*step->n_rel_neig[ck];
	  un_nb = un_ind - step->u_node*step->n_rel_neig[ck];
	  fn_nb = fn_ind - step->f_rel_neig[ck];
	  
	  //rho_u_mean = rho[rn_ind]*u[un_ind] + rho[rn_nb]*u[un_nb];
	  //g[fn_ind + ck] = g[fn_ind - step->f_rel_neig[ck] + sys->k_bb[ck]] + sys->w[ck]*conc*(2*rho_u_mean + fluid->gravity[0])/sys->cs_2;
	  g[fn_ind + ck] = g[fn_nb + sys->k_bb[ck]] + sys->w[ck]*conc*(rho[rn_ind]*u[un_ind] + rho[rn_nb]*u[un_nb] + fluid->gravity[0])/sys->cs_2;	  
#else
          u_mean = 0.5*( u[un_ind] + u[un_ind - step->u_node*step->n_rel_neig[ck]] ); // only x-direction
	  //if (u_mean<0) nneg++;//printf("u_mean: %g\n", u_mean);// u_mean = 0.0;
          g[fn_ind + ck] = g[fn_ind - step->f_rel_neig[ck] + sys->k_bb[ck]] + 2.0*sys->w[ck]*conc*u_mean/sys->cs_2;
#endif
        } else {
          g[fn_ind + ck] = 0;
        }
      }
    }        
  }
  communicate(species->f    , species->n_c, sys);

  //if (ncall%100==0) {
  //  broadcast_sum_of_int_to_all_procs(&nneg, sys->mpi);
  //  if (sys->mpi->my_rank==0) 
  //    printf("\nneg: %d\n",nneg);
  //  nneg = 0;
  //}

}


/* ---------------------------------------------------------------------------------  inlet_bc_g */
void inlet_leaky_bc_g(struct Field *species, struct Field *fluid, struct System *sys, struct Boundary *bndry)
/*
	inlet_bc_g:
	sets the inlet boundary condition to a predefined density

	INPUT  : void
	OUTPUT : void

 */
{
  int cc, ci, ck;         /* Counter phase, inlet node, direction */
  int in_ind;             /* inlet node index (in g) */
  int un_ind;             /* veloctiy index */
  real *g;                /* pointer to the first entry of a phase in species->f*/
  real *u;                /* */
  real rho ;              /* Initial density */

  struct Step *step = sys->step;
  struct Node *node = sys->node;

  u = &fluid->u[0];

  for( cc = 0; cc < species->n_c; ++cc) /* Phases */
  {
    g = &species->f[step->f_phase*cc];
	  rho = species->rho_inlet[cc];
	  
	  for( ci = 0; ci < node->nr_in; ++ci ) /* Inlet nodes */
	  {
	    /* Set the index in fluid->f for the inlet node*/
	    in_ind = step->f_node*node->list_in[ci];
	    un_ind = step->u_node*node->list_in[ci];

	    for( ck = 0; ck < sys->n_Q; ++ck ) { /* Inlet directions */
	      g[in_ind + ck] = calc_geq(ck, rho, &u[un_ind], sys);
	    }
	  }
	  
  }
  communicate(species->f    , species->n_c, sys);

}

//-----------------------------------
//
//
//-----------------------------------
void inlet_bc_3D_copy_fg(struct Field *field, struct Boundary *bndry, struct System *sys, int direction)
{
  int cc,ci,ck,in,nb;
  struct Node *node = sys->node;
  struct Step *step = sys->step;
  real *f, *f_tmp; /* Pointer to the first species->f for a phase  */
  
  for( cc = 0; cc < field->n_c; ++cc) { /* species */
    f = &field->f[step->f_phase*cc];
    f_tmp = &field->f_tmp[step->f_phase*cc];
    for( ci = 0; ci < node->nr_in; ++ci ) { /* outlet nodes */
      in = node->list_in[ci];
      nb = in - step->n_rel_neig[direction];
      //node->mode[nb] = SOLID;
      //printf("in, nb: %d, %d\n", get_global_xcoord(in, sys), get_global_xcoord(nb, sys));
      for( ck = 0; ck < sys->n_Q; ++ck )  { /* direction */
	//f[step->f_node*in + ck] = f[step->f_node*nb + ck];
	//f_tmp[step->f_node*in + ck] = f_tmp[step->f_node*nb + ck];
	f[step->f_node*nb + ck] = f[step->f_node*in + ck];
	f_tmp[step->f_node*nb + ck] = f_tmp[step->f_node*in + ck];
      }
    }
  }
}

