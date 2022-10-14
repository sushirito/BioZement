#include "global.h"

//-----------------------------------
// |0|1|2 |3 |4 |f|f|f....f|f|n1|-4|-3|-2|-1| N_MAX[0]
// |p|w|nb|n0|n1|f|f|f....f|f|n1|n0|nb|w |p |
//                              | 3| 2| 1| 0| outlet_x
//-----------------------------------
void bc_3D_copy_gradient(int nr_nodes, int *node_list, struct Field *field, struct System *sys, int nb_dir, int n1_dir)
{
  struct Step *step = sys->step;
  real *f, *f_tmp; /* Pointer to the first species->f for a phase  */

  for( int cc = 0; cc < field->n_c; ++cc) { /* loop species */
    f = &field->f[step->f_phase*cc];
    f_tmp = &field->f_tmp[step->f_phase*cc];
    for( int cn = 0; cn < nr_nodes; ++cn ) { /* loop nodes */
      int n0 = node_list[cn];
      int nb = n0 - step->n_rel_neig[nb_dir];
      int n1 = n0 - step->n_rel_neig[n1_dir];
      for( int ck = 0; ck < sys->n_Q; ++ck )  { /* loop directions */
        f[step->f_node*nb + ck]     = 2.0*f[step->f_node*n0 + ck]     - f[step->f_node*n1 + ck];
        f_tmp[step->f_node*nb + ck] = 2.0*f_tmp[step->f_node*n0 + ck] - f_tmp[step->f_node*n1 + ck];
      }
    }
  }
}

#ifdef _USE_COLOR_GRAD_
void bc_press_red_blue_fancy(int pos, double rho_boundary_red, double rho_boundary_blue, int nx, struct Node *node, struct Field *fluid, struct System *sys)
{
  int cc, ck, cd;
  int *ev_ind = sys->ev_ind;
  double *ev = sys->ev;
  double *w = sys->w;

  struct Step *step = sys->step;
  double *f_red, *f_blue;//, *f_green;
  double *f_tmp_red, *f_tmp_blue;
  double *rho_red, *rho_blue;
  double f_ck;
  //  double u[sys->n_D];
  double u[3];
  double rho;
  double w_tot;
  int pos_neig; /* Position of the nearest neighbor in the normal direction */
  int is_wall_node;
  int cn_neig;

  pos_neig = pos + nx;


  cc = 0;
  f_red = &fluid->f[step->f_phase*cc];
  f_tmp_red = &fluid->f_tmp[step->f_phase*cc];
  rho_red = &fluid->rho[step->rho_phase*cc];
  
  cc = 1;
  f_blue = &fluid->f[step->f_phase*cc];
  f_tmp_blue = &fluid->f_tmp[step->f_phase*cc];
  rho_blue = &fluid->rho[step->rho_phase*cc];

  //cc = 2;
  //f_green = &fluid->f[step->f_phase*cc];


  rho = 0;
  for (cd = 0; cd <sys->n_D; ++cd)
    u[cd] = 0;

  rho = rho_boundary_red + rho_boundary_blue;
  w_tot = 0;


  // Check if the node is a wall node
  is_wall_node = 0;
  for (ck = 0; ck < sys->n_Q; ++ck)
  {
    if (ev[ev_ind[ck]] == 0) {
      cn_neig = pos - step->n_rel_neig[ck];
      if (node->mode[cn_neig] == WALL)
        is_wall_node = 1;
    }
  }
  
  if (is_wall_node) {
    /* Calculate the velocity */
    for (ck = 0; ck < sys->n_Q; ++ck)
    {
      f_ck = f_tmp_red[step->f_node*pos_neig + ck] + f_tmp_blue[step->f_node*pos_neig + ck];
      rho += f_ck;
      for (cd = 0; cd < sys->n_D; ++cd) {
        u[cd] += ev[ev_ind[ck]+cd]*f_ck;
      }
    }
    for (cd = 0; cd <sys->n_D; ++cd) {
      u[cd] /= (rho_boundary_red + rho_boundary_blue);
      fluid->u[step->u_node*pos + cd] = u[cd];            
    }

    rho = 0;
    for (ck = 0; ck < sys->n_Q; ++ck)
    {
      //f_ck = f_tmp_red[step->f_node*pos_neig + ck] + f_tmp_blue[step->f_node*pos_neig + ck];
      f_ck = f_tmp_red[step->f_node*pos + ck] + f_tmp_blue[step->f_node*pos + ck];

      if ( ev[ev_ind[ck]]*nx < 0 ) {
        for (cd = 0; cd < sys->n_D; ++cd) {
          w_tot += 6*w[ck]*ev[ev_ind[ck]+cd]*u[cd];
        }
        rho += 2*f_ck;
      }
      if (ev[ev_ind[ck]] == 0) {
        for (cd = 0; cd < sys->n_D; ++cd) {
          w_tot += 6*w[ck]*ev[ev_ind[ck]+cd]*u[cd];
        }
        rho += f_ck;
      }
    }
    if (w_tot < -0.2) {
      rho = 1;
    } else {
      rho /= (1 + w_tot);
    }

    if (rho_boundary_red > rho_boundary_blue) {
      if (rho > rho_boundary_red)
        rho = rho_boundary_red;
      rho_red[pos] = rho;
      rho_blue[pos] = 0.0;
    } else {
      if (rho > rho_boundary_blue)
        rho = rho_boundary_blue;
      rho_red[pos] = 0.0;
      rho_blue[pos] = rho;
    }    
  }
  
  if (!is_wall_node) {
    for (ck = 0; ck < sys->n_Q; ++ck)
    {
      //f_ck = f_tmp_red[step->f_node*pos_neig + ck] + f_tmp_blue[step->f_node*pos_neig + ck];
      f_ck = f_tmp_red[step->f_node*pos + ck] + f_tmp_blue[step->f_node*pos + ck];

      if ( ev[ev_ind[ck]]*nx < 0 ) {
        w_tot += w[ck];
        rho -= 2*f_ck;
      }
      if (ev[ev_ind[ck]] == 0) {
        rho -= f_ck;
      }
    }
    
    u[0] = nx*rho/(6.0*(rho_boundary_red + rho_boundary_blue)*w_tot);
    
    for (cd = 0; cd <sys->n_D; ++cd) {
      fluid->u[step->u_node*pos + cd] = u[cd];
    }
    
    rho_red[pos] = rho_boundary_red;
    rho_blue[pos] = rho_boundary_blue;
    
  }
  for (ck = 0; ck < sys->n_Q; ++ck)
  {
    if (rho_boundary_red > rho_boundary_blue) {
      f_red[step->f_node*pos + ck] = calc_feq(ck, rho_red[pos], u, sys);
      f_blue[step->f_node*pos + ck] = f_tmp_blue[step->f_node*pos + sys->k_bb[ck]]; //calc_feq(ck, rho_blue[pos], u, sys);
    } else {
      f_red[step->f_node*pos + ck] = f_tmp_red[step->f_node*pos + sys->k_bb[ck]];
      f_blue[step->f_node*pos + ck] = calc_feq(ck, rho_blue[pos], u, sys);
    }
    //if ( (nx > 0) && (sys->t > 100000)) {
    //f_green[step->f_node*pos + ck] = 0.1*f_red[step->f_node*pos + ck];
    //}
  }
}



void bc_press_red_blue(int pos, double rho_boundary_red, double rho_boundary_blue, int nx, struct Node *node, struct Field *fluid, struct System *sys)
{
  int cc, ck, cd;
  int *ev_ind = sys->ev_ind;
  double *ev = sys->ev;
  struct Step *step = sys->step;
  double *f_red, *f_blue;//, *f_green;
  double *f_tmp_red, *f_tmp_blue;
  double *rho_red, *rho_blue;
  double f_ck;
  double u[3];
  double rho;
  int pos_neig; /* Position of the nearest neighbor in the normal direction */

  pos_neig = pos + nx;

  cc = 0;
  f_red = &fluid->f[step->f_phase*cc];
  f_tmp_red = &fluid->f_tmp[step->f_phase*cc];
  rho_red = &fluid->rho[step->rho_phase*cc];
  
  cc = 1;
  f_blue = &fluid->f[step->f_phase*cc];
  f_tmp_blue = &fluid->f_tmp[step->f_phase*cc];
  rho_blue = &fluid->rho[step->rho_phase*cc];

  //cc = 2;
  //f_green = &fluid->f[step->f_phase*cc];


  rho = 0;
  for (cd = 0; cd <sys->n_D; ++cd)
    u[cd] = 0;

  for (ck = 0; ck < sys->n_Q; ++ck)
    {
      f_ck = f_tmp_red[step->f_node*pos_neig + ck] + f_tmp_blue[step->f_node*pos_neig + ck];
      rho += f_ck;
      for (cd = 0; cd < sys->n_D; ++cd) {
	u[cd] += ev[ev_ind[ck]+cd]*f_ck;
      } 
    }
  for (cd = 0; cd <sys->n_D; ++cd) {
    u[cd] /= (rho_boundary_red + rho_boundary_blue);
    fluid->u[step->u_node*pos + cd] = u[cd];
  }

  rho_red[pos] = rho_boundary_red;
  rho_blue[pos] = rho_boundary_blue;

  for (ck = 0; ck < sys->n_Q; ++ck)
    {
      f_red[step->f_node*pos + ck] = calc_feq(ck, rho_boundary_red, u, sys);
      f_blue[step->f_node*pos + ck] = calc_feq(ck, rho_boundary_blue, u, sys);
      //if ( 1 ) 
      //  f_green[step->f_node*pos + ck] = 0.1*f_red[step->f_node*pos + ck];
    }

}

void bc_press_red_blue_z(int pos, double rho_boundary_red, double rho_boundary_blue, int nz, struct Field *fluid, struct System *sys)
{
  int cc, ck, cd;
  int *ev_ind = sys->ev_ind;
  double *ev = sys->ev;
  struct Step *step = sys->step;
  double *f_red, *f_blue;
  double *f_tmp_red, *f_tmp_blue;
  double f_ck;
  double u[3];
  double rho;
  int pos_neig = pos + nz*sys->max_n[0]*sys->max_n[1];

  cc = 0;
  f_red = &fluid->f[step->f_phase*cc];
  f_tmp_red = &fluid->f_tmp[step->f_phase*cc];
  cc = 1;
  f_blue = &fluid->f[step->f_phase*cc];
  f_tmp_blue = &fluid->f_tmp[step->f_phase*cc];

  rho = 0;
  for (cd = 0; cd <sys->n_D; ++cd)
    u[cd] = 0;

  for (ck = 0; ck < sys->n_D; ++ck)
    {
      f_ck = f_tmp_red[step->f_node*pos_neig] + f_tmp_blue[step->f_node*pos_neig];
      rho += f_ck;
      for (cd = 0; cd < sys->n_D; ++cd) {
	u[cd] += ev[ev_ind[ck]+cd]*f_ck;
      } 
      f_red[step->f_node*pos] = calc_feq(ck, rho_boundary_red, u, sys);
      f_blue[step->f_node*pos] = calc_feq(ck, rho_boundary_blue, u, sys);
    }

}

/* ---------------------------------------------------------------------------------  bc_press_f */
void bc_press_f_hack(int pos, double rho_boundary, int nx, struct Field *fluid, struct System *sys)
/*
 * bc_press_f : First versiion of the pressure boundary condtions. We will assume that the
 *              boundary normal is in either x or -x direction.
 *              NB: The neghboring nodes to the inlet/outlet in the x-direction and 
 *                  the inlet/outlet node must have the same geoemtry.
 *              Here we assume a single phase simulation
 *
 * pos : position of the boundary node
 * rho_boundary : Set boundary pressure
 * nx : boundary normal (left hand boundary nx = 1 right hand nx = -1)   
 * fluid : 
 * sys :
 */
{
  int pos_neig = pos + nx; /* Position of the nearest neighbor in the normal direction */
  //double rho_neig_normal, vel_neig_normal;
  real *f_tmp, *f_tmp_neig, *rho, *rho_neig,  *u, *u_neig; /* Pointers to the distribution, temprorar distribution, density and velocity */
  struct Step *step = sys->step;
  int alpha, alpha_rev, cd; /* counters */
  double df;
  int *ev_ind = sys->ev_ind;
  double *ev = sys->ev;
  double *w = sys->w;

  f_tmp = &fluid->f_tmp[step->f_node*pos]; 
  f_tmp_neig = &fluid->f_tmp[step->f_node*pos_neig]; 
  rho = &fluid->rho[pos];
  rho_neig = &fluid->rho[pos_neig];
  u = &fluid->u[step->u_node*pos];
  u_neig = &fluid->u[step->u_node*pos_neig];

  /* Calculate rho and vel for neig_normal */
  /* -- initiate to zero */
  (*rho_neig) = 0;
  for (cd = 0; cd < sys->n_D; ++cd)
    u_neig[cd] = 0;

  /* Calcualte the macroscopic values at the neighboring site */
  for (alpha = 0; alpha < sys->n_Q; ++alpha) 
    {
      (*rho_neig) += f_tmp_neig[alpha];
      for (cd = 0; cd < sys->n_D; ++cd)
	u_neig[cd] += ev[ev_ind[alpha] + cd]*f_tmp_neig[alpha];
    }
  for (cd = 0; cd < sys->n_D; ++cd) 
    {
      u[cd] = u_neig[cd]/rho_boundary;
      u_neig[cd] /= (*rho_neig);
    }

  /* Boundary condition implementation */
  (*rho) = 0;
  for (alpha = 0; alpha < sys->n_Q; ++alpha)
    {
      if (ev[ev_ind[alpha]]*nx > 0 ) { 
	/*if alpha is in direction of boundary normal */
	alpha_rev = sys->k_bb[alpha];
	df = (f_tmp_neig[alpha] - calc_feq(alpha, *rho_neig, u_neig, sys)) - (f_tmp_neig[alpha_rev] - calc_feq(alpha_rev, *rho_neig, u_neig, sys));
	f_tmp[alpha] = calc_feq(alpha, rho_boundary, u, sys) + f_tmp_neig[alpha_rev] - calc_feq(alpha_rev, rho_boundary, u, sys) + df;
	
      }
      (*rho) += f_tmp[alpha];      
    }

  /* Ensure that the value at the node is 'rho_boundary' */
  for (alpha = 0; alpha < sys->n_Q; ++alpha)
    {
      f_tmp[alpha] += w[alpha]*(rho_boundary - (*rho));
    }
  (*rho) = rho_boundary;
}

/* ---------------------------------------------------------------------------------  bc_press_f Z */
void bc_press_f_hack_z(int pos, double rho_boundary, int nz, struct Field *fluid, struct System *sys)
/*
 * bc_press_f : First versiion of the pressure boundary condtions. We will assume that the
 *              boundary normal is in either x or -x direction.
 *              NB: The neghboring nodes to the inlet/outlet in the x-direction and 
 *                  the inlet/outlet node must have the same geoemtry.
 *              Here we assume a single phase simulation
 *
 * pos : position of the boundary node
 * rho_boundary : Set boundary pressure
 * nx : boundary normal (left hand boundary nx = 1 right hand nx = -1)   
 * fluid : 
 * sys :
 */
{
  int pos_neig = pos + nz*sys->max_n[0]*sys->max_n[1]; /* Position of the nearest neighbor in the normal direction */
  //  double rho_neig_normal, vel_neig_normal;
  real *f_tmp, *f_tmp_neig, *rho, *rho_neig,  *u, *u_neig; /* Pointers to the distribution, temprorar distribution, density and velocity */
  struct Step *step = sys->step;
  int alpha, alpha_rev, cd; /* counters */
  double df;
  int *ev_ind = sys->ev_ind;
  double *ev = sys->ev;
  double *w = sys->w;

  f_tmp = &fluid->f_tmp[step->f_node*pos]; 
  f_tmp_neig = &fluid->f_tmp[step->f_node*pos_neig]; 
  rho = &fluid->rho[pos];
  rho_neig = &fluid->rho[pos_neig];
  u = &fluid->u[step->u_node*pos];
  u_neig = &fluid->u[step->u_node*pos_neig];

  /* Calculate rho and vel for neig_normal */
  /* -- initiate to zero */
  (*rho_neig) = 0;
  for (cd = 0; cd < sys->n_D; ++cd)
    u_neig[cd] = 0;

  /* Calcualte the macroscopic values at the neighboring site */
  for (alpha = 0; alpha < sys->n_Q; ++alpha) 
    {
      (*rho_neig) += f_tmp_neig[alpha];
      for (cd = 0; cd < sys->n_D; ++cd)
	u_neig[cd] += ev[ev_ind[alpha] + cd]*f_tmp_neig[alpha];
    }
  for (cd = 0; cd < sys->n_D; ++cd) 
    {
      u[cd] = u_neig[cd]/rho_boundary;
      u_neig[cd] /= (*rho_neig);
    }

  /* Boundary condition implementation */
  (*rho) = 0;
  for (alpha = 0; alpha < sys->n_Q; ++alpha)
    {
      if (ev[ev_ind[alpha] + 2]*nz > 0 ) { 
	/*if alpha is in direction of boundary normal */
	alpha_rev = sys->k_bb[alpha];
	df = (f_tmp_neig[alpha] - calc_feq(alpha, *rho_neig, u_neig, sys)) - (f_tmp_neig[alpha_rev] - calc_feq(alpha_rev, *rho_neig, u_neig, sys));
	f_tmp[alpha] = calc_feq(alpha, rho_boundary, u, sys) + f_tmp_neig[alpha_rev] - calc_feq(alpha_rev, rho_boundary, u, sys) + df;
	
      }
      (*rho) += f_tmp[alpha];      
    }

  /* Ensure that the value at the node is 'rho_boundary' */
  for (alpha = 0; alpha < sys->n_Q; ++alpha)
    {
      f_tmp[alpha] += w[alpha]*(rho_boundary - (*rho));
    }
  (*rho) = rho_boundary;
}
#endif  // _USE_COLOR_GRAD_

/* ---------------------------------------------------------------------------------  bc_init_run_g */
void bc_init_run_g_newconc(struct Field *species, struct Field *fluid, struct System *sys)
/*
bc_init_run_g :
	boudary conditions for the chemical components set before 
	the main time loop.

  INPUT  : void
  OUTPUT : void
*/
{
  periodic_bc(sys->node, species, sys);
  wall_bc_bounce_back(sys->node, species, sys);
  
}



/* ---------------------------------------------------------------------------------  bc_init_run_g */
void bc_init_run_g(struct InitChem *ICS, struct Field *species, struct Field *fluid,
    struct System *sys, struct Minerals *minerals, struct Boundary *bndry, struct splayTree **st_lst)
/*
bc_init_run_g :
	boudary conditions for the chemical components set before 
	the main time loop.

  INPUT  : void
  OUTPUT : void
*/
{
#ifndef BB_INLET
  periodic_bc(sys->node, species, sys);
#endif
#ifndef FLUID_OFF
  inlet_bc_g(species, fluid, sys, bndry);
#endif
  wall_bc_mid_grid_gca_ej_g_v3(ICS, species, fluid, sys, minerals, bndry, st_lst);
}


/* ---------------------------------------------------------------------------------  bc_init_run_f */
void bc_init_run_f(struct Node *node, struct Field *fluid, struct System *sys)
/*
  bc_init_run_f :
  boundary conditions for the fluid components set before
  the main time loop.
  
  INPUT  : void
  OUTPUT : void
*/
{

  periodic_bc(node, fluid, sys);
  //  inlet_bc_f(fluid, sys, bndry);
  wall_bc_bounce_back(node, fluid, sys);
}
	

/* ---------------------------------------------------------------------------------  bc_run_g */
void bc_run_g(struct InitChem *ICS, struct Field *species, struct Field *fluid,
    struct System *sys, struct Minerals *minerals, struct Boundary *bndry, struct splayTree **st_lst)
/*
bc_run_g :
	in time loop boundary conditions for chemical components.

  INPUT  : void
  OUTPUT : void
*/
{
  /* -- 1st boundary conditions TEST EJE */
  outlet_bc_3D_copy_g(species, fluid, bndry, sys);
  /* -- 2nd boundary conditions */
#ifndef BB_INLET
  periodic_bc(sys->node, species, sys);
#endif
#ifndef FLUID_OFF
  inlet_bc_g(species, fluid, sys, bndry);
#endif
  wall_bc_mid_grid_gca_ej_g_v3(ICS, species, fluid, sys, minerals, bndry, st_lst);
}


/* ---------------------------------------------------------------------------------  bc_run_f */
void bc_run_f(struct Node *node, struct Field *fluid, struct System *sys)
/*
bc_run_f :
	in time loop boundary conditions for the fluid components.

  INPUT  : void
  OUTPUT : void
*/
{
  periodic_bc(node, fluid, sys);
  wall_bc_bounce_back(node, fluid, sys);
}
