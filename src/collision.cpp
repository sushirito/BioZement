#include "global.h"

#ifdef _USE_COLOR_GRAD_


void collision_f_color_grad(struct Node *node, struct Field *fluid, struct System *sys)
/*
  collision_f :
  Accomplish a collision for the fluid part of the system
  NB: - Make a seperate force calculation that modifies the macroscopic veolocty (called after
        calc_rho_vel but before collision)
      - Not included mass-source 
      - No pressure boundaries 
      - B is hard-coded in this function

      In "Three-dimensional lattice Boltzmann model for immiscible two-phase flow simulations"
          Liu, Valocchi and Kang, Phys. Rev. E 85 (4) 1-14 (2012)

       *) Equation (1) and (3) describing the color fields collison is not correct. We need 
	  to assign a effective tau (see eq. (13)) making the collision color-less.

       *) Note that we also caluclate the effecite tau to assign the correct A value in the
	  surface tension pertubation.

       *) Add for corrections to u_tot: += 0.5*F/rho_tot

     A general multiphase solver including solute diffusion and transport is described in 

     "Modelingmass transfer and reaction of dilute solutes in a ternary phase system by the lattice Boltzmannmethod"
     Fu, Bai, Luo, Jin and Cheng, Phys. Rev. E 95, 043304 (2017). 

  INPUT  : void
  OUTPUT : void
*/
{
  int u_ind, f_ind; /* index in for the velocity and velocity distribution */
  struct Step *step = sys->step;
  double A;
  double B[19] = {-1.0/3.0, 1.0/18.0, 1.0/18.0, 1.0/18.0, 1.0/18.0, 1.0/18.0, 1.0/18.0,
      1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0,
      1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0};
  int *ev_ind = sys->ev_ind;
  double *ev = sys->ev;
  double *w = sys->w;
  double C2, C4;
  double *f_cc, * f_tmp_cc;
  double *f_cc1, *f_cc2;
  double *rho_cc, * u_cc;
  double *rho_cc1, *rho_cc2;
  double *u_tot, rho_tot;
  // double rho, u[3];
  double grad[3];
  double beta;
  double cos_theta, ev_norm;
  double omega_surf;

  double rho_n, tau_eff, tau_eff_surface_tension;

  /* TRT-variables */
  double f_sym, f_ant, feq_sym, feq_ant;
  double tau_eff_sym, tau_eff_ant;
  double Lambda = 3./16.;

  /* Forcing  */
  /* NB no source term */
  double *F;  
  double cu_tmp, cF_tmp, Fu_tmp;


  C2 = sys->cs_2;
  C4 = C2*C2;

  /* Set local paramters */
  // A = fluid->surface_tension; /* Not the actuall surface tension. Needs to be scaled by a lattice dependent factor*/
  beta = fluid->beta;
  F = fluid->gravity;


  /* Calculated rho_sigma and \sum_\alpha c_\alpha*f^\sigma_\sigma */
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
      

      /* Standard collision: 
         NB only tested for two phases 
       */
      /* --  calculation of effective tau and set velocity */
      tau_eff = 0.0;
      for (int cc = 0; cc < fluid->n_c; ++cc) {
        rho_cc = &fluid->rho[step->rho_phase*cc];
        tau_eff += rho_cc[cn]*fluid->tau[cc];
      }
      tau_eff /= rho_tot;      

      /* TRT - collision times */
      tau_eff_sym = tau_eff;
      tau_eff_ant = Lambda/(tau_eff_sym - 0.5) + 0.5;

      for (int cc = 0; cc < fluid->n_c; ++cc) {
        rho_cc = &fluid->rho[step->rho_phase*cc];
        f_tmp_cc = &fluid->f_tmp[step->f_phase*cc + f_ind];
        f_cc = &fluid->f[step->f_phase*cc + f_ind];
        for (int ck = 0; ck < sys->n_Q; ++ck) {
	  f_sym = 0.5*(f_tmp_cc[ck] + f_tmp_cc[sys->k_bb[ck]]);
	  f_ant = 0.5*(f_tmp_cc[ck] - f_tmp_cc[sys->k_bb[ck]]);

	  feq_sym = 0.5*(calc_feq(ck, rho_cc[cn], u_tot, sys) + calc_feq(sys->k_bb[ck], rho_cc[cn], u_tot, sys));
	  feq_ant = 0.5*(calc_feq(ck, rho_cc[cn], u_tot, sys) - calc_feq(sys->k_bb[ck], rho_cc[cn], u_tot, sys));


          f_cc[ck] = f_tmp_cc[ck] - (f_sym - feq_sym)/tau_eff_sym - (f_ant - feq_ant)/tau_eff_ant;
        }
      }
      
      /* Color gradient implementation of surface tension */
      for (int cc1 = 0; cc1 < fluid->n_c - 1; cc1++) {
        rho_cc1 = &fluid->rho[step->rho_phase*cc1];
        f_cc1 = &fluid->f[step->f_phase*cc1 + f_ind];
        for (int cc2 = cc1+1; cc2 < fluid->n_c; cc2++) {
          rho_cc2 = &fluid->rho[step->rho_phase*cc2];
          f_cc2 = &fluid->f[step->f_phase*cc2 + f_ind];
          color_gradient(cn, grad, rho_cc1, rho_cc2, sys);

          /* Calculate surface tension NB 9/4 scaling is a lattice dependent factor*/
          rho_n = (rho_cc1[cn] - rho_cc2[cn])/(rho_cc1[cn] + rho_cc2[cn]);
          tau_eff_surface_tension = 0.5*( (1+rho_n)*fluid->tau[cc1] + (1-rho_n)*fluid->tau[cc2]);
          A = 9*fluid->surface_tension_multiphase_matrix[cc1][cc2]/(4*tau_eff_surface_tension);

          for (int ck=0; ck < sys->n_Q; ++ck) {
            omega_surf = Omega_surface_tension(ck, grad, A, B, sys);
            f_cc1[ck] += omega_surf;
            f_cc2[ck] += omega_surf;
          }
        }
      }
      
      /* Latva-Kokko Rotham recoloring scheme (D'orthona) */      
      double f_ck_tot;
      double beta_term;
      for (int ck=0; ck < sys->n_Q; ++ck) {
        /* find total f in the ck-direction */
        f_ck_tot = 0.0;
        for (int cc=0; cc < fluid->n_c; ++cc) {
          f_ck_tot += fluid->f[step->f_phase*cc + f_ind + ck];
        }
        /* Adding a body force */
	/* Becomes SRT when tau_eff_sym = tau_eff_ant*/
	f_ck_tot += DeltaOmegaF(ck, tau_eff_sym,  tau_eff_ant, u_tot, F, sys);

	for (int cc=0; cc < fluid->n_c; ++cc) {
          rho_cc = &fluid->rho[step->rho_phase*cc];
          fluid->f[step->f_phase*cc + f_ind + ck] =  rho_cc[cn]*f_ck_tot/rho_tot;
        }

        /* Redistribution */
        for (int cc1 = 0; cc1 < fluid->n_c - 1; cc1++) {
          rho_cc1 = &fluid->rho[step->rho_phase*cc1];
          f_cc1 = &fluid->f[step->f_phase*cc1 + f_ind];
          for (int cc2 = cc1+1; cc2 < fluid->n_c; cc2++) {
            rho_cc2 = &fluid->rho[step->rho_phase*cc2];
            f_cc2 = &fluid->f[step->f_phase*cc2 + f_ind];
            color_gradient(cn, grad, rho_cc1, rho_cc2, sys);
            beta_term = sys->w[ck]*beta*cos_phi(ck, grad, sys)*rho_cc1[cn]*rho_cc2[cn]/rho_tot;
            f_cc1[ck] += beta_term;
            f_cc2[ck] -= beta_term;
          }
        }
      }
    } /* If cn is FLUID*/

    /* Update indices */
    u_ind += step->u_node;
    f_ind += step->f_node;
  } /* END for each node */ 
}


void correction_rho_vel(struct Node *node, struct Field *fluid, struct System *sys)
/* Add corrections to macroscopic varaibles.
   NB: Should be called after the calc_rho_vel procedures for the fluid and speceis.
       But before the collision.
 */
{
  int u_ind, f_ind; /* index in for the velocity and velocity distribution */
  struct Step *step = sys->step;
  double *u_tot, rho_tot;
  double *F;  

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
  } /* END for each node */ 
}


double DeltaOmegaF(int ck, double tau_eff_sym, double tau_eff_ant, double *u, double *F, struct System *sys)
/* Becomes SRT when tau_eff_sym = tau_eff_ant*/
{
  double C2 = sys->cs_2;
  double C4 = sys->cs_2*sys->cs_2;
  int *ev_ind = sys->ev_ind;
  double *ev = sys->ev;
  double *w = sys->w;
  double cu_tmp, cF_tmp, Fu_tmp;
    
  /* Adding a body force */
  Fu_tmp = 0.0;
  cF_tmp = 0.0;
  cu_tmp = 0.0;
  for (int nd = 0; nd < sys->n_D; nd++) {
    Fu_tmp += F[nd]*u[nd];
    cF_tmp += ev[ev_ind[ck]+nd]*F[nd];
    cu_tmp += ev[ev_ind[ck]+nd]*u[nd];
  }

  return w[ck]*( (1.0 - 0.5/tau_eff_ant)*cF_tmp/C2 + 
		 (1.0 - 0.5/tau_eff_sym)*(cF_tmp*cu_tmp - C2*Fu_tmp)/C4 );
}

double cos_phi(int ck, double * grad, struct System *sys)
{
  double ret = 0;
  double *ev;
  double grad_norm = 0;
  double ev_norm = 0;

  ev = &sys->ev[sys->ev_ind[ck]];

  for (int cd = 0; cd < sys->n_D; cd++) {
    ret += ev[cd]*grad[cd];
    grad_norm += grad[cd]*grad[cd];
    ev_norm += ev[cd]*ev[cd];
  }

  if ( (grad_norm == 0) || (ev_norm == 0) )
    return 0;

  ret /= sqrt(grad_norm)*sqrt(ev_norm);

  return ret;
}

double Omega_surface_tension(int ck, double * grad, double A, double * B, struct System *sys)
{
  double ret;
  double *ev;
  double grad_norm = 0;
  double ev_dot_grad = 0;

  ev = &sys->ev[sys->ev_ind[ck]];
  for (int cd=0; cd < sys->n_D; ++cd) {
    grad_norm += grad[cd]*grad[cd];
    ev_dot_grad += ev[cd]*grad[cd];
  }

  if (grad_norm == 0) 
    return 0;

  ret = sys->w[ck]*ev_dot_grad*ev_dot_grad/grad_norm - B[ck];
  ret *= 0.5*A*sqrt(grad_norm);

  return ret;
}



/* ---------------------------------------------------------------------------------  color_gradient */
void color_gradient(int cn, double *grad, double *rho_cc1, double *rho_cc2, struct System *sys)
/* Color gradient used in phys.rev.E 85, 046309 (2012) */
{
  struct Step *step = sys->step;
  int *ev_ind = sys->ev_ind;
  double *ev = sys->ev;
  double *w = sys->w;
  int cd, ck ;
  double rho_n;
  int cn_neig;
  double rho_sum, rho_diff;

  for (cd = 0; cd < sys->n_D; ++cd) {
    grad[cd] = 0.0;
    for (ck = 1; ck < sys->n_Q; ++ck) {
      cn_neig = cn - step->n_rel_neig[ck];
      rho_sum = (rho_cc1[cn_neig] + rho_cc2[cn_neig]);
      rho_diff = (rho_cc1[cn_neig] - rho_cc2[cn_neig]);
      if (rho_sum == 0.0)
        rho_n = 0.0;
      else
        rho_n =rho_diff/ rho_sum;
      grad[cd] += w[ck]*ev[ev_ind[ck]+cd]*rho_n;
    }
    grad[cd] *= 3;  
  }
}



#endif // _USE_COLOR_GRAD_


/* ---------------------------------------------------------------------------------  collision_f */
void collision_f(struct Node *node, struct Field *fluid, struct System *sys)
/* 
   collision_f :
   Accomplish a collision for the fluid part of the system 

   INPUT  : void
   OUTPUT : void
*/
{
  int cn, cc, cd, ck; /* Counter: node, phase, direction, dimension */
  int u_ind, f_ind; /* index in for the velocity and velocity distribution */
  real tau, tau_inv; /* relaxation time, invers of the relaxation time */ 
  real *f, *f_tmp, *rho, *u; /* Pointers to the distribution, temprorar distribution, density and velocity */

  struct Step *step = sys->step;
  //struct Node *node = sys->node;

  for( cc = 0; cc < fluid->n_c; ++cc ) { /* phases */
    f     = &fluid->f    [step->f_phase*cc];
    f_tmp = &fluid->f_tmp[step->f_phase*cc];
    rho   = &fluid->rho  [step->rho_phase*cc];
    u     = &fluid->u    [step->u_phase*cc];

    tau     = fluid->tau[cc];
    tau_inv = 1./tau;

    u_ind = 0;
    f_ind = 0;

    for( cn = 0; cn < sys->max_n_tot; ++cn) { /* nodes */
      if( node->mode[cn]==FLUID || node->mode[cn]==INLET || node->mode[cn]==OUTLET ) { /* part of the fluid */

        /* Calculate the density and the velocity */
        calc_rho_vel(&f_tmp[f_ind], &rho[cn], &u[u_ind], sys);

        /* Add gravity */
        for( cd = 0; cd < sys->n_D; ++cd )
          fluid->u_eq[cd] = u[u_ind+cd] + tau*fluid->gravity[cd]/rho[cn];

        /* Collision */
        for( ck = 0; ck < sys->n_Q; ++ck) { /* directions */
          f[f_ind + ck] = fluid->f_tmp[f_ind + ck] - (fluid->f_tmp[f_ind + ck] - calc_feq(ck, rho[cn], fluid->u_eq, sys))*tau_inv;
        }
      }
      u_ind += step->u_node;
      f_ind += step->f_node;
    }

  }
}



/* ---------------------------------------------------------------------------------  collision_f_newconc */
void collision_f_newconc(struct Node *node, struct Field *fluid, struct System *sys)
/*
   collision_f :
   Accomplish a collision for the fluid part of the system

   INPUT  : void
   OUTPUT : void
*/
{
  int cn, cc, cd, ck; /* Counter: node, phase, direction, dimension */
  int u_ind, f_ind; /* index in for the velocity and velocity distribution */
  real tau, tau_inv; /* relaxation time, invers of the relaxation time */
  real *f, *f_tmp, *rho, *u; /* Pointers to the distribution, temprorar distribution, density and velocity */
  double C2, C4;
  real delta_F;
  //double delta_Omega; /* Correction term for force and source terms */
  double q; /* NB! local source term. Should be a seperate scalar field */
  double *F; /* Force */
        
  struct Step *step = sys->step;
  //struct Node *node = sys->node;

  int *ev_ind = sys->ev_ind;
  double *ev = sys->ev;
  double *w = sys->w;

  // obsolete color-grad code?
  //  int nx, ny, nz; /* hold the y-positon */
  //  int nxyz[3];
  //  int max_nz = 0;
  //  int min_nz = 100;

  C2 = sys->cs_2;
  C4 = C2*C2;

  F = fluid->gravity;
  q = 0.0;

  for( cc = 0; cc < fluid->n_c; ++cc ) { /* phases */
    f     = &fluid->f    [step->f_phase*cc];
    f_tmp = &fluid->f_tmp[step->f_phase*cc];
    rho   = &fluid->rho  [step->rho_phase*cc];
    u     = &fluid->u    [step->u_phase*cc];
    
    tau     = fluid->tau[cc];
    tau_inv = 1./tau;
    
    u_ind = 0;
    f_ind = 0;

    //printf("Q_source: %g\n", fluid->q_source);
    
    for( cn = 0; cn < sys->max_n_tot; ++cn) { /* nodes */
      if( node->mode[cn]==FLUID || node->mode[cn]==INLET || node->mode[cn]==OUTLET ) { /* part of the fluid */

#ifdef _FLUID_SOURCE_
        q = 0;
        int x = get_global_xcoord(cn, sys);
        //int y = get_global_ycoord(cn, sys);
        if (x==3)
          q = fluid->q_source; //0.5e-2;
        if (x==sys->MAX_N[0]-1-3)
          q = -fluid->q_source;
#endif
#ifdef _FLUID_BDRY_PRESS_
        q = 0;
        int x = get_global_xcoord(cn, sys);
        if ( (x==node->in_x) || (x==node->out_x) ) {
          double f_sum = 0;
          for( ck = 0; ck < sys->n_Q; ++ck)
            f_sum += f_tmp[f_ind + ck];
          if (x==node->in_x)
            q = 2.0*(fluid->rho_inlet[cc] - f_sum);
          if (x==node->out_x)
            q = 2.0*(fluid->rho_outlet[cc] - f_sum);
        }
#endif

        /* Calculate the density and the velocity */
        calc_rho_vel_newconc(&f_tmp[f_ind], &rho[cn], &u[u_ind], &q, F, sys);
        /* Add gravity */
        for( cd = 0; cd < sys->n_D; ++cd )
          fluid->u_eq[cd] = u[u_ind+cd];


        /* Collision */
        for( ck = 0; ck < sys->n_Q; ++ck) { /* directions */

#if defined(_FLUID_SOURCE_) or defined(_FLUID_BDRY_PRESS_)

          double delta_Omega = 0.0;
          /* calculate the dot products */
          double e_dot_u = 0.0;
          double e_dot_F = 0.0;
          double u_dot_u = 0.0;
          double u_dot_F = 0.0;

          for (cd = 0; cd < sys->n_D; ++cd) {
            e_dot_u += ev[ev_ind[ck]+cd]*u[u_ind+cd];
            e_dot_F += ev[ev_ind[ck]+cd]*F[cd];
            u_dot_u += u[u_ind+cd]*u[u_ind+cd];
            u_dot_F += u[u_ind+cd]*F[cd];
          }
          /* -- force corrections */
          delta_Omega += e_dot_F/C2 + (e_dot_u*e_dot_F - C2*u_dot_F)/C4;
          /* -- mass source corrections */
          //delta_Omega += q - (e_dot_u*e_dot_u - C2*u_dot_u)/(2*C4);
          delta_Omega += q*(1.0 + e_dot_u/C2 + (e_dot_u*e_dot_u - C2*u_dot_u)/(2*C4));
          /* -- multiplicative factor */
          delta_Omega *= w[ck]*(1.0 - 0.5*tau_inv);
          f[f_ind + ck] = fluid->f_tmp[f_ind + ck] - (fluid->f_tmp[f_ind + ck] - calc_feq(ck, rho[cn], fluid->u_eq, sys))*tau_inv + delta_Omega;
#else
          delta_F = 0.0;
          double e_dot_u = 0.0;
          for (cd = 0; cd < sys->n_D; ++cd) {
            e_dot_u += ev[ev_ind[ck]+cd]*u[u_ind+cd];
          }
          e_dot_u *= 1.0/sys->cs_2;

          for (cd = 0; cd < sys->n_D; ++cd) {
            delta_F += ((ev[ev_ind[ck]+cd] - u[u_ind+cd]) + e_dot_u*ev[ev_ind[ck]+cd])*fluid->gravity[cd];
          }
          delta_F *= w[ck]*(1.0 - 0.5*tau_inv)/sys->cs_2;
          f[f_ind + ck] = fluid->f_tmp[f_ind + ck] - (fluid->f_tmp[f_ind + ck] - calc_feq(ck, rho[cn], fluid->u_eq, sys))*tau_inv + delta_F;
#endif
        }
      }
      u_ind += step->u_node;
      f_ind += step->f_node;
    }
    
  }
}

///* ---------------------------------------------------------------------------------  collision_shanchen_f */
//void collision_shanchen_f(struct Field *fluid, struct System *sys)
///*
//collision_shanchen_f :
//        Collision for the fluid part of the system including the phase forces
//
//  INPUT  : void
//  OUTPUT : void
//  */
//{
//        int cn, cc, ccl, cd, ck; /* Counter: node, phase, phase-local , direction, dimension */
//        int u_ind, f_ind; /* index in for the velocity and velocity distribution */
//        real tau, tau_inv; /* relaxation time, invers of the relaxation time */
//        real *f, *f_tmp, *rho, *u; /* Pointers to the distribution, temprorar distribution, density and velocity */
//        real u_norm; /* Normalization factor */
//
//        struct Step *step = sys->step;
//        struct Node *node = sys->node;
//
//        for( cc = 0; cc < fluid->n_c; ++cc )
//        {
//                f = &fluid->f[step->f_phase*cc];
//                f_tmp = &fluid->f_tmp[step->f_phase*cc];
//                rho = &fluid->rho[step->rho_phase*cc];
//                u   = &fluid->u[step->u_phase*cc];
//
//                tau     = fluid->tau[cc];
//                tau_inv = 1./tau;
//
//                u_ind = 0;
//                f_ind = 0;
//                for( cn = 0; cn < sys->max_n_tot; ++cn) /* nodes */
//                {
//                        /* Initiate the barisentric veloctiy */
//                        for( cd = 0; cd < sys->n_D; ++cd )
//                                fluid->u_tot[u_ind + cd] = 0.0;
//
//                        if( node->mode[cn]==FLUID ) /* part of the fluid */
//                        {
//                                /* Calculate the density and the velocity */
//                          calc_rho_vel_cn_f(cn, &u_norm, fluid, sys);
//
//
//                                for( cd = 0; cd < sys->n_D; ++cd )
//                                {
//                                        /* Update the baricentric velcotiy */
//                                        fluid->u_tot[u_ind + cd] = fluid->u_eq[cd]/u_norm;
//
//                                        /* Add gravity */
//                                        fluid->u_eq[cd] = fluid->u_eq[cd]/u_norm + tau*fluid->gravity[cd]/u_norm;
//
//                                        /* Multiphase */
//                                        /* -- rock fluid interaction */
//                                        fluid->u_eq[cd] -= fluid->WGB[cc]*tau*fluid->glforce[step->u_phase*cc + step->u_node*cn + cd];
//                                        /* -- fluid-fluid interaction */
//                                        for( ccl = 0; ccl < fluid->n_c; ++ccl )
//                                                fluid->u_eq[cd] -= fluid->WFB[fluid->n_c*cc + ccl]*tau*fluid->flforce[step->u_phase*ccl + step->u_node*cn + cd];
//                                }
//
//                                /* Collision */
//                                for( ck = 0; ck < sys->n_Q; ++ck) /* directions */
//                                  f[f_ind + ck] = fluid->f_tmp[f_ind + ck] - (fluid->f_tmp[f_ind + ck] - calc_feq(ck, rho[cn], fluid->u_eq, sys))*tau_inv;
//
//
//                        }
//
//                        u_ind += step->u_node;
//                        f_ind += step->f_node;
//                }
//        }
//}
//
        


/* ---------------------------------------------------------------------------------  collision_g */
void collision_g(struct Field *species, struct Field *fluid, struct System *sys)
/* 
collision_g :
Accomplish a collision for the advection diffusion part of the system 

  INPUT  : void
  OUTPUT : void
*/
{
  int cn, cc, ck; /* Counter: node, phase, direction, dimension */
  int u_ind, g_ind; /* index in for the velocity and velocity distribution */
  real tau_inv; /* invers of the relaxation time */ 
  real *g, *g_tmp, *rho, *u; /* Pointers to the distribution, temprorar distribution, density and velocity */

  struct Step *step = sys->step;
  struct Node *node = sys->node;
  
  u = fluid->u;
#ifdef BB_INLET
  int *n = sys->max_n;
  struct Mpi *mpi = sys->mpi;
  int *irank = sys->mpi->ind_rank;
  int *np = sys->mpi->np;
#endif

  for( cc = 0; cc < species->n_c; ++cc ) {
    g     = &species->f[step->f_phase*cc];
    g_tmp = &species->f_tmp[step->f_phase*cc];
    rho   = &species->rho[step->rho_phase*cc];
    
    tau_inv = 1./species->tau[cc];
    
    u_ind = 0;
    g_ind = 0;

    for( cn = 0; cn < sys->max_n_tot; ++cn) { /* nodes */
      if( node->mode[cn]==FLUID ) { /* part of the fluid */

        /* Calculate the density and the velocity */ 
        calc_rho(&g_tmp[g_ind], &rho[cn], sys);

#ifdef BB_INLET
        // select rho_inlet if this is an inlet node
        if ( (cn%n[0]+mpi->lower_bound[0]-1)<4 && sys->t>10) {
          rho[cn] = species->rho_inlet[cc];
        }
#endif // BB_INLET

        /* Collision */
        for( ck = 0; ck < sys->n_Q; ++ck) { /* directions */
          g[g_ind + ck] = g_tmp[g_ind + ck] - (g_tmp[g_ind + ck] - calc_geq(ck, rho[cn], &u[u_ind], sys))*tau_inv;
        }
      }
      
      u_ind += step->u_node;
      g_ind += step->f_node;
    }
    
  }
}


/* ---------------------------------------------------------------------------------  collision_g */
void collision_g_newconc(struct Field *species, struct Field *fluid, struct System *sys)
/* 
collision_g :
Accomplish a collision for the advection diffusion part of the system 

  INPUT  : void
  OUTPUT : void
*/
{
  int cn, cc, ck; /* Counter: node, phase, direction, dimension */
  int u_ind, g_ind; /* index in for the velocity and velocity distribution */
  real tau_inv; /* invers of the relaxation time */ 
  real *g, *g_tmp, *rho, *u, *psi; /* Pointers to the distribution, temprorar distribution, density and velocity */

  double *w = sys->w;
  int *ev_ind = sys->ev_ind;
  double *ev = sys->ev;
  struct Step *step = sys->step;
  struct Node *node = sys->node;

  real geq;
  real delta_F = 0.0;
  double A, Apsi;
  double delta_src = 0.0;
  /* double R_Cl = 0.0;
  double R_NH3 = 2*R_Cl;
  double R_HCO3 = R_Cl;
  int pos_cl[2] = {10, 20}; */
  /* int pos_urea[2] = {46, 46}; AMIN
  int pos_Lac[2] = { 10, 20 }; */

  /* int pos_Lac_x [] = { 10, 55, 10, 40, 70  };
  int pos_Lac_y []= { 20, 20, 40, 40, 70 };
  double R_Lac = 1.0e-2; */
  double* R_cc;

#ifdef _SPECIES_SOURCE_
  double q = fluid->q_source;
  //printf("fluid->q_source = %.3e\n", q);
#endif
#ifdef _SPECIES_SOURCE_BULK_
  double R = 0.0, Q_Vp, rho_c, rho_mean, phi;
  get_mean_g(sys, species);
  phi = get_porosity(sys);
  Q_Vp = fluid->u_Darcy*species->dt/(phi*sys->L_lb*sys->dx);
#endif

  u = fluid->u;
  
  /*  if (!sys->files["species_rho_in"].empty()) {
    set_rho_inlet_from_file(species, sys);
    } */

  for( cc = 0; cc < species->n_c; ++cc ) {
    g     = &species->f[step->f_phase*cc];
    g_tmp = &species->f_tmp[step->f_phase*cc];
    rho   = &species->rho[step->rho_phase*cc];
    psi   = &species->psi[step->rho_phase*cc];
	R_cc  = &species->R_local[step->rho_phase*cc]; /* AMIN EJE */
    //printf("<%d> %d: inlet: %g\n", sys->mpi->my_rank, cc, species->rho_inlet[cc]);
	//std::cout << "\tcc:"<<cc<<"\t"<< species->R_local[step->rho_phase*cc] <<std::endl;/* AMIN */
    
    tau_inv = 1./species->tau[cc];
#ifdef _SPECIES_SOURCE_BULK_
    rho_c = species->rho_inlet[cc];
    rho_mean = (species->mean_g[cc] + 0.5*Q_Vp*rho_c)/(1.0 + 0.5*Q_Vp);
    R =  Q_Vp*(rho_c - rho_mean);
    //printf("%d: R = %.5e\n", cc, R);
#endif
    
    u_ind = 0;
    g_ind = 0;

    A = (1.0 - 0.5*tau_inv); // /sys->cs_2;
    /* for each !FLUID! node DO*/
    for( cn = 0; cn < sys->max_n_tot; ++cn) { /* nodes */
      if( node->mode[cn]==FLUID) { /* part of the fluid */
      //if( node->mode[cn]==FLUID || node->mode[cn]==INLET || node->mode[cn]==OUTLET ) { // NBNBNBNBNB
	
        /*    calculate rho */
        rho[cn] = 0;
        for( ck = 0; ck < sys->n_Q; ++ck ) {
          rho[cn] += g_tmp[g_ind+ck];
        }

		/* AMIN EJE: add HCl and UREA source */
		delta_src = R_cc[cn];
		//std::cout << "\tcn:" << cn << "\t" << delta_src << std::endl;/* AMIN */
		rho[cn] += 0.5 * delta_src;


#ifdef _SPECIES_SOURCE_BULK_
        // R = Q/V_p*(rho_c-rho[cn]) --> rho[cn]:= sum g_tmp + 0.5*R
        // rho = (sum g_tmp + 0.5*rho_c*Q/Vp)/(1.0 + 0.5*Q/Vp)
        //rho[cn] += rho_c*0.5*Q_Vp;
        //rho[cn] /= (1.0 + 0.5*Q_Vp);
        //R = Q_Vp*(rho_c - rho[cn]);
        rho[cn] += 0.5*R;
#endif
#ifdef _SPECIES_SOURCE_
	int x = get_global_xcoord(cn, sys);
	double R = 0;
	if (x==node->in_x) {
	  if (sys->t > 1e4)
	    R = q*species->rho_inlet[cc];
	  else
	    R = q*species->rho_init[cc];
	}
	if (x==node->out_x) {
	  R = -q*rho[cn]/(fluid->rho[cn]+0.5*q);
	}
	rho[cn] += 0.5*R;
#endif




	/* AMIN EJE if (cc == 1)  { 
	  int x = get_global_xcoord(cn, sys);
	  int y = get_global_ycoord(cn, sys);
	  if (x == pos_cl[0] && y == pos_cl[1]) {

		  rho[cn] += 0.5*R_Cl;
	    delta_src = R_Cl;
	  } 
	} else if (cc == 2) { // Lac
	  int x = get_global_xcoord(cn, sys);
	  int y = get_global_ycoord(cn, sys);
	  
	  for (int i_1 = 0; i_1 < sizeof(pos_Lac_x) / sizeof(pos_Lac_x[0]); i_1++) {

		  if (x == pos_Lac_x[i_1] && y == pos_Lac_y[i_1]) {
			  rho[cn] += 0.5*R_Lac;
			  delta_src = R_Lac;
		  }
	  }
	} else if (cc == 3) { // NH3 
	  int x = get_global_xcoord(cn, sys);
	  int y = get_global_ycoord(cn, sys);

	  if (x == pos_urea[0] && y == pos_urea[1]) {
	    rho[cn] += 0.5*R_NH3;
	    delta_src = R_NH3;
	  } 
	} */

#ifdef _FIND_PERC_NODES_
	set_rho_in_percolation_nodes(cn, rho, species->rho_inlet[cc], sys);
#endif
	
        /*    calculate psi */
#ifdef FLUID_OFF
        psi[cn] = rho[cn];
#else
        psi[cn] = rho[cn]/fluid->rho[cn];
#endif

        Apsi = A*psi[cn];
        for( ck = 0; ck < sys->n_Q; ++ck) {

          /*    calculate the force corrections */
#ifndef FLUID_OFF // true if fluid is on
          delta_F = 0.0;
          for (int cd = 0; cd < sys->n_D; ++cd) {
            delta_F += ev[ev_ind[ck]+cd]*fluid->gravity[cd];
          }
          //delta_F *= w[ck]*(1.0 - 0.5*tau_inv)*psi[cn]/sys->cs_2;
          delta_F *= w[ck]*Apsi/sys->cs_2;
#endif
#ifdef _SPECIES_SOURCE_BULK_
          delta_src = w[ck]*A*R;  // only valid for u=0
#endif
#ifdef _SPECIES_SOURCE_
          delta_src = 0.0;
          for (int cd = 0; cd < sys->n_D; ++cd) {
            delta_src += ev[ev_ind[ck]+cd]*u[u_ind + cd];
          }
          delta_src = (1.0 + delta_src/sys->cs_2)*w[ck]*A*R;
#endif
          /*    calculate the equilibrium distribution (see fluid field) */
#ifdef FLUID_OFF
          geq = calc_geq(ck, rho[cn], &u[u_ind], sys);
#else // fluid is on
          geq = calc_feq(ck, rho[cn], &u[u_ind], sys);
#endif
	  /*    update g*/
	  /* AMIN EJE: set the source */
	  delta_src *= w[ck];
	  delta_F = 0.0;
          g[g_ind + ck] = g_tmp[g_ind + ck] - (g_tmp[g_ind + ck] - geq)*tau_inv + delta_F + delta_src;
        }
      }
      
      u_ind += step->u_node;
      g_ind += step->f_node;
    }
  }  
}

//-----------------------------------------------------------
//
//
//-----------------------------------------------------------
void set_rho_in_percolation_nodes(int cn, double *rho, double rho_set, struct System *sys)
{
  if (get_global_xcoord(cn, sys)==sys->node->perc_inlet) {
    rho[cn] = rho_set;
  }
  if (sys->t%50==0) {
    for (std::vector<int>::iterator it=sys->node->perc_cluster.begin(); it != sys->node->perc_cluster.end(); ++it) {
      //if (cn==sys->node->perc_cluster[i]) {
      if (cn == *it) {
        rho[cn] = rho_set;
        break;
      }
    }
  }
}

//-----------------------------------------------------------
//
//
//-----------------------------------------------------------
void get_mean_g(struct System *sys, struct Field *species)
{
  int nn, cc, ck;
  double fluid_nodes = 0.0;
  int ind;
  struct Node *node = sys->node;
  struct Step *step = sys->step;
  struct Mpi *mpi = sys->mpi;
  
  // reset
  for (cc = 0; cc<species->n_c; ++cc) {
    species->mean_g[cc] = 0.0;
  }

  for( nn = 0; nn < sys->max_n_tot; ++nn ) {
    if (node->mode[nn]!=FLUID || node->is_ghost[nn])
      continue;
    fluid_nodes += 1.0;
    for (cc = 0; cc<species->n_c; ++cc) {
      ind = step->f_phase*cc + step->f_node*nn;
      for( ck = 0; ck < sys->n_Q; ++ck ) {
        species->mean_g[cc] += species->f_tmp[ind + ck];
      }
    }
  }
  
  mpi->sbuf_eff[0] = fluid_nodes;
  for (cc = 0; cc<species->n_c; ++cc) {
    mpi->sbuf_eff[1+cc] = species->mean_g[cc];
  }
  MPI_Reduce(mpi->sbuf_eff, mpi->rbuf_eff, 1+species->n_c, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  if (mpi->my_rank == 0) {
    for (cc = 0; cc<species->n_c; ++cc) {
      species->mean_g[cc] = mpi->rbuf_eff[1+cc]/mpi->rbuf_eff[0];
    }
  }
  MPI_Bcast(species->mean_g, species->n_c, MPI_DOUBLE, 0, MPI_COMM_WORLD);
}

//-----------------------------------------------------------
//
//
//-----------------------------------------------------------
double get_porosity(struct System *sys)
{
  int nfluid=0, nfluid_mpi, cn;
  double phi;
  struct Node *node = sys->node;
  
  for( cn = 0; cn < sys->max_n_tot; ++cn ) { /* node */
    if (node->mode[cn] == FLUID && !node->is_ghost[cn] ) {
      nfluid++;
    }
  }
  MPI_Reduce(&nfluid, &nfluid_mpi, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  if (sys->mpi->my_rank == 0) {
    phi = nfluid_mpi/(double)(sys->V_lb);
  }
  MPI_Bcast(&phi, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  return(phi);
}


//-----------------------------------------------------------
//
//
//-----------------------------------------------------------
/* void set_rho_inlet_from_file(Field *species, System *sys)
{
  static int lines_read = 0;

  if (species->rho_inlet_file.size()==0) {
    // read values from rho_effluent.dat into array
    std::ifstream infile(sys->files["species_rho_in"]);
    if (infile.fail()) {
      std::cerr << "ERROR! Unable to open " << sys->files["species_rho_in"] << std::endl;
      MPI_Finalize();
      exit(1);
    }
    std::string line;
    // read header names
    std::string name;
    std::vector<std::string> variable;
    getline(infile, line);
    std::istringstream ss(line);
    ss >> name; // get rid of # at the beginning
    while (ss >> name) {
      variable.push_back(name);
    }

    int num_pos = 1+species->n_c;
    std::vector<int> col_pos(num_pos);

    // find position of 'realtime' column
    col_pos[0] = std::find(variable.begin(), variable.end(), "realtime") - variable.begin();
    // find species column positions
    for (int cc = 0; cc < species->n_c; ++cc) {
      std::string var_name = "mean_" + std::string(species->name[cc]);
      col_pos[cc+1] = std::find(variable.begin(), variable.end(), var_name) - variable.begin();
      //std::cout << var_name << ": ind[" << cc+1 << "] = " << col_pos[cc+1] << std::endl;
    }

    // populate 2D array with file values
    int nline = 0;
    while (getline(infile, line)) {
      ++nline;
      if (nline <= lines_read)
        continue;
      std::vector<double> col;
      double value;
      std::istringstream ss2(line);
      while (ss2 >> value) {
        col.push_back(value);
      }
      std::vector<double> row(num_pos);
      int i = 0;
      for (auto &val : row) {
        val = col[col_pos[i++]];
      }
      species->rho_inlet_file.push_back(row);
    }

    if (species->rho_inlet_file.size()<1) {
      if (sys->mpi->my_rank==0) {
        std::cout << std::endl << std::endl << "*** Simulation stopped due to missing species concentrations in " << sys->files["species_rho_in"] << " ***" << std::endl << std::endl;
      }
      MPI_Finalize();
      exit(0);
    }

    if (sys->mpi->my_rank==0) {
      std::cout << std::endl << "Inlet species concentrations are dynamically set by values in file " << sys->files["species_rho_in"]
                             << " (" << lines_read+1 << " - " << lines_read+species->rho_inlet_file.size() << " lines) "<< std::endl;
    }

    lines_read += species->rho_inlet_file.size();
  }


  for (int cc = 0; cc < species->n_c; ++cc) {
    species->rho_inlet[cc] = species->rho_inlet_file[0][cc+1];
  }

  if ( (sys->time+sys->dt > species->rho_inlet_file[0][0]) && (species->rho_inlet_file.size()>0))
    species->rho_inlet_file.erase(species->rho_inlet_file.begin());

  //  if (sys->mpi->my_rank==0)
  //    std::cout << std::endl << "time: " << species->rho_inlet_file[0][0] << ", Ca: " << species->rho_inlet_file[0][1] << ", size: " << species->rho_inlet_file.size() << std::endl;
  } */

