#include "global.h"
#include "chem_global.h"
#include "biozementfunctions.h"

void wall_bc_bounce_back(struct Node *node, struct Field *fluid, struct System *sys)
/*
  wall_bc_bounce_back :
  Mid grid bounce back algorithm


  OUTPUT : void
*/
{
  int nl, cc; /* counters: links, phases  */
  int dir, dir_bb; /* Direction,  bounce back direction*/
  real *fg; /* pointer to the first entry of a phase in fluid->f */

  struct Step *step = sys->step;
  struct Links *links = node->links;
  struct Links::Link link;

  /* loop phases */
  for( cc = 0; cc < fluid->n_c; ++cc ) {
    fg     = &fluid->f[step->f_phase*cc];
    /* loop links */
    for (nl = 0; nl < links->nlinks; ++nl) {
      link = links->list[nl];
      dir    = link.dir;  /* direction from solid to fluid */
      dir_bb = sys->k_bb[dir];  /* direction from fluid to solid */
      /* Copy fluid value to solid and change direction */
      fg[step->f_node*link.wall + dir] = fg[step->f_node*link.fluid + dir_bb];
    }
  }
  
  
  /* Dette skal bli gjort bedre (kun for Ã¥ teste) */
  periodic_bc(node, fluid, sys);
}


/* ---------------------------------------------------------------------------------*/
void wall_bc_mid_grid_gca_ej_g_v3(struct InitChem *ICS, struct Field *species, struct Field *fluid,
    struct System *sys, struct Minerals *minerals, struct Boundary *bndry, struct splayTree **st_lst)
{
#ifdef ALL_INERT
  wall_bc_bounce_back(sys->node, species, sys);
  return;
#endif
  
  /* NB Her har jeg oversatt cw med nl */
  int nl, cc; /* counters: links, phases, groups  */
  int dir, dir_bb; /* Direction,  bounce back direction*/
  real *fg; /* pointer to the first distribution entry of the solid node */

  struct Step *step = sys->step;
  struct Node *node = sys->node;
  struct Links *links = node->links;
  struct Links::Link *link;
  struct BasVec *Vchem;
  int cw, min, i, j;
  double *min_conc;
  real delta[20];
  int key, pos;
  double fluid_f;
  
  // prepare for output of pH
  struct treeNode *tn;
  struct ChemTable *SM_basis;
  SM_basis = ICS->SM_basis;

#ifdef _MAGN_ON_MAGN_
  min_conc = (double *) malloc((minerals->n_tot+1)*sizeof(double));
#else
  min_conc = (double *) malloc(minerals->n_tot*sizeof(double));
#endif
  
  for (nl = 0; nl < links->nlinks; ++nl) { /* loop links */
    link = &links->list[nl];
    dir = link->dir; /* direction from solid to fluid */
    dir_bb = sys->k_bb[dir]; /* direction from fluid to solid */
    
    if (link->inert) {
      // only do bounce-back, and move to next link
      for (cc = 0; cc < species->n_c; ++cc) {/* loop phases */
        fg = &species->f[step->f_phase*cc];
        fg[step->f_node*link->wall + dir] = fg[step->f_node*link->fluid + dir_bb];
      }
      continue;
    }
    
#ifdef FAKE_CHEM

    applyFakeChem(link, species, step, sys);

/*    double *g_in, *g_out;
    double c_eq = 0.0, k=0.1, c2_inv = 3.0;
    double kc2 = k*c2_inv;
    g_in = &(species->f[step->f_node*link->fluid + dir_bb]);
    g_out = &(species->f[step->f_node*link->wall + dir]);
    for (cc = 0; cc < species->n_c; ++cc) {
        g_out[step->f_phase*cc] = ( g_in[step->f_phase*cc]*(1.0-kc2) + kc2*2.0*sys->w[dir]*c_eq )/(1.0+kc2);
    } */

#endif

#ifndef FAKE_CHEM
    /* Copy distribution */
    fg = &(species->f[step->f_node*link->fluid + dir_bb]);

    for (cc = 0; cc < species->n_c; ++cc) /* loop phases */
        bndry->flux_l_wall[sys->step->wall_comp*cc] =
          bndry->rho_l_wall_tmp[cc] =
              fg[step->f_phase*cc];
    
    cw = link->wall;

    for (min = 0; min < minerals->n_tot; ++min) {
      min_conc[min] = minerals->list[min].sfrac[cw];
    }

#ifdef _MAGN_ON_MAGN_
    min_conc[minerals->n_tot] = (double) minerals->nucleation_site[cw];
#endif
    
    /* returns the key to the correct Vchem at that node */
#ifdef NEWCONC

#ifdef FLUID_OFF
    fluid_f = sys->w[dir];
#else // fluid on
    fluid_f = fluid->f[step->f_node*link->fluid + dir_bb];
#endif

    key = calc_eq_rho_rate_newconc(sys, bndry, ICS, minerals, fluid_f, species->n_c, 1.0, dir, min_conc, cw, &tn, st_lst);
#else // newconc off
    key = calc_eq_rho_rate_jlv(sys, bndry, ICS, minerals, species->n_c, 1.0, dir, min_conc, cw, &tn, st_lst);
#endif // NEWCONC


    Vchem = st_lst[key]->Vchem_key;

    // Store pH values for output
    // update pH of link
    if (ICS->INTERPOLATE) {
      // with splaytree
      link->pH = tn->val[st_lst[key]->num_val];
    } else {
      // only if SplayTree is OFF
      link->pH = Vchem->log_a[Vchem->pos_pH];
    }

    /* Assume that the new distributions are given in the bndry->rho_eq_wall_tmp variable */
    /* Complete the bounce back */
    for( cc = 0; cc < species->n_c; ++cc )  /* for each phase */
      bndry->flux_l_wall[step->wall_comp*cc] -=
          species->f[step->f_phase*cc + step->f_node*link->wall + dir] =
              bndry->rho_eq_wall_tmp[cc];

    // calculate basis species difference
    for (cc = 0; cc < species->n_c; ++cc) {
      delta[ICS->pos[cc]] = bndry->flux_l_wall[step->wall_comp*cc];
    }

    // calculate mineral change from species change
    //for(min=0;min<Vchem->size_rock;++min) {
    for (min = 0; min < minerals->n_tot; ++min) {
      // set to zero because size_rock might be < minerals->n_tot
      link->delta[min] = 0.;
    }

    // obtain mineral change from stoichimetric matrix 
    // and basis species differences
    for(j=0;j<Vchem->size_sup_min;++j) {
      pos = Vchem->pos_rel_sup_min[j];
      for(i=0;i<Vchem->size_sup_min;++i) {
        link->delta[pos] += Vchem->sm_sup_buf[j][i]*delta[Vchem->pos_sup_bas[i]];
      }
    } 

#ifdef _MAGN_ON_MAGN_
    if (min_conc[1]<sys->opt->magn_limit && minerals->nucleation_site[cw]==0)
      link->delta[1] = 0.;
#endif

    // delta_m_LB * dC * dx^3 = [mol]
    // convert delta from mol/L to volume fraction
    for(i=0;i<ICS->size_sup_min;++i) {
      pos = ICS->pos_sup_min[i];
      link->delta[i] *= 1.0e-3*ICS->SM_mineral->mol_volume[pos];
    }
    
    if (sys->convergence) {
      for (min = 0; min < minerals->n_tot; ++min) {
	link->delta_sum[min] += link->delta[min];
      }
    }

  #endif
  } /* end-block for number if links (nl) */

  free(min_conc);
  
  update_solid_v5(sys, minerals, links);
#ifndef NOT_MOVE // if not NOT_MOVE, i.e. if moving boundaries
  update_nodes(sys, bndry, minerals);
  update_linklist(node, sys);
  reset_rho_f_in_new_solid_nodes(species, sys);
#endif

#ifdef CONVERGE
  if (sys->t > sys->steps_to_wait_before_converge)
      check_convergence_v2(sys, minerals, species, bndry, links);
#endif // end CONVERGENCE

  /* Dette skal bli gjort bedre (kun for test) */
  periodic_bc(node, species, sys);

}



/* ---------------------------------------------------------------------------------*/
int calc_eq_rho_rate_jlv(struct System *sys, struct Boundary *bndry,
    struct InitChem *ICS, struct Minerals *minerals, int n_c, real norm_fac,
			  int ck, real *min_conc, int nn, struct treeNode **tn, struct splayTree **st_lst)
/*
calc_eq_rho_rate_gca :
             Calculates the bounceback distribution for the ck direction.

  INPUT  : 1) g_out_wall  : the bounce back distribution (to be calculated) (old  rho_eq_wall)
	   2) g_in_wall   : known distribution directed into the wall  (old rho_l_wall)
	   3) norm_fac    : scaling of the rate function
	   4) cw          : wall number
	   5) ck          : distribution direction (g_in_wall direction)

  OUTPUT : void
 */
{
  int cc, key;
  real *g_out_wall = bndry->rho_eq_wall_tmp;
  real *g_in_wall = bndry->rho_l_wall_tmp;
  
  /* Rescale the distribution for the chem solver */
  for( cc = 0; cc < n_c; ++cc ) {
    g_in_wall[cc] /= sys->w[ck];
  }
  /* g_out_wall holds the consentrations given by the chem_solver*/

  key = chem_calc_eq_rho_lb_moving(g_out_wall, g_in_wall, ICS, min_conc, tn, st_lst, sys);
  
  /* Calculate the bounce back distributions */
  for( cc = 0; cc < n_c; ++cc )
    g_out_wall[cc] = (2.*g_out_wall[cc] - g_in_wall[cc])*sys->w[ck];

  return key;
}


/* ---------------------------------------------------------------------------------*/
int calc_eq_rho_rate_newconc(struct System *sys, struct Boundary *bndry,
    struct InitChem *ICS, struct Minerals *minerals,  real fluid_f,
    int n_c, real norm_fac, int ck, real *min_conc, int nn, struct treeNode **tn, struct splayTree **st_lst)
/*
calc_eq_rho_rate_gca :
             Calculates the bounceback distribution for the ck direction.

  INPUT  : 1) g_out_wall  : the bounce back distribution (to be calculated) (old  rho_eq_wall)
	   2) g_in_wall   : known distribution directed into the wall  (old rho_l_wall)
	   3) norm_fac    : scaling of the rate function
	   4) cw          : wall number
	   5) ck          : distribution direction (g_in_wall direction)

  OUTPUT : void
 */
{
  int cc, key;
  real *g_out_wall = bndry->rho_eq_wall_tmp;
  real *g_in_wall = bndry->rho_l_wall_tmp;

  /* Rescale the distribution for the chem solver */
  for( cc = 0; cc < n_c; ++cc ) {
    g_in_wall[cc] /= fluid_f;
#ifdef SET_NEG_RHO_ZERO
    if (g_in_wall[cc]<0) {
      g_in_wall[cc] = 1e-8;
      sys->num_negative++;
    }
#endif

  }
  /* g_out_wall holds the consentrations given by the chem_solver*/
  key = chem_calc_eq_rho_lb_moving(g_out_wall, g_in_wall, ICS, min_conc, tn, st_lst, sys);

  /* Calculate the bounce back distributions */
  for( cc = 0; cc < n_c; ++cc )
    g_out_wall[cc] = (2.*g_out_wall[cc] - g_in_wall[cc])*fluid_f;
  
  return key;
} 


