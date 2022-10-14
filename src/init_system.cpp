#include "global.h"
#include "output.h"

void init_structs(Input &o, struct Options *opt, struct Step **step_ptr,
  struct Minerals **minerals_ptr, struct Links **links_ptr, struct Node **node_ptr,
  struct System *sys, struct Field **fluid_ptr, struct Field **species_ptr, struct Field **dif_ptr, struct Boundary **bndry_ptr)

{

  set_options(opt, o);

  // allocations
  *step_ptr     = new Step();
  *minerals_ptr = new Minerals();
  *links_ptr    = new Links();
  *node_ptr     = new Node();
  *fluid_ptr    = new Field();
  *species_ptr  = new Field();
  *bndry_ptr    = new Boundary();
  // New diffusive field
  *dif_ptr      = new Field();


  // system
  sys->mpi              = new Mpi();
  sys->max_n            = (int *) calloc(3, sizeof(int));
  sys->MAX_N            = (int *) calloc(3, sizeof(int));
  sys->output           = NULL; //*output;
  sys->opt              = opt;
  sys->step             = *step_ptr;
  sys->node             = *node_ptr;
  sys->n_reinit         = 0;
  sys->last_convergence = 0;
  sys->convergence = 0;
  sys->t = 0;
  sys->vel_reinit_steps = 0;
  sys->species_not_converged = sys->opt->find_steady_state;
  sys->velreinit = sys->rho_reinit = 0;
  sys->solid2fluid = sys->fluid2solid = 0;
  sys->num_new_nodes = sys->num_new_solid_nodes = 0;
  sys->t_extra = 0.;
  sys->time = 0.;
  sys->dt = 0;
  sys->nwrite = 0;
  sys->num_negative = 0;
  sys->steps_to_wait_before_converge = 0;
  sys->fluid                 = 1;
  sys->chem                  = o["switch"]["chem"]; 
  sys->mpi->pbc_sides        = o["switch"]["periodic_mpi"];
  sys->splaytree_res         = o["splaytree"]["resolution"];  
  sys->conv_jump             = o["convergence"]["jump"];
  sys->conv_crit_vel         = o["convergence"]["crit_vel"];
  sys->conv_crit_vel_target  = o["convergence"]["crit_vel_target"];
  sys->conv_crit_chem        = o["convergence"]["crit_chem"];
  sys->conv_crit_ss          = o["convergence"]["crit_ss"];
  sys->conv_crit_ss_min      = o["convergence"]["crit_ss_min"];
  sys->conv_length           = o["convergence"]["length"];
  sys->steps_until_check     = o["convergence"]["steps_until_check"];
  sys->D                     = 2.3e-9; // self-diffusion of water is the default value
  sys->D_b                   = 0.0; // Biomass doesn't grow /* AMIN */
  //sys->dt_ratio              = 1.0;

  // minerals
  struct Minerals *minerals = *minerals_ptr;
  minerals->n_tot = 0;
  minerals->SA_tot = 0.0;
  minerals->init_molvol = 0.0;
  minerals->no_comb = 0;

  // node
  struct Node *node = *node_ptr;
  node->links = *links_ptr;
  node->nr_per = node->nr_in = node->nr_out = node->nr_ghost = node->nr_wall = 0;
  node->in_x = o["geometry"]["inlet_x"];
  node->out_x = o["geometry"]["outlet_x"];
  //
  node->is_fluid_in_out[FLUID] = node->is_fluid_in_out[INLET] = node->is_fluid_in_out[OUTLET] = 1;
  node->is_fluid[FLUID] = 1;

  // fluid
  struct Field *fluid = *fluid_ptr;
  fluid->n_c = 0;
  fluid->step = *step_ptr;
  fluid->mean_velx = 0.;
  fluid->perm = 0.;
  fluid->flux = 0.;
  fluid->area = 0.;
  fluid->q_source = 0.; //o["LB"]["q_source"];
  fluid->u_Darcy = 0.;
  fluid->grad_P = 0.;

  // species
  struct Field *species = *species_ptr;
  species->n_c = 0;
  species->step = *step_ptr;
  species->mean_velx = 0.;
  species->perm = 0.;
  species->flux = 0.;
  species->area = 0.;
  species->q_source = 0.;
  species->u_Darcy = 0.;
  species->grad_P = 0.;

  // diffusive field
  struct Field *dif = *dif_ptr;
  dif->n_c = 0;
  dif->step = *step_ptr;
  dif->mean_velx = 0.;
  dif->perm = 0.;
  dif->flux = 0.;
  dif->area = 0.;
  dif->q_source = 0.;
  dif->u_Darcy = 0.;
  dif->grad_P = 0.;
  
  // boundary
  struct Boundary *bndry = *bndry_ptr;
  bndry->n_liquid = 0;
  bndry->n_solid = 0; //inp["num_solid"];

  // ----------------------------
  // checks and warnings
  // ----------------------------
  if ((sys->chem) && (opt->closed_cell==0) && ((opt->inertx_inlet<opt->in_nodes) || (opt->inertx_outlet<opt->out_nodes))) {
    printf("OPTION WARNING  Non-inert nodes in inlet or outlet chamber\n");
  }

  if (opt->find_steady_state && opt->skip_init_velocity_run) {
    printf("OPTION WARNING  find_steady_state and skip_init_velocity_run are both set in options.c: find_steady_state is now turned OFF!\n");
    opt->find_steady_state = 0;
  }
  if (opt->eff_offset < 1) {
    printf("OPTION WARNING  eff_offset is less than 1 (%d), and is now set to 1\n", opt->eff_offset);
    opt->eff_offset = 1;
  }
  if (opt->restart_run) {
    opt->find_steady_state = 2;
    opt->skip_init_velocity_run = 1;
    opt->read_vel_from_file = 0;
    opt->read_species_from_file = 0;
  }
#ifdef _FLUID_BDRY_PRESS_
  if (opt->copy_in_out<3) {
    printf("OPTION WARNING  Pressure boundary conditions (_FLUID_BDRY_PRESS_) need at least 3 identical\n");
    printf("OPTION WARNING  slices at inlet and outlet. Option copy_in_out changed from %d to 3.\n", opt->copy_in_out);
    printf("OPTION WARNING  Update copy_in_out to at least 3 in options.dat to silence this warning\n");
    opt->copy_in_out = 3;
  }
#endif

  #ifdef _FIND_PERC_NODES_
  if (opt->find_isolated_nodes == 0) {
    printf("OPTION WARNING  find_isolated_nodes is off, but macro _FIND_PERC_NODES_ is on\n");
    printf("OPTION WARNING  Either both is off, or both is on.\n");
    printf("OPTION WARNING  Aborting....\n");
    exit(0);
  }
#endif
  if (opt->find_isolated_nodes > 0) {
    sys->chem = 1;            // to allocate species->rho and psi
    opt->closed_cell = 1;
  }
}


/* ---------------------------------------------------------------------------------  init_system */
void init_system(struct System *sys)
/*
init_system :
    sets the non distribution specifics variables and arrays to their default values

  INPUT  : void

  OUTPUT : void
 */
{
  int ck, cd;   /* counters: directions, dimensions */
  int mult_d;   /* Used in finding propagating neighbors */
  //char fn[200]; /* File name used to read the basis vectors */

  struct Step *step = sys->step;
  int *n = sys->max_n;
  int dir;

  /* *** NUMERICS *** */
  /* Set the step values */
  /* -- velocity distribution 'fluid->f' */
  step->f_phase = sys->max_n_tot*sys->n_Q;
  step->f_node = sys->n_Q;
  /* -- for nodes */
  step->n_phase = sys->max_n_tot;
  /* -- basis velocity/position vector */
  step->stp_e_k = sys->n_D;
  /* -- density */
  step->rho_phase = sys->max_n_tot;
  step->rho_node = 1;
  /* -- velocities */
  step->u_phase = sys->max_n_tot*sys->n_D;
  step->u_node = sys->n_D;

  /* Set propagation neighbor */
  for (ck = 0; ck < sys->n_Q; ++ck)
  {
    mult_d = 1;
    step->n_rel_neig[ck] = 0;
    for (cd = 0; cd < sys->n_D; ++cd)
    {
      step->n_rel_neig[ck] -= mult_d*sys->ei[step->stp_e_k*ck + cd];
      mult_d *= n[cd];
    }

    step->f_rel_neig[ck] = step->f_node*step->n_rel_neig[ck];
  }

  /* to speed up collision-routine */
  for (ck = 0; ck < sys->n_Q; ++ck) { /* directions */
    sys->ev_ind[ck] = sys->step->stp_e_k*ck;
  }

  // find negative and positive x-direction for outlet_bc
  for(ck = 0; ck < sys->n_Q; ++ck) { /* directions */
    dir = sys->n_D*ck;
    if (sys->n_D > 2) {
      // 3D
      if ((sys->ei[dir] == -1) && (sys->ei[dir + 1] == 0) && (sys->ei[dir + 2] == 0)) {
        step->neg_x_dir = ck;
        //break;
      }
      if ( (sys->ei[dir]==1) && (sys->ei[dir+1]==0) && (sys->ei[dir+2]==0) ) {
        step->pos_x_dir = ck;
      }
    }
    else {
      // 2D
      if ((sys->ei[dir] == -1) && (sys->ei[dir + 1] == 0)) {
        step->neg_x_dir = ck;
        //break;
      }
      if ( (sys->ei[dir]==1) && (sys->ei[dir+1]==0) ) {
        step->pos_x_dir = ck;
      }
    }
  }
  //printf("\n\n---------- neg_x_dir = %d\n\n", step->neg_x_dir);
  //printf("\n\n---------- pos_x_dir = %d\n\n", step->pos_x_dir);

}


/* ---------------------------------------------------------------------------------  init_system_fg */
void init_system_fg(struct Field *field, struct System *sys)
/*
init_system_fg:
    set distribution relevant values to a default value

  INTPUT: 1) fg_c     : pointer to distribution
          2) fg_tmp_c : pointer to temporary distribution
          3) rho_c    : pointer to densities
          4) u_c      : pointer to velocities
          5) n_c      : number of phases

  OUTPUT: void
 */
{
  int cc, cn, ck, cd; /* counter: phase, node, direction, dimension */
  int fg_ind, u_ind, rho_ind; /* index to node values for the distribution, velocity, density */
  real *fg, *fg_tmp, *rho, *u; /* Pointers to one phase for distr., tmp dist., density, velocity*/

  struct Step *step = sys->step;

  for (cc = 0; cc < field->n_c; ++cc) /* Phases */
  {
    /* Pointers to one phase */
    fg = &field->f[step->f_phase*cc];
    fg_tmp = &field->f_tmp[step->f_phase*cc];
    rho = &field->rho[step->rho_phase*cc];
    u = &field->u[step->u_phase*cc];
    /* index initiation */
    fg_ind = 0;
    u_ind = 0;
    rho_ind = 0;

    for (cn = 0; cn < sys->max_n_tot; ++cn) /* Nodes */
    {
      /* Initiate density */
      rho[rho_ind] = 0.0;
      /* Initiate velocity */
      for (cd = 0; cd < sys->n_D; ++cd)
        u[u_ind + cd] = 0.0;
      /* Initiate distribution */
      for (ck = 0; ck < sys->n_Q; ++ck)
      {
        fg[fg_ind + ck] = 0.0;
        fg_tmp[fg_ind + ck] = 0.0;
      }

      /* Update indeces */
      fg_ind += step->f_node;
      u_ind += step->u_node;
      rho_ind += step->rho_node;
    }
  }
}


/* ---------------------------------------------------------------------------------  init_system_wall */
void init_system_wall_attr_real(real *attr_c, int n_c, int stp_comp, int nr_wall)
/*
init_system_wall :
    sets a wall attribute to zero:

  INPUT  : 1) attr_c : attribute array
           2) n_c    : number of phases for the attribute for each node

  OUTPUT : void
 */
{
  int cc, cw; /*counter components, wall nodes */
  real *attr; /* Pointer the the attribute for on phase */

  for (cc = 0; cc < n_c; ++cc)
  {
    /* Pointer to one phase*/
    attr = &attr_c[stp_comp*cc];

    for (cw = 0; cw < nr_wall; ++cw)
    {
      attr[cw] = 0.0;
    }
  }
}


void init_system_wall_attr_int(int *attr_c, int n_c, int stp_comp, int nr_wall)
/*
init_system_wall :
    sets a wall attribute to zero:

  INPUT  : 1) attr_c : attribute array
           2) n_c    : number of phases for the attribute for each node

  OUTPUT : void
 */
{
  int cc, cw; /*counter components, wall nodes */
  int *attr; /* Pointer the the attribute for on phase */

  for (cc = 0; cc < n_c; ++cc)
  {
    /* Pointer to one phase*/
    attr = &attr_c[stp_comp*cc];

    for (cw = 0; cw < nr_wall; ++cw)
    {
      attr[cw] = 0;
    }
  }
}


void init_system_wall_tmp_real(real *attr, int n_c)
/*
init_system_wall :
    sets a wall attribute to zero:

  INPUT  : 1) attr_c : attribute array
           2) n_c    : number of phases for the attribute for each node

  OUTPUT : void
 */
{
  int cc;

  for (cc = 0; cc < n_c; ++cc)
  {
    attr[cc] = 0.0;
  }
}


void init_system_wall_tmp_int(int *attr, int n_c)
/*
init_system_wall :
    sets a wall attribute to zero:

  INPUT  : 1) attr_c : attribute array
           2) n_c    : number of phases for the attribute for each node

  OUTPUT : void
 */
{
  int cc;

  for (cc = 0; cc < n_c; ++cc)
  {
    attr[cc] = 0;
  }
}


/* ---------------------------------------------------------------------------------  init_phase_forces */
/* void init_phase_forces(real *flforce, real *glforce, int n_c) */
//void init_phase_forces(struct Field *fluid, struct System *sys)
///*
//init_phase_forces :
//	set the phase forces to zero
//
//  INPUT  : 1) flforce : fluid-fluid forces
//           2) glforce : rock-fluid forces
//		   3) n_c     : number of phases
//
//  OUTPUT : void
// */
//{
//  int cc, cn, cd; /* Counter phase, node, dimension */
//  int fgl_ind;    /* index to force node */
//  real *flforce_c, *glforce_c; /* forces for one phase */
//
//  struct Step *step = sys->step;
//
//  for( cc = 0; cc < fluid->n_c; ++cc )
//  {
//    flforce_c = &fluid->flforce[step->u_phase*cc];
//    glforce_c = &fluid->glforce[step->u_phase*cc];
//
//    fgl_ind = 0;
//
//    for( cn = 0; cn < sys->max_n_tot; ++cn )
//    {
//      for( cd = 0; cd < sys->n_D; ++cd )
//        flforce_c[fgl_ind + cd] = glforce_c[fgl_ind + cd] = 0.0;
//
//      fgl_ind += step->u_node;
//    }
//  }
//}

