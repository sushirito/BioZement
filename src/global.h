#ifndef GLOBAL_H
#define GLOBAL_H

#include <cstdio>
#include <cerrno>
#include <cstring>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <climits>
#include <cfloat>
#include <string>
#include <iomanip>
#include <map>
#include <sys/types.h>
#if defined (_WIN32)
#include <direct.h>
#else
#include <unistd.h>
#endif
#include <sys/stat.h>
#include <bitset>
#include "defines.h"
#include "options.h"
#include "Input.h"
#include "macros.h"
#include "mpi.h"

#define _2D_RUN_
#define DIFFUSIVE_RUN

struct InitChem;
struct splayTree;
struct treeNode;
struct BasVec;
struct Node;
struct Mpi;
struct Output;


/********************/
/* M I N E R A L S */
/********************/
struct Minerals {
  struct Mineral {
    char name[30];
    real rho;
    real *sfrac;
    real *delta;
    real *rest;
    real *old_sfrac;
    real *dm;
    double cm3_mol;
    double SA;
    double SA_weight;
  } *list;
  double SA_tot;
  int n_tot;
  double *dm_sum;
  double *sfrac_sum;
  double init_molvol;
  int no_comb;
  double *sum_sfrac, *sum_sfrac_old, *sum_sfrac_prev;
  double *buffer;
#ifdef _MAGN_ON_MAGN_
  char *nucleation_site;
#endif
};

/********************/
/*  L I N K L I S T */
/********************/
struct Links {
  struct Link {
    int group;
    int wall;
    int fluid;
    unsigned char dir;
    unsigned char inert;
    unsigned char ghost;
    double pH;
    real delta[10];
    real delta_sum[10];
  } *list;
  struct Group {
    int wall;
    int first;
    int last;
  } *group;
  int nlinks;
  int ngroups;
  int max_links;
  int max_groups;
};

/********************/
/*     S T E P      */
/********************/
struct Step {

  /* stepping */
  int f_phase  , f_node;   /* Steps in the velocity distribution */
  int rho_phase, rho_node; /* Steps in density */
  int u_phase  , u_node;     /* Steps in the velocity  */
  int *f_rel_neig;      /* relative neighbor index for distributions */
  int n_phase;               /* Steps in nodes  */
  int *n_rel_neig;       /* relative neighbor index for nodes */
  int stp_e_k;               /* Steps in the velocity/position basis */
  int neg_x_dir;
  int pos_x_dir;
  int wall_comp;       /* Steps in the wall components */
};

/********************/
/*     N O D E      */
/********************/
struct Node {
  /* nodes */
  char *mode;     /* the mode of the node:  */
  char *Rmask;    /* AMIN mask for source. 1 source on, 0: off*/
  char *Bmask; // mask for biomass field. 1 bacterial colony, 0:no bacteria /* AMIN */
  int nr_wall;        /* number of wall nodes */
  int max_nr_wall;        /* initial number of wall nodes */
  int in_x, out_x;  // x-position of inlet and outlet nodes

  std::bitset<NUM_MODE> is_fluid;
  std::bitset<NUM_MODE> is_fluid_in_out;

  std::vector<int> isolated;
  std::vector<int> perc_cluster;
  int perc_inlet;
  //int *perc_list;  // list over nodes in percolating cluster
  //int num_perc;    // length of perc_list
  
  int nr_per;         /* number of periodic nodes */
  int nr_in;       /* number of inlet nodes */
  int nr_out;      /* number of outlet nodes */
  //int nr_isolated;  /* number of isolated nodes */

  /* lists of special nodes types */
  int *list_per;    /* list of periodic points */
  int *list_wall;   /* list of wall nodes */
  int *list_in;  /* list of inlet nodes */
  int *list_out; /* list of outlet nodes */

  /* _MPI_ stuff */
  int nr_ghost;       /* number of ghost nodes */
  char *is_ghost;     /* 1 if node is ghost  */
  char *is_ghost_or_periodic;     /* 1 if node is ghost or periodic  */
  int *list_ghost;    /* list of ghost nodes */

  /* Volume of fluid */
  Links *links;  /* list of links connecting a wall- and a fluid-node */


  /* BIOZEMENT */
  char *inert;

};


/********************/
/*   S Y S T E M    */
/********************/
struct System {
  std::map<std::string, std::string> files;
  char drive_type[30];
  double *drive;
  double *vel_target;
  
  double rho_not_converged;
  int node_not_converged;
  
  Output *output;

  Options *opt;

  real rate_constant;
  //double rate_shift;
  //double plug_mass_shift;
  //double source_flux;  // used in collision_g_newconc()
  
  double dt, dt_b, dx;
  double Dlb, D, D_b, Dlb_b, flux, nu_lb, nu;
  //double dt_ratio;
  double porosity;
  double pore_volume_tot, pore_volume_reac; // total (incl. inlet/outlet), and reactive pore volume
  double gx;
  int L_lb, H_lb, W_lb, A_lb, V_lb; // length, height, width, cross-section, volume
                                    // of REACTIVE system (i.e. excluding inlet/outlet) in lb-units 
  
  int real_units_in_input; // the first line in inp.dat must read REAL UNITS if it
                           // contains settings in real units to cause set_parameters()
                           // to be called in main.c
  int *is_inlet; /* array with 1 at inlet nodes and 0 otherwise */
  int fluid;
  int chem;    
  int t;   /* timesteps */
  double time; // in seconds
  int nwrite; // counts output files
  int steps_to_wait_before_converge;
  double t_extra;
  int n_reinit, vel_reinit_steps;
  double max_vel, mean_vel;
  int const_flow;
  double conv_crit_chem, conv_crit_ss, conv_crit_ss_min;
  double conv_crit_vel, conv_crit_vel_target;
  int conv_length;
  double conv_jump;
  int steps_until_check;
  int splaytree_res;
  int convergence;
  int last_convergence;
  int species_not_converged;

  /* dimensions */
  int n_D;                 /* number of dimensions */
  int n_Q;                 /* number of grid directions */
  int *max_n;              /* number of nodes in each spatial direction */
  int *MAX_N;              /* num nodes */
  int max_n_tot;           /* total number of nodes */

  int start_itr, end_itr; /* starting and end iteration */
  int n_itr;               /* number of iterations */
  double write_int;      /* interval for writing data to file IN SECONDS */

  /* the basis */
  real *w;      /* Weights for the different directions */
  int *ei;         /* position basis */
  real *ev;        /* velocity basis */
  int *k_bb;       /* Bounce back directions */

  /* constants */
  real eq_cst[3]; /* constants in the equilibrium distribution */
  real cs_2;       /* square of the lattice velocity */

  /* stepping */
  Step *step;
  Node *node;
  
  int velreinit, rho_reinit;
  //int *g_reinit, n_g_reinit;
  int *new_nodes, *new_nodes_cpy, *new_nodes_rest;
  int num_new_nodes;
  int *new_solid_nodes, num_new_solid_nodes;
  int solid2fluid, fluid2solid;

//  /* Flooding rate protocol  */
//  real ** vel_max_time; /* velocity and time */
//  int  num_vel_changes; /* number of vel changes */
//  int current_vel;
  
  /* helpful variables */
  int *ev_ind;
  int num_negative;

  Mpi *mpi;

  // used by perm measurement (e.g. during diffusive simulation)
  Node *perm_node;

  int charge_balance;
  
  
};

/********************/
/*    F I E L D     */
/*                  */
/* (bulk properties)*/
/********************/
struct Field{
  
  const char *name[30];

  /* constants */
  int n_c;          /* number of phases */
  int ntot;
  double rho_init0, rho_inlet0;
  real *tau;        /* relaxion time for each phase */
  double dt;

  /* parameters */
  real *gravity;           /* gravity */
  real *D;         /* diffusion constant */
  real C;         /* concentration at boundary */
  double q_source;
  double u_Darcy; // used in _SPECIES_SOURCE_BULK_
  double grad_P; // imposed pressure drop over the system
  
  /* Microscopic */
  real *f;         /* distribution function */
  real *f_tmp;     /* temporary distribution function used in the propagation step */

  /* Macroscopic */
  real *rho;       /* density */
  real *rho_old;     /* previous value used to compute convergence */
  real *u;         /* velocity */
  real *u_old;         /* used in velocity convergence calculation */
  real *psi;      /* the relative density */

  real *u_tot;     /* barisentrick velocity */
  real *rho_tot;   /* total density */

  /* AMIN EJE */
  real * R_local;

  /* initial values */
  real *rho_init;    /* initial density for each phase */
//  real *u_init;     /* initial velocity for each phase */
  //real *u_inlet;     /* velocity at the inlet */
  real *rho_inlet;   /* density at the inlet for each phase */
#ifdef _FLUID_BDRY_PRESS_
  real *rho_outlet;   /* fluid density at the outlet */
#endif
#ifdef _USE_COLOR_GRAD_
  //double inlet_pressure;
  //  double surface_tension;
  double beta;
  std::vector< std::vector<double> > surface_tension_multiphase_matrix;
  std::vector<double> rho_wall;
#endif
  std::vector<double> u_mean;
  std::vector<double> phase_u_mean;
  std::vector<double> phase_perm;
  std::vector<double> nu_lb;

  /* temp variable  */
  real *u_eq;      /* holds the local velocity used in calculate the equilibrium distribution */
  real *nc_tmp; /* holds the distribution for a nodal phase */ 

  /* stepping */
  Step *step;

  // fluid specific
  double flux, area, mean_velx, perm;

  // species specific
  double *effluent_conc;

  double *sum_conc, *sum_conc_old, *sum_conc_prev;
  double *mean_conc, *mean_conc_old;
  double *mean_g; 
  //  double *saturated_nodes; // used to calculate saturation
  std::vector<double> saturation;
  int not_converged = 1;
  double rho_max = 0, rho_min = 0;
  //double gx_change = 0;

//#ifdef _SPECIES_SOURCE_BULK_
  std::vector< std::vector<double> > rho_inlet_file;
//#endif
};


/********************/
/* B O U N D A R Y  */
/********************/
struct Boundary{
  /* -- -- outlet */
  int *dir_in;   /* directions in to the boundary nodes */     
  int *dir_out;  /* directions out from the boundary nodes */
  int n_dir_inout;     /* number of direction in and out of the boundary */

  char *w_mod;     /* wall mode (not yet used) */

  /* wall chemistry */
  int *rho_s_wall;   /* wall components in the solid */
  real *rho_l_wall;  /* wall components in the liquid */
  real *rho_eq_wall; /* chemical equlibrium at the wall */
  real *flux_l_wall; /* wall fluxes for the liquid components */
  real *old_flux_l_wall; /* wall fluxes for the liquid components */
  /* NB! rho_g_ holds the total amount of compoonents lost or gained at the wall */

  /* temp variables */
  int *rho_s_wall_tmp;   /* wall components in the solid */
  real *rho_l_wall_tmp;  /* wall components in the liquid */
  real *rho_eq_wall_tmp; /* chemical equlibrium at the wall */
  real *rate_wall_tmp;   /**/

  /* constants */
  int n_liquid;      /* number of liquid wall species (same as n_c_g_?) */
  int n_solid;      /* number of solid wall species */

  real J_influx;   /* total influx */

  /* stepping */ 
  //int stp_comp;   MOVED TO STEPS         /* Steps in the wall components */

};


/********************/
/*      M P I       */
/********************/
struct Mpi {
  struct Ghost {
    int *dir;
    int ndir;
  } *ghost;

  int my_rank;   /* rank of all processes */
  int eff_rank;  /* rank of processes in effluent calculation, >= 0 or MPI_UNDEFINED */
  int nr_procs;
  int lower_bound[3];
  int upper_bound[3];
  double time;
  int np[3];  // number of processes in each direction
  int ind_rank[3]; // rank index
  int *flag_rbuf;
  int eff_root_rank;
  int pbc_sides;

  // effluent communication
  MPI_Comm EFF_COMM;
  MPI_Group world_group, eff_group;

  // send and receive buffers
  double *double_rbuf;
  int sbuf_size;
  double *sbuf_eff, *rbuf_eff;  // for effluent calculations

  int *cpy_slice_LR;
  int *cpy_slice_TB;
  int *cpy_slice_FB;
  int *cpy_slice[3];
  int *cpy_slice_char_LR;
  int *cpy_slice_char_TB;
  int *cpy_slice_char_FB;
  int *cpy_slice_char[3];
  int cpy_slice_size[3];
  int *cpy_line_X;
  int *cpy_line_Y;
  int *cpy_line_Z;
  int *cpy_line[3];
  int *cpy_line_char_X;
  int *cpy_line_char_Y;
  int *cpy_line_char_Z;
  int *cpy_line_char[3];
  int cpy_line_size[3];

  MPI_Datatype slice[3];
  MPI_Datatype line[4];
  MPI_Datatype slice_char[3];
  MPI_Datatype line_char[3];
  MPI_Datatype slice_double[3];
  MPI_Datatype line_double[3];
  //MPI_Datatype slice_vel[3];
  //MPI_Datatype line_vel[3];

  /* These are used in communiate_init() and probably
     don't need to be here... */
  int S[3];
  int S_ncom;                /* number of sides to communicate */
  int S_ncpy;                /* number of sides to copy */
  int S_com[3];              /* which sides to communicate */
  int S_cpy[3];              /* which sides to copy */
  int side_proc[2][3];
  int side_proc_nobc[2][3];
  int side     [3][2];
  int side_bc  [3][2];
  int side_str[3];
  int side_str_bc[3];
  int side_rcv[3][2];
  int side_snd[3][2];
  int side_rcv_scalar[3][2];
  int side_snd_scalar[3][2];
  //int side_snd_vel[3][2];

  int E[3];
  int E_ncom;
  int E_com[6];
  int E_ncpy;
  int E_cpy[6];
  int edge_proc[2][6];
  int edge_proc_nobc[2][6];
  int edge     [3][4];
  int edge_bc  [3][4];
  int edge_str[3][2];
  int edge_str_bc[3][2];
  int edge_rcv[3][4];
  int edge_snd[3][4];
  int edge_rcv_scalar[3][4];
  int edge_snd_scalar[3][4];


  int C;
  int C_ncom;
  int C_com[4];
  int C_ncpy;
  int C_cpy[4];
  int crnr_proc[2][4];
  int crnr_proc_nobc[2][4];
  int crnr     [4][2];
  int crnr_bc  [4][2];
  int crnr_str[4];
  int crnr_str_bc[4];
  int crnr_rcv[4][2];
  int crnr_snd[4][2];
  int crnr_rcv_scalar[4][2];
  int crnr_snd_scalar[4][2];

//  real *ghost_sfrac_diff;
//  real *ghost_sfrac_tmp;

};


//------------------------------------------------//
//                                                //
//  F U N C T I O N    D E C L A R A T I O N S    //
//                                                //
//------------------------------------------------//


//------------------------------------------------//
//                allocate.c                      //
//------------------------------------------------//
int allocate(System *sys);
int allocate_fg(Field *field, System *sys);
int allocate_wall_attr_real(real **attr, int n_c, int nr_wall);
int allocate_wall_attr_int(int **attr, int n_c, int nr_wall);
int allocate_wall_tmp_real(real **attr, int n_c);
int allocate_wall_tmp_int(int **attr, int n_c);
int allocate_vel(real **u, int max_n_tot, int n_D);
int allocate_rho(real **rho, int max_n_tot);
//int allocate_phase_forces( Field *fluid,  System *sys);
int allocate_outlet( Boundary*,  System*);


//------------------------------------------------//
//                basis_calc.c                    //
//------------------------------------------------//
void calc_cs_2(int*,  System*);
int check_basis(int*,  System*);
void init_k_bb( System*);
void calc_fg_eq_c( System*);


//------------------------------------------------//
//            boundary_conditions.c               //
//------------------------------------------------//
void bc_3D_copy_gradient(int nr_nodes, int *node_list,  Field *field,  System *sys, int nb_dir, int n1_dir);
void bc_press_red_blue_fancy(int pos, double rho_boundary_red, double rho_boundary_blue, int nx,  Node *node,  Field *fluid,  System *sys);
void bc_press_red_blue(int pos, double rho_boundary_red, double rho_boundary_blue, int nx,  Node *node,  Field *fluid,  System *sys);
void bc_press_red_blue_z(int pos, double rho_boundary_red, double rho_boundary_blue, int nz,  Field *fluid,  System *sys);
void bc_press_f_hack(int pos, double rho_boundary, int nx,  Field *fluid,  System *sys);
void bc_press_f_hack_z(int pos, double rho_boundary, int nz,  Field *fluid,  System *sys);
void bc_init_run_f( Node *node,  Field*,  System*);
void bc_run_f( Node *node,  Field*,  System*);
void bc_init_run_g( InitChem*,  Field*,  Field*,  System*,  Minerals*,
     Boundary*,  splayTree **);
void bc_init_run_g_newconc( Field*,  Field*,  System*);
void bc_run_g( InitChem*,  Field*,  Field*,  System*,  Minerals*,
     Boundary*,  splayTree **);


//------------------------------------------------//
//        calc_equilibrium_conditions.c           //
//------------------------------------------------//
double calc_feq(int ck, real rho, real *u,  System *sys);
double calc_geq(int ck, real rho, real *u,  System *sys);


//------------------------------------------------//
//                calc_rho_vel.c                  //
//------------------------------------------------//
void correction_rho_vel_color_grad_fluid( Node *node,  Field *fluid, System *sys);
void correction_rho_vel_color_grad( Node *node,  Field *fluid,  Field *species,  System *sys);
void calc_rho_vel_f_color_grad( Node *node,  Field *fluid,  System *sys);
void calc_rho_vel(real*, real*, real*,  System*);
void calc_rho_vel_newconc(real*, real*, real*, real*, real*,  System*);
void calc_rho_vel_loop(double *f_dist,  Field *fluid,  System *sys, double *q);
void calc_rho(real*, real*,  System*);
void calc_mom(real*, real*,  System*);
void calc_rho_mom_all( Field*,  System*);
//void shanchen_f( Field*,  System*);
void calc_rho_vel_cn_f(int, real*,  Field*,  System*);


//------------------------------------------------//
//                collision.c                     //
//------------------------------------------------//
double DeltaOmegaF(int ck, double tau_eff_sym, double tau_eff_ant, double *u, double *F, struct System *sys);
double cos_phi(int ck, double * grad, struct System *sys);
double Omega_surface_tension(int ck, double * grad, double A, double * B, struct System *sys);
void color_gradient(int cn, double *grad, double *rho_cc1, double *rho_cc2, struct System *sys);
void collision_f_color_grad( Node *node,  Field *fluid,  System *sys);
void collision_f( Node *node,  Field*,  System*);
void collision_f_newconc( Node *node,  Field*,  System*);
//void collision_shanchen_f( Field*,  System*);
void collision_g( Field*,  Field*,  System*);
void collision_g_newconc( Field*,  Field*,  System*);
void collision_g_newconc_dif(struct Field *dif, struct Field *species, struct Field *fluid, struct System *sys, struct InitChem *ICS);
void get_mean_g( System*,  Field*);
double get_porosity( System *sys);
void set_rho_in_percolation_nodes(int cn, double *rho, double rho_set,  System *sys);
void set_rho_inlet_from_file(Field *species, System *sys);


//------------------------------------------------//
//                 communicate.c                  //
//------------------------------------------------//
int init_communicate( System *sys);
void communicate(real*, int,  System*);
void communicate_node( Node*,  System*);
void communicate_char( char*,  System*);
void communicate_double(double*, System *);
void communicate_sfrac( System*,  Minerals*);
void communicate_delta( System*,  Minerals*);
void communicate_mineral_delta(int,  System*,  Minerals*);
void communicate_mineral_sfrac(int,  System*,  Minerals*);
void communicate_rho_tot(Field *, System *);
void communicate_rho(Field *, System *);
void communicate_rho_tot(Field *, System *);
void setup_ghost_nodes( Node*,  System*);
void add_ghost_node( Node *node, int node_nr, int *c);
void ghost_loop( Node *node,  System *sys, int *a, int *s, int *c, int start, const char *str);
void set_ghost_neighbors( System*);
int global_to_local_ind(int Nn,  System *sys);
int global_to_local_xind(int Nx,  System *sys);
int get_local_ind(int, int, int,  System*);
int get_global_ind(int, int, int,  System*);
int get_global_xcoord(int,  System*);
int get_global_ycoord(int,  System*);
int get_global_zcoord(int,  System*);
//void broadcast_flag( Mpi *mpi, int *flag);
void broadcast_global_max( Mpi *mpi, double *max);
void broadcast_global_min( Mpi *mpi, double *min);
int get_global_max( Mpi *mpi, double *max);
void get_local_coord(int nn, int *ind,  System *sys);
void get_global_coord(int nn, int *ind,  System *sys);


//------------------------------------------------//
//                file_functions.c                //
//------------------------------------------------//
void check_file(const char *filename, std::ifstream &file);
FILE* my_fopen(const char*, const char*);
int getword(char *s, char *t);
int isdbl(char *s);
double atodbl(char *s);
int atoint(char *s);
int isint(char *s);
int cleanline(char s[]);
int fgetline(FILE *fp, char s[]);


//------------------------------------------------//
//                init_run.c                      //
//------------------------------------------------//
void init_run_g( Field*,  Field*,  System*);
void init_run_g_newconc( Field*,  Field*,  System*);
void init_run_3D_g( Field*,  Field*,  System*);
void init_run_f( Node *node,  Field*,  System*);
void init_run_f_nodes( Field*,  System*, int*, int);
void init_mpi_run( System *sys, int *argc, char ***argv,  Field *species);
void init_chem_splay_mineral( System *sys,  Field *species,  InitChem *InitChemSolver,
     splayTree ***st_lst,  BasVec **Vchem_key_tmp,  Minerals *minerals);
void init_3D_run(System *sys,  Field *species,  Field *fluid,  Field *dif,  InitChem **InitChemSolver,
     splayTree ***st_lst,  Minerals *minerals,  Boundary *bndry,  BasVec **Vchem_key_tmp);
void init_3D_fluid_run(System *sys,  Field *fluid);
void init_3D_chemistry_run( System *sys,  Field *species,  Field *fluid,  InitChem *InitChemSolver,
     splayTree **st_lst,  Minerals *minerals,  Boundary *bndry);
//void init_2D_run( System *sys,  Field *species,  Field *fluid,  InitChem *InitChemSolver,
//     splayTree **st_lst,  Minerals *minerals,  Boundary *bndry);
void steady_state_specie_run( Field *,  Field *,  System *);
void run_steady_state_all_species_and_fluid( Field *species,  Field *fluid,  Minerals *minerals,  System *sys);
//void find_isolated_fluid_nodes( Field *,  Field *,  Boundary *,  System *,  InitChem *,  splayTree **,  Minerals *);
void find_isolated_fluid_nodes( Field *,  Field *,  System *);
void set_isolated_nodes_inert( Node *node);
void init_ss_species_run( Field *,  Field *,  System *);
void end_ss_species_run( Field *,  System *);
void init_solid_fraction( System*,  Field*,  Minerals*,  InitChem*);
double get_dist_sqrd(int n1, int n2,  System *sys);
#ifdef _RANDOM_MINERAL_
void find_wall_neighbours(int nn0, int nn, double *conc, double rad, int *ncalls,  System *sys,  Node *node,  Minerals *minerals);
#endif



//------------------------------------------------//
//                init_system.c                     //
//------------------------------------------------//
void init_structs(Input &,  Options *,  Step **,  Minerals **,  Links **,  Node **,
     System *,  Field **,  Field **,  Field**, Boundary **);
void init_system( System *);
void init_system_fg( Field *,  System *);
void init_system_wall_attr_real(real*, int, int, int);
void init_system_wall_attr_int(int*, int, int, int);
void init_system_wall_tmp_real(real*, int);
void init_system_wall_tmp_int(int*, int);
void init_phase_forces( Field*,  System*);


//------------------------------------------------//
//                inlet_bc.c                      //
//------------------------------------------------//
void inlet_bc_f( Field*,  System*,  Boundary*);
void inlet_bc_g( Field*,  Field*,  System*,  Boundary*);
void inlet_leaky_bc_g( Field *species,  Field *fluid,  System *sys,  Boundary *bndry);
void inlet_bc_3D_copy_fg( Field *field,  Boundary *bndry,  System *sys, int direction);


//------------------------------------------------//
//              main_functions.c                  //
//------------------------------------------------//
char *get_cmd_option(char **begin, char **end, const std::string &option);
void read_cmd_args(int argc, char *argv[], System *sys);
int init_variables_and_geometry( System *sys,  Boundary *bndry,  Field *fluid,  Field *species,  Field *dif);
int allocate_init_wall( Boundary*,  System*);
int allocate_init_f( Field*,  System*);
int allocate_init_g( Field *fluid,  System *sys);
real convergence(real * f, real * f_tmp,  System *sys);
void write_error_to_file(real * rho, real * psi, real * rho_fluid, real * u_fluid,  System * sys, char fn[]);
void set_bulk_conc_after_reinit( Field *species,  Field *fluid,  System *sys);
void write_effluent_file( System*,  Field*,  Field*,  Minerals*,  InitChem*, time_t);
void calc_effluent_and_perm( System *sys,  Node *node,  Field *species,  Field *fluid, int broadcast);
double get_pressure_gradient( System *sys,  Field *fluid);
void print_status( System *sys,  InitChem *ICS, time_t t_begin,  splayTree **st_lst, int no_comb);
void print_splaytree_size( System *sys, int no_comb,  splayTree **st_lst);
int species_has_converged(int t,  System *sys,  Field *species, double conv_crit, double *diff);
void set_init_bulk_conc(double* rho_init,  Field *species,  Field *fluid,  System *sys);
void set_conc_field( Field *species,  Field *fluid,  System *sys);
void set_species_in_new_fluid_nodes( Field *species,  Field *fluid,  System *sys);
//void set_init_bulk_conc_from_file( Field *species,  System *sys);
char *get_time_string(double time);
void print_time(double);
//void read_ss_species_from_file( InitChem *,  Field *,  Field *,  System *,  Minerals *,
//     Boundary *,  splayTree **);
void advance_one_step( InitChem *,  Field *,  Field *,  Field *,  System *,  Minerals *,
     Boundary *,  splayTree **);
void print_status_and_write_output(struct InitChem &ICS_pH, Output &output,  System *sys, time_t t_begin, int no_comb,  InitChem *ICS,
     Field *species,  Field *fluid,  Minerals *minerals,  splayTree **st_lst, int);
void set_parameters   ( System *sys,  Field *species, Field *dif,  Field *fluid,  InitChem *ICS);
void set_parameters_LB( System *sys,  Field *species,  Field *fluid,  InitChem *ICS);
void print_run_parameters( System *sys,  Field *species,  Field *fluid,  InitChem *ICS,  Minerals *minerals);
//void read_fields_from_file( System *sys,  Field *species,  Field *fluid);
double get_mean_velocity(Node *node, System *sys, Field *fluid);
void get_max_min_rho(Node *node, System *sys, Field *fluid);
double get_flux( Node*,  System *,  Field *);
double get_umax( Node*,  System *,  Field *);
double get_rhomax( Node*,  System *,  Field *);
void get_pore_volume( Node *node, int *nfluid_tot, int *nfluid_reac,  System *sys);
void get_saturation( Node *node,  Field *fluid,  System *sys);
void get_solid_massfraction( Minerals *minerals,  System *sys);
void get_fluid_conc( Field *species,  System *sys);
void get_reactive_surface_area( Node *node, double *sumW,  System *sys,  Minerals *minerals);
int node_is_reactive(int nn,  System *sys);
void get_slice_surface( System *sys, int xpos, double *surf, int *nlinks);
void save_fluid_restart_files( Field *fluid,  System *sys);
void save_chem_restart_files( Field *species,  Minerals *minerals,  System *sys);
void load_fluid_restart_files( Field *fluid,  System *sys);
void load_chem_restart_files( Field *species,  System *sys);
void load_mineral_restart_files( Minerals *minerals,  System *sys);
void reset_rho_f_in_new_fluid_nodes( Field *field,  System *sys);
void reset_rho_f_in_new_solid_nodes( Field *field,  System *sys);
void update_rho(int *node_list, int num_nodes,  Field *fluid,  Field *species,  System *sys);
int copy_weighted_f_from_fluid_neighbors(int nn,  Field *species,  System *sys);
void apply_boundary_conditions(InitChem *ICS, Field *species, Field *dif, Field *fluid, System *sys, Minerals *minerals,
    Boundary *bndry, splayTree **st_lst);


//------------------------------------------------//
//             moving_boundaries.c                //
//------------------------------------------------//
void check_convergence_v2( System*,  Minerals*,  Field*,  Boundary*,  Links*);
void check_convergence(int *not_converged,  System *sys,  Field *species, double conv_crit, double *diff);
void check_mean_convergence(int *not_converged,  System *sys,  Field *species, double conv_crit, double *diff);
int init_linklist( Node*,  System*);
void update_linklist( Node*,  System*);
void update_solid_v5( System*,  Minerals*,  Links*);
void update_nodes( System*,  Boundary*,  Minerals*);
void has_node_changed(int,  Minerals *,  Node *, int *, int *);
void reinit_velocity( Node*,  System*,  Field*, int);
void initialize_velocity( Node*,  System*,  Field*, int);
void broadcast_sum_of_int_to_all_procs(int *flag,  Mpi *mpi);
void update_rho_old( System *sys,  Field *species);
void reset_rho_vel_in_new_fluid_nodes( Field *fluid,  System *sys);
int check_velocity_convergence(Node *node, Field *fluid, System *sys, double conv_crit);
int check_velocity_convergence_nphase(Node *node, Field *fluid, System *sys);
void init_fluid( Node*,  Field*,  System*);
void advance_fluid(Node *node,  Field *fluid,  System *sys);
void apply_fluid_bc(Node *node, Field *fluid, System *sys);
void init_permeability_measurement( System *sys,  Field *fluid,  Field *species,  Minerals *minerals);
void measure_permeability( System *sys,  Field *fluid,  Field *species,  Minerals *minerals);
void write_fluid_distribution_and_vector_to_file( System *sys,  Field *fluid);
int advance_fluid_to_convergence(double conv_crit,  Node *node,  Field *fluid,  System *sys, double starttime, int print_status);
int advance_fluid_and_species_to_convergence(double conv_crit,  Node *node,  Field *fluid,
					      Field *species, Field *dif,  InitChem *ICS,  Minerals *minerals,
					      Boundary *bndry,  splayTree **st_lst,  System *sys,
					     double starttime, int print_status);
void update_u_old( Node *node,  Field *fluid,  System *sys);
//int check_for_correct_vel_or_flux_and_adjust_gravity( System *sys,  Field *fluid, double starttime, int i);
int check_value_and_adjust_drive(double value, double target,  System *sys,  Field *fluid);
void calculate_permeability( System *sys,  Node *node,  Field *fluid,  Field *species);
void write_permeability_to_file(int timestep,  System *sys,  Field *fluid,  Minerals *minerals);
void copy_node_mode( Node *from_node,  Node *to_node,  System *sys);


//------------------------------------------------//
//                   outlet_bc.c                  //
//------------------------------------------------//
void outlet_bc_g( Field *species,  Field *fluid,  System *sys,  Boundary *bndry);
void outlet_bc_3D_copy_g( Field*,  Field*,  Boundary*,  System*);
void outlet_bc_3D_copy_fg( Field *field,  Boundary *bndry,  System *sys, int direction);


//------------------------------------------------//
//                  periodic.c                    //
//------------------------------------------------//
void add_periodic_bc( System*);
void setup_periodic_bc( Node*,  System*);
int find_sn(int,  System*);
void periodic_bc( Node *node,  Field*,  System*);


//------------------------------------------------//
//                 propagate.c                    //
//------------------------------------------------//
void propagate( Node *node,  Field*,  System*);
void propagate_dif( Node *node,  Field*,  System*); /* AMIN */


//------------------------------------------------//
//                   read_geo.c                   //
//------------------------------------------------//
void init_geometry_and_nodes( System *sys);
void read_geofile( Node *node,  System *sys);
void check_size_of_file( System *sys);
void set_wall_nodes( Node *node,  System *sys);
void copy_inlet_outlet_slice(Node *,  System *);
void add_overlapping_slice_at_inlet_and_outlet(Node *, System *);
void set_inlet_outlet_nodes( Node *node,  System *sys);
void init_node_lists( Node *node,  System *sys);
void add_outer_walls( Node *node,  System *sys);
void add_raster( Node *node,  System *sys);
void skip_header_line(FILE *fp,  Mpi *mpi);
void eof_error( Mpi *mpi);
void set_periodic_nodes_2D( Node *node,  System *sys);
void close_inlet_outlet( Node *node,  System *sys);
void set_dimensions_from_geofile( System *sys);
void update_dim_for_rotations(int &nx, int &ny, int &nz,  System *sys);
void update_dim_for_extra_walls(int add_sub, int &nx, int &ny, int &nz,  System *sys);
void set_rho_in_percolation_nodes(int cn, double *rho, double rho_set,  System *sys);
void set_system_dimensions(System *sys);
void get_x_limits( Node *node, int *xmin, int *xmax,  System *sys);


//------------------------------------------------//
//             read_write_fields.c                //
//------------------------------------------------//
void write_dist_binary(char*, real*, int,  System *);
void read_dist_binary(char*, real*, int,  System *);
void write_1d_array(char *fname, double *data, int n);
void read_1d_array(char *fname, double *data, int n);


//------------------------------------------------//
//                   read.c                       //
//------------------------------------------------//
int is_empty(const char *s);
void read_input(char *,  InitChem **,  System *,  Field *,  Field *,  Boundary *, char **);
void read_input_v2(Input &,  InitChem **,  System *,  Field *,  Field *,  Field *,  Boundary *, char **);
void myprint_double_vec(const char *vec_name, double *vec, int vec_size);
void myprint_int_vec(const char *vec_name, int *vec, int vec_size);
void myread_basis(const char *fname,  System *sys);
void myread_double_vec(double **vec , int n_e, char *ct, FILE *fp, int maxline);
void myread_int_vec(int **vec , int n_e, char *ct, FILE *fp, int maxline);
long myreadline_indx(char *line, int maxline, FILE *input_file);
char *myreadline_loc(char *line, int maxline, FILE *input_file);
#ifdef _FLUID_BDRY_PRESS_
void set_rho_inlet_outlet(Node *node, Field *fluid);
#endif


//------------------------------------------------//
//             wall_bc_gca_ej.c                   //
//------------------------------------------------//
void wall_bc_bounce_back( Node *node,  Field *fluid,  System *sys);
void wall_bc_mid_grid_gca_ej_g_v3( InitChem *ICS,  Field *species,  Field *fluid,
     System *sys,  Minerals *minerals,  Boundary *bndry,  splayTree **st_lst);
int calc_eq_rho_rate_jlv( System *sys,  Boundary *bndry,
     InitChem *ICS,  Minerals *minerals, int n_c, real norm_fac,
			  int ck, real *min_conc, int nn,  treeNode **tn,  splayTree **st_lst);
int calc_eq_rho_rate_newconc( System *sys,  Boundary *bndry,
     InitChem *ICS,  Minerals *minerals,  real fluid_f,
    int n_c, real norm_fac, int ck, real *min_conc, int nn,  treeNode **tn,  splayTree **st_lst);


#endif
