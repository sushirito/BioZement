#include "global.h"
#include "chem_global.h"
#include "output.h"
#include "biozementfunctions.h"

char *get_cmd_option(char **begin, char **end, const std::string &option)
{
  char **itr = std::find(begin, end, option);
  if (itr != end && ++itr != end) {
    return *itr;
  }
  return 0;
}

void read_cmd_args(int argc, char *argv[], System *sys)
{
  // default path and name of input files and output directory
  sys->files["basis"] = "./program_files/basis";
  sys->files["data"] = "./program_files/data";
  sys->files["input"] = "input.dat";
  sys->files["options"] = "options.dat";
  sys->files["out"] = "out";
  sys->files["geo"] = "./program_files/geo/dot_diss.dat"; // geo.dat";
  sys->files["path"] = "/home/amin/Documents/Codes/Qt/biozement/run/biozement";
  sys->files["species_rho_in"] = "";
  sys->files["Rfile"] = "./program_files/geo/dot_diss_sink.dat";
  sys->files["Bfile"] = "./program_files/geo/dot_diss_Biomass.dat";

  //for (std::map<std::string, std::string>::iterator it = sys->files.begin(); it != sys->files.end(); ++it) {
  for (auto & map : sys->files) {
    char *value = get_cmd_option(argv, argv + argc, "-" + map.first);
    if (value)
      map.second = value;
  }

  sys->files["input"] = sys->files["path"] + "/" + sys->files["input"];
  sys->files["options"] = sys->files["path"] + "/" + sys->files["options"];
  sys->files["out"] = sys->files["path"] + "/" + sys->files["out"];


  //for (std::map<std::string, std::string>::iterator it=sys->files.begin(); it!=sys->files.end();++it) {
  //  std::cout << it->first << " : " << it->second << std::endl;
  //}
}

/* ---------------------------------------------------------------------------------  allocate_init_general */
int init_variables_and_geometry(struct System *sys, struct Boundary *bndry, struct Field *fluid, struct Field *species, struct Field *dif)
/*
allocate_init_parallel :
    allocates and initiate the general variables for a parallel run.
    Including reading the geometry.

  INPUT  : 1) input_file : file name to the input file

  OUTPUT : (boolean) GOOD: managed to allocate memory.
                     BAD : didn'n manage to allocate memory
*/
{
  struct Node *node = sys->node;
  //int measure_perm = (sys->opt->perm_interval>0);

  if (sys->mpi->my_rank == 0) {
    std::cout << std::endl << "Input from \'" << sys->files["input"] << "\', \'" << sys->files["options"]
	      << "\' and \'" << sys->files["geo"] << "\'" << std::endl;
    std::cout << "Output folder is \'" << sys->files["out"] << "\'" << std::endl;
#ifdef _2D_RUN_
	std::cout << "This is a 2D run" << std::endl << std::endl;
#else
	std::cout << "This is a 3D run" << std::endl << std::endl;
#endif
//#ifdef _WIN32
//	std::cout << "Running on Windows" << std::endl;
//#endif
  }

  /* ALLOCATE MEMORY */
  /* -- Add space for a periodic rim/ghost nodes */
  add_periodic_bc(sys); /* AMIN */

  init_communicate(sys);

  /* -- general variables */
  if (allocate(sys) == BAD) {
    printf("error in allocate\n");
    exit(1);
  }

  /* SETUP THE SYSTEM */
  /* -- general initiations */
  init_system(sys);

  // periodic boundary conditions
  setup_periodic_bc(node, sys);

  //  if (sys->n_D<3) {
  //    sys->max_n[2]=1;
  //  }

    // setup and initialize ghost-nodes for mpi using system dimensions
  setup_ghost_nodes(node, sys);
  set_ghost_neighbors(sys);

//#ifdef _FLUID_BDRY_PRESS_
  sys->node->out_x = sys->MAX_N[0] - 1 - sys->node->out_x;
//#endif

  // read geometry file geo.dat and initialize node-lists and node-counters
  init_geometry_and_nodes(sys);
  
  allocate_outlet(bndry, sys);

  /* -- wall */
  sys->step->wall_comp = node->nr_wall;

#ifdef _DEBUG_MPI_
  printf("DIM <%02d (%d,%d,%d)> max_n: (%d, %d, %d)(%d), lb: (%3d, %3d, %3d), ub: (%3d, %3d, %3d), MAX_N: (%d, %d, %d)\n",
    sys->mpi->my_rank,
    sys->mpi->ind_rank[0], sys->mpi->ind_rank[1], sys->mpi->ind_rank[2],
    sys->max_n[0], sys->max_n[1], sys->max_n[2], sys->max_n_tot,
    sys->mpi->lower_bound[0], sys->mpi->lower_bound[1], sys->mpi->lower_bound[2],
    sys->mpi->upper_bound[0], sys->mpi->upper_bound[1], sys->mpi->upper_bound[2],
    sys->MAX_N[0], sys->MAX_N[1], sys->MAX_N[2]);
  printf("NODE <%02d> nr_in: %d, nr_out: %d, nr_wall: %d, nr_per: %d, nr_ghost: %d\n",
    sys->mpi->my_rank, node->nr_in, node->nr_out, node->nr_wall, node->nr_per, node->nr_ghost);
#endif // _DEBUG_MPI_

  /* list over wall-fluid links */
  if (init_linklist(node, sys) == BAD) {
    printf("Unable to allocate memory in init_linklist()\n");
    exit(1);
  }


  // init array that will hold a list over previous WALL nodes with an undefined 
  // species-distribution. New dist copied from a neigboring FLUID node
  int num_new_nodes = std::max(10 * node->nr_wall, 2000);
  if (!(sys->new_nodes = (int *)malloc(num_new_nodes * sizeof(int)))) return BAD;
  if (!(sys->new_solid_nodes = (int *)malloc(num_new_nodes * sizeof(int)))) return BAD;
  if (!(sys->new_nodes_cpy = (int *)malloc(num_new_nodes * sizeof(int)))) return BAD;
  if (!(sys->new_nodes_rest = (int *)malloc(num_new_nodes * sizeof(int)))) return BAD;

#ifdef FLUID_OFF
  sys->flux = -1.0;
  sys->fluid = 0;
  //fluid->gravity[0] = 0.0;
#endif

  if (sys->mpi->my_rank == 0) {
    if (sys->fluid == 0)
      printf("\n****  Fluid is off\n");
    if (sys->flux < 0)
      printf("****  Pure diffusive simulation");
    //#ifdef _SOURCE_
    //    printf(" with source term");
    //#endif
    if ((sys->fluid == 0) || (sys->flux < 0))
      printf("\n\n");
  }

  // allocation and initializations
  if (sys->fluid) {
    allocate_init_f(fluid, sys);   // sets fg_, fg_tmp_, rho, u to zero
  }
  if (sys->chem) {
    allocate_init_g(dif, sys);  // AMIN EJE
    allocate_init_g(species, sys);
    allocate_init_wall(bndry, sys);
  }

  return GOOD;
}





/* ---------------------------------------------------------------------------------  allocate_init_f */
int allocate_init_f(struct Field *fluid, struct System *sys)
/*
  allocate_init_f :
  allocates and initiate the variables related to the fluid.

  INPUT  : void
  OUTPUT : (boolean) GOOD: managed to allocate memory.
  BAD : didn'n manage to allocate memory
*/
{
  if (allocate_fg(fluid, sys) == BAD)
  {
    printf("error in allocate memory for the fluid component\n");
    exit(1);
    return BAD;
  }

  if (fluid->n_c==1) {
    // make u_tot the same as u if only one fluid phase is present
    fluid->u_tot = &fluid->u[0];
  } else {
    /* Allocate memory for the baricentric velocity and total density */
    if (allocate_vel(&fluid->u_tot, sys->max_n_tot, sys->n_D) == BAD)
    {
      printf("error in allocate memory for the fluid component\n");
      exit(1);
      return BAD;
    }
  }

  if (allocate_rho(&fluid->rho_tot, sys->max_n_tot) == BAD)
  {
    printf("error in allocate memory for the fluid component\n");
    exit(1);
    return BAD;
  }

  /* Allocate memory for the fluid-fluid and fluid-rock forces */
//	if( allocate_phase_forces(fluid,sys) == BAD )
//	{
//		printf("error in allocate memory for the fluid component\n");
//		exit(1);
//		return BAD;
//	}

    /* Initiate system */
  init_system_fg(fluid, sys);

  /* Initiate phase forces
  init_phase_forces(fluid, sys); */

  return GOOD;
}


/* ---------------------------------------------------------------------------------  allocate_init_g */
int allocate_init_g(struct Field *species, struct System *sys)
/*
allocate_init_g :
    allocates and initiate the variables realted to the chemical species.

  INPUT  : void
  OUTPUT : (boolean) GOOD: managed to allocate memory.
                     BAD : didn'n mange to allocate memory
*/
{

  int cc, cn; /* counter: phase, node, direction, dimension */
  int  psi_ind; /* index to node values for the distribution, velocity, density */
  real *psi; /* Pointers to one phase for distr., tmp dist., density, velocity*/
  real *R_local;

  struct Step *step = sys->step;



  if (allocate_fg(species, sys) == BAD)
  {
    printf("error in allocate memory for the diffusive component\n");
    exit(1);
    return BAD;
  }

  /* Add the molality */
  if (!(species->psi = (real *)calloc(sys->max_n_tot * species->n_c, sizeof(real)))) return BAD;

  /* AMIN EJE */
  if (!(species->R_local = (real *)calloc(sys->max_n_tot * species->n_c, sizeof(real)))) return BAD;


  init_system_fg(species, sys);

  for (cc = 0; cc < species->n_c; ++cc) /* Phases */
  {
    /* Pointers to one phase */
    psi = &species->psi[step->rho_phase*cc];
	R_local = &species->R_local[step->rho_phase*cc]; /* AMIN EJE */
    /* index initiation */
    psi_ind = 0;

    for (cn = 0; cn < sys->max_n_tot; ++cn) /* Nodes */
    {
      /* Initiate density */
      psi[psi_ind] = 0.0;
      // if (cc == 2){
      //     if (sys->node->Rmask[psi_ind])
      //         R_local[psi_ind] = 0.0; /* AMIN EJE R_local (mol/dt)*/
      //     else{
      //         R_local[psi_ind] = 0.0;
      //     }
      // }

      /* Update indeces */
      psi_ind += step->rho_node;
    }
  }

  return GOOD;
}

/* ---------------------------------------------------------------------------------  allocate_init_wall */
int allocate_init_wall(struct Boundary *bndry, struct System *sys)
/*
allocate_init_wall :
    allocates and initiate the variables related to the wall.

  INPUT  : void
  OUTPUT : (boolean) GOOD: managed to allocate memory.
                     BAD : didn'n manage to allocate memory
*/
{

  struct Node *node = sys->node;
  /* Allocate */
  /* -- fluid components */
  /* -- -- current density */
  if (allocate_wall_attr_real(&bndry->rho_l_wall, bndry->n_liquid, node->max_nr_wall) == BAD) return BAD;
  if (allocate_wall_tmp_real(&bndry->rho_l_wall_tmp, bndry->n_liquid) == BAD) return BAD;
  /* -- -- equlibrium  density */
  if (allocate_wall_attr_real(&bndry->rho_eq_wall, bndry->n_liquid, node->max_nr_wall) == BAD) return BAD;
  if (allocate_wall_tmp_real(&bndry->rho_eq_wall_tmp, bndry->n_liquid) == BAD) return BAD;
  /* -- -- fluxes */
  if (allocate_wall_attr_real(&bndry->flux_l_wall, bndry->n_liquid, node->max_nr_wall) == BAD) return BAD;
  //if( allocate_wall_attr_real(&bndry->old_flux_l_wall, bndry->n_liquid, node->max_nr_wall) == BAD ) return BAD;
  /* -- rates */
  //if( allocate_wall_tmp_real(&bndry->rate_wall_tmp, bndry->n_liquid) == BAD ) return BAD;
  /* -- solid components */
  //if( allocate_wall_attr_int(&bndry->rho_s_wall, bndry->n_solid, node->max_nr_wall) == BAD ) return BAD;
  //if( allocate_wall_tmp_int(&bndry->rho_s_wall_tmp, bndry->n_solid) == BAD ) return BAD;

  /* Initiate */
  /* -- fluid components */
  /* -- -- current density */
  init_system_wall_attr_real(bndry->rho_l_wall, bndry->n_liquid, sys->step->wall_comp, node->max_nr_wall);
  init_system_wall_tmp_real(bndry->rho_l_wall_tmp, bndry->n_liquid);
  /* -- -- equlibrium  density */
  init_system_wall_attr_real(bndry->rho_eq_wall, bndry->n_liquid, sys->step->wall_comp, node->max_nr_wall);
  init_system_wall_tmp_real(bndry->rho_eq_wall_tmp, bndry->n_liquid);
  /* -- rates */
  //init_system_wall_tmp_real(bndry->rate_wall_tmp, bndry->n_liquid);
  /* -- solid components */
  //init_system_wall_attr_int(bndry->rho_s_wall, bndry->n_solid, sys->step->wall_comp, node->max_nr_wall);
  //init_system_wall_tmp_int(bndry->rho_s_wall_tmp, bndry->n_solid);

  return GOOD;
}

real convergence(real * f, real * f_tmp, struct System *sys)
{
  int cn, nn;  /* phase counter, current node, neighbor node*/
  int ck; /* direction */
  char *mode;
  real max_error, current_error, error_nevner, error_teller;
  //int max_error_node;
  real f_nn, f_cn_tilde;
  //real tau;
  //int ntmp, nx, ny, nz;


  /* init varaibles */
  struct Step *step = sys->step; /* Steps in array */
  mode = sys->node->mode; /* Node mode ie. fluid, solid et c. */
  max_error = 0;

  /* for all fluid nodes */
  for (cn = 0; cn < sys->max_n_tot; ++cn) {
    if (mode[cn] == FLUID) {
      /* check that both x and x + e_alpha node is a fluid node */
      for (ck = 0; ck < sys->n_Q; ++ck) {
        nn = cn - step->n_rel_neig[ck];
        if (mode[nn] == FLUID) {
          //f_cn = f_tmp[cn*step->f_node + ck];
          f_cn_tilde = f[cn*step->f_node + ck];
          f_nn = f_tmp[nn*step->f_node + ck];
          /* check distance from stationary state */
          error_teller = abs(f_nn - f_cn_tilde);
          error_nevner = 0.5*abs(f_nn + f_cn_tilde);
          if (error_nevner < 1e-10)
            current_error = 0;
          else
            current_error = error_teller / error_nevner;
          if (current_error > max_error) {
            max_error = current_error;
            //max_error_node = cn;
            //max_f = f_cn;
            //max_f_tilde = f_cn_tilde;
          }
        }
      }
    }
  }

  return max_error;
}


void write_error_to_file(real * rho, real * psi, real * rho_fluid, real * u_fluid, struct System * sys, char fn[])
{
  FILE * f_ptr;
  real psi_min, psi_max, rho_min, rho_max;
  real rho_fluid_min, rho_fluid_max, u_max;
  char *mode;
  int cn;
  int u_step;

  mode = sys->node->mode;
  u_step = sys->step->u_node;

  /* open file */
  f_ptr = fopen(fn, "a");
  if (f_ptr == NULL) {
    printf("ERROR: could not open file %s\n", fn);
    exit(1);
  }

  psi_max = rho_max = rho_fluid_max = u_max = 0.0;
  psi_min = rho_min = rho_fluid_min = DBL_MAX;

  for (cn = 0; cn < sys->max_n_tot; ++cn) {
    if (mode[cn] == FLUID) {
      rho_min = std::min(rho_min, rho[cn]);
      rho_max = std::max(rho_max, rho[cn]);
      psi_min = std::min(psi_min, psi[cn]);
      psi_max = std::max(psi_max, psi[cn]);
      rho_fluid_min = std::min(rho_fluid_min, rho_fluid[cn]);
      rho_fluid_max = std::max(rho_fluid_max, rho_fluid[cn]);
      u_max = std::max(u_max, u_fluid[cn*u_step]); /* maximum x-component */
    }
  }

  /* print to file */
  fprintf(f_ptr, "%d %d %g %g %g %g %g %g %g\n", sys->max_n[0], sys->max_n[1], rho_min - 1.0, rho_max - 1.0, psi_min - 1.0, psi_max - 1.0,
    rho_fluid_min - 1.0, rho_fluid_max - 1.0, u_max);

  /* close file */
  fclose(f_ptr);

}



/* void read_fields_from_file(struct System *sys, struct Field *species, struct Field *fluid) */
/* /\************************************ */
/*  ************************************\/ */
/* { */
/*   /\* Read fields (binary) *\/ */
/*   char fn[50]; */

/*   if( sys->chem ) { */
/*     sprintf(fn, "species->f"); */
/*     read_dist_binary(fn, species->f, species->n_c*sys->max_n_tot*sys->n_Q); */

/*     sprintf(fn, "species->f_tmp"); */
/*     read_dist_binary(fn, species->f_tmp, species->n_c*sys->max_n_tot*sys->n_Q); */
/*   } */

/*   if( sys->fluid ) { */
/*     sprintf(fn, "fluid->f"); */
/*     read_dist_binary(fn, fluid->f, fluid->n_c*sys->max_n_tot*sys->n_Q); */

/*     sprintf(fn, "fluid->f_tmp"); */
/*     read_dist_binary(fn, fluid->f_tmp, fluid->n_c*sys->max_n_tot*sys->n_Q); */
/*   } */
/* } */


void print_status_header(struct System *sys)
{
  if (sys->mpi->my_rank == sys->mpi->nr_procs - 1) {
#ifdef _DEBUG_
    printf("\nfilenr  step  jump-step: time = dd-hh:mm:ss (cs), cpu = dd-hh:mm:ss (cs), reinit: ( nr [ steps ], s2f, f2s ), not-conv: nr \n");
    printf("----------------------------------------------------------------------------------------------------------------------\n");
#else
    printf("\nfilenr   step                 time = dd-hh:mm:ss     cpu = dd-hh:mm:ss     chem = nr. solutions\n");
    printf("------------------------------------------------------------------------------------------------------\n");
#endif
    fflush(stdout);
  }
}

void print_status(struct System *sys, struct InitChem *ICS, time_t t_begin, struct splayTree **st_lst, int no_comb)
  /************************************
   ************************************/
{
  struct Mpi *mpi = sys->mpi;
  int solutions[20] = { 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 };
  int rbuffer[20];
  if (no_comb > 20) {
    printf("ERROR in write_status: no_comb > 20\n");
    MPI_Finalize();
    exit(1);
  }

#ifdef _DEBUG_
  // sum check values
  double mean_steps = 0.0;
  int num_not_conv = 0;//, s2f=0;
  int num_neg = 0;
  // total number of times charge balance did not converge
  if (sys->chem) {
    MPI_Reduce(&ICS->num_chrg_not_conv, &num_not_conv, 1, MPI_INT, MPI_SUM, mpi->nr_procs - 1, MPI_COMM_WORLD);
#ifdef SET_NEG_RHO_ZERO
    // and total number of times g-dist < 1e-8 was set to 1e-8
    MPI_Reduce(&(sys->num_negative), &num_neg, 1, MPI_INT, MPI_SUM, mpi->nr_procs - 1, MPI_COMM_WORLD);
  }
#endif
#endif // _DEBUG_
  if (sys->chem) {
    for (int i = 0; i < no_comb; ++i) {
      solutions[i] += (int)st_lst[i]->size;
    }
    MPI_Reduce(&solutions, &rbuffer, no_comb, MPI_INT, MPI_SUM, mpi->nr_procs - 1, MPI_COMM_WORLD);
  }

  // print status
  if (mpi->my_rank == mpi->nr_procs - 1) {
    // real (simulated) time
    double step_sum = (double)(sys->t + 1) + sys->t_extra;
    printf("\r%4d:    step = %9.3e     time = ", sys->nwrite, step_sum);
    print_time(sys->time + sys->dt);
    // cpu time
    printf("     cpu = ");
    print_time(difftime(time(NULL), t_begin));
#ifdef _DEBUG_
    if (sys->n_reinit > 0)
      mean_steps = (double)sys->vel_reinit_steps / (double)sys->n_reinit;
    else
      mean_steps = 0.0;
#ifdef SET_NEG_RHO_ZERO
    printf(", reinit: (%3d [%7.1e], %4d, %4d), not-conv: %2d, neg: %.0e", sys->n_reinit, mean_steps, sys->solid2fluid, sys->fluid2solid, num_not_conv, (double)num_neg);
#else
    printf(", reinit: (%3d [%7.1e], %4d, %4d), not-conv: %2d", sys->n_reinit, mean_steps, sys->solid2fluid, sys->fluid2solid, num_not_conv);
#endif
#endif //_DEBUG_
    if (sys->chem) {
      printf("     chem = ");
      for (int i = 0; i < no_comb; ++i) {
        printf("%d", (int)(rbuffer[i] / mpi->nr_procs));
        if (i < no_comb - 1)
          printf(", ");
      }
    }
    fflush(stdout);
  }

  // reset flags
  if (sys->chem) {
    ICS->num_chrg_not_conv = 0;
    sys->num_negative = 0;
  }
  sys->n_reinit = 0;
  sys->solid2fluid = sys->fluid2solid = 0;
  sys->vel_reinit_steps = 0;
}

char *get_time_string(double time)
{
  /****************************************************
   * INPUT: time as double of seconds size of splaytree for individual processes *
   * OUTPUT: string of days, hours, mins, secs, and centi-seconds *
   ****************************************************/
  int days, hours, mins, secs;
  double dtmp, sec_fractional;
  static char string[50];

  sec_fractional = modf(time, &dtmp);
  secs = (int)dtmp;

  days = secs / 86400;
  secs %= 86400;
  hours = secs / 3600;
  secs %= 3600;
  mins = secs / 60;
  secs %= 60;
  //  sprintf(string, "%02d-%02d:%02d:%02d (%02d)",
  // 	  days, hours, mins, secs, (int)(sec_fractional*100));
  sprintf(string, "%02d-%02d:%02d:%02d",
    days, hours, mins, secs);
  return string;
}

void print_time(double time) {
  /****************************************************
   * INPUT: time as double of seconds size of splaytree for individual processes *
   * OUTPUT: string of days, hours, mins, secs, and centi-seconds *
   ****************************************************/
  int days, hours, mins, secs;
  double dtmp, sec_fractional;
  //int unit = 100; // 100: centi-seconds, 1000: milli-seconds

  sec_fractional = modf(time, &dtmp);
  secs = (int)dtmp;

  days = secs / 86400;
  secs %= 86400;
  hours = secs / 3600;
  secs %= 3600;
  mins = secs / 60;
  secs %= 60;
  //  printf("%02d-%02d:%02d:%02d (%02d)",
  //	 days, hours, mins, secs, (int)(sec_fractional*unit));
  printf("%02d-%02d:%02d:%02d", days, hours, mins, secs);
}

void print_splaytree_size(struct System *sys, int no_comb, struct splayTree **st_lst)
/****************************************************
 * print size of splaytree for individual processes *
 ****************************************************/
{
  char fn[50];
  FILE *fp;
  int i;

  sprintf(fn, "%s/%04d_splaytree.dat", sys->files["out"].c_str(), sys->mpi->my_rank);
  fp = my_fopen(fn, "a+b");
  fprintf(fp, "%10d ", sys->t + 1);
  for (i = 0; i < no_comb; ++i) {
    fprintf(fp, "%7lu", st_lst[i]->size);
    if (i == no_comb - 1)
      fprintf(fp, "\n");
    else
      fprintf(fp, ", ");
  }
  fclose(fp);
}


void write_effluent_file(struct System *sys, struct Field *species, struct Field *fluid,
  struct Minerals *minerals, struct InitChem *ICS, time_t t_begin)
  /**************** WRITE EFFLUENT ********************* *
   *******************************************************/
{
  char fn[50];
  FILE *fp = NULL;
  struct Mpi *mpi = sys->mpi;
  static int write_header = 1;
  int nfluid_tot, nfluid_reac;
  double phi = 0.0, SA = 0.0;
  double dx3_m3 = sys->dx*sys->dx*sys->dx; // voxel volume in m3
  double dx3_cm3 = 1.0e6*dx3_m3; // 1 m3 = 1e6 cm3
  double dx3_L = 1.0e3*dx3_m3; // 1 m3 = 1e3 L
  double sf_2_mol;

  if (mpi->eff_rank == 0)
    sprintf(fn, "%s/rho_effluent.dat", sys->files["out"].c_str());

  // write header if first call
  if (write_header && (mpi->eff_rank == 0)) {
    write_header = 0;
    fp = fopen(fn, "w");
    if ((fp == NULL) && (sys->mpi->my_rank == 0)) {
      printf("ERROR in write_effluent_file: fopen returned %d\n", errno);
      MPI_Finalize();
      exit(0);
    }
    fprintf(fp, "# %-11s %-11s %-11s", "timestep", "realtime", "CPUtime");
    // write fluid measurement names
    if (sys->fluid) {
      fprintf(fp, " %-11s %-11s %-11s %-11s %-11s %-11s", "gx", "area", "vel[m/s]", "perm[D]", "SA[m2]", "phi");
    }

    // write species and mineral measurement names
    if (sys->chem) {
      for (int cc = 0; cc < species->n_c; ++cc) {
        fprintf(fp, " %-15s mean_%-10s sum_%-11s diff_%-10s ddt_%-11s",
            species->name[cc], species->name[cc], species->name[cc], species->name[cc], species->name[cc]);
      }
      for (int min = 0; min < minerals->n_tot; ++min) {
        fprintf(fp, " sum_%-11s diff_%-10s ddt_%-11s SA_%-12s SA_w_%-10s",
            minerals->list[min].name, minerals->list[min].name, minerals->list[min].name, minerals->list[min].name, minerals->list[min].name);
      }
    }
    if (sys->fluid) {
      for (int cc = 0; cc < fluid->n_c; ++cc) {
        fprintf(fp, " sat_fluid%d  ", cc);
      }
#ifdef _USE_COLOR_GRAD_
      fprintf(fp, " rho_min     ");
      fprintf(fp, " rho_max     ");
      fprintf(fp, " Ca          ");
#endif
      fprintf(fp, " velx_mean    ");
      for (int cc = 0; cc < fluid->n_c; ++cc) {
        fprintf(fp, " velx_fluid%d   ", cc);
        //char comp[] = "xyz";
        //for (int d = 0; d < sys->n_D; ++d) {
        //fprintf(fp, " vel%c_fluid%d ",comp[d], cc);
        //}
      }
      for (int cc = 0; cc < fluid->n_c; ++cc) {
        fprintf(fp, " perm_fluid%d  ", cc);
      }
    }
    fprintf(fp, "\n");
    fclose(fp);
  }

  calc_effluent_and_perm(sys, sys->node, species, fluid, 0);

  get_pore_volume(sys->node, &nfluid_tot, &nfluid_reac, sys);
  phi = (double)nfluid_reac / (double)(sys->V_lb);

  if (sys->chem)
    get_reactive_surface_area(sys->node, &SA, sys, minerals);
  SA *= sys->dx*sys->dx;

  get_fluid_conc(species, sys);
  if (sys->chem)
    get_solid_massfraction(minerals, sys);

  if (sys->fluid) {
    get_saturation(sys->node, fluid, sys);
	get_mean_velocity(sys->node, sys, fluid);
    get_max_min_rho(sys->node, sys, fluid);
    double lb_to_D = M2_TO_DARCY*sys->dx*sys->dx; // 1 m2 = 1.013249966e12 Darcy, 1 Darcy = 1000 mD
    for (int cc=0; cc<fluid->n_c; ++cc) {
      double u_Darcy = fluid->phase_u_mean[cc*sys->n_D]*fluid->saturation[cc]*double(nfluid_tot)/double(sys->V_lb);
      fluid->phase_perm[cc] = 0.0;
      if (fluid->gravity[0]>0.0) {
	fluid->phase_perm[cc] = lb_to_D*fluid->nu_lb[cc]*u_Darcy/fluid->gravity[0];
      }
    }
  }

  // write flow and effluent data
  if (mpi->eff_rank == 0) {
    fp = fopen(fn, "a");
    // time data
    fprintf(fp, "  %11d %11.5e %11.5e",
      sys->t + 1,
      sys->time + sys->dt,
      difftime(time(NULL), t_begin));
    // flow data
    if (sys->fluid) {
      fprintf(fp, " %11.5e %11.5e %11.5e %11.5e %11.5e %11.5e",
          fluid->gravity[0],
          std::max(fluid->area, species->area),
          fluid->mean_velx,
          fluid->perm,
          SA,
          phi);
    }
    // species concentrations
    if (sys->chem) {
      for (int cc = 0; cc < species->n_c; ++cc) {
        fprintf(fp, " %15.8e %15.8e %15.8e %15.8e %15.8e",
            species->effluent_conc[cc],
            species->mean_conc[cc],
            species->sum_conc[cc] * dx3_L,
            (species->sum_conc[cc] - species->sum_conc_old[cc])*dx3_L,
            (species->sum_conc[cc] - species->sum_conc_prev[cc])*dx3_L);
        // store current value
        species->sum_conc_old[cc] = species->sum_conc[cc];
      }
      // mineral concentrations
      for (int min = 0; min < minerals->n_tot; ++min) {
        sf_2_mol = dx3_cm3 / minerals->list[min].cm3_mol;
        fprintf(fp, " %15.8e %15.8e %15.8e %15.8e %15.8e",
            minerals->sum_sfrac[min] * sf_2_mol,
            (minerals->sum_sfrac[min] - minerals->sum_sfrac_old[min])*sf_2_mol,
            (minerals->sum_sfrac[min] - minerals->sum_sfrac_prev[min])*sf_2_mol,
            minerals->list[min].SA / minerals->SA_tot,
            minerals->list[min].SA_weight / minerals->SA_tot);
        // store current value
        minerals->sum_sfrac_old[min] = minerals->sum_sfrac[min];
      }
    }
    // fluid saturations
    if (sys->fluid) {
      for (int cc = 0; cc < fluid->n_c; ++cc) {
        fprintf(fp, " %11.5e", fluid->saturation[cc]);
      }
#ifdef _USE_COLOR_GRAD_
      fprintf(fp, " %13.6e", fluid->rho_min);
      fprintf(fp, " %13.6e", fluid->rho_max);
      fprintf(fp, " %13.6e", fluid->u_mean[0]*fluid->nu_lb[0]/fluid->surface_tension_multiphase_matrix[0][1]);
#endif
      fprintf(fp, " %14.6e", fluid->u_mean[0]);
      for (int cc = 0; cc < fluid->n_c; ++cc) {
	fprintf(fp, " %14.6e", fluid->phase_u_mean[cc*sys->n_D]);
	//if (cc==0) 
	//  printf("\nERROR: %.8e - %.8e =  %.8e\n", fluid->u_mean[cc], (fluid->dt*fluid->mean_velx/sys->dx),
	//	 fluid->u_mean[cc] - (fluid->dt*fluid->mean_velx/sys->dx));
	//printf("fluid->u_mean: %.8e\n", fluid->u_mean[cc*sys->n_D]);
	//for (int d = 0; d < sys->n_D; ++d) {
        //  fprintf(fp, " %13.6e", fluid->u_mean[cc*sys->n_D + d]);
        //}
      }
      for (int cc = 0; cc < fluid->n_c; ++cc) {
	fprintf(fp, " %14.6e", fluid->phase_perm[cc]);
      }
    }
    fprintf(fp, "\n");
    fclose(fp);
  }
  
//  if (!sys->chem) {
//    species->n_c = n_c;
//  }

}



void calc_effluent_and_perm(struct System *sys, struct Node *node, struct Field *species, struct Field *fluid, int broadcast)
// ---------------------------
//
//
// ---------------------------
{
  struct Mpi *mpi = sys->mpi;
  struct Step *step = sys->step;
  int *n = sys->max_n;
  int cc, j, k, nn;
  //double m2_to_D = 1.013250e12;
  double lb_to_D = M2_TO_DARCY*sys->dx*sys->dx*sys->nu_lb; // 1 m2 = 1.013249966e12 Darcy, 1 Darcy = 1000 mD
  int xpos = n[0] - 1 - sys->opt->eff_offset;
  int A_lb;
  double *conc;
  int area = 0, vel = 1, chem = 2;

  if (sys->chem && species != NULL) {
#ifdef NEWCONC
    conc = species->psi;
#else
    conc = species->rho;
#endif
  }

  // only for processes in the effluent communication group
  if (mpi->eff_rank != MPI_UNDEFINED) {

    // zero out send-buffer
    for (cc = 0; cc < mpi->sbuf_size; ++cc)
      mpi->sbuf_eff[cc] = 0.0;

    // loop nodes to sum velocity, area and concentration
    for (j = 1; j < n[1] - 1; j++) {
      for (k = 1; k < n[2] - 1; k++) {
        nn = xpos + j*n[0] + k*n[0] * n[1];

        //if (node->mode[nn]==FLUID || node->mode[nn]==OUTLET || node->mode[nn]==INLET) {
        if (node->is_fluid_in_out[node->mode[nn]]) {
          // area
          mpi->sbuf_eff[area] += 1.0;
          // velocity
          if (sys->fluid)
            mpi->sbuf_eff[vel] += fluid->u_tot[step->u_node*nn];
	  //mpi->sbuf_eff[vel] += fluid->u[step->u_node*nn];
          // concentration
          if (sys->chem && species != NULL) {
            for (cc = 0; cc < species->n_c; ++cc) {
              mpi->sbuf_eff[chem + cc] += conc[step->rho_phase*cc + step->rho_node*nn];
            }
          }
        }

      }
    }

    // sum values from all effluent processes
    MPI_Reduce(mpi->sbuf_eff, mpi->rbuf_eff, mpi->sbuf_size, MPI_DOUBLE, MPI_SUM, 0, mpi->EFF_COMM);
  }

#ifdef _FLUID_BDRY_PRESS_
  double grad_P = get_pressure_gradient(sys, fluid);
  //  if (sys->mpi->my_rank==0) {
  //    std::cout << "measured grad_P - imposed grad_P = " << grad_P << " - " << fluid->grad_P*sys->cs_2 << " = "
  //        << grad_P - fluid->grad_P*sys->cs_2 << std::endl;;
  //}
#endif

        // calculate properties
  if (mpi->eff_rank == 0) {

    // reset
    if (sys->chem && species != NULL) {
      for (cc = 0; cc < species->n_c; ++cc) {
        species->effluent_conc[cc] = 0.0;
      }
    }

    if (sys->fluid) {
      fluid->mean_velx = fluid->perm = 0.0;
      // cross section area in nodes**2
      fluid->area = mpi->rbuf_eff[area];
      if (fluid->area > 0.0) {
        // velocity in m/s
        if (sys->opt->straight_tube) {
          A_lb = fluid->area;
        } else {
          A_lb = sys->A_lb;
        }
        fluid->mean_velx = mpi->rbuf_eff[vel] * sys->dx / fluid->dt / A_lb;
        // permeability in D, k_lb = Q*mu/(A*f) = u_mean*nu*rho/f, where f = grad(P)
        if (fluid->gravity[0] > 0)
          fluid->perm = lb_to_D*mpi->rbuf_eff[vel] / A_lb / fluid->gravity[0];
#ifdef _FLUID_BDRY_PRESS_
        if (fluid->grad_P > 0) {
          //double grad_P = fluid->grad_P*sys->cs_2;
          fluid->perm = lb_to_D*mpi->rbuf_eff[vel]/ A_lb / grad_P;
          //printf("fluid->perm = lb_to_D*mpi->rbuf_eff[vel]/ A_lb / grad_P = %.3e*%.3e/%d/%.3e = %.3e\n",
          //    lb_to_D, mpi->rbuf_eff[vel], A_lb, grad_P, fluid->perm);
        }
#endif
        //fluid->perm = lb_to_mD*mpi->rbuf_eff[vel]/fluid->area/fluid->gravity[0];
    //double A_gradP = sys->pore_volume_tot*fluid->gravity[0]/(sys->MAX_N[0]-2);
        //fluid->perm = lb_to_mD*mpi->rbuf_eff[vel]/A_gradP;
      }
    }

    if (sys->chem && species != NULL) {
      // effluent concentrations
      species->area = mpi->rbuf_eff[area];
      if (species->area > 0.0) {
        for (cc = 0; cc < species->n_c; ++cc) {
          species->effluent_conc[cc] = mpi->rbuf_eff[chem + cc] / species->area;
        }
      }
    }
  }

  if (broadcast) {
    if (sys->fluid) {
      MPI_Bcast(&(fluid->area), 1, MPI_DOUBLE, mpi->eff_root_rank, MPI_COMM_WORLD);
      MPI_Bcast(&(fluid->mean_velx), 1, MPI_DOUBLE, mpi->eff_root_rank, MPI_COMM_WORLD);
      MPI_Bcast(&(fluid->perm), 1, MPI_DOUBLE, mpi->eff_root_rank, MPI_COMM_WORLD);
    }
    if (sys->chem && species != NULL) {
      MPI_Bcast(&(species->area), 1, MPI_DOUBLE, mpi->eff_root_rank, MPI_COMM_WORLD);
    }
  }

}

//----------------------------------
//
//----------------------------------
double get_pressure_gradient(struct System *sys, struct Field *fluid)
{
  int in_x = sys->opt->copy_in_out+1;  // node->in_x;
  int out_x = sys->MAX_N[0]-2-sys->opt->copy_in_out; //  node->out_x;
  int x_pos[2] = {in_x, out_x};
#ifdef _DEBUG_
  double P_min[2] = {10,10}, P_max[2] = {0,0};
#endif
  double sbuff[4] = {0,0,0,0}, rbuff[4] = {0,0,0,0};
  struct Node *node = sys->node;

  for (int cn = 0; cn < sys->max_n_tot; ++cn) {
    if (node->is_ghost_or_periodic[cn])
      continue;
    for (int i=0; i<2; ++i) {
      if ( (get_global_xcoord(cn, sys)==x_pos[i]) && node->is_fluid_in_out[node->mode[cn]]) {
        sbuff[2*i  ] += fluid->rho[cn];
        sbuff[2*i+1] += 1.0;
#ifdef _DEBUG_
        if (fluid->rho[cn]>P_max[i]) {
          P_max[i] = fluid->rho[cn];
        }
        if (fluid->rho[cn]<P_min[i]) {
          P_min[i] = fluid->rho[cn];
        }
#endif
      }
      //      if (get_global_xcoord(cn, sys)==x_pos[i]) {
      //        sys->node->mode[cn] = FLUID;
      //      }
    }
  }
  //sys->output->write_geo_file("test_geo", sys->node, sys);

#ifdef _DEBUG_
  for (int i=0; i<2; ++i) {
    broadcast_global_min(sys->mpi, &P_min[i]);
    broadcast_global_max(sys->mpi, &P_max[i]);
  }
#endif
  MPI_Reduce(sbuff, rbuff, 4, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Bcast(rbuff, 4, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  for (int i=0; i<2; ++i) {
    sbuff[i] = rbuff[2*i]/rbuff[2*i+1];
  }

#ifdef _DEBUG_
  double dPdx1 = (P_max[0]-P_min[1])*sys->cs_2/(out_x-in_x);
  double dPdx2 = (P_min[0]-P_max[1])*sys->cs_2/(out_x-in_x);
  double dPdx3 = (P_mean[0]-P_mean[1])*sys->cs_2/(out_x-in_x);
  std::cout << "dPdx1 = " << dPdx1 << ", dPdx2 = " << dPdx2 << ", dPdx3 = " << dPdx3 << std::endl;
#endif
  return((sbuff[0]-sbuff[1])*sys->cs_2/(out_x-in_x));
}


int species_has_converged(int t, struct System *sys, struct Field *species, double conv_crit, double *diff)
/**************************************
 **************************************/
{
  int not_converged;

  *diff = -1.0e10;
  if (t%sys->conv_length == 0) {
    not_converged = 1;

    //if (t > sys->steps_to_wait_before_converge) {
    if (t > sys->steps_until_check) {
      check_convergence(&not_converged, sys, species, conv_crit, diff);
    }

    if (not_converged) {
      update_rho_old(sys, species);
      return(0);
    }
    else {
      sys->last_convergence = t;
      return(1);
    }
  }

  return(0);
}


/* void set_init_bulk_conc_from_file(struct Field *species, struct System *sys) */
/* /\************************************** */
/*  **************************************\/ */
/* { */

/*   int cn, cc, ck; */
/*   int g_ind; */
/*   double *g; */
/*   char fn[50]; */

/*   struct Step *step = sys->step; */
/*   struct Node *node = sys->node; */

/*   double *ss_g = (double *) calloc(sys->max_n_tot*sys->n_Q, sizeof(double)); */

/*   if (sys->mpi->my_rank == 0) */
/*     printf("**** READING STEADY-STATE SPECIES FROM FILES out/0000_init_specie_f.bfd - out/%04d_init_specie_f.bfd **** \n\n", sys->mpi->nr_procs-1); */
/*   sprintf(fn, "out/%04d_init_specie_f",sys->mpi->my_rank); */
/*   read_dist_binary(fn, ss_g, sys->max_n_tot*sys->n_Q); */
/*   //  write_rho_vtk_binary_v2(0, ss_g, sys); */

/*   for( cc = 0; cc < species->n_c; ++cc ) {  /\* phase *\/ */
/*     g = &species->f[step->f_phase*cc];      /\* distribution in the current phase *\/ */
/*     g_ind = 0;                              /\* index for nodes distribtuion *\/ */
/*     for( cn = 0; cn < sys->max_n_tot; ++cn ) { /\* node *\/ */
/*       if( node->mode[cn] == FLUID) { */
/* 	if (ss_g[g_ind] > 0.0) { */
/* 	  for( ck = 0; ck < sys->n_Q; ++ck ) { */
/* 	    g[g_ind + ck] = species->rho_inlet[cc]*ss_g[g_ind + ck]; */
/* 	  } */
/* 	} else { */
/* 	  // steady-state dist is 0 because velocity is 0, i.e. isolated pore */
/* 	  for( ck = 0; ck < sys->n_Q; ++ck ) { */
/* 	    g[g_ind + ck] = species->rho_init[cc]*sys->w[ck]; // velocity is 0 */
/* 	  } */
/* 	} */
/*       } */
/*       g_ind   += step->f_node; */
/*     } */
/*   } */
/*   //  write_rho_vtk_binary_v2(1, species->f, sys); */

/*   free(ss_g); */
/* } */

void set_init_bulk_conc(double *rho_init, struct Field *species, struct Field *fluid, struct System *sys)
/**************************************
 **************************************/
{
  int cn, cc, ck;
  int g_ind;
  double *g, *rho, *psi;

  struct Step *step = sys->step;
  struct Node *node = sys->node;

  double *ss_g = (double *)calloc(sys->max_n_tot*sys->n_Q, sizeof(double));
  g = &species->f[0];
  g_ind = 0;
  // copy steady-state species concentration to temporary array ss_g
  for (cn = 0; cn < sys->max_n_tot; ++cn) {
    for (ck = 0; ck < sys->n_Q; ++ck) {
      if (sys->opt->find_steady_state) {
        // g is dist-func from recent steady-state run with rho_init == 1
        ss_g[g_ind + ck] = g[g_ind + ck];
      }
      else {
        // no steady-state run and rho_init != 1
        ss_g[g_ind + ck] = g[g_ind + ck] / species->rho_init[0];
      }
    }
    g_ind += step->f_node;
  }

  for (cc = 0; cc < species->n_c; ++cc) {  // loop phases
    g = &species->f[step->f_phase*cc];      // distribution in the current phase
    rho = &species->rho[step->rho_phase*cc];
#ifdef NEWCONC
    psi = &species->psi[step->rho_phase*cc];
#endif
    g_ind = 0;
    for (cn = 0; cn < sys->max_n_tot; ++cn) { // loop nodes
      for (ck = 0; ck < sys->n_Q; ++ck) {
        g[g_ind + ck] = rho_init[cc] * ss_g[g_ind + ck];
      }

      // update rho and psi
      if (node->mode[cn] == FLUID) {
        rho[cn] = 0;
        for (ck = 0; ck < sys->n_Q; ++ck) {
          rho[cn] += g[g_ind + ck];
        }
#ifdef NEWCONC
#ifdef FLUID_OFF
        psi[cn] = rho[cn];
#else
        psi[cn] = rho[cn] / fluid->rho[cn];
#endif
#endif
      }

      g_ind += step->f_node;
    }
  }

  // set isolated nodes to initial (equilibrium) value
  for (cc = 0; cc < species->n_c; ++cc) {  // loop phases
    g = &species->f[step->f_phase*cc];      // distribution in the current phase
    //for( i = 0; i < node->nr_isolated; ++i ) { 
    for (std::vector<int>::size_type i = 0; i < node->isolated.size(); ++i) { // loop isolated fluid nodes
      //cn = node->list_isolated[i];
      cn = node->isolated[i];
      g_ind = cn*step->f_node;
      for (ck = 0; ck < sys->n_Q; ++ck) {
        g[g_ind + ck] = species->rho_init[cc] * ss_g[g_ind + ck];
      }
      // update rho and psi
      rho[cn] = 0;
      for (ck = 0; ck < sys->n_Q; ++ck) {
        rho[cn] += g[g_ind + ck];
      }
#ifdef NEWCONC
      psi[cn] = rho[cn] / fluid->rho[cn];
#endif
    }
  }

  free(ss_g);

  communicate(species->f, species->n_c, sys);
}


void set_conc_field(struct Field *species, struct Field *fluid, struct System *sys)
/**************************************
 **************************************/
{
  int cn, cc, ck;
  int g_ind;
  double *g, *rho, *psi;

  struct Step *step = sys->step;
  struct Node *node = sys->node;

  double *ss_g = (double *)calloc(sys->max_n_tot*sys->n_Q, sizeof(double));
  g = &species->f[0];
  g_ind = 0;
  // copy steady-state species concentration to temporary array ss_g
  for (cn = 0; cn < sys->max_n_tot; ++cn) {
    for (ck = 0; ck < sys->n_Q; ++ck) {
      ss_g[g_ind + ck] = g[g_ind + ck];
    }
    g_ind += step->f_node;
  }

  for (cc = 0; cc < species->n_c; ++cc) {  // loop phases
    g = &species->f[step->f_phase*cc];      // distribution in the current phase
    rho = &species->rho[step->rho_phase*cc];
    psi = &species->psi[step->rho_phase*cc];
    g_ind = 0;
    for (cn = 0; cn < sys->max_n_tot; ++cn) { // loop nodes
      for (ck = 0; ck < sys->n_Q; ++ck) {
        if (rho[cn] < 0) rho[cn] = 1e-8;
        g[g_ind + ck] = rho[cn] * ss_g[g_ind + ck];
      }

      // update rho and psi
      if (node->mode[cn] == FLUID) {
        rho[cn] = 0;
        for (ck = 0; ck < sys->n_Q; ++ck) {
          rho[cn] += g[g_ind + ck];
        }
        psi[cn] = rho[cn] / fluid->rho[cn];

      }

      g_ind += step->f_node;
    }
  }

  free(ss_g);

  communicate(species->f, species->n_c, sys);

}

void set_species_in_new_fluid_nodes(struct Field *species, struct Field *fluid, struct System *sys)
/**************************************
 **************************************/
{
  int nc, success, new_node;
  int imax = 10, i = 0;
  int num_rest_nodes = 0;
  int remaining_nodes = 1;
  int num_new_nodes = sys->num_new_nodes;
  struct Node *node = sys->node;

  //if (sys->mpi->my_rank==0)
  //  printf("\n**** Setting g and rho in %d new nodes\n", sys->rho_reinit);

  // make copy of new_nodes-list
  for (nc = 0; nc < num_new_nodes; ++nc) {
    sys->new_nodes_cpy[nc] = sys->new_nodes[nc];
  }

  reset_rho_f_in_new_fluid_nodes(species, sys);

  update_rho(node->list_ghost, node->nr_ghost, fluid, species, sys);
  //update_rho(sys->new_nodes, num_new_nodes, species, sys);


  while (remaining_nodes && i < imax) {

    //printf("%d: new_nodes = %d\n", sys->mpi->my_rank, num_new_nodes);
    for (nc = 0; nc < num_new_nodes; ++nc) {
      new_node = sys->new_nodes_cpy[nc];
      success = copy_weighted_f_from_fluid_neighbors(new_node, species, sys);
      if (success == 0) {
        sys->new_nodes_rest[num_rest_nodes++] = new_node;
      }
    }

    communicate(species->f, species->n_c, sys);

    update_rho(node->list_ghost, node->nr_ghost, fluid, species, sys);
    update_rho(sys->new_nodes, num_new_nodes, fluid, species, sys);

    num_new_nodes = 0;

    if (num_rest_nodes) {
      num_new_nodes = num_rest_nodes;
      for (nc = 0; nc < num_new_nodes; ++nc) {
        sys->new_nodes_cpy[nc] = sys->new_nodes_rest[nc];
      }
    }

    MPI_Reduce(&num_rest_nodes, &remaining_nodes, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Bcast(&remaining_nodes, 1, MPI_INT, 0, MPI_COMM_WORLD);

    num_rest_nodes = 0;

#ifdef _DEBUG_
    if (sys->mpi->my_rank == 0 && remaining_nodes > 0)
      printf("*** %d: remaining nodes: %d \n", i, remaining_nodes);
#endif
    i++;
  }

  if (sys->mpi->my_rank == 0 && i >= imax) {
    printf("\n*** WARNING! Max number of iterations in set_species_in_new_fluid_nodes() exceeded! Remaining nodes is %d\n", remaining_nodes);
  }

}

// ---------------------------------------------------------
// ---------------------------------------------------------
void reset_rho_f_in_new_fluid_nodes(struct Field *field, struct System *sys)
{
  int nc, cc, nn, ck;
  struct Step *step = sys->step;
  double *f;

  // reset rho, f in new FLUID nodes
  for (nc = 0; nc < sys->num_new_nodes; ++nc) {
    nn = sys->new_nodes[nc]; // new fluid node
    for (cc = 0; cc < field->n_c; ++cc) {  // phase
      field->rho[step->rho_phase*cc + nn] = 0.0;
      f = &field->f[step->f_phase*cc + step->f_node*nn];
      for (ck = 0; ck < sys->n_Q; ++ck) {
        f[ck] = 0.0;
      }
    }
  }
  communicate(field->f, field->n_c, sys);
}

// ---------------------------------------------------------
// ---------------------------------------------------------
void reset_rho_f_in_new_solid_nodes(struct Field *field, struct System *sys)
{
  int nc, cc, nn, ck;
  struct Step *step = sys->step;
  double *f;

  // reset rho, u, f in new SOLID nodes
  for (nc = 0; nc < sys->num_new_solid_nodes; ++nc) {
    nn = sys->new_solid_nodes[nc]; // new fluid node
    for (cc = 0; cc < field->n_c; ++cc) {  // phase
      field->rho[step->rho_phase*cc + nn] = 0.0;
      f = &field->f[step->f_phase*cc + step->f_node*nn];
      for (ck = 0; ck < sys->n_Q; ++ck) {
        f[ck] = 0.0;
      }
    }
  }
  communicate(field->f, field->n_c, sys);
}


// ---------------------------------------------------------
// copy g from neighboring old fluid nodes weighted by link weight
// ---------------------------------------------------------
int copy_weighted_f_from_fluid_neighbors(int nn, struct Field *species, struct System *sys)
{
  char *node_mode = sys->node->mode;
  struct Step *step = sys->step;
  int ck, ck2, nbn, cc, node, neighbor;
  double *g;
  int success = 0;
  double tot_w = 0.0; // total sum of weight of old fluid neighbors

  // loop over neighbors to find (old, i.e rho>0) fluid nodes
  for (ck = 0; ck < sys->n_Q; ++ck) {  // neighbor direction loop
    nbn = nn - step->n_rel_neig[ck];  // neighbor node
    //if (sys->node->is_ghost[nbn])
    //  continue;
    if (node_mode[nbn] == FLUID && species->rho[nbn] > 0.0) {
      // (old) fluid neighbor found 
      //if (sys->node->is_ghost[nbn])
      //printf("\n*** <%d> copy g from ghost node %d\n", sys->mpi->my_rank, nbn);
      tot_w += sys->w[ck];
      node = nn*step->f_node;
      neighbor = nbn*step->f_node;
      // copy all species and all directions from fluid neighbor
      for (cc = 0; cc < species->n_c; ++cc) {  // specie loop
        g = &species->f[step->f_phase*cc];
        for (ck2 = 0; ck2 < sys->n_Q; ++ck2) { // copy direction loop
          g[node + ck2] += sys->w[ck] * g[neighbor + ck2];
        }
      }
    }
  } // end neighbor direction loop

  if (tot_w > 0.0) {
    // (at least one) fluid neighbor found
    success = 1;
    // set f in new node to weighted mean value
    for (cc = 0; cc < species->n_c; ++cc) {  // specie loop
      g = &species->f[step->f_phase*cc];
      node = nn*step->f_node;
      for (ck = 0; ck < sys->n_Q; ++ck) {
        g[node + ck] /= tot_w;
      }
    }
  }

  return(success);
}


// ---------------------------------------------------------
// 
// ---------------------------------------------------------
//void update_rho(struct Field *species, struct System *sys)
void update_rho(int *node_list, int num_nodes, struct Field *fluid, struct Field *species, struct System *sys)
{
  int nc, nn, cc, ck;
  double *rho, *g;
  struct Step *step = sys->step;
  struct Node *node = sys->node;
#ifdef NEWCONC
  double *psi;
#endif

  for (nc = 0; nc < num_nodes; ++nc) {
    nn = node_list[nc];
    if (node->mode[nn] == FLUID) {
      for (cc = 0; cc < species->n_c; ++cc) {  // phase
        g = &species->f[step->f_phase*cc + step->f_node*nn];
        rho = &species->rho[step->rho_phase*cc];
        rho[nn] = 0.0;
        for (ck = 0; ck < sys->n_Q; ++ck) {
          rho[nn] += g[ck];
        }
#ifdef NEWCONC
        psi = &species->psi[step->rho_phase*cc];
#ifdef FLUID_OFF
        psi[nn] = rho[nn];
#else
        psi[nn] = rho[nn] / fluid->rho[nn];
#endif
#endif      
      } // phase
    }
  }

}



//void set_bulk_conc_after_reinit(struct Field *species, struct Field *fluid, struct System *sys)
///**************************************
// **************************************/
//{
//  int cn, cc, ck, nc, nn, nbn;
//  int g_ind;
//  double *g, *rho, *psi, tot_w;
//
//  struct Step *step = sys->step;
//  struct Node *node = sys->node;
//
//  //  printf("Init species after velocity re-init\n");
//
//  // reset rho in new fluid nodes (should be zero, but just to be sure...)
//  for( nc = 0; nc < sys->num_new_nodes; ++nc ) {
//    nn = sys->new_nodes[nc]; // new fluid node
//    for( cc = 0; cc < species->n_c; ++cc ) {  // phase
//      rho = &species->rho[step->rho_phase*cc];
//      rho[step->rho_node*nn] = 0.0;
//    }
//  }
//
//  // summarize rho in neighboring fluid nodes and set
//  // rho in new nodes to mean rho value
//  for( nc = 0; nc < sys->num_new_nodes; ++nc ) {
//    nn = sys->new_nodes[nc]; // new fluid node
//    tot_w = 0.0;  // total sum of weight of old fluid neighbors
//    for(ck = 0; ck < sys->n_Q; ++ck) {  // direction-loop
//      nbn = nn - step->n_rel_neig[ck];  // neighbor node
//      if (node->mode[nbn] == FLUID && species->rho[step->rho_node*nbn]>0.0) {
//        tot_w += sys->w[ck];
//        for( cc = 0; cc < species->n_c; ++cc ) {  // phase-loop
//          rho = &species->rho[step->rho_phase*cc];
//          rho[step->rho_node*nn] += sys->w[ck]*rho[step->rho_node*nbn];
//        }
//      }
//    } // end ck-loop
//    // set rho in new node to weighted mean value
//    for( cc = 0; cc < species->n_c; ++cc ) {  // phase
//      rho = &species->rho[step->rho_phase*cc];
//      rho[step->rho_node*nn] /= tot_w;
//    }
//  }
//
//  // set g-dist to product of rho_g and f-dist
//  for( cc = 0; cc < species->n_c; ++cc ) {  /* phase */
//    g = &species->f[step->f_phase*cc];      /* distribution in the current phase */
//    rho   = &species->rho[step->rho_phase*cc];
//    psi   = &species->psi[step->rho_phase*cc];
//    g_ind = 0;                              /* index for nodes distribtuion */
//    for( cn = 0; cn < sys->max_n_tot; ++cn ) { /* node */
//      for( ck = 0; ck < sys->n_Q; ++ck ) {
//        g[g_ind + ck] = rho[cn]*fluid->f[g_ind + ck];
//      }
//
//      // update rho and psi
//      if( node->mode[cn] == FLUID) {
//        rho[cn] = 0;
//        for( ck = 0; ck < sys->n_Q; ++ck ) {
//          rho[cn] += g[g_ind + ck];
//        }
//        psi[cn] = rho[cn]/fluid->rho[cn];
//
//      }
//
//      g_ind   += step->f_node;
//    }
//  }
//  //sys->num_new_nodes = 0;
//
//  communicate(species->f, species->n_c, sys);
//
//}



void advance_one_step(struct InitChem *ICS, struct Field *species, struct Field *fluid, struct Field *dif,
  struct System *sys, struct Minerals *minerals, struct Boundary *bndry, struct splayTree **st_lst)
  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
{
  struct Node *node = sys->node;

  // streaming step
  if (sys->fluid) {
    propagate(node, fluid, sys);
  }
  if (sys->chem) {  // AMIN
    propagate(node, species, sys);
    int D_ratio = round(sys->D / sys->D_b); // /* AMIN */ 
    if (sys->t % D_ratio == 0 && sys->t != 0)
      propagate_dif(node, dif, sys); // dif field is updated with the same rate! /* AMIN */   
  }
  // Calculate distribution momenta
#ifdef _USE_COLOR_GRAD_
  if (sys->fluid) {
    calc_rho_vel_f_color_grad(node, fluid, sys);
    // NB: Must communicate  fluid->rho_tot used in correction_rho_vel_color_grad
    communicate_rho_tot(fluid, sys);
  }
//  if (sys->chem) {
//    calc_rho_g_color_grad(node, species, sys);
//    // NB: Must communicate species->rho used in correction_rho_vel_color_grad
//    communicate_rho(species, sys);
//  }

  if (sys->fluid) {
    if (sys->chem)
      correction_rho_vel_color_grad(node, fluid, species, sys);
    else
      correction_rho_vel_color_grad_fluid(node, fluid, sys);
    // NB: Communicate fluid->rho
    communicate_rho(fluid, sys);
  }
#endif

  // Calculate distribution momenta
#ifdef _USE_COLOR_GRAD_
  if (sys->fluid) {
    calc_rho_vel_f_color_grad(node, fluid, sys);
    // NB: Must communicate  fluid->rho_tot used in correction_rho_vel_color_grad
    communicate_rho_tot(fluid, sys);
  }
//  if (sys->chem) {
//    calc_rho_g_color_grad(node, species, sys);
//    // NB: Must communicate species->rho used in correction_rho_vel_color_grad
//    communicate_rho(species, sys);
//  }

  if (sys->fluid) {
    if (sys->chem)
      correction_rho_vel_color_grad(node, fluid, species, sys);
    else
      correction_rho_vel_color_grad_fluid(node, fluid, sys);
    // NB: Communicate fluid->rho
    communicate_rho(fluid, sys);
  }
#endif

  // collision step
  if (sys->fluid) {
#ifdef NEWCONC
#ifdef _USE_COLOR_GRAD_
//    calc_rho_vel_f_color_grad(node, fluid, sys);
    collision_f_color_grad(node, fluid, sys);
#else
    collision_f_newconc(node, fluid, sys);
#endif
#else
    collision_f(node, fluid, sys);
#endif
  }

  if (sys->chem) {
#ifdef NEWCONC
// #ifdef _USE_COLOR_GRAD_
//     collision_g_color_grad(node, species, fluid, sys);
// #else
	  // AMIN
    //collision_g_newconc(species, fluid, sys);
    //collision_g_newconc(dif, fluid, sys);
    collision_g_newconc_dif(dif, species, fluid, sys, ICS);
#else
    collision_g(species, fluid, sys);
    collision_g(dif, fluid, sys);
#endif
  }

}

//------------------------------------------
//
//------------------------------------------
void apply_boundary_conditions(InitChem *ICS, Field *species, Field *dif, Field *fluid, System *sys, Minerals *minerals,
    Boundary *bndry, splayTree **st_lst)
{
  // BOUNDARY CONDITIONS (NB! the order is important)
  if( sys->fluid ) {
#ifdef _FLUID_BDRY_PRESS_
    bc_3D_copy_gradient(sys->node->nr_in , sys->node->list_in , fluid, sys, sys->step->neg_x_dir, sys->step->pos_x_dir);
    bc_3D_copy_gradient(sys->node->nr_out, sys->node->list_out, fluid, sys, sys->step->pos_x_dir, sys->step->neg_x_dir);
#endif
    bc_run_f(sys->node, fluid, sys);
  }

  if( sys->chem ) {
#ifdef _PERIODIC_RUN_
    periodic_bc(sys->node, species, sys);
    periodic_bc(sys->node, dif, sys);
    //wall_bc_bounce_back(species, sys);
    wall_bc_mid_grid_gca_ej_g_v3(ICS, species, fluid, sys, minerals, bndry, st_lst);
    wall_bc_bounce_back(sys->node, dif, sys);
#else
    bc_run_g(ICS, species, fluid, sys, minerals, bndry, st_lst);
    wall_bc_bounce_back(sys->node, dif, sys);
#endif
    communicate(species->f, species->n_c, sys);
    communicate(species->f_tmp, species->n_c, sys);
    
    communicate(dif->f, dif->n_c, sys);
    communicate(dif->f_tmp, dif->n_c, sys);
  }
}

//------------------------------------------
//
//------------------------------------------
void print_status_and_write_output(struct InitChem &ICS_pH, Output &output, struct System *sys, time_t t_begin, int no_comb, struct InitChem *ICS,
  struct Field *species, struct Field *fluid, struct Minerals *minerals, struct splayTree **st_lst, int write_velocity)
  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
{
  static int print_header = 1;
  if (print_header) {
    print_status_header(sys);
    print_header = 0;
  }

  if ((sys->t + 1) % 100 == 0) {
    print_status(sys, ICS, t_begin, st_lst, no_comb);
  }

  if (sys->time + sys->dt > sys->nwrite*sys->write_int) {
    print_status(sys, ICS, t_begin, st_lst, no_comb);
    if (sys->chem) {
        update_pH_buffer(ICS_pH, &(output.chem), species, sys);
      // write species concentrations and solid fraction
      output.chem.write(sys);
      output.mineral.write(sys);
      output.dif.write(sys);
    }
    
    // write density and velocity
    if (sys->fluid && write_velocity) {
      output.fluid.write(sys);
    }
    sys->nwrite++;
    if (sys->mpi->my_rank == sys->mpi->nr_procs - 1) {
      std::cout << std::endl << &std::fflush;// ("\n");
    }
  }
#ifdef _DEBUG_SPLAY_
  print_splaytree_size(sys, no_comb, st_lst);
#endif
}

//--------------------------------------------------------
// Set parameters from inp.dat in LB units (i.e. inp.dat has no header # REAL UNITS )
//--------------------------------------------------------
void set_parameters_LB(struct System *sys, struct Field *species, struct Field *fluid, struct InitChem *ICS)
{
  double lb;
  double C0 = 1e3; // conversion from moles/L to moles/m^3
  int i;
  struct Options *opt = sys->opt;
  double dx = sys->dx;
  double dx2 = dx*dx;
  int nfluid_tot = 0, nfluid_reac = 0;

  //sys->D = sys->opt->D;
  double tau = 1.0;
  if (sys->chem) {
    tau = species->tau[0];  
  } 
  sys->Dlb = (1. / 3.)*(tau - 0.5);
  sys->dt = sys->Dlb*sys->dx*sys->dx / sys->D; // timestep

  // perform fluid run tau=1
  fluid->dt = sys->cs_2*(1.0 - 0.5)*sys->dx*sys->dx / sys->D;
  if (sys->chem) {
    species->dt = sys->dt;
  }
  
  // convert reporting timestep from steps to seconds
  sys->write_int *= sys->dt;
  sys->write_int -= 1e-20;

  // set effluent and perm write interval
  opt->effluent_interval *= sys->dt;
  opt->effluent_interval -= 1e-14;
  if (opt->perm_interval > 0.0) {
    opt->perm_interval *= sys->dt;
    opt->perm_interval -= 1e-14;
  }

  // calculate rate constants in LB units
  // rate constants are given in mol/m^2/s in inp.dat
  if (sys->chem) {
    lb = sys->dt / (C0*sys->dx);
    for (i = 0; i < ICS->size_sup_min; ++i) {
      ICS->rate[i][0] *= lb;
      ICS->rate[i][1] *= lb;
    }
  }
  
  // calculate pore volume and porosity (nfluid_flow = nfluid_reac)
  get_pore_volume(sys->node, &nfluid_tot, &nfluid_reac, sys);
  sys->porosity = (double)nfluid_reac / (double)sys->V_lb;
  sys->pore_volume_tot = (double)nfluid_tot;
  sys->pore_volume_reac = (double)nfluid_reac;

  // kinematic viscosity
  //sys->nu_lb = sys->cs_2*(fluid->tau[0] - 0.5);
  sys->nu = sys->nu_lb*dx2 / fluid->dt;

  
  sys->steps_to_wait_before_converge = 0;
  
//  if (sys->mpi->my_rank == 0) {
//    printf("\n\n");
//    printf("***************************************************************************\n");
//    printf("*  Input given in LB units, dt = Dlb*dx2/D                                *\n");
//    //printf("*  dt = %.3e s, dx = %.3e m, D = %.2e m/s2, Dlb = %.2e  *\n", sys->dt, sys->dx, sys->D, sys->Dlb);
//    printf("***************************************************************************\n\n");
//  }

}

//-------------------------------------------
// set parameters from inp.dat in REAL units (i.e. inp.dat has # REAL UNITS as header)
//-------------------------------------------
void set_parameters(struct System *sys, struct Field *species, struct Field *dif, struct Field *fluid, struct InitChem *ICS)
{
  int cc, i;
  double tPV_real = 0, tau;
  double C0 = 1e3; // conversion from moles/L to moles/m^3
  struct Options *opt = sys->opt;
  int nfluid_tot, nfluid_reac;
  double dx = sys->dx;
  double dx2 = dx*dx;
  double dx3 = dx2*dx;
  double dt_min, dt_tau1, dt_tau_inp, dt_b;
#ifndef _TAU_FROM_INP_
  int sum_inlet;
#endif

  // set limit on dt by lowest tau-value
  dt_min = sys->cs_2*(sys->opt->tau_min - 0.5)*dx2 / sys->D;
  dt_tau1 = sys->cs_2*(1.0 - 0.5)*dx2 / sys->D;
  dt_b = sys->cs_2*(1.0 - 0.5)*dx2 / sys->D_b; /* AMIN */
  if (sys->chem) {
    dt_tau_inp = sys->cs_2*(species->tau[0] - 0.5)*dx2 / sys->D;
  } else {
    dt_tau_inp = dt_tau1;
  }
    
  // perform fluid run at tau=1
  fluid->dt = sys->cs_2*(fluid->tau[0] - 0.5)*sys->dx*sys->dx / sys->D;

  // calculate pore volume and porosity (nfluid_flow = nfluid_reac)
  get_pore_volume(sys->node, &nfluid_tot, &nfluid_reac, sys);
  sys->pore_volume_tot = (double)nfluid_tot;
  sys->pore_volume_reac = (double)nfluid_reac;
  //sys->porosity = sys->pore_volume_tot / (double)sys->V_lb;
  sys->porosity = sys->pore_volume_reac / (double)sys->V_lb;  // 
  if (opt->straight_tube)
    sys->porosity = 1.0;

#ifdef _TAU_FROM_INP_
  sys->dt = dt_tau_inp;
  sys->flux = fluid->flux*dx3 / sys->dt;
  if (fluid->flux > 0)
    tPV_real = nfluid_reac*sys->dt / fluid->flux;  //tPV_real = nfluid_flow*dx3/sys->flux;
#else
  if ((sys->flux>0) && (fluid->flux>0)) {
    // sys->flux is given as mean fluid velocity (Darcy velocity) in inp.dat in m/s
    if (opt->straight_tube) {
      MPI_Reduce(&(sys->node->nr_in), &sum_inlet, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
      MPI_Bcast(&sum_inlet, 1, MPI_INT, 0, MPI_COMM_WORLD);
      sys->flux *= dx2*sum_inlet;
    } else {
      sys->flux *= dx2*sys->W_lb*sys->H_lb;
    }
    
    // fluid->flux is in LB units and calculated in reinit_velocity
    // The time step is given by the ratio of the two fluxes
    sys->dt = fluid->flux*dx3 / sys->flux; // fluid->flux is calculated in initialize_velocity and given in LB units
    
    // check if the timestep yields a too low tau-value
    if (sys->dt < dt_min) {
      sys->dt = dt_min;
      sys->flux = fluid->flux*dx3 / sys->dt;
    }
    tPV_real = nfluid_reac*dx3 / sys->flux;
  } else {
    // flux is zero. Set dt using tau=1 to a fixed value
    sys->dt = dt_tau1;
  }
#endif

  if (opt->boost_rate) {
    sys->dt *= opt->boost_rate;
  }

  if (sys->chem)
    species->dt = sys->dt;
  if (opt->closed_cell == 0)
    fluid->dt = sys->dt;


  // set effluent and perm write interval
  opt->effluent_interval *= sys->dt;
  opt->effluent_interval -= 1e-14;
  if (opt->perm_interval > 0.0) {
    opt->perm_interval *= sys->dt;
    opt->perm_interval -= 1e-14;
  }

  // used in perm calculation
  //sys->dt_ratio = dt_tau1 / sys->dt;

  // convert from m3/s to LB units
  //Q2Qlb = sys->dt/dx3;

  // calculate tau from diffusion relation and set species->tau
  sys->Dlb = sys->D*sys->dt / dx2;
  sys->Dlb_b = sys->D_b*sys->dt_b / dx2; // for biomass field/* AMIN */
  tau = 0.5 + sys->Dlb*(1. / sys->cs_2);

  if (sys->chem) {
    for (cc = 0; cc < species->n_c; ++cc)
      species->tau[cc] = tau;
    sys->dt_b = dt_b; /* AMIN */
    (*ICS).BiP->mu_max *= sys->dt / 3600 ;
    (*ICS).BiP->K_i *=C0;
    // calculate rate constants in LB units
    // rate constants k are given in mol/m^2/s in inp.dat. 
    // The rate constant in the code, k' = S*k/V, where [S] = m2, and [V] = m3
    double lb = sys->dt / (C0*dx);
    for (i = 0; i < ICS->size_sup_min; ++i) {
      ICS->rate[i][0] *= lb;
      ICS->rate[i][1] *= lb;
    }
  }

  // kinematic viscosity
  sys->nu = sys->nu_lb*dx2 / fluid->dt;

  //   sys->steps_to_wait_before_converge = (int)(tPV_real*1.5/sys->dt);
  sys->steps_to_wait_before_converge = 0;

  if (sys->opt->restart_run)
    sys->steps_to_wait_before_converge = 0;


  // //  calc_effluent_and_perm(sys, sys->node, species, fluid, 1);
  // if (sys->fluid && )
  //   umax = get_umax(sys->node, sys, fluid);

}

// -------------------------------------------------------
//
// -------------------------------------------------------
void print_run_parameters(struct System *sys, struct Field *species, struct Field *fluid, struct InitChem *ICS, struct Minerals *minerals)
{
  int p = 0, i;
  FILE *data_table;
  int buffer_size = 2700;
  double umax = 0.0;
  double dx = sys->dx;
  double dx2 = dx*dx;
  double dx3 = dx2*dx;
  double sumW = 0;
  double Q2Qlb = sys->dt / dx3;
  double L2m3 = 1e-3; // convert from liter to m3
  //double mD_to_lb = 1. / (sys->dx*sys->dx*1.013250e15); // 1 m2 = 1.013249966e12 Darcy, 1 Darcy = 1000 mD
  double D_to_lb = 1. / (sys->dx*sys->dx*1.013250e12); // 1 m2 = 1.013249966e12 Darcy, 1 Darcy = 1000 mD
  double C0 = 1e3; // conversion from moles/L to moles/m^3
  double lb = sys->dt / (C0*dx);

  if (sys->opt->skip_init_velocity_run == 0)
    calc_effluent_and_perm(sys, sys->node, species, fluid, 1);

  if (sys->fluid)
    umax = get_umax(sys->node, sys, fluid);

  if (sys->chem)
    get_reactive_surface_area(sys->node, &sumW, sys, minerals);
  else
    get_reactive_surface_area(sys->node, &sumW, sys, NULL);

  if (sys->mpi->my_rank == 0) {
    //char buffer[buffer_size];
    char *buffer = (char *)malloc(buffer_size * sizeof(char));
    p = sprintf(buffer, "\n   -----------------------------------------------------\n");
    p += sprintf(buffer + p, "   |   Parameter   |  Physical units |     LB units    |\n");
    p += sprintf(buffer + p, "   -----------------------------------------------------\n");
    p += sprintf(buffer + p, "   | Dimensions    |      -----      | %3d x %3d x %3d |\n", sys->MAX_N[0], sys->MAX_N[1], sys->MAX_N[2]);
    p += sprintf(buffer + p, "   | L [m]         |  %12.5e   |     %4d        |\n", dx*sys->L_lb, sys->L_lb);
    p += sprintf(buffer + p, "   | W [m]         |  %12.5e   |     %4d        |\n", dx*sys->W_lb, sys->W_lb);
    p += sprintf(buffer + p, "   | H [m]         |  %12.5e   |     %4d        |\n", dx*sys->H_lb, sys->H_lb);
    p += sprintf(buffer + p, "   | dx [m]        |  %12.5e   |        1        |\n", dx);
    p += sprintf(buffer + p, "   | dt [s]        |  %12.5e   |        1        |", sys->dt);
    if (sys->real_units_in_input==0) {
      p += sprintf(buffer + p, "  dt = D_lb*dx^2/D");
    }
    p += sprintf(buffer + p, "\n");
    p += sprintf(buffer + p, "   | dt_fluid  [s] |  %12.2e   |        1        |\n", fluid->dt);
    p += sprintf(buffer + p, "   | D [m2/s]      |  %12.5e   |   %10.3e    |\n", sys->D, sys->Dlb);
    for (int c=0; c<fluid->n_c; ++c) {
      p += sprintf(buffer + p, "   | tau_fluid_%d   |       ---       |       %4.2f      |\n", c, fluid->tau[c]); 
    }
    for (int c=0; c<fluid->n_c; ++c) {
      p += sprintf(buffer + p, "   | nu_%d [m2/s]   |  %12.5e   |   %10.3e    |\n", c, fluid->nu_lb[c]*dx2/fluid->dt, fluid->nu_lb[c]);
    }
    //for (int c=0; c<fluid->n_c; ++c) {
    //  p += sprintf(buffer + p, "   | Ca_%d          |  %12.5e   |   %10.3e    |\n", c, fluid->nu_lb[c]*dx2/fluid->dt, fluid->nu_lb[c]);
    //}
    if (sys->chem)
      p += sprintf(buffer + p, "   | T [deg C]     |     %6.2f      |       ---       |\n", ICS->Temp - 273.15);
    if (sys->chem) {
      p += sprintf(buffer + p, "   | tau_chem [-]  |       ---       |       %4.2f      |", species->tau[0]);
      if ((species->tau[0] > 1.0) || (species->tau[0] < 1.0)) {
	p += sprintf(buffer + p, "*");
      }
      p += sprintf(buffer + p, "\n");
    }
    //p += sprintf(buffer + p, "   | fluid phases  |       ---       |        %d        |\n", fluid->n_c);
    if (sys->chem)
      p += sprintf(buffer + p, "   | chem phases   |       ---       |        %d        |\n", species->n_c);
    p += sprintf(buffer + p, "   | Q [L/s]       |   %10.3e    |   %10.3e    |\n", fluid->flux / Q2Qlb / L2m3, fluid->flux);
    //p += sprintf(buffer+p, "   | Q [L/s]       |    %9.3e    |   %10.4e    |%s\n", sys->flux/L2m3, sys->flux*Q2Qlb, "  from inp.dat");
    p += sprintf(buffer + p, "   | PV_tot [m3]   |   %10.3e    |   %10d    |\n", sys->pore_volume_tot*dx3, (int)sys->pore_volume_tot);
    p += sprintf(buffer + p, "   | PV_reac [m3]  |   %10.3e    |   %10d    |\n", sys->pore_volume_reac*dx3, (int)sys->pore_volume_reac);
    //p += sprintf(buffer+p, "   | tPV_tot [s]   |    %9.2e    |   %10.2e    |\n" , tPV_real, tPV_real/sys->dt);
    p += sprintf(buffer + p, "   | SA [m2]       |   %10.3e    |   %10.3e    |\n", sumW*dx2, sumW);
    p += sprintf(buffer + p, "   | porosity      |   %10.3e    |   %10.3e    |\n", sys->porosity, sys->porosity);
    p += sprintf(buffer + p, "   | perm [D]      |   %10.3e    |   %10.3e    |\n", fluid->perm, fluid->perm*D_to_lb);
    p += sprintf(buffer + p, "   | ux_mean [m/s] |   %10.3e    |       ---       |\n", fluid->mean_velx);
    p += sprintf(buffer + p, "   | ux_max  [m/s] |   %10.3e    |   %10.3e    |\n", dx*umax / sys->dt, umax);
    //#ifdef _SOURCE_
    //    p += sprintf(buffer+p, "   |uD_source [m/s]|    %9.2e    |       ---       |\n", sys->opt->u_Darcy);
    //#endif
    if (sys->chem) {
      p += sprintf(buffer + p, "   -----------------------------------------------------\n");
      p += sprintf(buffer + p, "   -----------------------------------------------------\n");
      for (i = 0; i < ICS->size_sup_min; ++i) {
	p += sprintf(buffer + p, "   | %-13s |                 |                 |\n", ICS->sup_min_name[i]);
	p += sprintf(buffer + p, "   | k1 [mol/m2/s] |    %9.2e    |    %9.2e    |\n", ICS->rate[i][0] / lb, ICS->rate[i][0]);
	p += sprintf(buffer + p, "   | k2 [mol/m2/s] |    %9.2e    |    %9.2e    |\n", ICS->rate[i][1] / lb, ICS->rate[i][1]);
      }
    }
    p += sprintf(buffer + p, "   -----------------------------------------------------\n\n");
    
    p += sprintf(buffer + p, "   Number of processes: %d x %d x %d\n\n\n", sys->mpi->np[0], sys->mpi->np[1], sys->mpi->np[2]);

    if (p > buffer_size) {
      printf("\n!!ERROR!! Buffer overflow in set_parameters(): %d > %d. Aborting....\n\n", p, buffer_size);
      MPI_Finalize();
      exit(0);
    }

    printf("%s", buffer);
    fflush(stdout);
    std::string filename(sys->files["out"] + "/data_table.txt");
    data_table = fopen(filename.c_str(), "w");
    fprintf(data_table, "%s", buffer);
    fclose(data_table);
    free(buffer);
  }

}


// -------------------------------------------------------
// find number of fluid nodes in all processes
// 
// -------------------------------------------------------
void get_pore_volume(struct Node *node, int *nfluid_tot, int *nfluid_reac, struct System *sys)
{
  // find pore volume by summing fluid nodes (in reactive part of system)
  *nfluid_tot = *nfluid_reac = 0;
  for (int cn = 0; cn < sys->max_n_tot; ++cn) { /* node */
    if ((node->mode[cn] == FLUID) && !(node->is_ghost_or_periodic[cn])) {
      (*nfluid_tot)++;
      // check if fluid node is in the reactive part of the system
      if (node_is_reactive(cn, sys)) {
        (*nfluid_reac)++;
      }
    }
  }
  // sum nfluid_tot and nfluid_reac between all processes
  MPI_Allreduce(MPI_IN_PLACE, nfluid_tot, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, nfluid_reac, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  //  isbuf = *nfluid_tot;
  //  MPI_Reduce(&isbuf, nfluid_tot, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  //  MPI_Bcast(nfluid_tot, 1, MPI_INT, 0, MPI_COMM_WORLD);
  //  isbuf = *nfluid_reac;
  //  MPI_Reduce(&isbuf, nfluid_reac, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  //  MPI_Bcast(nfluid_reac, 1, MPI_INT, 0, MPI_COMM_WORLD);
}

// -------------------------------------------------------
// Get saturation, only different from pore volume for two-phase flows
// 
// -------------------------------------------------------
void get_saturation(struct Node *node, struct Field *fluid, struct System *sys)
{
  struct Step *step = sys->step;
  struct Mpi *mpi = sys->mpi;

  // reset
  for (int cc = 0; cc < fluid->n_c; ++cc) {
    fluid->saturation[cc] = 0.0;
  }
  
  //double rho_tot = 0.0;
  // loop over fluid nodes
  for (int cn = 0; cn < sys->max_n_tot; ++cn) {
    if (node->mode[cn] == FLUID && !node->is_ghost_or_periodic[cn]) {
      for (int cc = 0; cc < fluid->n_c; ++cc) {
	//int ind = step->rho_phase*cc + cn;
        fluid->saturation[cc] += fluid->rho[step->rho_phase*cc + cn];
        //fluid->saturated_nodes[cc] += fluid->rho[step->rho_phase*cc + cn];
	//rho_tot += fluid->rho[ind];
      }
    }
  }

  // sum, communicate and broadcast using allreduce
  MPI_Allreduce(MPI_IN_PLACE, &fluid->saturation[0], fluid->n_c, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  //MPI_Allreduce(MPI_IN_PLACE, &rho_tot             , 1         , MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  double sat_sum = 0;
  for (auto &sat : fluid->saturation)
    sat_sum += sat;

  for (auto &sat : fluid->saturation)
    if (sat_sum>0) {
      sat /= sat_sum;
    } else {
      sat = 0.0;
    }

}

// -------------------------------------------------------
// get total solid mass of system
// 
// -------------------------------------------------------
void get_solid_massfraction(struct Minerals *minerals, struct System *sys)
{
  int nn, min, X, Y, Z;
  struct Node *node = sys->node;
  int Xlim[2] = { sys->opt->in_nodes + 1, sys->MAX_N[0] - 2 - sys->opt->out_nodes };
  int Ylim[2] = { sys->opt->wall_nodes + 1, sys->MAX_N[1] - 2 - sys->opt->wall_nodes };
  int Zlim[2] = { sys->opt->wall_nodes + 1, sys->MAX_N[2] - 2 - sys->opt->wall_nodes };

#ifdef _2D_RUN_
  Zlim[0] = Zlim[1] = 0;
#endif

  // reset
  for (min = 0; min < minerals->n_tot; ++min) {
    minerals->sum_sfrac[min] = 0.0;
  }

  for (nn = 0; nn < sys->max_n_tot; ++nn) {
    if (node->is_ghost_or_periodic[nn])
      continue;

    X = get_global_xcoord(nn, sys);
    Y = get_global_ycoord(nn, sys);
    Z = get_global_zcoord(nn, sys);

    if (X<Xlim[0] || Y<Ylim[0] || Z<Zlim[0] || X>Xlim[1] || Y>Ylim[1] || Z>Zlim[1])
      continue;

    for (min = 0; min < minerals->n_tot; ++min) {
      minerals->sum_sfrac[min] += minerals->list[min].sfrac[nn];
    }
  }

  for (min = 0; min < minerals->n_tot; ++min) {
    minerals->buffer[min] = minerals->sum_sfrac[min];
  }
  MPI_Reduce(minerals->buffer, minerals->sum_sfrac, minerals->n_tot, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Bcast(minerals->sum_sfrac, minerals->n_tot, MPI_DOUBLE, 0, MPI_COMM_WORLD);
}

/* ---------------------------------------------------------------------------------------------- */
double get_rhomax(struct Node *node, struct System *sys, struct Field *fluid)
/* ---------------------------------------------------------------------------------------------- */
/* Calculate max rho      */
/*                        */
{
  int nn;
  double rhomax = -9.e9, rhomax_all, rho;

  for (nn = 0; nn < sys->max_n_tot; ++nn) {
    if (node->mode[nn] == FLUID) {
      if ((rho = fluid->rho[sys->step->rho_node*nn]) > rhomax)
        rhomax = rho;
    }
  }
#ifdef _MPI_
  // communicate max rho
  MPI_Reduce(&rhomax, &rhomax_all, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  if (sys->mpi->my_rank == 0) rhomax = rhomax_all;
  // root-process broadcasts to all processes
  MPI_Bcast(&rhomax, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif

  return rhomax;
}


/* ---------------------------------------------------------------------------------------------- */
double get_umax(struct Node *node, struct System *sys, struct Field *fluid)
/* ---------------------------------------------------------------------------------------------- */
/* Calculate max velocity      */
/*                             */
{
  int nn;
  double umax = -9.e9, umax_all, *u, u2;

  for (nn = 0; nn < sys->max_n_tot; ++nn) {
    if (node->mode[nn] == FLUID) {
      //u = &fluid->u[sys->step->u_node*nn];
      u = &fluid->u_tot[sys->step->u_node*nn];
      u2 = u[0]*u[0] + u[1]*u[1] + u[2]*u[2];
      if (u2 > umax)
        umax = u2;
    }
  }
  umax = sqrt(umax);
  // communicate max velocity
  MPI_Allreduce(MPI_IN_PLACE, &umax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  //if (sys->mpi->my_rank == 0) umax = umax_all;
  // root-process broadcasts to all processes
  //MPI_Bcast(&umax, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  
  return umax;
}


/* ---------------------------------------------------------------------------------------------- */
double get_flux(struct Node *node, struct System *sys, struct Field *fluid)
/* ---------------------------------------------------------------------------------------------- */
/* Calculate effluent flux */
{
  double flux = 0., *u, area;
  int cc, j, k, nn;
  int *n = sys->max_n;
  int A = 0, F = 1;
  //double u2_min = sys->opt->u_min_for_flow*sys->opt->u_min_for_flow;

  struct Mpi *mpi = sys->mpi;

  if (mpi->eff_rank != MPI_UNDEFINED) {

    // zero out the buffer
    for (cc = 0; cc < 2; ++cc)
      mpi->sbuf_eff[cc] = 0.0;

    // loop over effluent nodes
    for (j = 1; j < n[1] - 1; j++) {
      for (k = 1; k < n[2] - 1; k++) {
        nn = (n[0] - 1 - sys->opt->eff_offset) + j*n[0] + k*n[0] * n[1];
        if (node->mode[nn]==FLUID || node->mode[nn]==INLET || node->mode[nn]==OUTLET) {
          //u = &fluid->u[sys->step->u_node*nn];
          u = &fluid->u_tot[sys->step->u_node*nn];
          mpi->sbuf_eff[A] += 1.0;
          mpi->sbuf_eff[F] += u[0];
        }
      }
    }
    // communicate
    MPI_Reduce(mpi->sbuf_eff, mpi->rbuf_eff, 2, MPI_DOUBLE, MPI_SUM, 0, mpi->EFF_COMM);

    // set effluent flux
    if (mpi->eff_rank == 0) {
      flux = mpi->rbuf_eff[F];
      area = mpi->rbuf_eff[A];
    }

  } // end if (mpi->eff_rank

  // broadcast effluent flux
  MPI_Bcast(&flux, 1, MPI_DOUBLE, mpi->eff_root_rank, MPI_COMM_WORLD);
  MPI_Bcast(&area, 1, MPI_DOUBLE, mpi->eff_root_rank, MPI_COMM_WORLD);

  //fluid->flux = flux;
  fluid->area = area;

  return flux;
}

// --------------------------------------------------------------------------
void get_max_min_rho(Node *node, System *sys, Field *fluid)
// --------------------------------------------------------------------------
{
  fluid->rho_max = 0;
  fluid->rho_min = 1e2;
  for (int cn = 0; cn < sys->max_n_tot; ++cn) { /* nodes */
    if ( (node->mode[cn] == FLUID) && !(node->is_ghost_or_periodic[cn])) {
      double sum = 0;
      for (int cc = 0; cc < fluid->n_c; ++cc) { /* Phases */
        sum += fluid->rho[sys->step->rho_phase*cc + cn];
      }
      if (sum > fluid->rho_max)
        fluid->rho_max = sum;
      if (sum < fluid->rho_min)
        fluid->rho_min = sum;
    }
  }

}
// --------------------------------------------------------------------------
double get_mean_velocity(Node *node, System *sys, Field *fluid)
// --------------------------------------------------------------------------
{
  Mpi *mpi = sys->mpi;
  Step *step = sys->step;

  // zero out
  for (auto & u_mean : fluid->phase_u_mean)
    u_mean = 0;
  
  // calculate sum
  std::vector<double> rho_sum(fluid->n_c, 0.0);
  for (int cn = 0; cn < sys->max_n_tot; ++cn) { /* nodes */
    if ( (node->mode[cn] == FLUID) && !(node->is_ghost_or_periodic[cn])) {
      // sum rho at node cn
      double A = 0.0;
      for (int cc = 0; cc < fluid->n_c; ++cc) { /* Phases */
        A += fluid->rho[step->rho_phase*cc + cn];
      }
      double *u = &fluid->u_tot[step->u_node*cn];
      for (int cc = 0; cc < fluid->n_c; ++cc) { /* Phases */
        double rho = fluid->rho[step->rho_phase*cc + cn]/A;
        rho_sum[cc] += rho;
        for (int d=0; d < sys->n_D; ++d) {
          fluid->phase_u_mean[cc*sys->n_D + d] += u[d]*rho;
        } 
      }
    }
  }

  // communicate (MPI_IN_PLACE replaces send and receive buffers)
  MPI_Allreduce(MPI_IN_PLACE, &fluid->phase_u_mean[0], fluid->phase_u_mean.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &rho_sum[0]      , rho_sum.size()      , MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  double rho_tot = 0.0;
  for (auto & rho : rho_sum )
    rho_tot += rho;
  
  for (auto & u_mean : fluid->u_mean)
    u_mean = 0;

  // calculate mean velocity and perm
  for (int cc=0; cc<fluid->n_c; ++cc) {
    double saturation = rho_sum[cc]/rho_tot;
    for (int d=0; d<sys->n_D; ++d) {
      int ind = cc*sys->n_D + d;
      if (rho_sum[cc]>0) {
        fluid->phase_u_mean[ind] /= rho_sum[cc];
      } else {
        fluid->phase_u_mean[ind] = 0.0;
      }
      fluid->u_mean[d] += saturation*fluid->phase_u_mean[ind];
    }
  }
  
  return (fluid->u_mean[0]);
}


// -------------------------------------------------------
// get total amount of dissolved species (in moles)
// 
// -------------------------------------------------------
void get_fluid_conc(struct Field *species, struct System *sys)
{
  int nn, cc, X, Y, Z;
  struct Node *node = sys->node;
  struct Step *step = sys->step;
  struct Mpi *mpi = sys->mpi;
  int Xlim[2] = { sys->opt->in_nodes + 1, sys->MAX_N[0] - 2 - sys->opt->out_nodes };
  int Ylim[2] = { sys->opt->wall_nodes + 1, sys->MAX_N[1] - 2 - sys->opt->wall_nodes };
  int Zlim[2] = { sys->opt->wall_nodes + 1, sys->MAX_N[2] - 2 - sys->opt->wall_nodes };
  double fluid_nodes = 0.0;

#ifdef NEWCONC
  double *conc = species->psi;
#else
  double *conc = species->rho;
#endif
#ifdef _2D_RUN_
  Zlim[0] = Zlim[1] = 0;
#endif

  // reset
  for (cc = 0; cc < species->n_c; ++cc) {
    species->sum_conc[cc] = 0.0;
  }

  for (nn = 0; nn < sys->max_n_tot; ++nn) {
    if (node->is_ghost_or_periodic[nn] || node->mode[nn] != FLUID)
      continue;

    X = get_global_xcoord(nn, sys);
    Y = get_global_ycoord(nn, sys);
    Z = get_global_zcoord(nn, sys);

    if (X<Xlim[0] || Y<Ylim[0] || Z<Zlim[0] || X>Xlim[1] || Y>Ylim[1] || Z>Zlim[1])
      continue;

    fluid_nodes += 1.0;
    for (cc = 0; cc < species->n_c; ++cc) {
      species->sum_conc[cc] += conc[step->rho_phase*cc + nn];
    }
  }

  mpi->sbuf_eff[0] = fluid_nodes;
  for (cc = 0; cc < species->n_c; ++cc) {
    mpi->sbuf_eff[1 + cc] = species->sum_conc[cc];
  }

  //MPI_Reduce(mpi->sbuf_eff, species->total_mole, species->n_c, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(mpi->sbuf_eff, mpi->rbuf_eff, 1 + species->n_c, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Bcast(mpi->rbuf_eff, 1 + species->n_c, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  //if (mpi->my_rank == 0) {
  fluid_nodes = mpi->rbuf_eff[0];
  for (cc = 0; cc < species->n_c; ++cc) {
    species->sum_conc[cc] = mpi->rbuf_eff[1 + cc];
    species->mean_conc[cc] = species->sum_conc[cc] / fluid_nodes;
  }
  //}
  //MPI_Bcast(species->sum_conc , species->n_c, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  //MPI_Bcast(species->mean_conc, species->n_c, MPI_DOUBLE, 0, MPI_COMM_WORLD);
}


// -------------------------------------------------------
//
// -------------------------------------------------------
void get_slice_surface(struct System *sys, int xpos, double *surf, int *nlinks)
{
  int gn, ln;
  struct Links *links = sys->node->links;
  struct Links::Group *group;
  struct Links::Link *link;
  double dsbuf;
  int isbuf;
  int *n = sys->max_n;
  (*surf) = 0.0;
  (*nlinks) = 0;

  for (gn = 0; gn < links->ngroups; ++gn) { // loop link groups
    group = &(links->group[gn]);
    if (group->wall%n[0] == xpos) {
      for (ln = group->first; ln < group->last; ++ln) { // loop link-number in group
        link = &(links->list[ln]);
        if (link->fluid%n[0] == xpos) {
          (*surf) += sys->w[link->dir];
          (*nlinks)++;
        }
      }
    }
  }
  (*surf) *= 6.0;
  dsbuf = *surf;
  MPI_Reduce(&dsbuf, surf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Bcast(surf, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  isbuf = *nlinks;
  MPI_Reduce(&isbuf, nlinks, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Bcast(nlinks, 1, MPI_INT, 0, MPI_COMM_WORLD);
}

// -------------------------------------------------------
// calculate reactive surface area of geometry
// 
// -------------------------------------------------------
void get_reactive_surface_area(struct Node *node, double *sumW, struct System *sys, struct Minerals *minerals)
{
  struct Links::Link *link;
  struct Links *links = node->links;
  struct Mpi *mpi = sys->mpi;
  int nlink = 0;
  double sfrac_sum;
  int num_comm = 0;
  
  // find number of active links and sum of their weights
  (*sumW) = 0.;
  if (minerals != NULL) {
    minerals->SA_tot = 0.;
    for (int i = 0; i < minerals->n_tot; ++i) {
      minerals->list[i].SA = 0.0;
      minerals->list[i].SA_weight = 0.0;
    }
  }
  
  for (int nl = 0; nl < links->nlinks; ++nl) { /* loop links */
    link = &links->list[nl];
    if (link->inert || link->ghost)
      continue;
    nlink++;
    (*sumW) += sys->w[link->dir];

    if (minerals != NULL) {
      // weight sum for individual minerals
      // Find which mineral present in this node that has the highest solid fraction
      int wall = link->wall;
      double max = minerals->list[0].sfrac[wall];
      int imax = 0;
      sfrac_sum = 0;
      for (int i = 0; i < minerals->n_tot; ++i) {
        if (minerals->list[i].sfrac[wall] > max) {
          max = minerals->list[i].sfrac[wall];
          imax = i;
        }
        sfrac_sum += minerals->list[i].sfrac[wall];
      }
#ifdef _RANDOM_MINERAL_
      // add SA to max mineral
      minerals->list[imax].SA += sys->w[link->dir];
      // If this is an active site: distribute SA to minerals according to fraction of total mineral content
      for (int i = 0; i < minerals->n_tot; ++i) {
        if (minerals->list[i].sfrac[wall] > sys->opt->magn_limit || minerals->nucleation_site[wall] == 1) {
          minerals->list[i].SA_weight += sys->w[link->dir] * minerals->list[i].sfrac[wall] / sfrac_sum;
        }
      }
#endif
    }
    
  } // end of link loop
  
  // communicate sumW and SA
  (*sumW) *= 6.0;
  mpi->sbuf_eff[0] = (*sumW);
  num_comm = 1;
  if (minerals != NULL) {
    for (int i = 0; i < minerals->n_tot; ++i) {
      minerals->list[i].SA *= 6.0;
      mpi->sbuf_eff[1 + i] = minerals->list[i].SA;
    }
    num_comm = 1 + minerals->n_tot;
  }
  MPI_Reduce(mpi->sbuf_eff, mpi->rbuf_eff, num_comm, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Bcast(mpi->rbuf_eff, num_comm, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  (*sumW) = mpi->rbuf_eff[0];
  if (minerals != NULL) {
    minerals->SA_tot = (*sumW);
    for (int i = 0; i < minerals->n_tot; ++i) {
      minerals->list[i].SA = mpi->rbuf_eff[1 + i];
    }
  }
  
#ifdef _RANDOM_MINERAL_
  // communicate SA_weight
  if (minerals != NULL) {
    for (int i = 0; i < minerals->n_tot; ++i) {
      minerals->list[i].SA_weight *= 6.0;
      mpi->sbuf_eff[i] = minerals->list[i].SA_weight;
    }
    MPI_Reduce(mpi->sbuf_eff, mpi->rbuf_eff, minerals->n_tot, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Bcast(mpi->rbuf_eff, minerals->n_tot, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    for (int i = 0; i < minerals->n_tot; ++i) {
      minerals->list[i].SA_weight = mpi->rbuf_eff[i];
    }
  }
#endif
}


// -------------------------------------------------------
// checks if a given node is in the inert (non-reactive)
// part of the system based on the inertx_inlet and inertx_outlet variable
// -------------------------------------------------------
int node_is_reactive(int nn, struct System *sys)
{
  int gnx = get_global_xcoord(nn, sys);
  if ((gnx < sys->opt->inertx_inlet) || (gnx > sys->MAX_N[0] - sys->opt->inertx_outlet)) {
    return (0);
  }
  else {
    return (1);
  }
}


// -------------------------------------------------------
// 
// 
// -------------------------------------------------------
void save_fluid_restart_files(struct Field *fluid, struct System *sys)
{

  char fn[50];
  FILE *fp;

  char path[100];
  strcpy(path, sys->files["out"].c_str());
  strcat(path, "/restart");

  // fluid
  sprintf(fn, "%s/%04d_restart_f.bfd", path, sys->mpi->my_rank);
  write_1d_array(fn, fluid->f, fluid->n_c*sys->max_n_tot*sys->n_Q);

  if (sys->mpi->my_rank == 0) {
    sprintf(fn, "%s/restart_data.dat", path);
    fp = my_fopen(fn, "w");
    fprintf(fp, "%.12e %d %.12e %.12e %.12e %d",
      fluid->gravity[0], sys->t, sys->time, sys->dt, sys->t_extra, sys->nwrite);
    fclose(fp);
  }

}

// -------------------------------------------------------
// 
// 
// -------------------------------------------------------
void save_chem_restart_files(struct Field *species, struct Minerals *minerals, struct System *sys)
{

  char fn[50];
  FILE *fp;
  int cc;

  char path[100];
  strcpy(path, sys->files["out"].c_str());
  strcat(path, "/restart");

  // species
  sprintf(fn, "%s/%04d_restart_g.bfd", path, sys->mpi->my_rank);
  write_1d_array(fn, species->f, species->n_c*sys->max_n_tot*sys->n_Q);

  // minerals
  sprintf(fn, "%s/%04d_restart_min.bfd", path, sys->mpi->my_rank);
  fp = fopen(fn, "wb");
  if (fp == NULL) {
    printf("Error opening %s\n", fn);
  }
  else {
    for (cc = 0; cc < minerals->n_tot; ++cc) {
      fwrite(minerals->list[cc].sfrac, sizeof(double), sys->max_n_tot, fp);
    }
  }
  fclose(fp);

}

// -------------------------------------------------------
// 
// 
// -------------------------------------------------------
void load_fluid_restart_files(struct Field *fluid, struct System *sys)
{

  char fn[50];
  FILE *fp;

  char path[100];
  strcpy(path, sys->files["out"].c_str());
  strcat(path, "/restart");

  // fluid
  sprintf(fn, "%s/%04d_restart_f.bfd", path, sys->mpi->my_rank);
  read_1d_array(fn, fluid->f, fluid->n_c*sys->max_n_tot*sys->n_Q);
  communicate(fluid->f, fluid->n_c, sys);

  sprintf(fn, "%s/restart_data.dat", path);
  fp = my_fopen(fn, "r");
  fscanf(fp, "%lf %d %lf %lf %lf %d",
    &(fluid->gravity[0]), &(sys->start_itr), &(sys->time), &(sys->dt), &(sys->t_extra), &(sys->nwrite));
  fclose(fp);

  if (sys->mpi->my_rank == 0) {
    printf("From restart-files: gx = %.3e, t: %d, time: %.3e, dt: %.3e, t_extra: %.3e, nwrite: %d\n", fluid->gravity[0], sys->start_itr, sys->time, sys->dt, sys->t_extra, sys->nwrite);
  }
}

// -------------------------------------------------------
// 
// 
// -------------------------------------------------------
//void load_chem_restart_files(struct Field *species, struct Minerals *minerals, struct System *sys)
void load_chem_restart_files(struct Field *species, struct System *sys)
{

  char fn[50];

  char path[100];
  strcpy(path, sys->files["out"].c_str());
  strcat(path, "/restart");

  // species
  sprintf(fn, "%s/%04d_restart_g.bfd", path, sys->mpi->my_rank);
  read_1d_array(fn, species->f, species->n_c*sys->max_n_tot*sys->n_Q);
  communicate(species->f, species->n_c, sys);

  /* // minerals */
  /* sprintf(fn, "%s%04d_restart_min.bfd", path, sys->mpi->my_rank); */
  /* fp = fopen(fn,"rb"); */
  /* if ( fp == NULL ) { */
  /*   printf("Error opening %s\n",fn); */
  /* } else { */
  /*   for( cc = 0; cc < minerals->n_tot; ++cc ) {  */
  /*     fread(minerals->list[cc].sfrac, sizeof(double), sys->max_n_tot, fp); */
  /*   } */
  /* } */
  /* fclose(fp); */
  /* communicate_sfrac(sys, minerals); */

}

// -------------------------------------------------------
// 
// 
// -------------------------------------------------------
//void load_mineral_restart_files(struct Field *species, struct Minerals *minerals, struct System *sys)
void load_mineral_restart_files(struct Minerals *minerals, struct System *sys)
{

  char fn[50];
  FILE *fp;
  int cc;

  char path[100];
  strcpy(path, sys->files["out"].c_str());
  strcat(path, "/restart");

  /* // species */
  /* sprintf(fn, "%s%04d_restart_g.bfd", path, sys->mpi->my_rank); */
  /* read_1d_array(fn, species->f, species->n_c*sys->max_n_tot*sys->n_Q); */
  /* communicate(species->f, species->n_c, sys); */

  // minerals
  sprintf(fn, "%s/%04d_restart_min.bfd", path, sys->mpi->my_rank);
  fp = fopen(fn, "rb");
  if (fp == NULL) {
    printf("Error opening %s\n", fn);
  }
  else {
    for (cc = 0; cc < minerals->n_tot; ++cc) {
      fread(minerals->list[cc].sfrac, sizeof(double), sys->max_n_tot, fp);
    }
  }
  fclose(fp);
  communicate_sfrac(sys, minerals);

}


