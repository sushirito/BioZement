#include "global.h"
#include "chem_global.h"
#include <ctype.h>  // isspace()
//#ifdef _WIN32
//#include <io.h>
//#else
//#include <unistd.h> // dup()
//#endif

void read_input_v2(Input &inp, struct InitChem **ICS_ptr, struct System *sys,
  struct Field *fluid, struct Field *species, struct Field *dif, struct Boundary *bndry, char *argv[])
{
  for (int i = 0; i < 3; ++i) {
    sys->mpi->np[i] = inp["mpi"]["num_proc"][i];
  }

  sys->real_units_in_input = inp["real_units"];
  sys->n_D = inp["dim"]["nD"];
  if (sys->n_D > 3 || sys->n_D < 2) {
    printf("ERROR in input.dat: sys->n_D = %d [2,3]", sys->n_D);
    MPI_Finalize();
    exit(1);
  }
  sys->n_Q = inp["dim"]["nQ"];

  // Read basis
  std::stringstream basis_file;
  basis_file << sys->files["basis"] << "/d" << sys->n_D << "q" << sys->n_Q << ".basis";
  myread_basis(basis_file.str().c_str(), sys);

  set_dimensions_from_geofile(sys);

  sys->start_itr = 0; //inp["iterations"]["start"];
  sys->n_itr = inp["iterations"]["max"];
  sys->end_itr = sys->start_itr + sys->n_itr;
  sys->write_int = inp["iterations"]["write"];
  fluid->u_Darcy = inp["fluid"]["u_Darcy"];
  fluid->gravity = (double *)calloc(sys->n_D, sizeof(double));
  fluid->gravity[0] = inp["fluid"]["force_x"];
  sys->gx             = fluid->gravity[0];
  sys->const_flow     = inp["fluid"]["const_flow"];
  sys->max_vel        = inp["fluid"]["max_vel"];
  sys->mean_vel       = inp["fluid"]["mean_vel"];
  sys->flux           = sys->mean_vel;
  sys->rate_constant  = 0.0; 

  strcpy(sys->drive_type,"gx");
  sys->drive = &fluid->gravity[0];
  if (sys->max_vel>0.0) {
    sys->vel_target = &sys->max_vel;
  } else {
    sys->vel_target = &sys->mean_vel;    
  }
  //fluid->gx_change = fluid->gravity[0];
#ifdef _FLUID_SOURCE_
  strcpy(sys->drive_type,"q_source");
  sys->drive = &fluid->q_source;
#endif
#ifdef _FLUID_BDRY_PRESS_
  fluid->gravity[0] = 0.0;
  fluid->grad_P = inp["fluid"]["grad_P"];
  strcpy(sys->drive_type,"grad_P");
  sys->drive = &fluid->grad_P;
#endif
  
  fluid->n_c         = inp["fluid"]["phases"]["num"];
  fluid->q_source    = inp["fluid"]["q_source"];
  fluid->tau         = (double *) calloc(fluid->n_c, sizeof(double));
  fluid->rho_init    = (double *) calloc(fluid->n_c, sizeof(double));            // initial density
  fluid->rho_inlet   = (double *) calloc(fluid->n_c, sizeof(double));            // inlet density
#ifdef _FLUID_BDRY_PRESS_
  fluid->rho_outlet  = (double *) calloc(fluid->n_c, sizeof(double));            // outlet density
#endif
  //  fluid->u_init      = (double *) calloc(sys->n_D*fluid->n_c, sizeof(double));   // initial velocity
  for (int i=0; i<fluid->n_c; ++i) {
    fluid->tau[i]       = inp["fluid"]["phases"]["tau"][i];
    fluid->rho_init[i]  = inp["fluid"]["rho"]["init"][i];
    fluid->rho_inlet[i] = inp["fluid"]["rho"]["inlet"][i];
  }
  fluid->phase_perm = std::vector<double>(fluid->n_c, 0);
  fluid->saturation = std::vector<double>(fluid->n_c, 0);

  // kinematic viscosity
  fluid->nu_lb = std::vector<double>(fluid->n_c);
  for (int c=0; c<fluid->n_c; ++c)
    fluid->nu_lb[c] = sys->cs_2*(fluid->tau[c] - 0.5);
  sys->nu_lb = sys->cs_2*(fluid->tau[0] - 0.5);


#ifdef _USE_COLOR_GRAD_
  //fluid->inlet_pressure = inp["fluid"]["P_in_color_grad"];
  fluid->beta = inp["fluid"]["beta"];
  fluid->surface_tension_multiphase_matrix = inp["fluid"]["surface_tension"].get_symmetric_matrix();
  sys->mean_vel = double(inp["fluid"]["Ca"])*fluid->surface_tension_multiphase_matrix[0][1]/fluid->nu_lb[0];
  fluid->rho_wall = inp["fluid"]["rho"]["wall"].get_vector();
  fluid->u_mean = std::vector<double>(sys->n_D, 0);
  fluid->phase_u_mean = std::vector<double>(fluid->n_c*sys->n_D, 0);
  //  for (auto &val : fluid->rho_wall)
  //    std::cout << val << " ";
  //  std::cout << std::endl;
  //  std::cout << std::endl << "SURFACE_TENSION_MULTIPHASE_MATRIX:" << std::endl;
  //  for (auto &row : fluid->surface_tension_multiphase_matrix) {
  //    for (double &val : row) {
  //      std::cout << val << " ";
  //    }
  //    std::cout << std::endl;
  //  }
#endif

  // DIFFUSIVE FIELD [Amin]
  dif->n_c = inp["dif"]["num"];
  dif->tau = (double *)calloc(dif->n_c, sizeof(double));
  for (int i = 0; i < dif->n_c; i++) {
      dif->tau[i] = inp["dif"]["tau"][i];
      std::cout << "tau[" << i << "] = " << dif->tau[i] << std::endl;
  }
  dif->rho_inlet = (double *)malloc(sizeof(double)*dif->n_c);
  dif->rho_init = (double *)malloc(sizeof(double)*dif->n_c);
  dif->sum_conc = (double *)malloc(sizeof(double)*dif->n_c);
  dif->sum_conc_old = (double *)malloc(sizeof(double)*dif->n_c);
  dif->sum_conc_prev = (double *)malloc(sizeof(double)*dif->n_c);
  dif->mean_conc = (double *)malloc(sizeof(double)*dif->n_c);
  dif->mean_conc_old = (double *)malloc(sizeof(double)*dif->n_c); // used in convergence routine
  dif->mean_g = (double *)malloc(sizeof(double)*dif->n_c); // used to calculate source term in collision.c      
      
  // AMIN EJE
  std::cout << "dif->n_c = " << dif->n_c << std::endl;
  for (int i=0; i < dif->n_c; ++i)
      dif->rho_init[i] = inp["dif"]["rho"]["init"][i];  
  std::cout << "end n_c" << std::endl; 
      
  // SPECIES AND CHEMISTRY
  if (sys->chem == 0)
    return;

  if (inp["chem"]["inlet"].nrows() != inp["chem"]["bulk"].nrows()) {
    std::cerr << std::endl << "ERROR! Number of inlet- and outlet-species are different! Check input.dat" << std::endl << std::endl;
    exit(1);
  }

  sys->D = inp["chem"]["D"];
  sys->D_b = inp["dif"]["D"]; // Diffusion coefficient for biomass field /* AMIN */

  species->n_c = inp["chem"]["inlet"].nrows();
  species->tau = (double *)calloc(species->n_c, sizeof(double));
  for (int i = 0; i < species->n_c; ++i)
    species->tau[i] = inp["chem"]["tau"];

  bndry->n_liquid = species->n_c;

  std::vector<std::string> name;
  //*ICS_ptr = (struct InitChem *) calloc(3, sizeof(struct InitChem));
  *ICS_ptr = new InitChem[3]();
  struct InitChem *ICS = *ICS_ptr;

  // set inlet values
  InitICS(&(ICS[0]), species->n_c);
  ICS[0].Temp = double(inp["chem"]["temp"]) + 273.15; // Kelvin
  name = inp["chem"]["inlet"].get_names();
  for (int i = 0; i < species->n_c; ++i) {
    strcpy(ICS[0].species_name[i], name[i].c_str());
    ICS[0].c_vchem[i] = inp["chem"]["inlet"][name[i]];
    ICS[0].pos_mass[ICS[0].size_mass] = i;
    ICS[0].size_mass++;
  }

  // set bulk values
  InitICS(&(ICS[1]), species->n_c);
  ICS[1].Temp = double(inp["chem"]["temp"]) + 273.15; // Kelvin
  name = inp["chem"]["bulk"].get_names();
  for (int i = 0; i < species->n_c; ++i) {
    strcpy(ICS[1].species_name[i], name[i].c_str());
    ICS[1].c_vchem[i] = inp["chem"]["bulk"][name[i]];
    ICS[1].pos_mass[ICS[1].size_mass] = i;
    ICS[1].size_mass++;

  }

  //set biomass parameters /* AMIN */
  // BiP.yield = inp["biomass"]["dissolution"]["yield"];
// AMIN
  name = inp["biomass"]["precipitation"].get_names();
  ICS[1].BiP->yield = inp["biomass"]["precipitation"][name[0]];
  ICS[1].BiP->mu_max = inp["biomass"]["precipitation"][name[1]];
  ICS[1].BiP->K_i = inp["biomass"]["precipitation"][name[2]];
  dif->rho_max = inp["biomass"]["precipitation"][name[3]];

  // set reaction rates
  name = inp["chem"]["rate"].get_names();
  ICS[1].size_sup_min = name.size();
  std::vector<double> value;
  for (int i = 0; i < ICS[1].size_sup_min; ++i) {
    strcpy(ICS[1].sup_min_name[i], name[i].c_str());
    value = inp["chem"]["rate"][name[i]];
    ICS[1].c_sup_min[i] = value[0];
    ICS[1].Sg[i] = value[1];
    ICS[1].log_af_sup[i] = value[2];
    ICS[1].rate[i][0] = value[3];
    ICS[1].rate[i][1] = value[4];
    ICS[1].rate[i][2] = value[5];
    ICS[1].rate[i][3] = value[6];
    ICS[1].rate[i][4] = value[7]; // for variable rate /* AMIN */
    ICS[1].rate[i][5] = value[8]; // for variable rate /* AMIN */
  }
  ICS[1].nlinr = 1;

  if (ICS[0].size != species->n_c) {
    printf("NB! Du har satt antall kjemi arter %d og antall lb_spec %d\n", ICS[0].size, species->n_c);
  }
  if (ICS[0].size != ICS[1].size) {
    printf("Feil! Du har satt antall inlet kjemi arter %d og antall bulk arter %d\n", ICS[0].size, ICS[1].size);
    MPI_Finalize();
    exit(1);
  }

  ICS[0].INTERPOLATE = ICS[1].INTERPOLATE = sys->opt->splaytree_ON;
  ICS[0].num_chrg_not_conv = ICS[1].num_chrg_not_conv = 0;

  species->rho_inlet = (double *)malloc(sizeof(double)*(ICS[0].size + 1));
  species->rho_init = (double *)malloc(sizeof(double)*(ICS[1].size + 1));
  species->sum_conc = (double *)malloc(sizeof(double)*(ICS[1].size + 1));
  species->sum_conc_old = (double *)malloc(sizeof(double)*(ICS[1].size + 1));
  species->sum_conc_prev = (double *)malloc(sizeof(double)*(ICS[1].size + 1));
  species->mean_conc = (double *)malloc(sizeof(double)*(ICS[1].size + 1));
  species->mean_conc_old = (double *)malloc(sizeof(double)*(ICS[1].size + 1)); // used in convergence routine
  species->mean_g = (double *)malloc(sizeof(double)*(ICS[1].size + 1)); // used to calculate source term in collision.c


  // AMIN EJE
  for (int i=0; i < species->n_c; ++i)
      species->rho_init[i] = ICS[1].c_vchem[i];
}



/* --------------------------------------------------------------------------------- myread_basis */
void myread_basis(const char *fname, struct System *sys)
/*
myread_basis :
    reads the d#q# basis from file. Sets the basis related constants

  INPUT  : 1) fname : name to basis file

  OUTPUT :	 void
*/
{
  FILE *input_file;  /* input file pointer */
  char ct[] = " \t \0"; /* tokens used in reading a list of numbers from a string */
  int maxline = 200; /* maximum line length */
  int cd, ck;        /* Counter dimension counter direction */
  int *tmp_w, *tmp_ei;      /* Temporal vectors */

  /* Open input file */
  input_file = my_fopen(fname, "r");

  /* Read weights (double) */
  /* NB! the first entry is the largest commen divisior for the denominator */
  /*     so the weights are given by tmp_w[ck+1]/tmp_w[0]                   */
  myread_int_vec(&tmp_w, sys->n_Q + 1, ct, input_file, maxline);
  sys->w = (real *) calloc(sys->n_Q, sizeof(real));
  for (ck = 0; ck < sys->n_Q; ++ck) {
    sys->w[ck] = (1.*tmp_w[ck + 1]) / tmp_w[0];
  }


  /* Read basis vectors (int) */
  sys->ei = (int *) calloc(sys->n_Q*sys->n_D, sizeof(int));
  for (cd = 0; cd < sys->n_D; ++cd) {
    myread_int_vec(&tmp_ei, sys->n_Q, ct, input_file, maxline);
    for (ck = 0; ck < sys->n_Q; ++ck) {
      sys->ei[sys->n_D*ck + cd] = tmp_ei[ck];
    }
    free(tmp_ei);
  }

  /* Set velocity vector */
  sys->ev = (real *) calloc(sys->n_Q*sys->n_D, sizeof(real));
  for (ck = 0; ck < sys->n_Q; ++ck)
    for (cd = 0; cd < sys->n_D; ++cd)
      sys->ev[sys->n_D*ck + cd] = (real)sys->ei[sys->n_D*ck + cd];


  /* Check symmetries */
  if (sys->chem) {
    if (check_basis(tmp_w, sys) < 2) {
      printf("ERROR in the basis symmetry\n");
      exit(1);
    }
  }
  /* else if ( sys->fluid ) { */
  if (sys->fluid) {
    if (check_basis(tmp_w, sys) < 6) {
      printf("ERROR in the basis symmetry\n");
      exit(1);
    }
  }

  /* Check and set constants */
  calc_cs_2(tmp_w, sys);

  /* Set bounce back directions */
  init_k_bb(sys);

  /* Set constants used in the equlibrium distributions */
  calc_fg_eq_c(sys);

  /* Free memory */
  free(tmp_w);

}




/* --------------------------------------------------------------------------------- myread_int_vec */
void myread_int_vec(int **vec, int n_e, char *ct, FILE *fp, int maxline)
/*
myread_int_vec :
    reads a line with multiple integer entries, and creates an
    vector containing the entries.

  INPUT  : 1) vec : pointer to the vector pointer
           2) n_e : number of elements in the vector
           3) ct  : deliminater tokens
           4) fp  : file pointer
           5) maxline : maximum number of character in a file

  OUTPUT : void
*/
{
  int ce; /* Counter vector elements*/
  char line[200];  /* holds an uncommented line from the inputfile */
  char *s_tmp; /* temporary string */

  /* allocate memory */
  if (!(*vec = (int *)calloc(n_e, sizeof(int)))) {
    printf("Error allocating memory\n");
    exit(1);
  }

  /* read line */
  myreadline_loc(line, maxline, fp);
  /* -- first entry */
  if (!(s_tmp = strtok(line, ct))) {
    printf("error reading inputfile");
    exit(1);
  }
  sscanf(s_tmp, "%d", &((*vec)[0]));

  /* consecutive entries */
  for (ce = 1; ce < n_e; ++ce) {
    if (!(s_tmp = strtok(NULL, ct))) {
      printf("error reading inputfile");
      exit(1);
    }
    sscanf(s_tmp, "%d", &((*vec)[ce]));
  }
}


// /* --------------------------------------------------------------------------------- myread_double_vec */
// void myread_double_vec(double **vec , int n_e, char *ct, FILE *fp, int maxline)
// /*
// myread_double_vec :
// 	reads a line with multiple double entries, and creates an
// 	vector containing the entries.

//   INPUT  : 1) vec : pointer to the vector pointer
// 		   2) n_e : number of elements in the vector
// 		   3) ct  : deliminater tokens
// 		   4) fp  : file pointer
// 		   5) maxline : maximum number of character in a file

//   OUTPUT : void
// */
// {
// 	int ce; /* Counter vector elements*/
// 	char line[200];  /* holds an uncommented line from the inputfile */
// 	char *s_tmp; /* temporary string */

// 	/* allocate memory */
// 	if( !(*vec = (double *) calloc(n_e, sizeof(double))) )
// 	{
// 		printf("Error allocating memory\n");
// 		exit(1);
// 	}

// 	/* read line */
// 	myreadline_loc(line, maxline, fp);
// 	/* -- first entry */
// 	if( !(s_tmp = strtok(line, ct)) )
// 	{
// 		printf("error reading inputfile");
// 		exit(1);
// 	}
//  	sscanf(s_tmp, "%lf", &((*vec)[0]));

// 	/* consecutive entries */
// 	for( ce = 1; ce < n_e; ++ce ) 
// 	{
// 		if( !(s_tmp = strtok(NULL,ct)) )
// 		{
// 			printf("error reading inputfile");
// 			exit(1);
// 		}
// 		sscanf(s_tmp,"%lf",&((*vec)[ce]));
// 	}	
// }

/* ---------------------------------------------------------------------------------  myreadline_loc */
char *myreadline_loc(char *line, int maxline, FILE *input_file)
/* reads the first line without '#' */
{
  char *ff;
  int i;
  for (i = 0; i < maxline; ++i)  line[i] = ' ';

  ff = fgets(line, maxline, input_file);
  while (ff != NULL) {

    if (line[0] != '#' && line[0] != '\n' && line[0] != '\r') {
      return ff;
      break;
    }
    for (i = 0; i < maxline; ++i)  line[i] = ' ';
    ff = fgets(line, maxline, input_file);
  }
  return ff;
}

// long myreadline_indx(char *line, int maxline, FILE *input_file)
// /* reads the first line without '#' */
// /* returns the position in the file */
// /* it is possible to go back to this position with fseek() */
// {
// 	char *ff;
// 	long indx = 0;
// 	int i, Nblanks;

// 	ff = (char *) calloc(maxline,sizeof(char));

// 	for( i=0; i<maxline; ++i) line[i]=ff[i]= ' ';

// 	fgets(ff, maxline, input_file);
// 	while(ff != NULL){

//     if( ff[0] != '#' && ff[0] != '\n' )
// 	{ /* remove blanks from the end of the line */
// 		Nblanks=strlen(ff);
// 		for(i=0;i< (int) strlen(ff);++i)
// 		{
// 			if(ff[strlen(ff)-1-i]== ' ' || ff[strlen(ff)-1-i]== '\n' 
// 				|| ff[strlen(ff)-1-i]== '\r')
// 			{
// 				Nblanks--;
// 			} else break; 
// 		}
// 			strncpy(line,ff,Nblanks);
// 			indx = ftell(input_file);
// 			free(ff);
// 		return indx;
// 		break;
//     }	
// 	for( i=0; i<maxline; ++i)  ff[i]= ' ';
//     fgets(ff, maxline, input_file);
// 	indx = ftell(input_file);
//   }
// 	free(ff);
//   return indx; 
// }


// /* ---------------------------------------------------------------------------------  myprint_int_vec */
// void myprint_int_vec(const char *vec_name, int *vec, int vec_size)
// /*
// myprint_int_vec :
// 	print the vector to the screen

//   INPUT  : 1) vec_name : name of the vector
// 		   2) vec      : pointer to the vector
// 		   3) vec_size : number of elements in the vector

//   OUTPUT : void
// */
// {
// 	int vn; /* coutner vector elements */

// 	if( vec_size == 1 )
// 	{
// 		printf("%-20s = %d\n", vec_name, vec[0]);
// 	} else {
// 		printf("%-20s = [%d", vec_name, vec[0]);
// 		for( vn = 1; vn < vec_size; ++vn )
// 			printf(", %d", vec[vn]);
// 		printf("]\n");
// 	}
// }


// /* ---------------------------------------------------------------------------------  myprint_double_vec */
// void myprint_double_vec(const char *vec_name, double *vec, int vec_size)
// /*
// myprint_int_vec :
// 	print the vector to the screen

//   INPUT  : 1) vec_name : name of the vector
// 		   2) vec      : pointer to the vector
// 		   3) vec_size : number of elements in the vector

//   OUTPUT : void
// */
// {
// 	int vn; /* coutner vector elements */
// 	if( vec_size == 1 )
// 	{
// 		printf("%-20s = %g\n", vec_name, vec[0]);
// 	} else {
// 		printf("%-20s = [%g", vec_name, vec[0]);
// 		for( vn = 1; vn < vec_size; ++vn )
// 			printf(", %g", vec[vn]);
// 		printf("]\n");
// 	}
// }

//----------------------------------------
// return 1 if string has only whitespace
//----------------------------------------
int is_empty(const char *s)
{
  while (*s != '\0') {
    if (!isspace((unsigned char)*s))
      return 0;
    s++;
  }
  return 1;
}

