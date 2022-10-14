#ifndef CHEM_GLOBAL_H
#define CHEM_GLOBAL_H
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <string>
#include "defines.h"
#include "chem_defines.h"
#include "splaytree.h"

struct treeNode;
struct splayTree;
struct System;

struct ChemTable {
	char **row_name;  /* name of complexes or rock buffers */
	char **col_name;  /* name of col names */
	real *log_m;      /* vector of log concentrations */
	real *log_a;      /* vector of log activities */
	real *log_g;      /* vector of log activity const */
	real *log_af;     /* vector of log activity const-buffers */
	real *delta;	  /* amount of precipitated mineral */
	real *T;          /* Temperature row vector */
	real *log_QK;     /* log_Q/K value for mineral, usually 0 */
	real *mol_volume; /* mol weight of mineral = molecular weight/density */
	real *mol_weight; /* mol weight */
	real *a0;         /* Debue Huckel constant */
	real *charge;     /* charge of complex */
	real *scharge;    /* charge of surface complex */
	real **M;         /* stoichometric matrix */
	real *logK;       /* logK values for each complex or mineral phase */
	int *abs_pos;     /* absolute position in full database */
	int *abs_pos_bas; /* absolute position in full basis vector */
	int size[2];      /* dimension of matrix */
	int *type;        /* type = 1 surface complex type = 0 aqueous complex */
} ;
/* Struct holding physical/chemical parameters */
struct ChemParam{
	real kB; /*Boltzmanns konstant : J/K*/
	real F;  /*Faradays constant C/mol */
	real e0; /*Permittivity of vacuum A s /(V m)*/
	real ew;  /*Relative permittivity of water */
	real Rg; /*Ideal gas constant*/
	real Na; /*Avogadros number*/
	real rho_w; /*density of water */
	real Adh;
	real Bdh;
	real beta_chem;
	real SA;  /*surface area */
};

/* memory space for calculations */
struct MEMC{
	int sizeV, sizeSM;
	real *dVma;
	real *dVmb;
	real *dVmc;
	real *dSMa;
	real **ddVma;
	real **ddVmb;
	real **ddVmc;
	int *iVma;
};
/* structure for initializing the chemical solver */
struct InitChem{
  int INTERPOLATE;
  int num_phases;
  int num_chrg_not_conv;  // keeps track of how often charge balance do not converge
	real Temp;
	int *pos_rock; 
	int  size_rock;
	int *pos_buffer;
	int  size_mass;
	int *pos_mass;
	int size_sup_min;
	int *pos_sup_min;
	int *pos_sup_bas;
	char **species_name;
	char **buffer_name;
	char **sup_min_name;
	int  *pos;
	int size;
	int Nm; /* No of minerals           */
	int Na; /* No of adsorption sites   */
	int Nx; /* No of ion exchange sites */
	int nlinr; /* = 1 non-linear rate equations used = 0 linear rate used */ 
	int PRINT_DEBUG_CHEM;/* print more info for chem solver */
	real *c_vchem;
	real *c_buffer;
	real *c_sup_min;
	real *c_ads;
	real *log_af;
	real *log_af_sup;
	real *Sg;
	real **rate; /* list of rate constants on set for each mineral */
	struct ChemTable *SM_mineral; /* mineral database */
	struct ChemTable *SM_all;     /* aqueous complexes database */
	struct ChemTable *SM_basis;   /* basis species database */
	struct ChemTable *SM_aq_logK; /* logK aqueous species */
	struct ChemTable *SM_M_logK;  /* logK minerals */
	struct MEMC *mem;             /* memory space for calculations */
	struct ChemParam *CP; /* physical/chemical param */
	struct BioParam *BiP; /* AMIN */
};
struct RateParam{
	real dt;
	real porosity;
	real *mol;
};

// AMIN 
struct BioParam{
	real mu_max;
	real lambda;
	real K_i;
	real K_O2;
	real alpha;
	real beta;
	real yield;
	real M;
};


struct BasVec{
	int surface_flag; /* surface_calc ? */
	int chbal;     /* chbal = 1 use charge balance, chbal = 0 skip */
	int DL;        /* Diffusive Layer calculation = 1 */
	int nlinr;      /* = 1 if non linear rate equations used, 0 othervise */
	int size;      /* number of basis species */
	int size_rock; /* size of rock buffer */
	int size_mass; /* size of mass balance  */
	int *pos_buffer; /* pos rock buffer basis */
	int *pos_rock; /* pos to the rock buffered species */
	int *pos_mass; /* pos to the species determined by mass balance */
	int *pos_sup_min; /* pos to supersaturated minerals */
	int *pos_sup_bas; /* pos to basis species for supersaturated minerals */
	int *pos_rel_sup_min; /* pos rel to input file */
	int size_sup_min; /* size of supersaturated minerals */
	int pos_water; /* pos to H20 */
	int pos_pH;    /* pos to H+ - always determined by charge balance */
	int pos_exp;   /* pos to exp - to be used for surface charge */
	int pos_X;	   /* pos to X - to be used for ion exchange     */
	int equilibrate; /* equilebrate with ionexchanger */
	real Temp;     /* temperature */
	real Io;       /* Ionic strength */
	real sch;      /* surface charge */
	real psi;      /* surface potential */
	real WHTO;     /* conc of water */
	real **beta_bas_inv; /* inverse of basis switching matrix */
	real **sm_buf; /* reduced stoichiometric matrix of basis buffers */
	real **sm_sup_buf; /* reduced stoichiometric matrix of basis sup sat min */
	real *ctot;	   /* total concentrations of basis species */
	real *ctot_calc; /* total calc concentrations of basis species */
	real *ctot_mineral; /* total concentration of mineral phase */
	real *delta_mineral; /* amount of mineral precipitated */
	real *ctot_ads;   /* tot concentration at the surface */
	real *ctot_ads_calc;   /* calc concentration at the surface */
	real *ctot_dl;  /*concentration in the diffuse layer */
	real *log_m;   /* log10-consentration of basis specise */
	real *log_a;   /* log10-activity of basis species */
	real *log_g;   /* log10-activity coefficient */
	real **rate;   /* rate konstants for each mineral phase */
	
	struct InitChem *ICS;
	struct BasVec *Vsurf; /*surface node */
	struct BasVec *Vbulk; /* only bulk species */
	struct RateParam *RP;
};

real *pick_row(struct ChemTable  *, char *);
real *pick_col(struct ChemTable  *, const char *);
void calc_logK(struct ChemTable *);
void get_logK(real Temp, real *logK, struct ChemTable *tab_all);
void set_chem_param(struct ChemTable *);
real *remove_vec_element(int, real *, real *);
//double zeroin(double x_init, int indx, double ax, double bx, double (*f)( double , int, struct BasVec *V), double tol, struct BasVec *V);
void ludcmp(double **, int , int *, double *);
void lubksb(double **, int, int *, double b[]);
double f_pH(double x, int no_call, struct BasVec *Vchem, void (*mbal)(int, struct BasVec *));
int charge_balance(int, struct BasVec *V);
void mass_balance(int no_call, struct BasVec *Vchem);
void calc_complex(struct ChemTable *, struct BasVec *V);
void newton_ah_log(real *x, real *, real **,struct BasVec *V);
double grahame(double x, int no_call, struct BasVec *Vchem);
//double surface_charge(struct BasVec *V);

double f_pH_newton(double, double *, struct BasVec *V);
int charge_balance_newton(int indx, struct BasVec *V);
double grahame_newton(double, double *, struct BasVec *V);
double surface_charge_newton(struct BasVec *V);
void init_chemistry(struct InitChem *ICS, std::string &data_folder);
/* double zeroin2( ax, bx,f,tol,Vchem)	; */
void calc_ads(struct BasVec *Vchem);

/* UTIL functions */
int in_list(int a, int *list, int size_list);
void invert_matrix(real **a, int dim, real **a_inv, real *d);
void in_fmatrix(real **fmatrix, int *row, int dim_row, int *coloum, int dim_col, real **f_red);
char *myreadline(char *, int, FILE *);
void change_basis_database(char *old_db, char *old_bname_t, char *new_bname_t, struct ChemTable *SM_Old_db);
int *find_elements_db(int col, struct ChemTable *db, char *list_t, int *pos, int *size);
void get_mol_volume(struct ChemTable *tab, char *fname);
void free_db(struct ChemTable *tab);
int call_chem_solver(real *ctot, struct InitChem *ICS, struct BasVec *Vchem);
void free_mem(struct MEMC *mem);
void InitMem(struct MEMC *mem);

// copied from global.h
int atoint(char *s);
int isdbl(char *s);
int getword(char *s, char *t);
int fgetline(FILE *fp, char s[]);
int cleanline(char s[]);
int isint(char *s);

/* IO-STUFF */
void read_database(std::string &folder, const std::string &filename, struct ChemTable *tab);
void remove_element(struct ChemTable *tab, real *vec, int start_row, int start_col, struct ChemTable *CM);
void remove_element_rc (struct ChemTable *tab, real *vec_row, real *vec_col, int start_row, int start_col, struct ChemTable *CM);
FILE* my_fopen(const char *filename, const char *mode);
void write_BasVec_struct(struct BasVec *Vchem, const char *name);
void write_chem_struct(struct ChemTable *tab, const char *name);
long myreadline_indx(char *line, int maxline, FILE *input_file);
char *myreadline_loc(char *line, int maxline, FILE *input_file);
struct InitChem myread_chem(char *line, int maxline,FILE *input_file, int no_phases, real Temp, long indx);
struct InitChem chem_read_input(char *fname);
/* struct InitChem myread_chem_surf(line,maxline,input_file, n_c_g_, Temp,indx); */
int myread_rate(struct InitChem *ICS, char *line, int maxline, FILE *input_file, long indx);
void copy_ChemTable_vector(real **vec_out, struct ChemTable *CM, real *vec, int n);
void icopy_ChemTable_vector(int **vec_out, struct ChemTable *CM, int *vec, int n);
int *ipick_col(struct ChemTable  *tab, const char *name);
void writeVchem(struct BasVec *Vp, const char *name);
/* INITIALIZATION */
void init_ctot(struct BasVec *Vchem);
void make_basvec_struct( struct InitChem *ICS, struct BasVec *Vchem);
void make_surface_bulk_struct(struct InitChem *ICS, struct BasVec *Vchem);
struct BasVec make_basvec_struct_old( char *rock_t, char *buffer_t, char *mass_t, char *mineral_sup_t);
struct BasVec chem_init_boundary_node( struct InitChem *ICS );
void  init_const_chem_param(struct ChemParam *CP);
void InitICS(struct InitChem *ICS, int size);
real chem_set_eq_conc_nlin(real *c_vchem_eq, struct InitChem *ICS, int call);
void trans_nlin_eq(struct BasVec *Vchem);
/* STRUCTURE and ORGANIZE INFORMATION */
void add_new_buffer_mineral(int ICSrock, int ICSbuffer,   struct BasVec *Vchem);
void remove_buffer_mineral(int length, struct BasVec *Vchem);
int check_if_mineral_supersaturated(int *mineral_list, int size, int *old_list, int old_size, struct BasVec *Vchem);
int check_if_mineral_conc_neg(int *list, int size, struct BasVec *Vchem);
void remove_sup_mineral(int length,   struct BasVec *Vchem);
int chem_string_to_pos(struct ChemTable *sm, char *list_t, int **pos);
void update_sm_buf(struct BasVec *Vchem);
void update_sm_buf_new(struct BasVec *Vchem, int size_r, int *pos_b, int *pos_r, real **sm);
void update_beta_inv( struct BasVec *Vchem);
real chem_set_eq_conc(real *c_vchem_eq, struct InitChem *ICS, int call);
void remove_species(struct InitChem *ICS_new, char *name, struct InitChem *ICS_old);
void set_ctot(real *ctot, struct BasVec *Vchem);
void chem_get_pos(struct InitChem *ICS);
void remove_mineral_buffer(int cw, struct BasVec *Vp);
void add_new_sup_mineral(int length, struct InitChem *ICS,   struct BasVec *Vchem);
void db_set_activity(struct InitChem *ICS);
void make_surface_bulk_struct(struct InitChem *ICS, struct BasVec *Vchem);
void make_bulk_struct(struct InitChem *ICS, struct BasVec *Vchem);
void make_surface_bulk_struct_split(struct InitChem *ICS, struct InitChem *ICS_s, struct InitChem *ICS_b, struct BasVec *Vchem);
void splitICS(struct InitChem *ICS, struct InitChem **ICS_s, struct InitChem **ICS_b);
/* CALCULATIONAL */
int chem_calc_eq_rho_lb_moving(real *rho_eq_wall, real *rho_wall, struct InitChem *ICS,
                               real *min_conc, struct treeNode **tn_tmp, struct splayTree **st_lst, struct System *sys);
void calc_act_DB(real A, real B, real Io,struct BasVec *Vchem);
real calc_Io(struct BasVec *Vchem);
void calc_rock_spec_conc(struct BasVec *Vchem);
void jacobi_num(real *x, real **j_num, struct BasVec *Vchem);
void calc_ctot(struct BasVec *Vchem);
void calc_ctot_aq(struct BasVec *Vchem);
void calc_mineral_precip(struct BasVec *Vchem);
int charge_balance_newton(int indx, struct BasVec *Vchem);
double f_pH_newton(double x, double *df, struct BasVec *Vchem);
int charge_balance_secant(int indx, struct BasVec *Vchem);
void solve_chem_bulk_surface_new(real *ctot_aq, int call, struct BasVec *Vchem);
void solve_chem_bulk_surface(real *ctot_aq, int call, struct BasVec *Vchem);
/* ----- RATE ------*/
void calc_F_sup_min(real *Fchem, struct BasVec *Vchem);
void calc_ctot_sup_min(real *Fchem, struct BasVec *Vchem);
void calc_jacobi_sup_min(real **j, struct BasVec *Vchem);
void set_rate_val(int min_no, double *da, struct InitChem *ICS);
/*-----------------*/
int  solve_chem_boundary_node(real *ctot, int call, struct BasVec *Vchem);
void solve_chemistry( struct BasVec *Vchem);
void calc_surface_reaction(real *ctot_aq, int call, struct BasVec *Vchem);
void calc_surface_reaction_nlin(real *ctot_aq, int call, struct BasVec *Vchem);

void solve_chemistry_splaytree(real *in, real *out, struct BasVec *Vchem);
int my_floor(real x);
int mycomp(int * key1, int * key2);

void newton_ah_log_surf(real *x, real *Fchem, real **jacobi, struct BasVec *Vchem);
void mass_balance_surf(int no_call, struct BasVec *Vchem);
void newton_ah_log_io(real *x, real *Fchem, real **jacobi, struct BasVec *Vchem);
void get_logK(real Temp, real *logK, struct ChemTable *tab_all);
void set_temperature_db(struct BasVec *Vchem);

/* RATE EQUATIONS */
void rk4(double t, double *y, int dimension, double h,
	 double (*func)(double , double *, int , struct BasVec *), struct BasVec *Vchem);
double rate_eq(double t, double *y, int i, struct BasVec *Vchem);
void rate_to_eq(struct BasVec *Vchem);
int reaction_rate(real *c_updated, struct BasVec *Vchem);
void reaction_rate_lb(real *c, struct BasVec *Vchem);
void reaction_rate_julius(real *c_updated, struct BasVec *Vchem);
void update_conc(real *c_updated, int pos, real rate, struct BasVec *Vchem);

/* HELP FUNCTIONS */
void jacobi_num2(real *x, real *dx, int dim, struct BasVec *Vchem, real *chem_F, real **jc);
void mass_balance_num(int no_call, struct BasVec *Vchem);
void print_fpH(real x1, real x2, struct BasVec *Vchem);
void mass_balance_num(int no_call, struct BasVec *Vchem);
void newton_ah_log_num(real *x, real *Fchem, struct BasVec *Vchem);
int get_int_from_string(char *s);
void compare_names(struct ChemTable *Ma, struct ChemTable *Mb);
int find_no_basis_species(int pos_buf, int *pos_bas, int *list, int size_list, struct ChemTable *db);
void write_jacobi(real **j, int size, const char *fname);

/*SURFACE CHEMISTRY */
real qgaus(real a, real b, int pos, struct ChemTable * Tab, struct BasVec *Vchem);
real trapez(real (*func)(real, int, struct ChemTable *, struct BasVec *), real a, real b, int pos, struct ChemTable * Tab, struct BasVec *Vchem);
real grahame_DL(real xi, int pos, struct ChemTable *Tab, struct BasVec *Vchem);
void calc_diffuse_layer_conc(real *c_dl, struct BasVec *Vchem);

/* Adsorption - Ion exchange */
void calc_adsorption(struct BasVec *Vchem_surf, struct BasVec *Vchem);

/* MEMORY */
void free_BasVec(struct BasVec *Vchem);
void free_ICS(struct InitChem *ICS);

/*INTERPOLATE */
void gen_binary(int *b, int num, int length);
void gen_vchem_vector(int no_comb, struct BasVec *Vchem, struct InitChem *ICS);
int gen_new_key(real *ctot_min, int size);
int solve_chem_boundary_node_interpolate(real *ctot, real *ctot_min, int key, struct BasVec *Vchem);
void copy_BasVec(struct BasVec *V_out, struct BasVec *V_in);
 void calc_min(real *dc, real *dc_min, struct BasVec *Vchem);
 void update_BasVec(real *c_new, real *dc_min, struct BasVec *Vchem);

//------------------------------------------------//
//                julius_chem_util.c              //
//------------------------------------------------//
//real calc_rock_density(real *dx, struct InitChem *ICS);
//real calc_surface_area(real *dx, real por, struct InitChem *ICS);
//real calc_permeability(real *dx, real tau, real por, struct InitChem *ICS);
//void set_mineral_conc(real rho_r, real por, struct BasVec *Vchem);
//void calc_wt_change(real *dx, real rho_r, real por,  struct BasVec *Vchem);
//real calc_por_change(real *dphi, real *dx, real rho_r, real por,  struct BasVec *Vchem);
void calculate_mol_weight_mineral(struct ChemTable *Min, struct ChemTable *Bas);
//real equilibrate_pore_water(real *ctot, struct BasVec *Vchem);
void add_H_to_mbal(real c_h, struct BasVec *Vchem);
void add_H_to_ICS(struct InitChem *ICS, int posH);
//real anal_perm_mod(real time, real b, real m, real c);
void calc_mineral_change_julius(struct BasVec *Vchem);


#endif
