/*
make_basvec_struct :
  pick out the rock buffered species, species buffered by the rock, and minerals in the rock 

  INPUT:  1) *rock_t        list of names of the species buffered by a mineral
          2) *buffer_t      list of minerals corresponding to a species given in rock_t
		  3) *mass_t        list of names of the basis species that are to be determined by mass balance
		  4) *mineral_sup_t list of supersaturated minerals that are to be used in rate dependent calculations,
		                    i.e. minerals that equilibrate slow
  RETURN: BasVec struct 
*/
#include "chem_global.h"

void make_basvec_struct( struct InitChem *ICS, struct BasVec *Vchem)
{ 
	//char ct[] = "\n \t \0"; /* MAC LINUX USE 0*/
	int full_basis_size,pos;
	//real nu;
	int full_buffer_size;
	int i,j, length, length_i, length_f, length_r, length_m;
	int *list, *list_b,size_list, size_list_b, no_i, pos_i, pos_b;
	struct ChemTable *SM_all, *SM_mineral, *SM_basis;
	//real c_dum[4], c_dum_min[4], por, rho_r;

	Vchem->ICS = ICS;

	SM_all     = Vchem->ICS->SM_all;
	SM_mineral = Vchem->ICS->SM_mineral;
	SM_basis   = Vchem->ICS->SM_basis;



	/* SM_mineral : Subset of full mineral database                */
	/* SM_all_     : Subset of full secondary (complexes) data base */
	/* SM_basis   : Subset of full basis database                  */ 
	full_basis_size=SM_all->size[1];
	full_buffer_size=SM_mineral->size[0];

	/* allocate space */
	Vchem->size           =  full_basis_size;
	Vchem->log_a          = (real *)  calloc(Vchem->size,sizeof(real ));
	Vchem->log_g          = (real *)  calloc(Vchem->size,sizeof(real ));
	Vchem->log_m          = (real *)  calloc(Vchem->size,sizeof(real ));
	Vchem->ctot           = (real *)  calloc(Vchem->size,sizeof(real ));
	Vchem->ctot_ads       = (real *)  calloc(Vchem->size,sizeof(real ));
	Vchem->ctot_dl        = (real *)  calloc(Vchem->size,sizeof(real ));
	Vchem->ctot_ads_calc  = (real *)  calloc(Vchem->size,sizeof(real ));
	Vchem->ctot_calc      = (real *)  calloc(Vchem->size,sizeof(real ));
	Vchem->ctot_mineral	  = (real *)  calloc(full_buffer_size,sizeof(real ));
	Vchem->delta_mineral  = (real *)  calloc(Vchem->size,sizeof(real ));
	Vchem->pos_rock       = (int *)   calloc(Vchem->size,sizeof(int));  
	Vchem->pos_buffer     = (int *)   calloc(Vchem->size,sizeof(int));
	Vchem->pos_mass       = (int *)   calloc(Vchem->size,sizeof(int));
	Vchem->pos_sup_min    = (int *)   calloc(Vchem->size,sizeof(int));
	Vchem->pos_sup_bas    = (int *)   calloc(Vchem->size,sizeof(int));
	Vchem->pos_rel_sup_min    = (int *)   calloc(Vchem->size,sizeof(int));
	Vchem->beta_bas_inv   = (real **) calloc(Vchem->size,sizeof(real *));
	Vchem->Vsurf          = NULL;
	Vchem->rate           = NULL;


	Vchem->RP = (struct RateParam *) calloc(1,sizeof(struct RateParam ));
	Vchem->RP->dt = 1.;
	Vchem->RP->porosity = 1.;
	Vchem->RP->mol = NULL;

	Vchem->WHTO = 1.;
	
			
	Vchem->pos_X   = -1;
	Vchem->pos_exp = -1;
	Vchem->pos_water = -1;
	Vchem->pos_pH = -1;
	Vchem->equilibrate = 0;
	Vchem->surface_flag = 0;
	Vchem->DL = 0;

	if(ICS->size_sup_min >0)
	{
		Vchem->RP->mol = (real *) calloc(ICS->size_sup_min,sizeof(real));
		for(i=0;i<ICS->size_sup_min;++i) Vchem->RP->mol[i] = 1.;
		Vchem->rate=(real **) calloc(4*ICS->size_sup_min, sizeof(real *));
		for(i=0;i<ICS->size_sup_min;++i) Vchem->rate[i] = (real *) calloc(4,sizeof(real));
	}
	for(i=0;i<Vchem->size;++i)
	{
		Vchem->beta_bas_inv[i] = (real *) calloc(Vchem->size, sizeof(real ));
	}
	/* allocate finnished */

	/* initialize */
	Vchem->sm_buf = NULL;

	for(i=0;i<Vchem->size;++i)
	{
		Vchem->log_a[i]=Vchem->log_m[i]=log10(CHEM_LOWER_);
		Vchem->log_g[i]=0.;
		Vchem->delta_mineral[i]=Vchem->ctot_calc[i]=Vchem->ctot_ads[i]=Vchem->ctot_ads_calc[i]= Vchem->ctot_dl[i]= 0.;
		Vchem->ctot[i] = 1.; /* we are using the log values in the calculations */
		Vchem->pos_rock[i]=Vchem->pos_buffer[i]=Vchem->pos_mass[i]=0;
	}
	
	/* special species : H20, E, H+ */
	length = 0;
	for(i=0;i<Vchem->size; ++i)
	{
		if(!strcmp(SM_basis->row_name[i],"H2O")) Vchem->pos_water = length;
		if(!strcmp(SM_basis->row_name[i],"E"))   Vchem->pos_exp   = length;
		if(!strcmp(SM_basis->row_name[i],"H"))   Vchem->pos_pH    = length;
		if(!strcmp(SM_basis->row_name[i],"X"))   Vchem->pos_X    = length;
		if(SM_basis->type[i]==1 && (Vchem->pos_exp != length)) Vchem->surface_flag = 1;

		length++;
	}
	/* init values */
	Vchem->log_a[Vchem->pos_pH]=Vchem->log_m[Vchem->pos_pH]=-7.;
	if(Vchem->pos_exp>0) Vchem->log_a[Vchem->pos_exp]=Vchem->log_m[Vchem->pos_exp]=0.;
	Vchem->log_a[Vchem->pos_water]=Vchem->log_m[Vchem->pos_water]=0.;
	Vchem->Io = 0.;
	Vchem->sch = 0.;
	Vchem->psi = 0.;
	Vchem->Temp = ICS->Temp;
	Vchem->chbal = 1; /* use charge balance - default */
	Vchem->nlinr = ICS->nlinr; /* non-linear rate equations? */

	for(i=0;i<ICS->size;++i)
	{
		pos = ICS->pos[i];
		Vchem->ctot[pos] = ICS->c_vchem[i];
		Vchem->log_a[pos] = Vchem->log_m[pos] = log10(ICS->c_vchem[i]);
		Vchem->ctot_ads[pos] = ICS->c_ads[i]; 
	}
	/* X-species is unphysical and the correct value is low */
	if(Vchem->pos_X > -1) Vchem->log_a[Vchem->pos_X ] = Vchem->log_m[Vchem->pos_X ] = -20.;
	/* buffers */
	length = 0;

	for(j=0; j<ICS->size;++j)
	{
		for(i=0; i<SM_mineral->size[0];++i)
		{
			if(!strcmp(SM_mineral->row_name[i],ICS->buffer_name[j]))
			{
				if(ICS->c_buffer[j] > 0.) /* not part */
				{
					Vchem->pos_buffer[length]=i;
					length++;
				}
			}
		}
	}

/* check if all the basis species are valid species */
	
	length_m = length_r = 0;
	for(j=0; j<ICS->size;++j)
	{
		for(i=0; i<Vchem->size;++i)
		{
			if(!strcmp(SM_basis->row_name[i],ICS->species_name[j]))
			{
				if(ICS->c_buffer[j] <= 0.) /* mass balance is used */
				{
					Vchem->pos_mass[length_m]=i;
					length_m++;
				} else /* basis switching */
				{
					Vchem->ctot_mineral[Vchem->pos_buffer[length_r]] = ICS->c_buffer[j];
					Vchem->pos_rock[length_r] = i;
					length_r++;
				}
			}
		}
	}
	Vchem->size_mass = length_m;
	Vchem->size_rock = length_r;

	
		/*CHECK*/
	/* for each rock spesies there has to be a corresponding buffer */
	if(length != Vchem->size_rock)
	{
		printf("Wrong dimensions on rock and buffer basis species\n");
		printf("Rock basis species: %d",   Vchem->size_rock);
		printf("Buffer basis species: %d", length);
		exit(1);
	}


	/* sup saturated minerals*/
	/* find number of minerals, then allocate space */
	/* Trick: find the minerals with least species and run through the list in that order */
	list    = (int *) malloc(Vchem->size*sizeof(int));
	list_b  = (int *) malloc(Vchem->size*sizeof(int));
	list[0]=Vchem->pos_exp;   /* E */ 
	list[1]=Vchem->pos_pH;    /* H+ */ 
	list[2]=Vchem->pos_water; /* H2O */ 
	list[3]=Vchem->pos_X;     /* X */ 
	/*list[4]=4;*/ /*HCO3*/
	size_list   = 4;
	size_list_b = 0;

	for(i=0;i<Vchem->size_rock;++i) /* do not count rock buffered species */
	{
		list[size_list]=Vchem->pos_rock[i];
		size_list++;
	}
	length = 0;
	while(length<ICS->size_sup_min)
	{
		length_i=1000;
		no_i=0;
		pos = -1;
		for(i=0;i<ICS->size_sup_min;++i)
		{
			if(0<=in_list(ICS->pos_sup_min[i], list_b,size_list_b))
				length_f=1000;
			else
				length_f=find_no_basis_species(ICS->pos_sup_min[i],&pos, list, size_list,ICS->SM_mineral);
			
			if(length_f<length_i && length_f >0)
			{
				pos_i = pos;
				no_i  = i;
				length_i=length_f;
				pos_b = ICS->pos_sup_min[i];
			}
		}
		if(ICS->PRINT_DEBUG_CHEM) printf("Number of basis species in mineral %s: %d choose basis species %s\n", 
				ICS->SM_mineral->row_name[ICS->pos_sup_min[no_i]], length_i, ICS->SM_mineral->col_name[pos_i]);
		ICS->pos_sup_bas[no_i] = pos_i; 
		list[size_list]=pos_i; /* remove chosen basis species from list */
		list_b[size_list_b]=pos_b; /* remove buffer mineral from list */
		length++;
		size_list++;
		size_list_b++;
	}

	for(i=0;i<ICS->size_sup_min;++i)
		if(ICS->PRINT_DEBUG_CHEM) printf("For mineral %s: choose basis species %s\n", 
		ICS->SM_mineral->row_name[ICS->pos_sup_min[i]], ICS->SM_mineral->col_name[ICS->pos_sup_bas[i]]);

	/* The following assumes that minerals are linearly independent */
	Vchem->size_sup_min = 0;
	for(i=0;i<Vchem->ICS->size_sup_min;++i)
	{
	  //printf("Rate-constants for %20s (LB units): k1 = %.4e, k2 = %.4e\n",
	  //    Vchem->ICS->SM_mineral->row_name[Vchem->ICS->pos_sup_min[i]], ICS->rate[i][0], ICS->rate[i][1]);
	  if(Vchem->ICS->c_sup_min[i] >0)
		{	
			Vchem->pos_rel_sup_min[Vchem->size_sup_min]     = i;
			Vchem->pos_sup_bas[Vchem->size_sup_min]         = Vchem->ICS->pos_sup_bas[i];
			Vchem->pos_sup_min[Vchem->size_sup_min]         = Vchem->ICS->pos_sup_min[i];
			Vchem->ctot_mineral[Vchem->ICS->pos_sup_min[i]] = Vchem->ICS->c_sup_min[i];
			if(ICS->PRINT_DEBUG_CHEM) printf("Sup min: %s with %s\n",
				Vchem->ICS->SM_basis->row_name[Vchem->pos_sup_bas[Vchem->size_sup_min]],Vchem->ICS->SM_mineral->row_name[Vchem->pos_sup_min[Vchem->size_sup_min]]);
			
			/* default values for generic form  rate = (k_1+k_2*a_h)*(1-SI^k_3)^k_4 */
			Vchem->rate[Vchem->size_sup_min][0] = ICS->rate[i][0];
			Vchem->rate[Vchem->size_sup_min][1] = ICS->rate[i][1];
			Vchem->rate[Vchem->size_sup_min][2] = ICS->rate[i][2];
			Vchem->rate[Vchem->size_sup_min][3] = ICS->rate[i][3];
			Vchem->size_sup_min++;
		}
		Vchem->ICS->SM_mineral->log_af[Vchem->ICS->pos_sup_min[i]] = Vchem->ICS->log_af_sup[i];
	}


	/* find reduced stoichiometric matrix �for sup sat minerals*/
	Vchem->sm_sup_buf = (real **) calloc(Vchem->size_sup_min,sizeof(real *));
	for(i=0; i< Vchem->size_sup_min; ++i) Vchem->sm_sup_buf[i] = (real *) calloc(Vchem->size_sup_min,sizeof(real ));

	update_sm_buf_new(Vchem,Vchem->size_sup_min, Vchem->pos_sup_min, Vchem->pos_sup_bas, Vchem->sm_sup_buf);
	if(ICS->PRINT_DEBUG_CHEM){
		for(i=0;i<Vchem->size_sup_min;++i)
		{
			for(j=0;j<Vchem->size_sup_min;++j)
				printf("%g\t",Vchem->sm_sup_buf[i][j]);
			printf("\n");
		}
	}


	
	/* basis transformation matrix */
	update_beta_inv(Vchem);
	
	/* find reduced stoichiometric matrix */
	Vchem->sm_buf = (real **) calloc(Vchem->size_rock,sizeof(real *));
	for(i=0; i< Vchem->size_rock; ++i) Vchem->sm_buf[i] = (real *) calloc(Vchem->size_rock,sizeof(real ));

	update_sm_buf_new(Vchem,Vchem->size_rock, Vchem->pos_buffer, Vchem->pos_rock, Vchem->sm_buf);
	
	/* HACK */
/*	c_dum[0]=-0.05;  */ /* Ca   */
/*	c_dum[1]=0.01;   */ /* Mg   */
/*	c_dum[2]=0.005;   */ /* SO4  */
/*	c_dum[3]=-1.5e-3; */ /* SiO2 */
/*	por = 0.26; 
	rho_r = 2710;
	for(i=0;i<Vchem->size_sup_min;++i)
	{
		c_dum_min[i] =0.;
		for(j=0;j<Vchem->size_sup_min;++j)
		{
			c_dum_min[i]+=Vchem->sm_sup_buf[i][j]*c_dum[j];
		}
		pos = Vchem->ICS->pos_sup_min[i];
		c_dum_min[i] *= Vchem->ICS->SM_mineral->mol_weight[pos]*por/(1.-por)/rho_r;
		printf("%s %g \%wt\n",ICS->SM_mineral->row_name[pos],c_dum_min[i]*100);
	}

*/
	free(list);
	free(list_b);
}

void free_BasVec(struct BasVec *Vchem)
{ 
	int i;

	free(Vchem->log_a);   
	free(Vchem->log_g);    
	free(Vchem->log_m);    
	free(Vchem->ctot );    
	free(Vchem->ctot_ads);     
	free(Vchem->ctot_dl  );   
	free(Vchem->ctot_ads_calc );    
	free(Vchem->ctot_calc );
	free(Vchem->ctot_mineral);
	free(Vchem->delta_mineral );
	free(Vchem->pos_rock);  
	free(Vchem->pos_buffer);
	free(Vchem->pos_mass);
	free(Vchem->pos_sup_min);
	free(Vchem->pos_sup_bas);
	free(Vchem->pos_rel_sup_min);

	
	//Vchem->beta_bas_inv   = (real **) calloc(Vchem->size,sizeof(real *));
	for(i=0;i<Vchem->size;++i) free(Vchem->beta_bas_inv[i]);
	free(Vchem->beta_bas_inv);

	for(i=0; i< Vchem->size_rock; ++i) free(Vchem->sm_buf[i]);

	free(Vchem->sm_buf);
	for(i=0;i<Vchem->ICS->size_sup_min;++i)	free(Vchem->rate[i]);
			
	for(i=0;i<Vchem->size_sup_min;++i) free(Vchem->sm_sup_buf[i]);
	free(Vchem->rate);
	free(Vchem->sm_sup_buf);
	

	/* more than one Vchem points to a single ICS, should be freed seperately */
	/* if(Vchem->ICS != NULL) freeICS(Vchem->ICS); */
	Vchem->ICS=NULL;
	
	/* struct RP */
	free(Vchem->RP->mol);
	free(Vchem->RP);

}
/* pos_buf = pos to buffer mineral */
/* returns number of basis species in mineral not counting elements in list*/
int find_no_basis_species(int pos_buf, int *pos_bas, int *list, int size_list, struct ChemTable *db)
{
	int i, N;

	N=0;

	for(i=0;i<db->size[1];++i)
	{
		if(db->M[pos_buf][i] != 0 &&  (0>in_list(i, list,size_list)) )
		{
			N++;
			*pos_bas = i;
		}
	}
	return N;
}


	
/*
init_ctot:
  initialise total concentrations of basis species, and basis buffers
  INPUT:  1) *c_rock        list of the aqueous concentration of species buffered by the rock
		  2) *c_buffer      list of the concentration of the minerals at the rock node
		                    first elements are for the rock buffers, and the last for the sup sat
							minerals
		  3) *c_mass        list of the aqueous concentration of species determined by mass balance
		  4) *Vchem
  OUTPUT: updated Vchem->ctot and Vchem->ctot_mineral
  RETURN: void
  */

void init_ctot(struct BasVec *Vchem)
{

	int i, pos;

	for(i=0;i<Vchem->ICS->size;++i)
	{
		pos = Vchem->ICS->pos[i];
		if(Vchem->ICS->c_vchem[i]<CHEM_LOWER_) Vchem->ctot[pos] = CHEM_LOWER_;
		else Vchem->ctot[pos] = Vchem->ICS->c_vchem[i];
		Vchem->log_a[pos] = Vchem->log_m[pos] = log10(Vchem->ICS->c_vchem[i]);

		Vchem->ctot_mineral[Vchem->pos_buffer[i]] = Vchem->ICS->c_buffer[i];
	}


/*
	write_BasVec_struct((*Vchem), "Vchem->out");
	*/
	return;
}
/* 
in_list:
   checks if an element is in the list

INPUT:  1) element : element in list
		2) list : list of integers
		3) size : of list

RETURN: position of element if element was found -1 otherwise
*/
			
int in_list(int element, int *list, int size_list)
{
	int flag=-1;
	int i;


	for(i=0;i<size_list;++i)
	{
		if(element==list[i])
		{
			flag=i;
			break;
		}
	}
	return flag;
}

/*
invert_matrix:
   calculates the inverse of a square matrix by the use of LU decomposition
   the original matrix is not destroyed in the process
INPUT:  1) a - the matrix to be inverted index 0,..,dim-1
        2) dim - the dimension of the matrix
OUTPUT: 1) a_inv - the inverse of matrix a
        2) d - the determinant of matrix a
		   NB: should not be used to calculate determinant of a large matrix
*/
void invert_matrix(real **a, int dim, real **a_inv, real *d)
{
	int i, j;
	int *indx;
	real *col, **a_tmp;

	a_tmp = (real **) calloc(dim+1, sizeof(real *));
	col   = (real *)  calloc(dim+1, sizeof(real ));
	indx  = (int *)   calloc(dim+1, sizeof(int));

	for(i=0;i<=dim;++i) a_tmp[i] = (real *) calloc(dim+1,sizeof(real));
	
	/* make a copy of the original matrix */
	for(i=0;i<dim;++i)
		for(j=0;j<dim;++j)
			a_tmp[i+1][j+1]=a[i][j];


	ludcmp(a_tmp,dim,indx,d);
	for(j=1;j<=dim;j++) (*d) *= a_tmp[j][j];
	if( fabs(*d) < CHEM_EPSILON_)
	{
		printf("singular matrix, no inverse found\n");
		exit(1);
	}

	for(j=1;j<=dim;j++) 
	{ 
		for(i=1;i<=dim;i++) col[i]=0.0;
		col[j]=1.0;
		lubksb(a_tmp,dim,indx,col);
		for(i=1;i<=dim;i++) a_inv[i-1][j-1]=col[i];
	}

	/* free space */
	free(col); 
	free(indx); 
	for(i=0;i<=dim;++i) free(a_tmp[i]);
	free(a_tmp); 
	return;
}

/*
in_fmatrix:
  from a matrix extract the rows and coloumns given in the list
  INPUT: 1) fmatrix: input matrix to be shaved 
         2) row: vector containing a list of integers, row number 2, 4, etc.
            the integers don't have to be sorted
         3) dim_row: dimension of the row vector
		 4) col: vector containing a list of integers, col number 1, 2, etc.
		 5) dim_col: dimension of the col vector
  OUTPUT: 1) f_red: the reduced matrix
  RETURN: void
  */
void in_fmatrix(real **fmatrix, int *row, int dim_row, int *col, int dim_col, real **f_red)
{
	int i, j;
	int r_element, c_element;

	for(i=0; i<dim_row; ++i)
	{
		r_element = row[i];
		for(j=0;j<dim_col;++j)
		{
			c_element = col[j];
			f_red[i][j]=fmatrix[r_element][c_element];
		}
	}
	return;
}

/*
update_sm_buf:
 from information in Vchem, find the new reduced stoichiometric matrix,
 this matrix is transposed and inverted
 INPUT: Vchem
 OUTPUT: sm_buf matrix in Vchem is updated
 RETURN: void
*/
void update_sm_buf(struct BasVec *Vchem)
{
	real **m_dum,d;
	int i,j;
	struct ChemTable *SM_all, *SM_mineral;

	SM_all     = Vchem->ICS->SM_all;
	SM_mineral = Vchem->ICS->SM_mineral;

	m_dum = (real **) calloc(Vchem->size_rock, sizeof(real *));
	for(i=0; i<Vchem->size_rock; ++i) m_dum[i] = (real *) calloc(Vchem->size_rock, sizeof(real));

	in_fmatrix(SM_mineral->M, Vchem->pos_buffer, Vchem->size_rock, Vchem->pos_rock, Vchem->size_rock, Vchem->sm_buf);
	/* we need the inverse of the transposed reduced stoichiometric matrix */
	for(i=0; i<Vchem->size_rock;++i)
	{
		for(j=0; j<Vchem->size_rock;++j)
		{
			m_dum[i][j] = Vchem->sm_buf[j][i];
		}
	}
	invert_matrix(m_dum,Vchem->size_rock, Vchem->sm_buf,&d);
	
	for(i=0;i<Vchem->size_rock;++i) free(m_dum[i]);
	free(m_dum);
	return;
}

/*
update_sm_buf_new:
 from information in Vchem, find the new reduced stoichiometric matrix,
 this matrix is transposed and inverted
 INPUT: Vchem
 OUTPUT: sm_buf matrix in Vchem is updated
 RETURN: void
*/
void update_sm_buf_new(struct BasVec *Vchem, int size_r, int *pos_b, int *pos_r, real **sm)
{
	real **m_dum,d;
	int i,j;
	struct ChemTable *SM_all, *SM_mineral;

	SM_all     = Vchem->ICS->SM_all;
	SM_mineral = Vchem->ICS->SM_mineral;

	m_dum = (real **) calloc(size_r, sizeof(real *));
	for(i=0; i<size_r; ++i) m_dum[i] = (real *) calloc(size_r, sizeof(real));

	in_fmatrix(SM_mineral->M, pos_b, size_r, pos_r, size_r, sm);
	/* we need the inverse of the transposed reduced stoichiometric matrix */
	for(i=0; i<size_r;++i)
	{
		for(j=0; j<size_r;++j)
		{
			m_dum[i][j] = sm[j][i];
			if(Vchem->ICS->PRINT_DEBUG_CHEM) printf("%g\t",m_dum[i][j]);
		}
		//if(Vchem->ICS->PRINT_DEBUG_CHEM) printf("%\n");
	}
	invert_matrix(m_dum,size_r, sm,&d);
	
	for(i=0;i<size_r;++i) free(m_dum[i]);
	free(m_dum);
	return;
}

/* 
update_beta_inv:
  find the basis transformation matrix, calculates its inverse and store it 
  in Vchem->beta_inv
  INPUT: *Vchem
  OUTPUT: updatet Vchem->beta_inv
  RETURN: void
*/
void update_beta_inv(struct BasVec *Vchem) 
{
	int i, j, length;
	real **beta_bas, d;
	//FILE *fp;
	struct ChemTable *SM_all, *SM_mineral;

	SM_all     = Vchem->ICS->SM_all;
	SM_mineral = Vchem->ICS->SM_mineral;

	beta_bas        = (real **) calloc(Vchem->size,sizeof(real *));
	for(i=0;i<Vchem->size; ++i)
		beta_bas[i] = (real *) calloc(Vchem->size, sizeof(real ));
	/* basis transformation matrix */
	for(j=0;j<Vchem->size;++j)
	{
		/* check if basis species j is one of the rock buffered species */
		length=in_list(j,Vchem->pos_rock,Vchem->size_rock);
		if( length >= 0)
		{
			for(i=0;i<Vchem->size;++i) beta_bas[j][i]=SM_mineral->M[Vchem->pos_buffer[length]][i];
		} else
		{
			for(i=0;i<Vchem->size;++i)
			{
				if(i==j) 
					beta_bas[j][i]=1.; 
				else 
					beta_bas[j][i]=0.;
			}
		}
	}

	invert_matrix(beta_bas,Vchem->size, Vchem->beta_bas_inv, &d);
	
	/*DEBUG*/
	/*
	fp = my_fopen("debug.out","w");
	for(j=0;j<Vchem->size; ++j)
	{
		for(i=0;i<Vchem->size;++i)
		{
			fprintf(fp,"%lf\t", beta_bas[j][i]);

		}
		fprintf(fp,"\n");
	}
	fprintf(fp,"\n");

	for(j=0;j<Vchem->size; ++j)
	{
		for(i=0;i<Vchem->size;++i)
		{
			fprintf(fp,"%lf\t", Vchem->beta_bas_inv[j][i]);

		}
		fprintf(fp,"\n");
	}

	fprintf(fp,"\n");

	fclose(fp);
*/
	for(i=0;i<Vchem->size;++i) free(beta_bas[i]);
	free(beta_bas);
return;
}
/*
calc_rock_spec_conc:
  calculates the concentration of the rock buffered species, from the
  basis transformation matrix, i.e.:
  log_10 a = beta^-1 *( log_10 a' + log_10 K)
INPUT: Vchem
OUTPUT: concentrations of the rock buffered species, stored in Vchem.log_a
RETURN: void
*/

void calc_rock_spec_conc(struct BasVec *Vchem)
{
	real *log_bas;
	int i, j, pos, posb;

	log_bas  = (real *) calloc(Vchem->size, sizeof(real));

	for(i=0;i<Vchem->size;++i) log_bas[i]=(Vchem->log_a[i]+ Vchem->ICS->SM_basis->logK[i]); /* initialize */

	/* basis species -> basis mineral */
	for(i=0; i< Vchem->size_rock; ++i)
	{
		pos  = Vchem->pos_rock[i];
		posb = Vchem->pos_buffer[i]; 
		log_bas[pos] = Vchem->ICS->SM_mineral->log_af[posb] + Vchem->ICS->SM_mineral->logK[posb]; /* log activity of basis mineral */
	}

	/* do the matrix multiplication */

	for(i=0; i< Vchem->size_rock; ++i)
	{
		pos=Vchem->pos_rock[i];
		Vchem->log_a[pos] = 0.;
		for(j=0; j<Vchem->size; ++j) Vchem->log_a[pos] += Vchem->beta_bas_inv[pos][j]*log_bas[j];

		Vchem->log_m[pos] = Vchem->log_a[pos] - Vchem->log_g[pos];
	}

	free(log_bas); log_bas = NULL;

	return;	
}


/*
chem_string_to_pos:
 finds the positons of the row names in list_t sm in table, and returns the positions
 INPUT: 1) ChemTable *sm: table
		2) list_t: a subset of row names in table
 OUTPUT: 1) pos: the position of the names in list_t relative to table sm
 RETUNR: size of list
*/

//int chem_string_to_pos(struct ChemTable *sm, char *list_t, int **pos)
//{
//	int length, i;
//	char list[100];
//	char ct[] = "\n \t \0";/* MAC/LINUX USE 0 */
//	char *s_tmp;
//
//	/* make copy - strtok destroys original string */
//	if(sprintf(list, "%s", list_t)<0)
//	{
//		printf("allocate more char space in routine : chem_string_to_pos \n");
//		exit(1);
//	}
//	/* species determined by mass balance */
//	s_tmp = strtok(list,ct);
//	length=0;
//	while(s_tmp != NULL)
//	{
//		for(i=0; i<sm->size[0];++i)
//		{
//			if(!strcmp(sm->row_name[i],s_tmp)) length++;
//		}
//		s_tmp = strtok(NULL,ct);
//	}
//	(*pos)  = (int *) calloc(length,sizeof(int));
//
//	sprintf(list, "%s", list_t);
//	s_tmp = strtok(list,ct);
//	length=0;
//	while(s_tmp != NULL)
//	{
//		for(i=0; i<sm->size[0];++i)
//		{
//			if(!strcmp(sm->row_name[i],s_tmp))
//			{
//				(*pos)[length]=i;
//				length++;
//			}
//		}
//		s_tmp = strtok(NULL,ct);
//	}
//
//	return length;
//}




void db_set_activity(struct InitChem *ICS)
{
	int i;
	for(i=0;i<ICS->size_rock;++i)
	{
		ICS->SM_mineral->log_af[ICS->pos_buffer[i]] = ICS->log_af[ICS->pos_rock[i]];
	}

	for(i=0;i<ICS->size_sup_min;++i)
	{
		ICS->SM_mineral->log_af[ICS->pos_sup_min[i]] = ICS->log_af_sup[i];
	}
}

void trans_nlin_eq(struct BasVec *Vchem)
{
	int i;

	for(i=0;i<Vchem->size_sup_min;++i)
		add_new_buffer_mineral(Vchem->pos_sup_bas[i], Vchem->pos_sup_min[i], Vchem);

	
	if(Vchem->ICS->PRINT_DEBUG_CHEM)
		write_BasVec_struct(Vchem,"qqqV2.out");

	for(i=0;i<Vchem->size_sup_min;++i) free(Vchem->sm_sup_buf[i]);
	Vchem->size_sup_min = 0;

}
