#include "chem_global.h"

/* reads a database, the first coloumn and row are names */
/* allocate space and returns a matrix M of everything but */
/* the first coloumn and row. The first row and coloumn are */
/* returned in row_name and col_name, size returns the dimension */
/* of the matrix */
void read_database(std::string &folder, const std::string &filename, struct ChemTable *tab)
{
	FILE *fp;
	int i, j,k,Nblanks;
	char ct[] = "\n \t \r \0 \v"; /* tokens used in reading a list of numbers from a string */
	char *s_tmp,*fc; /* temporary string */
	int maxline = 400, test=GOOD; /* maximum line length */
	char line[400];  /* holds an uncommented line from the inputfile */
	char *dum_char;

	/* Open input file */
	std::string datafile(folder + filename);
	fp = my_fopen(datafile.c_str(), "r");
	
	fc=myreadline(line, maxline, fp);	
	s_tmp = strtok(line,ct);
	j=0;
	while( s_tmp != NULL)
	{
		s_tmp =strtok(NULL,ct);
		++j;
		/*printf("%s\n", s_tmp);*/
	}
	
	//printf("No col = %d\n", j);
	 
	i=0;
	
	while(!(fc == NULL))
	{
		fc=myreadline(line, maxline, fp);
	    ++i; 
	}

	//printf("No rows = %d\n",i);

	tab->size[0]=i-1;
	tab->size[1]=j-1;
/* allocate memory */

	/* allocate row names*/
  if(!(tab->row_name =(char **) calloc(tab->size[0],sizeof(char *)))) test=BAD;
  for(i=0;i<tab->size[0];++i) if(!(tab->row_name[i] =(char *) calloc(COMPLEXNAME,sizeof(char )))) test=BAD;
	/* allocate col names*/
  if(!(tab->col_name =(char **) calloc(tab->size[1],sizeof(char *)))) test=BAD;
  for(j=0;j<tab->size[1];++j) if(!(tab->col_name[j] =(char *) calloc(COMPLEXNAME,sizeof(char )))) test=BAD;

	/* allocate stoichiometric matrix */

  if(!(tab->M =(real **) calloc(tab->size[0],sizeof(real *)))) test=BAD;
  for(i=0;i<tab->size[0];++i) if(!(tab->M[i] =(real *) calloc(tab->size[1],sizeof(real )))) test=BAD;

  

	if(test==BAD)
    {
      printf("OUT OF MEMORY\n");
      exit(1);
    }

	/* read file again */
	rewind(fp);


	fc=myreadline(line, maxline, fp);	
	s_tmp = strtok(line,ct);

	for(j=0;j<tab->size[1]; ++j){
		dum_char=strtok(NULL,ct);
		/* find number of blanks at the end */
		Nblanks=strlen(dum_char);
			for(i=0;i<strlen(dum_char);++i)
		{
			if(dum_char[strlen(dum_char)-1-i]== ' ' || dum_char[strlen(dum_char)-1-i]== '\n')
			{
				Nblanks--;
			} else break; 
		}
			strncat(tab->col_name[j],dum_char,Nblanks);
	}

	for(i=0;i<tab->size[0];++i)
	{
		myreadline(line, maxline, fp);
		s_tmp=strtok(line,ct);
		Nblanks=strlen(s_tmp);
			for(k=0;k<strlen(s_tmp);++k)
		{
			if(s_tmp[strlen(s_tmp)-1-k]== ' ' || s_tmp[strlen(dum_char)-1-k]== '\n')
			{
				Nblanks--;
			} else break; 
		}
			strncat(tab->row_name[i],s_tmp,Nblanks);		
		for(j=0;j<tab->size[1]; ++j)
		{
			s_tmp=strtok(NULL,ct);
			sscanf(s_tmp,"%lf",&(tab->M[i][j]));
			
		}
	}
	
	fclose(fp);
	tab->a0      = NULL;
	tab->logK    = NULL;
	tab->charge  = NULL;
	tab->scharge = NULL;
	tab->log_m   = NULL;
	tab->log_a   = NULL;
	tab->log_g   = NULL;
	tab->log_af  = NULL;
	tab->delta   = NULL;
	tab->log_QK  = NULL;
	tab->mol_volume = NULL;
	tab->mol_weight = NULL;
	tab->abs_pos_bas = NULL;
	tab->abs_pos = NULL;
	tab->type = NULL;
	tab->T = NULL;
}
	


/* shaves a table of elements that are not needed in the calculations, based on         */
/* a row vector, if one of the row vectors has a 1 in a place where (vec) has a value 0  */
/* it is to be taken out of the calculation. Start row and start col must be given, to  */
/* allow for the possibility to store more information in the M matrix                  */
/* In the exmaple below start_row = 1 and start_col=1 if all the M values is to be      */
/* considered.                                                                          */ 
/*            <col_name> <col_name> <col_name> <col_name> <col_name> <col_name>         */
/* <row_name>  M00        M01        M02         M03       M04         M05              */
/* <row_name>  M10        M11        M12         M13       M14         M15              */
/* <row_name> etc ...                                                                   */
 

void remove_element(struct ChemTable *tab, real *vec, int start_row, int start_col, struct ChemTable *CM)
{
  int N_row, N_col,i,j,Npos_i, Npos_j;
  int test=GOOD, flag=1;
  

  /* number of  coloumns */
  N_col = tab->size[1];
  N_col -= (start_col-1);
  for(i=0;i<tab->size[1]-start_col+1;++i) if((int) vec[i] == 0) N_col--;
    
  /* remove complexes that is not needed */
  /* nuber of rows */
	
  N_row = tab->size[0];
  N_row -= (start_row-1);

  for(i=0; i < tab->size[0]-start_row+1; ++i)
    {
      for(j=0; j< tab->size[1]-start_col+1 ; ++j)
	{
	  if( (int) vec[j] == 0 && ((int) tab->M[start_row+i-1][start_col+j-1] )!= 0) 
	    {
	      N_row --;
	      break;
	    }
	}
    }
  
  CM->size[0]=N_row;
  CM->size[1]=N_col;
  /* allocate sufficient space */

  /*row names */
  if(!(CM->row_name =(char **) calloc(N_row,sizeof(char *)))) test=BAD;
  for(i=0;i<CM->size[0];++i) if(!(CM->row_name[i] =(char *) calloc(COMPLEXNAME,sizeof(char )))) test=BAD;
  /*col names */
  if(!(CM->col_name =(char **) calloc(N_col,sizeof(char *)))) test=BAD;
  for(j=0;j<CM->size[1];++j) if(!(CM->col_name[j] =(char *) calloc(COMPLEXNAME,sizeof(char )))) test=BAD;
  /* position */
  if(!(CM->abs_pos     = (int *) calloc(N_row,sizeof(int )))) test=BAD;
  if(!(CM->abs_pos_bas = (int *) calloc(N_col,sizeof(int )))) test=BAD;
  /*matrix */
  if(!(CM->M =(real **) calloc(CM->size[0],sizeof(real *)))) test=BAD;
  for(i=0;i<CM->size[0];++i) if(!(CM->M[i] =(real *) calloc(CM->size[1],sizeof(real )))) test=BAD;
  /* concentration of basis specise */
  if(!(CM->log_m =(real *) calloc(CM->size[0],sizeof(real )))) test=BAD;
  if(!(CM->log_a =(real *) calloc(CM->size[0],sizeof(real )))) test=BAD;
  if(!(CM->log_g =(real *) calloc(CM->size[0],sizeof(real )))) test=BAD;
  if(!(CM->log_QK =(real *) calloc(CM->size[0],sizeof(real )))) test=BAD;
  if(!(CM->delta =(real *) calloc(CM->size[0],sizeof(real )))) test=BAD;
  
  CM->logK       = NULL;
  CM->a0         = NULL;
  CM->scharge    = NULL;
  CM->charge     = NULL;
  CM->mol_volume = NULL;
  CM->mol_weight = NULL;
  CM->type       = NULL;
  CM->log_af     = NULL;
  CM->T          = NULL;

  /* initialize concentrations */
  for(i=0;i<CM->size[0];++i)
  { 
	  CM->log_m[i]=CM->log_a[i]=CM->log_g[i]=CM->log_QK[i]= CM->delta[i]=0.;
  }
  
  	if(test==BAD)
    {
      printf("OUT OF MEMORY\n");
      exit(1);
    }

  Npos_j=0;

 for(j=0; j <tab->size[1]-start_col+1 ; ++j)
 {
	 if((int) vec[j] == 1 )
	 {
		 CM->col_name[Npos_j]=tab->col_name[start_col+j-1];
		 CM->abs_pos_bas[Npos_j]=j;
		 Npos_j++;
	 }
 }
Npos_j=Npos_i=0;

 for(i=0; i < tab->size[0]-start_row+1; ++i)
 {
	 Npos_j=0;
	 
	 for(j=0; j< tab->size[1]-start_col+1 ; ++j)
	{
		if( (int) vec[j] == 0 && ((int) tab->M[start_row+i-1][start_col+j-1] )!= 0) 
		{
			Npos_i--;
			flag=0;
			break;
		} 
	 }
	 
	 if(flag) 
	 {/* the complex is in the calculation */
		 
		 for(j=0; j< tab->size[1]-start_col+1 ; ++j)
		 
		 {
			 if( (int) vec[j] == 1 )
			 { 
				 CM->abs_pos[Npos_i]=i;
				 CM->M[Npos_i][Npos_j]=tab->M[start_row+i-1][start_col+j-1];
				 CM->row_name[Npos_i]=tab->row_name[start_row+i-1];
				 Npos_j++;
			 }
		 }
	 }
	  flag=1;
	  Npos_i++;
 }
 
 /* if set, update other elements in the struct */

 copy_ChemTable_vector(&(CM->logK),CM, tab->logK,0);
  if(CM->logK == NULL) CM->logK =(real *) calloc(CM->size[0],sizeof(real ));
  copy_ChemTable_vector(&(CM->a0),CM, tab->a0,0);
  copy_ChemTable_vector(&(CM->charge),CM, tab->charge,0);
  copy_ChemTable_vector(&(CM->scharge),CM, tab->scharge,0);
  copy_ChemTable_vector(&(CM->mol_volume),CM, tab->mol_volume,0);
  copy_ChemTable_vector(&(CM->mol_weight),CM, tab->mol_weight,0);
  copy_ChemTable_vector(&(CM->log_af),CM, tab->log_af,0); 
  icopy_ChemTable_vector(&(CM->type) ,CM, tab->type,0); 
  copy_ChemTable_vector(&(CM->delta),CM, tab->delta,0); 
  copy_ChemTable_vector(&(CM->T),CM, tab->T,1); 
}
void remove_element_rc (struct ChemTable *tab, real *vec_row, real *vec_col, int start_row, int start_col, struct ChemTable *CM)
{
  int N_row, N_col,i,j,Npos_i, Npos_j;
  int test=GOOD, flag=0, flag_cc=0, flag_cr=0;
  
  
  if( vec_row == NULL )
  { /* assume all rows are kept */
	  vec_row = (real *) calloc(tab->size[0], sizeof(real) );
	  for(i=0;i<tab->size[0];++i) vec_row[i]=1.;
	  flag_cr = 1;
  }

  if(vec_col == NULL)
  { /* assume all col are kept */
	  vec_col = (real *) calloc(tab->size[1], sizeof(real) );
	  for(i=0;i<tab->size[1];++i) vec_col[i]=1.;
	  flag_cc=1;
  }


  /* number of  coloumns */
  N_col = tab->size[1];
  N_col -= (start_col-1);
  
  for(i=0;i<tab->size[1]-start_col+1;++i) if((int) vec_col[i] == 0) N_col--;
    
  /* remove rows that is not needed */
	
  N_row = tab->size[0];
  N_row -= (start_row-1);

  for(i=0; i < tab->size[0]-start_row+1; ++i) if((int) vec_row[i] == 0) N_row--;
  
  if(N_row != 0 || N_col != 0)
  {
  CM->size[0]=N_row;
  CM->size[1]=N_col;
  

  /* allocate sufficient space */

  /*row names */
  if(!(CM->row_name =(char **) calloc(N_row,sizeof(char *)))) test=BAD;
  for(i=0;i<CM->size[0];++i) if(!(CM->row_name[i] =(char *) calloc(COMPLEXNAME,sizeof(char )))) test=BAD;
  /*col names */
  if(!(CM->col_name =(char **) calloc(N_col,sizeof(char *)))) test=BAD;
  for(j=0;j<CM->size[1];++j) if(!(CM->col_name[j] =(char *) calloc(COMPLEXNAME,sizeof(char )))) test=BAD;
  /* position */
  if(!(CM->abs_pos     = (int *) calloc(N_row,sizeof(int )))) test=BAD;
  if(!(CM->abs_pos_bas = (int *) calloc(N_col,sizeof(int )))) test=BAD;
  /*matrix */
  if(!(CM->M =(real **) calloc(CM->size[0],sizeof(real *)))) test=BAD;
  for(i=0;i<CM->size[0];++i) if(!(CM->M[i] =(real *) calloc(CM->size[1],sizeof(real )))) test=BAD;
  /* concentration of basis specise */
  if(!(CM->log_m =(real *) calloc(CM->size[0],sizeof(real )))) test=BAD;
  if(!(CM->log_a =(real *) calloc(CM->size[0],sizeof(real )))) test=BAD;
  if(!(CM->log_g =(real *) calloc(CM->size[0],sizeof(real )))) test=BAD;
  if(!(CM->log_QK =(real *) calloc(CM->size[0],sizeof(real )))) test=BAD;
  if(!(CM->delta =(real *) calloc(CM->size[0],sizeof(real )))) test=BAD;
  CM->logK       = NULL;
  CM->a0         = NULL;
  CM->scharge    = NULL;
  CM->charge     = NULL;
  CM->mol_volume = NULL;
  CM->mol_weight = NULL;
  CM->type       = NULL;
  CM->log_af     = NULL;
  CM->T          = NULL;
  /* initialize concentrations */

  for(i=0;i<CM->size[0];++i)
  { 
	  CM->log_m[i]=CM->log_a[i]=CM->log_g[i]=CM->log_QK[i]=CM->delta[i]= 0.;

  }
  
  
 
  	if(test==BAD)
    {
      printf("OUT OF MEMORY\n");
      exit(1);
    }

  Npos_j=0;
 for(j=0; j <tab->size[1]-start_col+1 ; ++j)
 {
	 if((int) vec_col[j] == 1 )
	 {
		 CM->col_name[Npos_j]=tab->col_name[start_col+j-1];
		 CM->abs_pos_bas[Npos_j]=j;
		 Npos_j++;
	 }
 }
 Npos_j=Npos_i=0;

 for(i=0; i < tab->size[0]-start_row+1; ++i)
 {
	 Npos_j=0;
	 for(j=0; j< tab->size[1]-start_col+1 ; ++j)
	{
		 
		 if( ((int) vec_col[j] == 1) && ((int) vec_row[i] == 1 ) )
		 {/* add coloumn and row element */
			 CM->abs_pos[Npos_i]=i;
			 CM->M[Npos_i][Npos_j]=tab->M[start_row+i-1][start_col+j-1];
			 CM->row_name[Npos_i]=tab->row_name[start_row+i-1];
			 Npos_j++;
			 flag=1;
	    } 
	 }
	 if(flag) Npos_i++;
	 flag=0;
    }
  }
 /* if set, update other elements in the struct */
  copy_ChemTable_vector(&(CM->logK),CM, tab->logK,0);
  if(CM->logK == NULL) CM->logK =(real *) calloc(CM->size[0],sizeof(real ));
  copy_ChemTable_vector(&(CM->a0),CM, tab->a0,0);
  copy_ChemTable_vector(&(CM->charge),CM, tab->charge,0);
  copy_ChemTable_vector(&(CM->scharge),CM, tab->scharge,0);
  copy_ChemTable_vector(&(CM->mol_volume),CM, tab->mol_volume,0);
  copy_ChemTable_vector(&(CM->mol_weight),CM, tab->mol_weight,0);
  copy_ChemTable_vector(&(CM->log_af),CM, tab->log_af,0); 
  icopy_ChemTable_vector(&(CM->type) ,CM, tab->type,0); 
  copy_ChemTable_vector(&(CM->delta),CM, tab->delta,0);
  copy_ChemTable_vector(&(CM->T),CM, tab->T,1); 
 
	 
 
 /* remember to free space */
 if(flag_cc) free(vec_col);
 if(flag_cr) free(vec_row);

}



char *myreadline(char *line, int maxline, FILE *input_file)
/* reads the first line without '#' or '\n' */
{
  char *ff;
  int i;
  for( i=0; i<maxline; ++i)  line[i]= ' ';

  ff = fgets(line, maxline, input_file);
  while(!(ff == NULL) ){
    if( line[0] != '#' && line[0] != '\n' && line[0] != '\r' )
    {
      return ff;
      break;
    }	
    for( i=0; i<maxline; ++i)  line[i]= ' ';
    ff = fgets(line, maxline, input_file);
  }
  return ff; 
}

real *pick_col(struct ChemTable  *tab, const char *name)
{
  int j,i, flag=1;
  int test =GOOD;
  real *col=NULL;

	for(j = 0; j<tab->size[1]; ++j)
	{
		if (!strcmp(tab->col_name[j],name) )
		{
			if(!(col =(real *) calloc(tab->size[0],sizeof(real )))) test=BAD;
		  	if(test==BAD)
			{
				printf("OUT OF MEMORY\n");
				exit(1);
			}
		
		  for(i=0;i<tab->size[0];++i) col[i]= (real) tab->M[i][j];
		  
		  return col;
		  flag=0;
		  break;
		}
	}

	if (flag){
	  printf("WARNING: No coloumn named %s\n",name);
	}
return NULL;
	
}

int *ipick_col(struct ChemTable  *tab, const char *name)
{
  int j,i, flag=1;
  int test =GOOD;
  int *col=NULL;

	for(j = 0; j<tab->size[1]; ++j)
	{
		if (!strcmp(tab->col_name[j],name) )
		{
			if(!(col =(int *) calloc(tab->size[0],sizeof(int )))) test=BAD;
		  	if(test==BAD)
			{
				printf("OUT OF MEMORY\n");
				exit(1);
			}
		
		  for(i=0;i<tab->size[0];++i) col[i]= (int) tab->M[i][j];
		  
		  return col;
		  flag=0;
		  break;
		}
	}

	if (flag){
	  printf("WARNING: No coloumn named %s\n",name);
	}
return NULL;
	
}

/* copy_ChemTable_vector :
copy vector from a struct; n = 0 coloumn vector n = 1 row vector 
*/
void copy_ChemTable_vector(real **vec_out, struct ChemTable *CM, real *vec, int n)
{
	int i, *pos;
	pos = ( n==0 ? (CM->abs_pos) : (CM->abs_pos_bas));
	if(vec != NULL)
	{
		if( (*vec_out) == NULL) (*vec_out)     = (real *) calloc(CM->size[n],sizeof(real ));
		for(i=0;i<CM->size[n];++i) (*vec_out)[i]    = vec[pos[i]];
	} 
}

/* copy_ChemTable_vector :
copy vector from a struct; n = 0 coloumn vector n = 1 row vector 
*/
void icopy_ChemTable_vector(int **vec_out, struct ChemTable *CM, int *vec, int n)
{
	int i, *pos;
	pos = ( n == 0 ? (CM->abs_pos) : (CM->abs_pos_bas));

	if(vec != NULL)
	{
		if((*vec_out) == NULL) (*vec_out)     = (int *) calloc(CM->size[n],sizeof(int ));
		for(i=0;i<CM->size[n];++i) (*vec_out)[i]    = vec[pos[i]];
	} 
}

real *pick_row(struct ChemTable  *tab, char *name)
{
	int i, flag=1;
	
	for(i = 0; i<tab->size[0]; ++i)
	{
		if (!strcmp(tab->row_name[i],name) )
		{
			return tab->M[i];
			flag=0;
			break;
		}
	}

	if (flag){
		printf("WARNING : No row named %s\n",name);
		
	}

	return NULL;
}

void free_db(struct ChemTable *tab)
{
	int i;

	if(tab->a0 != NULL) free(tab->a0);
	if(tab->abs_pos != NULL) free(tab->abs_pos);
	if(tab->abs_pos_bas != NULL) free(tab->abs_pos_bas);
	if(tab->charge != NULL) free(tab->charge);
	if(tab->scharge != NULL) free(tab->scharge);
	if(tab->col_name != NULL) free(tab->col_name);
	if(tab->row_name != NULL) free(tab->row_name);
	if(tab->delta != NULL) free(tab->delta);
	if(tab->log_a != NULL) free(tab->log_a);
	if(tab->log_af != NULL) free(tab->log_af);
	if(tab->log_m != NULL) free(tab->log_m);
	if(tab->log_g != NULL) free(tab->log_g);
	if(tab->log_QK != NULL) free(tab->log_QK);
	if(tab->logK != NULL) free(tab->logK);
	if(tab->mol_volume != NULL) free(tab->mol_volume);
	if(tab->type != NULL) free(tab->type);
	if(tab->M != NULL)
	{
		for(i=0;i<tab->size[0];++i) free(tab->M[i]);
		free(tab->M);
	}
}



	


