#include "chem_global.h"

void init_chemistry(struct InitChem *ICS, std::string &data_folder)
{
	
	int i;
	char s_bas_tmp[200];
	char *s_tmp;
	char ct[] = "\n \t \0";/* MAC LINUX 0*/
	struct ChemTable *Mineral_db, *Aq_db, *Basis_db, *Aq_logK_db, *M_logK_db;
/*
	struct InitChem ICS_new;
	char name_tt[100];
*/
	real *bs;
	ICS->PRINT_DEBUG_CHEM = 0; /* 1=debug info 0=no output */
	

	/*  memory leak ............... */
	

	
	init_const_chem_param(ICS->CP); /* set constant parameters k_B, Na, etc.*/

	/* ...........................*/

	Mineral_db = (struct ChemTable *) malloc( 1 * sizeof(struct ChemTable )); 
	Aq_db      = (struct ChemTable *) malloc( 1 * sizeof(struct ChemTable ));
	Basis_db   = (struct ChemTable *) malloc( 1 * sizeof(struct ChemTable ));
	Aq_logK_db = (struct ChemTable *) malloc( 1 * sizeof(struct ChemTable ));
	M_logK_db  = (struct ChemTable *) malloc( 1 * sizeof(struct ChemTable ));


	/* read databases */
	read_database(data_folder, "/aq_all.dat",       Aq_db);
	read_database(data_folder, "/buffer_all.dat",   Mineral_db);
	read_database(data_folder, "/basis_large.dat",  Basis_db);
	read_database(data_folder, "/aq_logKT.dat",     Aq_logK_db);
	read_database(data_folder, "/buffer_logKT.dat", M_logK_db);

	set_chem_param(Aq_db);
	set_chem_param(Basis_db);

    Mineral_db->log_af     = pick_col(Mineral_db,"logaf");
	Mineral_db->mol_volume = pick_col(Mineral_db,"mol_volume");
	Aq_db->type            = ipick_col(Aq_db,"type");
	Basis_db->mol_weight   = pick_col(Basis_db,"mol_weight");
	Basis_db->type         = ipick_col(Basis_db,"type");
	/* ----------DEBUG-------------------------------------------------*/
	if(ICS->PRINT_DEBUG_CHEM) write_chem_struct(Aq_db, "Aq_db.out");
	if(ICS->PRINT_DEBUG_CHEM) write_chem_struct(Basis_db, "Basis_db.out");
	if(ICS->PRINT_DEBUG_CHEM) write_chem_struct(Mineral_db, "Mineral_db.out");
	if(ICS->PRINT_DEBUG_CHEM) write_chem_struct(Aq_logK_db, "Aq_logK_db.out");
	if(ICS->PRINT_DEBUG_CHEM) write_chem_struct(M_logK_db, "M_logK_db.out");
	/*-----------------------------------------------------------------*/
	/* make a logical vector,  bs has 1 if the basis  */
	/* species is in the calculation and zero othervise */
	bs = (real *) calloc(Basis_db->size[0],sizeof(real));

	for(i=0;i<Basis_db->size[0];++i) bs[i]=0.;

	sprintf(s_bas_tmp," H E H2O ");
	for(i=0;i<ICS->size;++i)
	{
		strcat(s_bas_tmp, ICS->species_name[i]);
		strcat(s_bas_tmp, " ");
	}

	s_tmp = strtok(s_bas_tmp,ct);
	while(s_tmp != NULL)
	{
		for(i=0; i<Basis_db->size[0];++i)
		{
			if(!strcmp(Basis_db->row_name[i],s_tmp))
			{
				bs[i]=1.;
			}
		}
			
		s_tmp = strtok(NULL,ct);
	}


	/* make reduced databases */
		remove_element(Aq_db, bs,1,5,ICS->SM_all);
		remove_element(Mineral_db,bs,1,4,ICS->SM_mineral);
		remove_element_rc(Basis_db, bs,NULL,1,1,ICS->SM_basis);
		free(bs);
		bs  = (real *) calloc(Aq_db->size[0], sizeof(real));
		for(i=0;i<Aq_db->size[0];++i) bs[i] = 0.;
		/* 1 for all elements in SM_all_ */
		for(i=0;i<ICS->SM_all->size[0];++i) bs[ICS->SM_all->abs_pos[i]] = 1.;

		remove_element_rc(Aq_logK_db,bs,NULL,1,1,ICS->SM_aq_logK); 


		free(bs);
		bs  = (real *) calloc(Mineral_db->size[0], sizeof(real));
		for(i=0;i<Mineral_db->size[0];++i) bs[i] = 0.;
		/* 1 for all elements in SM_all_ */
		for(i=0;i<ICS->SM_mineral->size[0];++i) bs[ICS->SM_mineral->abs_pos[i]] = 1.;

		remove_element_rc(M_logK_db,bs,NULL,1,1,ICS->SM_M_logK);


	    /* check if something is wrong in the database */
		compare_names(ICS->SM_aq_logK,ICS->SM_all); 
		compare_names(ICS->SM_M_logK,ICS->SM_mineral); 

		/* get temperature from col names */
		ICS->SM_aq_logK->T = (real *) calloc(ICS->SM_aq_logK->size[1],sizeof(real));
        ICS->SM_M_logK->T = (real *) calloc(ICS->SM_M_logK->size[1],sizeof(real));
		for(i=0;i<ICS->SM_aq_logK->size[1];++i) ICS->SM_aq_logK->T[i] = ((real) get_int_from_string(ICS->SM_aq_logK->col_name[i]))+273.15;
		for(i=0;i<ICS->SM_M_logK->size[1];++i)  ICS->SM_M_logK->T[i] = ((real) get_int_from_string(ICS->SM_M_logK->col_name[i]))+273.15;

		if(ICS->SM_basis->mol_weight != NULL) calculate_mol_weight_mineral(ICS->SM_mineral, ICS->SM_basis);
		/* define types */
		if(ICS->SM_mineral->type == NULL) ICS->SM_mineral->type = (int *) calloc(ICS->SM_mineral->size[0],sizeof(int));
		for(i=0;i<ICS->SM_mineral->size[0];++i) ICS->SM_mineral->type[i] = -1;


	/* only using the reduced basis in the rest of the program */
		
	free(bs);
	free_db(Aq_db);
	free_db(Mineral_db);
	free_db(Basis_db);
	free_db(Aq_logK_db);
	free_db(M_logK_db);

	if(ICS->PRINT_DEBUG_CHEM) write_chem_struct(ICS->SM_all, "SM_all.out");
	if(ICS->PRINT_DEBUG_CHEM) write_chem_struct(ICS->SM_mineral, "SM_mineral.out");
	if(ICS->PRINT_DEBUG_CHEM) write_chem_struct(ICS->SM_aq_logK, "SM_aq_logK.out");
	if(ICS->PRINT_DEBUG_CHEM) write_chem_struct(ICS->SM_M_logK, "SM_M_logK.out");
	if(ICS->PRINT_DEBUG_CHEM) write_chem_struct(ICS->SM_basis, "SM_basis.out");

	ICS->mem->sizeSM = ICS->SM_all->size[0]+2;
	ICS->mem->sizeV  = ICS->SM_all->size[1]+2;
	InitMem(ICS->mem);

	/* test - remove Cl */
/*
	sprintf(name_tt,"Cl");
	InitICS(&ICS_new, ICS->size);
	remove_species(&ICS_new, name_tt, ICS);
	if(Vchem->ICS->PRINT_DEBUG_CHEM) write_chem_struct(ICS_new.SM_all, "SM_all2.out");
	if(Vchem->ICS->PRINT_DEBUG_CHEM) write_chem_struct(ICS_new.SM_mineral, "SM_mineral2.out");
	if(Vchem->ICS->PRINT_DEBUG_CHEM) write_chem_struct(ICS_new.SM_aq_logK, "SM_aq_logK2.out");
	if(Vchem->ICS->PRINT_DEBUG_CHEM) write_chem_struct(ICS_new.SM_M_logK, "SM_M_logK2.out");
	if(Vchem->ICS->PRINT_DEBUG_CHEM) write_chem_struct(ICS_new.SM_basis, "SM_basis2.out");
*/


	return;

}

/* 
copy all information from bulk species to surface species
but we put all species as determined by mass balance
*/
void chem_copy_bulk(struct InitChem *ICS_s, struct InitChem *ICS_b)
{
	int i;

	ICS_s->Temp = ICS_b->Temp;
	ICS_s->PRINT_DEBUG_CHEM = ICS_b->PRINT_DEBUG_CHEM;
	for(i=0;i<ICS_b->size;++i)
	{
		ICS_s->c_buffer[i] = 0.;/* no buffer! */
		ICS_s->c_vchem[i]     = ICS_b->c_vchem[i];
		ICS_s->log_af[i]   = ICS_b->log_af[i];
		ICS_s->pos[i]      = ICS_b->pos[i];
		ICS_s->c_ads[i]    = ICS_b->c_ads[i];
		strcpy(ICS_s->species_name[i],ICS_b->species_name[i]);
	}

	for(i=0;i<ICS_b->size_rock;++i)
	{
		ICS_s->pos_mass[i]=ICS_b->pos_rock[i];
	}

	for(i=0;i<ICS_b->size_mass;++i)
	{
		ICS_s->pos_mass[i+ICS_b->size_rock]=ICS_b->pos_mass[i];
	}
	ICS_s->size_mass = ICS_b->size_mass + ICS_b->size_rock;
	ICS_s->size_rock = 0;

	if(ICS_b->size_sup_min >0)
	{
		ICS_s->size_sup_min = ICS_b->size_sup_min;

		for(i=0;i<ICS_b->size;++i)
		{
			ICS_s->rate[i][0]=ICS_b->rate[i][0];
			ICS_s->rate[i][1]=ICS_b->rate[i][1];
			ICS_s->rate[i][2]=ICS_b->rate[i][2];
			ICS_s->rate[i][3]=ICS_b->rate[i][3];
		}

		for(i=0;i<ICS_b->size_sup_min;++i)
		{
			strcpy(ICS_s->sup_min_name[i], ICS_b->sup_min_name[i]);
			ICS_s->c_sup_min[i]  = ICS_b->c_sup_min[i];
			ICS_s->log_af_sup[i] = ICS_b->log_af_sup[i];
			ICS_s->pos_sup_min[i] = ICS_b->pos_sup_min[i];
		}
	}

	ICS_s->mem->sizeSM = ICS_b->mem->sizeSM;
	ICS_s->mem->sizeV  = ICS_b->mem->sizeV;
	InitMem(ICS_s->mem);

	ICS_s->CP = ICS_b->CP;

}

/* 
copy all information from bulk species to surface species
but we put all species as determined by mass balance
*/
void chem_remove_surface(struct InitChem *ICS_b, struct InitChem *ICS)
{
	int i, pos,j;

	ICS_b->Temp = ICS->Temp;
	ICS_b->PRINT_DEBUG_CHEM = ICS->PRINT_DEBUG_CHEM;
	ICS_b->size=0;
	ICS_b->size_mass = 0;
	ICS_b->size_rock = 0;
	for(i=0;i<ICS->size;++i)
	{
		pos = ICS->pos[i];
		
		if( ICS->SM_basis->type[pos] == 0 )
		{
			strcpy(ICS_b->species_name[ICS_b->size],ICS->species_name[i]);
			
			ICS_b->c_ads[ICS_b->size]    = ICS->c_ads[i];
			ICS_b->c_buffer[ICS_b->size] = ICS->c_buffer[i];
			ICS_b->log_af[ICS_b->size]   = ICS->log_af[i];
			ICS_b->pos[ICS_b->size]      = ICS->pos[i];
			ICS_b->c_vchem[ICS_b->size]     = ICS->c_vchem[i];
			
			for(j=0;j<ICS->size_mass;++j)
			{
				if(i == ICS->pos_mass[j])
				{
					ICS_b->pos_mass[ICS_b->size_mass] = ICS_b->size;
					ICS_b->size_mass++;
				}
			} 

			for(j=0;j<ICS->size_rock;++j) /*AH*/
			{
				if(i == ICS->pos_rock[j])
				{
					ICS_b->pos_rock[ICS_b->size_rock] = ICS_b->size;
					ICS_b->pos_buffer[ICS_b->size_rock] =ICS->pos_buffer[j];
					strcpy(ICS_b->buffer_name[ICS_b->size],ICS->buffer_name[i]);
					ICS_b->size_rock++;
				}
			} 
		ICS_b->size++;
		}
	}

	ICS_b->mem->sizeSM = ICS->mem->sizeSM;
	ICS_b->mem->sizeV  = ICS->mem->sizeV;
	InitMem(ICS_b->mem);

	ICS_b->CP = ICS->CP;

}


//void make_surface_bulk_struct(struct InitChem *ICS, struct BasVec *Vchem)
//{
//	int i, surface, pos;
//	char spec_tmp[200];
//	struct InitChem *ICS_s, *ICS_b;
//
//	surface = 0;
//	/* check if surface species are present */
//	for(i=0;i<ICS->size;++i)
//	{
//		pos = ICS->pos[i];
//		if(ICS->SM_basis->type[pos] != 0)
//		{
//			if(surface == 0)
//			{
//				sprintf(spec_tmp,"%s ",ICS->species_name[i]);
//				surface++;
//			} else
//			{
//				strcat(spec_tmp, ICS->species_name[i]);
//				strcat(spec_tmp," ");
//				surface++;
//			}
//		}
//	}
//
//	if(surface==0)
//	{
//		make_basvec_struct(ICS,Vchem);
//		if(Vchem->ICS->PRINT_DEBUG_CHEM) write_BasVec_struct(Vchem, "vchem.out");
//	}
//	else /* split in surface and bulk part */
//	{
//		ICS_s = (struct InitChem *) malloc(sizeof(struct InitChem) * 1);
//		ICS_b = (struct InitChem *) malloc(sizeof(struct InitChem) * 1);
//		InitICS(ICS_s, ICS->size);
//		InitICS(ICS_b, ICS->size-surface);
//
//		chem_copy_bulk(ICS_s,ICS);      /* remove species buffered by the rock */
//
//		ICS_s->SM_all      = ICS->SM_all;
//		ICS_s->SM_aq_logK  = ICS->SM_aq_logK;
//		ICS_s->SM_basis    = ICS->SM_basis;
//		ICS_s->SM_M_logK   = ICS->SM_M_logK;
//		ICS_s->SM_mineral  = ICS->SM_mineral;
//
//
//		chem_remove_surface(ICS_b,ICS);         /* remove surface species from ICS databases */
//		remove_species(ICS_b, spec_tmp, ICS_s); /* remove surface species from the database - stored in ICS_b*/
//		make_basvec_struct(ICS_b,Vchem);
//
//		Vchem->Vsurf = (struct BasVec *) calloc(1, sizeof( struct BasVec ));
//		make_basvec_struct(ICS_s,Vchem->Vsurf);
//
//
//		if(Vchem->ICS->PRINT_DEBUG_CHEM) write_BasVec_struct(Vchem, "vchem_b.out");
//		if(Vchem->ICS->PRINT_DEBUG_CHEM) write_BasVec_struct(Vchem->Vsurf, "vchem_s.out");
//
//	}
//
//}

void make_surface_bulk_struct_split(struct InitChem *ICS, struct InitChem *ICS_s, struct InitChem *ICS_b, struct BasVec *Vchem)
{

	if(ICS_s == NULL && ICS_b == NULL)
	{
		make_basvec_struct(ICS,Vchem);
		if(Vchem->ICS->PRINT_DEBUG_CHEM) write_BasVec_struct(Vchem, "vchem.out");
	}
	else /* split in surface and bulk part */
	{
		make_basvec_struct(ICS_b,Vchem);
		
		Vchem->Vsurf = (struct BasVec *) calloc(1, sizeof( struct BasVec ));
		make_basvec_struct(ICS_s,Vchem->Vsurf);


		if(Vchem->ICS->PRINT_DEBUG_CHEM) write_BasVec_struct(Vchem, "vchem_b.out");
		if(Vchem->ICS->PRINT_DEBUG_CHEM) write_BasVec_struct(Vchem->Vsurf, "vchem_s.out");

	}

}

//void splitICS(struct InitChem *ICS, struct InitChem **ICS_s, struct InitChem **ICS_b)
//{
//	int i, surface, pos;
//	char spec_tmp[200];
//
//	surface = 0;
//	/* check if surface species are present */
//	for(i=0;i<ICS->size;++i)
//	{
//		pos = ICS->pos[i];
//		if(ICS->SM_basis->type[pos] != 0)
//		{
//			if(surface == 0)
//			{
//				sprintf(spec_tmp,"%s ",ICS->species_name[i]);
//				surface++;
//			} else
//			{
//				strcat(spec_tmp, ICS->species_name[i]);
//				strcat(spec_tmp," ");
//				surface++;
//			}
//		}
//	}
//
//	if(surface==0)
//	{
//		/* no need to split */
//		*ICS_s = NULL;
//		*ICS_b = NULL;
//		return;
//	}
//	else /* split in surface and bulk part */
//	{
//		(*ICS_s) = (struct InitChem *) malloc(sizeof(struct InitChem ) * 1);
//		(*ICS_b) = (struct InitChem *) malloc(sizeof(struct InitChem ) * 1);
//		InitICS((*ICS_s), ICS->size);
//		InitICS(*(ICS_b), ICS->size-surface);
//
//		chem_copy_bulk((*ICS_s),ICS);      /* remove species buffered by the rock */
//
//		(*ICS_s)->SM_all      = ICS->SM_all;
//		(*ICS_s)->SM_aq_logK  = ICS->SM_aq_logK;
//		(*ICS_s)->SM_basis    = ICS->SM_basis;
//		(*ICS_s)->SM_M_logK   = ICS->SM_M_logK;
//		(*ICS_s)->SM_mineral  = ICS->SM_mineral;
//
//
//		chem_remove_surface((*ICS_b),ICS);         /* remove surface species from ICS databases */
//		remove_species((*ICS_b), spec_tmp, (*ICS_s)); /* remove surface species from the database - stored in ICS_b*/
//	}
//
//}


/* remove a species from the databases */
/* name is a string of species to be removed */
//void remove_species(struct InitChem *ICS_new, char *name, struct InitChem *ICS_old)
//{
//	int i;
//	real *bs;
//	char *s_tmp;
//	char s_bas_tmp[200];
//	char ct[] = "\n \t \0";
//
//	/* make copy of string */
//	if(sprintf(s_bas_tmp,name)<0)
//	{
//		printf("allocate more char space in routine : remove_species\n");
//		exit(1);
//	}
//
//	bs = (real *) malloc(sizeof(real)*ICS_old->SM_basis->size[0]);
//
//	/* bs = 1 for species in the calculation 0 otherwise */
//	for(i=0;i<ICS_old->SM_basis->size[0];++i) bs[i] = 1.;
//
//	s_tmp = strtok(s_bas_tmp,ct);
//	while(s_tmp != NULL)
//	{
//		for(i=0; i<ICS_old->SM_basis->size[0];++i)
//		{
//			if(!strcmp(ICS_old->SM_basis->row_name[i],s_tmp))
//			{
//				bs[i]=0.;
//			}
//		}
//
//		s_tmp = strtok(NULL,ct);
//	}
//	/* remove species */
//	remove_element(ICS_old->SM_all     ,bs,1,1,ICS_new->SM_all);
//	remove_element(ICS_old->SM_mineral ,bs,1,1,ICS_new->SM_mineral);
//	remove_element_rc(ICS_old->SM_basis,bs,NULL,1,1,ICS_new->SM_basis);
//
//	free(bs);
//	bs  = (real *) malloc(ICS_old->SM_all->size[0]*sizeof(real));
//	for(i=0;i<ICS_old->SM_all->size[0];++i) bs[i] = 0.;
//	/* 1 for all elements in SM_all_ */
//	for(i=0;i<ICS_new->SM_all->size[0];++i) bs[ICS_new->SM_all->abs_pos[i]] = 1.;
//	remove_element_rc(ICS_old->SM_aq_logK,bs,NULL,1,1,ICS_new->SM_aq_logK);
//
//	free(bs);
//	bs  = (real *) malloc(ICS_old->SM_mineral->size[0]*sizeof(real));
//	for(i=0;i<ICS_old->SM_mineral->size[0];++i) bs[i] = 0.;
//	/* 1 for all elements in SM_all_ */
//	for(i=0;i<ICS_new->SM_mineral->size[0];++i) bs[ICS_new->SM_mineral->abs_pos[i]] = 1.;
//	remove_element_rc(ICS_old->SM_M_logK,bs,NULL,1,1,ICS_new->SM_M_logK);
//
//	/* check if something is wrong in the database */
//	compare_names(ICS_new->SM_aq_logK,ICS_new->SM_all);
//	compare_names(ICS_new->SM_M_logK,ICS_new->SM_mineral);
//
//	free(bs);
//}

void chem_get_pos(struct InitChem *ICS)
{
	int i, j, cc, length; 
		/* check if everything went ok */
	length = 0;
	for(i=0;i<ICS->size;++i)
	{
		for(j=0;j<ICS->SM_basis->size[0];++j)
		{
			if(!strcmp(ICS->SM_basis->row_name[j],ICS->species_name[i]))
			{
				ICS->pos[i] = j;
				length++;
			}
		}
		if(i != length-1) /* species not found in the database */
		{
			printf("No species named: %s\n", ICS->species_name[i]);
			printf("VALID BASIS:\n");
			for(cc=0;cc<ICS->SM_basis->size[0];++cc) printf("%s\t",ICS->SM_basis->row_name[cc]);
			printf("\n");
			exit(0);
		}
	}

	/* check if all the mineral species are valid species */
	length = 0;
	for(i=0;i<ICS->size_rock;++i)
	{
		for(j=0;j<ICS->SM_mineral->size[0];++j)
		{
			if(!strcmp(ICS->SM_mineral->row_name[j],ICS->buffer_name[ICS->pos_rock[i]]))
			{
				ICS->pos_buffer[length] = j;
				length++;

			}
		}
		if(i != length-1) /* species not found in the database */
		{
			printf("No buffer named: %s\n", ICS->buffer_name[ICS->pos_rock[i]]);
			printf("VALID BUFFERS:\n");
			for(cc=0;cc<ICS->SM_mineral->size[0];++cc) printf("%s\t",ICS->SM_mineral->row_name[cc]);
			printf("\n");
			exit(0);
		}
	}
	
	/* check if all the mineral species are valid species */
	length = 0;
	for(i=0;i<ICS->size_sup_min;++i)
	{
		for(j=0;j<ICS->SM_mineral->size[0];++j)
		{
			if(!strcmp(ICS->SM_mineral->row_name[j],ICS->sup_min_name[i]))
			{
				ICS->pos_sup_min[length] = j;
				length++;

			}
		}
		if(i != length-1) /* species not found in the database */
		{
			printf("No mineral named: %s\n", ICS->sup_min_name[i]);
			printf("VALID BUFFERS:\n");
			for(cc=0;cc<ICS->SM_mineral->size[0];++cc) printf("%s\t",ICS->SM_mineral->row_name[cc]);
			printf("\n");
			exit(0);
		}
	}
}
/*
chem_init_boundary_node:
  make boundary struct and initialise the concentrations
  INPUT:  1) *rock_t        list of names of the species buffered by a mineral
          2) *buffer_t      list of minerals corresponding to a species given in rock_t
		  3) *mass_t        list of names of the basis species that are to be determined by mass balance
		  4) *mineral_sup_t list of supersaturated minerals that are to be used in rate dependent calculations,
		                    i.e. minerals that equilibrate slow
          5) *c_rock        list of the aqueous concentration of species buffered by the rock
		  6) *c_buffer      list of the concentration of the minerals at the rock node
		  7) *c_mass        list of the aqueous concentration of species determined by mass balance
  RETURN: BasVec struct 
*/
//struct BasVec chem_init_boundary_node( struct InitChem *ICS)
//{
//	struct BasVec Vp;
//
//	make_surface_bulk_struct(ICS,&Vp);
//	init_ctot(&Vp);
//
//	return Vp;
//}

/* use chemical solver to redefine the LB concentrations */
//real chem_set_eq_conc(real *c_vchem_eq, struct InitChem *ICS, int call)
//{
//	struct BasVec V;
//	int c,pos;
//	real c_h;
//
//
//	make_surface_bulk_struct(ICS, &V);
//	V.equilibrate = 1;
//	if(V.Vsurf != NULL) V.Vsurf->equilibrate = 1;
//
//	for(c=0;c<ICS->size;++c) c_vchem_eq[c] = ICS->c_vchem[c];
//
//	solve_chem_bulk_surface(c_vchem_eq,call,&V);
//
//	for(c = 0; c < V.ICS->size; ++c)
//	{
//		pos = ICS->pos[c];
//		if(V.size_rock != 0 && V.Vsurf != NULL) c_vchem_eq[c] = V.ctot_calc[pos];
//		ICS->c_ads[c] = V.ctot_ads[pos];
//	}
//	c_h = V.ctot_calc[V.pos_pH];
//	/*NB!!!!!!!!!!!!!!!*/
//	/*REMEMBER TO FREE STRUCT */
//	if(V.Vsurf != NULL) free_BasVec(V.Vsurf);
//	free_BasVec(&V);
//
//return c_h;
//}

/* use chemical solver to redefine the LB concentrations */
/* switch rate minerals with equilibrium minerals */
real chem_set_eq_conc_nlin(real *c_vchem_eq, struct InitChem *ICS, int call)
{
	struct BasVec V;
	int c,pos;
	real c_h;


	make_basvec_struct(ICS, &V);
	trans_nlin_eq(&V); /* add rate minerals to equilibrium minerals */

	V.equilibrate = 1;
	V.chbal = 1;
	for(c=0;c<ICS->size;++c) c_vchem_eq[c] = ICS->c_vchem[c];

	calc_surface_reaction_nlin(c_vchem_eq,1,&V);
	if(ICS->PRINT_DEBUG_CHEM) writeVchem(&V,"qqqV.out");
	for(c = 0; c < V.ICS->size; ++c)
	{
		pos = ICS->pos[c];		 
		ICS->c_ads[c] = V.ctot_ads_calc[pos];
	}
	c_h = V.ctot_calc[V.pos_pH];
	/*NB!!!!!!!!!!!!!!!*/
	/*REMEMBER TO FREE STRUCT */
	if(V.Vsurf != NULL) free_BasVec(V.Vsurf);
	free_BasVec(&V);

return c_h;
}

/* use chemical solver to redefine the LB concentrations */
real *chem_set_eq_conc_old(struct InitChem ICS)
{
	struct BasVec V;
	int c,pos;
	real *c_vchem_eq;


	make_basvec_struct(&ICS, &V);
	set_temperature_db(&V);
	if(V.pos_X > -1) V.equilibrate = 1;
	solve_chemistry(&V);
	/* new aqueous concentrations */
	calc_ctot_aq(&V);

	c_vchem_eq = (real *) calloc(ICS.size,sizeof(real));
	for(c = 0; c < ICS.size; ++c)
	{
		pos = ICS.pos[c];
		if(ICS.SM_basis->type[pos] == 0) c_vchem_eq[c] = V.ctot_calc[pos];
	}
	/*NB!!!!!!!!!!!!!!!*/
	/*REMEMBER TO FREE STRUCT */
return c_vchem_eq;
}
/* sets constants */
void  init_const_chem_param(struct ChemParam *CP)
{

      real a[] = {295.68, -1.2283,2.094e-3,-1.41e-6};
      real b[] = {5321,233.76,-0.9397,1.417e-3,-8.292e-7};
      real Temp = 273.15+130;

      if(Temp < 373) CP->ew = a[0] + a[1]*Temp + a[2]*Temp*Temp + a[3]*Temp*Temp*Temp; /* from Revil et al. 1999*/ 
      else CP->ew = b[0]/Temp + b[1] + b[2]*Temp + b[3]*Temp*Temp + b[4]*Temp*Temp*Temp;

      CP->e0 = 8.854187817e-12; /* permitivity of vacuum-Unit: A s V^-1 m^-1*/
   /* CP->ew = 80;  Dielectric  constant of water - This is temperature and pressure dependent*/
      CP->F  = 9.648456e+4; /*Faradays constant-Unit: C mol^-1 */
      CP->kB = 1.3806503e-23; /*Boltzmanns constant-Unit: m^2 kg s^-2 K^-1*/
      CP->Na = 6.0221415e-23; /* Avogadros number-Unit: mol-1 */
      CP->Rg = 8.31441;  /* Ideal gas constant-Unit: J K^-1 mol^-1 */
	  CP->SA = 5200;
}




/* check if row names are equal, and prints out those who are not */
void compare_names(struct ChemTable *Ma, struct ChemTable *Mb)
{
	int i;
	if(Ma->size[0] != Mb->size[0])
	{
		printf("Wrong dimensions in the database\n");
	}

	for(i=0;i<Ma->size[0];++i)
	{
		if(strcmp(Ma->row_name[i],Mb->row_name[i]))
		{
			printf("Something is wrong in the database\n");
			printf("Row name: %s does not match row name: %s\n", Ma->row_name[i], Mb->row_name[i]);
		}

	}

}

void InitICS(struct InitChem *ICS, int size)
{
	int i;

	ICS->num_phases = size;
	ICS->size = size;
	ICS->nlinr = 0;
	ICS->species_name =(char **) calloc(size,sizeof(char *));
	ICS->buffer_name =(char **) calloc(size,sizeof(char *));
	ICS->sup_min_name =(char **) calloc(size,sizeof(char *));

	for(i=0;i<size;++i)
	{
		ICS->species_name[i] = (char *) calloc(COMPLEXNAME,sizeof(char ));
		ICS->buffer_name[i]  = (char *) calloc(COMPLEXNAME,sizeof(char ));
		ICS->sup_min_name[i]  = (char *) calloc(COMPLEXNAME,sizeof(char ));


	}
	ICS->c_vchem      = (real *) calloc(size, sizeof(real));
	ICS->c_buffer     = (real *) calloc(size, sizeof(real));
	ICS->c_sup_min    = (real *) calloc(size, sizeof(real));
	ICS->c_ads        = (real *) calloc(size, sizeof(real));
	ICS->log_af       = (real *) calloc(size, sizeof(real));
	ICS->log_af_sup   = (real *) calloc(size, sizeof(real));
	ICS->Sg           = (real *) calloc(size, sizeof(real));
	ICS->pos_rock     = (int *) calloc(size, sizeof(int));
	ICS->pos_mass     = (int *) calloc(size, sizeof(int));
	ICS->pos_buffer   = (int *) calloc(size, sizeof(int));
	ICS->pos          = (int *) calloc(size, sizeof(int));
	ICS->pos_sup_min  = (int *) calloc(size, sizeof(int));
	ICS->pos_sup_bas  = (int *) calloc(size, sizeof(int));
	ICS->size_rock    = 0;
	ICS->size_mass    = 0;	
	ICS->size_sup_min = 0;

	ICS->rate     = (real **) calloc(ICS->size,sizeof(real));
	for(i=0;i<ICS->size;++i) ICS->rate[i] = (real *) calloc(4,sizeof(real));

	ICS->SM_all      = (struct ChemTable *) malloc(1 * sizeof(struct ChemTable));
	ICS->SM_basis    = (struct ChemTable *) malloc(1 * sizeof(struct ChemTable));
	ICS->SM_aq_logK  = (struct ChemTable *) malloc(1 * sizeof(struct ChemTable));
	ICS->SM_M_logK   = (struct ChemTable *) malloc(1 * sizeof(struct ChemTable));
	ICS->SM_mineral  = (struct ChemTable *) malloc(1 * sizeof(struct ChemTable));
	ICS->mem         = (struct MEMC      *) malloc(1 * sizeof(struct MEMC     ));
	ICS->CP          = (struct ChemParam *) malloc(1 * sizeof(struct ChemParam)); 
	ICS->BiP          = (struct BioParam *) malloc(1 * sizeof(struct BioParam)); /* AMIN */

}

void freeICS(struct InitChem *ICS)
{
	int i;

	for(i=0;i<ICS->size;++i)
	{
		free(ICS->species_name[i]);
		free(ICS->buffer_name[i]);
		free(ICS->rate[i]);
	}
	free(ICS->species_name);
	free(ICS->buffer_name);
	free(ICS->rate);
	free(ICS->c_vchem);    
	free(ICS->c_buffer);
	free(ICS->c_ads);  
	free(ICS->log_af);  
	free(ICS->Sg);
	free(ICS->pos_rock);
	free(ICS->pos_mass); 
	free(ICS->pos_buffer); 
	free(ICS->pos);   


	free_db(ICS->SM_all);
	free_db(ICS->SM_basis);    
	free_db(ICS->SM_aq_logK);  
	free_db(ICS->SM_M_logK);   
	free_db(ICS->SM_mineral); 
	free_mem(ICS->mem);

	free(ICS->SM_all);
	free(ICS->SM_basis);    
	free(ICS->SM_aq_logK);  
	free(ICS->SM_M_logK);   
	free(ICS->SM_mineral); 

	free_mem(ICS->mem);
	free(ICS->mem);
}


int get_int_from_string(char *s)
{
	while(*s != '\0' && (*s > '9' || *s < '0') ) s++;
	return atoi(s);
}

void free_mem(struct MEMC *mem)
{
	int i;  
	free(mem->dVma);
	free(mem->dVmb);
	free(mem->dVmc);
	free(mem->iVma);
	free(mem->dSMa);
	for(i=0;i<mem->sizeV;++i)
	{
		free(mem->ddVma[i]);
		free(mem->ddVmb[i]);
		free(mem->ddVmc[i]);
	}

	free(mem->ddVma);
	free(mem->ddVmb);
	free(mem->ddVmc);

}

void InitMem(struct MEMC *mem)
{
	int i; 

	mem->dVma  = (real *) malloc(mem->sizeV*sizeof(real));
	mem->dVmb  = (real *) malloc(mem->sizeV*sizeof(real));
	mem->dVmc  = (real *) malloc(mem->sizeV*sizeof(real));
	mem->iVma  = (int  *) malloc(mem->sizeV*sizeof(int));
	mem->dSMa  = (real *) malloc(mem->sizeSM*sizeof(real ));

	
	mem->ddVma  = (real **) malloc(mem->sizeV*sizeof(real *));
	mem->ddVmb  = (real **) malloc(mem->sizeV*sizeof(real *));
	mem->ddVmc  = (real **) malloc(mem->sizeV*sizeof(real *));
	

	for(i=0;i<mem->sizeV;++i)
	{
		mem->ddVma[i]  = (real *) malloc(mem->sizeV*sizeof(real ));
		mem->ddVmb[i]  = (real *) malloc(mem->sizeV*sizeof(real ));
		mem->ddVmc[i]  = (real *) malloc(mem->sizeV*sizeof(real ));
	}

	
}


