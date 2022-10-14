#include "chem_global.h"

void solve_chemistry(struct BasVec *Vchem)
{
  int iter,i,j;
  real psi_old;
  real cPH, Io_old, *c_dl;
  real x_dum, x_next, f_next, f_old, x_old, deltax;
  real delta_dl;

  struct ChemTable *SM_mineral, *SM_all;
  SM_mineral = Vchem->ICS->SM_mineral;
  SM_all     = Vchem->ICS->SM_all;


  if(Vchem->DL)
    {
      c_dl = (real *) calloc(Vchem->size,sizeof(real)); 
      for(i=0;i<Vchem->size;++i) c_dl[i]=0.;
    }
  /* initialize aqueous complexes */
  for(i=0;i<SM_all->size[0];++i)
    {
      SM_all->log_a[i]=SM_all->log_m[i]=-20.;
      SM_all->log_g[i]=0.;
    }

  /*--------------------- Ionic strength ------------------------*/
  if(Vchem->Io<1e-5) Vchem->Io=calc_Io(Vchem);

  /*--------------------- activity  ------------------------*/
  calc_act_DB(Vchem->ICS->CP->Adh, Vchem->ICS->CP->Bdh, Vchem->Io,Vchem);

  Io_old = Vchem->Io;

  if(Vchem->chbal) /* TRICK IS TO RUN 1 TIME WITHOUT SURFACE CHARGE */
    {
      x_old = Vchem->log_a[Vchem->pos_pH];
      x_next = x_old*1.01;
      //if (_DEBUG_FLAG_) printf("\nf_pH 1\n");
      f_old  = f_pH(x_old,1,Vchem, mass_balance);
      //if (_DEBUG_FLAG_) printf("\nf_pH 2\n");
      f_next = f_pH(x_next,1,Vchem,mass_balance);
      /* secant step */
      deltax = (x_next-x_old)/(f_next-f_old)*f_next;
      if(fabs(deltax) > PH_MAX_STEP_)
	{
	  if(deltax<0) deltax = -PH_MAX_STEP_;
	  else deltax = PH_MAX_STEP_;
	}
      x_dum  = x_next - deltax;
      x_old  = x_next;
      f_old  = f_next;
      x_next = x_dum;
      //if (_DEBUG_FLAG_) printf("\nf_pH 3\n");
      f_next = f_pH(x_next,1,Vchem,mass_balance);
      cPH = x_old;
    }
  else 
    {/* pH fixed */
      //if (_DEBUG_FLAG_) printf("\nf_pH 4\n");
      f_pH(Vchem->log_a[Vchem->pos_pH],0,Vchem,mass_balance);
      cPH = Vchem->log_a[Vchem->pos_pH];

    }

	

	
  iter=0;
  psi_old=Vchem->psi;
  for(j=0; j<= CHEM_MAX_PH_ITER_; ++j)
    {
      psi_old=Vchem->psi;
      iter++;

		
      /*-----------------------------STEP 2---------------------------------*/
      /*--------------------- calculate Io----------------------------------*/
      Io_old = Vchem->Io;
      Vchem->Io=calc_Io(Vchem);

      /*-----------------------------STEP 3---------------------------------*/
      /*--------------------- calculate gamma-------------------------------*/

      calc_act_DB(Vchem->ICS->CP->Adh, Vchem->ICS->CP->Bdh, Vchem->Io,Vchem);

      /*-----------------------------STEP 1---------------------------------*/
      /*--------------------- charge balance -------------------------------*/

      if(Vchem->chbal)
	{
	  /* secant step */
	  deltax = (x_next-x_old)/(f_next-f_old)*f_next;
	  if(fabs(deltax) > PH_MAX_STEP_)
	    {
	      if(deltax<0) deltax = -PH_MAX_STEP_;
	      else deltax = PH_MAX_STEP_;
	    }
	  x_dum  = x_next - deltax;
	  x_old  = x_next;
	  f_old  = f_next;
	  x_next = x_dum;
	  if(Vchem->surface_flag == 1)
	    {
	      //if (_DEBUG_FLAG_) printf("\nf_pH 5\n");
	      f_next = f_pH(x_next,iter,Vchem,mass_balance_surf);
	      if(Vchem->DL)
		{
		  calc_diffuse_layer_conc(c_dl,Vchem);
		  /* update */
		  delta_dl =0.;
		  for(i=0;i<Vchem->size;++i)
		    {
		      delta_dl += fabs(Vchem->ctot_dl[i]-c_dl[i]);
		      Vchem->ctot_dl[i] = c_dl[i];
		    }
		}
	    }
	  else { 
	    //if (_DEBUG_FLAG_) printf("\nf_pH 6\n");
	    f_next = f_pH(x_next,iter,Vchem,mass_balance);
	  }
	  cPH = x_old;
	}
      else 
	{
	  if(Vchem->surface_flag == 1)
	    {
	      //if (_DEBUG_FLAG_) printf("\nf_pH 7\n");
	      f_pH(Vchem->log_a[Vchem->pos_pH],1,Vchem,mass_balance_surf);
	      if(Vchem->DL)
		{
		  calc_diffuse_layer_conc(c_dl,Vchem);
		  /* update */
		  delta_dl =0.;
		  for(i=0;i<Vchem->size;++i)
		    {
		      delta_dl += fabs(Vchem->ctot_dl[i]-c_dl[i]);
		      Vchem->ctot_dl[i] = c_dl[i];
		    }
		}
	    }
	  else {
	    //if (_DEBUG_FLAG_) printf("\nf_pH 8\n");
	    f_pH(Vchem->log_a[Vchem->pos_pH],1,Vchem,mass_balance);
	  }

	  cPH=Vchem->log_a[Vchem->pos_pH];
	  f_next = 0.; /* charge balance is not a constraint */
	}
		

      if(fabs(cPH-Vchem->log_a[Vchem->pos_pH]) < 1e-6  && fabs(f_next) < 1e-8
	 && fabs(Vchem->psi-psi_old) < 1e-8 && fabs(Vchem->Io-Io_old) < 1e-5 )
	{
	  break;
	}

      /*
	printf("%lf \t %lf \n", Vchem->Io, Io_old);
      */
	
	

      /*-----------------------------STEP 5---------------------------------*/
      /*--------------------- surface charge ---------------------------------*/

      if(Vchem->ICS->PRINT_DEBUG_CHEM)printf("Sch_: %12.8e psi_: %12.8e pH %12.9e\n", Vchem->sch, Vchem->psi,-Vchem->log_a[Vchem->pos_pH]);
      if(Vchem->ICS->PRINT_DEBUG_CHEM)printf("Io_: %lf\n", Vchem->Io);
    }

  /* saturation index of minerals */
  calc_complex(SM_mineral,Vchem);

  if(Vchem->ICS->PRINT_DEBUG_CHEM)
    printf("No global iter %d and final pH %12.8e\n", iter, -Vchem->log_a[Vchem->pos_pH]);

  if(iter==CHEM_MAX_PH_ITER_) {
    //printf("WARNING charge balance not converged !!!!!!\n");/* No pH found*/
    Vchem->ICS->num_chrg_not_conv++;
  }
  
  if(Vchem->ICS->PRINT_DEBUG_CHEM) calc_ctot_aq(Vchem);

  if(Vchem->DL) free(c_dl);


}
/* calculates ionic strength */
real calc_Io(struct BasVec *Vchem)
{
	real Io;
	int i;
	struct ChemTable *SM_all, *SM_basis;
	SM_all     = Vchem->ICS->SM_all;
	SM_basis   = Vchem->ICS->SM_basis;

	Io = 0.;
	/* basis species */
	for(i=0;i<Vchem->size; ++i) if(SM_basis->type[i] == 0) Io +=.5*SM_basis->charge[i]*SM_basis->charge[i]*pow10(Vchem->log_m[i]);
	/* secondary species */
	for(i=0;i<SM_all->size[0];++i) if(SM_all->type[i] == 0)Io +=.5*SM_all->charge[i]*SM_all->charge[i]*pow10(SM_all->log_m[i]);
return Io;
}

/* calculate activity of Water */
/* use same algorithm as PHREEQC */

/* calculates the activity based on the Debye Huckel formulation */
void calc_act_DB(real A, real B, real Io,struct BasVec *Vchem)
{
	real Io_sqrt;
	int i;
	struct ChemTable *SM_all, *SM_basis;
	SM_all     = Vchem->ICS->SM_all;
	SM_basis   = Vchem->ICS->SM_basis;

	Io_sqrt = sqrt(Io);

	/* basis species */
	for(i=0;i<Vchem->size; ++i)
	{
		Vchem->log_g[i] =-A*SM_basis->charge[i]*SM_basis->charge[i]*Io_sqrt/(1.+B*SM_basis->a0[i]*Io_sqrt);
	}
	/* secondary species */
	for(i=0;i<SM_all->size[0];++i)
	{
		SM_all->log_g[i]=-A*SM_all->charge[i]*SM_all->charge[i]*Io_sqrt/(1+B*SM_all->a0[i]*Io_sqrt);
	}
}


/*
solve_chem_boundary_node:
  solve equilibrium chemistry at the boundary node, calculates amount of mineral precipitated and
  dissolved
  INPUT: 1) c_tot : total concentration of aqueous species 
         2) basis_list: list of the position to basis species that are to be used in basis transformation
		 3) size_basis_list
		 4) buffer_list: position of the buffer minerals that are to be used in basis transformation
		 5) Vchem
  OUTPUT: 1) ctot : updated total concentrations of aqueous species
          mineral concentration at the surface is stored in Vchem->ctot_mineral
  RETURNS 1 if chem solver was called 0 othervise
*/

int solve_chem_boundary_node(real *ctot, int call, struct BasVec *Vchem)
{
	int i, pos, length;
	int chem_solver_call;
	real *c_updated;
	//real CAdh[]={0.51 , -1.154, 35.697, -182.023, 346.528},
	//real CBdh[]={0.325, 0.08  , 1.441 , -6.541  , 11.655};
	struct InitChem *ICS;
	struct ChemTable *SM_mineral, *SM_all, *SM_aq_logK, *SM_M_logK;
	//struct splayTreeList *TL_dum;

	SM_mineral = Vchem->ICS->SM_mineral;
	SM_all     = Vchem->ICS->SM_all;
	SM_aq_logK = Vchem->ICS->SM_aq_logK;
	SM_M_logK  = Vchem->ICS->SM_M_logK;

	ICS = Vchem->ICS;

	
	/* get correct logK values */
	set_temperature_db(Vchem);

	/* update total concentrations at the surface */
	
	chem_solver_call = 0;
	/*
	if(call || call_chem_solver(ctot, ICS, Vchem) )
	{
	*/
	length = -1;
	chem_solver_call = 1;
	for(i=0;i<ICS->size;++i)
	{
		pos = ICS->pos[i];
		Vchem->ctot[pos] = ctot[i];
	}
	do
	{
		/* solve equilibrium chemistry */
		solve_chemistry(Vchem);
		if(Vchem->ICS->PRINT_DEBUG_CHEM) writeVchem(Vchem,"aaajulius.out");

		/* 
		check if the mineral is one of the basis minerals, and returns the position
		to the basis mineral
		*/
			
		if(!Vchem->equilibrate) /* rock buffered */
		{
			length = check_if_mineral_supersaturated(ICS->pos_buffer, ICS->size_rock, 
				Vchem->pos_buffer, Vchem->size_rock, Vchem);
			if(length>=0) /* length now holds the position to the basis mineral in the buffer_list */
			{
				printf("Add buffer mineral %s\n", SM_mineral->row_name[ICS->pos_buffer[length]]);
				add_new_buffer_mineral(ICS->pos[ICS->pos_rock[length]], ICS->pos_buffer[length], Vchem);
			}
		}
		/* non-linear rate equations */				
		if(ICS->size_sup_min > 0 && (Vchem->size_sup_min != ICS->size_sup_min) && (!Vchem->equilibrate))
		{
			length = check_if_mineral_supersaturated(ICS->pos_sup_min, ICS->size_sup_min, 
				Vchem->pos_sup_min, Vchem->size_sup_min, Vchem);
			if(Vchem->ICS->PRINT_DEBUG_CHEM) writeVchem(Vchem,"debug.out");
			if(length>=0)
			{
				printf("Add buffer mineral %s\n", SM_mineral->row_name[ICS->pos_sup_min[length]]);
				add_new_sup_mineral(length, ICS, Vchem);
			}
		}

			
	} while(length >= 0);
		/* new aqueous concentrations */
	calc_ctot_aq(Vchem);

	if(!Vchem->nlinr && Vchem->size_rock >0)
	{
		c_updated = (real *) calloc(Vchem->size,sizeof(real));
		for(i=0;i<ICS->size;++i) c_updated[ICS->pos[i]] = ctot[i];
		/* reaction_rate_lb(c_updated, Vchem);*/
		reaction_rate_julius(c_updated, Vchem);
		for(i=0;i<ICS->size;++i) ctot[i]= c_updated[ICS->pos[i]];
		free(c_updated);
	}
	if(Vchem->nlinr)
	{
		c_updated = (real *) malloc(Vchem->size*sizeof(real));
		for(i=0;i<Vchem->size;++i) c_updated[i] = 0.;
		calc_ctot_sup_min(c_updated, Vchem);
		calc_mineral_change_julius(Vchem);
		/*
		calc_ads(Vchem);
		*/
		for(i=0;i<ICS->size;++i)
		{
			/*ctot[i] = ctot[i]-Vchem->ctot_ads_calc[ICS->pos[i]]-c_updated[ICS->pos[i]];*/
			ctot[i] = Vchem->ctot_calc[ICS->pos[i]];
			/*
			Vchem->ctot_calc[ICS->pos[i]] -= c_updated[ICS->pos[i]];
			ctot[i] = Vchem->ctot_calc[ICS->pos[i]];
			*/	
		}
						
		free(c_updated);
	}
	if(Vchem->size_sup_min == 0 && Vchem->size_rock == 0)
	{
		for(i=0;i<ICS->size;++i) ctot[i] = Vchem->ctot_calc[ICS->pos[i]];
	}

	do
	{
		length = check_if_mineral_conc_neg(Vchem->pos_buffer,Vchem->size_rock, Vchem);

		if(length >= 0)
		{
			printf("Buffer mineral %s exhausted\n", SM_mineral->row_name[Vchem->pos_buffer[length]]);
			remove_buffer_mineral(length, Vchem);
		}

		length = check_if_mineral_conc_neg(Vchem->pos_sup_min,Vchem->size_sup_min, Vchem);

		if(length >= 0)
		{
			printf("Sup Buffer mineral %s exhausted (solve_chem_boundary_node)\n", SM_mineral->row_name[Vchem->pos_sup_min[length]]);
			remove_sup_mineral(length, Vchem);
		}

	}while(length>0);

	return chem_solver_call;
}

int solve_chem_boundary_node_interpolate(real *ctot, real *ctot_min, int key, struct BasVec *Vchem)
{
  int i, pos, length, new_key;
  //real dum, fact;
  //real CAdh[]={0.51 , -1.154, 35.697, -182.023, 346.528},
  //  CBdh[]={0.325, 0.08  , 1.441 , -6.541  , 11.655};
  struct InitChem *ICS;
  struct ChemTable *SM_mineral, *SM_all, *SM_aq_logK, *SM_M_logK;

  SM_mineral = Vchem->ICS->SM_mineral;
  SM_all     = Vchem->ICS->SM_all;
  SM_aq_logK = Vchem->ICS->SM_aq_logK;
  SM_M_logK  = Vchem->ICS->SM_M_logK;

  ICS = Vchem->ICS;

  new_key = key;
  length = -1;
  
		

  for(i=0;i<Vchem->ICS->size_sup_min;++i) Vchem->ctot_mineral[Vchem->ICS->pos_sup_min[i]] = ctot_min[i];

  length = check_if_mineral_conc_neg(Vchem->pos_sup_min,Vchem->size_sup_min, Vchem);
  
  if(length >= 0)
  {
	  printf("Sup Buffer mineral %s exhausted (solve_chem_boundary_node_interpolate)\n", SM_mineral->row_name[Vchem->pos_sup_min[length]]);
	  new_key = gen_new_key(ctot_min, Vchem->ICS->size_sup_min);
	  return new_key;		
  }



  /* get correct logK values */
 if(JULIUS) set_temperature_db(Vchem);
  /* update total concentrations at the surface */
	
  for(i=0;i<ICS->size;++i)
  {
	  pos = ICS->pos[i];
	  Vchem->ctot[pos] = ctot[i];
  }
       
  /* solve equilibrium chemistry */
  //if (_DEBUG_FLAG_) printf("\nsolve_chemistry\n");
  solve_chemistry(Vchem);
  if(Vchem->ICS->PRINT_DEBUG_CHEM) writeVchem(Vchem,"aaajulius.out");
  /* 
  check if the mineral is one of the basis minerals, and returns the position
  to the basis mineral
  */
  	
   /* non-linear rate equations */				
  if(ICS->size_sup_min > 0 && (Vchem->size_sup_min != ICS->size_sup_min) && (!Vchem->equilibrate))
  {
	  length = check_if_mineral_supersaturated(ICS->pos_sup_min, ICS->size_sup_min, 
		  Vchem->pos_sup_min, Vchem->size_sup_min, Vchem);
  }
  if(length >=0)
  {
	  if(Vchem->ICS->PRINT_DEBUG_CHEM) printf("Add buffer mineral %s\n", SM_mineral->row_name[ICS->pos_sup_min[length]]);
	  if(Vchem->ICS->PRINT_DEBUG_CHEM) writeVchem(Vchem,"debugi.out");
	  ctot_min[length]=1.e-16;
	  Vchem->ctot_mineral[Vchem->ICS->pos_sup_min[length]] = 1.e-16;
	  new_key = gen_new_key(ctot_min, Vchem->ICS->size_sup_min);
	  return new_key;
  }

  /* new aqueous concentrations */
  calc_ctot_aq(Vchem);

  /* @ah we store the change in concentrations IN-OUT */
  for(i=0;i<ICS->size;++i)  ctot[i]     -= Vchem->ctot_calc[ICS->pos[i]];		
 

  	
  return new_key;

}

void calc_ads(struct BasVec *Vchem)
{
	int i,j, pos;

	for(i=0;i<Vchem->ICS->size;++i)
		 {
			 pos = Vchem->ICS->pos[i];

			 Vchem->ctot_ads_calc[pos]=0.;
			 for(j=0;j<Vchem->ICS->SM_all->size[0];++j)
			 {
				 if(Vchem->ICS->SM_all->type[j] != 0)
					 Vchem->ctot_ads_calc[pos] += Vchem->ICS->SM_all->M[j][pos]*pow10(Vchem->ICS->SM_all->log_m[j]);
			 }
			 
		 }
}




int my_floor(real x)
{
  int ret = (int) x;
  
  if ( x < 0)
    return ret-1;
  else
    return ret;
}

/*
 *  0 : same key value
 *  1 : key1 is larger than key2
 *  2 : key2 is larger than key1
 */

/* INPUT: ctot_aq: is the aqueous concentration of species i, at exit it holds the
                   aqueous concentration after adsorption and precip/diss reactions 
*/
void solve_chem_bulk_surface(real *ctot_aq, int call, struct BasVec *Vchem)
{
	int i;
	real *ctot_aq_in;

	ctot_aq_in = (real *) malloc(sizeof(real)* Vchem->ICS->size);

	for(i=0;i<Vchem->ICS->size;++i) ctot_aq_in[i] = ctot_aq[i];


	if(Vchem->equilibrate == 1) /* first equilibrate with minerals and then with surface */
	{
	
		if(Vchem->Vsurf != NULL)
		{
			calc_surface_reaction(ctot_aq, call, Vchem);
			
		} else
		{
			solve_chem_boundary_node(ctot_aq,call,Vchem);
			for(i=0; i<Vchem->ICS->size;++i) ctot_aq[i] = Vchem->ctot_calc[Vchem->ICS->pos[i]];
		}
	}
	else /* rate equation - surface reactions are faster */ 
	{
		if(Vchem->Vsurf != NULL)
		{
			calc_surface_reaction(ctot_aq, call, Vchem);
		} else if(Vchem->size_rock > 0) solve_chem_boundary_node(ctot_aq,call,Vchem);
		else if (Vchem->size_sup_min >0 && Vchem->Vsurf == NULL) solve_chem_boundary_node(ctot_aq,call,Vchem);
		if(Vchem->ICS->PRINT_DEBUG_CHEM) writeVchem(Vchem, "Vchem_init.dat");
	}

		
/*
if(fabs(Vchem->ctot_aq_delta[1]) > 1e-8 ) printf("%s too large: %g\n",Vchem->ICS->SM_basis->row_name[1],Vchem->ctot_aq_delta[1]);  
*/
	free(ctot_aq_in);

}

void solve_chem_bulk_surface_new(real *ctot_aq, int call, struct BasVec *Vchem)
{
	int i, j,pos;
	real *ctot_aq_in, *ctot;

	ctot_aq_in = (real *) malloc(sizeof(real)* Vchem->ICS->size);
	ctot       = (real *) malloc(sizeof(real)* Vchem->ICS->size);

	for(i=0;i<Vchem->ICS->size;++i) ctot_aq_in[i] = ctot_aq[i];
	for(i=0;i<Vchem->ICS->size;++i)
	{
		pos = Vchem->ICS->pos[i];
		ctot[i] = ctot_aq[i] + Vchem->ctot_ads[pos];     /*release surface species in bulk, and recalculate */
		if(pos == Vchem->pos_X) ctot[i] = Vchem->ctot[pos];
	}

	solve_chem_boundary_node(ctot,call,Vchem);

	for(i=0;i<Vchem->size;++i)/* calc adsorption */
	{
		Vchem->ctot_ads_calc[i]=0.;
		for(j=0;j<Vchem->ICS->SM_all->size[0];++j)
		{
			if(Vchem->ICS->SM_all->type[j] != 0)
			{/*surface species*/
				Vchem->ctot_ads_calc[i] += Vchem->ICS->SM_all->M[j][i]*pow10(Vchem->ICS->SM_all->log_m[j]);
			}
		}
	}

	for(i=0;i<Vchem->ICS->size;++i) 
	{
		pos = Vchem->ICS->pos[i];
		if(pos != Vchem->pos_X) ctot_aq[i] = ctot[i]-Vchem->ctot_ads_calc[pos];
	}
		
	free(ctot_aq_in);
	free(ctot);

}
void calc_surface_reaction_nlin(real *ctot_aq, int call, struct BasVec *Vchem)
{	
	int i, pos, chem_calc;
	real *ctot;

	/* calculate adsorption first - assume adsorption is a fast process*/
	/* add adsorbed rock species to total concentration*/

	/* NOTE ctot_aq is relative to Vchem not Vp !*/

	ctot = (real *) malloc(Vchem->ICS->size*sizeof(real));
	for(i=0;i<Vchem->ICS->size;++i)
	{
		pos = Vchem->ICS->pos[i];
		ctot[i] = ctot_aq[i] + Vchem->ctot_ads[pos];     /*release surface species in bulk, and recalculate */
		if(pos == Vchem->pos_X) ctot[i] = Vchem->ctot[pos];
	}

	 chem_calc = solve_chem_boundary_node(ctot,call, Vchem);

	 if(chem_calc) calc_ads(Vchem);
	 
	 /* else use previous calculated ctot_ads, stored in ctot_ads */

	 for(i=0;i<Vchem->ICS->size;++i)
	 {
		 pos = Vchem->ICS->pos[i];
		 if(pos != Vchem->pos_X) ctot_aq[i]=ctot[i]; /* leave exchange sites unaltered */
	 }


	
	if(Vchem->ICS->PRINT_DEBUG_CHEM) writeVchem(Vchem, "Vchem_init.dat");
	free(ctot);

	
	
}
void calc_surface_reaction(real *ctot_aq, int call, struct BasVec *Vchem)
{	
	int i, j, pos, chem_calc, posVp;
	struct BasVec *Vp;
	real *ctot, *ctotVp;

	Vp = Vchem->Vsurf;
	/* calculate adsorption first - assume adsorption is a fast process*/
	/* add adsorbed rock species to total concentration*/

	/* NOTE ctot_aq is relative to Vchem not Vp !*/

	ctot = (real *) malloc(Vp->ICS->size*sizeof(real));
	ctotVp = (real *) malloc(Vp->size*sizeof(real));
	for(i=0;i<Vp->size;++i) ctotVp[i] = 0.;
	for(i=0;i<Vchem->ICS->size;++i) ctotVp[Vchem->ICS->pos[i]] = ctot_aq[i];
	for(i=0;i<Vp->ICS->size;++i)
	{
		pos = Vp->ICS->pos[i];
		ctot[i] = ctotVp[pos] + Vp->ctot_ads[pos];     /*release surface species in bulk, and recalculate */
		if(pos == Vp->pos_X) ctot[i] = Vp->ctot[pos];
	}

	 chem_calc = solve_chem_boundary_node(ctot,call, Vp);
	 

	 if(chem_calc)
	 {
		 for(i=0;i<Vchem->ICS->size;++i)
		 {
			 pos = Vchem->ICS->pos[i];
			 posVp = Vp->ICS->pos[i];


			 Vp->ctot_ads_calc[posVp]=0.;
			 for(j=0;j<Vp->ICS->SM_all->size[0];++j)
			 {
				 if(Vp->ICS->SM_all->type[j] != 0)
					 Vp->ctot_ads_calc[posVp] += Vp->ICS->SM_all->M[j][posVp]*pow10(Vp->ICS->SM_all->log_m[j]);
			 }
			 
			 Vchem->ctot_ads[pos] = Vp->ctot_ads[posVp]=Vp->ctot_ads_calc[posVp];
		 }
	 }
	 /* else use previous calculated ctot_ads, stored in ctot_ads */


	 /* we assume that surface species are always at the end */
	 /*
	for(i=0;i<Vchem->ICS->size;++i) 
	{
		pos   = Vchem->ICS->pos[i];
		ctot_aq[i] = Vchem->ctot_calc[pos] = Vp->ctot_calc[pos]; 
	}
	*/
	 for(i=0;i<Vp->ICS->size;++i) ctotVp[Vp->ICS->pos[i]] = ctot[i];
	 for(i=0;i<Vchem->ICS->size;++i) ctot_aq[i]=ctotVp[Vchem->ICS->pos[i]];

	Vchem->ctot_calc[Vchem->pos_pH] = Vp->ctot_calc[Vp->pos_pH];
	Vchem->ctot[Vchem->pos_pH] = Vp->ctot[Vp->pos_pH];
	Vchem->log_a[Vchem->pos_pH] = Vp->log_a[Vp->pos_pH];
	Vchem->log_m[Vchem->pos_pH] = Vp->log_m[Vp->pos_pH];
	Vchem->log_g[Vchem->pos_pH] = Vp->log_g[Vp->pos_pH];
	if(Vchem->ICS->PRINT_DEBUG_CHEM) writeVchem(Vp, "Vp_init.dat");
	free(ctot);
	free(ctotVp);

	/* STORE SURFACE CHARGE & POTENTIAL */
	Vchem->sch = Vp->sch;
	Vchem->psi = Vp->psi;

	if(Vchem->size_rock>0) solve_chem_boundary_node(ctot_aq,call, Vchem);
}

void set_ctot(real *ctot, struct BasVec *Vchem)
{
	int i;

	for(i=0;i<Vchem->ICS->size;++i) ctot[i] = Vchem->ctot[Vchem->ICS->pos[i]];
}

/*
check_if_mineral_supersaturated:
 checks if mineral is supersaturated, returns the index to the most supersaturated mineral
 OR -1 if no superasturated minerals
 INPUT: 1) list of positions to basis minerals
        2) size of mineral list
 OUTPUT:
 RETURN: index to supersaturated mineral -1 if no minerals were supersaturated
 */
int check_if_mineral_supersaturated(int *mineral_list, int size, int *old_list, int old_size, struct BasVec *Vchem)
{
	int i, pos, pos_b;
	real a,b;


	a=1.e-8;
	pos_b = -1;
	for(i=0;i<size; ++i) /* check all the mineral phases */ 
	{
		pos = mineral_list[i];
		b = Vchem->ICS->SM_mineral->log_a[pos]-Vchem->ICS->SM_mineral->log_af[pos];
		if( b > a && (0>in_list(pos,old_list,old_size)) )   /* also check if mineral not added already */
		{
			pos_b = i;
			a = b;
		}
	}

	return pos_b;
}



/*
add_new_buffer_mineral:
   add new buffer mineral
   INPUT: 
          1) ICSrock: pos to new rock species rel to Vchem
		  2) ICSbuffer: pos to new buffer
		  3) Vchem
   OUTPUT: 1) Updated Vchem
*/

void add_new_buffer_mineral(int ICSrock, int ICSbuffer,   struct BasVec *Vchem)
{
	int *pos_mass;
	int i, idum, flag;

	flag = 0;
	/* remove the basis species from the mass balance */
	pos_mass = (int *) calloc(Vchem->size_mass-1,sizeof(int));
	idum=0;
	for(i=0;i<Vchem->size_mass;i++)
	{
		if(Vchem->pos_mass[i] != ICSrock)
		{
			pos_mass[idum]=Vchem->pos_mass[i];
			idum++;
		}
		else flag = 1;
	}

	if(flag == 1)
	{
		Vchem->size_mass--;
		for(i=0;i<Vchem->size_mass;++i) Vchem->pos_mass[i]=pos_mass[i];
	}
	else
	{
		printf("something is wrong in add_buffer_mineral ....\n");
		printf("basis specie %s is not part of the mass balance\n", Vchem->ICS->SM_basis->row_name[ICSrock]);
		exit(1);
	}
	/* add the basis species to the rock buffered ones */
	Vchem->pos_rock[Vchem->size_rock]   = ICSrock;
	Vchem->pos_buffer[Vchem->size_rock] = ICSbuffer;
	Vchem->size_rock++;

	/* basis transformation matrix */
	update_beta_inv(Vchem);

	/* find reduced stoichiometric matrix */
	for(i=0;i<Vchem->size_rock-1;++i) free(Vchem->sm_buf[i]);
	free(Vchem->sm_buf);

	Vchem->sm_buf = (real **) calloc(Vchem->size_rock,sizeof(real *));
	for(i=0; i< Vchem->size_rock; ++i) Vchem->sm_buf[i] = (real *) calloc(Vchem->size_rock,sizeof(real ));

	update_sm_buf_new(Vchem,Vchem->size_rock, Vchem->pos_buffer, Vchem->pos_rock, Vchem->sm_buf);

	free(pos_mass);
return;
}

/*
remove_buffer_mineral:
   remove buffer mineral
   INPUT: 1) length : position to the buffer mineral that is to be removed
          2) basis_list: list of the position to basis species that are to be used in basis transformation
		  3) size_basis_list
		  4) buffer_list: position of the buffer minerals that are to be used in basis transformation
   OUTPUT: 1) Updated Vchem
*/

void remove_buffer_mineral(int length, struct BasVec *Vchem)
{
	int *pos_rock, *pos_buffer;
	int i, idum, flag, pos_m;

	flag = 0;
	/* remove the basis species from the rock species */
	pos_rock = (int *) calloc(Vchem->size_rock-1,sizeof(int));
	pos_buffer = (int *) calloc(Vchem->size_rock-1,sizeof(int));
	idum=0;
	for(i=0;i<Vchem->size_rock;i++)
	{
		if(i != length)
		{/* copy vector */
			pos_rock[idum]=Vchem->pos_rock[i];
			pos_buffer[idum] = Vchem->pos_buffer[i];
			idum++;
		}
		else
		{
			flag = 1;
			pos_m = Vchem->pos_rock[i];
		}
	}

	if(flag == 1)
	{
		Vchem->size_rock--;
		for(i=0;i<Vchem->size_rock;++i) Vchem->pos_rock[i]=pos_rock[i];
		for(i=0;i<Vchem->size_rock;++i) Vchem->pos_buffer[i]=pos_buffer[i];
	}
	else
	{
		printf("something is wrong in remove_buffer_mineral ....\n");
		printf("basis specie %s is not part of the rock species\n", Vchem->ICS->SM_basis->row_name[Vchem->pos_buffer[length]]);
		exit(1);
	}

	/* add the basis species to the ones determined by mass balance*/
	Vchem->pos_mass[Vchem->size_mass]   = pos_m;
	Vchem->size_mass++;

	/* basis transformation matrix */
	update_beta_inv(Vchem);

	/* find reduced stoichiometric matrix */
	for(i=0;i<Vchem->size_rock+1;++i) free(Vchem->sm_buf[i]);
	free(Vchem->sm_buf);

	Vchem->sm_buf = (real **) calloc(Vchem->size_rock,sizeof(real *));
	for(i=0; i< Vchem->size_rock; ++i) Vchem->sm_buf[i] = (real *) calloc(Vchem->size_rock,sizeof(real ));

	update_sm_buf_new(Vchem,Vchem->size_rock, Vchem->pos_buffer, Vchem->pos_rock, Vchem->sm_buf);

	free(pos_rock);free(pos_buffer);
return;
}

/*
check_if_mineral_conc_neg:
 checks if mineral has a negative conc, returns the index to the most neg mineral
 OR -1 if no minerals has a neg conc
 INPUT: 1) list of positions to basis minerals
        2) size of mineral list
 OUTPUT:
 RETURN: index to supersaturated mineral -1 if no minerals were supersaturated
 */
int check_if_mineral_conc_neg(int *list, int size, struct BasVec *Vchem)
{
	int i, pos, pos_b;
	real a,b;

	a=1.e-32;
	pos_b = -1;
	for(i=0;i<size; ++i) /* check all the mineral phases */ 
	{
		pos = list[i];
		b = Vchem->ctot_mineral[pos];
		if( b < a )
		{
			pos_b = i;
			a = b;
		}
	}

	return pos_b;
}

/*
add_new_sup_mineral:
   add new buffer mineral
   INPUT: 1) length : position to the new rock buffer mineral
          2) basis_list: list of the position to basis species that are to be used in basis transformation
		  3) size_basis_list
		  4) buffer_list: position of the buffer minerals that are to be used in basis transformation
   OUTPUT: 1) Updated Vchem
*/

void add_new_sup_mineral(int length, struct InitChem *ICS,   struct BasVec *Vchem)
{
	//int *pos_mass;
	int i, j, flag;
	//real nu;

	flag = 0;
	

	/* add the basis species to the rock buffered ones */
	Vchem->pos_sup_min[Vchem->size_sup_min]   = ICS->pos_sup_min[length];
	Vchem->pos_sup_bas[Vchem->size_sup_min]   = ICS->pos_sup_bas[length];

	/* default values for generic form  rate = (k_1+k_2*a_h)*(1-SI^k_3)^k_4 */
	Vchem->rate[Vchem->size_sup_min][0] = ICS->rate[length][0];
	Vchem->rate[Vchem->size_sup_min][1] = ICS->rate[length][1];
	Vchem->rate[Vchem->size_sup_min][2] = ICS->rate[length][2];
	Vchem->rate[Vchem->size_sup_min][3] = ICS->rate[length][3];

	printf("Sup min: %s with %s\n",
		ICS->SM_basis->row_name[ICS->pos_sup_bas[length]],
		ICS->SM_mineral->row_name[ICS->pos_sup_min[length]]);
							
	Vchem->size_sup_min++;

	


	/* find reduced stoichiometric matrix */
	for(j=0;j<Vchem->size_sup_min-1;++j)
	{
		free(Vchem->sm_sup_buf[j]);
	}

	free(Vchem->sm_sup_buf);

	Vchem->sm_sup_buf = (real **) calloc(Vchem->size_sup_min,sizeof(real *));
	for(j=0; j< Vchem->size_sup_min; ++j) Vchem->sm_sup_buf[j] = (real *) calloc(Vchem->size_sup_min,sizeof(real ));

	update_sm_buf_new(Vchem,Vchem->size_sup_min, Vchem->pos_sup_min, Vchem->pos_sup_bas, Vchem->sm_sup_buf);

	for(i=0;i<Vchem->size_sup_min;++i)
	{
		for(j=0;j<Vchem->size_sup_min;++j)
			printf("%g\t",Vchem->sm_sup_buf[i][j]);
		printf("\n");
	}
return;
}
/* length = pos in pos_sup_bas i.e. pos_sup_bas[length] has to be removed */
void remove_sup_mineral(int length,   struct BasVec *Vchem)
{
	int *pos_rock, *pos_buffer;
	int i, idum, flag, pos_m;

	flag = 0;
	/* remove the basis species from the rock species */
	pos_rock = (int *) calloc(Vchem->size_sup_min-1,sizeof(int));
	pos_buffer = (int *) calloc(Vchem->size_sup_min-1,sizeof(int));
	idum=0;
	for(i=0;i<Vchem->size_sup_min;i++)
	{
		if(i != length)
		{/* copy vector */
			pos_rock[idum]=Vchem->pos_sup_bas[i];
			pos_buffer[idum] = Vchem->pos_sup_min[i];
			idum++;
		}
		else
		{
			flag = 1;
			pos_m = Vchem->pos_sup_bas[i];
		}
	}

	if(flag == 1)
	{
		Vchem->size_sup_min--;
		for(i=0;i<Vchem->size_sup_min;++i) Vchem->pos_sup_bas[i]=pos_rock[i];
		for(i=0;i<Vchem->size_sup_min;++i) Vchem->pos_sup_min[i]=pos_buffer[i];
	}
	else
	{
		printf("something is wrong in remove_buffer_mineral ....\n");
		printf("basis specie %s is not part of the rock species\n", Vchem->ICS->SM_basis->row_name[Vchem->pos_sup_min[length]]);
		exit(1);
	}

	/* find reduced stoichiometric matrix */
	for(i=0;i<Vchem->size_sup_min+1;++i) free(Vchem->sm_sup_buf[i]);
	free(Vchem->sm_sup_buf);

	Vchem->sm_sup_buf = (real **) calloc(Vchem->size_sup_min,sizeof(real *));
	for(i=0; i< Vchem->size_sup_min; ++i) Vchem->sm_sup_buf[i] = (real *) calloc(Vchem->size_sup_min,sizeof(real ));

	update_sm_buf_new(Vchem,Vchem->size_sup_min, Vchem->pos_sup_min, Vchem->pos_sup_bas, Vchem->sm_sup_buf);

	free(pos_rock);free(pos_buffer);
return;
}


/*call_chem_solver:
INPUT  : ctot - new concentrations
RETURN : 0 if |log_10(ctot_i)-log_10(ctot_old_i)| < CHEM_NEW_CALC
*/
int call_chem_solver(real *ctot, struct InitChem *ICS, struct BasVec *Vchem)
{
	int i;
	int check, pos;
	check = 0;

	for(i=0;i<ICS->size;++i)
	{
		pos = ICS->pos[i];
		if( fabs(log10(ctot[i])-log10(Vchem->ctot[pos])) > CHEM_NEW_CALC_)
		{
			check = 1;
			break;
		}
	}
	
/*	if(check == 0)
	{
		printf("Skip chem solver\n");
	}*/

	return 1;
}


void calc_adsorption(struct BasVec *Vchem_surf, struct BasVec *Vchem)
{
	int i, j;

	for(i=0;i<Vchem->size;++i)
	{
		Vchem_surf->log_a[i] = Vchem->log_a[i];
		Vchem_surf->log_m[i] = Vchem->log_m[i];
		Vchem_surf->log_g[i] = Vchem->log_g[i];
		Vchem_surf->ctot[i]  = Vchem->ctot[i];
	}

	/* add adsorbed species to ctot: */
	for(i=0;i<Vchem->size;++i) Vchem_surf->ctot[i] += Vchem->ctot_ads[i];

	if(Vchem_surf->log_a[Vchem_surf->pos_X] > -15) Vchem_surf->log_a[Vchem_surf->pos_X] = Vchem_surf->log_m[Vchem_surf->pos_X] = -20.;
	Vchem_surf->Io = Vchem->Io;
	Vchem_surf->psi = Vchem->psi;
	Vchem_surf->sch = Vchem->sch;
	
	solve_chemistry(Vchem_surf);

	Vchem->psi = Vchem_surf->psi;
	Vchem->sch = Vchem_surf->sch;
	Vchem->log_a[Vchem->pos_exp] = Vchem_surf->log_a[Vchem_surf->pos_exp];
	Vchem->log_m[Vchem->pos_exp] = Vchem_surf->log_m[Vchem_surf->pos_exp];
	Vchem->log_a[Vchem->pos_X] = Vchem_surf->log_a[Vchem->pos_X];
	Vchem->log_m[Vchem->pos_X] = Vchem_surf->log_m[Vchem->pos_X];
		/* adsorbed species concentration */
	for(i=0;i<Vchem->size;++i)
	{
		Vchem->ctot_ads[i]=0.;
		for(j=0;j<Vchem->ICS->SM_all->size[0];++j)
		{
			if(Vchem->ICS->SM_all->type[j] != 0)/*surface species*/
				Vchem->ctot_ads[i] += Vchem->ICS->SM_all->M[j][i]*pow10(Vchem->ICS->SM_all->log_m[j]);
		}
	}

}

/* convert species determined by the rock to mass balance */ 
void convert_rock_to_mass(struct BasVec *Vchem_surf, struct BasVec *Vchem)
{
	int i, length;

	Vchem_surf->pos_pH    = Vchem->pos_pH;
	Vchem_surf->pos_exp   = Vchem->pos_exp;
	Vchem_surf->pos_water = Vchem->pos_water;
	Vchem_surf->pos_X     = Vchem->pos_X;
	Vchem_surf->pos_pH    = Vchem->pos_pH;

	length = 0;
	for(i=0;i<Vchem->size;++i)
	{
		if(i == Vchem->pos_mass[i])      Vchem_surf->pos_mass[length] = i;
		else if(i == Vchem->pos_rock[i]) Vchem_surf->pos_mass[length] = i;
		length++;
	}
	Vchem_surf->size_mass = length;
	Vchem_surf->size_rock = 0;
	Vchem_surf->size      = Vchem->size;
	if(length != Vchem->size_mass + Vchem->size_rock)
	{
		printf("Warning: something is when converting from rock species to mass species\n");
	}
}

void set_temperature_db(struct BasVec *Vchem)
{
	struct ChemTable *SM_mineral, *SM_all, *SM_aq_logK, *SM_M_logK;
	real CAdh[]={0.51 , -1.154, 35.697, -182.023, 346.528},
		 CBdh[]={0.325, 0.08  , 1.441 , -6.541  , 11.655};
	int j;

	SM_mineral = Vchem->ICS->SM_mineral;
	SM_all     = Vchem->ICS->SM_all;
	SM_aq_logK = Vchem->ICS->SM_aq_logK;
	SM_M_logK  = Vchem->ICS->SM_M_logK;


	/* get correct logK values */
	get_logK(Vchem->Temp, SM_all->logK,  SM_aq_logK);
	get_logK(Vchem->Temp, SM_mineral->logK,  SM_M_logK);
	/* Debye-Huckel */
	Vchem->ICS->CP->Adh=Vchem->ICS->CP->Bdh=0.;

	for(j=0;j<=4;++j)
	{
		Vchem->ICS->CP->Adh += pow(0.001*(Vchem->Temp-273.15),j)*CAdh[j];
		Vchem->ICS->CP->Bdh += pow(0.001*(Vchem->Temp-273.15),j)*CBdh[j];
	}
	Vchem->ICS->CP->beta_chem = Vchem->ICS->CP->F/(Vchem->ICS->CP->Rg*Vchem->Temp);
}

void remove_mineral_buffer(int cw, struct BasVec *Vp)
{
	int length;
	struct BasVec *Vchem;
	Vchem = &Vp[cw];
	
	do
	{
		length = check_if_mineral_conc_neg(Vchem->pos_buffer,Vchem->size_rock,Vchem);
		if(length >= 0)
		{
		  /*
		    printf("Buffer mineral %s exhausted\n", Vchem->ICS->SM_mineral->row_name[Vchem->pos_buffer[length]]);*/
			remove_buffer_mineral(length, Vchem);
		}

		/*solve_chemistry(&mass_balance, Vchem);*/

	}while(length>0);
}
