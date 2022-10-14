#include "chem_global.h"
#include "mpi.h"

/* log increments */
void newton_ah_log(real *x, real *Fchem, real **jacobi, struct BasVec *Vchem)
{
  
  int i,j,pos,q, pos_q, p, pos_i;
  real  g_tmp_j, g_tmp_q, g_tmp_p;
  real *delta_calc, **jac_dum, *mpow, *smpow;

  struct ChemTable *SM_mineral, *SM_all, *SM_basis;
  SM_mineral = Vchem->ICS->SM_mineral;
  SM_all     = Vchem->ICS->SM_all;
  SM_basis   = Vchem->ICS->SM_basis;
  
  delta_calc = Vchem->ICS->mem->dVmb;
  jac_dum    = Vchem->ICS->mem->ddVmb;
  mpow       = Vchem->ICS->mem->dVmc;
  smpow      = Vchem->ICS->mem->dSMa;
  
  // REMOVE
  /* if (_DEBUG_FLAG_) { */
  /*   printf("Jacobi before:\n"); */
  /*   for (i=0; i<Vchem->size_mass; i++) { */
  /*     for (j=0; j<Vchem->size_mass; j++) { */
  /* 	printf("%g, ", jacobi[i+1][j+1]); */
  /*     } */
  /*     printf("\n"); */
  /*   } */
  /*   printf("\n"); */
  /* } */
  // REMOVE

  for (i = 0;i < Vchem->size_mass;++i) /* construct basis vector */
    {
      pos = Vchem->pos_mass[i];
      Vchem->log_m[pos] = (real) x[i];
      Vchem->log_a[pos] = Vchem->log_m[pos]+Vchem->log_g[pos];
      if(SM_basis->type[pos] == 2) Vchem->log_a[pos] += log10(Vchem->ctot[Vchem->pos_X]);
    }
	
  if(Vchem->size_rock>0) calc_rock_spec_conc(Vchem); 
	
  /* calculate m = 10^log_m */
  for(i=0;i<Vchem->size;++i) mpow[i] = pow10(Vchem->log_m[i]);
  /* calculate concentration of complexes */
  calc_complex(SM_all, Vchem);

  /* calculate SM_all_m = 10^log_SM_all_m */
  for(i=0;i<SM_all->size[0];++i) smpow[i] = pow10(SM_all->log_m[i]);

  /* amount of rock species precipitated/dissolved                                  */
  /* delta_calc > 0 => species conc decreases in bulk and increases at the surface  */
  for(j=0; j < Vchem->size_rock ; ++j)
    {
      pos  = Vchem->pos_rock[j];
      delta_calc[j]= - mpow[pos];
      for(i=0; i< SM_all->size[0]; ++i)
	{
	  delta_calc[j] -= SM_all->M[i][pos]*smpow[i];
	}
      delta_calc[j] += Vchem->ctot[pos];
    }

  for(j=0;j<Vchem->size_rock;++j)
    {
      pos = Vchem->pos_buffer[j];
      SM_mineral->delta[pos]=0.;
      for(i=0;i<Vchem->size_rock;++i)
	{
	  SM_mineral->delta[pos]+= Vchem->sm_buf[j][i]*delta_calc[i];
	}
    }

  /* calculate F_j*/ 
  /* F_j = c_tot_j - m_j - sum_i^N_rb M_ij delta_i-sum_i^Nx SM_i,j x_i */
  Fchem[0]=0.;
  for(j=0; j < Vchem->size_mass ; ++j)
    {
      pos  = Vchem->pos_mass[j];
      Fchem[j+1] = Vchem->ctot[pos]-mpow[pos];
      if(pos==Vchem->pos_pH) Fchem[j+1] -= Vchem->WHTO;
      for(i=0; i< SM_all->size[0]; ++i)
	{
	  Fchem[j+1] -= SM_all->M[i][pos]*smpow[i];
	}
      /* add amount lost or gained due to precipitation/dissolution */
      for(i=0;i<Vchem->size_rock;++i)
	{
	  pos_q=Vchem->pos_buffer[i];
	  Fchem[j+1] -= SM_mineral->M[pos_q][pos]*SM_mineral->delta[pos_q];
	}

    }
  if(Vchem->size_sup_min>0)
    {
      calc_F_sup_min(Fchem,Vchem);
    }


  /* and the jacobi */
  /* jacobi_j,q = dF_j/dm_q = -delta_j,q m_j - ( sum_i^Nx SM_i,q * SM_i,j * x_i ) */ 
  for(j=0; j <Vchem->size_mass ; ++j)
    {
      pos  = Vchem->pos_mass[j];
      g_tmp_j =0.;

      for(q=0; q< Vchem->size_mass; ++q)
	{
	  g_tmp_q = 0.;
	  pos_q  = Vchem->pos_mass[q];
	  jacobi[j+1][q+1]=-Vchem->beta_bas_inv[pos][pos_q]*mpow[pos_q];
	  for(i=0; i < SM_all->size[0]; ++i)
	    {
	      g_tmp_p = 0.;
	      for(p=0;p<Vchem->size;++p)/* basis transformation */
		{
		  g_tmp_p += SM_all->M[i][p]*Vchem->beta_bas_inv[p][pos_q];
		}

	      g_tmp_q += g_tmp_p*SM_all->M[i][pos]*smpow[i];
	    }
	  jacobi[j+1][q+1] -= g_tmp_q;
	  jacobi[j+1][q+1]*= LNTEN;
	}

    }

  // REMOVE
  /* if (_DEBUG_FLAG_) { */
  /*   printf("Jacobi middle:\n"); */
  /*   for (i=0; i<Vchem->size_mass; i++) { */
  /*     for (j=0; j<Vchem->size_mass; j++) { */
  /* 	printf("%g, ", jacobi[i+1][j+1]); */
  /*     } */
  /*     printf("\n"); */
  /*   } */
  /*   printf("\n"); */
  /* } */
  // REMOVE


  /* the rock buffer part of the jacobian */ 
  for(j=0; j <Vchem->size_rock ; ++j)
    {
      pos  = Vchem->pos_rock[j];
      g_tmp_j =0.;

      for(q=0; q< Vchem->size_mass; ++q)
	{
	  g_tmp_q = 0.;
	  pos_q  = Vchem->pos_mass[q];
	  jac_dum[j+1][q+1] =-Vchem->beta_bas_inv[pos][pos_q]*mpow[pos];		 
	  for(i=0; i < SM_all->size[0]; ++i)
	    {
	      g_tmp_p = 0.;
	      for(p=0;p<Vchem->size;++p)/* basis transformation */
		{
		  g_tmp_p += SM_all->M[i][p]*Vchem->beta_bas_inv[p][pos_q];
		}

	      g_tmp_q += g_tmp_p*SM_all->M[i][pos]*smpow[i];
	    }
	  jac_dum[j+1][q+1] -= g_tmp_q;
	}

    }

  for(j=0; j <Vchem->size_mass ; ++j)
    {
      pos  = Vchem->pos_mass[j];
      for(q=0; q< Vchem->size_mass; ++q)
	{
	  g_tmp_q = 0.;
	  for(i=0;i<Vchem->size_rock;++i)
	    {
	      pos_i = Vchem->pos_buffer[i];
	      g_tmp_p = 0.;
	      for(p=0; p<Vchem->size_rock;++p)
		{
		  g_tmp_p += Vchem->sm_buf[i][p]*jac_dum[p+1][q+1];
		}
	      g_tmp_q += SM_mineral->M[pos_i][pos]*g_tmp_p;
	    }
	  jacobi[j+1][q+1] -= g_tmp_q*LNTEN;
	}
    }

  /* SUP SAT MINERALS */
  if(Vchem->size_sup_min>0)
    {
      calc_jacobi_sup_min(jacobi,Vchem);
    }
	

  // REMOVE
  /* if (_DEBUG_FLAG_) { */
  /*   printf("Jacobi end:\n"); */
  /*   for (i=0; i<Vchem->size_mass; i++) { */
  /*     for (j=0; j<Vchem->size_mass; j++) { */
  /* 	printf("%g, ", jacobi[i+1][j+1]); */
  /*     } */
  /*     printf("\n"); */
  /*   } */
  /*   printf("\n"); */
  /* } */
  // REMOVE

  /*-----------------------*/
  /*debug*/
  /*
    j_num = (real **) calloc(Vchem->size_mass+1, sizeof(real *));
    for(i=0;i<Vchem->size_mass+1;++i) j_num[i] = (real *) calloc(Vchem->size_mass+1, sizeof(real ));
    calc_complex(SM_mineral,Vchem);
    jacobi_num(x, j_num, Vchem);

    fp=my_fopen("debug_jacobi.out", "w");
    fprintf(fp,"jacobi num\n");
    for(i=0;i<Vchem->size_mass;++i)
    {
    for(j=0;j<Vchem->size_mass;++j)
    {
    fprintf(fp,"%4.8e\t",j_num[i+1][j+1]);
    }
    fprintf(fp,"\n");
    }
    fprintf(fp,"\n");
    fprintf(fp,"jacobi analytic\n");
    for(i=0;i<Vchem->size_mass;++i)
    {
    for(j=0;j<Vchem->size_mass;++j)
    {
    fprintf(fp,"%4.8e\t",jacobi[i+1][j+1]);
    }
    fprintf(fp,"\n");
    }
    fprintf(fp,"\n");
    fclose(fp);

    for(i=0;i<Vchem->size_mass;++i)
    for(j=0;j<Vchem->size_mass;++j) jacobi[i+1][j+1]=j_num[i+1][j+1];
		

    for(i=0;i<Vchem->size_mass;++i) free(j_num[i]);
    free(j_num);
  */
  /*-------------------------------------------------------*/	


	
  return;
}

/* log increments */
/* includes equilebration with ion exchanger*/
void newton_ah_log_io(real *x, real *Fchem, real **jacobi, struct BasVec *Vchem)
{
  
    int i,j,pos,q, pos_q, p, pos_i;
    real  g_tmp_j, g_tmp_q, g_tmp_p;
	real *delta_calc, **jac_dum, *mpow, *smpow;

	struct ChemTable *SM_mineral, *SM_all, *SM_basis;
	SM_mineral = Vchem->ICS->SM_mineral;
	SM_all     = Vchem->ICS->SM_all;
	SM_basis   = Vchem->ICS->SM_basis;

	delta_calc = Vchem->ICS->mem->dVmb;
	jac_dum    = Vchem->ICS->mem->ddVmb;
	mpow       = Vchem->ICS->mem->dVmc;
	smpow      = Vchem->ICS->mem->dSMa;

    for (i = 0;i < Vchem->size_mass;++i) /* construct basis vector */
	{
		pos = Vchem->pos_mass[i];
		Vchem->log_m[pos] = (real) x[i];
		Vchem->log_a[pos] = Vchem->log_m[pos]+Vchem->log_g[pos];
		if(SM_basis->type[pos] == 2) Vchem->log_a[pos] += log10(Vchem->ctot[Vchem->pos_X]);
	}
	
	if(Vchem->size_rock>0) calc_rock_spec_conc(Vchem); 
	
	/* calculate m = 10^log_m */
	for(i=0;i<Vchem->size;++i) mpow[i] = pow10(Vchem->log_m[i]);
	/* calculate concentration of complexes */
	calc_complex(SM_all, Vchem);

	/* calculate SM_all_m = 10^log_SM_all_m */
	for(i=0;i<SM_all->size[0];++i) smpow[i] = pow10(SM_all->log_m[i]);

	/* amount of rock species precipitated/dissolved                                  */
	/* delta_calc > 0 => species conc decreases in bulk and increases at the surface  */
	for(j=0; j < Vchem->size_rock ; ++j)
	{
		pos  = Vchem->pos_rock[j];
		delta_calc[j]= - mpow[pos];
		for(i=0; i< SM_all->size[0]; ++i)
		{
			if(SM_all->type[i] != 2)
			delta_calc[j] -= SM_all->M[i][pos]*smpow[i];
		}
		delta_calc[j] += Vchem->ctot[pos];
	}
	

	for(j=0;j<Vchem->size_rock;++j)
	{
		pos = Vchem->pos_buffer[j];
		SM_mineral->delta[pos]=0.;
		for(i=0;i<Vchem->size_rock;++i)
		{
			SM_mineral->delta[pos]+= Vchem->sm_buf[j][i]*delta_calc[i];
		}
	}

	/* calculate F_j*/ 
	/* F_j = c_tot_j - m_j - sum_i^N_rb M_ij delta_i-sum_i^Nx SM_i,j x_i */
	Fchem[0]=0.;
	for(j=0; j < Vchem->size_mass ; ++j)
	{
		pos  = Vchem->pos_mass[j];
		Fchem[j+1] = Vchem->ctot[pos]-mpow[pos];
		if(pos==Vchem->pos_pH) Fchem[j+1] -= Vchem->WHTO;
		for(i=0; i< SM_all->size[0]; ++i)
		{
			/* assume that c_tot is the bulk concentration and does not include surface species */
			/* mass balance still valid for surface species */
			if(pos == Vchem->pos_X || SM_all->type[i] != 2) Fchem[j+1] -= SM_all->M[i][pos]*smpow[i];
		}
			/* add amount lost or gained due to precipitation/dissolution */
		for(i=0;i<Vchem->size_rock;++i)
		{
			pos_q=Vchem->pos_buffer[i];
		    Fchem[j+1] -= SM_mineral->M[pos_q][pos]*SM_mineral->delta[pos_q];
		}

	}
	
	if(Vchem->size_sup_min>0)
	{
		calc_F_sup_min(Fchem,Vchem);
	}


	/* and the jacobi */
	/* jacobi_j,q = dF_j/dm_q = -delta_j,q m_j - ( sum_i^Nx SM_i,q * SM_i,j * x_i ) */ 
	for(j=0; j <Vchem->size_mass ; ++j)
	{
		pos  = Vchem->pos_mass[j];
		g_tmp_j =0.;

		for(q=0; q< Vchem->size_mass; ++q)
		{
			g_tmp_q = 0.;
			pos_q  = Vchem->pos_mass[q];
			jacobi[j+1][q+1]=-Vchem->beta_bas_inv[pos][pos_q]*mpow[pos_q];
			for(i=0; i < SM_all->size[0]; ++i)
			{
				g_tmp_p = 0.;
				for(p=0;p<Vchem->size;++p)/* basis transformation */
				{
					g_tmp_p += SM_all->M[i][p]*Vchem->beta_bas_inv[p][pos_q];
				}

				if(pos == Vchem->pos_X || SM_all->type[i] != 2) g_tmp_q += g_tmp_p*SM_all->M[i][pos]*smpow[i];
			}
			jacobi[j+1][q+1] -= g_tmp_q;
			jacobi[j+1][q+1]*= LNTEN;
		}

	}
/* the rock buffer part of the jacobian */ 
	for(j=0; j <Vchem->size_rock ; ++j)
	{
		pos  = Vchem->pos_rock[j];
		g_tmp_j =0.;

		for(q=0; q< Vchem->size_mass; ++q)
		{
			g_tmp_q = 0.;
			pos_q  = Vchem->pos_mass[q];
			jac_dum[j+1][q+1] =-Vchem->beta_bas_inv[pos][pos_q]*mpow[pos];		 
			for(i=0; i < SM_all->size[0]; ++i)
			{
				g_tmp_p = 0.;
				for(p=0;p<Vchem->size;++p)/* basis transformation */
				{
					if(SM_all->type[i] != 2)
					g_tmp_p += SM_all->M[i][p]*Vchem->beta_bas_inv[p][pos_q];
				}

				g_tmp_q += g_tmp_p*SM_all->M[i][pos]*smpow[i];
			}
			jac_dum[j+1][q+1] -= g_tmp_q;
		}

	}

	for(j=0; j <Vchem->size_mass ; ++j)
	{
		pos  = Vchem->pos_mass[j];
		for(q=0; q< Vchem->size_mass; ++q)
		{
			g_tmp_q = 0.;
			for(i=0;i<Vchem->size_rock;++i)
			{
				pos_i = Vchem->pos_buffer[i];
				g_tmp_p = 0.;
				for(p=0; p<Vchem->size_rock;++p)
				{
					g_tmp_p += Vchem->sm_buf[i][p]*jac_dum[p+1][q+1];
				}
				g_tmp_q += SM_mineral->M[pos_i][pos]*g_tmp_p;
			}
			jacobi[j+1][q+1] -= g_tmp_q*LNTEN;
		}
	}

	/* SUP SAT MINERALS */
	if(Vchem->size_sup_min>0)
	{
		calc_jacobi_sup_min(jacobi,Vchem);
	}
	

	/*-----------------------*/
	/*debug*/
	/*
	j_num = (real **) calloc(Vchem->size_mass+1, sizeof(real *));
	for(i=0;i<Vchem->size_mass+1;++i) j_num[i] = (real *) calloc(Vchem->size_mass+1, sizeof(real ));
	calc_complex(&SM_mineral,Vchem);
	jacobi_num(x, j_num, Vchem);

	fp=my_fopen("debug_jacobi.out", "w");
	fprintf(fp,"jacobi num\n");
	for(i=0;i<Vchem->size_mass;++i)
	{
		for(j=0;j<Vchem->size_mass;++j)
		{
			fprintf(fp,"%4.8e\t%",j_num[i+1][j+1]);
		}
		fprintf(fp,"\n");
	}
	fprintf(fp,"\n");
	fprintf(fp,"jacobi analytic\n");
	for(i=0;i<Vchem->size_mass;++i)
	{
		for(j=0;j<Vchem->size_mass;++j)
		{
			fprintf(fp,"%4.8e\t%",jacobi[i+1][j+1]);
		}
		fprintf(fp,"\n");
	}
	fprintf(fp,"\n");
	fclose(fp);

	for(i=0;i<Vchem->size_mass;++i)
		for(j=0;j<Vchem->size_mass;++j) jacobi[i+1][j+1]=j_num[i+1][j+1];

	for(i=0;i<Vchem->size_mass;++i) free(j_num[i]);
	free(j_num);
	*/
/*-------------------------------------------------------*/	

	
return;
}
/* log increments */
/* add surface charge to the system, solves the Grahame equation */
/* uses log_10 E = F psi/(R T) as the variable */
void newton_ah_log_surf(real *x, real *Fchem, real **jacobi, struct BasVec *Vchem)
{
  
    int i,j,pos,q, pos_q, p, pos_i;
    real  g_tmp_j, g_tmp_q, g_tmp_p;
	real *delta_calc, **jac_dum, *mpow, *smpow;
	real kappa_s, E, g_tmp_1, g_tmp_2, g_tmp_3;

	struct ChemTable *SM_mineral, *SM_all, *SM_basis;
	SM_mineral = Vchem->ICS->SM_mineral;
	SM_all     = Vchem->ICS->SM_all;
	SM_basis   = Vchem->ICS->SM_basis;

	delta_calc = Vchem->ICS->mem->dVmb;
	jac_dum    = Vchem->ICS->mem->ddVmb;
	mpow       = Vchem->ICS->mem->dVmc;
	smpow      = Vchem->ICS->mem->dSMa;

	

/*	kappa_s = 1.778e-5*Vchem->Temp; BUG!!!!!!!!*/
	kappa_s = Vchem->ICS->CP->Rg*Vchem->ICS->CP->e0*Vchem->ICS->CP->ew*2000.;
	kappa_s *= Vchem->Temp;

    for (i = 0;i < Vchem->size_mass;++i) /* construct basis vector */
	{
		pos = Vchem->pos_mass[i];
		Vchem->log_m[pos] = (real) x[i];
		Vchem->log_a[pos] = Vchem->log_m[pos]+Vchem->log_g[pos];
	}

	Vchem->log_a[Vchem->pos_exp]=Vchem->log_m[Vchem->pos_exp] = (real) x[Vchem->size_mass];
	Vchem->psi = Vchem->log_a[Vchem->pos_exp]*LNTEN/Vchem->ICS->CP->beta_chem;
	E = pow10(Vchem->log_m[Vchem->pos_exp]);
	/*
	Vchem->psi = Vchem->log_a[Vchem->pos_exp]*LNTEN/Vchem->ICS->CP->beta_chem;
	E = pow(10,Vchem->log_m[Vchem->pos_exp]);
	printf("%g %g\n",E, exp(Vchem->ICS->CP->beta_chem*Vchem->psi));
	*/
	
	if(Vchem->size_rock>0) calc_rock_spec_conc(Vchem); 
	
	/* calculate m = 10^log_m */
	for(i=0;i<Vchem->size;++i) mpow[i] = pow10(Vchem->log_m[i]);
	/* calculate concentration of complexes */
	calc_complex(SM_all, Vchem);

	/* calculate SM_all_m = 10^log_SM_all_m */
	for(i=0;i<SM_all->size[0];++i) smpow[i] = pow10(SM_all->log_m[i]);

	/* amount of rock species precipitated/dissolved                                  */
	/* delta_calc > 0 => species conc decreases in bulk and increases at the surface  */
	for(j=0; j < Vchem->size_rock ; ++j)
	{
		pos  = Vchem->pos_rock[j];
		delta_calc[j]= - mpow[pos];
		for(i=0; i< SM_all->size[0]; ++i)
		{
			delta_calc[j] -= SM_all->M[i][pos]*smpow[i];
		}
		delta_calc[j] += Vchem->ctot[pos];
	}

	for(j=0;j<Vchem->size_rock;++j)
	{
		pos = Vchem->pos_buffer[j];
		SM_mineral->delta[pos]=0.;
		for(i=0;i<Vchem->size_rock;++i)
		{
			SM_mineral->delta[pos]+= Vchem->sm_buf[j][i]*delta_calc[i];
		}
	}

	/* calculate F_j*/ 
	/* F_j = c_tot_j - m_j - sum_i^N_rb M_ij delta_i-sum_i^Nx SM_i,j x_i */
	Fchem[0]=0.;
	for(j=0; j < Vchem->size_mass ; ++j)
	{
		pos  = Vchem->pos_mass[j];
		Fchem[j+1] = Vchem->ctot[pos]-mpow[pos];
		if(Vchem->DL) Fchem[j+1] -= Vchem->ctot_dl[pos]; /*diffuse layer concentration*/
		for(i=0; i< SM_all->size[0]; ++i)
		{
			Fchem[j+1] -= SM_all->M[i][pos]*smpow[i];
		}
			/* add amount lost or gained due to precipitation/dissolution */
		for(i=0;i<Vchem->size_rock;++i)
		{
			pos_q=Vchem->pos_buffer[i];
			Fchem[j+1] -= SM_mineral->M[pos_q][pos]*SM_mineral->delta[pos_q];
		}

	}
	
	if(Vchem->size_sup_min>0)
	{
		calc_F_sup_min(Fchem,Vchem);
	}


	/* and the jacobi */
	/* jacobi_j,q = dF_j/dm_q = -delta_j,q m_j - ( sum_i^Nx SM_i,q * SM_i,j * x_i ) */ 
	for(j=0; j <Vchem->size_mass ; ++j)
	{
		pos  = Vchem->pos_mass[j];
		g_tmp_j =0.;

		for(q=0; q< Vchem->size_mass; ++q)
		{
			g_tmp_q = 0.;
			pos_q  = Vchem->pos_mass[q];
			jacobi[j+1][q+1]=-Vchem->beta_bas_inv[pos][pos_q]*mpow[pos_q];
			for(i=0; i < SM_all->size[0]; ++i)
			{
				g_tmp_p = 0.;
				for(p=0;p<Vchem->size;++p)/* basis transformation */
				{
					g_tmp_p += SM_all->M[i][p]*Vchem->beta_bas_inv[p][pos_q];
				}

				g_tmp_q += g_tmp_p*SM_all->M[i][pos]*smpow[i];
			}
			jacobi[j+1][q+1] -= g_tmp_q;
			jacobi[j+1][q+1]*= LNTEN;
		}

	}

	/* ---------- add E - species ----------*/

	for(j=0; j <Vchem->size_mass ; ++j)
	{
		pos  = Vchem->pos_mass[j];
		g_tmp_j =0.;

		q = Vchem->size_mass;
		
		g_tmp_q = 0.;
		pos_q  = Vchem->pos_exp;
		jacobi[j+1][q+1]=-Vchem->beta_bas_inv[pos][pos_q]*mpow[pos_q];
		for(i=0; i < SM_all->size[0]; ++i)
		{
			g_tmp_p = 0.;
			for(p=0;p<Vchem->size;++p)/* basis transformation */
			{
				g_tmp_p += SM_all->M[i][p]*Vchem->beta_bas_inv[p][pos_q];
			}

			g_tmp_q += g_tmp_p*SM_all->M[i][pos]*smpow[i];
		}
		jacobi[j+1][q+1] -= g_tmp_q;
		jacobi[j+1][q+1]*= LNTEN;
	}
/* ----------------------------------------------------*/
/* the rock buffer part of the jacobian */ 
	for(j=0; j <Vchem->size_rock ; ++j)
	{
		pos  = Vchem->pos_rock[j];
		g_tmp_j =0.;

		for(q=0; q< Vchem->size_mass; ++q)
		{
			g_tmp_q = 0.;
			pos_q  = Vchem->pos_mass[q];
			jac_dum[j+1][q+1] =-Vchem->beta_bas_inv[pos][pos_q]*mpow[pos];		 
			for(i=0; i < SM_all->size[0]; ++i)
			{
				g_tmp_p = 0.;
				for(p=0;p<Vchem->size;++p)/* basis transformation */
				{
					g_tmp_p += SM_all->M[i][p]*Vchem->beta_bas_inv[p][pos_q];
				}

				g_tmp_q += g_tmp_p*SM_all->M[i][pos]*smpow[i];
			}
			jac_dum[j+1][q+1] -= g_tmp_q;
		}

	}

	for(j=0; j <Vchem->size_mass ; ++j)
	{
		pos  = Vchem->pos_mass[j];
		for(q=0; q< Vchem->size_mass; ++q)
		{
			g_tmp_q = 0.;
			for(i=0;i<Vchem->size_rock;++i)
			{
				pos_i = Vchem->pos_buffer[i];
				g_tmp_p = 0.;
				for(p=0; p<Vchem->size_rock;++p)
				{
					g_tmp_p += Vchem->sm_buf[i][p]*jac_dum[p+1][q+1];
				}
				g_tmp_q += SM_mineral->M[pos_i][pos]*g_tmp_p;
			}
			jacobi[j+1][q+1] -= g_tmp_q*LNTEN;
		}
	}

	/* Surface charge contribution */
	Vchem->sch=0.;
	/* basis species */
	for(i=0;i<Vchem->size;++i) Vchem->sch += SM_basis->scharge[i]*mpow[i];
	/* surface complexes */
	for(i=0;i<SM_all->size[0];++i) Vchem->sch += SM_all->scharge[i]*smpow[i];
	Vchem->sch *= Vchem->ICS->CP->F/Vchem->ICS->CP->SA;

	g_tmp_1 = g_tmp_2 = 0.;
	for(i=0;i<Vchem->size;++i) g_tmp_1 += mpow[i]*(pow(E,-SM_basis->charge[i])-1.);
	for(i=0;i<SM_all->size[0];++i) g_tmp_2 += smpow[i]*(pow(E,-SM_all->charge[i])-1.);
	Fchem[Vchem->size_mass+1] = kappa_s*g_tmp_1;
	Fchem[Vchem->size_mass+1] += kappa_s*g_tmp_2;
	Fchem[Vchem->size_mass+1] -= Vchem->sch*Vchem->sch;

	g_tmp_1=g_tmp_2=g_tmp_3=0.;
	for(i=0;i<SM_all->size[0];++i)
	{
		g_tmp_1 += SM_all->scharge[i]*smpow[i]*SM_all->M[i][Vchem->pos_exp];
		g_tmp_2 += smpow[i]*SM_all->M[i][Vchem->pos_exp]*(pow(E,-SM_all->charge[i])-E);
		g_tmp_3 += smpow[i]*SM_all->charge[i]*pow(E,-SM_all->charge[i]);

	}
	jacobi[Vchem->size_mass+1][Vchem->size_mass+1] = -2.*Vchem->ICS->CP->F/Vchem->ICS->CP->SA*Vchem->sch*g_tmp_1 + 
		g_tmp_2*kappa_s - g_tmp_3*kappa_s;
	
	g_tmp_1 = 0.;
	for(i=0;i<Vchem->size;++i) g_tmp_1 += SM_basis->charge[i]*mpow[i]*pow(E,-SM_basis->charge[i]);

	jacobi[Vchem->size_mass+1][Vchem->size_mass+1] -= kappa_s*g_tmp_1;
	jacobi[Vchem->size_mass+1][Vchem->size_mass+1]*= LNTEN;

	/* d(Grahame)/d log_10 m_q */

	for(q=0;q<Vchem->size_mass;++q)
	{
		pos_q = Vchem->pos_mass[q];
		jacobi[Vchem->size_mass+1][q+1] = kappa_s*mpow[pos_q]*(pow(E,-SM_basis->charge[pos_q])-1.);

		g_tmp_1 = SM_basis->scharge[pos_q]*mpow[pos_q];
		g_tmp_2=0.;

		for(i=0;i<SM_all->size[0];++i)
		{
			
			g_tmp_1 += SM_all->scharge[i]*smpow[i]*SM_all->M[i][pos_q];
			g_tmp_2 += smpow[i]*SM_all->M[i][pos_q]*(pow(E,-SM_all->charge[i])-1.);
		}
		jacobi[Vchem->size_mass+1][q+1] = jacobi[Vchem->size_mass+1][q+1] -2.*Vchem->ICS->CP->F/Vchem->ICS->CP->SA*Vchem->sch*g_tmp_1 
			+ g_tmp_2*kappa_s;
		jacobi[Vchem->size_mass+1][q+1] *= LNTEN;

	}


	/* SUP SAT MINERALS */
	if(Vchem->size_sup_min>0)
	{
		calc_jacobi_sup_min(jacobi,Vchem);
	}
	
	/*-----------------------*/
	/*debug*/
	/*
	j_num = (real **) calloc(Vchem->size_mass+1, sizeof(real *));
	for(i=0;i<Vchem->size_mass+1;++i) j_num[i] = (real *) calloc(Vchem->size_mass+1, sizeof(real ));
	calc_complex(&SM_mineral,Vchem);
	jacobi_num(x, j_num, Vchem);

	fp=my_fopen("debug_jacobi.out", "w");
	fprintf(fp,"jacobi num\n");
	for(i=0;i<Vchem->size_mass;++i)
	{
		for(j=0;j<Vchem->size_mass;++j)
		{
			fprintf(fp,"%4.8e\t%",j_num[i+1][j+1]);
		}
		fprintf(fp,"\n");
	}
	fprintf(fp,"\n");
	fprintf(fp,"jacobi analytic\n");
	for(i=0;i<Vchem->size_mass;++i)
	{
		for(j=0;j<Vchem->size_mass;++j)
		{
			fprintf(fp,"%4.8e\t%",jacobi[i+1][j+1]);
		}
		fprintf(fp,"\n");
	}
	fprintf(fp,"\n");
	fclose(fp);

	for(i=0;i<Vchem->size_mass;++i)
		for(j=0;j<Vchem->size_mass;++j) jacobi[i+1][j+1]=j_num[i+1][j+1];

	for(i=0;i<Vchem->size_mass;++i) free(j_num[i]);
	free(j_num);
	*/
/*-------------------------------------------------------*/	
	
return;
}

void newton_ah_log_num(real *x, real *Fchem, struct BasVec *Vchem)
{
  
    int i,j,pos, pos_q;
    //real  g_tmp_j, g_tmp_q, g_tmp_p;
	real *delta_calc, **jac_dum, *mpow, *smpow;

	struct ChemTable *SM_mineral, *SM_all, *SM_basis;
	SM_mineral = Vchem->ICS->SM_mineral;
	SM_all     = Vchem->ICS->SM_all;
	SM_basis   = Vchem->ICS->SM_basis;

	
	delta_calc = Vchem->ICS->mem->dVmb;
	jac_dum    = Vchem->ICS->mem->ddVmb;
	mpow       = Vchem->ICS->mem->dVmc;
	smpow      = Vchem->ICS->mem->dSMa;

    for (i = 0;i < Vchem->size_mass;++i) /* construct basis vector */
	{
		pos = Vchem->pos_mass[i];
		Vchem->log_m[pos] = (real) x[i];
		Vchem->log_a[pos] = Vchem->log_m[pos]+Vchem->log_g[pos];
	}

	
	/* Surface charge contribution */
	/* calculate surface charge */
	if(Vchem->surface_flag == 1)
	{	
		Vchem->log_a[Vchem->pos_exp]=Vchem->log_m[Vchem->pos_exp]=x[Vchem->size_mass];
		Vchem->psi = Vchem->log_a[Vchem->pos_exp]*LNTEN/Vchem->ICS->CP->beta_chem;
		
	}

	
/*-----------------------------------------------------------------------------------------------*/

	if(Vchem->size_rock>0) calc_rock_spec_conc(Vchem); 
	
	/* calculate m = 10^log_m */
	for(i=0;i<Vchem->size;++i) mpow[i] = pow10(Vchem->log_m[i]);
	/* calculate concentration of complexes */
	calc_complex(SM_all, Vchem);

	/* calculate SM_all_m = 10^log_SM_all_m */
	for(i=0;i<SM_all->size[0];++i) smpow[i] = pow10(SM_all->log_m[i]);

	/* amount of rock species precipitated/dissolved                                  */
	/* delta_calc > 0 => species conc decreases in bulk and increases at the surface  */
	for(j=0; j < Vchem->size_rock ; ++j)
	{
		pos  = Vchem->pos_rock[j];
		delta_calc[j]= - mpow[pos];
		for(i=0; i< SM_all->size[0]; ++i)
		{
			delta_calc[j] -= SM_all->M[i][pos]*smpow[i];
		}
		delta_calc[j] += Vchem->ctot[pos];
	}

	for(j=0;j<Vchem->size_rock;++j)
	{
		pos = Vchem->pos_buffer[j];
		SM_mineral->delta[pos]=0.;
		for(i=0;i<Vchem->size_rock;++i)
		{
			SM_mineral->delta[pos]+= Vchem->sm_buf[j][i]*delta_calc[i];
		}
	}

	/* calculate F_j*/ 
	/* F_j = c_tot_j - m_j - sum_i^N_rb M_ij delta_i-sum_i^Nx SM_i,j x_i */
	Fchem[0]=0.;
	for(j=0; j < Vchem->size_mass ; ++j)
	{
		pos  = Vchem->pos_mass[j];
		Fchem[j+1] = Vchem->ctot[pos]-mpow[pos];
		
		if(Vchem->DL) 
		for(i=0; i< SM_all->size[0]; ++i)
		{
			Fchem[j+1] -= SM_all->M[i][pos]*smpow[i];
		}
			/* add amount lost or gained due to precipitation/dissolution */
		for(i=0;i<Vchem->size_rock;++i)
		{
			pos_q=Vchem->pos_buffer[i];
			Fchem[j+1] -= SM_mineral->M[pos_q][pos]*SM_mineral->delta[pos_q];
		}

	}


	/* SURFACE CHARGE */
	if(Vchem->surface_flag == 1)
	{
		Vchem->sch=0.;
		/* basis species */
		for(i=0;i<Vchem->size;++i) Vchem->sch += SM_basis->scharge[i]*pow10(Vchem->log_m[i]);
		/* surface complexes */
		for(i=0;i<SM_all->size[0];++i) Vchem->sch += SM_all->scharge[i]*pow10(SM_all->log_m[i]);

		Vchem->sch *= Vchem->ICS->CP->F/Vchem->ICS->CP->SA;
		Fchem[Vchem->size_mass+1]=0.;
		/* sum up basis species */
		for(i=0;i<Vchem->size;++i) Fchem[Vchem->size_mass+1] += pow10(Vchem->log_m[i])*(exp(-SM_basis->charge[i]*Vchem->ICS->CP->beta_chem*Vchem->psi)-1.);
	
		/* sum up complexes */
		for(i=0;i<SM_all->size[0];++i) Fchem[Vchem->size_mass+1] += pow10(SM_all->log_m[i])*(exp(-SM_all->charge[i]*Vchem->ICS->CP->beta_chem*Vchem->psi)-1.);
		/* factor 2*e*e_0*k_B*1000*N_a = 1.1776E-05*/
		Fchem[Vchem->size_mass+1] *= 1.778e-5*Vchem->Temp;
		Fchem[Vchem->size_mass+1] -= Vchem->sch*Vchem->sch;
		/*------------------------*/
	}

	
return;
}



void mass_balance(int no_call, struct BasVec *Vchem)
{
  int i, iter, pos;
  real criterion, delta;
  real *x_init;
  real *chem_F, **chem_Jacobi;
  int *chem_indx;

  //static int count = 0;
  //count++;

  for(i=0;i<Vchem->size;++i) {
    if (Vchem->ctot_calc[i] < 0) {
      printf("WARNING! Neg. ctot for %s: %g\n", Vchem->ICS->SM_basis->row_name[i], Vchem->ctot_calc[i]);
      exit(1);
      MPI_Finalize();
    }
  }
  

  chem_F    = Vchem->ICS->mem->dVma;
  chem_indx = Vchem->ICS->mem->iVma;
  chem_Jacobi = Vchem->ICS->mem->ddVma;
  
  
  x_init = (real *) calloc(Vchem->size_mass, sizeof(real));
  for(i=0; i < Vchem->size_mass;++i) {
    pos=Vchem->pos_mass[i];
    x_init[i] = Vchem->log_m[pos];
  }
  
  iter=0;
  do {
    // REMOVE
    /* if (_DEBUG_FLAG_) */
    /*   printf("iter: %d\n",iter); */
    // REMOVE
    criterion = 0.;
    /* calculate jacobi_ and Fchem_ */
    /* equilibrate with ion exchanger! */
    if (Vchem->pos_X > -1 && Vchem->equilibrate) {
      newton_ah_log_io(x_init, chem_F, chem_Jacobi, Vchem);
      //if (count==10000) printf("log_io\n");
    } else { 
      //if (_DEBUG_FLAG_) printf("      newton_ah_log(x_init, chem_F, chem_Jacobi, Vchem)\n");
      newton_ah_log(x_init, chem_F, chem_Jacobi, Vchem);
      //if (count==10000) printf("ah_log\n");
    }
    
    for (i=0;i<Vchem->size_mass;++i) 
      criterion += fabs(chem_F[i+1]);
    if (criterion <1e-8) 
      break; 
    
    if (Vchem->ICS->PRINT_DEBUG_CHEM) 
      write_jacobi(chem_Jacobi, Vchem->size_mass,"debug_j.out");


    /* solve system */
    ludcmp(chem_Jacobi, Vchem->size_mass , chem_indx, chem_F);
    lubksb(chem_Jacobi, Vchem->size_mass , chem_indx, chem_F);
    
    for (i=0;i<Vchem->size_mass;++i) {
      if (fabs(chem_F[i+1]) > MBAL_MAX_STEP_) {
	if (chem_F[i+1] < 0) 
	  delta = -MBAL_MAX_STEP_;
	else 
	  delta = MBAL_MAX_STEP_;
      } else 
	delta = chem_F[i+1];
      x_init[i] -= delta;
      
      if(x_init[i]>10) { /* 10 mol/l ..... */
	pos=Vchem->pos_mass[i];
	x_init[i] = log10(Vchem->ctot[pos]);
      }
      
      
    }
    
    ++iter;
    if (iter > CHEM_MAXITER_MASS_) 
      break;
    
  } while ( criterion > 1e-8);
  
  //if (count > 20000) {
  //  writeVchem(Vchem, "vchemOK.out");
  //  exit(1);
  //}
  if (iter > CHEM_MAXITER_MASS_) {
    printf("1 WARNING MASSBALANCE DID NOT CONVERGE FINAL VALUE %g, %d\n", criterion, no_call);
    writeVchem(Vchem, "vchem.out");
    //exit(1);
  }
  free(x_init); 
  x_init=NULL;
  if (Vchem->ICS->PRINT_DEBUG_CHEM)
    printf("No iter %d and final value %12.8e\n", iter, criterion);
  
  
}
void mass_balance_surf(int no_call, struct BasVec *Vchem)
{
	int i, iter, pos, dim;
	real criterion, delta;
	real *x_init, E_init;
	real *chem_F, **chem_Jacobi;
	int *chem_indx;

	chem_F    = Vchem->ICS->mem->dVma;
	chem_indx = Vchem->ICS->mem->iVma;
	chem_Jacobi = Vchem->ICS->mem->ddVma;

	dim = Vchem->size_mass+1;

	x_init = (real *) calloc(dim, sizeof(real));
	
	for(i=0; i < dim-1;++i)
	{
		pos=Vchem->pos_mass[i];
		x_init[i] = Vchem->log_m[pos];
	}
	
	E_init = x_init[Vchem->size_mass] = Vchem->log_a[Vchem->pos_exp];	
	
	iter=0;
	do 
	{
		/* 
		Make sure that surface charge and surface potential has the same sign,
		solution of the euquations are not uniqe
		*/
		if(Vchem->psi*Vchem->sch < 0)
		{
		/*	printf("MASSBALANCE: CHANGE SIGN OF SURFACE POTENTIAL\n");*/
			x_init[Vchem->size_mass] = -x_init[Vchem->size_mass];
		}
	
		/* calculate jacobi_ and Fchem_ */
		newton_ah_log_surf(x_init, chem_F, chem_Jacobi, Vchem);
		
		criterion = 0.;		
		for(i=0;i<dim;++i) criterion += fabs(chem_F[i+1]);

		/* solve system */
		ludcmp(chem_Jacobi, dim , chem_indx, chem_F);
		lubksb(chem_Jacobi, dim , chem_indx, chem_F);
		for(i=0;i<Vchem->size_mass;++i)
		{
			if(fabs(chem_F[i+1]) > MBAL_MAX_STEP_)
			{
				if(chem_F[i+1] < 0) delta = -MBAL_MAX_STEP_;
				else delta = MBAL_MAX_STEP_;
			} else delta = chem_F[i+1];
			x_init[i] -= delta;
	
			if(x_init[i]>10) /* 10 mol/l ..... */
			{
				pos=Vchem->pos_mass[i];
				x_init[i] = log10(Vchem->ctot[pos]);
			}
			
			
		}

		if(fabs(chem_F[Vchem->size_mass+1]) > MBAL_SURF_MAX_STEP_)
		{
			if(chem_F[Vchem->size_mass+1] < 0) delta = -MBAL_SURF_MAX_STEP_;
			else delta = MBAL_SURF_MAX_STEP_;
		} 
		else delta = chem_F[Vchem->size_mass+1];

		x_init[Vchem->size_mass] -= delta;
		if(fabs(x_init[Vchem->size_mass]*LNTEN/Vchem->ICS->CP->beta_chem) >1) /* potental of 1 V*/
		{
			x_init[Vchem->size_mass] = E_init;
		}


		
		++iter;
		if(iter > CHEM_MAXITER_MASS_) break;

	} while( criterion > 1e-8);

	if(iter > CHEM_MAXITER_MASS_)
	{
		printf("2 WARNING MASSBALANCE DID NOT CONVERGE FINAL VALUE %g\n", criterion);
	}
	free(x_init); x_init=NULL;
	if(Vchem->ICS->PRINT_DEBUG_CHEM)printf("No iter %d and final value %12.8e\n", iter, criterion);
	if(Vchem->ICS->PRINT_DEBUG_CHEM)printf("Surface charge: %g\t Surface pot: %g\n",Vchem->sch, Vchem->psi);



}

void mass_balance_num(int no_call, struct BasVec *Vchem)
{
	int i, j,iter, pos, dim;
	real criterion, delta;
	real *x_init;

	real *dx;
	FILE *fp;
	
	real *chem_F, **chem_Jacobi;
	int *chem_indx;

	chem_F    = Vchem->ICS->mem->dVma;
	chem_indx = Vchem->ICS->mem->iVma;
	chem_Jacobi = Vchem->ICS->mem->ddVma;

	if(Vchem->surface_flag == 1) dim = Vchem->size_mass+1;
	else dim = Vchem->size_mass;
	x_init = (real *) calloc(dim, sizeof(real));
	
	for(i=0; i < dim;++i)
	{
		pos=Vchem->pos_mass[i];
		x_init[i] = Vchem->log_m[pos];
	}
	dx = (real *) calloc(dim,sizeof(real));
	for(i=0;i<dim;++i) dx[i] = 1.e-3;
	
		/* DEBUG */
	if(Vchem->surface_flag == 1)
	{
		x_init[Vchem->size_mass] = Vchem->log_a[Vchem->pos_exp];
		dx[Vchem->size_mass] = 1e-3;
	}
		
	
	

	iter=0;

	do 
	{
		/* DEBUG */
		if(Vchem->surface_flag == 1)
		{
			if(Vchem->psi*Vchem->sch < 0)
			{
				printf("MASSBALANCE: CHANGE SIGN OF SURFACE POTENTIAL\n");
				x_init[Vchem->size_mass] = -x_init[Vchem->size_mass];
			}
		}
		
		fp = fopen("jacobi.out","w");
		for(i=0;i<dim;++i)
			for(j=0; j<dim;++j) chem_Jacobi[i+1][j+1]=0.;

		newton_ah_log_surf(x_init, chem_F, chem_Jacobi, Vchem);
 		fprintf(fp,"analytical:\n");
		for(i=0;i<dim;++i)
		{
			for(j=0;j<dim;++j)
			{
				fprintf(fp,"%g\t",chem_Jacobi[i+1][j+1]);
			}
			fprintf(fp,"\n");
		}
		for(i=0;i<dim;++i)
			for(j=0; j<dim;++j) chem_Jacobi[i+1][j+1]=0.;
		jacobi_num2(x_init, dx, dim, Vchem, chem_F,chem_Jacobi);
		
		fprintf(fp,"numerical:\n");
		for(i=0;i<dim;++i)
		{
			for(j=0;j<dim;++j)
			{
				fprintf(fp,"%g\t",chem_Jacobi[i+1][j+1]);
			}
			fprintf(fp,"\n");
		}
		fclose(fp);
		
		criterion = 0.;
		/* calculate jacobi_ and Fchem_ */
		
		for(i=0;i<dim;++i) criterion += fabs(chem_F[i+1]);

	/* solve system */
		ludcmp(chem_Jacobi, dim , chem_indx, chem_F);
		lubksb(chem_Jacobi, dim , chem_indx, chem_F);
		for(i=0;i<dim;++i)
		{
			if(fabs(chem_F[i+1]) > MBAL_MAX_STEP_)
			{
				if(chem_F[i+1] < 0) delta = -MBAL_MAX_STEP_;
				else delta = MBAL_MAX_STEP_;
			} else delta = chem_F[i+1];
			x_init[i] -= delta;
	
			if(x_init[i]>10) /* 10 mol/l ..... */
			{
				pos=Vchem->pos_mass[i];
				x_init[i] = log10(Vchem->ctot[pos]);
			}
			
			
		}

		
		++iter;
		/*********TEST************/
		/*-----------------------------STEP 2---------------------------------*/
	/*--------------------- calculate Io----------------------------------*/
/*	Vchem->Io=calc_Io(Vchem);*/

    /*-----------------------------STEP 3---------------------------------*/
	/*--------------------- calculate gamma-------------------------------*/

		if(iter > CHEM_MAXITER_MASS_) break;

	} while( criterion > 1e-8);

	if(iter > CHEM_MAXITER_MASS_) printf("3 WARNING MASSBALANCE DID NOT CONVERGE FINAL VALUE %g\n", criterion);
	free(x_init); x_init=NULL;
	if(Vchem->ICS->PRINT_DEBUG_CHEM)printf("No iter %d and final value %12.8e\n", iter, criterion);
	if(Vchem->ICS->PRINT_DEBUG_CHEM)printf("Surface charge: %g\t Surface pot: %g\n",Vchem->sch, Vchem->psi);



}

/*
jacobi_num:
  calculates the jacobian of the species determined by mass balance numerically
  INPUT: 1) *x: the point which the Jacobian is to be calculated
         2) Vchem:  vector struct 
  OUTPUT: 1) j_num: the numerical jacobian indexed from 1,..,dim_mass
  RETURN: void
*/
void jacobi_num(real *x, real **j_num, struct BasVec *Vchem)
{
	real dx=1e-3, *F_init;
	int i, j, pos,q, pos_j;

	/* The function value at x */
	/* F = c_tot - c_tot_calc  */
	F_init = (real *) calloc(Vchem->size_mass,sizeof(real));
	for (i = 0;i < Vchem->size_mass;++i) 
	{
		pos = Vchem->pos_mass[i];
		Vchem->log_m[pos] = (real) x[i];
		Vchem->log_a[pos] = Vchem->log_m[pos]+Vchem->log_g[pos];
	}
	calc_ctot(Vchem);
	for(i=0;i<Vchem->size_mass;++i)
	{
		pos = Vchem->pos_mass[i];
		F_init[i]=Vchem->ctot[pos]-Vchem->ctot_calc[pos];
	}
	

	for (q = 0; q < Vchem->size_mass;++q) /* construct basis vector */
	{ 
		pos = Vchem->pos_mass[q];
		Vchem->log_m[pos] = (real) x[q]-dx;
		Vchem->log_a[pos] = Vchem->log_m[pos]+Vchem->log_g[pos];
		/* calculate c_tot */
		calc_ctot(Vchem);
		for(j=0;j<Vchem->size_mass;++j)
		{
			pos_j=Vchem->pos_mass[j];
			j_num[j+1][q+1] = -(Vchem->ctot[pos_j]-(Vchem->ctot_calc[pos_j]+F_init[j]))/dx;
		}
		/* set back */
		Vchem->log_m[pos] = (real) x[q];
		Vchem->log_a[pos] = Vchem->log_m[pos]+Vchem->log_g[pos];
		if(Vchem->size_rock>0) calc_rock_spec_conc(Vchem);
	}
	free(F_init);

	return;
}
/*
calc_ctot:
  calculates the total basis species concentration aqueous + surface
  INPUT: Vchem
  OUTPUT: Vchem.ctot_calc
  RETURN: void
*/
void calc_ctot(struct BasVec *Vchem)
{
	int pos, i,j;
	real *delta_calc;
	struct ChemTable *SM_mineral, *SM_all;
	SM_mineral = Vchem->ICS->SM_mineral;
	SM_all     = Vchem->ICS->SM_all;

	delta_calc = (real *) calloc(Vchem->size_rock, sizeof(real));
	
	if(Vchem->size_rock>0) calc_rock_spec_conc(Vchem);
	calc_complex(SM_all, Vchem); 

	/* add amount lost or gained */
	for(j=0; j < Vchem->size_rock ; ++j)
	{
		pos  = Vchem->pos_rock[j];
		delta_calc[j]=-pow10(Vchem->log_m[pos]);
		for(i=0; i< SM_all->size[0]; ++i)
		{
			delta_calc[j] -= SM_all->M[i][pos]*pow10(SM_all->log_m[i]);
		}
		delta_calc[j] += Vchem->ctot[pos];
	}

	for(j=0;j<Vchem->size_rock;++j)
	{
		pos = Vchem->pos_buffer[j];
		SM_mineral->delta[pos]=0.;
		for(i=0;i<Vchem->size_rock;++i)
		{
			SM_mineral->delta[pos]+= Vchem->sm_buf[j][i]*delta_calc[i];
		}
	}

	/* c_tot_j = m_j + sum_i^Nx SM_i,j x_i */
	for(j=0; j < Vchem->size; ++j)
	{
		Vchem->ctot_calc[j] = pow10(Vchem->log_m[j]);
		for(i=0; i< SM_all->size[0]; ++i)
		{
			Vchem->ctot_calc[j] += SM_all->M[i][j]*pow10(SM_all->log_m[i]);
		}
		for(i=0;i<Vchem->size_rock;++i)
		{
			pos=Vchem->pos_buffer[i];
			Vchem->ctot_calc[j] += SM_mineral->M[pos][j]*SM_mineral->delta[pos];
		}
	}
	if(Vchem->size_sup_min>0)
	{
		calc_ctot_sup_min(Vchem->ctot_calc,Vchem);
	}

	free(delta_calc);
	
	return;
}

/*
calc_ctot_aq:
  calculates the total basis species concentration aqueous
  INPUT: Vchem
  OUTPUT: Vchem.ctot_calc
  RETURN: void
*/
void calc_ctot_aq(struct BasVec *Vchem)
{
	int pos, i,j,posb;
	struct ChemTable *SM_all, *SM_mineral;
	
	SM_all     = Vchem->ICS->SM_all;
	SM_mineral = Vchem->ICS->SM_mineral;

	if(Vchem->size_rock>0) calc_rock_spec_conc(Vchem);
	calc_complex(SM_all, Vchem); 

	/* c_tot_j = m_j + sum_i^Nx SM_i,j x_i */
	for(j=0; j < Vchem->size; ++j)
	{
		Vchem->ctot_calc[j] = 0.;
		if(SM_all->type[j] == 0) Vchem->ctot_calc[j] = pow10(Vchem->log_m[j]);
		if(j==Vchem->pos_pH) Vchem->ctot_calc[j] += Vchem->WHTO;

		for(i=0; i< SM_all->size[0]; ++i)
		{
			/*only aqueous species, not surface species - treated seperately*/
			if(SM_all->type[i] == 0) Vchem->ctot_calc[j] += SM_all->M[i][j]*pow10(SM_all->log_m[i]);
		}
	}
	
	/* store amount precip */
	for(i=0;i<Vchem->size_rock;++i)
	{
		pos = Vchem->pos_rock[i];
		posb = Vchem->pos_buffer[i];
		Vchem->delta_mineral[pos] = SM_mineral->delta[posb];
	}

	return;
}

void jacobi_num2(real *x, real *dx, int dim, struct BasVec *Vchem, real *chem_F, real **jc)
{
	real *x_init, *xp, *xm;

	int i, j,k;
	x_init = (real *) calloc(dim, sizeof(real));
	xp = (real *) calloc(dim, sizeof(real));
	xm = (real *) calloc(dim, sizeof(real));


	for(i=0;i<dim;++i) x_init[i] = x[i];

	for(i=0;i<dim;++i)
	{
		for(k=0;k<dim;++k) xp[k] = xm[k]=x_init[k] = x[k];
		xp[i] = x_init[i]+ dx[i];
		xm[i] = x_init[i]- dx[i];

		newton_ah_log_num(xp, chem_F, Vchem);
		for(j=0;j<dim;++j) jc[j+1][i+1] = .5*chem_F[j+1]/dx[i];

		newton_ah_log_num(xm, chem_F, Vchem);
		for(j=0;j<dim;++j) jc[j+1][i+1] -= .5*chem_F[j+1]/dx[i];

	}

	for(k=0;k<dim;++k) x_init[k] = x[k];
	newton_ah_log_num(x_init, chem_F, Vchem);

	free(xp);
	free(xm);
	free(x_init);

}


/*NB: IF RATE EQUATIONS ARE CHANGED - REMEMEBER TO CHANGE IN "calc_ctot_sup_min" ASWELL */
/* c_tot_calc = m_i + sum_j mu_ij n_i + e.q. minerals - (k_1+k_2*a_H)*(1-SI^k_3)^k_4 */
/* F_i = c_tot-c_tot_calc */
void calc_F_sup_min(real *Fchem, struct BasVec *Vchem)
{
  real SI, sgn, fact,mol;
  int i, j, k, posi;
  /* 1J = 0.238902957619 cal conversion between kcal J */
  /*real Ea = 15e+3/0.238902957619; 
    real T0 = 273.15+130;

    fact = exp(-Ea/Vchem->ICS->CP->Rg*(1./Vchem->Temp-1./T0));*/
  fact = 1.;
  fact *= Vchem->RP->dt;

  for(j=0;j<Vchem->size_mass;++j) {
    for(i=0;i<Vchem->size_sup_min;++i) {
      posi = Vchem->pos_sup_min[i];
      mol = Vchem->RP->mol[i];
			
      SI =0.;
      for(k=0;k<Vchem->ICS->SM_mineral->size[1];++k) SI += Vchem->ICS->SM_mineral->M[posi][k]*Vchem->log_a[k];
      SI -= Vchem->ICS->SM_mineral->logK[posi];
      SI -= Vchem->ICS->SM_mineral->log_af[posi];
      Vchem->ICS->SM_mineral->log_a[posi] = SI;
      Vchem->ICS->SM_mineral->log_m[posi] = SI-Vchem->ICS->SM_mineral->log_g[posi];
      sgn = ( SI < 0 ? 1. : -1.);

      //printf("Vchem->rate[%d][1] = %.3e\n",i, Vchem->rate[i][1]);
      
      Fchem[j+1] += sgn*mol*fact*Vchem->ICS->SM_mineral->M[posi][Vchem->pos_mass[j]]
	*(Vchem->rate[i][0]+Vchem->rate[i][1]*pow10(Vchem->log_a[Vchem->pos_pH]))
	*pow(fabs(1.-pow10(SI*Vchem->rate[i][2])),Vchem->rate[i][3]); 
    }
  }
}

void calc_ctot_sup_min(real *Fchem, struct BasVec *Vchem)
{
  real SI, sgn, fact, mol;
  int i, j, k, posi;
  /* 1J = 0.238902957619 cal conversion between kcal J */
  /*real Ea = 15e+3/0.238902957619; 
    real T0 = 273.15+130;

    fact = exp(-Ea/Vchem->ICS->CP->Rg*(1./Vchem->Temp-1./T0));*/
  fact = 1.;
  fact *= Vchem->RP->dt;

  
  for(j=0;j<Vchem->size;++j) {
    for(i=0;i<Vchem->size_sup_min;++i) {
      posi = Vchem->pos_sup_min[i];
      mol = Vchem->RP->mol[i];
      
      SI =0.;
      for(k=0;k<Vchem->ICS->SM_mineral->size[1];++k)
	SI += Vchem->ICS->SM_mineral->M[posi][k]*Vchem->log_a[k];
      SI -= Vchem->ICS->SM_mineral->logK[posi];
      SI -= Vchem->ICS->SM_mineral->log_af[posi];
      sgn = ( SI < 0 ? 1. : -1.);
      
      Vchem->ICS->SM_mineral->log_a[posi] = SI;
      Vchem->ICS->SM_mineral->log_m[posi] = SI-Vchem->ICS->SM_mineral->log_g[posi];

      Fchem[j] -= sgn*mol*fact*Vchem->ICS->SM_mineral->M[posi][j]
	*(Vchem->rate[i][0]+Vchem->rate[i][1]*pow10(Vchem->log_a[Vchem->pos_pH]))
	*pow(fabs(1.-pow10(SI*Vchem->rate[i][2])),Vchem->rate[i][3]); 
    }
  }
}


/* j_i,k += sum_a (k_1+k_2*a_h)*(1-SI^k_3)^(k_4-1)*k_3*k_4*SI^k_3 */
void calc_jacobi_sup_min(real **jacobi, struct BasVec *Vchem)
{
	int j, q, i, pos_q, pos_i, pos_ib, pos;
	real g_tmp_q, SI,SIp, k1, k2, k3, k4, aH, sgn, fact,mol;

	aH = pow10(Vchem->log_a[Vchem->pos_pH]);
	fact = Vchem->RP->dt;

	for(j=0; j <Vchem->size_mass ; ++j)
	{
		pos  = Vchem->pos_mass[j];

		for(q=0; q<Vchem->size_mass; ++q)
		{
			g_tmp_q = 0.;
			pos_q  = Vchem->pos_mass[q];
			
			for(i=0; i < Vchem->size_sup_min; ++i)
			{
				pos_i  = Vchem->pos_sup_min[i];
				pos_ib = Vchem->pos_sup_bas[i];
			
				mol = Vchem->RP->mol[i];

				k1  = Vchem->rate[i][0];
				k2  = Vchem->rate[i][1];
				k3  = Vchem->rate[i][2];
				k4  = Vchem->rate[i][3];
				SI  = pow10(Vchem->ICS->SM_mineral->log_a[pos_i]);
				SIp = pow(SI,k3);
				sgn = ( SI < 1 ? 1. : -1.);
				
				if(SIp < 1)
					g_tmp_q += sgn*mol*fact*(k1+k2*aH)*pow(fabs(1.-SIp),k4-1.)*k3*k4*SIp*Vchem->ICS->SM_mineral->M[pos_i][pos]
				*Vchem->ICS->SM_mineral->M[pos_i][pos_q];
				else 
					g_tmp_q -= sgn*mol*fact*(k1+k2*aH)*pow(fabs(1.-SIp),k4-1.)*k3*k4*SIp*Vchem->ICS->SM_mineral->M[pos_i][pos]
				*Vchem->ICS->SM_mineral->M[pos_i][pos_q];
			}

	
			jacobi[j+1][q+1] -= g_tmp_q*LNTEN;
		}
	}
}
		

		

		



