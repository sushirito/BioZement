#include "chem_global.h"
void write_chem_struct(struct ChemTable *tab, const char *name)
{
  FILE *fpd;
  int i, j;

	fpd=my_fopen(name , "w");
	
	
	/*
	if(tab->abs_pos_bas != NULL)for(i=0;i<tab->size[1];++i)fprintf(fpd,"%d\t", tab->abs_pos_bas[i]);
	fprintf(fpd,"\n");
	*/
	fprintf(fpd,"row_name\t");
	if(tab->logK != NULL) fprintf(fpd,"logK\t");
	if(tab->mol_volume != NULL) fprintf(fpd,"mol_volume\t");
	if(tab->mol_volume != NULL) fprintf(fpd,"mol_weight\t");
	if(tab->log_a != NULL) fprintf(fpd,"log_a\t");
	if(tab->log_m != NULL) fprintf(fpd,"log_m\t");
	for(i=0;i<tab->size[1];++i)
	{
		fprintf(fpd,"%s\t", tab->col_name[i]);
	}
	fprintf(fpd,"\n");
	for(i=0;i<tab->size[0];++i)
	{
		fprintf(fpd,"%s\t", tab->row_name[i]);
		if(tab->logK != NULL) fprintf(fpd,"%lf\t", tab->logK[i]);
		if(tab->mol_volume != NULL) fprintf(fpd,"%lf\t", tab->mol_volume[i]);
		if(tab->mol_weight != NULL) fprintf(fpd,"%lf\t", tab->mol_weight[i]);
		if(tab->log_a != NULL) fprintf(fpd,"%lf\t", tab->log_a[i]);
		if(tab->log_m != NULL) fprintf(fpd,"%lf\t", tab->log_m[i]);
		for(j=0;j<tab->size[1];++j) fprintf(fpd,"%lf\t", tab->M[i][j]);
		fprintf(fpd,"\n");
	}
	
	if(tab->charge != NULL)
	{
		fprintf(fpd,"charge\n");
		for(i=0;i<tab->size[0];++i) fprintf(fpd,"%lf\t", tab->charge[i]);
		fprintf(fpd,"\n");
	}
	if(tab->scharge != NULL)
	{
		fprintf(fpd,"scharge\n");
		for(i=0;i<tab->size[0];++i) fprintf(fpd,"%lf\t", tab->scharge[i]);
		fprintf(fpd,"\n");
	}
	if(tab->a0 != NULL)
	{
		fprintf(fpd,"a0\n");
		for(i=0;i<tab->size[0];++i) fprintf(fpd,"%lf\t", tab->a0[i]);
		fprintf(fpd,"\n");
	}


	fclose(fpd);

}



void write_BasVec_struct(struct BasVec *Vchem, const char *name)
{
  FILE *fpd;
  int i, j;
  struct ChemTable *SM_mineral, *SM_all, *SM_basis;
	SM_mineral = Vchem->ICS->SM_mineral;
	SM_all     = Vchem->ICS->SM_all;
	SM_basis   = Vchem->ICS->SM_basis;

	fpd=my_fopen(name , "w");
	
	fprintf(fpd,"name\t c_tot\t log_m\t log_a\t log_g\t c_ads\n");
	for(i=0;i<Vchem->size;++i)
		fprintf(fpd,"%s\t%lf\t%lf\t%lf\t%lf\t%lf\n", SM_basis->row_name[i], Vchem->ctot[i],Vchem->log_m[i],Vchem->log_a[i],Vchem->log_g[i], Vchem->ctot_ads[i]);
	fprintf(fpd,"Mass balance:\n");
	for(i=0;i<Vchem->size_mass;++i)
		fprintf(fpd,"%s ",SM_basis->row_name[Vchem->pos_mass[i]]);
	fprintf(fpd,"\n");
	fprintf(fpd,"Rock basis:\n");
	for(i=0;i<Vchem->size_rock;++i)
		fprintf(fpd,"%s ",SM_basis->row_name[Vchem->pos_rock[i]]);
	fprintf(fpd,"\n");
	fprintf(fpd,"Sup basis:\n");
	for(i=0;i<Vchem->size_sup_min;++i)
		fprintf(fpd,"%s ",SM_basis->row_name[Vchem->pos_sup_bas[i]]);
	fprintf(fpd,"\n");
	fprintf(fpd,"Rock buffer:\n");
	for(i=0;i<Vchem->size_rock;++i)
		fprintf(fpd,"%s ",SM_mineral->row_name[Vchem->pos_buffer[i]]);
	fprintf(fpd,"\n");
	fprintf(fpd,"Supersaturated minerals:\n");
	for(i=0;i<Vchem->size_sup_min;++i)
		fprintf(fpd,"%s ",SM_mineral->row_name[Vchem->pos_sup_min[i]]);
	fprintf(fpd,"\n");

	fprintf(fpd,"Basis transformation matrix:\n");
	for(i=0;i<Vchem->size;++i)
	{
		for(j=0;j<Vchem->size;++j)
			fprintf(fpd,"%lf\t",Vchem->beta_bas_inv[i][j]);

		fprintf(fpd,"\n ");
	}

	fclose(fpd);

}

void writeVchem(struct BasVec *Vp, const char *name)
{
	int i;
	real ch_dum;
	FILE *fp;
	struct ChemTable *SM_all, *SM_mineral, *SM_basis;

	fp = my_fopen(name,"w");
	
	SM_all     = Vp->ICS->SM_all;
	SM_mineral = Vp->ICS->SM_mineral;
	SM_basis   = Vp->ICS->SM_basis;

	fprintf(fp,"pH\t Sch\t psi[mV]\t F\t S\n");
	fprintf(fp,"%g\t %g\t %g\t %g\t %g\n",-Vp->log_a[Vp->pos_pH],Vp->sch, Vp->psi*1000, Vp->ICS->CP->F, Vp->ICS->CP->SA);
	calc_complex(SM_mineral,Vp);
	fprintf(fp,"Temp \t Io_\n");
	fprintf(fp,"%lf\t%lf\n",Vp->Temp-273.15, Vp->Io);

	fprintf(fp,"\t");
	for(i=0; i< Vp->size;++i)
		fprintf(fp,"%s\t", SM_basis->row_name[i]);
	for(i=0; i< SM_all->size[0];++i)
		fprintf(fp,"%s\t",SM_all->row_name[i]);
	
	fprintf(fp,"\n");
	fprintf(fp,"type\t");
	for(i=0; i< Vp->size;++i)
		fprintf(fp,"%d\t",(int) SM_basis->type[i]);
	for(i=0; i< SM_all->size[0];++i)
		fprintf(fp,"%d \t",(int) SM_all->type[i]);
	fprintf(fp,"\n");
	fprintf(fp,"log_a\t");
	for(i=0; i< Vp->size;++i)
		fprintf(fp,"%lf\t",Vp->log_a[i]);
	for(i=0; i< SM_all->size[0];++i)
		fprintf(fp,"%12.8e \t",SM_all->log_a[i]);
	fprintf(fp,"\n");
	fprintf(fp,"log_m\t");
	for(i=0; i< Vp->size;++i)
		fprintf(fp,"%lf\t",Vp->log_m[i]);
	for(i=0; i< SM_all->size[0];++i)
		fprintf(fp,"%lf \t",SM_all->log_m[i]);
	fprintf(fp,"\n");
	fprintf(fp,"log_g\t");
	for(i=0; i< Vp->size;++i)
		fprintf(fp,"%lf\t",Vp->log_g[i]);
	for(i=0; i< SM_all->size[0];++i)
		fprintf(fp,"%lf \t",SM_all->log_g[i]);
	fprintf(fp,"\n");
	fprintf(fp,"charge\t");
	for(i=0; i< Vp->size;++i)
		fprintf(fp,"%lf\t",SM_basis->charge[i]);
	for(i=0; i< SM_all->size[0];++i)
		fprintf(fp,"%lf \t",SM_all->charge[i]);
	fprintf(fp,"\n");

	fprintf(fp,"Species\t ch\t c_tot\t c_dl:\n");
	write_BasVec_struct(Vp,"vchem.out");
	calc_ctot_aq(Vp);
	for(i=0;i<Vp->size;++i)
	{
		fprintf(fp,"%s\t %g\t %12.8e\t %12.8e\n", SM_basis->row_name[i], SM_basis->charge[i],Vp->ctot_calc[i], Vp->ctot_dl[i] );
	}

	fprintf(fp,"Saturation index:\n");
	for(i=0; i< SM_mineral->size[0];++i)
		fprintf(fp,"%s\t%lf\n",SM_mineral->row_name[i], SM_mineral->log_a[i]);

	fprintf(fp,"Charge Balance:\t");
	ch_dum = 0.;
	for(i=0;i<Vp->size;++i)  if(SM_basis->type[i] == 0) ch_dum += SM_basis->charge[i]*pow10(Vp->log_m[i]);
	for(i=0;i<SM_all->size[0];++i)  if(SM_all->type[i] == 0)ch_dum += SM_all->charge[i]*pow10(SM_all->log_m[i]);
	fprintf(fp,"%12.8e\n",ch_dum);

	if(Vp->pos_X >= 0)
	{
		fprintf(fp,"\n");
		fprintf(fp,"\t concen-\t Equiv-\t Equivalent\t Log\n");
		fprintf(fp,"species\t tration\t alents\t fraction\t gamma\n");
		for(i=0;i<SM_all->size[0];++i)
		{
			if(SM_all->type[i] == 2)
			{
				fprintf(fp,"%s\t %10.6e\t%10.6e\t%10.6e\t%g",SM_all->row_name[i],
					pow10(SM_all->log_m[i]),pow10(SM_all->log_m[i])*SM_all->charge[i],pow10(SM_all->log_a[i]-SM_all->log_g[i]), SM_all->log_g[i]);
				fprintf(fp,"\n");
			}
		}
	}

	fprintf(fp, "\n");
	fclose(fp);
}

void write_jacobi(real **jc, int size, const char *fname)
{
	FILE *fp;
	int i, j;

	fp=my_fopen(fname, "w");
	for(i=0;i<size;++i)
	{
		for(j=0;j<size;++j)
		{
			fprintf(fp,"%4.8e\t",jc[i+1][j+1]);
		}
		fprintf(fp,"\n");
	}
	fclose(fp);
}
	
