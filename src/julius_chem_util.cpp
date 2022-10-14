#include "chem_global.h"
//#include "global.h"


/* Mol Volume = Molecular Weight/Density */

//real calc_rock_density(real *dx, struct InitChem *ICS)
//{
//	int i;
//	real r_denst;
//	r_denst = 0.;
//	for(i=0;i<ICS->size_rock;++i)
//		r_denst += (ICS->c_buffer[ICS->pos_rock[i]]+dx[i])*ICS->SM_mineral->mol_volume[ICS->pos_buffer[i]]/ICS->SM_mineral->mol_weight[ICS->pos_buffer[i]];
//	for(i=0;i<ICS->size_sup_min;++i)
//		r_denst += (ICS->c_sup_min[i]+dx[i])*ICS->SM_mineral->mol_volume[ICS->pos_sup_min[i]]/ICS->SM_mineral->mol_weight[ICS->pos_sup_min[i]];
//
//	if(r_denst < 1e-16) return 0; else
//	return 1.e+3/r_denst; /* note mole_volume is given in L/mol & mole_weight in gram/L density in g/L*/
//}

/* calculates permeability in mD from a Carman-Kozeny relation, 
input: tortuosity, porosity, wt-fraction mineral change*/
//real calc_permeability(real *dx, real tau, real por, struct InitChem *ICS)
//{
//	real S, alpha, perm;
//	alpha = 1.01325;
//
//	S = calc_surface_area(dx, por, ICS); /* m^2/L = 10^3 1/micron */
//
//	perm = alpha*por/(2*tau*S*S);
//	return perm;
//
//}

/* Surface area for one node -  1/micron, rho = kg/L,  --> S_p = (1-por)/por * S_g * rho */
//real calc_surface_area(real *dx, real por, struct InitChem *ICS)
//{
//	int i, pos;
//	real S;
//	S=0.;
//	for(i=0;i<ICS->size_rock;++i)
//	{
//		pos = ICS->pos_rock[i];
//		S += (ICS->c_buffer[pos]+dx[i])*ICS->Sg[pos]*ICS->SM_mineral->mol_weight[ICS->pos_buffer[i]]/ICS->SM_mineral->mol_volume[ICS->pos_buffer[i]];
//	}
//	for(i=0;i<ICS->size_sup_min;++i)
//		S += (ICS->c_sup_min[i]+dx[i])*ICS->Sg[i]*ICS->SM_mineral->mol_weight[ICS->pos_sup_min[i]]/ICS->SM_mineral->mol_volume[ICS->pos_sup_min[i]];
//	S *= (1-por)/por; /* rock to pore volume */
//	return S;
//}

//void set_mineral_conc(real rho_r, real por, struct BasVec *Vchem)
//{
//	int i, pos;
//
//	for(i=0;i<Vchem->ICS->size_rock;++i)
//	{
//		pos   = Vchem->ICS->pos_buffer[i];
//		Vchem->ctot_mineral[pos] *= rho_r*(1.-por)/por/Vchem->ICS->SM_mineral->mol_weight[pos];
//	}
//
//	for(i=0;i<Vchem->ICS->size_sup_min;++i)
//	{
//		pos   = Vchem->ICS->pos_sup_min[i];
//		Vchem->ctot_mineral[pos] *= rho_r*(1.-por)/por/Vchem->ICS->SM_mineral->mol_weight[pos];
//	}
//
//	if(Vchem->Vsurf != NULL)
//	{
//		for(i=0;i<Vchem->Vsurf->ICS->size_rock;++i)
//		{
//			pos   = Vchem->Vsurf->ICS->pos_buffer[i];
//			Vchem->Vsurf->ctot_mineral[pos] *= rho_r*(1.-por)/por/Vchem->Vsurf->ICS->SM_mineral->mol_weight[pos];
//		}
//		for(i=0;i<Vchem->Vsurf->ICS->size_sup_min;++i)
//		{
//			pos   = Vchem->Vsurf->ICS->pos_sup_min[i];
//			Vchem->Vsurf->ctot_mineral[pos] *= rho_r*(1.-por)/por/Vchem->Vsurf->ICS->SM_mineral->mol_weight[pos];
//		}
//	}
//
//}

/* note c_vchem[i] is in wt% and ctot_mineral in mol/l */
/* dx relative to ICS numbering                        */
//void calc_wt_change(real *dx, real rho_r, real por,  struct BasVec *Vchem)
//{
//	int i, pos, pos_b;
//
//	for(i=0;i<Vchem->ICS->size_rock;++i)
//	{
//		pos = Vchem->ICS->pos_buffer[i];
//		pos_b = Vchem->ICS->pos_rock[i];
//		dx[i] = -Vchem->ICS->c_buffer[pos_b] + Vchem->ctot_mineral[pos]*Vchem->ICS->SM_mineral->mol_weight[pos]*por/(1.-por)/rho_r;
//	}
//
//	for(i=Vchem->ICS->size_rock;i<Vchem->ICS->size_sup_min+Vchem->ICS->size_rock;++i)
//	{
//		pos = Vchem->ICS->pos_sup_min[i];
//		dx[i] = -Vchem->ICS->c_sup_min[i] + Vchem->ctot_mineral[pos]*Vchem->ICS->SM_mineral->mol_weight[pos]*por/(1.-por)/rho_r;
//	}
//
//
//}

/* returns total porosity change, and contribution from each mineral phase */
//real calc_por_change(real *dphi, real *dx, real rho_r, real por,  struct BasVec *Vchem)
//{
//	int i, pos;
//	real dphi_tot;
//
//	dphi_tot = 0.;
//	for(i=0;i<Vchem->ICS->size_rock;++i)
//	{
//		pos = Vchem->ICS->pos_buffer[i];
//		dphi[i] = -por*(1.-por)*rho_r*dx[i]*1.e-3*Vchem->ICS->SM_mineral->mol_volume[pos]/Vchem->ICS->SM_mineral->mol_weight[pos];
//		dphi_tot += dphi[i];
//	}
//	for(i=Vchem->ICS->size_rock;i<Vchem->ICS->size_sup_min;++i)
//	{
//		pos = Vchem->ICS->pos_buffer[i];
//		dphi[i] = -por*(1.-por)*rho_r*dx[i]*1.e-3*Vchem->ICS->SM_mineral->mol_volume[pos]/Vchem->ICS->SM_mineral->mol_weight[pos];
//		dphi_tot += dphi[i];
//	}
//
//	return dphi_tot;
//
//}

void calculate_mol_weight_mineral(struct ChemTable *Min, struct ChemTable *Bas)
{
	int i,j;

	if(Min->mol_weight == NULL) Min->mol_weight = (real *) calloc(Min->size[0],sizeof(real));

	for(i=0;i<Min->size[0];++i)
	{
		Min->mol_weight[i] = 0.;
		for(j=0;j<Min->size[1];++j) Min->mol_weight[i] += Min->M[i][j]*Bas->mol_weight[j];
	}

}

void calc_mineral_change_julius(struct BasVec *Vchem)
{
	real SI, sgn, fact, mol;
	int i, k, posi;
	fact = Vchem->RP->dt;
	for(i=0;i<Vchem->size_sup_min;++i)
	{
		posi = Vchem->pos_sup_min[i];
		mol = Vchem->RP->mol[i];
		
		SI =0.;
		for(k=0;k<Vchem->ICS->SM_mineral->size[1];++k) SI += Vchem->ICS->SM_mineral->M[posi][k]*Vchem->log_a[k];
		SI -= Vchem->ICS->SM_mineral->logK[posi];
		sgn = ( SI < 0 ? 1. : -1.);			

		Vchem->ICS->SM_mineral->log_a[posi] = SI;
		Vchem->ICS->SM_mineral->log_m[posi] = SI-Vchem->ICS->SM_mineral->log_g[posi];
		Vchem->ctot_mineral[posi] -= sgn*mol*fact*(Vchem->rate[i][0]+Vchem->rate[i][1]*pow10(Vchem->log_a[Vchem->pos_pH]))
			*pow(fabs(1.-pow10(SI*Vchem->rate[i][2])),Vchem->rate[i][3]); 
	}
}

/* 
make ICS that have all species as equilibrium species
*/


//real equilibrate_pore_water(real *ctot, struct BasVec *Vchem)
//{
//	struct BasVec *Vchem_new=NULL;
//	int size_sup,i,j;
//	real c_h;
//
//	Vchem_new = (struct BasVec *) calloc(1,sizeof(struct BasVec));
//	make_basvec_struct(Vchem->ICS,Vchem_new);
//	if(Vchem_new->pos_X > -1) Vchem_new->equilibrate = 1;
//
//	 /* initialize */
//	Vchem_new->Temp = Vchem->Temp;
//	for(i=0;i<Vchem->size;++i)
//	{
//		Vchem_new->ctot[i]=Vchem->ctot[i];
//		Vchem_new->ctot_ads[i]=Vchem->ctot_ads[i];
//		Vchem_new->log_a[i]=Vchem->log_a[i];
//		Vchem_new->log_m[i]=Vchem->log_m[i];
//		Vchem_new->log_g[i]=Vchem->log_g[i];
//	}
//
//
//
//	for(i=0;i<Vchem->size_sup_min;++i)
//	{
//		add_new_buffer_mineral(Vchem_new->pos_sup_bas[i], Vchem_new->pos_sup_min[i],Vchem_new);
//		printf("Add buffer mineral %s with basis %s\n",
//			Vchem_new->ICS->SM_mineral->row_name[Vchem_new->pos_sup_min[i]],
//			Vchem_new->ICS->SM_basis->row_name[Vchem_new->pos_sup_bas[i]]);
//	}
//
//	/* add sup sat min to rock buffer for fast equilibration calc */
//	size_sup = Vchem_new->size_sup_min;
//
//	Vchem_new->size_sup_min = 0;
//		/* basis transformation matrix */
//
//	/* calculate equilibrium */
//	set_temperature_db(Vchem_new);
//	Vchem_new->equilibrate = 1; /* important for ion exchange */
//	solve_chemistry(Vchem_new);
//	calc_ctot_aq(Vchem_new);
//	if(Vchem->pos_X > -1) /* calc adsorbes species */
//	{
//		for(i=0;i<Vchem->size;++i)
//		{
//			Vchem->ctot_ads[i]=0.;
//			for(j=0;j<Vchem->ICS->SM_all->size[0];++j)
//			{
//				if(Vchem->ICS->SM_all->type[j] != 0)/*surface species*/
//					Vchem->ctot_ads[i] += Vchem->ICS->SM_all->M[j][i]*pow10(Vchem->ICS->SM_all->log_m[j]);
//			}
//		}
//	}
//
//	c_h = Vchem_new->ctot_calc[Vchem_new->pos_pH];
//
//	for(i=0;i<Vchem->ICS->size;++i) ctot[i] = Vchem_new->ctot_calc[Vchem_new->ICS->pos[i]];
//	/* set back and free struct*/
//	Vchem_new->size_sup_min=size_sup;
//	free_BasVec(Vchem_new);
//	free(Vchem_new);
//
//	return c_h;
//}

void add_H_to_mbal(real c_h, struct BasVec *Vchem)
{
	Vchem->pos_mass[Vchem->size_mass] = Vchem->pos_pH;
	Vchem->ctot[Vchem->pos_pH]        = c_h;
	Vchem->size_mass++;
	if(Vchem->Vsurf != NULL)
	{
		Vchem->Vsurf->pos_mass[Vchem->Vsurf->size_mass] = Vchem->Vsurf->pos_pH;
		Vchem->Vsurf->ctot[Vchem->Vsurf->pos_pH]        = c_h;
		Vchem->Vsurf->size_mass++;
	}
}

void add_H_to_ICS(struct InitChem *ICS, int posH)
{
	real *ctmp;
	int *pos, *posm, i, flag;
	flag=0;
	/* check if H+ is member */
	for(i=0;i<ICS->size;++i) if(ICS->pos[i] == posH) flag = 1;

	if(!flag)
	{
		ctmp = (real *) calloc(ICS->size+1,sizeof(real));
		pos  = (int  *) calloc(ICS->size+1,sizeof(int));
		posm = (int  *) calloc(ICS->size_mass+1,sizeof(int));

		for(i=0;i<ICS->size;++i) pos[i]=ICS->pos[i];
		for(i=0;i<ICS->size_mass;++i) posm[i]=ICS->pos_mass[i];

		pos[ICS->size] = posH;
		posm[ICS->size_mass] = posH;
		ICS->size++;
		ICS->size_mass++;
		free(ICS->pos);free(ICS->pos_mass);
		ICS->pos = pos;
		ICS->pos_mass = posm;
	}
}

/* returns volumetric strain in % */
/* eps = 10**b/(m+1)*time**(m+1) + c */
/* t in timer */
//real anal_perm_mod(real time, real b, real m, real c)
//{
//	return pow10(b)/(m+1)*pow(time,m+1)+c;
//}
//
