#include "chem_global.h"


/* dc/dt = -k*(c-c_eq) */
/* d (c-c_eq)/dt = -k(c-c_eq) => c(t)-c_eq = (c(0)-c_eq)*exp(-k*t) */
void reaction_rate_lb(real *c_updated, struct BasVec *Vchem)
{
	int pos_j,j;
	real k;

	printf("REACTION RATE LB\n");
	k =0.1;
	k=3.4e-3;
	k=999999;

	for(j=0;j<Vchem->size_rock;++j)
	  {
	    pos_j = Vchem->pos_rock[j];
	    c_updated[pos_j] += (Vchem->ctot_calc[pos_j]-Vchem->ctot[pos_j]) - exp(-k)*(Vchem->ctot_calc[pos_j]-Vchem->ctot[pos_j]);
	  }
	
	for(j=0;j<Vchem->size_mass;++j)
	  {
	    pos_j = Vchem->pos_mass[j];
	    c_updated[pos_j] += (Vchem->ctot_calc[pos_j]-Vchem->ctot[pos_j]) - exp(-k)*(Vchem->ctot_calc[pos_j]-Vchem->ctot[pos_j]);
	  }
}

void reaction_rate_julius(real *c_updated, struct BasVec *Vchem)
	{
		int j,pos_j, pos_k;
		real Ea = 15e+3/0.238902957619; /* 1J = 0.238902957619 cal conversion between kcal J */
		real T0 = 273.15+130;
		real k, k1, k2, k3, ctot_mag;
		real eq_flag;

		/* Y2K acid*/
		k1 = 3e-3;

		k1 = 1e-5;

		k1 = 5e-7;


	printf("REACTION RATE JULIUS\n");

		if(Vchem->equilibrate == 1) eq_flag = 0.;/* No rate equation */
		else eq_flag = 1.; 

		/*		eq_flag = 0.; */

		k1 = 1E-3; k2 = 0.4; k3 =15;
/*		ctot_mag = Vchem->ctot_mineral[9]+Vchem->ctot_mineral[16] ;*/
		ctot_mag = Vchem->ctot_mineral[60];
		k = 2e-5*exp(-Ea/Vchem->ICS->CP->Rg*(1./Vchem->Temp-1./T0))*(k2+k1*k3*ctot_mag)/(1+k3*ctot_mag)*Vchem->RP->dt/Vchem->RP->porosity;
		k = 2e-5*exp(-Ea/Vchem->ICS->CP->Rg*(1./Vchem->Temp-1./T0))*(k2/1.5+40*k1*k3*ctot_mag)/(1+1.2*k3*ctot_mag)*Vchem->RP->dt/Vchem->RP->porosity;
       
/*		k = 1e-6*exp(-Ea/CP_->Rg*(1./Vchem->Temp-1./T0))*Vchem->RP->dt/Vchem->RP->porosity;*/

/*		k = 1./(9.49*24*60*60)*Vchem->RP->dt;*/
/*		k*=100;*/
/*		k=1e-6*Vchem->RP->dt/Vchem->RP->porosity;*/
		/* RATE 3PV day */
		/*	      	k = .8e-6*exp(-Ea/Vchem->ICS->CP->Rg*(1./Vchem->Temp-1./T0))*Vchem->RP->dt/Vchem->RP->porosity;*/


		for(j=0;j<Vchem->size_rock;++j)
		{
			pos_j = Vchem->pos_rock[j];
			pos_k = Vchem->pos_buffer[j];
			c_updated[pos_j] += Vchem->ctot_calc[pos_j]-Vchem->ctot[pos_j]-eq_flag*(Vchem->ctot_calc[pos_j]-Vchem->ctot[pos_j])*exp(-k);
/*
			c_updated[pos_j] += k*(Vchem->ctot_calc[pos_j]-Vchem->ctot[pos_j]);
			*/
			Vchem->ctot_mineral[pos_k] += Vchem->delta_mineral[pos_j] - eq_flag*Vchem->delta_mineral[pos_j]*exp(-k);
		}

	
		for(j=0;j<Vchem->size_mass;++j)
		{
			pos_j = Vchem->pos_mass[j];
			c_updated[pos_j] += Vchem->ctot_calc[pos_j]-Vchem->ctot[pos_j]-eq_flag*(Vchem->ctot_calc[pos_j]-Vchem->ctot[pos_j])*exp(-k);
		}

	}

