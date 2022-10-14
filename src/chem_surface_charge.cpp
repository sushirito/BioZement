#include "chem_global.h"

//double surface_charge(struct BasVec *Vchem)
//{
//	double lower, upper, tol, value, surf_pot;
//	int i;
///* calculate surface charge */
//
//	Vchem->sch=0.;
//
//	/* basis species */
//	for(i=0;i<Vchem->size;++i) Vchem->sch += Vchem->ICS->SM_basis->scharge[i]*pow10(Vchem->log_m[i]);
//	/* surface complexes */
//	for(i=0;i<Vchem->ICS->SM_all->size[0];++i) Vchem->sch += Vchem->ICS->SM_all->scharge[i]*pow10(Vchem->ICS->SM_all->log_m[i]);
//
//	Vchem->sch *= Vchem->ICS->CP->F/Vchem->ICS->CP->SA;
//	if(Vchem->sch < 0.)
//	{
//		lower = -1.;
//		upper = 0.;
//	} else
//	{
//		lower = 0.;
//		upper = 1.;
//	}
//	tol = 1e-8;
//	if(fabs(Vchem->sch) >= tol){
//		surf_pot = zeroin( 0,0,lower, upper,&grahame,tol, Vchem);
//		value=grahame(surf_pot,0, Vchem);
//	} else
//		surf_pot = 0.;
//	return surf_pot;
//
//}
/* returns the Grahame equation for a given surface potential - x */
double grahame(double x, int no_call, struct BasVec *Vchem)
{
	double gr;
	int i;

	gr=0.;
	/* sum up basis species */
	for(i=0;i<Vchem->size;++i) gr+= pow10(Vchem->log_m[i])*(exp(-Vchem->ICS->SM_basis->charge[i]*Vchem->ICS->CP->beta_chem*x)-1.);
	
	/* sum up complexes */
	for(i=0;i<Vchem->ICS->SM_all->size[0];++i) gr += pow10(Vchem->ICS->SM_all->log_m[i])*(exp(-Vchem->ICS->SM_all->charge[i]*Vchem->ICS->CP->beta_chem*x)-1.);
/* factor 2*e*e_0*k_B*1000*N_a = 1.1776E-05 when concentrations are given in mol/l*/
	/* 2*e*e_0*k_B for concentrations in number/m^3 */
	/* k_B = R/Na --> 2*e*e_0*k_B*1000*N_a = 2e*e0*R*1000*/
	gr *= 1.1778e-5*Vchem->Temp;
	gr -= Vchem->sch*Vchem->sch;

	return gr;
}


double surface_charge_newton(struct BasVec *Vchem)
{
	double tol, df_value, fvalue;
	double x_init, x_old, x_new;
	int i, iter;
/* calculate surface charge */
	Vchem->sch=0.;

	/* basis species */
	for(i=0;i<Vchem->size;++i) Vchem->sch += Vchem->ICS->SM_basis->scharge[i]*pow10(Vchem->log_m[i]);
	/* surface complexes */
	for(i=0;i<Vchem->ICS->SM_all->size[0];++i) Vchem->sch += Vchem->ICS->SM_all->scharge[i]*pow10(Vchem->ICS->SM_all->log_m[i]);

	Vchem->sch *= Vchem->ICS->CP->F/Vchem->ICS->CP->SA;

	x_init = Vchem->sch/10.; /* initial value */
	iter = 0;
	
	fvalue = grahame_newton(x_init, &df_value, Vchem);
	x_new=x_old = x_init;
	tol = 1e-16;
	while(fabs(fvalue) > tol)
	{
		x_new = x_old - fvalue/df_value;
		fvalue = grahame_newton(x_new, &df_value, Vchem);
		x_old = x_new;
		iter++;
	}

	printf("surface charge newtons method: \n");
	printf("No iter %d and value %g : \n", iter, fvalue);
	return x_new;
}
double grahame_newton(double x, double *df, struct BasVec *Vchem)
{
	double gr, df_tmp, exp_tmp, mc_tmp;
	int i;

	gr=0.;
	df_tmp=0.;
	/* sum up basis species */
	for(i=0;i<Vchem->size;++i)
	{
		mc_tmp  = pow10(Vchem->log_m[i]);
		exp_tmp = exp(-Vchem->ICS->SM_basis->charge[i]*Vchem->ICS->CP->beta_chem*x); 
		gr     += mc_tmp*(exp_tmp-1.);
		df_tmp -= mc_tmp*Vchem->ICS->CP->beta_chem*Vchem->ICS->SM_basis->charge[i]*exp_tmp; 
	}
	
	/* sum up complexes */
	for(i=0;i<Vchem->ICS->SM_all->size[0];++i)
	{
		mc_tmp  = pow10(Vchem->ICS->SM_all->log_m[i]);
		exp_tmp = exp(-Vchem->ICS->SM_all->charge[i]*Vchem->ICS->CP->beta_chem*x); 
		gr += mc_tmp*(exp_tmp-1.);
		df_tmp -= mc_tmp*Vchem->ICS->CP->beta_chem*Vchem->ICS->SM_all->charge[i]*exp_tmp;
	}

	df_tmp *= 1.1778e-5*Vchem->Temp;
	gr *= 1.1778e-5*Vchem->Temp;
	gr -= Vchem->sch*Vchem->sch;
	
	*df = df_tmp;
	return gr;
}

real trapez(real (*func)(real, int, struct ChemTable *, struct BasVec *), real a, real b, int pos, struct ChemTable * Tab, struct BasVec *Vchem)
/*
Returns the integral of the function func between a and b, by trapez method:
int_a^b f(x) = (b-a)/(2 N) ( f(a) +f(b) + 2*sum_i=1^N-1 f(a+i*(b-a)/N))
*/
{
	int j, N=10000;
	real xm,xr,s;

	xr=(b-a)/((real) N);
	xm=a;
	s =  (*func)(a,pos,Tab,Vchem)+(*func)(b, pos, Tab,Vchem);
	for (j=1;j<=N;j++) 
	{
		xm += xr;
		s += 2*(*func)(xm,pos,Tab,Vchem);
	}
	return s *= xr*.5; /*Scale the answer to the range of integration.*/
}

real qgaus(real a, real b, int pos, struct ChemTable * Tab, struct BasVec *Vchem)
/*
Returns the integral of the function func between a and b, by ten-point Gauss-Legendre integration:
the function is evaluated exactly ten times at interior points in the range of integration.
*/
{
	int j;
	real xr,xm,dx,s;
	/*The abscissas and weights. First value of each array not used.*/
	real x[]={0.0,0.1488743389,0.4333953941, 
		0.6794095682,0.8650633666,0.9739065285};
	real w[]={0.0,0.2955242247,0.2692667193,
		0.2190863625,0.1494513491,0.0666713443};
	xm=0.5*(b+a);
	xr=0.5*(b-a);
	s=0; /*Will be twice the average value of the function, since the
			ten weights (five numbers above each used twice) sum to 2. */
	for (j=1;j<=5;j++) 
	{
		dx=xr*x[j];
		s += w[j]*(grahame_DL(xm+dx,pos,Tab,Vchem)+grahame_DL(xm-dx, pos, Tab,Vchem));
		//s += w[j]*((*func)(xm+dx,pos,Tab,Vchem)+(*func)(xm-dx, pos, Tab,Vchem));
	}
	return s *= xr; /*Scale the answer to the range of integration.*/
}
/* 
Estimate of the concentration in the diffusive layer, from:
Borkovec M, Westall J (1983) Solution of the poisson-boltzmann equation for surface
excesses of ions in the diffuse layer at the oxide-electrolyte interface. J Electroanal
Chem 150:325�337

xi = exp(-Vchem->ICS->CP->beta_chem*psi)
*/
real grahame_DL(real xi, int pos, struct ChemTable *Tab, struct BasVec *Vchem)
{
	real gr;
	int i;
	struct ChemTable *SM_mineral, *SM_all, *SM_basis;
	SM_mineral = Vchem->ICS->SM_mineral;
	SM_all     = Vchem->ICS->SM_all;
	SM_basis   = Vchem->ICS->SM_basis;

	if(xi == 1) return 0;
	else
	{
		gr=0.;
		/* sum up basis species */
		for(i=0;i<Vchem->size;++i) gr+= pow10(Vchem->log_m[i])*(pow(xi,SM_basis->charge[i])-1.);
		/* sum up complexes */
		for(i=0;i<SM_all->size[0];++i) gr += pow10(SM_all->log_m[i])*(pow(xi,SM_all->charge[i])-1.);
		gr *= xi*xi;
		if(gr<0) 
		{
			gr = 0.;
			if(Vchem->ICS->PRINT_DEBUG_CHEM) printf("WARNING: Negative value in routine: grahame_DL\n");
		}
		else
			gr = (pow(xi,Tab->charge[pos])-1.)/sqrt(gr);

		return gr;
	}
}

/*
calculates diffuse layer concentration, at exit Vchem->ctot_dl contains diffuse
layer concentration
*/
void calc_diffuse_layer_conc(real *c_dl, struct BasVec *Vchem)
{
	int i, j;
	real x_i, x_f, sign_xi, E, int_gr;
	real kappa_s, *cc_DL, g_tmp;

	struct ChemTable *SM_all, *SM_basis;
	SM_all     = Vchem->ICS->SM_all;
	SM_basis   = Vchem->ICS->SM_basis;


	cc_DL = (real *) calloc(SM_all->size[0],sizeof(real));
	


	/*fp = my_fopen("dl_debug.out","w");*/
	
	E = pow10(Vchem->log_m[Vchem->pos_exp]);

	x_i = 1;
	x_f = 1./E;

	sign_xi = (x_f > 1 ? (1.):(-1.));
	/* 1.1778515e-08 = sqrt(e*e_0*R/2)/F */
	kappa_s = sign_xi*Vchem->ICS->CP->SA*sqrt(.5*Vchem->ICS->CP->e0*Vchem->ICS->CP->ew*Vchem->ICS->CP->Rg*Vchem->Temp*1e+3)/Vchem->ICS->CP->F;
	
	/* surface excess of the complexes */
	for(j=0;j<SM_all->size[0];++j)
	{
		cc_DL[j] = 0.;
		if(SM_all->charge[j] != 0)
		{
		  cc_DL[j] = kappa_s*qgaus(x_i, x_f, j, SM_all, Vchem)*pow10(SM_all->log_m[j]);
		  //cc_DL[j] = kappa_s*qgaus(grahame_DL, x_i, x_f, j, SM_all,Vchem)*pow10(SM_all->log_m[j]);
		}
	}

	for(i=0;i<Vchem->size;++i)
	{
	
		c_dl[i] = 0.;
		if(SM_basis->charge[i]!= 0)
		{
		  int_gr = qgaus(x_i, x_f, i, SM_basis, Vchem)*pow10(Vchem->log_m[i]);
		  //int_gr = qgaus(grahame_DL, x_i, x_f, i, SM_basis,Vchem)*pow10(Vchem->log_m[i]);
		  /*		int_gr = trapez(grahame_DL, x_i, x_f, i, &SM_basis,Vchem)*pow(10,Vchem->log_m[i]);*/
		  c_dl[i] = kappa_s*int_gr;
		}
		g_tmp =0.;
		for(j=0;j<SM_all->size[0];++j) g_tmp += SM_all->M[j][i]*cc_DL[j];
		c_dl[i] += g_tmp;
	}


	/* calculate complex contribution */


	
/*
	delta = (x_f-x_i)/1000;
	x = x_i;
	fprintf(fp,"x\t %s\t %s\n",SM_basis->row_name[Vchem->pos_mass[0]],SM_basis->row_name[Vchem->pos_mass[1]]);
	for(i=0;i<1001;++i)
	{
		fprintf(fp,"%g\t%g\t%g\n",x, grahame_DL(x, Vchem->pos_mass[0], &SM_basis, Vchem)*pow(10,Vchem->log_m[Vchem->pos_mass[0]]),
			grahame_DL(x, Vchem->pos_mass[1], &SM_basis,Vchem)*pow(10,Vchem->log_m[Vchem->pos_mass[1]]));
		x += delta;
	}
	fclose(fp);
*/
	free(cc_DL);
}

