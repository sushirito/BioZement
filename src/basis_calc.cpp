#include "global.h"


/* ---------------------------------------------------------------------------------  check_basis */
int check_basis(int *w, struct System *sys)
/*
check_basis:
	checks that the symmetry demands are fulfilled

  INPUT  : 1) w: the weights where w[0] is the grates common divisor for the denominator
                 and w[1] to w[n_Q] are the numerators
	       2) ei : relative position of neighbors

  OUTPUT : (int) returns the lowest order of the symmetry that is broken ([0,..5])
                 returns 6 if OK
*/
{
	int od, cd[5], ck; /* counter order, dimension, direction */
	int cd_lim, cdl, cod;
	int cs_2, cs_4,w_tot, cdl_tmp;
	int cn_tmp, sum_tmp;
	int *ei = sys->ei;
	int n_D = sys->n_D;
	int n_Q = sys->n_Q;

	/* Find the numerator for the lattice velocity cs_2 */
	cs_2 = 0;
	for( ck = 0; ck < n_Q; ++ck )
		for( cdl = 0; cdl < n_D; ++cdl )
		cs_2 += w[ck+1] * ei[n_D*ck + cdl] * ei[n_D*ck + cdl];
	
	cs_2 /= n_D;
	cs_4 = (cs_2*cs_2)/w[0];

	/* Check weights */
	w_tot = 0;
	for( ck = 0; ck < n_Q; ++ck )
		w_tot += w[ck+1];

	if( w_tot != w[0] )
	{
		printf("ERROR: sum_i w_i = %lf ( = 1.0)\n", (double) (1.*w_tot)/w[0] );
		return BAD;
	}
	/* Check symmetries */
	cd_lim = 1;
	for( od = 1; od < 6; ++od )
	{
		cd_lim *= n_D;
		for( cdl = 0; cdl < cd_lim; ++cdl)
		{
			cdl_tmp = cdl;
			cd[0]   = cdl_tmp % n_D;
			for( cod = 1; cod < od; ++cod )
			{
				cdl_tmp -= cd[cod-1];
				cdl_tmp /= n_D;
				cd[cod] = cdl_tmp % n_D;
			}

			/* Sum */
			sum_tmp = 0;
			for( ck = 0; ck < n_Q; ++ck )
			{
				cn_tmp  = 1;
				for( cod = 0; cod < od; ++cod)
					cn_tmp *= ei[n_D*ck + cd[cod]];

				sum_tmp += w[ck+1]*cn_tmp;
			}

			switch( od ) 
			{
			case 2:
				if( sum_tmp != cs_2*dirac_delta(cd[0],cd[1]) )
				{
					printf("ERROR in %d order symmetry for the basis\n", od);
					return od;
				} else {
					//printf("%d order symmetry OK\n", od);
				}
				break;
			case 4:
				if( sum_tmp != cs_4*( dirac_delta(cd[0],cd[1])*dirac_delta(cd[2],cd[3]) + 
					                  dirac_delta(cd[0],cd[2])*dirac_delta(cd[1],cd[3]) +
								      dirac_delta(cd[0],cd[3])*dirac_delta(cd[1],cd[2]) ) )
				{
					printf("ERROR in %d order symmetry for the basis\n", od);
					return od;
				} else {
					//printf("%d order symmetry OK\n", od);
				} 
				break;
			default :
				if( sum_tmp != 0 )
				{
					printf("ERROR in %d order symmetry for the basis\n", od);
					return od;
				} else {
					//printf("%d order symmetry OK\n", od);
				}
			}
		}
	}

	
	return od;
}


/* ---------------------------------------------------------------------------------  init_k_bb */
void init_k_bb(struct System *sys)
/*
init_k_bb:
	set the bounce back directions

  INPUT  : 1) ei : relative postion of neighbors
		   2) ck_bb :bounce back directions
		   
  OUTPUT : void
*/
{
	int ck, ck1; /*vel. direction and bounce back direction */
	int cd;        /* direction counter */
	int test;      /* boolean test */
	
	sys->k_bb = (int *) calloc(sys->n_Q, sizeof(int));

	int *k_bb = sys->k_bb;
	int *ei = sys->ei;
	int n_D = sys->n_D;
	int n_Q = sys->n_Q;

	for( ck = 0; ck < n_Q; ++ck )
	{
		for( ck1 = 0; ck1 < n_Q; ++ck1 )
		{
			test = TRUE;
			for( cd = 0; cd < n_D; ++cd )
			{
				if( ei[n_D*ck + cd] != -1*ei[n_D*ck1 + cd] )
				{
					test = FALSE;
					break;
				}
			}

			if( test )
			{
				k_bb[ck] = ck1;
				break;
			}
		}
	}
}


/* ---------------------------------------------------------------------------------  calc_cs_2 */
/* real calc_cs_2(int *w, struct System *sys) */
void calc_cs_2(int *w, struct System *sys)
/*
calc_cs_2 :
	calculates the square of the lattice velocity

  INPUT  : 1) w: the weights where w[0] is the grates commen divisor for the denominator
                 and w[1] to w[n_Q] are the numenators
	       2) ei : relative postion of neighbors

  OUTPUT : (real) value of c_s^2
*/
{
	int cd, ck; /* counter dimension, direction*/
	real cs_2;

	int *ei = sys->ei;
	int n_D = sys->n_D;
	int n_Q = sys->n_Q;

	/* Find the lattice velocity sys->cs_2 */
	cs_2 = 0.0;
	for( ck = 0; ck < n_Q; ++ck ) 
		for( cd = 0; cd < n_D; ++cd )
		cs_2 += w[ck+1]*ei[n_D*ck + cd]*ei[n_D*ck + cd];


	cs_2 /= w[0]*n_D;
	sys->cs_2 = cs_2;
		
	return;
}

/* ---------------------------------------------------------------------------------  calc_fg_eq_c */
void calc_fg_eq_c(struct System *sys)
/*
calc_fg_eq_c :
	calculates the constants in the equlibirum distributions

  INPUT  : 1) cs_2 : the square of the lattice velocity

  OUTPUT : void
  
*/
{
        real cs_2 = sys->cs_2;
	
	sys->eq_cst[0] = 1./cs_2;
	sys->eq_cst[1] = 0.5/(cs_2*cs_2);
	sys->eq_cst[2] = -0.5/cs_2; 
	return;
}
