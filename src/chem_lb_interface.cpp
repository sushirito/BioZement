#include "chem_global.h"
#include "global.h"

void chem_calc_eq_rho_lb(real *rho_eq_wall, real *rho_wall , int cs, struct BasVec *Vchem,
    real rate_constant)
/*
chem_calc_eq_rho_lb :
	Calculates the equilibrium concentrations based on the wall concentrations

  INPUT : 1) rho_eq_wall : equilibrium concentration (to be calculated)
		  2) rho_wall    : Concentration (supplied)

  OUTPUT : void
*/
{
	/* int cc, num_itr; /\* Counter species *\/ */
	/* int pos,i; */
        int i;
	real *c_tot;

	c_tot = (real *) calloc(Vchem[cs].ICS->size,sizeof(real));

	for(i=0;i<Vchem[cs].ICS->size;++i) c_tot[i]=rho_wall[i] ;

	solve_chem_bulk_surface(c_tot,0,&Vchem[cs]);

	for(i=0;i<Vchem[cs].ICS->size;++i)
	{
		rho_eq_wall[i]=c_tot[i];
	/*	rho_eq_wall[i]=1e-9; */

	}


/*printf("Wall node : %d and pH %lf and surface charge: %lf\n", cs, -Vchem_[cs].log_a[0], psi_ );*/


	free(c_tot);
}

/* ---------------------------------------------------------------------------------  chem_calc_eq_rho_lb */
int chem_calc_eq_rho_lb_moving(real *rho_eq_wall, real *rho_wall, struct InitChem *ICS,
                               real *min_conc, struct treeNode **tn_tmp, struct splayTree **st_lst, struct System *sys)
/* 
chem_calc_eq_rho_lb :
	Calculates the equilibrium concentrations based on the wall concentrations

  INPUT : 1) rho_eq_wall : equilibrium concentration (to be calculated)
		  2) rho_wall    : Concentration (supplied)

  OUTPUT : void
*/
{
  int i, new_key, key;
  struct splayTree * st_tmp;
  //struct treeNode * tn_tmp;
  
  
  new_key = gen_new_key(min_conc,ICS->size_sup_min);
  
  
  do {
    key = new_key;
    st_tmp = st_lst[key];
    st_tmp->Vchem_key->RP->dt = 1.;
    st_tmp->Vchem_key->RP->porosity=1.;
    st_tmp->Vchem_key->Temp = ICS->Temp;
      
#ifdef _MAGN_ON_MAGN_
    //st_tmp->Vchem_key->RP->mol[0] = 1.0;
    st_tmp->Vchem_key->RP->mol[1] = 0.0;
    if ((min_conc[1]>sys->opt->magn_limit) || (min_conc[2]>0.0)) {
      //if (key==4) {
      st_tmp->Vchem_key->RP->mol[1] = 1.0;
      //printf("min_conc = %.15e, %.15e\n", min_conc[0], min_conc[1]);
    }
#endif
    
   // if (ICS->INTERPOLATE) { /* AMIN no inperpolation*/
    if (0) {
      *tn_tmp = FindInsert(rho_wall, NULL, st_tmp, min_conc, st_tmp->Vchem_key->ICS->size_sup_min);
    } else {
      //if (_DEBUG_FLAG_) printf("\nFindInsertNoInterp\n");
      *tn_tmp = FindInsertNoInterp(rho_wall, NULL, st_tmp, min_conc, st_tmp->Vchem_key->ICS->size_sup_min);
    }
    new_key = (*tn_tmp)->min_key;
  } while(new_key != key);
  /* We store the change in concentration, to avoid changes in conservative species such as Cl-*/
  /* C_IN-C_OUT is stored in tn_tmp->val */
  
  /* check if some conc are negative */
  for(i=0;i<ICS->num_phases ;++i) {
    //if (rho_wall[i]-(*tn_tmp)->val[i]<-1e-5)
    if (rho_wall[i]-(*tn_tmp)->val[i]<0) {
      printf("WARNING: Negative %s: %g, rho_wall: %g\n", ICS->SM_basis->row_name[ICS->pos[i]], rho_wall[i]-(*tn_tmp)->val[i], rho_wall[i]);
    }
    rho_eq_wall[i]= ( (rho_wall[i]-(*tn_tmp)->val[i]) < 0 ? (1e-8):(rho_wall[i]-(*tn_tmp)->val[i]) );
  }
  
  return new_key;
}

