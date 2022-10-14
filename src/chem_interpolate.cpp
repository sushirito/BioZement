#include "chem_global.h"
#include "macros.h"

/* generate binary number */

void gen_binary(int *b, int num, int length)
{
  int i;
  for(i=0;i<length;++i)
    {
      b[length-i-1]=num % 2;
      num /= 2;
    }
}

void gen_vchem_vector(int no_comb, struct BasVec *Vchem, struct InitChem *ICS)
{
  //int no_comb;
  int i, j,*b;
  real *cmin;

  b = (int *) malloc(sizeof(int)*ICS->size_sup_min);

  cmin = (real *) malloc(sizeof(real)*ICS->size_sup_min);
  for(i=0;i<ICS->size_sup_min;++i) cmin[i]=ICS->c_sup_min[i];

  //no_comb = (int) pow(2,ICS->size_sup_min);

  for(i=0;i<no_comb;++i)
    {
      gen_binary(b,i,ICS->size_sup_min);
      for(j=0;j<ICS->size_sup_min;++j) ICS->c_sup_min[j] = (real) b[j];
      make_basvec_struct(ICS,&(Vchem[i]));
    }

  for(i=0;i<ICS->size_sup_min;++i) ICS->c_sup_min[i]=cmin[i];
  free(cmin);
  free(b);
}

int gen_new_key(real *ctot_min, int size)
{
  int i, key;
  
  //#ifdef _MAGN_ON_MAGN_
  //if ( (ctot_min[1]>opt->magn_limit) || (ctot_min[2]>0.0))
  //  return 4;
  //#endif

  key=0;
  for(i=0;i<size;++i)
    key = key*2 + (ctot_min[i]> 0. ? (1):(0) );
  
  return key;
}

 void copy_BasVec(struct BasVec *V_out, struct BasVec *V_in)
 {
	 int i;

	 for(i=0;i<V_in->size;++i)
	 {
		 V_out->ctot[i] = V_in->ctot[i];
		 V_out->ctot_calc[i] = V_in->ctot_calc[i];;
	 }
	 for(i=0;i<V_in->ICS->SM_mineral->size[0];++i)
	 {
		 V_out->ctot_mineral[i] = V_in->ctot_mineral[i];
	 }
 }

 /* 
 input dc : conc flux and Vchem
 output dc:_min : mineral flux 
 */
 void calc_min(real *dc, real *dc_min, struct BasVec *Vchem)
 {
	 int i,j,pos;
	 real c_tmp;

	 for(i=0;i<Vchem->ICS->size;++i) Vchem->ctot[Vchem->ICS->pos[i]]=dc[i];
	 for(i=0;i<Vchem->ICS->size_sup_min;++i) Vchem->ctot_mineral[Vchem->ICS->pos_sup_min[i]]=0.;
	 

	 for(i=0;i<Vchem->size_sup_min;++i)
	 {
		 pos = Vchem->pos_sup_min[i];
		 c_tmp=0.;
		 for(j=0;j<Vchem->size_sup_min;++j) c_tmp += Vchem->sm_sup_buf[i][j]*Vchem->ctot[Vchem->pos_sup_bas[j]] ;
		
		 Vchem->ctot_mineral[Vchem->pos_sup_min[i]] = c_tmp;
	 }

	 for(i=0;i<Vchem->ICS->size_sup_min;++i) dc_min[i] = Vchem->ctot_mineral[Vchem->ICS->pos_sup_min[i]];
 
 }

 void update_BasVec(real *c_new, real *dc_min, struct BasVec *Vchem)
 {
	 int i;
	 for(i=0;i<Vchem->ICS->size;++i) Vchem->ctot[Vchem->ICS->pos[i]] = c_new[i];
	 for(i=0;i<Vchem->ICS->size_sup_min;++i)	 Vchem->ctot_mineral[Vchem->ICS->pos_sup_min[i]] += dc_min[i];
	 
 }
