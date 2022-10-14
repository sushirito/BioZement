#include "chem_global.h"

void calc_complex(struct ChemTable *sm, struct BasVec *Vchem)
{
  int i,j;
  /* calculate concentration of complexes */
  for(i=0;i<sm->size[0];++i)
    {
      sm->log_a[i]=sm->log_m[i]=0.;

      for(j=0;j<sm->size[1];++j) sm->log_a[i] += sm->M[i][j]*Vchem->log_a[j];
	
      sm->log_a[i] -= sm->logK[i];
	  sm->log_m[i] =  sm->log_a[i]-sm->log_g[i];
	 
	  /* special treatment of exchange species */
	  if(sm->type[i] == 2) sm->log_m[i] += -log10(sm->charge[i])+log10(Vchem->ctot[Vchem->pos_X]);
    }
  
}


