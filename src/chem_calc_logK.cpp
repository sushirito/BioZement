#include "chem_global.h"

void calc_logK(struct ChemTable *tab_all)
{

  tab_all->logK = pick_col(tab_all,"logK");
}

/* temperature points are stored in tab_all->delta */
void get_logK(real Temp, real *logK, struct ChemTable *tab_all)
{
	int x, i;
	real left, right;
	/* check if Temp is equal to tabulated values */
	for(i=0;i<tab_all->size[1];++i)
	{
		if(Temp == tab_all->T[i])
		{
			for(x=0;x<tab_all->size[0];++x) logK[x] = tab_all->M[x][i];
			return; /* no interpolation */
		}
	}
	/* start interpolation */
	right = left = tab_all->T[0];
	x=0;
	while( Temp > right)
	{
		x++;
		right = tab_all->T[x];
		left  = tab_all->T[x-1];
	}

	if(i > tab_all->size[1])
	{
		printf("Temperature outside range for logK values\n");
		printf("Temp: %g Range: %g-%g\n", Temp, tab_all->T[0],tab_all->T[tab_all->size[1]-1]);
		exit(0);
	}

	for(i=0;i<tab_all->size[0];++i) 
		logK[i] = tab_all->M[i][x-1] + (tab_all->M[i][x]-tab_all->M[i][x-1])/(right-left)*(Temp-left);
	return;
}
