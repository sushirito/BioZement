#include "chem_global.h"

void set_chem_param(struct ChemTable *tab_all)
{
	tab_all->a0      = pick_col(tab_all,"a0");
	tab_all->charge  = pick_col(tab_all,"charge");
	tab_all->scharge = pick_col(tab_all,"scharge");
}
