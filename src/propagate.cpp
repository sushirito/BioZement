#include "global.h"

/* ---------------------------------------------------------------------------------  propagate */
/* void propagate(real *fg, real *fg_tmp, int n_c) */
void propagate(struct Node *node, struct Field *field, struct System *sys)
/* 
propagate :
	propagates the distribution 	

  INPUT : 1) fg     : pointer to the distribution
          2) fg_tmp : pointer to the temporal distribution
          3) n_c    : number of phases

  OUTPUT : void
 */
{
  int cn, cc, ck; /* counter: node, phase, direction */
  int n_ind; /* Node index */
  real *fg_c, *fg_tmp_c; /* pointer to the first entry of a phase in fg_ and fg_tmp_ */
  
  //struct Node *node = sys->node;
  struct Step *step = sys->step;
  
  for( cc = 0; cc < field->n_c; ++cc) { /* Phase */
    fg_c = &field->f[step->f_phase*cc];
    fg_tmp_c = &field->f_tmp[step->f_phase*cc];
    
    n_ind = 0;
    
    for( cn = 0; cn < sys->max_n_tot; ++cn ) { /* node */
      //#ifdef _MPI_ 
      if( (node->mode[cn]<PERIODIC) && !(node->is_ghost[cn]) ) {
        //#else // NOT _MPI_ 
        //	if(node->mode[cn]<PERIODIC) {
        //#endif // _MPI_ 
        for( ck = 0; ck < sys->n_Q; ++ck ) { /* direction */
          /* Propagate */
          fg_tmp_c[n_ind + ck] = fg_c[n_ind + step->f_rel_neig[ck] + ck];
          // fg_c[n_ind + step->f_rel_neig[ck] + ck] = 0;
        }
      }
      n_ind += step->f_node;
    }
  }
}

void propagate_dif(struct Node *node, struct Field *dif, struct System *sys)
/* 
propagate :
	propagates the distribution 	

  INPUT : 1) fg     : pointer to the distribution
          2) fg_tmp : pointer to the temporal distribution
          3) n_c    : number of phases

  OUTPUT : void
 */
{
  int cn, cc, ck; /* counter: node, phase, direction */
  int n_ind; /* Node index */
  real *fg_c, *fg_tmp_c; /* pointer to the first entry of a phase in fg_ and fg_tmp_ */
  
  //struct Node *node = sys->node;
  struct Step *step = sys->step;
  
  for( cc = 0; cc < dif->n_c; ++cc) { /* Phase */
    fg_c = &dif->f[step->f_phase*cc];
    fg_tmp_c = &dif->f_tmp[step->f_phase*cc];
    
    n_ind = 0;
    
    for( cn = 0; cn < sys->max_n_tot; ++cn ) { /* node */
      //#ifdef _MPI_ 
      if( ((node->mode[cn]<PERIODIC) && !(node->is_ghost[cn])) && (dif->rho[cn] > dif->rho_max) ) {
        //#else // NOT _MPI_ 
        //	if(node->mode[cn]<PERIODIC) {
        //#endif // _MPI_ 
        for( ck = 0; ck < sys->n_Q; ++ck ) { /* direction */
          /* Propagate */
          fg_tmp_c[n_ind + ck] = fg_c[n_ind + step->f_rel_neig[ck] + ck];
          node->Bmask[(n_ind + step->f_rel_neig[ck])/step->f_node] = 1;
          // fg_c[n_ind + step->f_rel_neig[ck] + ck] = 0;
        }
      }
      else{
        for( ck = 0; ck < sys->n_Q; ++ck )
        fg_tmp_c[n_ind + ck] = fg_c[n_ind + ck];
      }
      n_ind += step->f_node;
    }
  }
}