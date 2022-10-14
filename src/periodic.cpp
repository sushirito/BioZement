#include "global.h"

/* ---------------------------------------------------------------------------------  add_periodic_bc */
void add_periodic_bc(struct System *sys)
/* 
add_periodic_bc :
	adds space for a rim of periodic nodes. Updates
	'sys->max_n' and 'sys->max_n_tot'.


  INPUT  : void
  OUTPUT : void 
 */
{
  int cd; /* Counter for dimensions */

  sys->max_n_tot = 1;
  for( cd = 0; cd < sys->n_D; ++cd )
  {
    sys->max_n[cd] += 2;
    sys->max_n_tot *= sys->max_n[cd];
  }
}


/* ---------------------------------------------------------------------------------  setup_periodic_bc */
void setup_periodic_bc(struct Node *node, struct System *sys)
/*
setup_periodic_bc :
	initiates 1) node->nr_per : number of periodic nodes
	          2) node->list_per : list over periodic nodes and which nodes they represent.
					    The structure is [pn_0, sn_0, pn_1, sn_1, ..., pn_n_per_-1, sn_n_per_-1],
						where pn_0 representes sn_0 and so on.

  INPUT  : void

  OUTPUT : void
 */
{
  int cd, ccd, cp, ccp, cccp; /* Counter dimension, dimension(2), periodic, periodic(2), periodic(3) */
  int n_per, n_per_tmp; /* Help for counting number of periodic nodes */

  //struct Node *node = sys->node;

  /* Find number of periodic nodes */
  node->nr_per = 2;

  for( cd = 1; cd < sys->n_D; ++cd )
  {
    node->nr_per *= (sys->max_n[cd] - 2);

    n_per_tmp = 1;
    for( ccd = 0; ccd < cd; ++ccd )
      n_per_tmp *= sys->max_n[ccd];

    node->nr_per += 2*n_per_tmp;
  }

  /* Allocate memory for the periodic node array */
  /* Note per_nodes [periodic node, system node, ...] */
  if ( !( node->list_per = (int *) calloc(2*node->nr_per, sizeof(int)) )  )
  {
    printf("error allocating memory for 'node->list_per'\n");
    exit(1);
  }

  /* Initiate the list of periodic boundary nodes */
  n_per         = 2;
  node->list_per[0] = 0;
  node->list_per[2] = sys->max_n[0]-1;

  for( cd = 1; cd < sys->n_D; ++cd )
  {
    /* Calculate the multiplicative factor */
    n_per_tmp = 1;
    for( ccd = 0; ccd < cd; ++ccd )
      n_per_tmp *= sys->max_n[ccd];

    /* Initiate counter */
    cp = 0;

    /* Add confining "side" boundary */
    for( cccp = 0; cccp < n_per; ++cccp )
    {
      node->list_per[cp] += n_per_tmp;
      cp += 2;
    }

    for( ccp = 1; ccp < sys->max_n[cd] - 2; ++ccp )
    {
      for( cccp = 0; cccp < n_per; ++cccp )
      {
        node->list_per[cp] = node->list_per[cp-2*n_per] + n_per_tmp;
        cp += 2;
      }
    }

    n_per  *= (sys->max_n[cd] - 2);

    /* Add "top" and "bottom" - layer */
    for( ccp = 0; ccp < n_per_tmp; ++ccp )
    {
      node->list_per[cp] = ccp;
      cp            += 2;
      node->list_per[cp] = ccp + n_per_tmp*(sys->max_n[cd] - 1);
      cp            += 2;
    }

    n_per  += 2*n_per_tmp;
  }

  /* Find the system node represented by the periodic nodes */
  // Due to MPI this is only valid for 1 process. Otherwise set to 0 for compatibility.
  if (sys->mpi->nr_procs<2) {
    for( cp = 0; cp < node->nr_per; ++cp ) {
      node->list_per[2*cp + 1] = find_sn(node->list_per[2*cp], sys);
    }
  } else {
    for( cp = 0; cp < node->nr_per; ++cp ) {
      node->list_per[2*cp + 1] = 0;
    }
  }

}


/* ---------------------------------------------------------------------------------  find_sn */
int find_sn(int pn, struct System *sys)
/*
find_sn :
	finds the system node that is represented a periodic node

  INPUT  : 1) pn : index to a periodic node

  OUTPUT : (int) system nodes represented by 'pn'
 */
{
  int cd; /* Counter dimension */
  int mult_f, sn, nd; /* Multiplicative factor, system node, value node direction */

  /* Initiations */
  mult_f  = 1;
  sn = pn;

  for( cd = 0; cd < sys->n_D; ++cd )
  {
    nd = pn % sys->MAX_N[cd];

    if( nd == 0 )
      sn += (sys->MAX_N[cd] - 2)*mult_f;
    else if( nd == sys->MAX_N[cd] - 1 )
      sn -= (sys->MAX_N[cd] - 2)*mult_f;

    pn = (pn - nd)/sys->MAX_N[cd];
    mult_f *= sys->MAX_N[cd];
  }

  return sn;
}



/* ---------------------------------------------------------------------------------  periodic_bc */
void periodic_bc(struct Node *node, struct Field *field, struct System *sys)
/*
periodic_bc :
	Enforces the periodic boundary conditions

  INPUT : 1) fg  : pointer to a vel. dist.
		  2) fg_tmp : pinter to the temp velocity distribution
          2) n_c : number of phases

  OUTPUT : void
 */
{

  /* No need to update the node (yet) since  
     node-modes dont change (yet) during simulation */
  /*   communicate_node(sys); */

  /* printf("<%d> PERIODIC_BC\n",sys->mpi->my_rank);fflush(stdout); */

#ifndef _2D_RUN_  // true for 3D runs
  communicate(field->f    , field->n_c, sys);
  communicate(field->f_tmp, field->n_c, sys); // uncommenting this line leads to problems
                                              // when we have only one _MPI_ process in y- or z-dir
  return;
#else           // true for 2D runs

  int cn, cp, cc, ck; /* Counters periodic nodes, node->list_per, phases, directions */
  int pn, sn; /* periodic and system node */
  int pn_ind, sn_ind; /* index for pn and sn in f */
  real *fg_c;    /* Pointer to the first index in vel. dist. for a phase */
  real *fg_tmp_c;    /* Pointer to the first index in temp vel. dist. for a phase */

  //struct Node *node = sys->node;
  struct Step *step = sys->step;

  for( cc = 0; cc < field->n_c; ++cc )
  {
    fg_c     = &(field->f)[step->f_phase*cc];
    fg_tmp_c = &(field->f_tmp)[step->f_phase*cc];

    cp = 0;

    for( cn = 0; cn < node->nr_per; ++cn )
    {
      pn = node->list_per[cp++]; /* periodic node */
      sn = node->list_per[cp++]; /* system node (represented by the periodic node) */
      node->mode[pn] = PERIODIC + node->mode[sn]; /* update the node mode such that periodic nodes can
  								   be recognized as either SOLID LIQUID WALL*/
      pn_ind = step->f_node*pn;
      sn_ind = step->f_node*sn;

      /* copy information to the periodic node */
      for( ck = 0; ck < sys->n_Q; ++ck)
      {
        fg_c[pn_ind + ck] = fg_c[sn_ind + ck];
        fg_tmp_c[pn_ind + ck] = fg_tmp_c[sn_ind + ck];
      }
    }
  }
#endif
}

