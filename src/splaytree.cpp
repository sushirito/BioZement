#include "splaytree.h"


void SetNewVal(struct treeNode * tn, struct splayTree * st, real * ctot, real * ctot_min, int ctot_min_size)
{
  int n;
  VALTYP * val;

  
  // additional value added to hold pH
  if ( (tn->val = (VALTYP *) malloc((st->num_val+1)*sizeof(VALTYP))) == NULL ) {
    printf("ERROR (SetNewVal): could not allocate space for tn->val\n");
    exit(1);
  }

  /* Set input values values based on the key */
  /* key = tn->key; */
  val = st->key_to_val;
  for ( n=0; n<st->num_val_log; ++n ) 
    *val++ = ctot[n];

  tn->min_key = solve_chem_boundary_node_interpolate(st->key_to_val, ctot_min, st->min_key, st->Vchem_key);  

  if (tn->min_key == st->min_key) {
    /* Copy consentrations */
    for ( n=0; n<st->num_val; ++n) 
      tn->val[n] = st->key_to_val[n];
    // write H+ activity (pH)
    tn->val[st->num_val] = st->Vchem_key->log_a[st->Vchem_key->pos_pH];
  }
}

void SetNewKey(KEYTYP *key, VALTYP * val_log, VALTYP * val_lin, struct splayTree * st)
{
  int n;

  for ( n=0; n<st->num_val_log; ++n )
    *key++ = NEGFLR( st->log_a_inv*log10(st->val_scale_log_inv[n]*std::max(val_log[n],DBL_MIN) ) );
  for ( n=0; n<st->num_val_lin; ++n ) {
    *key++ = NEGFLR(st->val_scale_lin_inv[n]*val_lin[n] + 0.5); 
    /*    printf("ERROR (eje)\n st->num_val_lin = %d ()\n", st->num_val_lin); */
  }
}

void CopyKey(KEYTYP *dest, KEYTYP *orig, struct splayTree * st) 
{
  int n;
  
  for ( n=0; n<st->num_key; ++n )
    dest[n] = orig[n];
}

int CompKey(KEYTYP * k1, KEYTYP * k2, struct splayTree * st)
{
  int n;

  for ( n=0; n<st->num_key; ++n ) {
    if (k1[n]>k2[n])  return 1;
    if (k2[n]>k1[n])  return 2;
  }

  return 0;
}

void FreeTreeNode(struct treeNode *tn)
{
  if (tn) {
  /* FREE KEY */
    if (tn->key)
      free(tn->key);
  /* FREE VAL */    
    if (tn->val)
      free(tn->val);
    free(tn);
  }
}

struct splayTree * InitializeSplayTree(int num_val, int num_val_log, int num_val_lin, int num_interval_log, double interval_lin, VALTYP * val_scale_log, VALTYP * val_scale_lin)
/*  num_val          : number of return values 
 *  num_val_log      : number of key values with logarithmic binning
 *  num_val_lin      : number of key values with linear binning
 *  num_interval_log : number of intervals for the logarithmic binning (base 10)
 *  interval_lin     : length of a intevall for the linear binning
 *  val_scale_log    : rescale factor for the logarithim binned values, if NULL => set to 1.0
 *  val_scale_lin    : rescale factor for the linear binned values, if NULL => set to 1.0
*/
{
  int n; 
  struct splayTree * ret;

  if ( (ret = (struct splayTree *) malloc(sizeof(struct splayTree))) == NULL )
    return NULL;
  ret->t = ret->nil_node = ret->header = NULL;
  ret->size = 0;
  ret->num_val = num_val;
  ret->num_val_log = num_val_log;
  ret->num_val_lin = num_val_lin;
  ret->num_key = num_val_log + num_val_lin;

  /* Special splay tree nodes */
  ret->nil_node = MakeNewNode(ret);
  ret->nil_node->left = ret->nil_node->right = ret->nil_node;
  ret->header = MakeNewNode(ret);
  ret->t = ret->nil_node; 
  
  /* Set key related constants */
  ret->val_log = NULL;
  if (ret->num_val_log)
    if ( (ret->val_log = (VALTYP *) malloc(ret->num_val_log*sizeof(VALTYP))) == NULL )
      return NULL;  
  /* Initialize value */
  for (n=0; n<ret->num_val_log; ++n)
    ret->val_log[n] = 1;


  ret->val_lin = NULL;
  if (ret->num_val_lin)
    if ( (ret->val_lin = (VALTYP *) malloc(ret->num_val_lin*sizeof(VALTYP))) == NULL )
      return NULL;  
  /* Initialize value */
  for (n=0; n<ret->num_val_lin; ++n)
  ret->val_lin[n] = 1; 

  if ( (ret->new_key = (KEYTYP *) malloc(ret->num_key*sizeof(KEYTYP))) == NULL )
    return NULL;  

  if ( (ret->val_scale_log = (VALTYP *) malloc(ret->num_val_log*sizeof(VALTYP))) == NULL )
    return NULL;

  if ( (ret->val_scale_lin = (VALTYP *) malloc(ret->num_val_lin*sizeof(VALTYP))) == NULL )
    return NULL;


  if ( (ret->key_to_val = (VALTYP *) malloc(ret->num_key*sizeof(VALTYP))) == NULL )
    return NULL;
  
  if ( val_scale_log != NULL )
    for ( n=0; n<ret->num_val_log; ++n )
      ret->val_scale_log[n] = val_scale_log[n];
  else
    for ( n=0; n<ret->num_val_log; ++n )
      ret->val_scale_log[n] = 1.0;

  if ( val_scale_lin != NULL )
    for ( n=0; n<num_val_lin; ++n )
      ret->val_scale_lin[n] = val_scale_lin[n];
  else
    for ( n=0; n<num_val_lin; ++n )
      ret->val_scale_lin[n] = 1.0;

  if ( (ret->val_scale_log_inv = (VALTYP *) malloc(ret->num_val_log*sizeof(VALTYP))) == NULL )
    return NULL;

  if ( (ret->val_scale_lin_inv = (VALTYP *) malloc(ret->num_val_lin*sizeof(VALTYP))) == NULL )
    return NULL;


  /* -- calculate the log base */
  ret->a = pow10(1./num_interval_log);
  ret->log_a_inv = 1.*num_interval_log;
  /* -- -- set scaling values to find key */
  for ( n=0; n<ret->num_val_log; ++n )
    ret->val_scale_log_inv[n] = 2.0*ret->a/(ret->val_scale_log[n]*(ret->a + 1.0));

  /* -- calculate the lin base */
  ret->b = interval_lin;
  for ( n=0; n<ret->num_val_lin; ++n ) {
    ret->val_scale_lin[n] *= ret->b;
    ret->val_scale_lin_inv[n] = 1./ret->val_scale_lin[n];
  }

  /* Make new node */
  ret->new_node = MakeNewNode(ret);
  ret->size++;
  
  /*** Variable used for chemical sollutions ***/
  ret->Vchem_key = NULL;
  ret->min_key   = -1;

  return ret;
}

void Splay(KEYTYP * key, struct splayTree * st) 
{
  struct treeNode *left_tree_max, *right_tree_min; /* Make static ?*/
  struct treeNode *tmp;

  st->header->left = st->header->right = st->nil_node;
  left_tree_max = right_tree_min = st->header;

  CopyKey(st->nil_node->key, key, st);

  while ( CompKey(key,st->t->key, st) ) { /* while key not found */

    /* EJ: TEST 
       st->num_while++; */

    if ( CompKey(key, st->t->key, st) == 2 ) { 
      if ( CompKey(key, st->t->left->key, st) == 2 ) {
      /* Rotate with left child */
	tmp = st->t;
	st->t = tmp->left;
	tmp->left = st->t->right;
	st->t->right = tmp;
      }
      if (st->t->left == st->nil_node )
	break;
      /* Link right */
      right_tree_min->left = st->t;
      right_tree_min = st->t;
      st->t = st->t->left;
    } else {
      if ( CompKey(key, st->t->right->key, st) == 1 ) {
	/* Rotate with right child */
	tmp = st->t;
	st->t = tmp->right;
	tmp->right = st->t->left;
	st->t->left = tmp;
      }
      if (st->t->right == st->nil_node)
	break;
      /* Link left */
      left_tree_max->right = st->t;
      left_tree_max = st->t;
      st->t = st->t->right;
    }    

  } /* end wile key not found */

  /* Reassmeble */
  left_tree_max->right = st->t->left;
  right_tree_min->left = st->t->right;
  st->t->left = st->header->right;
  st->t->right = st->header->left;



}

struct treeNode * FindInsertNoInterp(VALTYP * val_log, VALTYP * val_lin, struct splayTree *st, real * ctot_min, int ctot_min_size)
{
  int n;
  VALTYP * val;
  struct treeNode * tn;

  /* BEGIN The part from FindInsert */
  SetNewKey(st->new_key, val_log, st->val_lin, st);
  CopyKey(st->new_node->key, st->new_key, st);

  st->t = st->new_node;
  /* END part from FindInsert */

  tn = st->t;
  /* BEGIN Part from SetNewVal */
  if ( (tn->val == NULL) && ((tn->val = (VALTYP *) malloc(st->num_val*sizeof(VALTYP))) == NULL) ) {
    printf("ERROR (SetNewVal): could not allocate space for tn->val\n");
    exit(1);
  }

  /* Set input values values based on the key */
  /* key = tn->key; */
  val = st->key_to_val;
  for ( n=0; n<st->num_val_log; ++n ) 
    *val++ = val_log[n];

  //if (_DEBUG_FLAG_) printf("\nsolve_chem_boundary_node_interpolate\n");
  tn->min_key = solve_chem_boundary_node_interpolate(st->key_to_val, ctot_min, st->min_key, st->Vchem_key);  

  if (tn->min_key == st->min_key) {
    /* Copy consentrations */
    for ( n=0; n<st->num_val; ++n) 
      tn->val[n] = st->key_to_val[n];
  }
  /* END Part from SetNewVal */

  return st->t;
}


struct treeNode * FindInsert(VALTYP * val_log, VALTYP * val_lin, struct splayTree *st, real * ctot_min, int ctot_min_size)
/* val_log : values logarithmic binning
 * val_lin : values linear binning
 * st      : pointer to splay tree
 */
{ 
  SetNewKey(st->new_key, val_log, st->val_lin, st);
  CopyKey(st->new_node->key, st->new_key, st);

  if (st->t == st->nil_node) {  /* If empty tree */   /* Denne kan vi ta bort, ved aa initialisere st->t */
    st->new_node->left = st->new_node->right = st->nil_node;
    st->t = st->new_node;
  } else { /* else not empty tree */
    Splay(st->new_key, st);
    
    if ( CompKey(st->new_key, st->t->key, st) == 2 ) {
      st->new_node->left = st->t->left;
      st->new_node->right = st->t;
      st->t->left = st->nil_node;
      st->t = st->new_node;
    } else if ( CompKey(st->new_key, st->t->key, st) == 1 ) {
      st->new_node->right = st->t->right;
      st->new_node->left = st->t;
      st->t->right = st->nil_node;
      st->t = st->new_node;
    }
  }

  if (st->t == st->new_node) {
    /* Set new val */
    SetNewVal(st->t, st, val_log, ctot_min, ctot_min_size);
    /* Make a new new_node */
    st->new_node = MakeNewNode(st);
    st->size++;
  }
  return st->t;
}

struct treeNode * MakeNewNode(struct splayTree * st)
{
  struct treeNode * ret;

  if ( (ret = (struct treeNode *) malloc(sizeof(struct treeNode))) == NULL )
    return NULL;
  /* Allocate memory to the key */
  if ( (ret->key = (KEYTYP *) malloc(st->num_key*sizeof(KEYTYP))) == NULL )
    return NULL;
 
  ret->val = NULL; /* Set to default (no value calculated) */
  
  return ret;
}


void DeleteKey(KEYTYP * key, struct splayTree *st)
{
  struct treeNode  * tmp;

  if ( st->t != st->nil_node ) {
    Splay(key, st);
    if ( CompKey(st->t->key, key, st) == 0 ) {
      tmp = st->t;
      if ( st->t->left == st->nil_node ) {
        st->t = tmp->right;
      } else {
        st->t = tmp->left;
        Splay(key, st);
        st->t->right = tmp->right;
      }
      FreeTreeNode(tmp);
      st->size--;
    }
  }
}

void EmptySplayTree(struct splayTree *st)
{
  while (st->t != st->nil_node )
    DeleteKey( st->t->key, st );
}

void FreeSplayTree(struct splayTree *st)
{
  if (st) {
    EmptySplayTree(st);

    if (st->new_node)
      FreeTreeNode(st->new_node);
    if (st->header)
      FreeTreeNode(st->header);
    if (st->nil_node)
      FreeTreeNode(st->nil_node);
    if (st->val_scale_log)
      free(st->val_scale_log);
    if (st->val_scale_log_inv)
      free(st->val_scale_log_inv);
    if (st->val_scale_lin)
      free(st->val_scale_lin);
    if (st->val_scale_lin_inv)
      free(st->val_scale_lin_inv);
    if (st->key_to_val)
      free(st->key_to_val);

    if (st->val_log)
      free(st->val_log);
    if (st->val_lin)
      free(st->val_lin);

    if (st->new_key)
      free(st->new_key);

    free(st);
  }
}
