#ifndef SPLAYTREE_H
#define SPLAYTREE_H

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include "chem_global.h"
#include <climits>
#include <cfloat>
#include <algorithm>
#include "defines.h"

#define KEYTYP int
#define VALTYP double  

/* floor function for negative values */
#define NEGFLR(x) ( ((x)<0) ? (  ((int) (x))-1) : ( ((int) (x)) ) )

/* STRUCTURES                                                   */
/* splayTree : main splaytree structure                         */
/*             Contains T, the root and the information header  */
/*             and nil node                                     */
/* treeNode  : Holds the key, value and left and right pointers */

struct treeNode {
  KEYTYP  * key;   /* list of key values used by the search algorithm */
  VALTYP  * val;   /* values for the given key */
  int min_key;     /* mineral combination kode */
  struct treeNode *left, *right;  /* pointer to the left- and right-hand nodes (ST) */
};

struct splayTree {
  /* size      : number of tree nodes 
   * num_while : use to test the efficiency
   */
  long int size, num_while; 
  /* Special splay tree nodes 
   */
  struct treeNode *t, *nil_node, *header, *new_node; 
  /* val_scale_log     : rescale factor for the logarithim binned values
   * val_scale_log_inv : factor for scaling logarithmic binned key values
   * val_scale_lin     : rescale factor for the linear binned values
   * val_scale_lin_inv : factor for scaling linear binned key values 
   * key_to_val        : used to translate 'key values' to physical values
   */
  VALTYP * val_scale_log, * val_scale_log_inv, * val_scale_lin, * val_scale_lin_inv, * key_to_val;
  
  VALTYP * val_log, * val_lin;
  /* a         : logarithmic bin width factor (physcical value)
   * log_a_inv : invers of the logartihm (base 10) of 'a'
   * b         : linear intevall length
   */
  VALTYP a, log_a_inv, b;
  KEYTYP * new_key;
  /* num_val     : number of return values
   * num_val_log : number of logarithmic binned key values
   * num_val_lin : number of lineare binned key values
   * num_key     : total number of key values (num_val_log + num_val_lin)
   */
  int  num_val, num_val_log, num_val_lin, num_key;
  /*** chemistry relevant paramters ***
   * Vchem_key : holds the vchem structure for the given mineral combination
   * min_key   : mineral key. Binary representation of the mineral combination
   */
  struct BasVec * Vchem_key;
  int min_key;
};

/* FUNCTIONS  */

void               SetNewVal(struct treeNode *, struct splayTree *, real * ctot, real * ctot_min, int ctot_min_size);
void               CopyKey(KEYTYP *dest, KEYTYP *orig, struct splayTree * st); 
int                CompKey(KEYTYP * k1, KEYTYP * k2, struct splayTree * st);
void               FreeTreeNode(struct treeNode *tn);
struct splayTree * InitializeSplayTree(int num_val, int num_val_log, int num_val_lin, int num_interval_log, double interval_lin, VALTYP * val_scale, VALTYP * val_scale_lin);
void               Splay(KEYTYP * key, struct splayTree * st);
struct treeNode *  FindInsert(VALTYP * val_log, VALTYP * val_lin, struct splayTree *st, real * ctot_min, int ctot_min_size);
struct treeNode * FindInsertNoInterp(VALTYP * val_log, VALTYP * val_lin, struct splayTree *st, real * ctot_min, int ctot_min_size);
struct treeNode *  MakeNewNode(struct splayTree *);
void               DeleteKey(KEYTYP * key, struct splayTree *st);
void               EmptySplayTree(struct splayTree *st);
void               FreeSplayTree(struct splayTree *st); 
void               SetNewKey(KEYTYP *key, VALTYP * val_log, VALTYP * val_lin, struct splayTree * st);


#endif
