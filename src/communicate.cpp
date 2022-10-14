/*
 * communicate.c
 *
 *  Created on: 7. mars 2011
 *      Author: janlv
 */
#include "global.h"

//---------------------------------
// Return: local node index or -1 if node is not in this process
// Input: global node index
//---------------------------------
int global_to_local_ind(int Nn, struct System *sys)
{
  int *N = sys->MAX_N;

  int Nx = Nn%N[0];
  int Ny = ((Nn-(Nx))/N[0])%N[1];
  int Nz = (int)(Nn/(N[0]*N[1]));

  // add 1 to avoid the local ghost/periodic node
  int nx, ny, nz;
  nx = Nx - sys->mpi->lower_bound[0] + 1;
  if ( (nx<sys->mpi->lower_bound[0]) || (nx>sys->mpi->upper_bound[0]) )
    return(-1);
  ny = Ny - sys->mpi->lower_bound[1] + 1;
  if ( (ny<sys->mpi->lower_bound[1]) || (ny>sys->mpi->upper_bound[1]) )
    return(-1);
  nz = Nz - sys->mpi->lower_bound[2] + 1;
#ifdef _2D_RUN_
  nz = 0;
#else
  if ( (nz<sys->mpi->lower_bound[2]) || (nz>sys->mpi->upper_bound[2]) )
    return(-1);
#endif
  int n_ind = sys->max_n[1]*sys->max_n[0]*nz + sys->max_n[0]*ny + nx;
  return(n_ind);
}

//---------------------------------
// Return: local node index
// Input: global node coordinate
//---------------------------------
int get_local_ind(int Nx,int Ny,int Nz, struct System *sys)
{
  // add 1 to avoid the local ghost/periodic node
  int nx = Nx - sys->mpi->lower_bound[0] + 1;
  int ny = Ny - sys->mpi->lower_bound[1] + 1;
  int nz = 0;
  if (sys->n_D>2)
    nz = Nz - sys->mpi->lower_bound[2] + 1;
  int n_ind = sys->max_n[1]*sys->max_n[0]*nz + sys->max_n[0]*ny + nx;
  return(n_ind);
}

//---------------------------------
// 
// 
//---------------------------------
int global_to_local_xind(int Nx, struct System *sys)
{
  // add 1 to avoid the local ghost/periodic node
  int nx = Nx - sys->mpi->lower_bound[0] + 1;
  if (nx>sys->max_n[0]) {
    nx = -1;
  }
  return(nx);
}

//---------------------------------
// Return: global node index
// Input:  local node coordinate
//---------------------------------
int get_global_ind(int nx,int ny,int nz, struct System *sys)
{
  // add 1 to avoid the local ghost/periodic node
  int Nx = nx + sys->mpi->lower_bound[0] - 1;
  int Ny = ny + sys->mpi->lower_bound[1] - 1;
  int Nz = 0;
  if (sys->n_D>2)
    Nz = nz + sys->mpi->lower_bound[2] - 1;
  int n_ind = sys->MAX_N[1]*sys->MAX_N[0]*Nz + sys->MAX_N[0]*Ny + Nx;
  return(n_ind);
}

//---------------------------------
// Output: global node x-coordinate
// Input: local node index
//---------------------------------
int get_global_xcoord(int nn, struct System *sys)
{
  int nx = nn%sys->max_n[0];
  // subtract 1 to avoid the local ghost/periodic node
  return(nx + sys->mpi->lower_bound[0] - 1);
}

//---------------------------------
// Output: global node y-coordinate
// Input: local node index
//---------------------------------
int get_global_ycoord(int nn, struct System *sys)
{
  int *n = sys->max_n;
  int ny = ((nn-(nn%n[0]))/n[0])%n[1];
  // subtract 1 to avoid the local ghost/periodic node
  return(ny + sys->mpi->lower_bound[1] - 1);
}

//---------------------------------
// Output: global node z-coordinate
// Input: local node index
//---------------------------------
int get_global_zcoord(int nn, struct System *sys)
{
#ifdef _2D_RUN_
  return(0);
#endif
  int *n = sys->max_n;
  int nz = (int)(nn/(n[0]*n[1]));
  // subtract 1 to avoid the local ghost/periodic node
  return(nz + sys->mpi->lower_bound[2] - 1);
}

//---------------------------------
//---------------------------------
void get_local_coord(int nn, int *ind, struct System *sys)
{
  int *n = sys->max_n;
  ind[0] = nn%n[0];
  ind[1] = ((nn-(ind[0]))/n[0])%n[1];
  ind[2] = (int)(nn/(n[0]*n[1]));
#ifdef _2D_RUN_
  ind[2] = 0;
#endif
}

//---------------------------------
//---------------------------------
void get_global_coord(int nn, int *ind, struct System *sys)
{
  get_local_coord(nn, ind, sys);
  ind[0] += sys->mpi->lower_bound[0] - 1;
  ind[1] += sys->mpi->lower_bound[1] - 1;
  ind[2] += sys->mpi->lower_bound[2] - 1;
#ifdef _2D_RUN_
  ind[2] = 0;
#endif
}




int init_communicate(struct System *sys)

/* Defines new MPI_Datatypes, and sets
 * send/receive arrays used in communicate() */

{
  struct Mpi *mpi = sys->mpi;

  
  mpi->S_ncom = 0;
  mpi->S_ncpy = 0;
  mpi->E_ncom = 0;
  mpi->E_ncpy = 0;
  mpi->C_ncom = 0;
  mpi->C_ncpy = 0;

  // receive buffer for _MPI_ to communicate flags between processes, e.g
  // velreinit or converge flags
  mpi->flag_rbuf = (int *) malloc(mpi->nr_procs*sizeof(int));

  // receive buffer for double value
  mpi->double_rbuf = (double *) malloc(mpi->nr_procs*sizeof(double));

#ifdef _2D_RUN_
  return GOOD;
#endif

  int *ind = mpi->ind_rank;
  int *np = mpi->np;
  int *n = sys->max_n;
  //int *N = sys->MAX_N;
  int nQ = sys->n_Q;
  int nD = sys->n_D;
  int i, j;
  
  /* define _MPI_ datatypes used by communicate() */

  //
  //  datatypes for distribution-array (n_Q doubles)
  //
  int X=0, Y=1, Z=2;
  // x-line
  int count       = n[0]-2;    // num blocks to send (don't incl. edges!)
  int blocklength = sys->n_Q;  // num elements in block
  int stride      = sys->n_Q;  // num elements between start of blocks
  MPI_Type_vector(count,blocklength,stride,MPI_DOUBLE,&mpi->line[X]);
  MPI_Type_commit(&mpi->line[X]);

  // y-line
  count       = n[1]-2;
  blocklength = sys->n_Q;
  stride      = n[0]*sys->n_Q;
  MPI_Type_vector(count,blocklength,stride,MPI_DOUBLE,&mpi->line[Y]);
  MPI_Type_commit(&mpi->line[Y]);

  // z-line
  count       = n[2]-2;
  blocklength = sys->n_Q;
  stride      = n[0]*n[1]*sys->n_Q;
  MPI_Type_vector(count,blocklength,stride,MPI_DOUBLE,&mpi->line[Z]);
  MPI_Type_commit(&mpi->line[Z]);

  int LR=0, TB=1, FB=2;

  // right/left slice
  count       = n[2]-2;
  blocklength = 1;
  MPI_Aint byte_stride = n[0]*n[1]*sys->n_Q*sizeof(double);
  MPI_Type_create_hvector(count,blocklength,byte_stride,mpi->line[Y],&mpi->slice[LR]);
  MPI_Type_commit(&mpi->slice[LR]);

  // top/bottom slice
  count       = n[2]-2;
  blocklength = (n[0]-2)*sys->n_Q;
  stride      = n[0]*n[1]*sys->n_Q;
  MPI_Type_vector(count,blocklength,stride,MPI_DOUBLE,&mpi->slice[TB]);
  MPI_Type_commit(&mpi->slice[TB]);

  // front/back slice
  count       = n[1]-2;
  blocklength = (n[0]-2)*sys->n_Q;
  stride      = n[0]*sys->n_Q;
  MPI_Type_vector(count,blocklength,stride,MPI_DOUBLE,&mpi->slice[FB]);
  MPI_Type_commit(&mpi->slice[FB]);

  //
  //  datatypes for char-array (e.g. node mode array)
  //
  // x-line
  count       = n[0]-2; // num blocks to send (don't incl. edges!)
  blocklength = 1;      // num elements in block
  stride      = 1;      // num elements between start of blocks
  MPI_Type_vector(count,blocklength,stride,MPI_CHAR,&mpi->line_char[X]);
  MPI_Type_commit(&mpi->line_char[X]);

  // y-line
  count       = n[1]-2;
  blocklength = 1;
  stride      = n[0];
  MPI_Type_vector(count,blocklength,stride,MPI_CHAR,&mpi->line_char[Y]);
  MPI_Type_commit(&mpi->line_char[Y]);

  // z-line
  count       = n[2]-2;
  blocklength = 1;
  stride      = n[0]*n[1];
  MPI_Type_vector(count,blocklength,stride,MPI_CHAR,&mpi->line_char[Z]);
  MPI_Type_commit(&mpi->line_char[Z]);

  // right/left slice
  count       = n[2]-2;
  blocklength = 1;
  byte_stride = n[0]*n[1]*sizeof(char);
  MPI_Type_create_hvector(count,blocklength,byte_stride,mpi->line_char[Y],&mpi->slice_char[LR]);
  MPI_Type_commit(&mpi->slice_char[LR]);

  // top/bottom slice
  count       = n[2]-2;
  blocklength = (n[0]-2);
  stride      = n[0]*n[1];
  MPI_Type_vector(count,blocklength,stride,MPI_CHAR,&mpi->slice_char[TB]);
  MPI_Type_commit(&mpi->slice_char[TB]);

  // front/back slice
  count       = n[1]-2;
  blocklength = (n[0]-2);
  stride      = n[0];
  MPI_Type_vector(count,blocklength,stride,MPI_CHAR,&mpi->slice_char[FB]);
  MPI_Type_commit(&mpi->slice_char[FB]);

  //
  //  datatypes for double-array (e.g. sfrac-array)
  //
  // x-line
  count       = n[0]-2; // num blocks to send (don't incl. edges!)
  blocklength = 1;      // num elements in block
  stride      = 1;      // num elements between start of blocks
  MPI_Type_vector(count,blocklength,stride,MPI_DOUBLE,&mpi->line_double[X]);
  MPI_Type_commit(&mpi->line_double[X]);

  // y-line
  count       = n[1]-2;
  blocklength = 1;
  stride      = n[0];
  MPI_Type_vector(count,blocklength,stride,MPI_DOUBLE,&mpi->line_double[Y]);
  MPI_Type_commit(&mpi->line_double[Y]);

  // z-line
  count       = n[2]-2;
  blocklength = 1;
  stride      = n[0]*n[1];
  MPI_Type_vector(count,blocklength,stride,MPI_DOUBLE,&mpi->line_double[Z]);
  MPI_Type_commit(&mpi->line_double[Z]);

  // right/left slice
  count       = n[2]-2;
  blocklength = 1;
  byte_stride = n[0]*n[1]*sizeof(double);
  MPI_Type_create_hvector(count,blocklength,byte_stride,mpi->line_double[Y],&mpi->slice_double[LR]);
  MPI_Type_commit(&mpi->slice_double[LR]);

  // top/bottom slice
  count       = n[2]-2;
  blocklength = (n[0]-2);
  stride      = n[0]*n[1];
  MPI_Type_vector(count,blocklength,stride,MPI_DOUBLE,&mpi->slice_double[TB]);
  MPI_Type_commit(&mpi->slice_double[TB]);

  // front/back slice
  count       = n[1]-2;
  blocklength = (n[0]-2);
  stride      = n[0];
  MPI_Type_vector(count,blocklength,stride,MPI_DOUBLE,&mpi->slice_double[FB]);
  MPI_Type_commit(&mpi->slice_double[FB]);


//  //
//  //  datatypes for velocity-array (n_D doubles); only L and R
//  //  (used for effluent output)
//  //
//
//  // y-line
//  count       = n[1]-2;
//  blocklength = sys->n_D;
//  stride      = n[0]*sys->n_D;
//  MPI_Type_vector(count,blocklength,stride,MPI_DOUBLE,&mpi->line_vel[Y]);
//  MPI_Type_commit(&mpi->line_vel[Y]);
//
//  // right/left slice
//  count       = n[2]-2;
//  blocklength = 1;
//  byte_stride = n[0]*n[1]*sys->n_D*sizeof(double);
//  MPI_Type_create_hvector(count,blocklength,byte_stride,mpi->line_vel[Y],&mpi->slice_vel[LR]);
//  MPI_Type_commit(&mpi->slice_vel[LR]);

  /* Define send and receive arrays */

  // node stride in x,y,z dir
  int s[3] = {1,n[0],n[0]*n[1]};
  // outer edge lengths (x,y,z), corner stride
  int ol[3] = {n[0]-1, n[0]*(n[1]-1), n[0]*n[1]*(n[2]-1)};
  // inner edge lengths (x,y,z)
  int il[3] = {n[0]-3, n[0]*(n[1]-3), n[0]*n[1]*(n[2]-3)};
  // inner origo
  int io = n[0]*(n[1]+1)+1;

  // process stride in x,y,z dir
  int ps[3] = {1,np[0],np[0]*np[1]};
  // process outer edge lengths (x,y,z), process corner stride
  int pol[3] = {np[0]-1, np[0]*(np[1]-1), np[0]*np[1]*(np[2]-1)};

  // used to determine process locations
  int I[2][3],I_bc[2][3];
  for (i = 0; i < sys->n_D; ++i) {
    I[0][i] = (ind[i] < np[i]-1);
    I[1][i] = (ind[i] > 0      );
    I_bc[0][i] = (ind[i] == np[i]-1);
    I_bc[1][i] = (ind[i] == 0      );
  }

  int A[3]={1,0,0}, a;
  int B[3]={2,2,1}, b;


  /* SIDES */
  // side process stride
  int tmp_side_str[3] = {1, np[0], np[0]*np[1]};
  int tmp_side_str_bc[3] = {np[0]-1, np[0]*(np[1]-1), np[0]*np[1]*(np[2]-1)};

  // a[3][2] = { {L,R}, {BT,T}, {F,BK} };
  /* periodic boundary NOT included */
  /* original version */
  int tmp_side_rcv[3][2] =
  {
      { s[1]+s[2], ol[0]+s[1]+s[2] },
      { s[0]+s[2], ol[1]+s[0]+s[2] },
      { s[0]+s[1], ol[2]+s[0]+s[1] }
  };

  /* periodic boundary IS included */
  /* 	int tmp_side_rcv[3][2] = */
  /* 	  { */
  /* 	    { 0, ol[0] }, */
  /* 	    { 0, ol[1] }, */
  /* 	    { 0, ol[2] } */
  /* 	  }; */

  for (i = 0; i < 3; ++i) {
    mpi->side_str[i] = tmp_side_str[i];
    mpi->side_str_bc[i] = tmp_side_str_bc[i];

    mpi->side_rcv[i][0] = nQ*tmp_side_rcv[i][0];
    mpi->side_rcv[i][1] = nQ*tmp_side_rcv[i][1];
    mpi->side_snd[i][0] = nQ*(tmp_side_rcv[i][0]+s[i]);
    mpi->side_snd[i][1] = nQ*(tmp_side_rcv[i][1]-s[i]);

    mpi->side_rcv_scalar[i][0] = tmp_side_rcv[i][0];
    mpi->side_rcv_scalar[i][1] = tmp_side_rcv[i][1];
    mpi->side_snd_scalar[i][0] = tmp_side_rcv[i][0]+s[i];
    mpi->side_snd_scalar[i][1] = tmp_side_rcv[i][1]-s[i];

    //mpi->side_snd_vel[i][0] = nD*(tmp_side_rcv[i][0]+s[i]);
    //mpi->side_snd_vel[i][1] = nD*(tmp_side_rcv[i][1]-s[i]-sys->opt->eff_offset);
  }

  //   left: side_[0][0]==1
  //  right: side_[0][1]==1
  // bottom: side_[1][0]==1
  // etc...
  //I_bc[0][i] = (ind[i] == np[i]-1);
  //I_bc[1][i] = (ind[i] == 0      );
  for (i = 0; i < sys->n_D; ++i) {
    mpi->side[i][0] = I[1][i];
    mpi->side[i][1] = I[0][i];
    mpi->side_bc[i][0] = I_bc[1][i]; //= (ind[i] == 0      );
    mpi->side_bc[i][1] = I_bc[0][i]; //= (ind[i] == np[i]-1);
  }


  /* EDGES */

  // process stride for transferring edges between processes
  // first index loops over lines along x-, y-, and z-axis
  int tmp_edge_str[3][2] =
  {
      { ps[1]+ps[2], -ps[1]+ps[2] },
      { ps[0]+ps[2], -ps[0]+ps[2] },
      { ps[0]+ps[1], -ps[0]+ps[1] }
  };
  /* 	{ */
  /* 			{ np[0]*(np[1]+1), np[0]         }, */
  /* 			{ np[0]*np[1]+1  , np[0]*np[1]-1 }, */
  /* 			{ np[0]+1        , np[0]-1       } */
  /* 	}; */
  /* 	int tmp_edge_str_bc[3][2] = */
  /* 	  { */
  /* 	    { np[0]*(np[1]*np[2]-1)      , np[0]*(np[1]*(np[2]-2)+1) }, */
  /* 	    { np[0]*np[1]*(np[2]-1)+np[0]-1, np[0]*np[1]*(np[2]-1)-np[0]+1 }, */
  /* 	    { np[0]+1              , np[0]-1       } */
  /* 	  }; */
  int tmp_edge_str_bc[3][2] =
  {
      { pol[1]+pol[2], -pol[1]+pol[2] },
      { pol[0]+pol[2], -pol[0]+pol[2] },
      { pol[0]+pol[1], -pol[0]+pol[1] }
  };
  int tmp_edge_rcv[3][4] =
  {
      { ol[1]+ol[2], 0, ol[2], ol[1] },
      { ol[0]+ol[2], 0, ol[2], ol[0] },
      { ol[0]+ol[1], 0, ol[1], ol[0] }
  };
  int tmp_edge_snd[3][4] =
  {
      { il[1]+il[2], 0, il[2], il[1] },
      { il[0]+il[2], 0, il[2], il[0] },
      { il[0]+il[1], 0, il[1], il[0] }
  };
  for (i = 0; i < 3; ++i) {
    for (j = 0; j < 2; ++j) {
      mpi->edge_str[i][j] = tmp_edge_str[i][j];
      mpi->edge_str_bc[i][j] = tmp_edge_str_bc[i][j];
    }
    for (j = 0; j < 4; ++j) {
      mpi->edge_snd[i][j] = nQ*(io + tmp_edge_snd[i][j]);
      mpi->edge_rcv[i][j] = nQ*(s[i] + tmp_edge_rcv[i][j]);
      mpi->edge_snd_scalar[i][j] = io + tmp_edge_snd[i][j];
      mpi->edge_rcv_scalar[i][j] = s[i] + tmp_edge_rcv[i][j];
    }
  }
  // i=0: ignore x-index, i=1: ignore y-index, etc.
  for (i = 0; i < sys->n_D; ++i) {
    a = A[i]; // {1,0,0}
    b = B[i]; // {2,2,1}
    mpi->edge[i][0] = ( (ind[a] < np[a]-1) && (ind[b] < np[b]-1) );
    mpi->edge[i][1] = ( (ind[a] > 0      ) && (ind[b] > 0      ) );
    mpi->edge[i][2] = ( (ind[a] > 0      ) && (ind[b] < np[b]-1) );
    mpi->edge[i][3] = ( (ind[a] < np[a]-1) && (ind[b] > 0      ) );

    mpi->edge_bc[i][0] = ( (ind[a] == np[a]-1) && (ind[b] == np[b]-1) );
    mpi->edge_bc[i][1] = ( (ind[a] == 0      ) && (ind[b] == 0      ) );
    mpi->edge_bc[i][2] = ( (ind[a] == 0      ) && (ind[b] == np[b]-1) );
    mpi->edge_bc[i][3] = ( (ind[a] == np[a]-1) && (ind[b] == 0      ) );
  }


  /* CORNERS */
  // corner process stride

  int tmp_crnr_str[4] =
  {
      ps[2]+ps[1]+ps[0],
      ps[2]+ps[1]-ps[0],
      ps[2]-ps[1]+ps[0],
      ps[2]-ps[1]-ps[0]
  };
  int tmp_crnr_str_bc[4] =
  {
      pol[2]+pol[1]+pol[0],
      pol[2]+pol[1]-pol[0],
      pol[2]-pol[1]+pol[0],
      pol[2]-pol[1]-pol[0]
  };
  int tmp_crnr_snd[4][2] =
  {
      { il[0]+il[1]+il[2], 0           },
      { il[1]+il[2]      , il[0]       },
      { il[0]+il[2]      , il[1]       },
      { il[2]            , il[0]+il[1] }
  };
  int tmp_crnr_rcv[4][2] =
  {
      { ol[0]+ol[1]+ol[2], 0           },
      { ol[1]+ol[2]      , ol[0]       },
      { ol[0]+ol[2]      , ol[1]       },
      { ol[2]            , ol[0]+ol[1] }
  };
  for (i = 0; i < 4; ++i) {
    mpi->crnr_str[i] = tmp_crnr_str[i];
    mpi->crnr_str_bc[i] = tmp_crnr_str_bc[i];
    for (j = 0; j < 2; ++j) {
      mpi->crnr_snd[i][j] = nQ*(io+tmp_crnr_snd[i][j]);
      mpi->crnr_rcv[i][j] = nQ*tmp_crnr_rcv[i][j];
      mpi->crnr_snd_scalar[i][j] = io+tmp_crnr_snd[i][j];
      mpi->crnr_rcv_scalar[i][j] = tmp_crnr_rcv[i][j];
    }
  }

  mpi->crnr[0][0] = ( I[0][0] && I[0][1] && I[0][2] );
  mpi->crnr[0][1] = ( I[1][0] && I[1][1] && I[1][2] );

  mpi->crnr[1][0] = ( I[1][0] && I[0][1] && I[0][2] );
  mpi->crnr[1][1] = ( I[0][0] && I[1][1] && I[1][2] );

  mpi->crnr[2][0] = ( I[0][0] && I[1][1] && I[0][2] );
  mpi->crnr[2][1] = ( I[1][0] && I[0][1] && I[1][2] );

  mpi->crnr[3][0] = ( I[1][0] && I[1][1] && I[0][2] );
  mpi->crnr[3][1] = ( I[0][0] && I[0][1] && I[1][2] );


  mpi->crnr_bc[0][0] = ( I_bc[0][0] && I_bc[0][1] && I_bc[0][2] );
  mpi->crnr_bc[0][1] = ( I_bc[1][0] && I_bc[1][1] && I_bc[1][2] );

  mpi->crnr_bc[1][0] = ( I_bc[1][0] && I_bc[0][1] && I_bc[0][2] );
  mpi->crnr_bc[1][1] = ( I_bc[0][0] && I_bc[1][1] && I_bc[1][2] );

  mpi->crnr_bc[2][0] = ( I_bc[0][0] && I_bc[1][1] && I_bc[0][2] );
  mpi->crnr_bc[2][1] = ( I_bc[1][0] && I_bc[0][1] && I_bc[1][2] );

  mpi->crnr_bc[3][0] = ( I_bc[1][0] && I_bc[1][1] && I_bc[0][2] );
  mpi->crnr_bc[3][1] = ( I_bc[0][0] && I_bc[0][1] && I_bc[1][2] );


  /* Define process-arrays to ease communication */

  int p1,p2,id;
  int rank = mpi->my_rank;


  /*      */
  /* SIDE */
  /*      */
  int nsides = 3;
  int dir_side[3][3] = {{1,0,0},{0,1,0},{0,0,1}};
  int k,p,sign[2] = {-1,1};

  // set up process communication array
  // mpi->side_proc[][]
  for (i = 0; i < nsides; ++i) {

    // periodic communication
    for (k = 0; k < 2; ++k) {
      p = 0;
      for (j = 0; j < sys->n_D; ++j) {
        id = ind[j]+sign[k]*dir_side[i][j];
        p += ( (id+np[j])%np[j] )*mpi->side_str[j];
      }
      mpi->side_proc[k][i] = p;
    }

    // NON-periodic communication
    for (k = 0; k < 2; ++k) {
      p = 0;
      for (j = 0; j < sys->n_D; ++j) {
        id = ind[j]+sign[k]*dir_side[i][j];
        if (id>=0 && id<=np[j]-1 && p>=0)
          p += id*mpi->side_str[j];
        else
          p = -1;
      }
      mpi->side_proc_nobc[k][i] = p;
    }

    // in case of only one process in
    // direction i we have to copy instead
    // of using _MPI_ 
    p1 = mpi->side_proc[0][i];
    p2 = mpi->side_proc[1][i];
    if ( (p1!=rank) && (p2!=rank) ) {
      mpi->S_com[mpi->S_ncom++] = i;
    } else {
      mpi->S_cpy[mpi->S_ncpy++] = i;
    }
  }

  // odd procs sends first, while
  // even procs receives first (or vice versa...)
  for (i = 0; i < nsides; ++i) {
    mpi->S[i] = (ind[i]%2==0);
  }

  // test
  // for (i = 0; i < nsides; ++i) {
  //   printf("SIDE: <%02d> i: %d, %02d<- nobc ->%02d, %02d<- pbc ->%02d\n",mpi->my_rank,i, 
  // 	   mpi->side_proc_nobc[0][i],mpi->side_proc_nobc[1][i],mpi->side_proc[0][i],mpi->side_proc[1][i]); 
  // }
  // for (int ii = 0; ii < mpi->S_ncom; ++ii) {
  //   i = mpi->S_com[ii];
  //   if ( mpi->S[i] ) {
  //     printf("%d -> %d",rank, mpi->side_proc[1][i]);
  //     if (mpi->side_proc_nobc[1][i]>=0)
  // 	printf(", %d -> %d",rank, mpi->side_proc_nobc[1][i]);
  //     printf("\n");
      
  //     printf("%d <- %d",rank, mpi->side_proc[1][i]);
  //     if (mpi->side_proc_nobc[1][i]>=0)
  // 	printf(", %d <- %d",rank, mpi->side_proc_nobc[1][i]);
  //     printf("\n");

  //     printf("%d <- %d",rank, mpi->side_proc[0][i]);
  //     if (mpi->side_proc_nobc[0][i]>=0)
  // 	printf(", %d <- %d",rank, mpi->side_proc_nobc[0][i]);
  //     printf("\n");

  //     printf("%d -> %d",rank, mpi->side_proc[0][i]);
  //     if (mpi->side_proc_nobc[0][i]>=0)
  // 	printf(", %d -> %d",rank, mpi->side_proc_nobc[0][i]);
  //     printf("\n");
  //   } else {
  //     printf("%d <- %d",rank, mpi->side_proc[0][i]);
  //     if (mpi->side_proc_nobc[0][i]>=0)
  // 	printf(", %d <- %d",rank, mpi->side_proc_nobc[0][i]);
  //     printf("\n");

  //     printf("%d -> %d",rank, mpi->side_proc[0][i]);
  //     if (mpi->side_proc_nobc[0][i]>=0)
  // 	printf(", %d -> %d",rank, mpi->side_proc_nobc[0][i]);
  //     printf("\n");

  //     printf("%d -> %d",rank, mpi->side_proc[1][i]);
  //     if (mpi->side_proc_nobc[1][i]>=0)
  // 	printf(", %d -> %d",rank, mpi->side_proc_nobc[1][i]);
  //     printf("\n");

  //     printf("%d <- %d",rank, mpi->side_proc[1][i]);
  //     if (mpi->side_proc_nobc[1][i]>=0)
  // 	printf(", %d <- %d",rank, mpi->side_proc_nobc[1][i]);
  //     printf("\n");
  //   }
  // }
  
  /*      */
  /* EDGE */
  /*      */
  int dir_edge[6][3] = {{0,1,1},{0,-1,1},
			{1,0,1},{-1,0,1},
			{1,1,0},{-1,1,0}};
  // set up process communication array
  for (i = 0; i < 6; ++i) {

    // periodic communication
    for (k = 0; k < 2; ++k) {
      p = 0;
      for (j = 0; j < sys->n_D; ++j) {
        id = ind[j]+sign[k]*dir_edge[i][j];
        p += ( (id+np[j])%np[j] )*mpi->side_str[j];
      }
      mpi->edge_proc[k][i] = p;
    }

    // NON-periodic communication
    for (k = 0; k < 2; ++k) {
      p = 0;
      for (j = 0; j < sys->n_D; ++j) {
        id = ind[j]+sign[k]*dir_edge[i][j];
        if (id>=0 && id<=np[j]-1 && p>=0) {
          p += id*mpi->side_str[j];
        } else {
          p = -1;
        }
      }
      mpi->edge_proc_nobc[k][i] = p;
    }

    // in case of only one process in
    // direction i we have to copy instead
    // of using MPI
    p1 = mpi->edge_proc[0][i];
    p2 = mpi->edge_proc[1][i];
    if (i%2==0) {
      if ( (p1!=rank) && (p2!=rank) )
        mpi->E_com[mpi->E_ncom++] = (int)i*0.5;
      else
        mpi->E_cpy[mpi->E_ncpy++] = (int)i*0.5;
    }
  }

  // odd procs sends first, while
  // even procs receives first (or vice versa...)
  mpi->E[0] = (ind[2]%2==0);
  mpi->E[1] = (ind[0]%2==0);
  mpi->E[2] = (ind[1]%2==0);

  if (np[0]==1) {
    mpi->E[1] = (ind[2]%2==0);
  }
  if (np[1]==1) {
    mpi->E[2] = (ind[0]%2==0);
  }
  if (np[2]==1) {
    mpi->E[0] = (ind[1]%2==0);
  }

  // test
  // for (i = 0; i < 6; ++i) {
  //   printf("EDGE: <%02d> i: %d, %02d<- nobc ->%02d, %02d<- pbc ->%02d\n",mpi->my_rank,i, 
  // 	   mpi->edge_proc_nobc[0][i],mpi->edge_proc_nobc[1][i],mpi->edge_proc[0][i],mpi->edge_proc[1][i]); 
  // }


  /*        */
  /* CORNER */
  /*        */

  int dir_crnr[4][3] = {{1,1,1},{-1,1,1},{1,-1,1},{-1,-1,1}};
  for (i = 0; i < 4; ++i) {

    // periodic communication
    for (k = 0; k < 2; ++k) {
      p = 0;
      for (j = 0; j < sys->n_D; ++j) {
        id = ind[j]+sign[k]*dir_crnr[i][j];
        p += ( (id+np[j])%np[j] )*mpi->side_str[j];
      }
      mpi->crnr_proc[k][i] = p;
    }

    // NON-periodic communication
    for (k = 0; k < 2; ++k) {
      p = 0;
      for (j = 0; j < sys->n_D; ++j) {
        id = ind[j]+sign[k]*dir_crnr[i][j];
        if (id>=0 && id<=np[j]-1 && p>=0)
          p += id*mpi->side_str[j];
        else
          p = -1;
      }
      mpi->crnr_proc_nobc[k][i] = p;
    }
    
    p1 = mpi->crnr_proc[0][i];
    p2 = mpi->crnr_proc[1][i];
    if ( (p1!=rank) && (p2!=rank) ) {
      mpi->C_com[mpi->C_ncom++] = i;
    } else {
      mpi->C_cpy[mpi->C_ncpy++] = i;
    }
  }
  mpi->C = (ind[0]%2==0);
  for (i = 0; i < 2; ++i) {
    if (np[i]==1) {
      mpi->C = (ind[i+1]%2==0);
    }
  }

  for (i = 0; i < 3; ++i) {
    for (j = 1; j < 3; ++j) {
      k = (i+j)%3;
      if ( (np[i]==1) && (np[k]==1) ) {
        mpi->C = (ind[(j+k)%3]%2==0);
      }
    }
  }

  // test
  // for (i = 0; i < 4; ++i) {
  //   printf("CRNR: <%02d> i: %d, %02d<- nobc ->%02d, %02d<- pbc ->%02d\n",mpi->my_rank,i, 
  // 	   mpi->crnr_proc_nobc[0][i],mpi->crnr_proc_nobc[1][i],mpi->crnr_proc[0][i],mpi->crnr_proc[1][i]); 
  // }


  /* setup copy-arrays */

  /* SIDES */
  int str[3] = {1,n[0],n[0]*n[1]};
  mpi->cpy_slice_size[LR] = (n[1]-2)*(n[2]-2);
  mpi->cpy_slice_size[TB] = (n[0]-2)*(n[2]-2);
  mpi->cpy_slice_size[FB] = (n[0]-2)*(n[1]-2);

  if ( !( mpi->cpy_slice_LR      = (int *) calloc(mpi->cpy_slice_size[LR], sizeof(int)) )  ) return BAD;
  if ( !( mpi->cpy_slice_TB      = (int *) calloc(mpi->cpy_slice_size[TB], sizeof(int)) )  ) return BAD;
  if ( !( mpi->cpy_slice_FB      = (int *) calloc(mpi->cpy_slice_size[FB], sizeof(int)) )  ) return BAD;
  if ( !( mpi->cpy_slice_char_LR = (int *) calloc(mpi->cpy_slice_size[LR], sizeof(int)) )  ) return BAD;
  if ( !( mpi->cpy_slice_char_TB = (int *) calloc(mpi->cpy_slice_size[TB], sizeof(int)) )  ) return BAD;
  if ( !( mpi->cpy_slice_char_FB = (int *) calloc(mpi->cpy_slice_size[FB], sizeof(int)) )  ) return BAD;

  mpi->cpy_slice[LR] = &mpi->cpy_slice_LR[0];
  mpi->cpy_slice[TB] = &mpi->cpy_slice_TB[0];
  mpi->cpy_slice[FB] = &mpi->cpy_slice_FB[0];
  mpi->cpy_slice_char[LR] = &mpi->cpy_slice_char_LR[0];
  mpi->cpy_slice_char[TB] = &mpi->cpy_slice_char_TB[0];
  mpi->cpy_slice_char[FB] = &mpi->cpy_slice_char_FB[0];

  int x,y,z,d[3];
  int c[3]={0,0,0};
  for (i=0; i<3; ++i) {
    x=(i+0)%3;
    y=(i+1)%3;
    z=(i+2)%3;
    d[x] = 0;
    for (d[y] = 0; d[y] < n[y]-2; ++d[y]) {
      for (d[z] = 0; d[z] < n[z]-2; ++d[z]) {
        mpi->cpy_slice_char[i][c[i]] = d[x]*str[x] + d[y]*str[y] + d[z]*str[z];
        mpi->cpy_slice     [i][c[i]] = nQ*mpi->cpy_slice_char[i][c[i]];
        c[i]++;
      }
    }
  }

  for (i=0; i<3; ++i) {
    if (c[i] != mpi->cpy_slice_size[i]) {
      printf("<%02d> ERROR! Copy-slices have the wrong size!: %d != %d\n",
          mpi->my_rank, c[i], mpi->cpy_slice_size[i]);
      exit(1);
    }
  }

  /* 	if (mpi->my_rank==0) */
  /* 	  for (i=0; i<3; ++i) { */
  /* 	    for (j=0; j<mpi->cpy_slice_size[i]; ++j) { */
  /* 	      printf("[%d]",mpi->cpy_slice_node[i][j]); */
  /* 	    } */
  /* 	    printf("\n\n"); */
  /* 	  } */

  /* EDGES */
  mpi->cpy_line_size[X] = n[0]-2;
  mpi->cpy_line_size[Y] = n[1]-2;
  mpi->cpy_line_size[Z] = n[2]-2;

  if ( !( mpi->cpy_line_X      = (int *) calloc(mpi->cpy_line_size[X], sizeof(int)) )  ) return BAD;
  if ( !( mpi->cpy_line_Y      = (int *) calloc(mpi->cpy_line_size[Y], sizeof(int)) )  ) return BAD;
  if ( !( mpi->cpy_line_Z      = (int *) calloc(mpi->cpy_line_size[Z], sizeof(int)) )  ) return BAD;
  if ( !( mpi->cpy_line_char_X = (int *) calloc(mpi->cpy_line_size[X], sizeof(int)) )  ) return BAD;
  if ( !( mpi->cpy_line_char_Y = (int *) calloc(mpi->cpy_line_size[Y], sizeof(int)) )  ) return BAD;
  if ( !( mpi->cpy_line_char_Z = (int *) calloc(mpi->cpy_line_size[Z], sizeof(int)) )  ) return BAD;

  mpi->cpy_line[X] = &mpi->cpy_line_X[0];
  mpi->cpy_line[Y] = &mpi->cpy_line_Y[0];
  mpi->cpy_line[Z] = &mpi->cpy_line_Z[0];
  mpi->cpy_line_char[X] = &mpi->cpy_line_char_X[0];
  mpi->cpy_line_char[Y] = &mpi->cpy_line_char_Y[0];
  mpi->cpy_line_char[Z] = &mpi->cpy_line_char_Z[0];

  c[0]=0; c[1]=0; c[2]=0;
  for (i=0; i<3; ++i) {
    x=(i+2)%3;
    y=(i+1)%3;
    z=(i+0)%3;
    d[x] = 0;
    d[y] = 0;
    for (d[z] = 0; d[z] < n[z]-2; ++d[z]) {
      mpi->cpy_line_char[i][c[i]] = d[x]*str[x] + d[y]*str[y] + d[z]*str[z];
      mpi->cpy_line     [i][c[i]] = nQ*mpi->cpy_line_char[i][c[i]];
      c[i]++;
    }
  }
  for (i=0; i<3; ++i) {
    if (c[i] != mpi->cpy_line_size[i]) {
      printf("<%02d> ERROR! Copy-lines have the wrong size!: %d != %d\n",
          mpi->my_rank, c[i], mpi->cpy_line_size[i]);
      exit(1);
    }
  }

#ifdef _DEBUG_MPI_
  if (rank==0) {
    printf("\nCOMMUNICATION\n");
    printf("-------------\n");
    printf("SIDE: mpi: ");
    for (i = 0; i<mpi->S_ncom; ++i)
      printf("%d",mpi->S_com[i]);
    printf(" copy: ");
    for (i = 0; i<mpi->S_ncpy; ++i)
      printf("%d",mpi->S_cpy[i]);
    printf("\n");
    printf("EDGE: mpi: ");
    for (i = 0; i<mpi->E_ncom; ++i)
      printf("%d",mpi->E_com[i]);
    printf(" copy: ");
    for (i = 0; i<mpi->E_ncpy; ++i)
      printf("%d",mpi->E_cpy[i]);
    printf("\n");
    printf("CORNER: mpi: ");
    for (i = 0; i<mpi->C_ncom; ++i)
      printf("%d",mpi->C_com[i]);
    printf(" copy: ");
    for (i = 0; i<mpi->C_ncpy; ++i)
      printf("%d",mpi->C_cpy[i]);
    printf("\n\n");
  }
#endif //_DEBUG_MPI_

  return GOOD;
}


void communicate(real *fg, int n_c, struct System *sys)
{
#ifdef _2D_RUN_
  return;
#endif
  struct Mpi *mpi = sys->mpi;
  int step_f = sys->step->f_phase;
  int i,ii,j,k,tag,cc,node,v;
  MPI_Status status;
  real *vec;
  /* printf("<%d> COMMUNICATE\n",mpi->my_rank);fflush(stdout); */
  /* MPI_Finalize(); */
  /* exit(1); */
  //int p0 = mpi->my_rank;
  //  int *ind = mpi->ind_rank;
  real *snd0,*rcv0,*snd1,*rcv1;

  for( cc = 0; cc < n_c; ++cc ) { // loop over phases
    vec     = &fg[step_f*cc];

    /* SIDES */
    for (ii = 0; ii < mpi->S_ncom; ++ii) {
      i = mpi->S_com[ii];
      if ( mpi->S[i] ) {
        MPI_Ssend(&(vec[ mpi->side_snd[i][1] ]),1,mpi->slice[i],mpi->side_proc[1][i],i*2+1,MPI_COMM_WORLD);
        MPI_Recv (&(vec[ mpi->side_rcv[i][1] ]),1,mpi->slice[i],mpi->side_proc[1][i],i*2+2,MPI_COMM_WORLD,&status);
        MPI_Recv (&(vec[ mpi->side_rcv[i][0] ]),1,mpi->slice[i],mpi->side_proc[0][i],i*2+1,MPI_COMM_WORLD,&status);
        MPI_Ssend(&(vec[ mpi->side_snd[i][0] ]),1,mpi->slice[i],mpi->side_proc[0][i],i*2+2,MPI_COMM_WORLD);
      } else {
        MPI_Recv (&(vec[ mpi->side_rcv[i][0] ]),1,mpi->slice[i],mpi->side_proc[0][i],i*2+1,MPI_COMM_WORLD,&status);
        MPI_Ssend(&(vec[ mpi->side_snd[i][0] ]),1,mpi->slice[i],mpi->side_proc[0][i],i*2+2,MPI_COMM_WORLD);
        MPI_Ssend(&(vec[ mpi->side_snd[i][1] ]),1,mpi->slice[i],mpi->side_proc[1][i],i*2+1,MPI_COMM_WORLD);
        MPI_Recv (&(vec[ mpi->side_rcv[i][1] ]),1,mpi->slice[i],mpi->side_proc[1][i],i*2+2,MPI_COMM_WORLD,&status);
      }
    }
    for (ii = 0; ii < mpi->S_ncpy; ++ii) {
      i = mpi->S_cpy[ii];
      snd0 = &vec[ mpi->side_snd[i][0] ];
      rcv0 = &vec[ mpi->side_rcv[i][0] ];
      snd1 = &vec[ mpi->side_snd[i][1] ];
      rcv1 = &vec[ mpi->side_rcv[i][1] ];
      for (j=0; j<mpi->cpy_slice_size[i]; ++j) {
        node = mpi->cpy_slice[i][j];
        for (v=0; v<sys->n_Q; ++v) {
          rcv0[node] = snd1[node];
          rcv1[node] = snd0[node];
          node++;
        }
      }
    }
    /* #endif */

    /* EDGES */
    for (ii = 0; ii < mpi->E_ncom; ++ii) { // loop over x-,y-,z-lines
      i = mpi->E_com[ii];
      for (j = 0; j < 2; ++j) {
        tag = i*4+j*2;
        k = 2*i+j;
        if ( mpi->E[i] ) {
          MPI_Ssend(&(vec[ mpi->edge_snd[i][j*2  ] ]),1,mpi->line[i],mpi->edge_proc[1][k],tag+2,MPI_COMM_WORLD);
          MPI_Recv (&(vec[ mpi->edge_rcv[i][j*2  ] ]),1,mpi->line[i],mpi->edge_proc[1][k],tag+1,MPI_COMM_WORLD,&status);
          MPI_Recv (&(vec[ mpi->edge_rcv[i][j*2+1] ]),1,mpi->line[i],mpi->edge_proc[0][k],tag+2,MPI_COMM_WORLD,&status);
          MPI_Ssend(&(vec[ mpi->edge_snd[i][j*2+1] ]),1,mpi->line[i],mpi->edge_proc[0][k],tag+1,MPI_COMM_WORLD);
        } else {
          MPI_Recv (&(vec[ mpi->edge_rcv[i][j*2+1] ]),1,mpi->line[i],mpi->edge_proc[0][k],tag+2,MPI_COMM_WORLD,&status);
          MPI_Ssend(&(vec[ mpi->edge_snd[i][j*2+1] ]),1,mpi->line[i],mpi->edge_proc[0][k],tag+1,MPI_COMM_WORLD);
          MPI_Ssend(&(vec[ mpi->edge_snd[i][j*2  ] ]),1,mpi->line[i],mpi->edge_proc[1][k],tag+2,MPI_COMM_WORLD);
          MPI_Recv (&(vec[ mpi->edge_rcv[i][j*2  ] ]),1,mpi->line[i],mpi->edge_proc[1][k],tag+1,MPI_COMM_WORLD,&status);
        }
      }
    }
    for (ii = 0; ii < mpi->E_ncpy; ++ii) {
      i = mpi->E_cpy[ii];
      for (j = 0; j < 2; ++j) {
        snd0 = &vec[ mpi->edge_snd[i][j*2  ] ];
        rcv0 = &vec[ mpi->edge_rcv[i][j*2  ] ];
        snd1 = &vec[ mpi->edge_snd[i][j*2+1] ];
        rcv1 = &vec[ mpi->edge_rcv[i][j*2+1] ];
        for (k=0; k<mpi->cpy_line_size[i]; ++k) {
          node = mpi->cpy_line[i][k];
          for (v=0; v<sys->n_Q; ++v) {
            rcv0[node] = snd1[node];
            rcv1[node] = snd0[node];
            node++;
          }
        }
      }
    }

    /* CORNERS */
    for (i = 0; i < mpi->C_ncom; ++i) {
      if ( mpi->C ) {
        MPI_Ssend(&(vec[ mpi->crnr_snd[i][0] ]),sys->n_Q,MPI_DOUBLE,mpi->crnr_proc[1][i],i*2+2,MPI_COMM_WORLD);
        MPI_Recv (&(vec[ mpi->crnr_rcv[i][0] ]),sys->n_Q,MPI_DOUBLE,mpi->crnr_proc[1][i],i*2+1,MPI_COMM_WORLD,&status);
        MPI_Recv (&(vec[ mpi->crnr_rcv[i][1] ]),sys->n_Q,MPI_DOUBLE,mpi->crnr_proc[0][i],i*2+2,MPI_COMM_WORLD,&status);
        MPI_Ssend(&(vec[ mpi->crnr_snd[i][1] ]),sys->n_Q,MPI_DOUBLE,mpi->crnr_proc[0][i],i*2+1,MPI_COMM_WORLD);
      } else {
        MPI_Recv (&(vec[ mpi->crnr_rcv[i][1] ]),sys->n_Q,MPI_DOUBLE,mpi->crnr_proc[0][i],i*2+2,MPI_COMM_WORLD,&status);
        MPI_Ssend(&(vec[ mpi->crnr_snd[i][1] ]),sys->n_Q,MPI_DOUBLE,mpi->crnr_proc[0][i],i*2+1,MPI_COMM_WORLD);
        MPI_Ssend(&(vec[ mpi->crnr_snd[i][0] ]),sys->n_Q,MPI_DOUBLE,mpi->crnr_proc[1][i],i*2+2,MPI_COMM_WORLD);
        MPI_Recv (&(vec[ mpi->crnr_rcv[i][0] ]),sys->n_Q,MPI_DOUBLE,mpi->crnr_proc[1][i],i*2+1,MPI_COMM_WORLD,&status);
      }
    }
    for (ii = 0; ii < mpi->C_ncpy; ++ii) {
      i = mpi->C_cpy[ii];
      snd0 = &(vec[ mpi->crnr_snd[i][0] ]);
      rcv0 = &(vec[ mpi->crnr_rcv[i][0] ]);
      snd1 = &(vec[ mpi->crnr_snd[i][1] ]);
      rcv1 = &(vec[ mpi->crnr_rcv[i][1] ]);
      node = 0;
      for (v=0; v<sys->n_Q; ++v) {
        rcv0[node] = snd1[node];
        rcv1[node] = snd0[node];
        node++;
      }
    }

  }
}


//-----------------------------------------
//
//-----------------------------------------
void communicate_char(char *data, struct System *sys)
{
  //static int cnt = 0;
  struct Mpi *mpi = sys->mpi;
  //struct Node *node = sys->node;

  int i,ii,j,k,tag,cn,pn,nn;
  MPI_Status status;
  char *snd0,*rcv0,*snd1,*rcv1;
  
  /* SIDES */
  for (ii = 0; ii < mpi->S_ncom; ++ii) {
    i = mpi->S_com[ii];
    if ( mpi->S[i] ) {
      MPI_Ssend(&(data[ mpi->side_snd_scalar[i][1] ]),1,mpi->slice_char[i],mpi->side_proc[1][i],i*2+1,MPI_COMM_WORLD);
      MPI_Recv (&(data[ mpi->side_rcv_scalar[i][1] ]),1,mpi->slice_char[i],mpi->side_proc[1][i],i*2+2,MPI_COMM_WORLD,&status);
      MPI_Recv (&(data[ mpi->side_rcv_scalar[i][0] ]),1,mpi->slice_char[i],mpi->side_proc[0][i],i*2+1,MPI_COMM_WORLD,&status);
      MPI_Ssend(&(data[ mpi->side_snd_scalar[i][0] ]),1,mpi->slice_char[i],mpi->side_proc[0][i],i*2+2,MPI_COMM_WORLD);
    } else {
      MPI_Recv (&(data[ mpi->side_rcv_scalar[i][0] ]),1,mpi->slice_char[i],mpi->side_proc[0][i],i*2+1,MPI_COMM_WORLD,&status);
      MPI_Ssend(&(data[ mpi->side_snd_scalar[i][0] ]),1,mpi->slice_char[i],mpi->side_proc[0][i],i*2+2,MPI_COMM_WORLD);
      MPI_Ssend(&(data[ mpi->side_snd_scalar[i][1] ]),1,mpi->slice_char[i],mpi->side_proc[1][i],i*2+1,MPI_COMM_WORLD);
      MPI_Recv (&(data[ mpi->side_rcv_scalar[i][1] ]),1,mpi->slice_char[i],mpi->side_proc[1][i],i*2+2,MPI_COMM_WORLD,&status);
    }
  }
  for (ii = 0; ii < mpi->S_ncpy; ++ii) {
    i = mpi->S_cpy[ii];

    snd0 = &data[ mpi->side_snd_scalar[i][0] ];
    rcv0 = &data[ mpi->side_rcv_scalar[i][0] ];
    snd1 = &data[ mpi->side_snd_scalar[i][1] ];
    rcv1 = &data[ mpi->side_rcv_scalar[i][1] ];

    for (j=0; j<mpi->cpy_slice_size[i]; ++j) {
      nn = mpi->cpy_slice_char[i][j];
      rcv0[nn] = snd1[nn];
      rcv1[nn] = snd0[nn];
    }
  }

  /* EDGES */
  for (ii = 0; ii < mpi->E_ncom; ++ii) { // loop over x-,y-,z-lines
    i = mpi->E_com[ii];
    for (j = 0; j < 2; ++j) {
      tag = i*4+j*2;
      k = 2*i+j;
      if ( mpi->E[i] ) {
        //printf("A1 (%d) %d --> %d\n",i,mpi->my_rank,                                  mpi->edge_proc[1][k]);fflush(stdout);
        MPI_Ssend(&(data[ mpi->edge_snd_scalar[i][j*2  ] ]),1,mpi->line_char[i],mpi->edge_proc[1][k],tag+2,MPI_COMM_WORLD);
        //printf("A2 (%d) %d <-- %d\n",i,mpi->my_rank,                                  mpi->edge_proc[1][k]);fflush(stdout);
        MPI_Recv (&(data[ mpi->edge_rcv_scalar[i][j*2  ] ]),1,mpi->line_char[i],mpi->edge_proc[1][k],tag+1,MPI_COMM_WORLD,&status);
        //printf("A3 (%d) %d <-- %d\n",i,mpi->my_rank,                                  mpi->edge_proc[0][k]);fflush(stdout);
        MPI_Recv (&(data[ mpi->edge_rcv_scalar[i][j*2+1] ]),1,mpi->line_char[i],mpi->edge_proc[0][k],tag+2,MPI_COMM_WORLD,&status);
        //printf("A4 (%d) %d --> %d\n",i,mpi->my_rank,                                  mpi->edge_proc[0][k]);fflush(stdout);
        MPI_Ssend(&(data[ mpi->edge_snd_scalar[i][j*2+1] ]),1,mpi->line_char[i],mpi->edge_proc[0][k],tag+1,MPI_COMM_WORLD);
      } else {
        //printf("B1 (%d) %d <-- %d\n",i,mpi->my_rank,                                  mpi->edge_proc[0][k]);fflush(stdout);
        MPI_Recv (&(data[ mpi->edge_rcv_scalar[i][j*2+1] ]),1,mpi->line_char[i],mpi->edge_proc[0][k],tag+2,MPI_COMM_WORLD,&status);
        //printf("B2 (%d) %d --> %d\n",i,mpi->my_rank,                                  mpi->edge_proc[0][k]);fflush(stdout);
        MPI_Ssend(&(data[ mpi->edge_snd_scalar[i][j*2+1] ]),1,mpi->line_char[i],mpi->edge_proc[0][k],tag+1,MPI_COMM_WORLD);
        //printf("B3 (%d) %d --> %d\n",i,mpi->my_rank,                                  mpi->edge_proc[1][k]);fflush(stdout);
        MPI_Ssend(&(data[ mpi->edge_snd_scalar[i][j*2  ] ]),1,mpi->line_char[i],mpi->edge_proc[1][k],tag+2,MPI_COMM_WORLD);
        //printf("B4 (%d) %d <-- %d\n",i,mpi->my_rank,                                  mpi->edge_proc[1][k]);fflush(stdout);
        MPI_Recv (&(data[ mpi->edge_rcv_scalar[i][j*2  ] ]),1,mpi->line_char[i],mpi->edge_proc[1][k],tag+1,MPI_COMM_WORLD,&status);
      }
    }
  }
  for (ii = 0; ii < mpi->E_ncpy; ++ii) {
    i = mpi->E_cpy[ii];
    for (j = 0; j < 2; ++j) {
      snd0 = &data[ mpi->edge_snd_scalar[i][j*2  ] ];
      rcv0 = &data[ mpi->edge_rcv_scalar[i][j*2  ] ];
      snd1 = &data[ mpi->edge_snd_scalar[i][j*2+1] ];
      rcv1 = &data[ mpi->edge_rcv_scalar[i][j*2+1] ];
      for (k=0; k<mpi->cpy_line_size[i]; ++k) {
        nn = mpi->cpy_line_char[i][k];
        rcv0[nn] = snd1[nn];
        rcv1[nn] = snd0[nn];
      }
    }
  }

  /* CORNERS */
  for (i = 0; i < mpi->C_ncom; ++i) {
    if ( mpi->C ) {
      MPI_Ssend(&(data[ mpi->crnr_snd_scalar[i][0] ]),1,MPI_CHAR,mpi->crnr_proc[1][i],i*2+2,MPI_COMM_WORLD);
      MPI_Recv (&(data[ mpi->crnr_rcv_scalar[i][0] ]),1,MPI_CHAR,mpi->crnr_proc[1][i],i*2+1,MPI_COMM_WORLD,&status);
      MPI_Recv (&(data[ mpi->crnr_rcv_scalar[i][1] ]),1,MPI_CHAR,mpi->crnr_proc[0][i],i*2+2,MPI_COMM_WORLD,&status);
      MPI_Ssend(&(data[ mpi->crnr_snd_scalar[i][1] ]),1,MPI_CHAR,mpi->crnr_proc[0][i],i*2+1,MPI_COMM_WORLD);
    } else {
      MPI_Recv (&(data[ mpi->crnr_rcv_scalar[i][1] ]),1,MPI_CHAR,mpi->crnr_proc[0][i],i*2+2,MPI_COMM_WORLD,&status);
      MPI_Ssend(&(data[ mpi->crnr_snd_scalar[i][1] ]),1,MPI_CHAR,mpi->crnr_proc[0][i],i*2+1,MPI_COMM_WORLD);
      MPI_Ssend(&(data[ mpi->crnr_snd_scalar[i][0] ]),1,MPI_CHAR,mpi->crnr_proc[1][i],i*2+2,MPI_COMM_WORLD);
      MPI_Recv (&(data[ mpi->crnr_rcv_scalar[i][0] ]),1,MPI_CHAR,mpi->crnr_proc[1][i],i*2+1,MPI_COMM_WORLD,&status);
    }
  }
  //cnt++;
  for (ii = 0; ii < mpi->C_ncpy; ++ii) {
    i = mpi->C_cpy[ii];
    snd0 = &(data[ mpi->crnr_snd_scalar[i][0] ]);
    rcv0 = &(data[ mpi->crnr_rcv_scalar[i][0] ]);
    snd1 = &(data[ mpi->crnr_snd_scalar[i][1] ]);
    rcv1 = &(data[ mpi->crnr_rcv_scalar[i][1] ]);
    *rcv0 = *snd1;
    *rcv1 = *snd0;
    //printf("%d = %d, %d = %d\n", *rcv0, *snd1, *rcv1, *snd0);
    //printf("%d: %d = %d, %d = %d\n", cnt, mpi->crnr_rcv_scalar[i][0], mpi->crnr_snd_scalar[i][1], mpi->crnr_rcv_scalar[i][1], mpi->crnr_snd_scalar[i][0]);
  }
}

//-----------------------------------------
//
//-----------------------------------------
void communicate_node(struct Node *node, struct System *sys)
{
#ifdef _2D_RUN_
  return;
#endif
  communicate_char(node->mode, sys);

  // add 5 to the periodic nodes so that they
  // can be recognized as FLUID, WALL, or SOLID
  for(int cn = 0; cn < node->nr_per; ++cn ) {
    int pn = node->list_per[2*cn];
    node->mode[pn] += PERIODIC;
  }
}


//-----------------------------------------
//
//-----------------------------------------
void communicate_double(double *data, struct System *sys)
{
  struct Mpi *mpi = sys->mpi;
  MPI_Status status;

  /* SIDES */
  for (int ii = 0; ii < mpi->S_ncom; ++ii) {
    int i = mpi->S_com[ii];
    if ( mpi->S[i] ) {
      MPI_Ssend(&data[ mpi->side_snd_scalar[i][1] ],1,mpi->slice_double[i],mpi->side_proc[1][i],i*2+1,MPI_COMM_WORLD);
      MPI_Recv (&data[ mpi->side_rcv_scalar[i][1] ],1,mpi->slice_double[i],mpi->side_proc[1][i],i*2+2,MPI_COMM_WORLD,&status);
      MPI_Recv (&data[ mpi->side_rcv_scalar[i][0] ],1,mpi->slice_double[i],mpi->side_proc[0][i],i*2+1,MPI_COMM_WORLD,&status);
      MPI_Ssend(&data[ mpi->side_snd_scalar[i][0] ],1,mpi->slice_double[i],mpi->side_proc[0][i],i*2+2,MPI_COMM_WORLD);
    } else {
      MPI_Recv (&data[ mpi->side_rcv_scalar[i][0] ],1,mpi->slice_double[i],mpi->side_proc[0][i],i*2+1,MPI_COMM_WORLD,&status);
      MPI_Ssend(&data[ mpi->side_snd_scalar[i][0] ],1,mpi->slice_double[i],mpi->side_proc[0][i],i*2+2,MPI_COMM_WORLD);
      MPI_Ssend(&data[ mpi->side_snd_scalar[i][1] ],1,mpi->slice_double[i],mpi->side_proc[1][i],i*2+1,MPI_COMM_WORLD);
      MPI_Recv (&data[ mpi->side_rcv_scalar[i][1] ],1,mpi->slice_double[i],mpi->side_proc[1][i],i*2+2,MPI_COMM_WORLD,&status);
    }
  }
  for (int ii = 0; ii < mpi->S_ncpy; ++ii) {
    int i = mpi->S_cpy[ii];

    double *snd0 = &data[ mpi->side_snd_scalar[i][0] ];
    double *rcv0 = &data[ mpi->side_rcv_scalar[i][0] ];
    double *snd1 = &data[ mpi->side_snd_scalar[i][1] ];
    double *rcv1 = &data[ mpi->side_rcv_scalar[i][1] ];

    for (int j=0; j<mpi->cpy_slice_size[i]; ++j) {
      int nn = mpi->cpy_slice_char[i][j];
      rcv0[nn] = snd1[nn];
      rcv1[nn] = snd0[nn];
    }
  }

  /* EDGES */
  for (int ii = 0; ii < mpi->E_ncom; ++ii) { // loop over x-,y-,z-lines
    int i = mpi->E_com[ii];
    for (int j = 0; j < 2; ++j) {
      int tag = i*4+j*2;
      int k = 2*i+j;
      if ( mpi->E[i] ) {
        MPI_Ssend(&(data[ mpi->edge_snd_scalar[i][j*2  ] ]),1,mpi->line_double[i],mpi->edge_proc[1][k],tag+2,MPI_COMM_WORLD);
        MPI_Recv (&(data[ mpi->edge_rcv_scalar[i][j*2  ] ]),1,mpi->line_double[i],mpi->edge_proc[1][k],tag+1,MPI_COMM_WORLD,&status);
        MPI_Recv (&(data[ mpi->edge_rcv_scalar[i][j*2+1] ]),1,mpi->line_double[i],mpi->edge_proc[0][k],tag+2,MPI_COMM_WORLD,&status);
        MPI_Ssend(&(data[ mpi->edge_snd_scalar[i][j*2+1] ]),1,mpi->line_double[i],mpi->edge_proc[0][k],tag+1,MPI_COMM_WORLD);
      } else {
        MPI_Recv (&(data[ mpi->edge_rcv_scalar[i][j*2+1] ]),1,mpi->line_double[i],mpi->edge_proc[0][k],tag+2,MPI_COMM_WORLD,&status);
        MPI_Ssend(&(data[ mpi->edge_snd_scalar[i][j*2+1] ]),1,mpi->line_double[i],mpi->edge_proc[0][k],tag+1,MPI_COMM_WORLD);
        MPI_Ssend(&(data[ mpi->edge_snd_scalar[i][j*2  ] ]),1,mpi->line_double[i],mpi->edge_proc[1][k],tag+2,MPI_COMM_WORLD);
        MPI_Recv (&(data[ mpi->edge_rcv_scalar[i][j*2  ] ]),1,mpi->line_double[i],mpi->edge_proc[1][k],tag+1,MPI_COMM_WORLD,&status);
      }
    }
  }
  for (int ii = 0; ii < mpi->E_ncpy; ++ii) {
    int i = mpi->E_cpy[ii];
    for (int j = 0; j < 2; ++j) {
      double *snd0 = &data[ mpi->edge_snd_scalar[i][j*2  ] ];
      double *rcv0 = &data[ mpi->edge_rcv_scalar[i][j*2  ] ];
      double *snd1 = &data[ mpi->edge_snd_scalar[i][j*2+1] ];
      double *rcv1 = &data[ mpi->edge_rcv_scalar[i][j*2+1] ];
      for (int k=0; k<mpi->cpy_line_size[i]; ++k) {
        int nn = mpi->cpy_line_char[i][k];
        rcv0[nn] = snd1[nn];
        rcv1[nn] = snd0[nn];
      }
    }
  }

  /* CORNERS */
  for (int i = 0; i < mpi->C_ncom; ++i) {
    if ( mpi->C ) {
      MPI_Ssend(&(data[ mpi->crnr_snd_scalar[i][0] ]),1,MPI_DOUBLE,mpi->crnr_proc[1][i],i*2+2,MPI_COMM_WORLD);
      MPI_Recv (&(data[ mpi->crnr_rcv_scalar[i][0] ]),1,MPI_DOUBLE,mpi->crnr_proc[1][i],i*2+1,MPI_COMM_WORLD,&status);
      MPI_Recv (&(data[ mpi->crnr_rcv_scalar[i][1] ]),1,MPI_DOUBLE,mpi->crnr_proc[0][i],i*2+2,MPI_COMM_WORLD,&status);
      MPI_Ssend(&(data[ mpi->crnr_snd_scalar[i][1] ]),1,MPI_DOUBLE,mpi->crnr_proc[0][i],i*2+1,MPI_COMM_WORLD);
    } else {
      MPI_Recv (&(data[ mpi->crnr_rcv_scalar[i][1] ]),1,MPI_DOUBLE,mpi->crnr_proc[0][i],i*2+2,MPI_COMM_WORLD,&status);
      MPI_Ssend(&(data[ mpi->crnr_snd_scalar[i][1] ]),1,MPI_DOUBLE,mpi->crnr_proc[0][i],i*2+1,MPI_COMM_WORLD);
      MPI_Ssend(&(data[ mpi->crnr_snd_scalar[i][0] ]),1,MPI_DOUBLE,mpi->crnr_proc[1][i],i*2+2,MPI_COMM_WORLD);
      MPI_Recv (&(data[ mpi->crnr_rcv_scalar[i][0] ]),1,MPI_DOUBLE,mpi->crnr_proc[1][i],i*2+1,MPI_COMM_WORLD,&status);
    }
  }
  for (int ii = 0; ii < mpi->C_ncpy; ++ii) {
    int i = mpi->C_cpy[ii];
    double *snd0 = &(data[ mpi->crnr_snd_scalar[i][0] ]);
    double *rcv0 = &(data[ mpi->crnr_rcv_scalar[i][0] ]);
    double *snd1 = &(data[ mpi->crnr_snd_scalar[i][1] ]);
    double *rcv1 = &(data[ mpi->crnr_rcv_scalar[i][1] ]);
    *rcv0 = *snd1;
    *rcv1 = *snd0;
  }

}


//-----------------------------------------
//
//-----------------------------------------
void communicate_rho(Field *field, System *sys)
{
#ifdef _2D_RUN_
  return;
#endif
  for (int cc = 0; cc < field->n_c; ++cc) {
    double *data = &field->rho[sys->step->rho_phase*cc];
    communicate_double(data, sys);
  }
}


//-----------------------------------------
//
//-----------------------------------------
void communicate_rho_tot(Field *field, System *sys)
{
#ifdef _2D_RUN_
  return;
#endif
  communicate_double(field->rho_tot, sys);
}


//-----------------------------------------
//
//-----------------------------------------
void communicate_sfrac(struct System *sys, struct Minerals *minerals)
{
#ifdef _2D_RUN_
  return;
#endif
  for (int min = 0; min < minerals->n_tot; ++min) {
    double *data = minerals->list[min].sfrac;
    communicate_double(data, sys);
  }
}



//-----------------------------------------
//
//-----------------------------------------
void communicate_mineral_delta(int min, struct System *sys, struct Minerals *minerals)
{
#ifdef _2D_RUN_
  return;
#endif
  communicate_double(minerals->list[min].delta, sys);
}


//-----------------------------------------
//
//-----------------------------------------
void communicate_mineral_sfrac(int min, struct System *sys, struct Minerals *minerals)
{
#ifdef _2D_RUN_
  return;
#endif
  communicate_double(minerals->list[min].sfrac, sys);
}

//-----------------------------------------
//
//-----------------------------------------
void setup_ghost_nodes(struct Node *node, struct System *sys)
// The routine also fixes the list of periodic nodes
//
{
  int i;
  struct Mpi *mpi = sys->mpi;
  //struct Node *node = sys->node;
  int *ind = mpi->ind_rank;
  int *np = mpi->np;
  int *n = sys->max_n;

  node->nr_ghost=0;

  //  if (sys->n_D<3) {
  //    sys->max_n[2]=1;
  //  }

  // number of ghosts for different faces
  int ng[3] = {n[1]*n[2], n[0]*n[2], n[0]*n[1]};

  // get number of ghost nodes (n_ghost_)
  // for each sub-process
  for (i = 0; i < sys->n_D; ++i) {
    if (ind[i]<np[i]-1) {
      // sub-procs: NOT rightmost, topmost, backmost
      // faces: right(i=0) / top(i=1) / back(i=2)
      node->nr_ghost += ng[i];
    }
    if (ind[i]>0) {
      // sub-procs: NOT leftmost, bottomost, frontmost
      // faces: left(i=0) / bottom(i=1) / front(i=2)
      node->nr_ghost += ng[i];
    }
  }

  // nr_ghost is now too long since we have counted edges and corners twice,
  // but nr_ghost is updated to the correct value after ghost_loop
  if ( !( node->list_ghost = (int *) malloc(node->nr_ghost*sizeof(int)) )  ) {
    printf("error allocating memory for 'node->list_ghost'\n");
    exit(1);
  }
  // init array to -1 to otherwise node 0 may not be added
  // to node->list_ghost
  for (i = 0; i < node->nr_ghost; ++i) {
    node->list_ghost[i] = -1;
  }

  //int cc = 0;

  /*   int first_rtb[3] = {n[0]-1, n[1]-1, n[2]-1}; */
  /*   int last_rtb [3] = {n[0]  , n[1]  , n[2]}; */

  /*   int first_lbf[3] = {0,0,0}; */
  /*   int last_lbf [3] = {n[0],n[1],n[2]}; */
  int c=0,start,s[3],a[2];
  //char type[10];

  for (i = 0; i < sys->n_D; ++i) {
    if (ind[i]<np[i]-1) {

      // right
      if (i==0) {
        start = n[0]-1;
        s[0]=0; s[1]=1; s[2]=2;
        a[0]=1; a[1]=2;
        ghost_loop(node, sys,a,s,&c,start,"R");
        //if (mpi->my_rank==0)
        //  printf("R, c: %d\n",c);
      }
      // top
      if (i==1) {
        start = n[1]-1;
        s[0]=1; s[1]=0; s[2]=2;
        a[0]=0; a[1]=2;
        ghost_loop(node, sys,a,s,&c,start,"T");
        //if (mpi->my_rank==0)
        //  printf("T, c: %d\n",c);
      }

      // back
      if (i==2) {
        start = n[2]-1;
        s[0]=2; s[1]=0; s[2]=1;
        a[0]=0; a[1]=1;
        ghost_loop(node, sys,a,s,&c,start,"B");
        //if (mpi->my_rank==0)
        //  printf("B, c: %d\n",c);
      }
    }

    if (ind[i]>0) {

      // left
      if (i==0) {
        s[0]=0; s[1]=1; s[2]=2;
        a[0]=1; a[1]=2;
        ghost_loop(node, sys,a,s,&c,0,"L");
      }

      // bottom
      if (i==1) {
        s[0]=1; s[1]=0; s[2]=2;
        a[0]=0; a[1]=2;
        ghost_loop(node, sys,a,s,&c,0,"B");
      }

      // front
      if (i==2) {
        s[0]=2; s[1]=0; s[2]=1;
        a[0]=0; a[1]=1;
        ghost_loop(node, sys,a,s,&c,0,"F");
      }
    } 
  }

  node->nr_ghost = c;

  // reallocate per_nodes_
  int *per_nodes_new;
  int n_per_new = node->nr_per - node->nr_ghost;
  if ( !( per_nodes_new = (int *) calloc(2*n_per_new, sizeof(int)) )  ) {
    printf("ERROR! re-allocating memory for 'node->list_per' failed!\n");
    exit(1);
  }

  int ii,missed=1;
  c = 0;
  for( i = 0; i < node->nr_per; ++i ) {
    missed = 1;
    for( ii = 0; ii < node->nr_ghost; ++ii ) {
      if ( node->list_per[i*2] == node->list_ghost[ii] ) {
        missed = 0;
        break;
      }
    }
    if (missed) {
      per_nodes_new[c++] = node->list_per[i*2];
      per_nodes_new[c++] = node->list_per[i*2+1];
    }
  }


  int check = c - 2*n_per_new;
  if (check != 0) {
    fprintf(stderr,"<%d> ERROR! Index mismatch in setup_ghost_nodes()\n",mpi->my_rank);
    fprintf(stderr,"<%d> ERROR! c - 2*n_per_new = %d\n",mpi->my_rank,check);
    exit(1);
  }

  //printf("<%d> node->nr_per: %d, n_per_new: %d\n",mpi->my_rank,node->nr_per,n_per_new);
  node->nr_per = n_per_new;
  free(node->list_per);
  node->list_per = per_nodes_new;

  int cn, pn;
  for( cn = 0; cn < node->nr_ghost; ++cn ) {
    pn = node->list_ghost[cn];
    node->is_ghost[pn] = node->is_ghost_or_periodic[pn] = 1;
  }
  for( cn = 0; cn < node->nr_per; ++cn ) {
    pn = node->list_per[2*cn];
    node->is_ghost_or_periodic[pn] = 1;
  }
  return;
}


//-----------------------------------------
//
//-----------------------------------------
void add_ghost_node(struct Node *node, int nn, int *c)
/* Add node nn to the ghost-list *
 * only if it is new             */
{
  int i;
  for (i = 0; i < node->nr_ghost; ++i) {
    if (nn == node->list_ghost[i])
      return;
  }
  node->mode[nn] = GHOST;
  node->list_ghost[(*c)++] = nn;
  return;
}

//-----------------------------------------
//
//-----------------------------------------
void ghost_loop(struct Node *node, struct System *sys, int *a, int *s, int *c, int start, const char *str)
{
  struct Mpi *mpi = sys->mpi;
  //struct Node *node = sys->node;

  int *ind = mpi->ind_rank; // rank x,y,z-index
  int *np = mpi->np;        // number of procs in x,y,z dir
  int *n = sys->max_n;
  int imin[3] = {(ind[0]==0)      ,(ind[1]==0)      ,(ind[2]==0)};
  int imax[3] = {(ind[0]==np[0]-1),(ind[1]==np[1]-1),(ind[2]==np[2]-1)};
  int A;
  int first[3] = {0,0,0};
  int last [3] = {n[0],n[1],n[2]};
  int stride[3] = {1,n[0],n[0]*n[1]};
  int d[3];
  int i,nn;

  // uncomment this block to avoid ghost nodes
  // on the "surface" of the system. This block
  // change the values of first and last in the
  // for-loops below to 1 and n-1 instead of 0 and n
  for (i=0; i<2; ++i) {
    A = a[i];
    if (imin[A]) {
      first[A] += 1;
    }
    if (imax[A]) {
      last[A] -= 1;
    }
  }

  d[s[0]] = start;
  for (d[s[1]] = first[s[1]]; d[s[1]] < last[s[1]]; ++d[s[1]]) {
    for (d[s[2]] = first[s[2]]; d[s[2]] < last[s[2]]; ++d[s[2]]) {
      nn = 0;
      for (i=0; i<3; ++i) {
        nn += (d[s[i]]*stride[s[i]]);
      }
      add_ghost_node(node,nn,c);
    }
  }
  /*   printf("\n",str); */
  /*   fflush(stdout); */
  return;
}

//-----------------------------------------
//
//-----------------------------------------
void set_ghost_neighbors(struct System *sys)
{
  int ng, nn, nb, ix[3], ck, cd, ndir;
  int valid;
  int *n = sys->max_n;
  int *ei;
  struct Mpi *mpi = sys->mpi;
  struct Mpi::Ghost *ghost;
  struct Node *node = sys->node;

  // allocate
  //mpi->ghost = (struct Mpi::Ghost *) calloc(node->nr_ghost,sizeof(struct Mpi::Ghost));
  mpi->ghost = new Mpi::Ghost[node->nr_ghost]();
  for( ng = 0; ng < node->nr_ghost; ++ng ) {
    mpi->ghost[ng].dir = (int *) calloc(sys->n_Q,sizeof(int));
  }
  
  /*   FILE *fp; */
  /*   char name[50]; */
  /*   sprintf(name,"%03d_ghost-dir.txt",mpi->my_rank); */
  /*   fp = fopen(name,"w"); */
  
  for( ng = 0; ng < node->nr_ghost; ++ng ) {
    ndir = 0;
    nn = node->list_ghost[ng];
    ghost = &mpi->ghost[ng];
    ix[0] = nn%n[0];
    ix[1] = ((nn-ix[0])/n[0])%n[1];
    ix[2] = (int)(nn/(n[0]*n[1]));
    /*     fprintf(fp,"(%02d,%02d,%02d) -> ",ix[0],ix[1],ix[2]); */
    for( ck = 1; ck < sys->n_Q; ++ck ) {
      ei = &sys->ei[sys->n_D*ck];
      valid = 1;
      for( cd = 0; cd < sys->n_D; ++cd ) {
	//nb = ix[cd]+ei[cd];
	nb = ix[cd]+ei[cd];
	if (nb < 0 || nb >= n[cd] ) {
	  valid = 0;
	  break;
	}
      }
      if (valid) {
	ghost->dir[ndir] = ck;
	ndir++;
      }
    }
    ghost->ndir = ndir;
  }
}


//-----------------------------------------
//
//-----------------------------------------
void broadcast_global_max(struct Mpi *mpi, double *max)
{
  int nn;
  // All processes send double to root-process using double_rbuf
  MPI_Gather(max, 1, MPI_DOUBLE, mpi->double_rbuf, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  // root-process loops through double_rbuf to find global max
  if (mpi->my_rank == 0) {
    (*max) = 0.;
    for (nn = 0; nn < mpi->nr_procs; ++nn) {
      if (mpi->double_rbuf[nn] > *max) {
        (*max) = mpi->double_rbuf[nn];
      }
    }
  }
  
  // root-process broadcasts max-value to all processes
  MPI_Bcast(max, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
}

////////////////
// calculate global max at process 0
////////////////
int get_global_max(struct Mpi *mpi, double *max)
{
  int nn, proc;
  // All processes send double to root-process using double_rbuf
  MPI_Gather(max, 1, MPI_DOUBLE, mpi->double_rbuf, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  // root-process loops through double_rbuf to find global max
  if (mpi->my_rank == 0) {
    (*max) = -1.0e10;
    for (nn = 0; nn < mpi->nr_procs; ++nn) {
      if (mpi->double_rbuf[nn] > *max) {
        (*max) = mpi->double_rbuf[nn];
        proc = nn;
      }
    }
  }
  MPI_Bcast(&(proc), 1, MPI_INT, 0, MPI_COMM_WORLD);
  return(proc);
}

////////////////

////////////////
void broadcast_global_min(struct Mpi *mpi, double *min)
{
  int nn;
  // All processes send double to root-process using double_rbuf
  MPI_Gather(min, 1, MPI_DOUBLE, mpi->double_rbuf, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  // root-process loops through double_rbuf to find global min
  if (mpi->my_rank == 0) {
    *min = mpi->double_rbuf[0];
    for (nn = 1; nn < mpi->nr_procs; ++nn) {
      if (mpi->double_rbuf[nn] < *min) {
	*min = mpi->double_rbuf[nn];
      }
    }

#ifdef _DEBUG_
    FILE *fp;
    char fn[50];
    for (nn = 1; nn < mpi->nr_procs; ++nn) {
      if (mpi->double_rbuf[nn] == *min) {
	sprintf(fn, "%s/%04d_converge.dat", (*sys->files)["out"], nn);
	fp = my_fopen(fn, "a");
	fprintf(fp, "*********** GLOBAL MIN IN PROCESS <%03d> **********\n", nn);
	fclose(fp);
	break;
      }
    }
#endif
    
  }

  // root-process broadcasts min-value to all processes
  MPI_Bcast(min, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
}

