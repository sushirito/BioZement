#include "global.h"
#include "biozementfunctions.h"


//--------------------------------------------
// Read geo.dat, set up wall-nodes, inlet/outlet-nodes,
// and allocate and initialize node-lists
// geo.dat only needs to specify solid (1) and fluid (0) nodes
// and wall-nodes are set up by set_wall_nodes().
// 
// Note that communicate_node also updates the periodic nodes, but
// the periodic nodes should be left untoched at this stage. This is
// why each call to communicate_node() is followed by a "reset"
// of the periodic nodes.
//--------------------------------------------
void init_geometry_and_nodes(struct System *sys)
{
  struct Node *node = sys->node;
  struct Options *opt = sys->opt;

  if (sys->opt->measure_inlet_perm == 0) {
    // check if size of geofile match dimensions in inp.dat
    check_size_of_file(sys);

    read_geofile(node, sys);
    read_Rfile(node, sys); /* AMIN */
    read_Bfile(node, sys); /* AMIN */
  }

  init_node_lists(node, sys);

  // need to set periodic nodes for 2D runs explicitly since
  // mpi-communication is turned off for 2D runs
  if (sys->n_D < 3) {
    set_periodic_nodes_2D(node, sys);
  }

  if (opt->copy_in_out) {
    copy_inlet_outlet_slice(node, sys);
  }

  if (opt->add_union_in_out) {
    add_overlapping_slice_at_inlet_and_outlet(node, sys);
  }

  if (opt->wall_nodes) {
    add_outer_walls(node, sys);
  }

  if (opt->add_raster) {
    add_raster(node, sys);
  }

  set_wall_nodes(node, sys);
#ifndef FLUID_OFF // if fluid is on
  set_inlet_outlet_nodes(node, sys);
#endif

  if (opt->closed_cell) {
    close_inlet_outlet(node, sys);
  }

  set_system_dimensions(sys);

}

//--------------------------------------------
//
//--------------------------------------------
void set_system_dimensions(System *sys)
{
  Options *opt = sys->opt;
  int xmin, xmax;
  if (opt->closed_cell) {
    get_x_limits(sys->node, &xmin, &xmax, sys);
    opt->inertx_inlet = xmin;
    opt->inertx_outlet = sys->MAX_N[0]-1-xmax;
    //    printf("inertx_inlet, inertx_outlet = %d, %d\n",opt->inertx_inlet, opt->inertx_outlet );
  }

  // dimensions of reactive system
  if (opt->measure_inlet_perm) {
    sys->L_lb = opt->in_nodes + opt->out_nodes;
  } else {
    sys->L_lb = sys->MAX_N[0] - 2 - opt->inertx_inlet - opt->inertx_outlet;
  }
  //printf("sys->L_lb = sys->MAX_N[0] - opt->inertx_inlet - opt->inertx_outlet + 1 = %d - %d - %d + 1 = %d",
  //    sys->MAX_N[0], opt->inertx_inlet, opt->inertx_outlet, sys->L_lb);

  sys->W_lb = sys->MAX_N[1] - 2 - 2 * opt->wall_nodes;
  sys->H_lb = std::max(sys->MAX_N[2] - 2 - 2 * opt->wall_nodes, 0);
  if (sys->n_D > 2) {
    sys->A_lb = sys->W_lb*sys->H_lb;
  } else {
    sys->A_lb = sys->W_lb;
  }
  sys->V_lb = sys->A_lb*sys->L_lb;
}

// -------------------------------------------------------
//
//
// -------------------------------------------------------
void get_x_limits(struct Node *node, int *xmin, int *xmax, struct System *sys)
{
  int x, y, z, nn = 0;
  int nx = sys->max_n[0];
  int nxny = nx*sys->max_n[1];

  // xmin
  for (x=0; x<sys->max_n[0]; ++x) {
    for (y=0; y<sys->max_n[1]; ++y) {
      z = 0;
      do {
        nn = x + y*nx + z*nxny;
        if (!node->is_ghost_or_periodic[nn] && node->mode[nn] == FLUID) {
          goto END1;
        }
        ++z;
      } while (z<sys->max_n[2]);
    }
  }
END1:
  x = get_global_xcoord(nn, sys);
  MPI_Reduce(&x, xmin, 1, MPI_INT, MPI_MIN, 0, MPI_COMM_WORLD);
  MPI_Bcast(xmin, 1, MPI_INT, 0, MPI_COMM_WORLD);

  // xmax
  for (x=sys->max_n[0]-1; x>=0; --x) {
    for (y=0; y<sys->max_n[1]; ++y) {
      z = 0;
      do {
        nn = x + y*nx + z*nxny;
        if (!node->is_ghost_or_periodic[nn] && node->mode[nn] == FLUID) {
          goto END2;
        }
        ++z;
      } while (z<sys->max_n[2]);
    }
  }
END2:
  x = get_global_xcoord(nn, sys);
  MPI_Reduce(&x, xmax, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
  MPI_Bcast(xmax, 1, MPI_INT, 0, MPI_COMM_WORLD);
  //printf("x_min, x_max: %d, %d\n", *xmin, *xmax);
}



//--------------------------------------------
// Read geofile (geo.dat)
// Ignore SOLID, INLET, OUTLET nodes
// After completion node->mode is FLUID (0) or SOLID (4)
//--------------------------------------------
void read_geofile(struct Node *node, struct System *sys)
{
  FILE *fp;
  int Nx, Ny, Nz, x, y, z, xr, yr, zr, tmp, n_ind;
  struct Mpi *mpi = sys->mpi;
  struct Options *opt = sys->opt;
  int *lb = mpi->lower_bound;
  int *ub = mpi->upper_bound;
  int *N = sys->MAX_N;
  int geodim[3] = { N[0] - 1, N[1] - 1, N[2] - 1 };
  int in_nodes = opt->in_nodes;
  int out_nodes = opt->out_nodes;
  //int y_wall_nodes, z_wall_nodes, wall_nodes;
  int cpy;

  //y_wall_nodes = z_wall_nodes = wall_nodes = opt->wall_nodes;

  //std::ifstream geo_file("geo.dat");
  //if (!geo_file.good()) {
     // std::cout << "ERROR: Unable to open geo.dat" << std::endl;
     // MPI_Finalize();
     // exit(0);
  //}
  //open_readfile("geo.dat", geo_file);
  fp = my_fopen(sys->files["geo"].c_str(), "r"); // open file

  skip_header_line(fp, sys->mpi);

  if (opt->add_in_out_wall) {
    geodim[0] -= (in_nodes + out_nodes);
    geodim[1] -= 2 * opt->wall_nodes;
    if (sys->n_D > 2)
      geodim[2] -= 2 * opt->wall_nodes;
  }
  if (opt->copy_in_out) {
    geodim[0] -= 2 * opt->copy_in_out;
  }
  if (opt->flip_x) {
    in_nodes = -in_nodes;
    out_nodes = -out_nodes;
  }
  if (opt->rotate_y) {
    cpy = geodim[0];
    geodim[0] = geodim[2];
    geodim[2] = cpy;
    //z_wall_nodes = wall_nodes - opt->in_nodes - opt->out_nodes;
  }
  if (opt->rotate_z) {
    cpy = geodim[0];
    geodim[0] = geodim[1];
    geodim[1] = cpy;
    //y_wall_nodes = wall_nodes - opt->in_nodes - opt->out_nodes;
  }

  geodim[2] = std::max(geodim[2], 0);

  // read geo.dat, treat INLET/OUTLET as FLUID and WALL as SOLID
  // disregard the outer periodic rim
  z = 1;
  do {  // use do-while to ensure at least one execution (for 2D runs)
    for (y = 1; y < geodim[1]; ++y) {
      for (x = 1; x < geodim[0]; ++x) {

        if (fscanf(fp, "%1d", &tmp) == EOF) {
          eof_error(sys->mpi);
        }

        xr = x; yr = y; zr = z;

        if (opt->rotate_y) {
          xr = z;
          //zr = N[0]-1-x;
          zr = geodim[0] - x;
        }
        if (opt->rotate_z) {
          xr = y;
          yr = geodim[0] - x;
        }

        Nx = xr; Ny = yr; Nz = zr;

        if (opt->flip_x)
          Nx = N[0] - 1 - xr;

        if (opt->add_in_out_wall) {
          Nx += in_nodes;
          //Ny += y_wall_nodes;
          //Nz += z_wall_nodes;
          Ny += opt->wall_nodes;
          if (sys->n_D > 2)
            Nz += opt->wall_nodes;
        }

        if (opt->copy_in_out) {
          Nx += opt->copy_in_out;
        }

        if (sys->n_D > 2) {
          if (Nz<lb[2] || Nz>ub[2] ||
              Ny<lb[1] || Ny>ub[1] ||
              Nx<lb[0] || Nx>ub[0])
            continue;
        }

        n_ind = get_local_ind(Nx, Ny, Nz, sys);

        // if tmpR == 1 then set node->Rmask[n_ind] = 1;
        // else tmpR == 0 the set node->Rmask[n_ind] = 0;


        /* Default value for inert if off (0) */
        node->inert[n_ind] = 0;
        switch (tmp)
        {
        case FLUID:
          node->mode[n_ind] = FLUID;
          break;
        case WALL:
          node->mode[n_ind] = SOLID;
          break;
        case INLET:
          node->mode[n_ind] = FLUID;
          break;
        case OUTLET:
          node->mode[n_ind] = FLUID;
          break;
        case SOLID:
          node->mode[n_ind] = SOLID ;
          break;
        case SOLID_INERT:
          node->mode[n_ind] = SOLID;
          node->inert[n_ind] = 1;
          break;
        default:
          printf("ERROR in reading geometry\n");
          exit(1);
          break;
        }
      }
      if (fscanf(fp, "\n") == EOF) {
        eof_error(sys->mpi);
      }
    }
    if (fscanf(fp, "\n") == EOF) {
      eof_error(sys->mpi);
    }
    ++z;
  } while (z < geodim[2]);

  fclose(fp);
  // update ghost-nodes and set periodic nodes (for 3D runs)
  communicate_node(node, sys);
}


//--------------------------------------------
// Allocate memory and initialize lists for inlet, outlet, and wall nodes
//--------------------------------------------
void init_node_lists(struct Node *node, struct System *sys)
{
  int *n = sys->max_n;
  //int cn, gn;

  // reset ghost nodes (otherwise this routine will not work as expected)
  // (maybe check for ghost nodes in the loop below?)
  /* for( cn = 0; cn < node->nr_ghost; ++cn ) { */
  /*   gn = node->list_ghost[cn]; */
  /*   node->mode[gn] = GHOST; */
  /* } */

  node->max_nr_wall = 5 * (2 * n[0] * n[1] + 2 * n[0] * n[2] + 2 * n[1] * n[2]);

  // allocate memory
  if (!(node->list_wall = (int *)calloc(node->max_nr_wall, sizeof(int)))) {
    printf("error allocating memory for 'node->list_wall'\n");
    exit(1);
  }

  // moved to set_inlet_outlet_nodes
  /* if ( !( node->list_in = (int *) calloc(node->nr_in, sizeof(int)) )  ) { */
  /*   printf("error allocating memory for 'node->list_in'\n"); */
  /*   exit(1); */
  /* } */

  /* if ( !( node->list_out = (int *) calloc(node->nr_out, sizeof(int)) )  ) { */
  /*   printf("error allocating memory for 'node->list_out'\n"); */
  /*   exit(1); */
  /* } */

  /* // Initiate the lists of wall, inlet and outlet nodes */
  /* cw = ci = co = 0; */
  /* for( cn = 0; cn < sys->max_n_tot; ++cn ) { */
  /*   if ( node->mode[cn] == WALL ) */
  /*     node->list_wall[cw++] = cn; */
  /*   else if ( node->mode[cn] == INLET ) */
  /*     node->list_in[ci++] = cn; */
  /*   else if ( node->mode[cn] == OUTLET ) */
  /*     node->list_out[co++] = cn; */
  /* } */
}


//--------------------------------------------
// Check if size of geofile (geo.dat) matches dimensions
// given in header of geofile
//--------------------------------------------
void check_size_of_file(struct System *sys)
{
  if (sys->mpi->my_rank == 0) {
    FILE *fp;
    int *N = sys->MAX_N;
    int Nx, Ny, Nz;
    unsigned long size_file;
    int diff;

    fp = my_fopen(sys->files["geo"].c_str(), "rb");
    skip_header_line(fp, sys->mpi);
    size_t start = ftell(fp); // get current position
    fseek(fp, 0L, SEEK_END);  // go to end of file
    size_t size_read = ftell(fp) - start;
    //std::cout << "start : " << start << ", end : " << ftell(fp) << std::endl;
    fclose(fp);

    Nx = N[0] - 2;
    Ny = N[1] - 2;
    Nz = N[2] - 2;

    update_dim_for_extra_walls(-1, Nx, Ny, Nz, sys);

    update_dim_for_rotations(Nx, Ny, Nz, sys);

	int end_chars = 1;
#ifdef _WIN32
	end_chars = 2;
#endif

    size_file = Nz*((Nx + end_chars)*Ny + 1) - 1; // - 1 to avoid newline at the end of file
    if (sys->n_D < 3)
      size_file = (Nx + end_chars)*Ny;

    diff = size_read - size_file;

    if ((diff < 0) || (diff > 2)) {
      printf("******\nERROR! Wrong dimensions for %s: [%d, %d, %d] (size_file = %lu, size_read = %zu)\n******\n", 
		  sys->files["geo"].c_str(), Nx, Ny, Nz, size_file, size_read);
      if (diff < 2) {
        printf("HINT! Try to remove newline at the end of file!\n******\n");
      }
      MPI_Finalize();
      exit(1);
    }
  }

}

//--------------------------------------------
//  Set solid-fluid neighbors to WALL nodes and update counters
//--------------------------------------------
void set_wall_nodes(struct Node *node, struct System *sys)
{
  int cn, ck, nb;
  struct Step *step = sys->step;
  int nw = 0;

  // set all WALL nodes (if any) to SOLID
  for (cn = 0; cn < sys->max_n_tot; ++cn) {
    if (node->is_ghost_or_periodic[cn])
      continue;
    if (node->mode[cn] == WALL) {
      node->mode[cn] = SOLID;
    }
  }

  // set WALL nodes
  for (cn = 0; cn < sys->max_n_tot; ++cn) {
    if (node->is_ghost_or_periodic[cn])
      continue;
    if (node->mode[cn] == SOLID) {
      for (ck = 0; ck < sys->n_Q; ++ck) {  // loop directions
        nb = cn - step->n_rel_neig[ck];
        if (node->mode[nb] == FLUID ||
          node->mode[nb] == FLUID + PERIODIC) {  // check neighbor node
          node->mode[cn] = WALL;
          node->list_wall[nw++] = cn;
          break;
        }
      }
    }
  }
  node->nr_wall = nw;
  //sys->step->wall_comp = node->nr_wall;

  if (node->nr_wall > node->max_nr_wall) {
    printf("ERROR! Number of wall nodes exceeds list length!");
    exit(1);
  }

  communicate_node(node, sys);
}


//--------------------------------------------
//  Set INLET, OUTLET nodes and update counters
//--------------------------------------------
void set_inlet_outlet_nodes(struct Node *node, struct System *sys)
{
  //struct Node *node = sys->node;
  int cn, ci, co; //, ny, nz, yz;

  ci = co = 0;

  // count inlet/outlet nodes
  for (cn = 0; cn < sys->max_n_tot; ++cn) {
    if (node->is_ghost_or_periodic[cn])
      continue;
    // INLET
    //if (get_global_xcoord(cn, sys) == 1) {
    if (get_global_xcoord(cn, sys) == node->in_x) {
      if (node->mode[cn] == FLUID) {
        ci++;
      }
    }
    // OUTLET
    //if (get_global_xcoord(cn, sys) == sys->MAX_N[0]-2) {
    if (get_global_xcoord(cn, sys) == node->out_x) {
      if (node->mode[cn] == FLUID) {
        co++;
      }
    }
  }
  node->nr_in = ci;
  node->nr_out = co;

  // allocate
  if (!(node->list_in = (int *)calloc(node->nr_in, sizeof(int)))) {
    printf("error allocating memory for 'node->list_in'\n");
    exit(1);
  }

  if (!(node->list_out = (int *)calloc(node->nr_out, sizeof(int)))) {
    printf("error allocating memory for 'node->list_out'\n");
    exit(1);
  }

  // update lists
  ci = co = 0;
  for (cn = 0; cn < sys->max_n_tot; ++cn) {
    if (node->is_ghost_or_periodic[cn])
      continue;
    // set INLET
    if (get_global_xcoord(cn, sys) == node->in_x) {
      if (node->mode[cn] == FLUID) {
        node->mode[cn] = INLET;
        node->list_in[ci++] = cn;
      }
    }
    // set OUTLET
    if (get_global_xcoord(cn, sys) == node->out_x) {
      if (node->mode[cn] == FLUID) {
        node->mode[cn] = OUTLET;
        node->list_out[co++] = cn;
      }
    }
#ifdef _FLUID_BDRY_PRESS_
    if (get_global_xcoord(cn, sys) == node->in_x-1) {
      if (node->mode[cn] == FLUID) {
        node->mode[cn] = DISCARD;
      }
    }
    if (get_global_xcoord(cn, sys) == node->out_x+1) {
      if (node->mode[cn] == FLUID) {
        node->mode[cn] = DISCARD;
      }
    }
#endif

  }

  communicate_node(node, sys);
}

//--------------------------------------------
//  Add outer solid walls
//--------------------------------------------
void add_outer_walls(struct Node *node, struct System *sys)
{
  int cn, Y, Z;
  int Ylim[2] = { sys->opt->wall_nodes + 1, sys->MAX_N[1] - 2 - sys->opt->wall_nodes };
  int Zlim[2] = { sys->opt->wall_nodes + 1, sys->MAX_N[2] - 2 - sys->opt->wall_nodes };
#ifdef _2D_RUN_
  Zlim[0] = Zlim[1] = 0;
#endif
  for (cn = 0; cn < sys->max_n_tot; ++cn) {
    if (node->is_ghost_or_periodic[cn])
      continue;
    // make wall parallel to flow
    Y = get_global_ycoord(cn, sys);
    Z = get_global_zcoord(cn, sys);
    if (Y<Ylim[0] || Z<Zlim[0] || Y>Ylim[1] || Z>Zlim[1]) {
      node->mode[cn] = SOLID;
    }
  }
  communicate_node(node, sys);
}

//--------------------------------------------
//  Copy slice at x=1 to x=0, and slice at x=-2 to x=-1
//--------------------------------------------
void copy_inlet_outlet_slice(struct Node *node, struct System *sys)
{
  int nx, ny, nz, i, j;
  int cpy_xpos[2] = {1+sys->opt->copy_in_out, sys->MAX_N[0]-2-sys->opt->copy_in_out};
  int cpy_xdir[2] = {-1, 1}; // 
  int *n = sys->max_n;
  int cpy_from, cpy_to;

  for (i = 0; i < 2; ++i) {
    nx = global_to_local_xind(cpy_xpos[i], sys);
    if (nx > 0) {
      //printf("<%d> nx=%d\n",sys->mpi->my_rank, nx);
      for (ny=0; ny<n[1]; ++ny) {
        for (nz=0; nz<n[2]; ++nz) {
          cpy_from = nx + ny*n[0] + nz*n[0]*n[1];
          //node->mode[cpy_from] = 5;
          for (j=1; j<sys->opt->copy_in_out+1; ++j) {
            cpy_to = cpy_from + cpy_xdir[i]*j;
            node->mode[cpy_to] = node->mode[cpy_from];
          }
        }
      }
    }
  }
  communicate_node(node, sys);
}

//--------------------------------------------
//  Copy slice at x=1 to x=0, and slice at x=-2 to x=-1
//--------------------------------------------
void add_overlapping_slice_at_inlet_and_outlet(struct Node *node, struct System *sys)
{
  int cpy_xpos[2] = {1, sys->MAX_N[0]-2};
  int periodic_node_xdir[2] = {-1, 1}; // 
  int *n = sys->max_n;
  
  for (int i = 0; i < 2; ++i) {
    int nx = global_to_local_xind(cpy_xpos[i], sys);
    if (nx > 0) {
      //printf("<%d> nx=%d\n",sys->mpi->my_rank, nx);
      for (int ny=0; ny<n[1]; ++ny) {
        for (int nz=0; nz<n[2]; ++nz) {
          int node_to_update = nx + ny*n[0] + nz*n[0]*n[1];
          //node->mode[cpy_from] = 5;
	  int periodic_node = node_to_update + periodic_node_xdir[i];
	  node->mode[node_to_update] = ((node->mode[node_to_update]==0) || (node->mode[periodic_node]==5)) ? 0 : 1;
        }
      }
    }
  }
  communicate_node(node, sys);
}


//--------------------------------------------
//  Close inlet and outlet (for diffusive run)
//--------------------------------------------
void close_inlet_outlet(struct Node *node, struct System *sys)
{
  //struct Node *node = sys->node;
  int cn, X;
  int Xlim[2] = { sys->opt->in_nodes + 1, sys->MAX_N[0] - 2 - sys->opt->out_nodes };

  for (cn = 0; cn < sys->max_n_tot; ++cn) {
    if (node->is_ghost_or_periodic[cn])
      continue;
    X = get_global_xcoord(cn, sys);
    if (X <= Xlim[0] || X >= Xlim[1]) {
      node->mode[cn] = SOLID;
    }
  }
  communicate_node(node, sys);
}

//--------------------------------------------
//  Close outlet (to find percolating cluster)
//--------------------------------------------
void close_outlet(struct Node *node, struct System *sys)
{
  //struct Node *node = sys->node;
  int cn, X;

  for (cn = 0; cn < sys->max_n_tot; ++cn) {
    if (node->is_ghost_or_periodic[cn])
      continue;
    X = get_global_xcoord(cn, sys);
    if (X > sys->MAX_N[0] - 2 - sys->opt->out_nodes) {
      node->mode[cn] = SOLID;
    }
  }
  communicate_node(node, sys);
}

//--------------------------------------------
//  Add raster at inlet and outlet
//--------------------------------------------
void add_raster(struct Node *node, struct System *sys)
{
  //struct Node *node = sys->node;
  int cn, *cy, *cz, i, y, z;
  int mid_y, mid_z, dy, dz, y0, z0, Y, Z;
  int *N = sys->MAX_N;
  struct Options *opt = sys->opt;
  int num = opt->add_raster; // number of cubes in y-, z-dir
  int n;
  int w = 6;

  int shift = 1 + sys->opt->wall_nodes;
  int Ny = N[1] - 2 * shift;
  int Nz = std::max(N[2] - 2 * shift, 0);

  cy = (int *)malloc(num*num * sizeof(int));
  cz = (int *)malloc(num*num * sizeof(int));

  // center of inlet
  mid_y = (int)(0.5*Ny);
  mid_z = (int)(0.5*Nz);
  // distance between cubes
  dy = (int)((Ny) / num);
  dz = (int)((Nz) / num);
  // 
  y0 = mid_y - dy*(int)(0.5*(num - 1));
  z0 = mid_z - dz*(int)(0.5*(num - 1));

  //if (sys->mpi->my_rank==0)
  //  printf("%d, %d, %d, %d, %d, %d, %d, %d\n", N[1], N[2], mid_y, mid_z, dy, dz, y0, z0);
  n = 0;
  for (z = 0; z < num; ++z) {
    for (y = 0; y < num; ++y) {
      cy[n] = y0 + y*dy + shift;
      cz[n] = z0 + z*dz + shift;
      //if (sys->mpi->my_rank==0) {
      //printf("cy[%d] = %d, cz[%d] = %d\n", n, cy[n], n, cz[n]);
      //fflush(stdout);
      //}
      n++;
    }
  }

  for (cn = 0; cn < sys->max_n_tot; ++cn) {
    if (node->is_ghost_or_periodic[cn])
      continue;
    if (get_global_xcoord(cn, sys) < opt->raster_nodes + 1 ||
      get_global_xcoord(cn, sys) > sys->MAX_N[0] - 2 - opt->raster_nodes) {
      node->mode[cn] = SOLID;
      Y = get_global_ycoord(cn, sys);
      Z = get_global_zcoord(cn, sys);
      for (i = 0; i < num*num; ++i) {
        if (Y > (cy[i] - w) && Y<(cy[i] + w) && Z>(cz[i] - w) && Z < (cz[i] + w)) {
          node->mode[cn] = FLUID;
        }
      }
    }
  }
  free(cy);
  free(cz);
  communicate_node(node, sys);
}

//--------------------------------------------
//--------------------------------------------
void set_dimensions_from_geofile(struct System *sys)
{
  int *n = sys->max_n;
  // read grid size and resolution from geo-file
  // check for format: int int int double
  std::ifstream infile(sys->files["geo"].c_str());
  if (!infile.good()) {
    std::cerr << "ERROR! Unable to open " << sys->files["geo"] << std::endl;
    exit(0);
  }
  std::string line;
  std::getline(infile, line);
  infile.close();
  std::queue<std::string> tokens;
  std::istringstream stream(line);
  std::string word;
  do {
    stream >> word;
    tokens.push(word);
  } while (stream);
  if ((tokens.size() - 1) != unsigned(sys->n_D + 1)) {
    printf("\nERROR! Wrong heading in %s. First line must be formatted as 'nx ny dx' (2D) or 'nx ny nz dx' (3D)\n\n",
		sys->files["geo"].c_str());
    exit(1);
  }
  for (int i = 0; i < sys->n_D; i++) {
    std::istringstream(tokens.front()) >> sys->max_n[i];
    tokens.pop();
  }
  std::istringstream(tokens.front()) >> sys->dx;

  //for (int i=0; i<3; i++) {
  //  std::cout << sys->max_n[i] << ", ";
  //}
  //std::cout << sys->dx << std::endl;

  update_dim_for_rotations(n[0], n[1], n[2], sys);

  update_dim_for_extra_walls(+1, n[0], n[1], n[2], sys);

}

//--------------------------------------------
//--------------------------------------------
void update_dim_for_rotations(int &nx, int &ny, int &nz, struct System *sys)
{
  double tmp;
  if (sys->opt->rotate_y) {
    tmp = nx;
    nx = nz;
    nz = tmp;
  }
  if (sys->opt->rotate_z) {
    tmp = nx;
    nx = ny;
    ny = tmp;
  }
}

//--------------------------------------------
//--------------------------------------------
void update_dim_for_extra_walls(int pos_neg, int &nx, int &ny, int &nz, struct System *sys)
{
  if (sys->opt->add_in_out_wall) {
    nx += pos_neg*(sys->opt->in_nodes + sys->opt->out_nodes);
    ny += pos_neg * 2 * sys->opt->wall_nodes;
    if (nz > 0) {
      nz += pos_neg * 2 * sys->opt->wall_nodes;
    }
  }
  if (sys->opt->copy_in_out) {
    nx += pos_neg * 2 * sys->opt->copy_in_out;
  }
}


//--------------------------------------------
//--------------------------------------------
/*void remove_raster(struct System *sys)
{
  struct Node *node = sys->node;
  struct Options *opt = sys->opt;
  int *N = sys->MAX_N;
  int cn, X, Y, Z;
  int Xlim[2] = {opt->raster_nodes+1, N[0]-2-opt->raster_nodes};
  int Ylim[2] = {opt->wall_nodes+1  , N[1]-2-opt->wall_nodes};
  int Zlim[2] = {opt->wall_nodes+1  , N[2]-2-opt->wall_nodes};

  for( cn = 0; cn < sys->max_n_tot; ++cn ) {
    if (node->is_ghost[cn])
      continue;
    X = get_global_xcoord(cn, sys);
    Y = get_global_ycoord(cn, sys);
    Z = get_global_zcoord(cn, sys);
    // set all raster nodes to FLUID
    if (X<Xlim[0] || X>Xlim[2]) {
      node->mode[cn] = FLUID;
    }
    if (Y<Ylim[0] || Y>Ylim[1] || Z<Zlim[0] || Z>Zlim[1]) {
      node->mode[cn] = SOLID;
    }

  communicate_node(sys);
  set_wall_nodes(sys);
  set_inlet_outlet_nodes(sys);

  // update ghost-nodes
  communicate_node(sys);

  update_linklists(sys);


  }
} */

//--------------------------------------------
//  set periodic nodes (only called for 2D runs)
//--------------------------------------------
void set_periodic_nodes_2D(struct Node *node, struct System *sys)
{
  int cn, pn, sn;
  for (cn = 0; cn < sys->node->nr_per; ++cn) {
    pn = sys->node->list_per[2 * cn];
    sn = sys->node->list_per[2 * cn + 1];
    sys->node->mode[pn] = sys->node->mode[sn] + PERIODIC;
  }
}

//--------------------------------------------
//  skip header in geo.dat
//--------------------------------------------
void skip_header_line(FILE *fp, struct Mpi *mpi)
{
  char buffer[100];
  // skip header line
  if (fgets(buffer, sizeof(buffer), fp) == NULL) {
    if (mpi->my_rank == 0) {
      perror("Error reading geo-file\n");
      MPI_Finalize();
      exit(0);
    }
  }

}

//--------------------------------------------
//  Error output
//--------------------------------------------
void eof_error(struct Mpi *mpi) {
  if (mpi->my_rank == 0)
    printf("ERROR in read_geofile: End-of-file reached before reading finished!\n");
  MPI_Finalize();
  exit(0);
}
