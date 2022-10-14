/*
 * moving_bondaries.c
 *
 *  Created on: 1. june 2012
 *      Author: janlv
 *
 */
#include "global.h"
#include "chem_global.h"
#include "output.h"
#include "biozementfunctions.h"


 /* * * * * * * * * * * * *
  *                       *
  * * * * * * * * * * * * */
int init_linklist(struct Node *node, struct System *sys)
{
  //struct Node *node = sys->node;
  struct Links *links = node->links;
  struct Step *step = sys->step;
  int nlinks, nw, ck, nb, nn, nbm;

  // find number of links
  nlinks = 0;
  for (nw = 0; nw < node->nr_wall; ++nw) {  /* loop over wall-nodes */
    nn = node->list_wall[nw];
    for (ck = 0; ck < sys->n_Q; ++ck) { /* direction */
      nb = nn - step->n_rel_neig[ck];
      //if (node->is_ghost[nb])
        //continue;
      nbm = node->mode[nb];
      if (nbm != SOLID && nbm != SOLID + PERIODIC &&
        nbm != WALL  && nbm != WALL + PERIODIC) {
        nlinks++;
      }
    }
  }

  if (nlinks < 1000)
    nlinks = 1000;

  links->max_links = 8 * nlinks;
  links->max_groups = 8 * nlinks;

  //if (!(links->list = (struct Links::Link *)  calloc(links->max_links, sizeof(struct Links::Link)))) return BAD;
  //if (!(links->group = (struct Links::Group *) calloc(links->max_groups, sizeof(struct Links::Group)))) return BAD;
  links->list = new Links::Link[links->max_links]();
  links->group = new Links::Group[links->max_groups]();

  communicate_node(node, sys);
  update_linklist(node, sys);

#ifdef _DEBUG_MPI_
  printf("LINKLIST <%03d> ngroups: %d/%d, nlinks: %d/%d\n",
    sys->mpi->my_rank, node->links->ngroups, node->links->max_groups,
    node->links->nlinks, node->links->max_links);
#endif // _DEBUG_MPI_

  return GOOD;
}


/* * * * * * * * * * * * *
 *                       *
 * * * * * * * * * * * * */
void update_linklist(struct Node *node, struct System *sys)
{

  /* set up list over wall-fluid node links */
  int nn, ck, nb, nbm;
  //struct Node *node = sys->node;
  struct Step *step = sys->step;
  struct Links *links = node->links;

  int *N = sys->MAX_N;
  struct Mpi::Ghost ghost;
  struct Mpi *mpi = sys->mpi;
  int nck, nlg, gnx;
  int ln, lg, nw, first_link;/* , ghost_links = 0; */
  int inx = sys->opt->inertx_inlet;
  int outx = sys->opt->inertx_outlet;
  //int xcoord, ycoord, zcoord;

  ln = 0; /* link number */
  lg = 0; /* link group number */

  /* loop over wall-nodes */
  for (nw = 0; nw < node->nr_wall; ++nw) {
    nn = node->list_wall[nw];

    first_link = 1;
    for (ck = 0; ck < sys->n_Q; ++ck) { /* direction */
      nb = nn - step->n_rel_neig[ck];
      nbm = node->mode[nb];
      if (nbm != SOLID && nbm != SOLID + PERIODIC &&
        nbm != WALL  && nbm != WALL + PERIODIC) {

        /* new link */
        if (first_link) {
          /* make new group */
          links->group[lg].wall = nn;
          links->group[lg].first = ln;
          first_link = 0;
        }
        links->list[ln].group = lg;
        links->list[ln].wall = nn;
        links->list[ln].fluid = nb;
        links->list[ln].dir = ck;
        links->list[ln].ghost = 0;

        // set nodes around inlet and outlet inert (opt->inertx_inlet/outlet)
        gnx = get_global_xcoord(nn, sys);
        if ((gnx < inx) || (gnx > N[0] - outx)) {
          links->list[ln].inert = 1;
        }
        else {
          links->list[ln].inert = 0;
        }

        // check if wall links should be set inert (default is reactive)
        int xcoord, ycoord, zcoord;
        switch (sys->opt->inert_walls) {
        case (0):
          // all walls are reactive, do nothing
          break;
        case (1):
          // set outer walls inert
          ycoord = get_global_ycoord(nn, sys);
          zcoord = get_global_zcoord(nn, sys);
#ifdef _2D_RUN_
          if (ycoord < sys->opt->wall_nodes + 1 ||
            ycoord > N[1] - 2 - sys->opt->wall_nodes) {
            links->list[ln].inert = 1;
          }
#else
          if (ycoord < sys->opt->wall_nodes + 1 ||
            zcoord < sys->opt->wall_nodes + 1 ||
            ycoord > N[1] - 2 - sys->opt->wall_nodes ||
            zcoord > N[2] - 2 - sys->opt->wall_nodes) {
            links->list[ln].inert = 1;
          }
#endif
          // ... and also inlet/outlet walls if the cell is closed
          if (sys->opt->closed_cell) {
            xcoord = get_global_xcoord(nn, sys);
            if (xcoord < sys->opt->in_nodes + 1 ||
              xcoord > sys->MAX_N[0] - 2 - sys->opt->out_nodes) {
              links->list[ln].inert = 1;
            }
          }
          break;
        case (2):
          // set everything except the bottom plate inert
          if (get_global_ycoord(nn, sys) > sys->opt->wall_nodes) {
            links->list[ln].inert = 1;
          }
          break;
        }
        ln++;
      } // end of new link test 
    } // end of directions loop
    if (first_link == 0) {
      links->group[lg].last = ln;
      lg++;
    }
  } // end of wall-node loop

  /* loop over ghost-nodes */
  nlg = 0;
  for (nw = 0; nw < node->nr_ghost; ++nw) {
    nn = node->list_ghost[nw];

    if (node->mode[nn] == WALL) {
      ghost = mpi->ghost[nw];
      first_link = 1;
      for (nck = 0; nck < ghost.ndir; ++nck) { /* direction */
        ck = ghost.dir[nck];
        nb = nn - step->n_rel_neig[ck];
        nbm = node->mode[nb];
        if (nbm != SOLID && nbm != SOLID + PERIODIC &&
          nbm != WALL  && nbm != WALL + PERIODIC) {

          /* new link */
          if (first_link) {
            /* make new group */
            links->group[lg].wall = nn;
            links->group[lg].first = ln;
            first_link = 0;
          }
          links->list[ln].group = lg;
          links->list[ln].wall = nn;
          links->list[ln].fluid = nb;
          links->list[ln].dir = ck;
          links->list[ln].ghost = 1;
          gnx = get_global_xcoord(nn, sys);
          if ((gnx < inx) || (gnx > N[0] - outx)) {
            links->list[ln].inert = 1;
          }
          else {
            links->list[ln].inert = 0;
          }
          ln++;
          nlg++;
        }
      }
      if (first_link == 0) {
        links->group[lg].last = ln;
        lg++;
      }
    }
  }

  links->nlinks = ln;
  links->ngroups = lg;

  if (sys->opt->find_isolated_nodes && sys->opt->set_isolated_inert) {
    set_isolated_nodes_inert(node);
  }

  /* Biozement added inert solids */
    set_nodes_inert(node);

  if (links->ngroups > links->max_groups) {
    printf("<%03d> ERROR! Max number of groups exceeded: %d/%d !\n",
      mpi->my_rank, links->ngroups, links->max_groups);
    MPI_Finalize();
    exit(1);
  }
  if (links->nlinks > links->max_links) {
    printf("<%03d> ERROR! Max number of links exceeded: %d/%d !\n",
      mpi->my_rank, links->nlinks, links->max_links);
    MPI_Finalize();
    exit(1);
  }

  return;
}



/* * * * * * * * * * * * *
 *                       *
 * * * * * * * * * * * * */
void update_solid_v5(struct System *sys, struct Minerals *minerals, struct Links *links)
{
  /*
     dgi = dgi_out - dgi_in
     if dgi > 0: precipitation, dS_w = sum(dgi) put in delta wall-node
     if dgi < 0: dissolution,   dS_f = sum(dgi) put in delta fluid-node
     S(t+1) = S(t) + dS
     if S>1 => b_w = (S - 1)/dS_w
     if S<1 => b_f = S/dS_f
     -Loop over links
     -Distribute positive sfrac (b_w) to FLUID nodes via links and
      weighted by dgi for the i'th link
     -Distribute negative sfrac (b_f) to WALL nodes
   */

  real change, sum; /* tot_fluid; */
  int min, nl, nn;
  int i;
  struct Minerals::Mineral *mineral;
  struct Links::Link *link;
  struct Node *node = sys->node;

  real W, surplus;

  // zero out delta (dS) for all minerals
  for (nn = 0; nn < sys->max_n_tot; ++nn) {
    for (min = 0; min < minerals->n_tot; ++min) {
      minerals->list[min].delta[nn] = 0.;
    }
  }

  /* outermost loop over minerals; one mineral at a time */
  for (min = 0; min < minerals->n_tot; ++min) {
    mineral = &minerals->list[min];

    //printf("<%d> nlinks: %d\n",mpi->my_rank, links->nlinks);
    for (nl = 0; nl < links->nlinks; ++nl) { /* loop links */
      link = &links->list[nl];

      if (link->inert)
        continue;

      change = link->delta[min];

      if (change > 0.) {
        /*  PRECIPITATION: add to WALL node */
        mineral->delta[link->wall] += change;
      }
      if (change < 0.) {
        /*  DISSOLUTION: remove from FLUID node */
        mineral->delta[link->fluid] += change;
      }
    }


    /* add delta to sfrac and put surplus or *
     * deficit sfrac/delta in delta         */
    for (nn = 0; nn < sys->max_n_tot; ++nn) {  /* loop nodes */

      if ((mineral->delta[nn] == 0.) || node->is_ghost_or_periodic[nn])
        continue;

      if ((node->mode[nn] == WALL) || (node->mode[nn] == FLUID)) {
        // Add mineral change to node nn
        mineral->sfrac[nn] += mineral->delta[nn]; /* Si(t+1) = Si(t) + dSi */
        // JLV
        // To add mass to the system, like in a plug-run where a real plug
        // is much heavier than the numerical system, the mass-flux can be reduced
        // by uncommenting the following line (and commenting the above line)
        //mineral->sfrac[nn] += sys->plug_mass_shift*mineral->delta[nn]; /* Si(t+1) = Si(t) + dSi */
        // JLV

        if (node->mode[nn] == WALL) {
          /* check for sfrac_total > 1 */
          sum = 0.;
          for (i = 0; i < minerals->n_tot; ++i)               /* S = sum_i Si */
            sum += minerals->list[i].sfrac[nn];
          if (sum > 1.0) {
            surplus = sum - 1.0;
            //mineral->sfrac[nn] = mineral->sfrac[nn]-sum+1.;                      /* remove surplus and */
            //mineral->delta[nn] = sum/mineral->delta[nn] - 1./mineral->delta[nn]; /* put it in dS */
            mineral->sfrac[nn] -= surplus;                   /* remove surplus and */
            mineral->delta[nn] = surplus / mineral->delta[nn]; /* put it in dS */
          }
          else {
            mineral->delta[nn] = 0.;
          }
        }
        else if (node->mode[nn] == FLUID) {
          /* check for sfrac < 0. */
          if (mineral->sfrac[nn] < 0.) {
            mineral->delta[nn] = fabs(mineral->sfrac[nn] / mineral->delta[nn]); // NB!! delta<0 for FLUID
            mineral->sfrac[nn] = 0.;
          }
          else {
            mineral->delta[nn] = 0.;
          }
        }
      }
    } /* end node-loop */

    /* communicate both sfrac and the remainder in delta   *
     * which is normalized by the sum of dgi, dgi = gi_out - gi_in */
     // I think it is unnecessary to communicate sfrac here...
     //communicate_mineral_sfrac(min, sys, minerals); // or mineral_sfrac since we are inside mineral-loop...
    communicate_mineral_delta(min, sys, minerals);

    for (nl = 0; nl < links->nlinks; ++nl) { /* loop links */
      link = &links->list[nl];

      if (link->inert)
        continue;

      W = link->delta[min];

      if (W > 0.) {
        /* PRECIPITATION: sfrac surplus in wall added to fluid node */
        mineral->sfrac[link->fluid] += mineral->delta[link->wall] * W;
      }

      if (W < 0.) {
        /* DISSOLUTION: sfrac deficit in fluid subtracted from wall node */
        mineral->sfrac[link->wall] += mineral->delta[link->fluid] * W;
        if (mineral->sfrac[link->wall] < 0.) {
          //mineral->rest[link->wall] += mineral->sfrac[link->wall];
          mineral->sfrac[link->wall] = 0.;
        }
      }
    } /* end link-loop */


    communicate_mineral_sfrac(min, sys, minerals);
  } /* end mineral-loop */

}



//
//
//
void check_convergence_v2(struct System *sys, struct Minerals *minerals,
  struct Field *species, struct Boundary *bndry, struct Links *links)
{
  int nn, min;
  struct Minerals::Mineral *mineral;
  struct Node *node = sys->node;
  double sfrac;
  double tstep = -1.0, min_tstep = 1.0e25, jump_time;
#ifdef KEEP_INTERVAL
  double next_write;
#endif
  double dm;
  int time_since_conv = sys->t - sys->last_convergence;
  int not_converged;
  double diff = 0.0;
#ifdef _DEBUG_
  int n0 = sys->max_n[0];
  int n1 = sys->max_n[1];
  int n01 = sys->max_n[0] * sys->max_n[1];
  int type = -1, site = 0, mm = 0;
  double dms = 0, sf = 0, sfs = 0, dmss = 0;
#endif
  double offset = sys->conv_jump;
  double conv_crit;
  int nl;
  struct Links::Link *link;
  //  static int first_convergence = 1;

  // if (first_convergence) {
  //   conv_crit = sys->opt->init_conv_lim*sys->conv_length;
  // } else {
  conv_crit = sys->conv_crit_chem*sys->conv_length;
  //  }

  if (time_since_conv%sys->conv_length == 0) {
    not_converged = 1;

    if (time_since_conv > sys->steps_until_check) {
#ifdef _MEAN_CONVERGENCE_
      check_mean_convergence(&not_converged, sys, species, conv_crit, &diff);
#else
      check_convergence(&not_converged, sys, species, conv_crit, &diff);
#endif

#ifdef _DEBUG_CONVERGENCE_
      // screendump steps and difference every 100 steps
      int proc = get_global_max(sys->mpi, &diff);
      MPI_Bcast(&(sys->rho_not_converged), 1, MPI_DOUBLE, proc, MPI_COMM_WORLD);
      MPI_Bcast(&(sys->node_not_converged), 1, MPI_INT, proc, MPI_COMM_WORLD);
      int ind[3];
      get_global_coord(sys->node_not_converged, ind, sys);
      if (sys->mpi->my_rank == 0) {
        char fn[50];
        sprintf(fn, "%s/debug_converge.dat", (*sys->files)["out"].c_str());
        FILE *fp = my_fopen(fn, "a");
        fprintf(fp, "%6d : diff = %.2e (%.2e), rho = %.8e, (%d,%d,%d)\n",
          sys->t, diff, conv_crit, sys->rho_not_converged, ind[0], ind[1], ind[2]);
        //fflush(stdout);
        fclose(fp);
      }
#endif
    }

    if (not_converged) {
#ifdef _MEAN_CONVERGENCE_
      // do nothing
#else
      update_rho_old(sys, species);
#endif
    }
    else {
      // CONVERGENCE !!
      sys->convergence = 1;
      sys->last_convergence = sys->t;
      // if (first_convergence) {
      // 	first_convergence = 0;
      // }

      // save current mineral values for calculation  
      // of the change in mineral mass next convergence step
      for (min = 0; min < minerals->n_tot; ++min) { /* loop minerals */
        mineral = &minerals->list[min];
        for (nn = 0; nn < sys->max_n_tot; ++nn) {  /* loop nodes */
          mineral->old_sfrac[nn] = mineral->sfrac[nn];
        }
        // reset delta_sum
        for (nl = 0; nl < links->nlinks; ++nl) { /* loop links */
          link = &links->list[nl];
          link->delta_sum[min] = 0.0;
        }
      }
    } // end if (not_converged)

    if (sys->convergence && (time_since_conv == sys->conv_length)) {
      // calculate new dm and update sfrac (volume fraction)
      sys->convergence = 0;

      // find convergence dm
      for (min = 0; min < minerals->n_tot; ++min) { /* loop minerals */
        mineral = &minerals->list[min];
        for (nn = 0; nn < sys->max_n_tot; ++nn) {  /* loop nodes */
          mineral->dm[nn] = mineral->sfrac[nn] - mineral->old_sfrac[nn];
          //#ifdef _MAGN_ON_MAGN_
          //      	  if ((min==1) && (mineral->sfrac[nn]<sys->opt->magn_limit))
          //      	    mineral->dm[nn] = 0.0;
          //#endif	  
        }
      }

      // calculate next timestep
      for (nn = 0; nn < sys->max_n_tot; ++nn) {  /* loop nodes */

        if (node->is_ghost_or_periodic[nn])
          continue;

        minerals->sfrac_sum[nn] = 0.;
        minerals->dm_sum[nn] = 0.;
        for (min = 0; min < minerals->n_tot; ++min) { /* loop minerals */
          minerals->dm_sum[nn] += minerals->list[min].dm[nn];
          minerals->sfrac_sum[nn] += minerals->list[min].sfrac[nn];
        }

        if (minerals->dm_sum[nn] > 0.) {
          // more precipitation (than dissolution)
          if (minerals->sfrac_sum[nn] < 0.5) {
            // when will this FLUID node be SOLID (>=0.5)?
            //tstep = (0.5001 - minerals->sfrac_sum[nn])/minerals->dm_sum[nn];
            tstep = (0.5 + offset - minerals->sfrac_sum[nn]) / minerals->dm_sum[nn];
            if (tstep < min_tstep) {
              min_tstep = tstep;
#ifdef _DEBUG_
              type = 0;
              dms = -1.0;
              dmss = minerals->dm_sum[nn];
              site = nn;
              mm = tstep;
              sf = -1.0;
              sfs = minerals->sfrac_sum[nn];
#endif // DEBUG
            }
          }
          if (minerals->sfrac_sum[nn] >= 0.5) {
            // when will this SOLID node be full?
            tstep = (1.0 + offset - minerals->sfrac_sum[nn]) / minerals->dm_sum[nn];
            if (tstep < min_tstep) {
              min_tstep = tstep;
#ifdef _DEBUG_
              type = 1;
              dms = -1.0;
              dmss = minerals->dm_sum[nn];
              site = nn;
              mm = tstep;
              sf = -1.0;
              sfs = minerals->sfrac_sum[nn];
#endif // DEBUG
            }
          }
        }
        if (minerals->dm_sum[nn] < 0.) {
          // more dissolution (than precipitation)
          if (minerals->sfrac_sum[nn] >= 0.5) {
            // when will this SOLID node be FLUID?
            tstep = (minerals->sfrac_sum[nn] - (0.5 - offset)) / fabs(minerals->dm_sum[nn]);
            if (tstep < min_tstep) {
              min_tstep = tstep;
#ifdef _DEBUG_
              type = 2;
              dms = -1.0;
              dmss = minerals->dm_sum[nn];
              site = nn;
              mm = tstep;
              sf = -1.0;
              sfs = minerals->sfrac_sum[nn];
#endif // DEBUG
            }
          }
        }

        // check FLUID and WALL nodes for when they empty of minerals
        for (min = 0; min < minerals->n_tot; ++min) { /* loop minerals */
          mineral = &minerals->list[min];
          dm = mineral->dm[nn];
          sfrac = mineral->sfrac[nn];
          if (dm < 0.) {
            // dissolution
                //if (node->mode[nn] == FLUID || node->mode[nn] == WALL) {
            //if (minerals->sfrac_sum[nn] >= 0.5) {
            if (node->mode[nn] == WALL) {
              // when will WALL empty of this mineral?
              // don't care about FLUID nodes since only mineral
              // conc. in wall is passed to chem-solver
              tstep = (sfrac + offset) / fabs(dm);
              if (tstep < min_tstep) {
                min_tstep = tstep;
#ifdef _DEBUG_
                type = 3;
                dms = dm;
                dmss = minerals->dm_sum[nn];
                site = nn;
                mm = min;
                sf = sfrac;
                sfs = minerals->sfrac_sum[nn];
#endif // DEBUG
              }
            }
          }
        } /* end minerals-loop */
      } /* end nodes-loop */


#ifdef _DEBUG_
      sprintf(fn, "%s/%04d_converge.dat", (*sys->files)["out"].c_str(), sys->mpi->my_rank);
      fp = my_fopen(fn, "a");
      fprintf(fp, "tmin = %5.4e, type: %d, dm: % 5.4e, dm_sum: % 5.4e, site: %2d,%2d,%2d, min: %d, sf: %5.4e, sfs: %5.4e)\n",
        min_tstep, type, dms, dmss, site%n0, ((site - site%n0) / n0) % n1, (int)site / n01, mm, sf, sfs);
      fclose(fp);
#endif // DEBUG

      broadcast_global_min(sys->mpi, &min_tstep);

      jump_time = min_tstep*sys->conv_length*sys->dt;

      // check if jump will skip next reporting timestep
#ifdef KEEP_INTERVAL
      next_write = sys->nwrite*sys->write_int;
      if ((sys->time + jump_time > next_write) && (next_write > sys->time + sys->dt)) {
        jump_time = next_write - sys->time;
      }

#endif
      sys->time += jump_time;
      sys->t_extra += jump_time / sys->dt;


#ifdef _DEBUG_CONV_
      if (sys->mpi->my_rank == 0) {
        printf("*** Advance %4.3e timesteps at t = %7d ***\n", jump_time / sys->dt, sys->t);
      }
#endif // _MPI_ 

      // update sfrac (volume fraction) with convergence dm
      //      for( nn = 0; nn < sys->max_n_tot; ++nn ) {  /* loop nodes */
      //        minerals->sfrac_sum[nn] = 0.;
      //        for (min = 0; min < minerals->n_tot; ++min) { /* loop minerals */
      //          mineral = &minerals->list[min];
      //          mineral->sfrac[nn] += (min_tstep*mineral->dm[nn]);
      //        }
      //      }

      // set link->delta
      for (min = 0; min < minerals->n_tot; ++min) { /* loop minerals */
        for (nl = 0; nl < links->nlinks; ++nl) { /* loop links */
          link = &links->list[nl];
          link->delta[min] = link->delta_sum[min] * min_tstep;
          link->delta_sum[min] = 0.0;
        }
      }

      update_solid_v5(sys, minerals, links);

      update_nodes(sys, bndry, minerals);
      update_linklist(node, sys);

    }
  }

}



/* * * * * * * * * * * * *
 *  Update nodes         *
 * * * * * * * * * * * * */
void update_nodes(struct System *sys, struct Boundary *bndry, struct Minerals *minerals)
{
  struct Node *node = sys->node;
  struct Step *step = sys->step;
  struct Mpi *mpi = sys->mpi;
  int nn, int_receive_buffer;
  int ck, nbn, nw;
  int f2w, w2f;

  /* update node mode:     *
   * FLUID if sfrac < 0.5  *
   * SOLID if sfrac >= 0.5 */

  sys->num_new_nodes = sys->num_new_solid_nodes = 0;
  sys->velreinit = 0;
  sys->rho_reinit = 0;

  // update nodes that have changed to/from FLUID/SOLID
  // and set velocity reiteration flag to one if they have.
  // SOLID -> FLUID nodes must also be initialized with
  // species distributions from neighboring FLUID nodes 
  for (nn = 0; nn < sys->max_n_tot; ++nn) {
    if (node->mode[nn] == FLUID ||
      node->mode[nn] == WALL) {
      has_node_changed(nn, minerals, node, &f2w, &w2f);
      if (f2w || w2f)
        sys->velreinit = 1;
      if (w2f) {
        sys->new_nodes[sys->num_new_nodes++] = nn;
      }
      if (f2w) {
        sys->new_solid_nodes[sys->num_new_solid_nodes++] = nn;
      }
    }
  }

  // add and broadcast velreinit-flag
  MPI_Reduce(&(sys->velreinit), &int_receive_buffer, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  if (mpi->my_rank == 0) sys->velreinit = int_receive_buffer;
  MPI_Bcast(&(sys->velreinit), 1, MPI_INT, 0, MPI_COMM_WORLD);

  // add and broadcast rho_reinit-flag
  MPI_Reduce(&(sys->num_new_nodes), &int_receive_buffer, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  if (mpi->my_rank == 0) sys->rho_reinit = int_receive_buffer;
  MPI_Bcast(&(sys->rho_reinit), 1, MPI_INT, 0, MPI_COMM_WORLD);
  sys->solid2fluid += sys->rho_reinit;

  // add f2w
  MPI_Reduce(&(sys->num_new_solid_nodes), &int_receive_buffer, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  if (mpi->my_rank == 0) sys->fluid2solid += int_receive_buffer;
  MPI_Bcast(&(sys->fluid2solid), 1, MPI_INT, 0, MPI_COMM_WORLD);

  // set SOLID nodes to WALL if they have FLUID neighbours
  nw = 0;
  for (nn = 0; nn < sys->max_n_tot; ++nn) {
    // wall-list should not contain ghost-nodes
    if (node->is_ghost_or_periodic[nn])
      continue;
    if (node->mode[nn] == SOLID) {
      for (ck = 0; ck < sys->n_Q; ++ck) {
        nbn = node->mode[nn - step->n_rel_neig[ck]];
        if (nbn != SOLID && nbn != SOLID + PERIODIC &&
          nbn != WALL  && nbn != WALL + PERIODIC) {
          node->mode[nn] = WALL;
          node->list_wall[nw++] = nn;
          break;
        }
      }
    }
  }
  node->nr_wall = nw;
  sys->step->wall_comp = node->nr_wall;

  if (node->nr_wall > node->max_nr_wall) {
    printf("ERROR! Number of wall nodes exceeds list length!");
    exit(1);
  }

  communicate_node(node, sys);

  return;
}

/* * * * * * * * * * * * *
 *                       *
 * * * * * * * * * * * * */
void has_node_changed(int nn, struct Minerals *minerals, struct Node *node,
  int *f2w, int *w2f)
{
  int mn;
  real sfrac_tot = 0.;

  for (mn = 0; mn < minerals->n_tot; ++mn) {
    sfrac_tot += minerals->list[mn].sfrac[nn];
  }

  *f2w = *w2f = 0;
  if (sfrac_tot >= 0.5) {
    // FLUID --> WALL
    if (node->mode[nn] == FLUID)
      *f2w = 1;
    node->mode[nn] = SOLID;
  }
  else {
    // WALL --> FLUID
    if (node->mode[nn] == WALL)
      *w2f = 1;
    node->mode[nn] = FLUID;
  }

  return;
}


// ------------------------------------------------------------------------------------
// when geometry change a new velocity field is found.
//      constant flow: find a gx that gives same flux as before, stored in sys->max_vel
//
//      constant pressure: find steady-state solution using initial gx
// ------------------------------------------------------------------------------------
void reinit_velocity(struct Node *node, struct System *sys, struct Field *fluid, int print_status)
{
  int nsteps;
  double conv_crit = sys->conv_crit_vel*sys->conv_length;
  double starttime = MPI_Wtime();
  double value;
  double target = *(sys->vel_target); // pointer to sys->mean_vel or sys->max_vel
  int not_on_target = 1;
  int constant_pressure_run = (sys->const_flow == 0);
  //int i=0, maxiter = 5;
  //double conv_red = 0.5;

  sys->n_reinit++;
  sys->velreinit = 0;
  
  init_fluid(node, fluid, sys);
  do {
    nsteps = advance_fluid_to_convergence(conv_crit, node, fluid, sys, starttime, print_status);
    if (constant_pressure_run) {
      //value = get_umax(node, sys, fluid);
      //if (print_status && sys->mpi->my_rank == 0) printf(" ---> umax = %.3e (%.3e)\n", value, target);
      value = get_mean_velocity(node, sys, fluid);
      if (print_status && sys->mpi->my_rank == 0) printf(" ---> u_mean = %.3e (%.3e)\n", value, target);
      break;
    }
    if (sys->vel_target == &(sys->max_vel)) {
      value = get_flux(node, sys, fluid);
    } else {  // AMIN : HACK NB! NB!
      value = 0; //get_mean_velocity(node, sys, fluid);
    }
    not_on_target = check_value_and_adjust_drive(value, target, sys, fluid);
    if (print_status && sys->mpi->my_rank == 0) printf(" ---> flux = %.3e (%.3e)\n", value, target);
  } while (not_on_target);

  sys->vel_reinit_steps += nsteps;
}


//-------------------------------------------------
//
//-------------------------------------------------
void initialize_velocity(struct Node *node, struct System *sys, struct Field *fluid, int print_status)
{
  int nsteps, not_on_target = 1;
  struct Mpi *mpi = sys->mpi;
  double conv_crit = sys->conv_crit_vel*sys->conv_length;
  double starttime = MPI_Wtime();
  double value;
  double target = sys->max_vel; // from inp.dat
  //int const_pressure_run = (sys->const_flow==0);

  if (print_status && mpi->my_rank == 0) {
    printf("*** Initializing fluid to max velocity = %.2e (<%.1e between %d steps)\n", sys->max_vel, conv_crit, sys->conv_length);
  }

  if (sys->opt->perm_interval>0.0)
    init_fluid(node, fluid, sys);
  do {
    nsteps = advance_fluid_to_convergence(conv_crit, node, fluid, sys, starttime, print_status);
    calculate_permeability(sys, node, fluid, NULL);
    write_permeability_to_file(sys->t, sys, fluid, NULL);
    value = get_umax(node, sys, fluid);
    not_on_target = check_value_and_adjust_drive(value, target, sys, fluid);
    if (print_status && mpi->my_rank == 0) printf("\n");
  } while (not_on_target);

  fluid->flux = get_flux(node, sys, fluid);
  if (sys->const_flow) {
    // max_value holds the target value for re-init velocity runs
    sys->max_vel = fluid->flux;
  }

  //write_fluid_distribution_and_vector_to_file(sys, fluid);
  //write_to_file(sys->output->fluid, sys);
}


//-----------------------------------------------
//  Sum flags of all processes and broadcast.
//  Flag is 0 if all flags are zero, otherwise
//  number of processes with flag = 1
//-----------------------------------------------
void broadcast_sum_of_int_to_all_procs(int *flag, struct Mpi *mpi)
{
  int sum_flags;
  MPI_Reduce(flag, &sum_flags, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  if (mpi->my_rank == 0)
    *flag = sum_flags;
  MPI_Bcast(flag, 1, MPI_INT, 0, MPI_COMM_WORLD);
}


/* ---------------------------------------------------------------------------------------------- */
void check_convergence(int *not_converged, struct System *sys, struct Field *species, double conv_crit, double *diff)
/* ---------------------------------------------------------------------------------------------- */
/* check if species has converged      */
/*                                     */
{
  //std::cout << "Node convergence" << std::endl;
  real *rho, *rho_old;
  struct Step *step = sys->step;
  struct Node *node = sys->node;

  *not_converged = 0;
  sys->rho_not_converged = 0;

  // calculate rho and check for convergence 
  for (int cc = 0; cc < species->n_c; ++cc) { /* loop species */
#ifdef NEWCONC
    rho = &species->psi[step->rho_phase*cc];
#else
    rho = &species->rho[step->rho_phase*cc];
#endif
    rho_old = &species->rho_old[step->rho_phase*cc];

    for (int cn = 0; cn < sys->max_n_tot; ++cn) { /* loop nodes */
      if (node->is_ghost_or_periodic[cn] || node->mode[cn] != FLUID)
        continue;

      if ((rho[cn] > 0.0) && (fabs(rho[cn] / rho_old[cn] - 1.0) > conv_crit)) {
        if (rho[cn] > 0.0) {
          *diff = fabs(rho[cn] / rho_old[cn] - 1.0);
          //*diff = fabs(rho[cn]-rho_old[cn]);
        }
        else {
          *diff = rho[cn];
        }
        *not_converged = 1;
        sys->rho_not_converged = rho[cn];
        sys->node_not_converged = cn;
        break; // breaks node loop
      }
    }
    if (*not_converged)
      break; // breaks species loop
  } /* end species-loop */

#ifdef _MPI_ 
  // broadcast convergence flag to all processes
  //broadcast_flag(sys->mpi, not_converged);
  broadcast_sum_of_int_to_all_procs(not_converged, sys->mpi);
#endif

}

/* ---------------------------------------------------------------------------------------------- */
void check_mean_convergence(int *not_converged, struct System *sys, struct Field *species, double conv_crit, double *diff)
/* ---------------------------------------------------------------------------------------------- */
/* check if mean species concentration has converged      */
/*                                     */
{
  //std::cout << "Mean convergence" << std::endl;
  *not_converged = 0;
  sys->rho_not_converged = 0;
  sys->node_not_converged = 0;

  get_fluid_conc(species, sys);

  // calculate mean concentration
  if (sys->mpi->my_rank == 0) {
    for (int cc = 0; cc < species->n_c; ++cc) {
      //species->mean_conc[cc] = species->total_mole[cc]/species->total_fluid_nodes;
      *diff = fabs(species->mean_conc[cc] / species->mean_conc_old[cc] - 1.0);
      if ((*not_converged == 0) && (*diff > conv_crit)) {
        *not_converged = 1;
        sys->rho_not_converged = species->mean_conc[cc];
        sys->node_not_converged = cc;
      }
      // store current value
      species->mean_conc_old[cc] = species->mean_conc[cc];
    }
  }

#ifdef _MPI_ 
  MPI_Bcast(not_converged, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif

}


/* ---------------------------------------------------------------------------------------------- */
void update_rho_old(struct System *sys, struct Field *species)
/* ---------------------------------------------------------------------------------------------- */
/*  update rho_old with current rho-values to  *
 *  be used in the next convergence check      *
 *  Need to recompute rho since the            *
 *  check_convergence loop only calculates rho *
 *  up to the first non-convergent node        */
{
  int cc, g_ind, cn;
  real *rho, *rho_old;
  struct Step *step = sys->step;
  struct Node *node = sys->node;

  // save current species values for next convergence check
  for (cc = 0; cc < species->n_c; ++cc) { /* loop species */
    //g       = &species->f[step->f_phase*cc];
#ifdef NEWCONC
    rho = &species->psi[step->rho_phase*cc];
#else
    rho = &species->rho[step->rho_phase*cc];
#endif
    rho_old = &species->rho_old[step->rho_phase*cc];
    g_ind = 0;
    for (cn = 0; cn < sys->max_n_tot; ++cn) { /* loop nodes */
      if (node->mode[cn] == FLUID) {
        //calc_rho(&g[g_ind], &rho[cn], sys);
        rho_old[cn] = rho[cn];
      }
      g_ind += step->f_node;
    }
  }

}


// ---------------------------------------------------------
// ---------------------------------------------------------
void reset_rho_vel_in_new_fluid_nodes(struct Field *fluid, struct System *sys)
{
  int nc, nn, ck;
  struct Step *step = sys->step;

  // reset rho and vel ONLY in new fluid nodes
  for (nc = 0; nc < sys->num_new_nodes; ++nc) {
    nn = sys->new_nodes[nc]; // new fluid node
    fluid->rho[nn] = 1.0;
    fluid->u[step->u_node*nn] = fluid->u[step->u_node*nn + 1] = fluid->u[step->u_node*nn + 2] = 0.0;
    for (ck = 0; ck < sys->n_Q; ++ck) {  // direction-loop
      fluid->f[step->f_node*nn + ck] = calc_feq(ck, fluid->rho[nn], &fluid->u[step->u_node*nn], sys);
    }
  }
}


// ---------------------------------------------------------
// ---------------------------------------------------------
int check_velocity_convergence(struct Node *node, struct Field *fluid, struct System *sys, double conv_crit)
{
  int nn;
  //struct Node *node = sys->node;
  struct Step *step = sys->step;
  int not_converged = 1;
  double conv_crit2 = conv_crit*conv_crit;

  // compare ux with stored value; assume convergence until
  // first occurence of too large deviation
  not_converged = 0;
  for (nn = 0; nn < sys->max_n_tot; ++nn) {
    if (node->is_ghost_or_periodic[nn] || node->mode[nn] != FLUID)
      continue;
    //if (node->mode[nn] == FLUID) {
    //ux = fluid->u[step->u_node*nn];
    double *u = &fluid->u[step->u_node*nn];
    double u2 = u[0]*u[0] + u[1]*u[1] + u[2]*u[2];
    if (not_converged == 0) {
      if (fabs(u2 - fluid->u_old[nn]) > conv_crit2) {
        not_converged = 1;
      }
    }
    // store ux for next check
    fluid->u_old[nn] = u2;
    //}
  } // end node-loop

  // root-process broadcasts not_converged-flag to all processes
  broadcast_sum_of_int_to_all_procs(&not_converged, sys->mpi);

  return(not_converged);
}

// ---------------------------------------------------------
//
// ---------------------------------------------------------
int check_velocity_convergence_nphase(struct Node *node, struct Field *fluid, struct System *sys)
// compare ux with stored value; assume convergence until
// first occurrence of too large deviation
{
  static int last_convergence_step = 0;
  static std::vector<double> u_mean_old(sys->n_D, 0.0);
  
  if (sys->t%sys->conv_length == 0) {
    fluid->not_converged = 0;
    if (sys->t - last_convergence_step < 12*sys->conv_length) {
      // wait before making next check
      fluid->not_converged = 1; // skip check, just copy values
    }
    double conv_crit = sys->conv_crit_vel*sys->conv_length;
    double conv_crit2 = conv_crit*conv_crit;
    get_mean_velocity(node, sys, fluid);
    if (fluid->u_mean[0]<0.0)
      return(-1);
    if (fluid->not_converged == 0) {
      for (int d=0; d<sys->n_D; ++d) {
        if (fabs(fluid->u_mean[d]*fluid->u_mean[d] - u_mean_old[d]*u_mean_old[d]) > conv_crit2) {
          fluid->not_converged = 1;
        }
      }
    }
    u_mean_old = fluid->u_mean;

    //    for (int nn = 0; nn < sys->max_n_tot; ++nn) {
    //      if (node->mode[nn] == FLUID && !node->is_ghost_or_periodic[nn]) {
    //        double *u = &fluid->u_tot[sys->step->u_node*nn];
    //        double u2 = u[0]*u[0] + u[1]*u[1] + u[2]*u[2];
    //        if (fluid->not_converged == 0) {
    //          if (fabs(u2 - fluid->u_old[nn]) > conv_crit2) {
    //            fluid->not_converged = 1;
    //          }
    //        }
    //        // store ux for next check
    //        fluid->u_old[nn] = u2;
    //      }
    //    } // end node-loop

    // sum and broadcasts not_converged-flag to all processes
    MPI_Allreduce(MPI_IN_PLACE, &fluid->not_converged, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    if (fluid->not_converged==0) {
      // Convergence!
      last_convergence_step = sys->t;
      double value = fluid->u_mean[0];//get_mean_velocity(node, sys, fluid);
      double target = sys->mean_vel;
      check_value_and_adjust_drive(value, target, sys, fluid);
      if (sys->mpi->my_rank==0)
        printf("\nCONVERGENCE at %d, value = %13.6e, target = %13.6e, gx = %13.6e\n", sys->t, value, target, *(sys->drive));
    }
    return(fluid->not_converged);
  }
  return(-1);
}

// ---------------------------------------------------------
//
// ---------------------------------------------------------
void write_fluid_distribution_and_vector_to_file(struct System *sys, struct Field *fluid)
{
  FILE *fp;
  char fn[50];

  if (sys->mpi->my_rank == 0) {
    sprintf(fn, "%s/bfd/gx.dat", sys->files["out"].c_str());
    fp = my_fopen(fn, "w");
    fprintf(fp, "%.8e\n", fluid->gravity[0]);
    fclose(fp);
  }

  sprintf(fn, "%04d_init_f%s", sys->mpi->my_rank, sys->opt->flip_x ? "__flip_x" : "");
  write_dist_binary(fn, fluid->f, fluid->n_c*sys->max_n_tot*sys->n_Q, sys);
  sys->output->fluid.write(sys);
}

// ---------------------------------------------------------
//
// ---------------------------------------------------------
void init_fluid(struct Node *node, struct Field *fluid, struct System *sys)
{
  reset_rho_vel_in_new_fluid_nodes(fluid, sys); // should be zero, but just in case
  communicate(fluid->f, fluid->n_c, sys);
  communicate_node(node, sys);
  bc_init_run_f(node, fluid, sys);
}

// ---------------------------------------------------------
//
// ---------------------------------------------------------
void advance_fluid(struct Node *node, struct Field *fluid, struct System *sys)
{
  propagate(node, fluid, sys);
#ifdef NEWCONC
  collision_f_newconc(node, fluid, sys);
#else
  collision_f(node, fluid, sys);
#endif
}

// ---------------------------------------------------------
//
// ---------------------------------------------------------
void apply_fluid_bc(struct Node *node, struct Field *fluid, struct System *sys)
{
  //bc_run_f(node, fluid, sys);
#ifdef _FLUID_BDRY_PRESS_
  bc_3D_copy_gradient(node->nr_in , node->list_in , fluid, sys, sys->step->neg_x_dir, sys->step->pos_x_dir);
  bc_3D_copy_gradient(node->nr_out, node->list_out, fluid, sys, sys->step->pos_x_dir, sys->step->neg_x_dir);
#endif
  //periodic_bc(node, fluid, sys);
  bc_run_f(node, fluid, sys);
}

// ---------------------------------------------------------
//
// ---------------------------------------------------------
int advance_fluid_to_convergence(double conv_crit, struct Node *node, struct Field *fluid, struct System *sys, double starttime, int print_status)
{
  int i = 0;
  int not_converged = 1;
  double now, laptime = MPI_Wtime();
  int maxi = 20000 * sys->conv_length;
  int perm_int = 1000;
  double umax = 0.0;
  
  while (not_converged) {

    if (++i >= maxi) {
      if (sys->mpi->my_rank == 0) {
        printf("\n*** WARNING velocity calculation did not converge after %d steps, aborting...\n", maxi);
      }
      MPI_Finalize();
      exit(1);
    }

    advance_fluid(node, fluid, sys);

    if((sys->opt->perm_interval==0.0) && (i%perm_int==0)) {
      calculate_permeability(sys, node, fluid, NULL);
      write_permeability_to_file(i, sys, fluid, NULL);
    }

    if (i%sys->conv_length == 0) {
      // NB
      //sys->output->fluid.write(sys);
      //sys->nwrite++;
      // NB
      if (i == sys->conv_length) {
        update_u_old(node, fluid, sys);
      } else {
        not_converged = check_velocity_convergence(node, fluid, sys, conv_crit);
        umax = get_umax(node, sys, fluid);
        if (print_status && sys->mpi->my_rank == 0) {
          now = MPI_Wtime();
          printf("\r*** time: %s, %4d steps, dt (%d steps) = %4.1f s, not_converged = %3d, %s = %.3e, umax = %.6e",
              get_time_string(now - starttime), i, sys->conv_length, now - laptime, not_converged,
              sys->drive_type, *(sys->drive), umax);
          fflush(stdout);
        }
      }
      laptime = MPI_Wtime();
    }

    apply_fluid_bc(node, fluid, sys);
  }
  sys->output->fluid.write(sys);
  return(i);
}


// ---------------------------------------------------------
//
// ---------------------------------------------------------
int advance_fluid_and_species_to_convergence(double conv_crit, struct Node *node, struct Field *fluid,
  struct Field *species, struct Field *dif, struct InitChem *ICS, struct Minerals *minerals,
  struct Boundary *bndry, struct splayTree **st_lst, struct System *sys,
  double starttime, int print_status)
{
  int i = 0;
  int not_converged = 1;
  double now, laptime = MPI_Wtime();
  int maxi = 2000 * sys->conv_length;

  while (not_converged) {

    if (++i >= maxi) {
      if (sys->mpi->my_rank == 0) {
        printf("\n*** WARNING velocity calculation did not converge after %d steps, aborting...\n", maxi);
      }
      MPI_Finalize();
      exit(1);
    }

    //advance_fluid(node, fluid, sys);
    advance_one_step(ICS, species, fluid, dif, sys, minerals, bndry, st_lst);

    if (i%sys->conv_length == 0) {
      // NB NB
      //if ((sys->opt->perm_interval < 1.1*sys->dt) && (i % 5000 == 0)) {
      //  calc_effluent_and_perm(sys, node, NULL, fluid, 1);
      //  write_permeability_to_file(i, sys, fluid, NULL);
      //}
      // NB NB
      if (i == sys->conv_length) {
        update_u_old(node, fluid, sys);
      }
      else {
        not_converged = check_velocity_convergence(node, fluid, sys, conv_crit);
        if (print_status && sys->mpi->my_rank == 0) {
          now = MPI_Wtime();
          printf("\r*** time: %s, %4d steps, dt (%d steps) = %4.1f s, not_converged = %3d, %s = %.3e",
            get_time_string(now - starttime), i, sys->conv_length, now - laptime, not_converged,
            sys->drive_type, *(sys->drive));
          fflush(stdout);
        }
      }
      laptime = MPI_Wtime();
    }

  }

  return(i);
}



// ---------------------------------------------------------
// Update u_old with x-component of velocity vector
// ---------------------------------------------------------
void update_u_old(struct Node *node, struct Field *fluid, struct System *sys)
{
  int nn;
  for (nn = 0; nn < sys->max_n_tot; ++nn) {
    if (node->mode[nn] == FLUID && !node->is_ghost_or_periodic[nn]) {
      fluid->u_old[nn] = fluid->u[sys->step->u_node*nn];
    }
  }
}


//---------------------------------
//
//---------------------------------
int check_value_and_adjust_drive(double value, double target, struct System *sys, struct Field *fluid)
{
  //double target_conv_limit = 1e-3;
  int value_not_on_target;

  //  if ((value < *(sys->drive)*0.001) || (value < _EPS_)) {
  //  // value is too low
  //  if (sys->mpi->my_rank == 0) {
  //    printf("\n\n            !!!!! ERROR !!!!!\n\n");
  //    printf(" Velocity is too low (umax = %.3e, %s = %.3e)!\n", value, sys->drive_type, *(sys->drive));
  //    printf(" Possible reasons:\n");
  //    printf(" 1) The cell is closed: set closed_cell = 0 in options.dat\n");
  //    printf(" 2) No confining walls parallel to flow: set add_in_out_wall = 1 in options.dat)\n");
  //    printf(" 3) Wrong inlet/outlet conditions in geo.dat (e.g. non-periodic inlet/outlet)\n");
  //  }
  // MPI_Finalize();
  //  exit(1);
  //}

  // assume value on target and check if assumption holds (relative convergence criteria)
  value_not_on_target = 0;
  if (fabs((target - value) / target) > sys->conv_crit_vel_target) {
    // value too far from target
    value_not_on_target = 1;
    // adjust drive
    //fluid->gx_change = *(sys->drive)*((target/value) - 1.0);
    *(sys->drive) *= (target / value);
#ifdef _FLUID_BDRY_PRESS_
    init_run_f(sys->node, fluid, sys);
    //reset_rho_vel_in_new_fluid_nodes(fluid, sys); // should be zero, but just in case
    //communicate(fluid->f, fluid->n_c, sys);
    //communicate_node(node, sys);
    bc_init_run_f(sys->node, fluid, sys);
    //init_fluid(sys->node, fluid, sys);
#endif
  }

  return(value_not_on_target);
}

//--------------------------------------------
//
//--------------------------------------------
void init_permeability_measurement(struct System *sys, struct Field *fluid, struct Field *species, struct Minerals *minerals)
{
  struct Node *perm_node, *node, *pn;
  int n;
  int flow_cpy = sys->const_flow;
  int print_status = 1; //_DEBUG_PERM_;

  if (sys->mpi->my_rank == 0)
    printf("\n*** Initializing permeability measurement\n");

  if (sys->opt->closed_cell && fluid->gravity[0]>0.0) {
    // allocate
    //sys->perm_node = (struct Node *) malloc(1 * sizeof(struct Node));
    sys->perm_node = new Node();
    pn = sys->perm_node;
    //pn->links = (struct Links *) malloc(1 * sizeof(struct Links));
    pn->links = new Links();
    pn->mode = (char *)calloc(sys->max_n_tot, sizeof(char));
    pn->is_ghost = (char *)calloc(sys->max_n_tot, sizeof(char));
    pn->is_ghost_or_periodic = (char *)calloc(sys->max_n_tot, sizeof(char));
    pn->nr_per = pn->nr_in = pn->nr_out = pn->nr_ghost = pn->nr_wall = 0;
    pn->is_fluid_in_out[FLUID] = pn->is_fluid_in_out[INLET] = pn->is_fluid_in_out[OUTLET] = 1;
    pn->is_fluid[FLUID] = 1;
    init_node_lists(pn, sys);

    node = sys->node;
    perm_node = sys->perm_node;

    setup_periodic_bc(perm_node, sys);

    setup_ghost_nodes(perm_node, sys);

    for (n = 0; n < node->nr_ghost; ++n) {
      if (node->list_ghost[n] != perm_node->list_ghost[n]) {
        printf("         !!! ERROR !!!\n ");
        printf("Permeability nodes have different ghost-nodes!\n");
        MPI_Finalize();
        exit(1);
      }
    }

    add_outer_walls(perm_node, sys);

    copy_node_mode(node, perm_node, sys);

    set_wall_nodes(perm_node, sys);

    sys->output->write_geo_file("perm_geo", perm_node, sys);

    init_linklist(perm_node, sys);

    allocate_init_f(fluid, sys);

    sys->output->init_fluid(fluid, perm_node, sys);

    init_run_f(perm_node, fluid, sys);
  } else {
    perm_node = sys->node;
  }

  sys->const_flow = 1;
  initialize_velocity(perm_node, sys, fluid, print_status);
  sys->const_flow = flow_cpy;

  calculate_permeability(sys, perm_node, fluid, species);
  write_permeability_to_file(sys->t, sys, fluid, minerals);

  //if ((sys->opt->perm_interval < 1.1*sys->dt) && (sys->mpi->my_rank == 0)) {
  if (sys->mpi->my_rank == 0) {
    printf("*** fluid_flux = %.3e, mean_velx = %.3e, perm = %.3e D\n", fluid->flux, fluid->mean_velx, fluid->perm);
    fflush(stdout);
  }
  // set flux to use as source-term in collision
  //sys->source_flux = fluid->flux*sys->opt->perm_vel/fluid->mean_velx;

  //if (print_status && sys->mpi->my_rank==0) {
  //printf("*** source_flux = %.2e, fluid_flux = %.2e, mean_velx = %g\n", sys->source_flux, fluid->flux, fluid->mean_velx);
  //}
}


//--------------------------------------------
//
//--------------------------------------------
void measure_permeability(struct System *sys, struct Field *fluid, struct Field *species, struct Minerals *minerals)
{
  struct Node *perm_node = sys->perm_node;
  struct Node *node = sys->node;
  int flow_cpy = sys->const_flow;
  int print_status = _DEBUG_PERM_;

  if (print_status && sys->mpi->my_rank == 0)
    printf("\n*** Starting permeability measurement at step %d\n", sys->t + 1);

  copy_node_mode(node, perm_node, sys);

  set_wall_nodes(perm_node, sys);

  update_linklist(perm_node, sys);

  init_run_f(perm_node, fluid, sys);  // zero out

  sys->const_flow = 1; // if 0 the routine stops after first convergence (0 means a constant pressure run)
  reinit_velocity(perm_node, sys, fluid, print_status);
  sys->const_flow = flow_cpy;

  calculate_permeability(sys, perm_node, fluid, species);
  write_permeability_to_file(sys->t, sys, fluid, minerals);

  //write_to_file(sys->output->fluid, sys);

}

//--------------------------------------------
// calculate_permeability
// turn chem off (exclude effluent) and fluid on
//--------------------------------------------
void calculate_permeability(struct System *sys, struct Node *node, struct Field *fluid, struct Field *species)
{
  // copy values
  int chem_cpy = sys->chem;
  int fluid_cpy = sys->fluid;
  //double ratio_cpy = sys->dt_ratio;
  int broadcast = 1;

  //sys->dt_ratio = 1.0;  // perm measurement performed at tau=1
                        // and species->tau can be > 1 => dt_ratio != 1
  sys->chem = 0;        // chem off to avoid effluent values
  sys->fluid = 1;       // fluid on to get perm values

  calc_effluent_and_perm(sys, node, species, fluid, broadcast);

  // restore values
  //sys->dt_ratio = ratio_cpy;
  sys->chem = chem_cpy;
  sys->fluid = fluid_cpy;
}

//--------------------------------------------
//
//--------------------------------------------
void write_permeability_to_file(int timestep, struct System *sys, struct Field *fluid, struct Minerals *minerals)
{
  char fn[50];
  FILE *fp = NULL;
  int min;
  struct Mpi *mpi = sys->mpi;
  static int write_header = 1;
  double phi = 0.0, SA = 0.0;
  int nfluid_tot, nfluid_reac;
  //double gradP_mu;
  // filename
  //double *solid_mass_diff, *fluid_mole_diff;
  sprintf(fn, "%s/permeability.dat", sys->files["out"].c_str());

  // create file and write header if first call
  if (write_header && (mpi->my_rank == 0)) {
    write_header = 0;
    fp = fopen(fn, "w");
    if (fp == NULL) {
      printf("ERROR in write_permeability_to_file: fopen returned %d\n", errno);
      MPI_Finalize();
      exit(0);
    }
    fprintf(fp, "# %-11s %-11s %-11s %-11s %-11s %-11s %-11s %-11s %-11s",
      "timestep", "realtime", "gx", "area", "vel[m/s]", "gradP_mu", "perm[D]", "porosity", "SA[m2]");
    // mineral SA
    if (minerals != NULL) {
      for (min = 0; min < minerals->n_tot; ++min) {
        fprintf(fp, " SA_%-12s", minerals->list[min].name);
      }
    }
    fprintf(fp, "\n");
    fclose(fp);
  }

  // pore volume
  get_pore_volume(sys->node, &nfluid_tot, &nfluid_reac, sys);
  phi = (double)nfluid_reac / (double)(sys->V_lb);

  // surface area
  if (minerals != NULL) {
    get_reactive_surface_area(sys->node, &SA, sys, minerals);
    SA *= sys->dx*sys->dx;
  }

  // write flow and perm to file
  if (mpi->my_rank == 0) {
    fp = fopen(fn, "a");
    // flow data
    fprintf(fp, "%11d %11.5e %11.5e %11.5e %11.5e %11.5e %11.5e %11.5e %11.5e",
      timestep + 1, sys->time + sys->dt, fluid->gravity[0], fluid->area, fluid->mean_velx,
      fluid->gravity[0] * sys->dx / (sys->nu*fluid->dt*fluid->dt),
      fluid->perm, phi, SA);
    // mineral SA
    if (minerals != NULL) {
      for (min = 0; min < minerals->n_tot; ++min) {
        fprintf(fp, " %15.9e", minerals->list[min].SA / minerals->SA_tot);
      }
    }
    fprintf(fp, "\n");
    fclose(fp);
  }

}

//--------------------------------------------
// copy from from_node to to_node
//  - exclude nodes outside Xlim
//  - set WALL nodes to SOLID nodes
//--------------------------------------------
void copy_node_mode(struct Node *from_node, struct Node *to_node, struct System *sys)
{
  int cn, X, m;
  int Xlim[2] = {sys->opt->in_nodes+1, sys->MAX_N[0]-2-sys->opt->out_nodes};
  //int Xlim[2] = { sys->opt->in_nodes, sys->MAX_N[0] - 1 - sys->opt->out_nodes };
  for (cn = 0; cn < sys->max_n_tot; ++cn) {
    if (from_node->is_ghost_or_periodic[cn])
      continue;
    X = get_global_xcoord(cn, sys);
    if (X > Xlim[0] && X < Xlim[1]) {
      m = from_node->mode[cn];
      if (m == FLUID || m == INLET || m == OUTLET) {
        to_node->mode[cn] = FLUID;
      } else {
        to_node->mode[cn] = SOLID;
      }
    }
  }
  communicate_node(to_node, sys);
}
