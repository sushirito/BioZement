#include "global.h"
#include "chem_global.h"
#include "output.h"
#include "biozementfunctions.h"

#ifdef _FLUID_BDRY_PRESS_
//------------------------------------------------------
//
//------------------------------------------------------
void set_rho_inlet_outlet(Node *node, Field *fluid)
{
    for (int i=0; i<fluid->n_c; ++i) {
        double delta_P = fluid->grad_P*(node->out_x - node->in_x);
        fluid->rho_outlet[i] = fluid->rho_init[i] - 0.5*delta_P;
        fluid->rho_inlet[i]  = fluid->rho_init[i] + 0.5*delta_P;
//    std::cout << i << ": rho_inlet, rho_outlet = " << fluid->rho_inlet[i] << ", "<< fluid->rho_outlet[i]
//              << "grad_P: " << fluid->grad_P << "delta_P: "<< delta_P << std::endl;
    }
}
#endif

/* ---------------------------------------------------------------------------------  init_run_f */
void init_run_f(struct Node *node, struct Field *fluid, struct System *sys)
/*
init_run_f:
    initiate a run for fluid

  INPUT : void
  OUTPUT : void
 */
{
    int cn, cc, ck, cd; /* counter: node, phase, direction, dimension */
    int f_ind, u_ind, rho_ind; /* Node index, velocity index, density index */
    real *f; /* Pointer to the first node for a given phase */
    real *u; /* Pointer to the first node for the advective velocity for a phase */
    //  real *u_init; /* Pointer to the initial velocity for a phase */
    real *rho; /* Pointer to first node for the density for a phase */
    real rho_init; /* Initial value for the density of a phase */
#ifdef _USE_COLOR_GRAD_
    int nxyz[3];
    double rand_max = (double) RAND_MAX;
    double sum_rand;
    double pertubation_factor = 0.1;
    double rho_init_min;
    double *rand_list;

    rand_list = (double*) malloc(sizeof(double)*fluid->n_c);

    rho_init_min = rand_max;
    for (cc = 0; cc < fluid->n_c; cc++) {
        if (rho_init_min > fluid->rho_init[cc])
            rho_init_min = fluid->rho_init[cc];
    }

    pertubation_factor *= 2*rho_init_min;

    srand(time(NULL) + sys->mpi->my_rank);

#endif

    struct Step *step = sys->step;

#ifdef _FLUID_BDRY_PRESS_
    printf("USING _FLUID_BDRY_PRESS_ \n");
    set_rho_inlet_outlet(sys->node, fluid);
#endif

    for (cc = 0; cc < fluid->n_c; ++cc) { /* phase */
        //printf("cc = %d\n", cc);
        f = &fluid->f[step->f_phase*cc];       /* distribution in the current phase */
        u = &fluid->u[step->u_phase*cc];       /* velocities in the current phase */
        //u_init = &fluid->u_init[sys->n_D*cc];  /* initial velocity for the current phase */
        rho = &fluid->rho[step->rho_phase*cc]; /* density in the current phase */
        rho_init = fluid->rho_init[cc];        /* initial density for the current phase */
//#ifdef _FLUID_BDRY_PRESS_
//    double *f_tmp = &fluid->f_tmp[step->f_phase*cc];
//#endif
        f_ind = 0;                    /* index for nodes distribution */
        u_ind = 0;                    /* index for velocities */
        rho_ind = 0;                  /* index for densities */

        for (cn = 0; cn < sys->max_n_tot; ++cn) {
            if (node->mode[cn] != WALL && node->mode[cn] != SOLID) {
#ifdef _USE_COLOR_GRAD_
                get_global_coord(cn, nxyz, sys);
                if (cc == 0) {
                    sum_rand = 0;
                    for (int n=0; n < fluid->n_c; n++) {
                        rand_list[n] = rand()/rand_max;
                        sum_rand += rand_list[n];
                    }
                    sum_rand = sum_rand/fluid->n_c;

                    for (int n=0; n < fluid->n_c; n++) {
                        fluid->rho[step->rho_phase * n + cn] = fluid->rho_init[n];
                        fluid->rho[step->rho_phase * n + cn] += pertubation_factor*(rand_list[n] - sum_rand);
                    }

                    /*          if ( (nxyz[0] > 10) && (nxyz[0] < 30) &&
                     (nxyz[1] > 10) && (nxyz[1] < 30) &&
                     (nxyz[2] > 10) && (nxyz[2] < 30)) {
                      fluid->rho[step->rho_phase * 0 + cn] = 1;
                      fluid->rho[step->rho_phase * 1 + cn] = 0;
                      fluid->rho[step->rho_phase*2 + cn ] = 0;
                    } else {
                      fluid->rho[step->rho_phase * 0 + cn] = 0;
                      fluid->rho[step->rho_phase * 1 + cn] = 1;
                      fluid->rho[step->rho_phase*2 + cn ] = 0;
                    } */
                }
#else
                rho[rho_ind] = rho_init;
#ifdef _FLUID_BDRY_PRESS_
                int x = get_global_xcoord(cn, sys);
                if ((x>=node->in_x) && (x<=node->out_x)) {
                    rho[rho_ind] = fluid->rho_inlet[cc] - fluid->grad_P*(x-node->in_x);
                }
                //        if (x==node->in_x)
                //          rho[rho_ind] = fluid->rho_inlet[cc];
                //        if (x==node->out_x)
                //          rho[rho_ind] = fluid->rho_outlet[cc];
#endif
#endif
            } else {
                rho[rho_ind] = 0.0;
                rho[rho_ind] = rho_init;  // AMIN

            }

            /* init velocity */
            for( cd = 0; cd < sys->n_D; ++cd )
                u[u_ind + cd] = 0.0; //u_init[cd];

            /* init velocity distribution */
            if (node->mode[cn] != WALL && node->mode[cn] != SOLID) {
                for (ck = 0; ck < sys->n_Q; ++ck) {
                    f[f_ind + ck] = calc_feq(ck, rho[rho_ind], &u[u_ind], sys);
//#ifdef _FLUID_BDRY_PRESS_
//          f_tmp[f_ind + ck] = f[f_ind + ck];
//#endif
                }
            }

            f_ind += step->f_node;
            u_ind += step->u_node;
            rho_ind += step->rho_node;
        }
    }

#ifdef _USE_COLOR_GRAD_
    free(rand_list);
#endif
}


/* ---------------------------------------------------------------------------------  init_run_g */
void init_run_g(struct Field *species, struct Field *fluid, struct System *sys)
/*
init_run_g:
    initiate a run for diffusion

  INPUT : void
  OUTPUT : void
 */
{
    int cn, cc, ck; /* counter: node, phase, direction, dimension */
    int g_ind, u_ind, rho_ind; /* Node index, velocity index, density index */
    real *g; /* Pointer to the first node for a given phase */
    real *u; /* Pointer to the first node for the advective velocity for a phase*/
    real *rho; /* Pointer to first node for the density for a phase */
    real rho_init; /* Inital value for the density of a phase */

    struct Step *step = sys->step;
    struct Node *node = sys->node;

#ifdef FLUIDCHEM
    u = fluid->u;
#endif

    for (cc = 0; cc < species->n_c; ++cc) {  /* phase */
        g = &species->f[step->f_phase*cc];  /* distribution in the current phase */
#ifndef FLUIDCHEM
        u = &species->u[step->u_phase*cc]; /* velocities in the current phase */
#endif
        rho = &species->rho[step->rho_phase*cc]; /* Density for the current phase */
        rho_init = species->rho_init[cc];   /* initial density for the current phase */

        g_ind = 0;              /* index for nodes distribtuion */
        u_ind = 0;              /* index for velocities */
        rho_ind = 0;            /* index for densities */

        for (cn = 0; cn < sys->max_n_tot; ++cn) { /* node */
            /* init density */
            if (node->mode[cn] != WALL && node->mode[cn] != SOLID)
                rho[rho_ind] = rho_init;
            else
                rho[rho_ind] = 0.0;

            /* init velocity
            for( cd = 0; cd < sys->n_D; ++cd )
            u[u_ind + cd] = 0.0; */

            /* init velocity distribution */
            if (node->mode[cn] != WALL && node->mode[cn] != SOLID)
                for (ck = 0; ck < sys->n_Q; ++ck)
                    g[g_ind + ck] = calc_geq(ck, rho[rho_ind], &u[u_ind], sys);

            g_ind += step->f_node;
            u_ind += step->u_node;
            rho_ind += step->rho_node;
        }
    }
}


/* ---------------------------------------------------------------------------------  init_run_3D_g */
void init_run_3D_g(struct Field *species, struct Field *fluid, struct System *sys)
/*
init_run_3D_g:
    initiate a run for diffusion in three dimensions

  INPUT  : 1) fn : filename for the 'rho' distribution

  OUTPUT : void
 */
{
    int cn, cc, ck; /* counter: node, phase, direction, dimension */
    int g_ind, u_ind, rho_ind; /* Node index, velocity index, density index */
    real *g; /* Pointer to the first node for a given phase */
    real *u; /* Pointer to the first node for the advective velocity for a phase*/
    real *rho; /* Pointer to first node for the density for a phase */
    real rho_init; /* Inital value for the density of a phase */
    real* R_cc; /* AMIN EJE */
    int cn_neig, is_wall_node;


#ifdef NEWCONC
    real *psi;
#endif

    struct Step *step = sys->step;
    struct Node *node = sys->node;

    u = fluid->u;

    for (cc = 0; cc < species->n_c; ++cc) {  /* phase */
        g = &species->f[step->f_phase*cc];  /* distribution in the current phase */
        //std::cout << "step->f_phase*cc = " << step->f_phase*cc << std::endl; /* AMIN */
        rho = &species->rho[step->rho_phase*cc]; /* Density for the current phase */
#ifdef NEWCONC
        psi = &species->psi[step->rho_phase*cc]; /* Density for the current phase */
#endif
        rho_init = species->rho_init[cc];   /* initial density for the current phase */
        R_cc = &species->R_local[step->rho_phase*cc];

        g_ind = 0;              /* index for nodes distribtuion */
        u_ind = 0;              /* index for velocities */
        rho_ind = 0;            /* index for densities */

        for (cn = 0; cn < sys->max_n_tot; ++cn) { /* node */
            /* init density */
            //if( node->mode[cn] != WALL && node->mode[cn] != SOLID)
#ifdef NEWCONC
            //if( node->mode[cn] == FLUID) {
            // need to include INLET and OUTLET to avoid zero to be propagated in first step
            if (node->mode[cn] == FLUID || node->mode[cn] == INLET || node->mode[cn] == OUTLET) { /* Inside the fluid */
                psi[rho_ind] = rho_init;
                /* if (cc == 2) { AMIN
                    is_wall_node = 0;
                	for (ck = 0; ck < sys->n_Q; ++ck)
                	{
                		cn_neig = cn - step->n_rel_neig[ck];
                		if (node->mode[cn_neig] == WALL)
                			is_wall_node = 1;
                	}
                     if (!is_wall_node) { AMIN
                        R_cc[rho_ind] = 0.0;
                        std::cout << rho_ind << std::endl;
                    }
                } */
#ifdef FLUID_OFF
                rho[rho_ind] = psi[rho_ind];
#else
                rho[rho_ind] = psi[rho_ind] * fluid->rho[rho_ind];
#endif
            } else {
                psi[rho_ind] = rho[rho_ind] = 0.0;
            }
#else
            //if( node->mode[cn] == FLUID)
            // need to include INLET and OUTLET to avoid zero to be propagated in first step
            if (node->mode[cn] == FLUID || node->mode[cn] == INLET || node->mode[cn] == OUTLET)
                rho[rho_ind] = rho_init;
            else
                rho[rho_ind] = 0.0;
#endif

            /* init velocity distribtuion */
            if (node->mode[cn] != WALL && node->mode[cn] != SOLID) {
                for (ck = 0; ck < sys->n_Q; ++ck) {
#ifdef NEWCONC // true for NEWCONC on and FLUID_OFF off 
#ifdef FLUID_OFF
                    g[g_ind + ck] = calc_geq(ck, rho[rho_ind], &u[u_ind], sys);
#else // fluid is on
                    double ev_F = 0;
                    for (int cd = 0; cd < sys->n_D; ++cd) {
                        ev_F += sys->ev[sys->ev_ind[ck] + cd] * fluid->gravity[cd];
                    }
                    g[g_ind + ck] = calc_feq(ck, rho[rho_ind], &u[u_ind], sys) + 0.5*ev_F / sys->cs_2;
#endif
#else
                    g[g_ind + ck] = calc_geq(ck, rho[rho_ind], &u[u_ind], sys);
#endif
                }
            }
            g_ind += step->f_node;
            u_ind += step->u_node;
            rho_ind += step->rho_node;
        }
    }
}


/* --------------------------------------------------------------------------------- */
void init_run_g_newconc(struct Field *species, struct Field *fluid, struct System *sys)
/*
init_run_g:
    initiate a run for diffusion

  INPUT : void
  OUTPUT : void
 */
{
    int cn, cc, ck; /* counter: node, phase, direction, dimension */
    int g_ind, u_ind, rho_ind; /* Node index, velocity index, density index */
    //int nx, ny, nz;
    real *g; /* Pointer to the first node for a given phase */
    real *u; /* Pointer to the first node for the advective velocity for a phase*/
    real *rho, *psi; /* Pointer to first node for the density for a phase */
    //real rho_init; /* Inital value for the density of a phase */

    struct Step *step = sys->step;
    struct Node *node = sys->node;

    u = fluid->u;

    for (cc = 0; cc < species->n_c; ++cc) {  /* phase */
        g = &species->f[step->f_phase*cc];  /* distribution in the current phase */
        rho = &species->rho[step->rho_phase*cc]; /* Density for the current phase */
        psi = &species->psi[step->rho_phase*cc]; /* Density for the current phase */
        /*    rho_init = species->rho_init[cc];    initial density for the current phase */

        g_ind = 0;              /* index for nodes distribtuion */
        u_ind = 0;              /* index for velocities */
        rho_ind = 0;            /* index for densities */

        for (cn = 0; cn < sys->max_n_tot; ++cn) { /* node */
            /* init density */
            if (node->mode[cn] != WALL && node->mode[cn] != SOLID)
                psi[rho_ind] = rho[rho_ind] = 1.0;
            else
                psi[rho_ind] = rho[rho_ind] = 0.0;


            /* init velocity distribution */
            if (node->mode[cn] != WALL && node->mode[cn] != SOLID)
                for (ck = 0; ck < sys->n_Q; ++ck)
                    g[g_ind + ck] = calc_feq(ck, rho[rho_ind], &u[u_ind], sys);

            g_ind += step->f_node;
            u_ind += step->u_node;
            rho_ind += step->rho_node;
        }
    }
}

void init_mpi_run(struct System *sys, int *argc, char ***argv, struct Field *species)
/************************************
 ************************************/
{
    int i, c;
    struct Mpi *mpi = sys->mpi;

    MPI_Init(argc, argv);                        // start up _MPI_
    MPI_Comm_size(MPI_COMM_WORLD, &mpi->nr_procs); // number of processes
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi->my_rank);  // process rank
    if (mpi->np[0] * mpi->np[1] * mpi->np[2] != mpi->nr_procs) {
        fprintf(stderr, "ERROR! Wrong number of processes: %d*%d*%d != %d\n",
                mpi->np[0], mpi->np[1], mpi->np[2], mpi->nr_procs);
        exit(1);
    }

    /* DATA PARTITION AND LOAD BALANCING */
    int *np, *irank, rank;
    np = mpi->np;
    irank = mpi->ind_rank;
    rank = mpi->my_rank;

    irank[0] = rank%np[0];
    irank[1] = ((rank - irank[0]) / np[0]) % np[1];
    irank[2] = (int)(rank / (np[1] * np[0]));

    int n_rest, n_tmp;
    int *n = sys->max_n;
    int *lb = mpi->lower_bound;
    int *ub = mpi->upper_bound;
    lb[2] = ub[2] = 0;  // compatibility with 2D runs

    for (i = 0; i < sys->n_D; ++i) {
        sys->MAX_N[i] = 2 + n[i];
        n_tmp = n[i] / np[i];
        n_rest = n[i] % np[i];
        lb[i] = irank[i] * n_tmp;
        lb[i] += (irank[i] <= n_rest) ? irank[i] : n_rest;
        if (irank[i] + 1 <= n_rest)
            ++n_tmp;
        n[i] = n_tmp;
        lb[i] += 1;
        ub[i] = lb[i] + n[i] - 1;
    }

    //
    // set up communication group and new communicator
    // for the effluent calculations

    int nprocs = np[1] * np[2];
    int *eff_ranks = (int *)malloc(nprocs * sizeof(int));

    // only effluent from the outlet-end
    for (c = 0; c < nprocs; c++) {
        eff_ranks[c] = np[0] * (c + 1) - 1;
    }
    mpi->eff_root_rank = eff_ranks[0];

    // create new effluent group
    MPI_Comm_group(MPI_COMM_WORLD, &mpi->world_group);
    MPI_Group_incl(mpi->world_group, nprocs, eff_ranks, &mpi->eff_group);

    // create new effluent communicator
    MPI_Comm_create(MPI_COMM_WORLD, mpi->eff_group, &mpi->EFF_COMM);
    //MPI_Comm_create_group(MPI_COMM_WORLD, mpi->eff_group, 0, &mpi->EFF_COMM);
    MPI_Group_rank(mpi->eff_group, &mpi->eff_rank);

    // START alternative implementation
    //int color = (irank[0]==np[0]-1)? 1 : 0;
    // END alternative implementation

    // allocate send buffer for effluent calculation
    mpi->sbuf_size = 2 + species->n_c;
    if (!(mpi->sbuf_eff = (double *)malloc((mpi->sbuf_size)*nprocs * sizeof(double))))
        printf("ERROR! Unable to allocate memory\n");

    // allocate receive buffer for effluent calculation only at recvroot
    //if (mpi->eff_rank == 0) {
    if (!(mpi->rbuf_eff = (double *)malloc((3 + species->n_c)*nprocs * sizeof(double))))
        printf("ERROR! Unable to allocate memory\n");
    //}

    // cleanup
    free(eff_ranks);
}


void init_chem_splay_mineral(struct System *sys, struct Field *species, struct InitChem *InitChemSolver,
                             struct splayTree ***st_lst, struct BasVec **Vchem_key_tmp, struct Minerals *minerals)
/************************************
 ************************************/
{
    /* Read database etc. BULK and INLET CHEM can be used
         if SURFACE_CHEM = 1 we have to use SURFACE*/

    double c_h_inlet, c_h_bulk;
    int i, no_comb, key_tmp_size;

    init_chemistry(&InitChemSolver[0], sys->files["data"]);
    init_chemistry(&InitChemSolver[1], sys->files["data"]);

    chem_get_pos(&InitChemSolver[0]); /* INLET */
    chem_get_pos(&InitChemSolver[1]); /* BULK */

    db_set_activity(&InitChemSolver[0]);
    db_set_activity(&InitChemSolver[1]);

    for (i = 0; i < InitChemSolver[0].num_phases; ++i) {
        species->name[i] = InitChemSolver[0].SM_basis->row_name[InitChemSolver[0].pos[i]];
    }

    c_h_inlet = chem_set_eq_conc_nlin(species->rho_inlet, &InitChemSolver[0], 1);
    c_h_bulk = chem_set_eq_conc_nlin(species->rho_init, &InitChemSolver[1], 1);

    minerals->no_comb = (int)pow(2, InitChemSolver[1].size_sup_min);
#ifdef _MAGN_ON_MAGN_
    minerals->no_comb++; // solutions for precipitating magnesite nodes kept in separate splaytree
#endif
    no_comb = minerals->no_comb;
    //(*Vchem_key_tmp) = (struct BasVec *) malloc(no_comb * sizeof(struct BasVec));
    *Vchem_key_tmp = new BasVec[no_comb]();
    gen_vchem_vector(no_comb, *Vchem_key_tmp, &InitChemSolver[1]);

    //(*st_lst) = (struct splayTree**) malloc(no_comb * sizeof(struct splayTree *));
    *st_lst = new splayTree*[no_comb]();

    for (i = 0; i < no_comb; ++i) {
        (*Vchem_key_tmp)[i].log_m[(*Vchem_key_tmp)[i].pos_pH] = -7;
        (*Vchem_key_tmp)[i].chbal = sys->opt->charge_balance;
        key_tmp_size = (*Vchem_key_tmp)[i].ICS->size;
        if (!((*Vchem_key_tmp)[i].chbal)) {
            add_H_to_mbal(c_h_bulk, &((*Vchem_key_tmp)[i]));
            species->rho_inlet[key_tmp_size] = c_h_inlet;
            species->rho_init[key_tmp_size] = c_h_bulk;
            key_tmp_size++;
        }
        (*st_lst)[i] = InitializeSplayTree(key_tmp_size, key_tmp_size, 1, sys->splaytree_res, 1, NULL, NULL);
        (*st_lst)[i]->min_key = i;
        (*st_lst)[i]->Vchem_key = &((*Vchem_key_tmp)[i]);
    }

    if (sys->opt->charge_balance) {
        add_H_to_ICS(&InitChemSolver[0], (*Vchem_key_tmp)[0].pos_pH);
        add_H_to_ICS(&InitChemSolver[1], (*Vchem_key_tmp)[0].pos_pH);
    }

    init_solid_fraction(sys, species, minerals, &InitChemSolver[1]);

    get_solid_massfraction(minerals, sys);

    for (i = 0; i < minerals->n_tot; ++i) {
        minerals->sum_sfrac_old[i] = minerals->sum_sfrac_prev[i] = minerals->sum_sfrac[i];
    }

    if (sys->mpi->my_rank == 0) {
        printf("\nCONCENTRATION\n");
        printf("             inlet   (input.dat)         bulk    (input.dat)\n");
        for (i = 0; i < InitChemSolver[0].num_phases; ++i) {
            printf("  %-7s %10.3e (%10.3e) %14.3e (%10.3e)\n", species->name[i],
                   species->rho_inlet[i], InitChemSolver[0].c_vchem[i], species->rho_init[i], InitChemSolver[1].c_vchem[i]);
        }
        printf("\nMASS\n");
        for (i = 0; i < minerals->n_tot; ++i) {
            printf("  %-15s %15.4e [sf] %12.2e [cm3/mol]\n", minerals->list[i].name, minerals->sum_sfrac[i], minerals->list[i].cm3_mol);
        }
        printf("\n");
    }
}



void init_3D_run(struct System *sys, struct Field *species, struct Field *fluid, struct Field *dif, struct InitChem **InitChemSolver,
                 struct splayTree ***st_lst, struct Minerals *minerals, struct Boundary *bndry, struct BasVec **Vchem_key_tmp)
/************************************
 ************************************/
{

    if (sys->opt->find_isolated_nodes) {
        find_isolated_fluid_nodes(species, fluid, sys);
        MPI_Finalize();
        exit(0);
    }

    if (sys->fluid) {
        sys->output->init_fluid(fluid, sys->node, sys);
        init_3D_fluid_run(sys, fluid);
        communicate(fluid->f, fluid->n_c, sys);
        communicate(fluid->f_tmp, fluid->n_c, sys);
        //double q=0;
        //calc_rho_vel_loop(fluid->f, fluid, sys, &q);
        sys->output->fluid.write(sys);
    }

    communicate_node(sys->node, sys);
    if (sys->real_units_in_input) {
        set_parameters(sys, species, dif, fluid, &((*InitChemSolver)[1]));
    } else {
        set_parameters_LB(sys, species, fluid, &((*InitChemSolver)[1]));
    }

    if (sys->chem) {
        real c_init_original[species->n_c];
        //for (int cc = 0; cc < species->n_c; ++cc) // AMIN & EJE
        //    c_init_original[cc] = species->rho_init[cc];

        init_chem_splay_mineral(sys, species, *InitChemSolver, st_lst, Vchem_key_tmp, minerals);
        sys->output->init_chem(species, InitChemSolver, minerals, sys);

        //for (int cc = 0; cc < species->n_c; ++cc) // AMIN & EJE
        //    species->rho_init[cc] = c_init_original[cc];

        init_3D_chemistry_run(sys, species, fluid, *InitChemSolver, *st_lst, minerals, bndry);
        communicate(species->f, species->n_c, sys);
        communicate(species->f_tmp, species->n_c, sys);


        init_run_3D_dif(dif, fluid, sys);
 
        sys->output->init_dif(dif, sys->node, sys);

#ifndef BB_INLET
        periodic_bc(sys->node, dif, sys);
#endif
#ifndef FLUID_OFF
        inlet_bc_g(species, dif, sys, bndry);
#endif
        wall_bc_bounce_back(sys->node, dif, sys);
        communicate(dif->f, dif->n_c, sys);
        communicate(dif->f_tmp, dif->n_c, sys);

    }

    print_run_parameters(sys, species, fluid, &((*InitChemSolver)[1]), minerals);
}

void init_3D_fluid_run(struct System *sys, struct Field *fluid)
/************************************
 ************************************/
{
    char fn[50];
    int i;
    struct Mpi *mpi = sys->mpi;
    FILE *fp;
    //double gx_copy;
    int nsteps = 0;
    struct Options *opt = sys->opt;
    struct Node *node = sys->node;
    std::string out(sys->files["out"]);

    if (mpi->my_rank == 0)
        printf("***\n*** Starting constant %s simulation***\n", (sys->const_flow) ? "FLOW" : "PRESSURE");

    init_run_f(node, fluid, sys);
    communicate_node(node, sys);
    bc_init_run_f(node, fluid, sys);
    communicate(fluid->f, fluid->n_c, sys);

    //fluid->saturated_nodes = (double *)malloc(sizeof(double)*(fluid->n_c));

#ifdef _USE_COLOR_GRAD_
    return;
#endif

    if (opt->skip_init_velocity_run && (opt->restart_run == 0)) {
        return;
    }

    //if (sys->flux <= 0.0) {
    //  return;
    //}

    if (opt->read_vel_from_file) {
        if (mpi->my_rank == 0)
            printf("\n**** READING VELOCITY FROM FILES 0000_init_f.bfd - %04d_init_f.bfd from folder %s **** \n", mpi->nr_procs - 1, out.c_str());
        sprintf(fn, "%04d_init_f%s", mpi->my_rank, opt->flip_x ? "__flip_x" : "");
        read_dist_binary(fn, fluid->f, fluid->n_c*sys->max_n_tot*sys->n_Q, sys);
        // read gx from file
        sprintf(fn, "%s/bfd/gx.dat", out.c_str());
        fp = my_fopen(fn, "r");
        fscanf(fp, "%lf", &(fluid->gravity[0]));
        fclose(fp);
        if (mpi->my_rank == 0) {
            printf("**** fluid->gravity[0] = %.5e ****\n", fluid->gravity[0]);
        }
        nsteps = 10;
    }

    if (opt->restart_run) {
        load_fluid_restart_files(fluid, sys);
        nsteps = 10;
    }

    if (nsteps > 0) {
        // advance a few steps to set things right
        for (i = 0; i < nsteps; i++) {
            advance_fluid(node, fluid, sys);
        }
        fluid->flux = get_flux(sys->node, sys, fluid);
        if (mpi->my_rank == 0) {
            printf("**** fluid->flux       = %.5e ****\n\n", fluid->flux);
        }
        sys->output->fluid.write(sys);

    } else {
        // start steady state velocity run
#ifndef REINIT_OFF // if velocity reinit is on
        initialize_velocity(sys->node, sys, fluid, 1);
#endif
    }
}


void init_3D_chemistry_run(struct System *sys, struct Field *species, struct Field *fluid, struct InitChem *InitChemSolver,
                           struct splayTree **st_lst, struct Minerals *minerals, struct Boundary *bndry)
/************************************
 ************************************/
{
    int cc;

    // init_run_3D_g(species, fluid, sys);
    // bc_init_run_g(&InitChemSolver[1], species, fluid, sys, minerals, bndry, st_lst);
    // //communicate(species->f, species->n_c, sys);
    // //communicate(species->f_tmp, species->n_c, sys);


    // if (sys->opt->find_isolated_nodes) {
    //   find_isolated_fluid_nodes(species, fluid, bndry, sys, InitChemSolver, st_lst, minerals);
    //   if (sys->opt->set_isolated_inert) {
    //     if (sys->mpi->my_rank==0) printf("**** Setting isolated nodes inert\n");
    //     set_isolated_nodes_inert(sys->node);
    //   }
    // }
    init_run_3D_g(species, fluid, sys);
    bc_init_run_g(&InitChemSolver[1], species, fluid, sys, minerals, bndry, st_lst);

    // init to steady-state and set init concentration
    //if (sys->species_not_converged) {
    switch (sys->opt->find_steady_state) {

    case 0:
        // no nothing
        break;

    case 1:
        // read from file, or iterate until convergence is met
        steady_state_specie_run(species, fluid, sys);

        if (sys->opt->kick_start) {
            if (sys->mpi->my_rank == 0) printf("**** Setting inlet concentration in whole sample\n");
            // set inlet concentraion in whole system
            set_init_bulk_conc(species->rho_inlet, species, fluid, sys);
        } else {
            // set equilibrium init concentration
            set_init_bulk_conc(species->rho_init, species, fluid, sys);
        }
        break;

    case 2:
        run_steady_state_all_species_and_fluid(species, fluid, minerals, sys);
        break;

    default:
        if (sys->mpi->my_rank == 0) printf("ERROR!! Unrecognized case in init_3D_chemistry_run(): %d\n", sys->opt->find_steady_state);
        break;
    }

    if (sys->opt->kick_start > 1) {  // diffusive run
        switch (sys->opt->kick_start) {
        case 2:
            // set INLET concentraion in whole system
            if (sys->mpi->my_rank == 0) printf("**** Setting INLET concentration in whole sample\n");
            set_init_bulk_conc(species->rho_inlet, species, fluid, sys);
            break;
        case 3:
            // set INITIAL concentraion in whole system
            if (sys->mpi->my_rank == 0) printf("**** Setting INITIAL concentration in whole sample\n");
            set_init_bulk_conc(species->rho_init, species, fluid, sys);
            break;
        default:
            if (sys->mpi->my_rank == 0)
                printf("**** ERROR non-valid option in init_3D_chemistry_run(): sys->opt->kick_start = %d\n",
                       sys->opt->kick_start);
            break;
        }
    }

    // get initial sum of species concentration
    get_fluid_conc(species, sys);
    for (cc = 0; cc < species->n_c; ++cc) {
        species->sum_conc_old[cc] = species->sum_conc_prev[cc] = species->sum_conc[cc];
    }

    // write mineral and species to file
    sys->output->mineral.write(sys);
    sys->output->chem.write(sys);
    
}

//-----------------------------------------------------------------
//
//
//
//-----------------------------------------------------------------
void run_steady_state_all_species_and_fluid(struct Field *species, struct Field *fluid, struct Minerals *minerals, struct System *sys)
{
    int t_max = 1e6;
    int t, cc, ci;
    double diff = -9.9e9;
    double conv_crit = sys->conv_crit_ss*sys->conv_length;
    char fn[50];
    double max_conc = -9e9;
    struct Node *node = sys->node;
    FILE *fp;

    // SPECIES
    for (cc = 0; cc < species->n_c; ++cc) {
        if (sys->opt->kick_start)
            species->rho_init[cc] = species->rho_inlet[cc];
        if (species->rho_init[cc] > max_conc) {
            max_conc = species->rho_init[cc];
        }
    }
    conv_crit *= max_conc;
    if (conv_crit < sys->conv_crit_ss_min) conv_crit = sys->conv_crit_ss_min;
    conv_crit *= sys->conv_length;

    // change INLET and OUTLET nodes to FLUID
    // inlet nodes
    for (ci = 0; ci < node->nr_in; ++ci) {
        node->mode[node->list_in[ci]] = FLUID;
    }
    // outlet nodes
    for (ci = 0; ci < node->nr_out; ++ci) {
        node->mode[node->list_out[ci]] = FLUID;
    }

    // SPECIES
    if (sys->opt->read_species_from_file) {
        // read steady state from file
        sprintf(fn, "%04d_init_ss_g", sys->mpi->my_rank);
        read_dist_binary(fn, species->f, species->n_c*sys->max_n_tot*sys->n_Q, sys);
    } else if (sys->opt->restart_run) {
        // read restart files
        load_chem_restart_files(species, sys);
    } else {
        // initialize new run
        init_run_3D_g(species, fluid, sys);
    }
    // boundary conditions
    periodic_bc(node, species, sys);
    wall_bc_bounce_back(node, species, sys);
    communicate(species->f, species->n_c, sys);
    communicate(species->f_tmp, species->n_c, sys);

    /* // FLUID */
    /* if (sys->opt->read_species_from_file) { */
    /*   // read fluid */
    /*   sprintf(fn, "%04d_init_ss_f", sys->mpi->my_rank); */
    /*   read_dist_binary(fn, fluid->f, fluid->n_c*sys->max_n_tot*sys->n_Q); */
    /*   // read gx from file */
    /*   sprintf(fn, "out/init_ss_f_gx.dat"); */
    /*   fp = my_fopen(fn, "r"); */
    /*   fscanf(fp, "%lf", &(fluid->gravity[0])); */
    /*   fclose(fp); */
    /* } else if (sys->opt->restart_run==0) { */
    /*   // init fluid run only if NOT restart run */
    /*   // (fluid is restarted in init_fluid_run()) */
    /*   init_run_f(node, fluid, sys); */
    /* } */
    /* communicate(fluid->f, fluid->n_c, sys); */
    /* communicate_node(node, sys); */
    /* periodic_bc(node, species, sys); */
    /* wall_bc_bounce_back(node, species, sys); */
    /* communicate(fluid->f, fluid->n_c, sys); */
    /* communicate(fluid->f_tmp, fluid->n_c, sys); */

    if (sys->opt->read_species_from_file) {
        sys->opt->read_species_from_file = 0;
        if (sys->mpi->my_rank == 0) printf("**** Reading steady state solution from file ");
        // advance a few steps to set things right
        t_max = 3;
    }

    if (sys->opt->restart_run) {
        //sys->opt->restart_run = 0;
        if (sys->mpi->my_rank == 0) printf("**** Restarting at t = %d", sys->start_itr);
        // advance a few steps to set things right
        t_max = 3;
    }

    if (sys->t < 1 && t_max>10 && sys->mpi->my_rank == 0) {
        printf("**** Initializing to %s concentration for all species (convergence < %.2e after %d steps)...",
               (sys->opt->kick_start) ? "INLET" : "EQUILIBRIUM", conv_crit, sys->conv_length);
        fflush(stdout);
    }
    // run until convergence for all species is met, gradually increase gx
    sys->gx = fluid->gravity[0];
    for (t = 0; t < t_max; ++t) {
        sys->t = t;
        if (t % 100 == 0 && sys->mpi->my_rank == 0) {
            printf("%5d\b\b\b\b\b", t);
            fflush(stdout);
        }
        //  if (t%100==0) {
        //  sys->nwrite = t;
        //  write_to_file(sys->output->chem, sys);
        //  write_to_file(sys->output->fluid, sys);
        //}
        //write_effluent_file(sys, species, fluid, NULL, 0);

        // gradually increase gx
        //if (t_max>10 && sys->t<200) {
        //  a = ((double)(sys->t+1)/200) - 1.0;
        //  fluid->gravity[0] = sys->gx*(1.0+a*a*a);
        //}

        propagate(node, fluid, sys);
        propagate(node, species, sys);
#ifdef NEWCONC
        collision_f_newconc(node, fluid, sys);
        collision_g_newconc(species, fluid, sys);
#else
        collision_f(node, fluid, sys);
        collision_g(species, fluid, sys);
#endif

        communicate(species->f, species->n_c, sys);
        bc_run_f(node, fluid, sys);
        periodic_bc(node, species, sys);
        wall_bc_bounce_back(node, species, sys);
        communicate(species->f, species->n_c, sys);
        communicate(species->f_tmp, species->n_c, sys);

        // check convergence
        if (t_max > 10 && species_has_converged(t, sys, species, conv_crit, &diff)) {
            sys->species_not_converged = 0;
#ifdef _DEBUG_
            if (sys->mpi->my_rank == 0) printf("\n*******  Finished steady state specie run\n\n");
#endif
            break;
        }

        // screendump steps and difference every 100 steps
        //#ifdef _DEBUG_CONVERGENCE_
        if (t % 100 == 0) {
            get_global_max(sys->mpi, &diff);
            if (sys->mpi->my_rank == 0) {
                if (t < 1) printf("\n");
                printf("\r**** %6d : diff = %.2e (%.2e)", t, diff, conv_crit);
                fflush(stdout);
            }
        }
        //#endif
    }
    if (sys->t < 1 && t_max>10 && sys->mpi->my_rank == 0) printf("took %d steps with criteria %.1e", t, conv_crit);

    if (sys->opt->read_species_from_file == 0) {
        // save steady-state solution to file unless read from file
        // species
        sprintf(fn, "%04d_init_ss_g", sys->mpi->my_rank);
        write_dist_binary(fn, species->f, species->n_c*sys->max_n_tot*sys->n_Q, sys);
        // fluid
        sprintf(fn, "%04d_init_ss_f", sys->mpi->my_rank);
        write_dist_binary(fn, fluid->f, fluid->n_c*sys->max_n_tot*sys->n_Q, sys);
        if (sys->mpi->my_rank == 0) {
            sprintf(fn, "%s/init_ss_f_gx.dat", sys->files["out"].c_str());
            fp = my_fopen(fn, "w");
            fprintf(fp, "%.8e\n", fluid->gravity[0]);
            fclose(fp);
        }

    }

#ifndef _PERIODIC_RUN_
    // change INLET and OUTLET nodes back
    // inlet nodes
    for (ci = 0; ci < node->nr_in; ++ci) {
        node->mode[node->list_in[ci]] = INLET;
    }
    // outlet nodes
    for (ci = 0; ci < node->nr_out; ++ci) {
        node->mode[node->list_out[ci]] = OUTLET;
    }
    communicate_node(node, sys);
#endif

    //sys->t = 0;
    //sys->nwrite = 0;
}

/* * * * * * * * * * * * *
 *                       *
 * * * * * * * * * * * * */
void init_solid_fraction(struct System *sys, struct Field *species,
                         struct Minerals *minerals, struct InitChem *ICS)
{
    struct Node *node = sys->node;
    int i, nn, cc, pos;
    struct Options *opt = sys->opt;
#ifdef _RANDOM_MINERAL_
    //                  calc  magn
    //double conc[2] = {0.85, 0.1};
    double conc[2] = { 1.0, 0.0 };
    double mtot = conc[0] + conc[1];
    int magn_per_10000 = (int)(sys->opt->magn_percent*100.0);
    int nw, nsteps, site_sum, solid_nodes, site_sum_mpi, solid_nodes_mpi;
#endif

    minerals->n_tot = ICS->size_sup_min;

    // allocate
    //minerals->list = (struct Minerals::Mineral *) calloc(minerals->n_tot, sizeof(struct Minerals::Mineral));
    minerals->list = new Minerals::Mineral[minerals->n_tot]();
    for (i = 0; i < minerals->n_tot; ++i) {
        strcpy(minerals->list[i].name, ICS->sup_min_name[i]);
        minerals->list[i].sfrac = (real *)calloc(sys->max_n_tot, sizeof(real));
        minerals->list[i].delta = (real *)calloc(sys->max_n_tot, sizeof(real));
        minerals->list[i].rest = (real *)calloc(sys->max_n_tot, sizeof(real));
        minerals->list[i].old_sfrac = (real *)calloc(sys->max_n_tot, sizeof(real));
        minerals->list[i].dm = (real *)calloc(sys->max_n_tot, sizeof(real));
        minerals->list[i].cm3_mol = ICS->SM_mineral->mol_volume[ICS->pos_sup_min[i]];
        minerals->list[i].SA = 0.0;
        minerals->list[i].SA_weight = 0.0;
    }
    minerals->SA_tot = 0.0;
    minerals->dm_sum = (real *)calloc(sys->max_n_tot, sizeof(real));
    minerals->sfrac_sum = (real *)calloc(sys->max_n_tot, sizeof(real));
#ifdef _MAGN_ON_MAGN_
    minerals->nucleation_site = (char *)malloc(sys->max_n_tot * sizeof(char));
    for (nn = 0; nn < sys->max_n_tot; ++nn) {
        minerals->nucleation_site[nn] = 0;
    }
#endif

    minerals->sum_sfrac = (double *)calloc(minerals->n_tot, sizeof(double));
    minerals->sum_sfrac_old = (double *)calloc(minerals->n_tot, sizeof(double));
    minerals->sum_sfrac_prev = (double *)calloc(minerals->n_tot, sizeof(double));
    minerals->buffer = (double *)calloc(minerals->n_tot, sizeof(double));

    // initialize arrays
    if (opt->restart_run) {
        // read from restart files
        load_mineral_restart_files(minerals, sys);
    } else {
        // set initial value from inp.dat
        for (nn = 0; nn < sys->max_n_tot; ++nn) { /* nodes */
            if (node->mode[nn] == SOLID || node->mode[nn] == WALL) {
                for (cc = 0; cc < minerals->n_tot; ++cc) {
#ifdef _RANDOM_MINERAL_
                    minerals->list[cc].sfrac[nn] = mtot*ICS->c_sup_min[cc];
#else
                    minerals->list[cc].sfrac[nn] = ICS->c_sup_min[cc];
#endif
                }
            }
        }
#ifdef _RANDOM_MINERAL_
        if (minerals->n_tot < 2) {
            std::cerr << "ERROR! _RANDOM_MINERAL_ requires at least 2 minerals, only "
                      << minerals->n_tot << " given!" << std::endl << std::endl;
            MPI_Finalize();
            exit(0);
        }
        site_sum = solid_nodes = 0;
        //srand((unsigned) time(&t));
        srand((unsigned)14097234168 + sys->mpi->my_rank);
        if (opt->magn_surface) {
            nsteps = node->nr_wall;
        } else {
            nsteps = sys->max_n_tot;
        }
        for (nw = 0; nw < nsteps; ++nw) {
            if (opt->magn_surface) {
                nn = node->list_wall[nw];
            } else {
                nn = nw;
            }
            if (node->is_ghost_or_periodic[nn])
                continue;

            if (node->mode[nn] == WALL || node->mode[nn] == SOLID) {
                solid_nodes++;
                if (rand() % 10000 < magn_per_10000) {
                    minerals->list[0].sfrac[nn] = conc[0];
                    minerals->list[1].sfrac[nn] = conc[1];
                    //minerals->list[0].sfrac[nn] = 0.0;
                    //minerals->list[1].sfrac[nn] = 0.999;
                    minerals->nucleation_site[nn] = 1;
                    site_sum++;
                    //// set neighbours recursively
                    ////ncalls = 0;
                    ////find_wall_neighbours(nn, nn, conc, _MAGN_RADIUS_, &ncalls, sys, sys->node, minerals);
                }
                //else {
                //  minerals->list[0].sfrac[nn] = 0.999;
                //  minerals->list[1].sfrac[nn] = 0.0;
                //}
            }
        }
        MPI_Reduce(&site_sum, &site_sum_mpi, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(&solid_nodes, &solid_nodes_mpi, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
        if (sys->mpi->my_rank == 0) {
            printf("\nFraction of nucleation nodes: %d/%d = %.2f%% (target value is %.2f %%)\n\n",
                   site_sum_mpi, solid_nodes_mpi, 100.0*site_sum_mpi / solid_nodes_mpi, (double)(magn_per_10000 / 100.0));
        }
#endif
    }


    minerals->init_molvol = 0.;
    for (cc = 0; cc < minerals->n_tot; ++cc) {
        pos = ICS->pos_sup_min[cc];
        minerals->init_molvol += ICS->c_sup_min[cc] * 1000. / ICS->SM_mineral->mol_volume[pos];
    }
    if (sys->mpi->my_rank == 0) {
        printf("Initial mineral concentration [mol/L]: %.5e\n\n", minerals->init_molvol);
//#ifdef _RANDOM_MINERAL_
//    if (opt->magn_surface) {
//      printf("Nucleation sites cover %.2f %% of the SURFACE\n\n", (double)(magn_per_10000 / 100.0));
//    } else {
//      printf("Nucleation sites cover %.2f %% of the SOLID VOLUME\n\n", (double)(magn_per_10000 / 100.0));
//    }
//#endif
    }
    // allocate arrays to hold sfrac updates associated with links to ghost-nodes
//  sys->mpi->ghost_sfrac_diff = (real *)malloc(node->nr_ghost*minerals->n_tot * sizeof(real));
//  sys->mpi->ghost_sfrac_tmp = (real *)malloc(node->nr_ghost*minerals->n_tot * sizeof(real));

    communicate_sfrac(sys, minerals);

}


#ifdef _RANDOM_MINERAL_
//---------------------------------------------
//
//
//---------------------------------------------
void find_wall_neighbours(int nn0, int nn, double *conc, double rad,
                          int *ncalls, struct System *sys, struct Node *node, struct Minerals *minerals)
{
    int ck, nb;
    double d1, d2;
    //static int ncall = 0;
    (*ncalls)++;
    for (ck = 1; ck < sys->n_Q; ++ck) {
        nb = nn - sys->step->n_rel_neig[ck];
        if (nb >= 0 && nb < sys->max_n_tot && !node->is_ghost[nb] && node->mode[nb] == WALL) {
            minerals->list[0].sfrac[nb] = conc[0];
            minerals->list[1].sfrac[nb] = conc[1];
            minerals->nucleation_site[nb] = 1;
            d1 = get_dist_sqrd(nn0, nn, sys);
            d2 = get_dist_sqrd(nn0, nb, sys);
            if ((d2 > d1) && (d2 < (rad*rad) + 1e-14) && ((*ncalls) < 1e7)) {
                find_wall_neighbours(nn0, nb, conc, rad, ncalls, sys, sys->node, minerals);
            }
        }
    }

}
#endif

//---------------------------------------------
//
//---------------------------------------------
double get_dist_sqrd(int n1, int n2, struct System *sys)
{
    int i, ind1[3], ind2[3];
    double dist = 0;
    get_local_coord(n1, ind1, sys);
    get_local_coord(n2, ind2, sys);
    for (i = 0; i < 3; ++i) {
        ind1[i] -= ind2[i];
        dist += (double)(ind1[i] * ind1[i]);
    }
    return(dist);
}

//---------------------------------------------
//
//---------------------------------------------
//void find_isolated_fluid_nodes(struct Field *species, struct Field *fluid, struct Boundary *bndry, struct System *sys, struct InitChem *InitChemSolver,
//			       struct splayTree **st_lst, struct Minerals *minerals)
void find_isolated_fluid_nodes(struct Field *species, struct Field *fluid, struct System *sys)
{
    double rho_init = 1.0, rho_inlet = 1e10;
    double crit = (1.0 + 1e-8)*rho_init;
    int t_max = 1e6, check_int = 100;

    struct Node *node = sys->node;
    //double *u;
    double *conc;

    // only use phase 0
#ifdef NEWCONC
    conc = &species->psi[0];
#else
    conc = &species->rho[0];
#endif

    // init
    // only use specie 0 during s-s run
    species->ntot = species->n_c;
    species->n_c = 1;
    species->rho_init0 = species->rho_init[0];
    species->rho_init[0] = rho_init;
    species->rho_inlet0 = species->rho_inlet[0];
    species->rho_inlet[0] = rho_inlet;

    std::vector< std::vector<int> > fluid_nodes(2, std::vector<int>());
    fluid_nodes[0].clear();
    fluid_nodes[1].clear();
    node->isolated.clear();
    node->perc_cluster.clear();
    for (int a = 0; a < 2; ++a) {  // loop over positive and negative x-direction (flooding direction)
        // choose direction
        if (a == 0) {
            node->perc_inlet = sys->node->in_x + sys->opt->closed_cell;
        } else {
            node->perc_inlet = sys->node->out_x - sys->opt->closed_cell;
        }

        if (sys->mpi->my_rank == 0) {
            //printf("**** Looking for isolated fluid nodes in %s x-direction ... (crit = %.10e)", (a == 0) ? "positive" : "negative", crit);
            printf("**** Looking for isolated fluid nodes in %s x-direction ... ", (a == 0) ? "positive" : "negative");
            fflush(stdout);
        }

        // find fluid nodes
        for (int cn = 0; cn < sys->max_n_tot; ++cn) { /* node */
            if ( (node->mode[cn]==FLUID) && !node->is_ghost[cn]) {
                fluid_nodes[0].push_back(cn);
            }
        }
        int nfluid_tot = fluid_nodes[0].size();
        broadcast_sum_of_int_to_all_procs(&nfluid_tot, sys->mpi);

        init_run_3D_g(species, fluid, sys);
        periodic_bc(sys->node, species, sys);
        wall_bc_bounce_back(sys->node, species, sys);
        //bc_init_run_g(&InitChemSolver[1], species, fluid, sys, minerals, bndry, st_lst);
        wall_bc_bounce_back(node, species, sys);

        int ncheck = 0, nneg = 0;
        int not_converged = 1;

        for (sys->t = 0; (sys->t < t_max) && (not_converged > 0); ++sys->t) {
            propagate(node, species, sys);
#ifdef NEWCONC
            collision_g_newconc(species, fluid, sys);
#else
            collision_g(species, fluid, sys);
#endif
            communicate(species->f, species->n_c, sys);
            wall_bc_bounce_back(node, species, sys);
            communicate(species->f, species->n_c, sys);

            if (sys->t > 0 && sys->t%check_int == 0) {
                broadcast_sum_of_int_to_all_procs(&nneg, sys->mpi);
                if (sys->mpi->my_rank == 0) {
                    printf("%7d    (%.1e)\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b", sys->t, (float)nneg);
                    fflush(stdout);
                }
                nneg = 0;
                not_converged = 1;
                int n0 = ncheck % 2;
                int n1 = (ncheck + 1) % 2;
                fluid_nodes[n1].clear();
                for (std::vector<int>::iterator it = fluid_nodes[n0].begin(); it != fluid_nodes[n0].end(); ++it) {
                    int fn = *it; //fluid_nodes[n0][*it];
                    if (conc[fn] < crit) {
                        fluid_nodes[n1].push_back(fn);
                    } else {
                        node->perc_cluster.push_back(fn);
                    }
                    if (conc[fn] < 0.0) {
                        nneg++;
                    }
                }
                if ((fluid_nodes[0].size()==fluid_nodes[1].size()) && (nneg==0)) {
                    not_converged = 0;
                }
                broadcast_sum_of_int_to_all_procs(&not_converged, sys->mpi);
                ncheck++;
            }

            // if (a==1 && sys->t%1==0) {
            // 	//printf("\r%d",t);fflush(stdout);
            // 	sys->nwrite++;
            // 	write_to_file(sys->output->chem, sys);
            // }
        }

        // allocate and init lists
        for (std::vector<int>::iterator it = fluid_nodes[1].begin(); it != fluid_nodes[1].end(); ++it) {
            //node->isolated.push_back(fluid_nodes[1][i]);
            node->isolated.push_back(*it);
        }

        // output
        //int nfluid = fluid_nodes[0].size();
        int nfluid = node->isolated.size();// fluid_nodes[0].size();
        int nperc = node->perc_cluster.size();
        broadcast_sum_of_int_to_all_procs(&nfluid, sys->mpi);
        broadcast_sum_of_int_to_all_procs(&nperc, sys->mpi);
        if (sys->mpi->my_rank == 0) {
            printf("\nFound %d non-percolating nodes (%3.2f %%) in %d steps                             \n",
                   nfluid, 100.0*nfluid/((sys->MAX_N[0]-2)*(sys->MAX_N[1]-2)*(sys->MAX_N[2]-2)), sys->t - 1);
            if (nfluid + nperc != nfluid_tot) {
                printf("WARNING! Sum of non- and percolating nodes different from fluid nodes: %d - %d = %d\n",
                       nfluid_tot, nfluid + nperc, nfluid_tot - (nfluid + nperc));
            }
            fflush(stdout);
        }

        //for (std::vector<int>::iterator it = node->isolated.begin(); it != node->isolated.end(); ++it) {
        for (std::vector<int>::iterator it = fluid_nodes[0].begin(); it != fluid_nodes[0].end(); ++it) {
            node->mode[*it] = SOLID;
        }
        //    for (cn = 0; cn < sys->max_n_tot; ++cn) {
        //      if (get_global_ycoord(cn, sys)==20) {
        //        node->mode[cn] = SOLID;
        //      }
        //    }
        communicate_node(node, sys);
        set_wall_nodes(node, sys);
        sys->output->write_geo_file((a == 0) ? "geo_perc_in" : "geo_perc_final", node, sys);
        //write_geo_file("geo_perc_c", node, sys->output, sys);
        if (sys->mpi->my_rank == 0) {
            printf("Percolating cluster (from %s) saved in geo_perc_%s.pvti\n\n", (a == 0) ? "inlet" : "outlet", (a == 0) ? "in" : "final");
        }

        fluid_nodes[0].clear();
        fluid_nodes[1].clear();
        node->isolated.clear();
        node->perc_cluster.clear();
    }

    // reset to original values
    sys->t = 0;
    species->n_c = species->ntot;
    species->rho_init[0] = species->rho_init0;
    species->rho_inlet[0] = species->rho_inlet0;
}



//---------------------------------------------
//
//
//---------------------------------------------
void set_isolated_nodes_inert(struct Node *node)
{
    int nl;//, c=0;
    struct Links::Link *link;

    if (node->isolated.size() == 0)
        return;

    // set isolated nodes inert
    for (nl = 0; nl < node->links->nlinks; ++nl) { /* loop links */
        link = &node->links->list[nl];
        //for (i=0; i<node->nr_isolated; ++i) {
        for (std::vector<int>::size_type i = 0; i < node->isolated.size(); ++i) {
            if (link->fluid == node->isolated[i]) {
                //printf("%d: link %d with node %d set inert\n", ++c, nl, link->fluid);
                link->inert = 1;
                break;
            }
        }
    }
}


//---------------------------------------------
//
//  Only for specie 0: Set rho to 1 and loop until steady state
//                     Temporarily store original rho-field in rho_copy
//
//---------------------------------------------
void steady_state_specie_run(struct Field *species, struct Field *fluid, struct System *sys)
{
    struct Node *node = sys->node;
    int t_max = 1e6;
    int t;
    //double *rho_copy;
    double diff = 0.0;
    //double conv_crit = sys->conv_crit_chem * sys->conv_length;
    double conv_crit = sys->conv_crit_ss*sys->conv_length;
    char fn[50];

    //rho_copy = malloc(sys->max_n_tot*sizeof(double));

    // backup original rho-field for specie 0
    //for( cn = 0; cn < sys->max_n_tot; ++cn ) { // nodes
    //  rho_copy[cn] = species->rho[cn]; // only phase 0
    //}

    // backup variables, set INLET and OUTLET nodes to FLUID, call init routines
    init_ss_species_run(species, fluid, sys);

    // read ss solution from file?
    if (sys->opt->read_species_from_file) {
        // only read from file during initialization
        sys->opt->read_species_from_file = 0;
        if (sys->mpi->my_rank == 0) printf("**** Reading steady state specie solution from files ");
        // read species
        sprintf(fn, "%04d_init_ss_g", sys->mpi->my_rank);
        read_dist_binary(fn, species->f, species->n_c*sys->max_n_tot*sys->n_Q, sys);
        // advance a few steps to set things right
        t_max = 10;
    }

    if (sys->t < 1 && t_max>10 && sys->mpi->my_rank == 0) {
        printf("**** Finding steady state specie solution ... ");
        fflush(stdout);
    }
    // time-loop for steady-state run of 1 specie initialized to 1
    // apply periodic bc and bounce back
    for (t = 0; t < t_max; ++t) {
        sys->t = t;
        //write_effluent_file(sys, species, fluid, NULL, 0);
        if (t % 100 == 0 && sys->mpi->my_rank == 0) {
            printf("%5d\b\b\b\b\b", t);
            fflush(stdout);
        }
        //if (t%100==0) {
        //  sys->nwrite++;
        //  write_to_file(sys->output->chem, sys);
        //}

        propagate(node, species, sys);
#ifdef NEWCONC
        collision_g_newconc(species, fluid, sys);
#else
        collision_g(species, fluid, sys);
#endif

        communicate(species->f, species->n_c, sys);
        periodic_bc(node, species, sys);
        wall_bc_bounce_back(node, species, sys);
        communicate(species->f, species->n_c, sys);

        // check convergence
        if (t_max > 10 && species_has_converged(t, sys, species, conv_crit, &diff)) {
            sys->species_not_converged = 0;
#ifdef _DEBUG_
            if (sys->mpi->my_rank == 0) printf("\n*******  Finished steady state specie run\n\n");
#endif
            break;
        }

        // screendump steps and difference every 100 steps
#ifdef _DEBUG_CONVERGENCE_
        if (t % 100 == 0) {
            get_global_max(sys->mpi, &diff);
            if (sys->mpi->my_rank == 0) {
                if (t < 1) printf("\n");
                printf("\r**** %6d : diff = %.2e (%.2e)", t, diff, conv_crit);
                fflush(stdout);
            }
        }
#endif
    }
    if (sys->t < 1 && t_max>10 && sys->mpi->my_rank == 0) printf("took %d steps with criteria %.1e", t, conv_crit);

    if (sys->opt->read_species_from_file == 0) {
        // save steady-state solution to file unless read from file
        // species
        sprintf(fn, "%04d_init_ss_g", sys->mpi->my_rank);
        write_dist_binary(fn, species->f, species->n_c*sys->max_n_tot*sys->n_Q, sys);
    }

    // restore variables and set INLET and OUTLET nodes
    end_ss_species_run(species, sys);

    // restore original rho-field for specie 0
    //for( cn = 0; cn < sys->max_n_tot; ++cn ) { // nodes
    //  species->rho[cn] = rho_copy[cn]; // only phase 0
    //}
    //sys->t = 0;
    //sys->nwrite = 0;
    //free(rho_copy);
}



void init_ss_species_run(struct Field *species, struct Field *fluid, struct System *sys)
/////////////////////////////////////////////////////////////////////////////
// use periodic bc until steady state. Remember to change INLET and
// OUTLET nodes back to FLUID nodes, and vice versa when s-s is reached.
// Only use one chem phase (remember to change back after s-s)
/////////////////////////////////////////////////////////////////////////////
{
    struct Node *node = sys->node;
    int ci;
    // only use specie 0 during s-s run
    species->ntot = species->n_c;
    species->n_c = 1;
    // set init values to 1
    species->rho_init0 = species->rho_init[0];
    species->rho_init[0] = 1.0;
    // change nodes
    // loop inlet nodes...
    for (ci = 0; ci < node->nr_in; ++ci) {
        node->mode[node->list_in[ci]] = FLUID;
    }
    // ...and outlet nodes
    for (ci = 0; ci < node->nr_out; ++ci) {
        node->mode[node->list_out[ci]] = FLUID;
    }
    init_run_3D_g(species, fluid, sys);
    // init boundary conditions
    periodic_bc(node, species, sys);
    wall_bc_bounce_back(node, species, sys);
}


void end_ss_species_run(struct Field *species, struct System *sys)
/////////////////////////////////////////////////////////////////////////////
// steady-state for specie run is reached, and all species are set to
// init values. The inlet- and outlet-nodes are changed back to INLET and OUTLET
// (they were FLUID during s-s run due to periodic boundary conditions)
/////////////////////////////////////////////////////////////////////////////
{
    int ci;
    struct Node *node = sys->node;

    species->n_c = species->ntot;
    species->rho_init[0] = species->rho_init0;
    // change nodes back to INLET and OUTLET
    // loop inlet nodes...
    for (ci = 0; ci < node->nr_in; ++ci) {
        node->mode[node->list_in[ci]] = INLET;
    }
    // ...and outlet nodes
    for (ci = 0; ci < node->nr_out; ++ci) {
        node->mode[node->list_out[ci]] = OUTLET;
    }
}
