#include "biozementfunctions.h"

void init_run_3D_dif(struct Field *species, struct Field *fluid, struct System *sys)
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
#ifdef NEWCONC
            if (node->mode[cn] == FLUID || node->mode[cn] == INLET || node->mode[cn] == OUTLET) { /* Inside the fluid */
                if (get_global_xcoord(cn, sys) == 5) {
                    psi[rho_ind] = 1;
                } else {
                    if (node->Bmask[cn] == 1)
                        psi[rho_ind] = rho_init;
                    else
                        psi[rho_ind] = 0;
                }


#ifdef FLUID_OFF
                rho[rho_ind] = psi[rho_ind];
#else
                rho[rho_ind] = psi[rho_ind] * fluid->rho[rho_ind];
#endif
            } else {
                psi[rho_ind] = rho[rho_ind] = 0.0;
            }
#else
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


void collision_g_newconc_dif_calc_rho(struct Field *species, struct System *sys)
/*
 * Calculates rho = \sum_\alpha g_\alpha for all nodes and phases.
 *
 * After this function is run
 *   species->rho = \sum_\alpha species->t_tmp[\alpha]
 */
{
    int g_ind = 0; /* index in for the velocity and velocity distribution */
//    real *g_tmp, *rho; /* Pointers to the distribution, temprorar distribution, density and velocity */
    struct Step *step = sys->step;
    struct Node *node = sys->node;


    for(int cc = 0; cc < species->n_c; ++cc ) {
        real * g_tmp = &species->f_tmp[step->f_phase*cc];
        real * rho   = &species->rho[step->rho_phase*cc];


        /* for each !FLUID! node DO*/
        g_ind = 0;
        for(int cn = 0; cn < sys->max_n_tot; ++cn) { /* nodes */
            if( node->mode[cn]==FLUID) { /* part of the fluid */
                //if( node->mode[cn]==FLUID || node->mode[cn]==INLET || node->mode[cn]==OUTLET ) { // NBNBNBNBNB

                /*    calculate rho */
                rho[cn] = 0;
                for(int ck = 0; ck < sys->n_Q; ++ck ) {
                    rho[cn] += g_tmp[g_ind+ck];
                }

            }
            g_ind += step->f_node;
        }
    }
}


void collision_g_newconc_dif_source(struct Field *dif, struct Field *species, struct Field *fluid, struct System *sys, struct InitChem *ICS)
/*
 * Here we will use our kinetic relations to obtain the source terms.
 * Remember that:
 *  1) R_cc = (d/dt)c
 *  2) c    = \sum_\alpha g_\alpha + (1/2)R_cc
 *
 *  In 'collision_g_newconc_dif_calc_rho' we have already calculated the sums g and
 *  saved them in 'species->rho'
 *
 * After this function is run...
 *  - species->R_local is set
 *  - species->rho[output] = species->rho[input] + (1/2)*species->R_local
 */
{
    struct Step *step = sys->step;
    struct Node *node = sys->node;

    real yield = (*ICS).BiP->yield; 
    real mu_max = (*ICS).BiP->mu_max;
    real K_u = (*ICS).BiP->K_i;

    for( int cn = 0; cn < sys->max_n_tot; ++cn) { /* nodes */
        if( node->mode[cn]==FLUID) { /* part of the fluid */


            //  *** Initiate all rates to zero ***
            for( int cc = 0; cc < species->n_c; ++cc ) {
                int nodeNo = step->rho_phase*cc + cn;
                species->R_local[nodeNo] += 0;

            } // For each phase
        int D_ratio = round(sys->D / sys->D_b); // /* AMIN */ 
            if (sys->t % D_ratio == 0 ) {
                for (int cc = 0; cc < dif->n_c; ++cc) {
                    int nodeNo = step->rho_phase*cc + cn;
                    dif->R_local[nodeNo] = 0;
                } // For each phase
            }

            // END end initate



            // ******************************* ADD THE RATE LAW
            
            // ******************************* END ADD RATE LAW
            
            

            // *** Update consentrations

            int nodeNo_dif = step->rho_phase * 0 + cn; // node ID for biomass field /* AMIN */
            //for( int cc = 0; cc < species->n_c; ++cc ) {
             //   int nodeNo = step->rho_phase * cc + cn;
                int nodeNo = step->rho_phase * 3 + cn;
                if (node->Bmask[cn] == 1){ // for urea
                    real delta = pow((2*K_u-2*species->rho[nodeNo] - mu_max*dif->rho[nodeNo_dif]*(1+yield)), 2) + 8*species->rho[nodeNo]*K_u*(2+yield*mu_max);
                    real rho_urea_old = species->rho[nodeNo];
                   // species->rho[nodeNo] = (-2*K_u+2*species->rho[nodeNo] + mu_max*dif->rho[nodeNo_dif]*(1+yield)+sqrt(delta))/(4+2*yield*mu_max);
                    species->R_local[step->rho_phase*3 + cn] = -.1 * species->rho[nodeNo]; // for Urea
                    dif->R_local[nodeNo_dif] += -yield * species->R_local[step->rho_phase*3 + cn];
                    species->R_local[step->rho_phase + cn] = - species->R_local[step->rho_phase*3 + cn]; // for HCO3   
                    species->R_local[step->rho_phase*2 + cn] = -2*species->R_local[step->rho_phase*3 + cn]; // for NH3   
                    //std::cout<<"species->rho[step->rho_phase + cn] = "<< species->rho[step->rho_phase + cn] <<std::endl;
                }
                // species->rho[nodeNo] += 0.5 * species->R_local[nodeNo];
             // For each phase
            //for( int cc = 0; cc < dif->n_c; ++cc ) {
            //    int nodeNo = step->rho_phase*cc + cn;
            //    dif->rho[nodeNo] += 0.5*dif->R_local[nodeNo];
            //} // For each phase
            //** END update consentrations

        } // For each fluid node
    }
}

void collision_g_newconc_dif_collide(struct Field *species, struct Field *fluid, struct System *sys)
/*
 * Performs a collision of the species field.
 * We assume that species->rho and species->R_local are already calcualted.
 */
{
    struct Step *step = sys->step;
    real * u = fluid->u;

    for (int cc = 0; cc < species->n_c; ++cc) {
        real * g     = &species->f[step->f_phase*cc];
        real * g_tmp = &species->f_tmp[step->f_phase*cc];
        real * rho   = &species->rho[step->rho_phase*cc];
        real * psi   = &species->psi[step->rho_phase*cc];
        real * R_cc  = &species->R_local[step->rho_phase*cc]; /* AMIN EJE */

        real tau_inv = 1./species->tau[cc];

        int u_ind = 0; /* index in for the velocity and velocity distribution */
        int g_ind = 0;

        double A = (1.0 - 0.5*tau_inv); // /sys->cs_2;
        /* for each !FLUID! node DO*/
        for (int cn = 0; cn < sys->max_n_tot; ++cn) { /* nodes */
            if( sys->node->mode[cn]==FLUID) { /* part of the fluid */

                /*    calculate psi */
#ifdef FLUID_OFF
                psi[cn] = rho[cn];
#else
                psi[cn] = rho[cn]/fluid->rho[cn];
#endif

                // Collision
                for (int ck = 0; ck < sys->n_Q; ++ck) {

                    /*    calculate the force corrections */
                    /*    calculate the equilibrium distribution (see fluid field) */
#ifdef FLUID_OFF
                    real geq = calc_geq(ck, rho[cn], &u[u_ind], sys);
#else // fluid is on
                    real geq = calc_feq(ck, rho[cn], &u[u_ind], sys);
#endif
                    /*    update g*/
                    /* AMIN EJE: set the source */
                    real delta_src = sys->w[ck] * A * R_cc[cn];
                    // Propagation
                    g[g_ind + ck] = g_tmp[g_ind + ck] - (g_tmp[g_ind + ck] - geq)*tau_inv + delta_src;
                }
            }

            u_ind += step->u_node;
            g_ind += step->f_node;
        } // For each node
    } // For each phase
}


/* ---------------------------------------------------------------------------------  collision_g */
void collision_g_newconc_dif(struct Field *dif, struct Field *species, struct Field *fluid, struct System *sys, struct InitChem *ICS)
/*
 * collision_g :
 *
 *   Need to add source terms dependent on bacteria consentration, glucose, ....
 *          and the chemical species.
 *
 *   We need to first calculate the consentrations -> and then collide and propagate.
 *
 */
{
    // Find the zeroth moment
    collision_g_newconc_dif_calc_rho(species, sys);
    // density update rate for dif should be different from diffusive field! /* AMIN */
    int D_ratio = round(sys->D / sys->D_b); // /* AMIN */ 
    if (sys->t % D_ratio == 0 & sys->t != 0)
        collision_g_newconc_dif_calc_rho(dif, sys); 
    //for (int cn = 0; cn < sys->max_n_tot; ++cn)

        // Set the source term
        collision_g_newconc_dif_source(dif, species, fluid, sys, ICS);

    // Collition and propagation
    collision_g_newconc_dif_collide(species, fluid, sys);
    if (sys->t % D_ratio == 0  & sys->t != 0)
        collision_g_newconc_dif_collide(dif, fluid, sys);
}




void applyFakeChem(struct Links::Link *link, struct Field *species, struct Step *step, struct System *sys)
{
    int dir = link->dir; /* direction from solid to fluid */
    int dir_bb = sys->k_bb[dir]; /* direction from fluid to solid */

    double *g_in, *g_out;
    double c_eq = 0.1, k=0.1, c2_inv = 3.0;
    double kc2 = k*c2_inv;
    g_in = &(species->f[step->f_node*link->fluid + dir_bb]);
    g_out = &(species->f[step->f_node*link->wall + dir]);
    for (int cc = 0; cc < species->n_c; ++cc) {
        g_out[step->f_phase*cc] = ( g_in[step->f_phase*cc]*(1.0-kc2) + kc2*2.0*sys->w[dir]*c_eq )/(1.0+kc2);
    }
}

void printRho(struct Field *species, struct System *sys, FILE *outputRho)
{
    if ((sys->t % 1) == 0) {
        int numFluidNodes = 0;
        double rhoCaSum = 0.0;
        double rhoHCO3Sum = 0.0;
        double rhoLacSum = 0.0;
        for (int cn = 0; cn < sys->max_n_tot; ++cn) {
            if (sys->node->mode[cn] == FLUID) {
                numFluidNodes += 1;
                //std::cout << numFluidNodes << std::endl; /* AMIN */
                rhoCaSum += species->rho[cn + 0 * sys->step->rho_phase];
                rhoHCO3Sum += species->rho[cn + 1 * sys->step->rho_phase];
                rhoLacSum += species->rho[cn + 2 * sys->step->rho_phase];
            }
        }
        /* Write to file */
        fprintf(outputRho, "%i \t %lf \t %lf \t %lf \t %d\n",sys->t, rhoCaSum / numFluidNodes, rhoHCO3Sum / numFluidNodes, rhoLacSum / numFluidNodes, numFluidNodes);
    }
}

void printChem(struct Field *species, struct System *sys, FILE *ChemOutput)
{
    double rhoLac;
    for (int cn = 0; cn < sys->max_n_tot; ++cn) {
        rhoLac =  species->rho[cn + 2 * sys->step->rho_phase];
        fprintf(ChemOutput, "%lf \t", rhoLac);
        if ( (cn+1) % sys->max_n[0] == 0)
            fprintf(ChemOutput, "\n");
    }
    fprintf(ChemOutput, "\n");
}

void read_Rfile(struct Node *node, struct System *sys)
{
    FILE *fpR;
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
    fpR = my_fopen(sys->files["Rfile"].c_str(), "r"); // open file

    skip_header_line(fpR, sys->mpi);

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

                if (fscanf(fpR, "%1d", &tmp) == EOF) {
                    eof_error(sys->mpi);
                }

                xr = x;
                yr = y;
                zr = z;

                if (opt->rotate_y) {
                    xr = z;
                    //zr = N[0]-1-x;
                    zr = geodim[0] - x;
                }
                if (opt->rotate_z) {
                    xr = y;
                    yr = geodim[0] - x;
                }

                Nx = xr;
                Ny = yr;
                Nz = zr;

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

                if( tmp == 1)
                    sys->node->Rmask[n_ind] = 1;
                else {
                    sys->node->Rmask[n_ind] = 0;
                }
            }
            if (fscanf(fpR, "\n") == EOF) {
                eof_error(sys->mpi);
            }
        }
        if (fscanf(fpR, "\n") == EOF) {
            eof_error(sys->mpi);
        }
        ++z;
    } while (z < geodim[2]);

    fclose(fpR);
    // update ghost-nodes and set periodic nodes (for 3D runs)
    communicate_node(node, sys);
}


void update_pH_buffer(InitChem &ICS_pH, Outfile* chem, struct Field *species, struct System *sys)
{
    struct Step *step = sys->step;
    struct Node *node = sys->node;
    real c_vals[species->n_c];

    for (int cn = 0; cn<sys->max_n_tot; cn++) {
        if (node->mode[cn] == FLUID) {
            // Asign the basis specie consentrations
            for (int cc = 0; cc < species->n_c; cc++)
                c_vals[cc] = species->rho[step->rho_phase*cc + cn];
            chem->buffer[cn] = get_pH(ICS_pH, c_vals);
        } else {
            chem->buffer[cn] = 0.0;
        }
    }

}

void read_Bfile(struct Node *node, struct System *sys)
{
    FILE *fpR;
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
    fpR = my_fopen(sys->files["Bfile"].c_str(), "r"); // open file

    skip_header_line(fpR, sys->mpi);

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

                if (fscanf(fpR, "%1d", &tmp) == EOF) {
                    eof_error(sys->mpi);
                }

                xr = x;
                yr = y;
                zr = z;

                if (opt->rotate_y) {
                    xr = z;
                    //zr = N[0]-1-x;
                    zr = geodim[0] - x;
                }
                if (opt->rotate_z) {
                    xr = y;
                    yr = geodim[0] - x;
                }

                Nx = xr;
                Ny = yr;
                Nz = zr;

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

                if( tmp == 1)
                    sys->node->Bmask[n_ind] = 1;
                else {
                    sys->node->Bmask[n_ind] = 0;
                }
            }
            if (fscanf(fpR, "\n") == EOF) {
                eof_error(sys->mpi);
            }
        }
        if (fscanf(fpR, "\n") == EOF) {
            eof_error(sys->mpi);
        }
        ++z;
    } while (z < geodim[2]);

    fclose(fpR);
    // update ghost-nodes and set periodic nodes (for 3D runs)
    communicate_node(node, sys);
}

void init_pH_cal(InitChem *ICS, struct Field *species, Input &input, struct System *sys)
{
//    ICS = new InitChem();
    // Set inlet values
    InitICS(ICS, species->n_c);
    // Temperature in kelvin
    ICS->Temp = double(input["chem"]["temp"]) + 273.15;
    // Species names
    std::vector<std::string> name = input["chem"]["bulk"].get_names();
    for (int i = 0; i < species->n_c; ++i) {
        strcpy(ICS->species_name[i], name[i].c_str());
        ICS->c_vchem[i] = input["chem"]["inlet"][name[i]];
        ICS->pos_mass[ICS->size_mass] = i;
        ICS->size_mass++;
    }

    ICS->INTERPOLATE = sys->opt->splaytree_ON;
    ICS->num_chrg_not_conv = 0;

    init_chemistry(ICS, sys->files["data"]);
    chem_get_pos(ICS);
    db_set_activity(ICS);
}

real get_pH(InitChem &ICS, real * c_vchem_eq)
{
    for (int i = 0; i < ICS.size; ++i)
        ICS.c_vchem[i] = c_vchem_eq[i];


    struct BasVec V;
//     int pos;
//     real c_h;


    make_basvec_struct(&ICS, &V);

    trans_nlin_eq(&V); /* add rate minerals to equilibrium minerals */

    V.equilibrate = 1;
    V.chbal = 1;
    for(int c=0; c<ICS.size; ++c) c_vchem_eq[c] = ICS.c_vchem[c];

    calc_surface_reaction_nlin(c_vchem_eq,1,&V);
//     if(ICS.PRINT_DEBUG_CHEM) writeVchem(&V,"qqqV.out");
//     for(int c = 0; c < V.ICS->size; ++c)
//       {
//         pos = ICS.pos[c];
//         ICS.c_ads[c] = V.ctot_ads_calc[pos];
//       }


    /*NB!!!!!!!!!!!!!!!*/
    /*REMEMBER TO FREE STRUCT */
    if(V.Vsurf != NULL) free_BasVec(V.Vsurf);
    free_BasVec(&V);
    return -V.log_a[V.pos_pH];
}

void set_nodes_inert(struct Node *node)
{
    int nl;//, c=0;
    struct Links::Link *link;


    // set isolated nodes inert
    for (nl = 0; nl < node->links->nlinks; ++nl) { /* loop links */
        link = &node->links->list[nl];
        if (node->inert[link->wall]) {
            link->inert = 1;
        }
    }
}

void temporal_rate(struct InitChem *ICS_ptr, struct BasVec **Vchem, struct System *sys, Input &inp)
{
    double k_neut = inp["chem"]["rate"]["calcite"][3]*sys->dt/(1e3*sys->dx);
    double k_acid = inp["chem"]["rate"]["calcite"][4]*sys->dt/(1e3*sys->dx);
    double t_max = inp["chem"]["rate"]["calcite"][8]; // t_max in seconds
    double n = inp["chem"]["rate"]["calcite"][7]; // k_max = n * k_0
    (*Vchem)[1].rate[0][0] = n*k_neut*exp(-sys->time/t_max)+k_neut; // exponential rate decrease
    (*Vchem)[1].rate[0][1] = n*k_acid*exp(-sys->time/t_max)+k_acid;
    //(*Vchem)[1].rate[0][0] = n*k_neut - k_neut * (n-1)/t_max * sys->time; // linear rate decrease
    //(*Vchem)[1].rate[0][1] = n*k_acid - k_acid * (n-1)/t_max * sys->time;
    //ICS_ptr->rate[0][0] = 1.e-8;//n*k_neut - k_neut * (n-1)/t_max * sys->time;
    //ICS_ptr->rate[0][1] = 1.e-8;//n*k_acid - k_acid * (n-1)/t_max * sys->time;
}
void initial_temporal_rate(struct InitChem **ICS_ptr, struct System *sys, Input &inp)
{
    double k_neut = inp["chem"]["rate"]["calcite"][3];
    double k_acid = inp["chem"]["rate"]["calcite"][4];
    double n = inp["chem"]["rate"]["calcite"][7]; //k = n * k_0
    (*ICS_ptr)[1].rate[0][0] = n*k_neut;
    (*ICS_ptr)[1].rate[0][1] = n*k_acid;

}
