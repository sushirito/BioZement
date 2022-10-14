
#include "global.h"
#include "chem_global.h"
#include "output.h"
#include "biozementfunctions.h"

int main(int argc, char *argv[])
{
  // ----------------------------
  //  variable declarations
  // ----------------------------
  Input options, input;
  Options opt;
  Step *step = nullptr;
  Minerals *minerals = nullptr;
  Links *links = nullptr;
  Node *node = nullptr;
  System *sys = nullptr;
  Field *fluid = nullptr;
  Field *species = nullptr;
  Field *dif_field = nullptr;
  Boundary *bndry = nullptr;
  InitChem *InitChemSolver = nullptr;
  InitChem ICS_pH = InitChem();
  BasVec *Vchem_key_tmp = nullptr;
  splayTree **st_lst = nullptr;
  Output output;
  double last_eff_write = 0., last_perm_write = 0.; // timing counter for effluent and perm file
  time_t t_begin = time(nullptr);
  int i;

  // struct BioParam BiP; // /* AMIN */
  /* EJE AMIN*/
  FILE *rate_data;
  rate_data = fopen("rate_data.dat", "w");
  //FILE *RhoOutput;
  //FILE *ChemOutput;
  //ChemOutput = fopen("FakeChemOutput.dat", "w");
  //RhoOutput = fopen("RhoOutput.dat", "w");
  // ----------------------------
  // general initializations, read input, setup mpi-run
  // ----------------------------
  sys = new System();
  read_cmd_args(argc, argv, sys);
  options.read(sys->files["options"]);  //options.print();
  input.read(sys->files["input"]);  //input.print();
  init_structs(options, &opt, &step, &minerals, &links, &node, sys, &fluid, &species, &dif_field, &bndry);

  read_input_v2(input, &InitChemSolver, sys, fluid, species, dif_field, bndry, argv);

  init_mpi_run(sys, &argc, &argv, species);

  init_variables_and_geometry(sys, bndry, fluid, species, dif_field);
  init_pH_cal(&ICS_pH, species, input, sys);
//  printf("before ph\n");
//  double ctot[] = {1e-3, 1e-8, 1e-8};
//  for (int i = 0; i < 10; ++i)
//    std::cout << "pH = " << get_pH(ICS_pH, ctot) << std::endl;
//   printf("after ph\n");

  output.init(sys);
  sys->output = &output;
  output.write_geo_file("geo", sys->node, sys);

  initial_temporal_rate(&InitChemSolver, sys,  input); // Changing rate parameters for variable dissolution rate  /* AMIN */

  init_3D_run(sys, species, fluid, dif_field, &InitChemSolver, &st_lst, minerals, bndry, &Vchem_key_tmp); // also works for 2D run


    output.dif.write(sys);  //  EJ AMIN Print the initial dist of the diffusive fields
    
  if (opt.perm_interval >= 0.0) {
    init_permeability_measurement(sys, fluid, species, minerals);
    if (opt.perm_interval == 0.0) {
      MPI_Finalize();
      exit(0);
    }
  }

  //---------------------------------
  // time loop
  //---------------------------------
 

  sys->nwrite++;
  for (sys->t = sys->start_itr; sys->t < sys->end_itr; ++(sys->t), sys->time += sys->dt) {

     // Chaging rate constants /* AMIN */

   //
    temporal_rate(&InitChemSolver[1], &Vchem_key_tmp,  sys,  input); // applying variable dissolution rate /* AMIN */
    //if (sys->t % 1000 == 0)
         fprintf(rate_data, "%lf \t %lf \t %lf \n", sys->time, Vchem_key_tmp[1].rate[0][0]*1e6, Vchem_key_tmp[1].rate[0][1] );
    //std::cout << "Vchem_key_tmp[1].rate[0][0] = " << Vchem_key_tmp[1].rate[0][0] << std::endl;
    // Streaming and collision
    advance_one_step(&InitChemSolver[1], species, fluid, dif_field, sys, minerals, bndry, st_lst);

    //printChem(species, sys, ChemOutput);
    //printRho(species, sys, RhoOutput);


	// Apply boundary condtions
    apply_boundary_conditions(&InitChemSolver[1], species, dif_field, fluid, sys, minerals, bndry, st_lst);

    // re-initialize velocity if nodes change and
    // update species concentration in new fluid nodes
#ifndef REINIT_OFF // true if velocity reinit is on
    if (sys->velreinit) {
      if (sys->fluid)
        reinit_velocity(sys->node, sys, fluid, 0);
      set_species_in_new_fluid_nodes(species, fluid, sys);
      wall_bc_bounce_back(sys->node, species, sys);
      communicate(species->f, species->n_c, sys);
      sys->velreinit = 0;
    }
#else
    if (sys->rho_reinit) {
      set_species_in_new_fluid_nodes(species, fluid, sys);
      //set_conc_field(species, fluid, sys); // only if steady-state run first
      // apply boundary conditions
      wall_bc_bounce_back(sys->node, species, sys);
      //bc_run_g(&InitChemSolver[1], species, fluid, sys, minerals, bndry, st_lst);
      // mpi
      communicate(species->f, species->n_c, sys);
      sys->rho_reinit = 0;
    }
#endif

	// After this we should not change anyting EJE /AMIN
    // free SplayTree if total size exceeds the limit
    if (sys->chem) {
      for (i = 0; i < minerals->no_comb; ++i) {
        if (st_lst[i]->size > (long int)opt.splaytree_length) {
          EmptySplayTree(st_lst[i]);
        }
      }
    }

 
    //---------------------------------
    // save data
    //---------------------------------

    // make independent permeability measurement
    if ( (opt.perm_interval>0) && (sys->time+sys->dt>last_perm_write+opt.perm_interval) ) {
      last_perm_write = sys->time + sys->dt;
      measure_permeability(sys, fluid, species, minerals);
    }


    // write and print interval
    print_status_and_write_output(ICS_pH, output, sys, t_begin, minerals->no_comb, &InitChemSolver[1], species, fluid, minerals, st_lst, opt.write_velocity);
	// write effluent concentration to file
   /* if (sys->t == sys->start_itr || sys->time + sys->dt > last_eff_write + opt.effluent_interval) {
      last_eff_write = sys->time + sys->dt;
      if (sys->flux > 0 || sys->opt->skip_init_velocity_run)
        write_effluent_file(sys, species, fluid, minerals, InitChemSolver, t_begin);
	} */

    // save binary files every restart_interval
    if (opt.restart_interval > 0 && (sys->t + 1) % opt.restart_interval == 0) {
      if (sys->fluid)
        save_fluid_restart_files(fluid, sys);
      save_chem_restart_files(species, minerals, sys);
    }

  }

  /* EJE AMIN*/
  //fclose(ChemOutput);
  //fclose(RhoOutput);
  fclose(rate_data);

  MPI_Finalize();
  return 0;
}


