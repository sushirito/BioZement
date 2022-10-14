/*
 * options.c
 *
 *  Created on: 23. jun. 2016
 *      Author: janlv
 */

#include "options.h"
#include "macros.h"
#include <cstdio>
#include "Input.h"

void set_options(struct Options *opt, Input &o)
{

  // ----------------------------
  //  general switches
  // ----------------------------
  opt->charge_balance = 1; //o["switch"]["charge_balance"];
  opt->find_steady_state = o["switch"]["find_steady_state"];
  opt->kick_start = o["switch"]["kick_start"];
  opt->skip_init_velocity_run = o["switch"]["skip_init_velocity_run"];
  opt->write_velocity = o["switch"]["write_velocity"];
  opt->read_vel_from_file = o["switch"]["read_vel_from_file"];
  opt->read_species_from_file = o["switch"]["read_species_from_file"];
  opt->find_isolated_nodes = o["switch"]["find_isolated_nodes"];
  opt->set_isolated_inert = 0; //o["switch"]["set_isolated_inert"];
  opt->restart_run = 0; //o["switch"]["restart_run"];
  opt->boost_rate = o["switch"]["boost_rate"];
  opt->measure_inlet_perm = 0; //o["switch"]["measure_inlet_perm"];

  // ----------------------------
  //  geometry control parameters
  // ----------------------------
  opt->straight_tube = 0; //o["geometry"]["straight_tube"];
  opt->flip_x = o["geometry"]["flip_x"];
  opt->rotate_y = o["geometry"]["rotate_y"];
  opt->rotate_z = o["geometry"]["rotate_z"];
  opt->add_raster = 0; //o["geometry"]["add_raster"];
  opt->raster_nodes = 0; //o["geometry"]["raster_nodes"];
  opt->wall_nodes = o["geometry"]["wall_nodes"];
  opt->copy_in_out = o["geometry"]["copy_in_out"];
  opt->add_union_in_out = o["geometry"]["union_in_out"];
  opt->closed_cell = o["geometry"]["closed_cell"];
  opt->in_nodes = o["geometry"]["in_nodes"];
  opt->out_nodes = o["geometry"]["out_nodes"];
  opt->inertx_inlet = o["geometry"]["inertx_inlet"];
  opt->inertx_outlet = o["geometry"]["inertx_outlet"];
  opt->inert_walls = o["geometry"]["inert_walls"];
  #ifdef _RANDOM_MINERAL_
  opt->magn_surface = o["geometry"]["magn_surface"];
  opt->magn_percent = o["geometry"]["magn_percent"];
  opt->magn_limit = o["geometry"]["magn_limit"];
#else
  opt->magn_surface = 0; 
  opt->magn_percent = 0; 
  opt->magn_limit = 0; 
#endif
  //opt->add_in_out_wall = o["geometry"]["add_in_out_wall"];
  if (opt->in_nodes > 0 || opt->out_nodes > 0 || opt->wall_nodes > 0)
    opt->add_in_out_wall = 1;

  // ----------------------------
  //  output control parameters
  // ----------------------------
  opt->output_pH = o["output"]["pH"];
  opt->eff_offset = o["output"]["eff_offset"];
  opt->effluent_interval = o["output"]["effluent_interval"];
  opt->perm_interval = o["output"]["perm_interval"];
  opt->restart_interval = 0; //o["output"]["restart_interval"];

  // ----------------------------
  //  LB parameters
  // ----------------------------
  opt->tau_min = 0.53;// o["LB"]["tau_min"];

  // ----------------------------
  //  splaytree control parameters
  // ----------------------------
  opt->splaytree_ON = o["splaytree"]["on"];
  opt->splaytree_length = o["splaytree"]["length"];


  // ----------------------------
  // Dependent options
  // ----------------------------
  if (opt->add_union_in_out) {
    opt->copy_in_out = 1;
  }
}
