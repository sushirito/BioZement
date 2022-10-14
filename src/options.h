/*
 * options.h
 *
 *  Created on: 23. jun. 2016
 *      Author: janlv
 */

#ifndef SRC_OPTIONS_H_
#define SRC_OPTIONS_H_

//#include "Input.h"
#include <string>
class Input;

struct Options {
  // ----------------------------
  //  general switches
  // ----------------------------
  int charge_balance;         // 0: advect H+ (should be faster), number of tau's in inp.dat must be nspecies+1
  int find_steady_state;      // 1: find steady-state solution of species before chem-solver is turned on
  int kick_start;             // 1: inlet conc in all fluid nodes 2: inlet conc only in high vel nodes
  int skip_init_velocity_run; // 0: init velocity to steady-state before main loop, 1: velocity is zero
  int write_velocity;         // switch on or off velocity data to file (vel is unchanged if nodes are unchanged)
  int read_vel_from_file;     // read velocity and species from file
  int read_species_from_file; // (remember to set correct gx in inp.dat)
  int find_isolated_nodes;    //
  int set_isolated_inert;     //
  int restart_run;            //
  int boost_rate;
  int measure_inlet_perm;
  
  // ----------------------------
  //  geometry control parameters
  // ----------------------------
  int straight_tube;      // 1: porosity is set to 1, 0: porosity is calculated
  int flip_x;             // flood in opposite x-dir, i.e. rotate geometry 90 deg with axis in y-dir
  int rotate_y;           // rotate geometry 45 deg clockwise with axis in y-dir
  int rotate_z;           // rotate geometry 45 deg clockwise with axis in z-dir
  int add_in_out_wall;
  int add_raster;
  int raster_nodes;
  int wall_nodes;
  int copy_in_out;
  int add_union_in_out;
  int closed_cell;        // 1: close inlet and outlet to do pure diffusive simulation
  int in_nodes;
  int out_nodes;
  int inertx_outlet;      // this number of nodes at inlet
  int inertx_inlet;       // and outlet are excluded from chem reactions
  int inert_walls;        // 1: outer walls are inert
  double magn_percent;
  int magn_surface;
  double magn_limit;  // at which solid fraction non-nucleation sites can precipitate
  
  // ----------------------------
  //  output control parameters
  // ----------------------------
  int output_pH;
  int eff_offset;             // position from outlet where flux and effluent is measured (should be inert nodes)
  double effluent_interval;   // how often (in sec.) effluent data are written to file (neg. values = every step)
  double perm_interval;
  int restart_interval;
  
  // ----------------------------
  //  LB parameters
  // ----------------------------
  double tau_min; // minimum tau-value
  
  // ----------------------------
  //  splaytree control parameters
  // ----------------------------
  int splaytree_ON;
  int splaytree_length;

};

void set_options(struct Options *opt, Input &o);

#endif /* SRC_OPTIONS_H_ */
