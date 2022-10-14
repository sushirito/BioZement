/*
 * macros.h
 * 
 *  Created on: 1. jul. 2016
 *      Author: janlv
 */

#ifndef SRC_MACROS_H_
#define SRC_MACROS_H_

// Collection of pre-processor macros that control the code
// To turn a macro off: comment the #define statement (or add a #undef statement below)

//#define _RANDOM_MINERAL_
//#define _MAGN_ON_MAGN_

//#ifndef _USE_COLOR_GRAD_
//#define _USE_COLOR_GRAD_
//#endif

//#define _FIND_PERC_NODES_

//
// ***
//#define _PERIODIC_RUN_
// ***


// Uncomment this line to do a pure diffusive simulation
//#define DIFFUSIVE_RUN
//#define NOT_MOVE
//#define _SPECIES_SOURCE_BULK_
//#define _FLUID_BDRY_PRESS_

//#define _FLUID_SOURCE_
//#define _SPECIES_SOURCE_
// -----------------------------------

// If on: the new definition of concentration which accounts for compressibility is used.
// For pure diffuive simulations it should be undefined (its done automatically)
#define NEWCONC  // use psi = sum g / sum f as concentration
// -----------------------------------


// If on: force tau to be one.
// If off: tau is allowed to vary to increase the timesteps.
//         This can lead to small neg. rho-values, and the macro
//         SET_NEG_RHO_ZERO is auto-set if this macro is not set (i.e. if tau != 1)
//#define _TAU_FROM_INP_  // affects calculation of dt and hence comp. efficiency
#define SET_NEG_RHO_ZERO
// -----------------------------------

// If on: find steady-state solution and jump in time until nodes are empty or full
// #define CONVERGE
// #define _MEAN_CONVERGENCE_
// -----------------------------------

// Controls the outlet zero-flux boundary condition
// If on: dist-func at outlet copied to the 2. and 3.nd last node
// If off: dist-func at outlet only copied to the 2. last node
//#define OUTLET2
// -----------------------------------

// If on: the jump in time caused by a steady-state solution is split into fixed intervals
//        given by the reporting timestep in inp.dat
// WARNING: should be off if reporting timestep is small, i.e. < 100
#define KEEP_INTERVAL
// -----------------------------------

// If on: use linear rate law J = -k*(c-c_eq) and skip chem-solver
//#define FAKE_CHEM
// -----------------------------------

// If on: skip interaction with the mineral surface by calling bounce-back for all links
//#define ALL_INERT
// -----------------------------------

// If on: fluid is off, no advection in the system
//#define FLUID_OFF
// -----------------------------------

// If on: no reinitialization of velocity if nodes change to/from fluid/solid
//#define REINIT_OFF
// -----------------------------------

// If on: use bounce-back inlet boundary condition
//#define BB_INLET
// -----------------------------------

// Debugging macros that control output to screen
//#define _DEBUG_MPI_ 1
//#define _DEBUG_VEL_ 1
//#define _DEBUG_ 1
//#define _DEBUG_CONVERGENCE_ 1
#define _DEBUG_PERM_ 1

// some macros control other macros

#ifdef _FIND_PERC_NODES_
#define DIFFUSIVE_RUN
#endif

#ifdef _SPECIES_SOURCE_BULK_
#define DIFFUSIVE_RUN
#endif

#ifdef DIFFUSIVE_RUN
#define FLUID_OFF
#define REINIT_OFF
#endif

#ifdef NEWCONC
#undef REINIT_OFF
#endif

// ALWAYS ON
#define _MPI_ 1     // these macros are always on, but we still keep them here
#define FLUIDCHEM 1 // since the code is not completely cleaned yet

#endif /* SRC_MACROS_H_ */
