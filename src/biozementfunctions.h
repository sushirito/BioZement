#ifndef BIOZEMENTFUNCTIONS_H
#define BIOZEMENTFUNCTIONS_H

#include "global.h"
#include "output.h"
#include "chem_global.h"

void applyFakeChem(struct Links::Link *link, struct Field *species, struct Step *step, struct System *sys);
void printChem(struct Field *species, struct System *sys, FILE *FakeChemOutput);
void read_Rfile(struct Node *node, struct System *sys);
void read_Bfile(struct Node *node, struct System *sys);
void printRho(struct Field *species, struct System *sys, FILE *Rho);
real get_pH(InitChem &ICS, real * c_vchem_eq);
void init_pH_cal(InitChem *ICS, struct Field *species, Input &input, struct System *sys);
void set_nodes_inert(struct Node *node);
void temporal_rate(struct InitChem *ICS_ptr, struct BasVec **Vchem_key_tmp, System *sys, Input &inp);
void initial_temporal_rate(struct InitChem **ICS_ptr, struct System *sys, Input &inp);
void update_pH_buffer(InitChem &ICS_pH, Outfile* chem, struct Field *species, struct System *sys);
void init_run_3D_dif(struct Field *species, struct Field *fluid, struct System *sys);

#endif // BIOZEMENTFUNCTIONS_H
