#ifndef CHEM_DEFINES_H
#define CHEM_DEFINES_H

#define COMPLEXNAME 20
#define LNTEN 2.302585092994
#define CHEM_LOWER_ 1.0e-20    /* for lower limit on chemical concentrations */
#define CHEM_EPSILON_ 1.0e-20  /* for numerical calculations */
#define CHEM_NEW_CALC_ 1e-5     /* call chemical solver */
#define CHEM_MAXITER_MASS_ 1000  /* max iter in the massbalance routine */ 
#define CHEM_MAX_PH_ITER_  500  /* max iter in the charge balance solver */
#define PH_MAX_STEP_ .5        /* max length of pH step in each iteration */
#define MBAL_MAX_STEP_ 5.        /* max length of log10 step in each iteration for massbalance*/
#define MBAL_SURF_MAX_STEP_ 1.2        /* max length of log10 step in each iteration for surface potential*/
#define JULIUS 0

#endif // CHEM_DEFINES_H
