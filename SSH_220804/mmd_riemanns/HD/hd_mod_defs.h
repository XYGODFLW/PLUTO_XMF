/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Set labels, indexes and prototypes for the HD module.

  Contains variable names and prototypes for the HD module

  \author A. Mignone (mignone@ph.unito.it)
  \date   April, 2, 2015
*/
/* ///////////////////////////////////////////////////////////////////// */

/* *********************************************************
    Set flow variable indices.
    Extra vector components, when not needed, point to the
    last element (255) of the array stored by startup.c.  
   ********************************************************* */
    //shared with MHD
    
/* ***********************************************************
                   Prototyping goes here          
   *********************************************************** */

int  hd_ConsToPrim   (double **, double **, int, int, unsigned char *);
void hd_PrimToCons   (double **, double **, int, int);

void hd_Eigenvalues (double **, double *, double **, int, int);

void hd_PrimEigenvectors (double *, double, double, double *, double **, double **);
void hd_ConsEigenvectors (double *, double *, double, 
                       double **, double **, double *);

void hd_Flux      (double **, double **, double *, double **, double *, int, int);
void hd_HLL_Speed (double **, double **, double *, double *, 
                double *, double *, int, int);
void hd_MaxSignalSpeed (double **, double *, double *, double *, int, int);

void hd_PrimRHS    (double *, double *, double, double, double *);
void hd_PrimSource (const State_1D *, int, int, 
                 double *, double *, double **, Grid *);

void hd_GetAreaFlux (const State_1D *, double **, double **, int, int, Grid *);

void hd_RightHandSide (const State_1D *, Time_Step *, int, int, double, Grid *);
void hd_RightHandSideSource (const State_1D *, Time_Step *, int, int, double,
                          double *, Grid *);

Riemann_Solver *hd_SetSolver (const char *solver);

Riemann_Solver hd_LF_Solver, hd_Roe_Solver, hd_HLL_Solver, hd_HLLC_Solver;




















