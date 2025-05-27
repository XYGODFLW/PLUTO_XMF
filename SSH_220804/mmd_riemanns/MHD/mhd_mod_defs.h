/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Set labels, indexes and prototypes for the MHD module.

  Contains basic macro definitions, structure definitions and global
  variable declarations used by the MHD module.

  \author A. Mignone (mignone@ph.unito.it)
  \date   April, 2, 2015
*/
/* ///////////////////////////////////////////////////////////////////// */

/* *********************************************************
    Set flow variable indices.
    Extra vector components, when not needed, point to the
    last element (255) of the array stored by startup.c.  
   ********************************************************* */

#define  RHO 0

#define  MX1 1
#define  MX2 (COMPONENTS >= 2 ? 2: 255)
#define  MX3 (COMPONENTS == 3 ? 3: 255) 

#if HAVE_ENERGY
  #define ENG  (COMPONENTS + 1)
  #define PRS  ENG
#endif

#define  BX1 (HAVE_ENERGY + COMPONENTS + 1)
#define  BX2 (COMPONENTS >= 2 ? (BX1+1): 255)
#define  BX3 (COMPONENTS == 3 ? (BX1+2): 255)

#if DIVB_CONTROL == DIV_CLEANING
  #define PSI_GLM  (2*COMPONENTS + 1 + HAVE_ENERGY)
#endif

#define VX1   MX1
#define VX2   MX2
#define VX3   MX3

#define NFLX (1 + 2*COMPONENTS + HAVE_ENERGY + (DIVB_CONTROL == DIV_CLEANING))
#define NFLX_LOOP(n)     for ((n) = NFLX;   (n)--;  )
#define NVAR NFLX
#define NVAR_LOOP(n)     for ((n) = NVAR;   (n)--;  )
#define VAR_LOOP(n)   for ((n) = NVAR; (n)--;    )

#define MHD_NFLX (1 + 2*COMPONENTS + HAVE_ENERGY + (DIVB_CONTROL == DIV_CLEANING))
#define MHD_NFLX_LOOP(n) for ((n) = MHD_NFLX;   (n)--;  )
#define MHD_NVAR MHD_NFLX 
#define MHD_NVAR_LOOP(n)     for ((n) = MHD_NVAR;   (n)--;  )
#define MHD_VAR_LOOP(n)   for ((n) = MHD_NVAR; (n)--;    )

#define HD_NFLX (1 + COMPONENTS + HAVE_ENERGY)
#define HD_NFLX_LOOP(n) for ((n) = HD_NFLX;   (n)--;  )
#define HD_NVAR HD_NFLX
#define HD_NVAR_LOOP(n)     for ((n) = HD_NVAR;   (n)--;  )
#define HD_VAR_LOOP(n)   for ((n) = HD_NVAR; (n)--;    )



/* ****************************************************************
   **************************************************************** */

/* ********************************************************************* */
/*! Label the different waves in increasing order 
    following the number of vector components.

    \b IMPORTANT: the KPSI_GLMM & KPSI_GLMP modes are 
                  present only in the MHD-GLM formulation.
                  We keep them at the END of the enumeration
                  so we can skip them in unnecessary loops.
                  Please do NOT change them !
   ********************************************************************* */

enum KWAVES {
 KFASTM, KFASTP
 #if HAVE_ENERGY
  , KENTRP
 #endif

 #if DIVB_CONTROL != DIV_CLEANING
  , KDIVB
 #endif

 #if COMPONENTS >= 2
  , KSLOWM, KSLOWP
  #if COMPONENTS == 3
   , KALFVM, KALFVP
  #endif
 #endif

 #if DIVB_CONTROL == DIV_CLEANING  
  , KPSI_GLMM, KPSI_GLMP 
 #endif
};

/*! \name Vector Potential Labels 
    These may only be used in the STARTUP / INIT  functions.
    They're convenient in obtaining a discretization that preserve 
    the divergence-free condition (for staggered field) or if you simply
    wish to initialize the magnetic field from the vector potential.     */
/**@{ */
#define   AX1  (NVAR + 1)
#define   AX2  (NVAR + 2)
#define   AX3  (NVAR + 3)
/**@} */

#define AX  AX1  
#define AY  AX2
#define AZ  AX3

/* *************************************************
     Now define more convenient and user-friendly 
     pointer labels for geometry setting      
   ************************************************* */

#if GEOMETRY == CARTESIAN
  #define VX    VX1
  #define VY    VX2
  #define VZ    VX3

  #define MX    MX1
  #define MY    MX2
  #define MZ    MX3

  #define BX    BX1
  #define BY    BX2
  #define BZ    BX3
#endif

#if GEOMETRY == CYLINDRICAL 
 #define iVR    VX1
 #define iVZ    VX2
 #define iVPHI  VX3

 #define iMR    MX1
 #define iMZ    MX2
 #define iMPHI  MX3

 #define iBR    BX1
 #define iBZ    BX2
 #define iBPHI  BX3
#endif

#if GEOMETRY == POLAR 
 #define iVR    VX1
 #define iVPHI  VX2
 #define iVZ    VX3

 #define iMR    MX1
 #define iMPHI  MX2
 #define iMZ    MX3

 #define iBR    BX1
 #define iBPHI  BX2
 #define iBZ    BX3
#endif

#if GEOMETRY == SPHERICAL 
 #define iVR     VX1
 #define iVTH    VX2
 #define iVPHI   VX3

 #define iMR    MX1
 #define iMTH   MX2
 #define iMPHI  MX3

 #define iBR    BX1
 #define iBTH   BX2
 #define iBPHI  BX3
#endif

/* ***********************************************************
    \cond REPEAT_FUNCTION_DOCUMENTATION_IN_HEADER_FILES 
    Function prototyping
   *********************************************************** */

void mhd_BackgroundField (double x1, double x2, double x3, double *B0);

#if BACKGROUND_FIELD == YES
 double **mhd_GetBackgroundField (int, int, int, Grid *);
#endif 

int  mhd_ConsToPrim   (double **, real **, int , int, unsigned char *);
void mhd_PrimToCons   (double **, double **, int, int);

void mhd_Eigenvalues (double **, double *, double **, int, int);

void mhd_PrimEigenvectors (double *, double, double, double *, double **, double **);
void mhd_ConsEigenvectors (double *, double *, double, double **, double **, double *);

void mhd_Flux      (double **, double **, double *, double **, double **,
                double *, int, int);
void mhd_HLL_Speed (double **, double **, double *, double *, double **, 
                double *, double *, int, int);
void mhd_MaxSignalSpeed (double **, double *, double *, double *, double **, int, int);

void mhd_PrimRHS    (double *, double *, double, double, double *);
void mhd_PrimSource (const State_1D *, int, int, 
                 double *, double *, double **, Grid *);

void mhd_GetAreaFlux (const State_1D *, double **, double **, int, int, Grid *);
void mhd_PrimToChar (double **, double *, double *); 

void mhd_RightHandSide (const State_1D *, Time_Step *, int, int, double, Grid *);
void mhd_RightHandSideSource (const State_1D *, Time_Step *, int, int, double,
                          double *, Grid *);

Riemann_Solver *mhd_SetSolver (const char *);


#if DIVB_CONTROL == EIGHT_WAVES

 void mhd_Roe_DivBSource (const State_1D *, int, int, Grid *);
 void mhd_HLL_DivBSource (const State_1D *, double **, int, int, Grid *);

#elif DIVB_CONTROL == DIV_CLEANING

 #include "mmd_riemanns/MHD/GLM/glm.h"

#elif DIVB_CONTROL == CONSTRAINED_TRANSPORT

 #include "mmd_riemanns/MHD/CT/ct.h"

#endif

Riemann_Solver mhd_HLL_Solver, mhd_HLLC_Solver, mhd_HLLD_Solver;
Riemann_Solver mhd_LF_Solver, mhd_Roe_Solver;

#if RESISTIVITY != NO
 #include "mmd_riemanns/Resistivity/res.h"
#endif

#ifdef SHEARINGBOX
 #include "mmd_riemanns/MHD/ShearingBox/shearingbox.h"
#endif
/* \endcond */




