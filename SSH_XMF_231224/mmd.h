#ifndef MMD_H
#define MMD_H

//MMD_setting beg
#define NFLD   2 //number of fluids 
#define ELC_ON 0
#define MAG_ON 1

#define FLD_LOOP(f)   for ((f) = 0; (f) < NFLD; (f)++)
#define RK4_LOOP(k)   for ((k) = 1; (k) <= 4; (k)++)
#define RK3_LOOP(k)   for ((k) = 1; (k) <= 3; (k)++)

#define CPN_LOOP(nv)  for ((nv) = 0; (nv)<COMPONENTS ;(nv)++) //supplyment for origin PLUTO

enum{F000,F001,F002,F003,F004,
     F005,F006,F007,F008,F009,};
enum{NEU,ION,ELC};

//element
enum{ELECTRON, H, HE,LI,BE,B,C,N,O,F,NE,};

//fluid sequence
enum{HI, HII,};//HEI,HEII, ELCI

extern double g_dt_temp;

//In this mode base variable is number density, velocity and Tempture
#define NDT RHO
#define TMP PRS

#define BM1 0  //B in model
#define BM2 1
#define BM3 2 
//Variables looping tools for MHD
#if ELC_ON==1
#define RK_ELC (HD_NVAR*NFLD)
#endif

#if MAG_ON==1
#define RK_MAG (HD_NVAR*NFLD+ELC_ON)
#endif

#define MD_TOT (HD_NVAR*NFLD+ELC_ON+MAG_ON*COMPONENTS)//Total model variables loop.
#define RK_TOT MD_TOT
#define RK_VAR_LOOP(nv) for(nv=0;nv<RK_TOT;nv++)
#define IM_TOT (HD_NVAR*NFLD+ELC_ON)
#define IM_VAR_LOOP(nv) for(nv=0;nv<IM_TOT;nv++) //Loop for the variables which need implicit treatment.



//MMD_setting end
#include "mmd_structs.h"
#include "mmd_globals.h"

#include "mhd_mod_defs.h"  /* Include MHD physics header file (search path is set
                          in the makefile) */
#include "hd_mod_defs.h"   /* Include HD physics header file (search path is set
                          in the makefile) */

#include "mmd_prototypes.h"
#include "mhd_convection.h"
#include "hd_convection.h"

//MMD_state_machine
#include "mmd_sm.h"

//Common tools
#include "pd.h"
#include "matrix_tools.h"
#include "aly.h"

#endif

//EXT
#define RK_EXT(RK,nv,k,j,i)         ((RK)[nv][k][j][i]) 
#define FLD_EXT(model,f,nv,k,j,i)   ((model)->fluid[f].Vc[nv][k][j][i])
#define ELC_EXT(model,nv,k,j,i)     ((model)->E_fluid[nv][k][j][i])
#define MAG_EXT(model,nv,k,j,i)    ((model)->B_field[nv][k][j][i])
#define MGS_EXT(model,nv,k,j,i)     ((model)->Bs[nv][k][j][i])

#define FLUID_EXT(fluid,f,nv,k,j,i) ((fluid)[f].Vc[nv][k][j][i])
#define FLUX_EXT(model,b,k,j,i) ((model)->xuv.flux[b][k][j][i])

#define ATOM_MASS(model,f) ((model)->fluid[f].ptc.mass*CONST_amu)
#define ATOM_MASS_AM(model,f) ((model)->fluid[f].ptc.mass)

//loop_tools 
#define BEAM_LOOP(b,nbeam) for((b)=0;(b)<(nbeam);(b)++)

//other tools
#define G_DT g_dt*UNIT_TIME

#define RLT_NDT(model,f1,f2,k,j,i)  (MODEL_EXT(model,f1,NDT,k,j,i)/MODEL_EXT(model,f2,NDT,k,j,i)) // f1/f2 
#define RLT_RHO(model,f1,f2,k,j,i)  (RLT_NDT(model,f1,f2,k,j,i)*ATOM_MASS(model,f1)/ATOM_MASS(model,f2)) // f1/f2
#define DIFF_VX1(model,f1,f2,k,j,i) (MODEL_EXT(model,f1,VX1,k,j,i) - MODEL_EXT(model,f2,VX1,k,j,i) )// f1-f2
#define DIFF_TMP(model,f1,f2,k,j,i) (MODEL_EXT(model,f1,TMP,k,j,i) - MODEL_EXT(model,f2,TMP,k,j,i) ) // f1-f2

#define TMP_REDUCE(model,f1,f2,k,j,i) ( (MODEL_EXT(model,f1,PRS,k,j,i)/MODEL_EXT(model,f1,NDT,k,j,i) * ATOM_MASS_AM(model,f2) \
                                      + MODEL_EXT(model,f2,PRS,k,j,i)/MODEL_EXT(model,f2,NDT,k,j,i) * ATOM_MASS_AM(model,f1)) \
                                      /(ATOM_MASS_AM(model,f1)+ATOM_MASS_AM(model,f2))/CONST_kB )

#define MASS_REDUCE(model,f1,f2) ( ATOM_MASS(model,f1)*ATOM_MASS(model,f2)/(ATOM_MASS(model,f1)+ATOM_MASS(model,f2)) )

#define RHO_EXT(model,f,k,j,i) (MODEL_EXT(model,f,NDT,k,j,i)*ATOM_MASS(model,f))                                //RHO in c.g.s
#define MX1_EXT(model,f,k,j,i) (MODEL_EXT(model,f,NDT,k,j,i)*ATOM_MASS(model,f) * MODEL_EXT(model,f,VX1,k,j,i)) //fixed fraction of momenton
#define PRS_EXT(model,f,k,j,i) (MODEL_EXT(model,f,NDT,k,j,i) * MODEL_EXT(model,f,TMP,k,j,i) )                    //fixed fraction of pressure


#if DIMENSIONS == 1 //reversed dimension expand
 #define RD_EXPAND(a,b,c) c
#endif

#if DIMENSIONS == 2
 #define RD_EXPAND(a,b,c) b c
#endif

#if DIMENSIONS == 3
 #define RD_EXPAND(a,b,c) a b c
#endif























