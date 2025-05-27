#include "pluto.h"

/* ****************************************************************** */
void mhd_SetDefaultVarNames(Output *output)
/*
 *
 *  PURPOSE
 *
 *    Set file names for I/O
 *
 *
 ******************************************************************** */
{
  int nv;

/* ----------------------------------------------
    Physics module file names; 
    these pertain to the physics module ONLY
   ---------------------------------------------- */

  output->var_name[RHO] = "rho";
  EXPAND(output->var_name[VX1] = "vx1";  ,
         output->var_name[VX2] = "vx2";  ,
         output->var_name[VX3] = "vx3";)
#if HAVE_ENERGY
  output->var_name[PRS] = "prs";
#endif
  
#if PHYSICS == MHD || PHYSICS == RMHD
  EXPAND(output->var_name[BX1] = "bx1";  ,
         output->var_name[BX2] = "bx2";  ,
         output->var_name[BX3] = "bx3";)
#endif
  
  /* (staggered field names are set in SetOutput) */

#ifdef GLM_MHD
  output->var_name[PSI_GLM] = "psi_glm";
#endif
 
/* ------------------------------------------------
    Dust
   ------------------------------------------------ */

#if DUST == YES
  output->var_name[RHO_D] = "rho_d";
  EXPAND(output->var_name[VX1_D] = "vx1_d";  ,
         output->var_name[VX2_D] = "vx2_d";  ,
         output->var_name[VX3_D] = "vx3_d";)
#endif

/* ------------------------------------------------
                   Tracers 
   ------------------------------------------------ */

  NTRACER_LOOP(nv) sprintf (output->var_name[nv],"tr%d",nv - TRC + 1);

  #if ENTROPY_SWITCH
   sprintf (output->var_name[ENTR],"entropy");
  #endif

  return;
}




