#include "pluto.h"

/* ****************************************************************** */
void hd_SetDefaultVarNames(Output *output)
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
  
/* ------------------------------------------------
                   Tracers 
   ------------------------------------------------ */

  NTRACER_LOOP(nv) sprintf (output->var_name[nv],"tr%d",nv - TRC + 1);

  #if ENTROPY_SWITCH
   sprintf (output->var_name[ENTR],"entropy");
  #endif

  return;
}




