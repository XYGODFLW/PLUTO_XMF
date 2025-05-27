/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Return a pointer to a Riemann solver function.

  \author A. Mignone (mignone@ph.unito.it)
  \date   June 5, 2013
*/
/* ///////////////////////////////////////////////////////////////////// */
#include"pluto.h"

/* ********************************************************************* */
Riemann_Solver *mhd_SetSolver (const char *solver)
/*!
 *  Depending on the choice of the Riemann solver specified in
 *  pluto.ini, return a pointer to the corresponding Riemann solver
 *  function
 *
 *********************************************************************** */
{
  
/* ------------------------------------------------------
       Set Pointers for SOLVERS
   ------------------------------------------------------ */

  #if EOS == IDEAL
  
   if      (!strcmp(solver, "tvdlf"))   return (&mhd_LF_Solver);
   else if (!strcmp(solver, "roe"))     return (&mhd_Roe_Solver);
   else if (!strcmp(solver, "hlle") ||
            !strcmp(solver, "hll"))     return (&mhd_HLL_Solver);
   else if (!strcmp(solver, "hllc"))    return (&mhd_HLLC_Solver);
   else if (!strcmp(solver, "hlld"))    return (&mhd_HLLD_Solver);
/*
   else if (!strcmp(solver, "musta"))   return (&MUSTA_Solver);
*/

   #elif EOS == PVTE_LAW
   if (!strcmp(solver, "tvdlf"))       return (&mhd_LF_Solver);
   else if (!strcmp(solver, "hlle") || 
            !strcmp(solver, "hll"))    return (&mhd_HLL_Solver);
   else if (!strcmp(solver, "hllc"))   return (&mhd_HLLC_Solver);
   else if (!strcmp(solver, "hlld"))   return (&mhd_HLLD_Solver);

  #elif EOS == ISOTHERMAL

   if      (!strcmp(solver, "tvdlf"))   return (&mhd_LF_Solver);
   else if (!strcmp(solver, "roe"))     return (&mhd_Roe_Solver);
   else if (!strcmp(solver, "hlle") ||
            !strcmp(solver, "hll"))     return (&mhd_HLL_Solver);
   else if (!strcmp(solver, "hlld"))    return (&mhd_HLLD_Solver);
/*
   else if (!strcmp(solver, "musta"))   return (&MUSTA_Solver);
*/

  #elif EOS == BAROTROPIC

   if (!strcmp(solver, "tvdlf"))        return (&mhd_LF_Solver);
   else if (!strcmp(solver, "hlle") ||
             !strcmp(solver, "hll"))    return (&mhd_HLL_Solver);

  #endif

  print1 ("\n! SetSolver: '%s' is not available.\n", solver);
  QUIT_PLUTO(1);

}

