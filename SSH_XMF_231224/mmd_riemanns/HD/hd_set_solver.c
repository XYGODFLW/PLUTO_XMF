#include"pluto.h"

/* ********************************************************************* */
Riemann_Solver *hd_SetSolver (const char *solver)
/*!
 *
 * PURPOSE
 *
 *   return a pointer to a riemann solver function
 *
 *********************************************************************** */
{
  
/* ------------------------------------------------------
       Set Pointers for SOLVERS 
   ------------------------------------------------------ */

  #if EOS == IDEAL 
   if (!strcmp(solver, "tvdlf"))  return (&hd_LF_Solver);
   else if (!strcmp(solver, "roe"))    return (&hd_Roe_Solver);
   else if (!strcmp(solver, "hlle") || 
            !strcmp(solver, "hll"))    return (&hd_HLL_Solver);
   else if (!strcmp(solver, "hllc"))   return (&hd_HLLC_Solver);
/*
   else if (!strcmp(solver, "rusanov_dw")) return (&RusanovDW_Solver);
*/
  #elif EOS == PVTE_LAW
   if (!strcmp(solver, "tvdlf"))       return (&hd_LF_Solver);
   else if (!strcmp(solver, "hlle") || 
            !strcmp(solver, "hll"))    return (&hd_HLL_Solver);
   else if (!strcmp(solver, "hllc"))   return (&hd_HLLC_Solver);
/*
   else if (!strcmp(solver, "rusanov_dw")) return (&RusanovDW_Solver);
*/
  #elif EOS == ISOTHERMAL
   if (!strcmp(solver, "tvdlf"))       return (&hd_LF_Solver);
   else if (!strcmp(solver, "roe"))    return (&hd_Roe_Solver);
   else if (!strcmp(solver, "hlle") || 
            !strcmp(solver, "hll"))    return (&hd_HLL_Solver);
   else if (!strcmp(solver, "hllc"))   return (&hd_HLLC_Solver);
/*
   else if (!strcmp(solver, "rusanov_dw")) return (&RusanovDW_Solver);
*/
  #endif

  print1 ("\n! SetSolver: '%s' not available with this configuration.\n", 
          solver);
  QUIT_PLUTO(1);
}
