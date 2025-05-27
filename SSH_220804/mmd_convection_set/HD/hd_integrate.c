#include"pluto.h"

/* ******************************************************************** */
int hd_Integrate (Data* d, Riemann_Solver *Solver, Time_Step *Dts, Grid *grid)
/*!
 * Advance equations by a single time-step.

 * \param  d      pointer to PLUTO Data structure;
 * \param  Solver pointer to a Riemann solver function;
 * \param  Dts    pointer to time Step structure;
 * \param  grid   pointer to grid structure.
 * 
 * \return An integer giving success / failure (development).
 * 
 ********************************************************************** */
{
  int idim, err = 0;
  int k,j,i,nv;
   
  g_maxMach = 0.0;
  g_maxRiemannIter = 0;
  g_maxRootIter    = 0;

/* -------------------------------------------------------
    Initialize max propagation speed in Dedner's approach
   ------------------------------------------------------- */


  /* ---------------------------------------------
        perform Strang Splitting on directions 
        (if necessary) and sources 
     --------------------------------------------- */

  #ifdef FARGO
   FARGO_ComputeVelocity(d, grid);
  #endif

  g_operatorStep = HYPERBOLIC_STEP;
  #if DIMENSIONAL_SPLITTING == YES
   for (g_dir = 0; g_dir < DIMENSIONS; g_dir++){
     if (hd_AdvanceStep (d, Solver, Dts, grid) != 0) return (1);
   }
  #else
  if (hd_AdvanceStep (d, Solver, Dts, grid) != 0) return(1);
  #endif

  //du = new_u 
  sm_dt = hd_NextTimeStep (Dts, grid);
  
  HD_NVAR_LOOP(nv)DOM_LOOP(k,j,i){
    sm_du[nv][k][j][i] = d->Vc[nv][k][j][i];
  }

  return (0); /* -- ok, step achieved -- */
}










