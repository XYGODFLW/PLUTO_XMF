#include"pluto.h"

/* ******************************************************************** */
int mhd_Integrate (Data *d, Riemann_Solver *Solver, Time_Step *Dts, Grid *grid)//RK_var
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
  int i,j,k;
  
  g_maxMach = 0.0;
  g_maxRiemannIter = 0;
  g_maxRootIter    = 0;

  g_operatorStep = HYPERBOLIC_STEP;
  #if DIMENSIONAL_SPLITTING == YES
    for (g_dir = 0; g_dir < DIMENSIONS; g_dir++){
      if (mhd_AdvanceStep (d, Solver, Dts, grid) != 0) return (1);
    }
  #else
  if (mhd_AdvanceStep (d, Solver, Dts, grid) != 0) return(1);
  #endif

  return (0); 
}








