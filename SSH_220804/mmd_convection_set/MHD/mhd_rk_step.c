/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Advance equations with Runge Kutta time integrators.

  Main driver for RK split/unsplit integrations and finite difference
  methods (RK3).
  Time stepping include Euler, RK2 and RK3.

  \authors A. Mignone (mignone@ph.unito.it)\n
           P. Tzeferacos (petros.tzeferacos@ph.unito.it)
  \date    Dec 18, 2014
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ********************************************************************* */
int mhd_AdvanceStep (const Data *d, Riemann_Solver *Riemann, 
                 Time_Step *Dts, Grid *grid) //output is needed which also could be set as a global variable
/*!
 * Advance the equations by a single time step using unsplit 
 * integrators based on the method of lines.
 *
 * \param [in,out]      d  pointer to Data structure
 * \param [in]    Riemann  pointer to a Riemann solver function
 * \param [in,out]    Dts  pointer to time step structure
 * \param [in]       grid  pointer to array of Grid structures
 *    
 *********************************************************************** */
{
  int  i, j, k, nv;
  static first_call = 1;
   
  //updating

  /* -------------------------------------------------------
    Initialize max propagation speed in Dedner's approach
   ------------------------------------------------------- */

  #ifdef GLM_MHD  /* -- initialize glm_ch -- */
   GLM_Init (d, Dts, grid);   
   GLM_Source (d->Vc, 0.5*g_dt, grid);
  #endif

  /* ---------------------------------------------
        perform Strang Splitting on directions 
        (if necessary) and sources 
     --------------------------------------------- */

  TOT_LOOP(k,j,i) d->flag[k][j][i] = 0;

  #ifdef FARGO
   FARGO_ComputeVelocity(d, grid);
  #endif

  mhd_UpdateStage(d, d->Uc, NULL, Riemann, g_dt, Dts, grid);
  #ifdef GLM_MHD  /* -- GLM source for dt/2 -- */
   GLM_Source (d->Vc, 0.5*g_dt, grid);
  #endif

  first_call = 0;
  return 0; /* -- step has been achieved, return success -- */
}


























