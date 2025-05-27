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

/* Weight factor for 2nd stage of RK integrators */

#if TIME_STEPPING == RK2
 #define w0  0.5
 #define wc  0.5
#elif TIME_STEPPING == RK3
 #define w0 0.75
 #define wc 0.25
#endif

/* ********************************************************************* */
int hd_AdvanceStep (const Data *d, Riemann_Solver *Riemann, 
                 Time_Step *Dts, Grid *grid)
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
  static int first_call = 1;
  RBox *box = GetRBox (DOM, CENTER);

  TOT_LOOP(k,j,i) d->flag[k][j][i] = 0;
  hd_UpdateStage(d, d->Uc, NULL, Riemann, g_dt, Dts, grid);

  hd_ConsToPrim3D (d->Uc, d->Vc, d->flag, box);
  return 0; /* -- step has been achieved, return success -- */
}




















