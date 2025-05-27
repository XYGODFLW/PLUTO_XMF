#include"pluto.h"

/* ********************************************************************* */
double mhd_NextTimeStep (Time_Step *Dts, Runtime *ini, Grid *grid)
/*!
 * Compute and return the time step for the next time level
 * using the information from the previous integration
 * (Dts->inv_dta and Dts->inv_dp).
 *
 * \param [in] Dts    pointer to the Time_Step structure
 * \param [in] ini    pointer to the Runtime structure
 * \param [in] grid   pointer to array of Grid structures
 *
 * \return The time step for next time level
 *********************************************************************** */
{
  int idim;
  double dt_adv, dt_par, dtnext;
  double dxmin;
  double xloc, xglob;

/* ---------------------------------------------------
   1. Take the maximum of inv_dt across all processors
   --------------------------------------------------- */

  #ifdef PARALLEL
   xloc = Dts->inv_dta;
   MPI_Allreduce (&xloc, &xglob, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
   Dts->inv_dta = xglob;
   #if (PARABOLIC_FLUX != NO)
    xloc = Dts->inv_dtp;
    MPI_Allreduce (&xloc, &xglob, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    Dts->inv_dtp = xglob;
   #endif
   #if COOLING != NO
    xloc = Dts->dt_cool;
    MPI_Allreduce (&xloc, &xglob, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    Dts->dt_cool = xglob;
   #endif
  #endif

/* ----------------------------------
   2. Compute time step
   ---------------------------------- */

  #if (PARABOLIC_FLUX & EXPLICIT)
   dt_adv  = 1.0/(Dts->inv_dta + 2.0*Dts->inv_dtp);
  #else
   dt_adv  = 1.0/Dts->inv_dta;
  #endif
  dt_adv *= ini->cfl;
  dtnext  = dt_adv;

/* -------------------------------------------------------
   3. Maximum propagation speed for the local processor.
      Global glm_ch will be computed later in GLM_Init.
   ------------------------------------------------------- */

  #ifdef GLM_MHD
   dxmin = grid[IDIR].dl_min;
   for (idim = 1; idim < DIMENSIONS; idim++){ /*  Min cell length   */
     dxmin = MIN(dxmin, grid[idim].dl_min);
   }
   glm_ch = ini->cfl*dxmin/dtnext;
  #endif

/* ---------------------------------------------------------
   4. With STS, the ratio between advection (full) and 
      parabolic time steps should not exceed ini->rmax_par.
   --------------------------------------------------------- */
      
  #if (PARABOLIC_FLUX & SUPER_TIME_STEPPING) || (PARABOLIC_FLUX & RK_CHEBYSHEV)
   dt_par  = ini->cfl_par/(2.0*Dts->inv_dtp);
   dtnext *= MIN(1.0, ini->rmax_par/(dt_adv/dt_par));
  #endif

/* ----------------------------------
   5. Compute Cooling time step
   ---------------------------------- */

  #if COOLING != NO
   dtnext = MIN(dtnext, Dts->dt_cool);
  #endif
   
/* --------------------------------------------------------------
    6. Allow time step to vary at most by a factor 
       ini->cfl_max_var.
       Quit if dt gets too small, issue a warning if first_dt has
       been overestimated.
   -------------------------------------------------------------- */


/* --------------------------------------------
   7. Reset time step coefficients
   -------------------------------------------- */

  DIM_LOOP(idim) Dts->cmax[idim] = 0.0;
  Dts->inv_dta = 0.0;
  Dts->inv_dtp = 0.5e-38;
  Dts->dt_cool = 1.e38;

  return(dtnext);
}

