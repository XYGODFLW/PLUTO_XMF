/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief Perform index permutation and set domain integration indexes.

  The function SetIndexes() performs two different tasks, depending on
  the value of ::g_dir and the time-stepping algorithm:
  - perform a cyclic permutation of the array indexes corresponding to
    vector components (velocity, momentum and magnetic field);
  - set the starting and final grid indexes (IBEG, IEND,...) of the 
    computational domain before commencing integration.

  The indexes coincide with the usual ones (IBEG, IEND,...) most of the time,
  but can be expanded one further zone either in the transverse or normal
  directions or both. This depends on the integration algorithm.

  With cell-centered fields, the following table holds:
  \verbatim
                        RK        CTU
                   +-------------------------------------------
    Transverse++   |    NO        YES, at predictor step
                   |
    Normal++       |    NO        NO
  \endverbatim
 
  If constrained transport is enabled, then
  \verbatim
                        RK        CTU
                   +-------------------------------------------
    Transverse++   |    YES       YES
                   |
    Normal++       |    NO        YES, at predictor step
  \endverbatim
  
  Predictor step (g_intStage = 1) in CTU requires extra transverse loop
  to obtain transverse predictor in any case.
  Also, with CTU+CT, one needs to expand the grid of one zone in the
  \e normal direction as well.
  This allows to computed fully corner coupled states in the boundary to
  get electric field components during the constrained transport algorithm. 

  \author A. Mignone (mignone@ph.unito.it)\n
  \date   Sep 17, 2012
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ********************************************************************* */
void hd_SetIndexes (Index *indx, Grid *grid)
/*!
 * Set vector indices and integration index range.
 *
 * \param [out] indx pointer to an Index structure
 * \param [in]  grid pointer to an array of Grid structures
 *
 *********************************************************************** */
{
  IBEG = grid[IDIR].lbeg; IEND = grid[IDIR].lend;
  JBEG = grid[JDIR].lbeg; JEND = grid[JDIR].lend;
  KBEG = grid[KDIR].lbeg; KEND = grid[KDIR].lend;

  if (g_dir == IDIR) {   /* -- Order: X-Y-Z  {in,t1,t2 = i,j,k} -- */

    EXPAND(VXn = MXn = VX1; , 
           VXt = MXt = VX2; , 
           VXb = MXb = VX3;)


    indx->ntot   = grid[IDIR].np_tot;
    indx->beg    = IBEG; indx->end    = IEND;
    indx->t1_beg = JBEG; indx->t1_end = JEND;
    indx->t2_beg = KBEG; indx->t2_end = KEND;

  }else if (g_dir == JDIR){ /* -- Order: Y-X-Z  {in,t1,t2 = j,i,k} -- */

    EXPAND(VXn = MXn = VX2;  , 
           VXt = MXt = VX1;  , 
           VXb = MXb = VX3;)

    
    indx->ntot   = grid[JDIR].np_tot;
    indx->beg    = JBEG; indx->end    = JEND;
    indx->t1_beg = IBEG; indx->t1_end = IEND;
    indx->t2_beg = KBEG; indx->t2_end = KEND;

  }else if (g_dir == KDIR){ /* -- Order: Z-X-Y  {in,t1,t2 = k,i,j} -- */

    VXn = MXn = VX3;
    VXt = MXt = VX1;
    VXb = MXb = VX2;
    
    indx->ntot   = grid[KDIR].np_tot;
    indx->beg    = KBEG; indx->end    = KEND;
    indx->t1_beg = IBEG; indx->t1_end = IEND;
    indx->t2_beg = JBEG; indx->t2_end = JEND;

  }

/* -------------------------------------------------------
    Expand grid one further zone to account for proper 
    flux computation. This is necessary to obtain the EMF 
    in the boundary zones and to get transverse rhs for 
    corner coupled states.
   ------------------------------------------------------- */

   #ifdef CTU 
    if (g_intStage == 1){
      #if (PARABOLIC_FLUX & EXPLICIT)
       indx->beg--;
       indx->end++; 
      #endif
      D_EXPAND(                                 ,
               indx->t1_beg--; indx->t1_end++;  ,
               indx->t2_beg--; indx->t2_end++;) 
    }
   #endif
}








