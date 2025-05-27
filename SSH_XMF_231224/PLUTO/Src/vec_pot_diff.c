/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file 
  \brief Compute magnetic field from vector potential.
  
  The function VectorPotentialDiff() computes either staggered or 
  cell-center magnetic field components by suitable central differencing
  of the vector potential, \f$ \vec{B} = \nabla\times\vec{A} \f$.
  - In Cartesian geometry:
    \f[
    B_x = \pd{A_z}{y} - \pd{A_y}{z} \,,\quad
    B_y = \pd{A_x}{z} - \pd{A_z}{x} \,,\quad
    B_z = \pd{A_y}{x} - \pd{A_x}{y}
    \f]
 
  - In cylindrical geometry:
    \f[
     B_R = \frac{1}{R}\pd{A_z}{\phi}    - \pd{A_\phi}{z}      \,,\quad
     B_z = \frac{1}{R}\left(\pd{(RA_\phi)}{R} - \pd{A_R}{\phi}\right)
    \f]
 
  - In polar geometry:
    \f[
    B_R    = \frac{1}{R}\pd{A_z}{\phi} - \pd{A_\phi}{z}  \,,\quad
    B_\phi = \pd{A_R}{z} - \pd{A_z}{R}             \,,\quad
    B_z    = \frac{1}{R}\pd{(RA_\phi)}{R} - \frac{1}{R}\pd{A_R}{\phi}
    \f]
 
  - In spherical geometry:
    \f[ 
    B_r      =   \frac{1}{r\sin\theta}\pd{(\sin\theta A_\phi)}{\theta} 
               - \frac{1}{r\sin\theta}\pd{A_\theta}{\phi}           \,,\quad
    B_\theta =   \frac{1}{r\sin\theta}\pd{A_r}{\phi} 
               - \frac{1}{r}\pd{(rA_\phi)}{r}                         \,,\quad
    B_\phi   =   \frac{1}{r}\pd{(rA_\theta)}{r}                  
               - \frac{1}{r}\pd{A_r}{\theta}
    \f]
 
  For cell-centered MHD vector potential is compute at the cell-center.
  In the case of staggered MHD, the position of A is edge-centered and
  it is shown below:
  \verbatim
            ______________________
           /                     /|  
          /                     / |  
         /       z face        /  |
        /                     Ax  |
       /                     /    |
      /                     /     Az
     ----------Ay-----------      |
     |                     |      |
     |                     |   y  |
     |                     | face |
     |                     |     / 
     Az       x face      Az    /  
     |                     |   Ax  
     |                     |  /     
     |                     | /
     |                     |/
     ----------Ay-----------
  \endverbatim
 
 
  \authors A. Mignone (mignone@ph.unito.it)\n
           P. Tzeferacos (petros.tzeferacos@ph.unito.it)
  \date   Sep 24, 2012
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

#if PHYSICS == MHD || PHYSICS == RMHD
/* ********************************************************************* */
void VectorPotentialDiff (double *b, int i, int j, int k, Grid *grid)
/*!
 * Assign face- or cell-centered magnetic field by differentiating
 * the vector potential.
 *
 * \param [out] b  array of magnetic field starting at 0, 
 *               \f$ B_{x_1} = b[0]\,, B_{x_2} = b[1]\,, B_{x_3} = b[2] \f$
 * \param [in]  i  the cell index in the first coordinate direction
 * \param [in]  j  the cell index in the first coordinate direction
 * \param [in]  k  the cell index in the first coordinate direction
 * \param [in] grid pointer to an array of Grid structures
 *
 *********************************************************************** */
{
//OFF in multi-fluid module  
}
#endif /* PHYSICS == MHD || PHYSICS == RMHD */











