#include "pluto.h"

void mmd_sm_initialize(Runtime* ini)
{
  //initialize for sm_dts
  sm_dts.cmax     = ARRAY_1D(NMAX_POINT, double);
  sm_dts.inv_dta  = 0.0;
  sm_dts.inv_dtp  = 0.5e-38; 
  sm_dts.dt_cool  = 1.e38;
  sm_dts.cfl      = ini->cfl;
  sm_dts.cfl_par  = ini->cfl_par;
  sm_dts.rmax_par = ini->rmax_par;
  sm_dts.Nsts     = sm_dts.Nrkc = 0;
  
  //initialize for sm_d
  sm_data.Vc = ARRAY_4D(HD_NVAR, NX3_TOT, NX2_TOT, NX1_TOT, double);
  sm_data.Uc = ARRAY_4D(NX3_TOT, NX2_TOT, NX1_TOT, HD_NVAR, double); 
 
  #ifdef STAGGERED_MHD
   sm_data.Vs = ARRAY_1D(DIMENSIONS, double ***);
   D_EXPAND(
     sm_data.Vs[BX1s] = ArrayBox( 0, NX3_TOT-1, 0, NX2_TOT-1,-1, NX1_TOT-1); ,
     sm_data.Vs[BX2s] = ArrayBox( 0, NX3_TOT-1,-1, NX2_TOT-1, 0, NX1_TOT-1); ,
     sm_data.Vs[BX3s] = ArrayBox(-1, NX3_TOT-1, 0, NX2_TOT-1, 0, NX1_TOT-1);)
  #endif  

  #if UPDATE_VECTOR_POTENTIAL == YES 
   D_EXPAND(                                                  ,
     sm_data.Ax3 = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);  , 
     sm_data.Ax1 = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);  
     sm_data.Ax2 = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);  
   )
  #endif

  #if RESISTIVITY != NO
   sm_data.J = ARRAY_4D(3,NX3_TOT, NX2_TOT, NX1_TOT, double);
  #endif

  sm_data.flag = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, unsigned char);

  sm_du = ARRAY_4D(3, NX3_TOT, NX2_TOT, NX1_TOT, double);
  return;
}
























