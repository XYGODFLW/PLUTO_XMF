#include "pluto.h"
/* **************************
   notes for advection terms:

   fluids   :HI, HE, HII, HE+, Electric
   HDsolver :hd_LF_Solver, hd_Roe_Solver, hd_HLL_Solver, hd_HLLC_Solver;
   MHDsolver:mhd_HLL_Solver, mhd_HLLC_Solver, mhd_HLLD_Solver, mhd_LF_Solver, mhd_Roe_Solver
   sm       : state machine

 **************************** */

void mmd_get_lf(Model* model, Grid* grid, double ****RK_lf)
{

  int f,nv,k,j,i;
  double ****Vc, ****Vs;
  double t_HII, t_HEII, t_TOT, scrh, scrh1, dt_next, k2_b2, *U, wt_B;
  Fluid *fluid;

  static Model VAR_out;
  static int first_call = 1; 

  if(first_call){
  //allocate memoray for VAR_out (with rho, momentum, prs, B)
  VAR_out.fluid = ARRAY_1D(NFLD, Fluid);
  fluid = VAR_out.fluid;
  FLD_LOOP(f){
    fluid[f].Vc = ARRAY_4D(HD_NVAR, NX3_TOT, NX2_TOT, NX1_TOT, double);
    fluid[f].f  = model->fluid[f].f;
    fluid[f].ptc.charge = model->fluid[f].ptc.charge;
    fluid[f].ptc.mass   = model->fluid[f].ptc.mass;
    fluid[f].ptc.q_m    = model->fluid[f].ptc.q_m;
  }
  VAR_out.E_fluid = ARRAY_4D(HD_NVAR, NX3_TOT, NX2_TOT, NX1_TOT, double);
  VAR_out.B_field = ARRAY_4D(COMPONENTS, NX3_TOT, NX2_TOT, NX1_TOT, double);

  #ifdef STAGGERED_MHD
  VAR_out.Bs = ARRAY_1D(DIMENSIONS, double ***);
  EXPAND(
    VAR_out.Bs[BX1s] = ArrayBox( 0, NX3_TOT-1, 0, NX2_TOT-1,-1, NX1_TOT-1); ,
    VAR_out.Bs[BX2s] = ArrayBox( 0, NX3_TOT-1,-1, NX2_TOT-1, 0, NX1_TOT-1); ,
    VAR_out.Bs[BX3s] = ArrayBox(-1, NX3_TOT-1, 0, NX2_TOT-1, 0, NX1_TOT-1);)
  #endif
  }

  g_dt_temp = 100;



  //Advection BEG
  //Adv for HI 
  f = HI; wt_B = 0.0;

  TOT_LOOP(k,j,i){
    //load Uc
    U = sm_data.Uc[k][j][i];

    U[RHO] = FLD_EXT(model,f,NDT,k,j,i)*ATOM_MASS(model,f)/UNIT_DENSITY;
    k2_b2 = 0;
    CPN_LOOP(nv){
      U[MX1+nv] = FLD_EXT(model,f,MX1+nv,k,j,i)/UNIT_MMT;
      U[BX1+nv] = MAG_EXT(model,nv,k,j,i)/UNIT_B;
      k2_b2 += 0.5*(U[MX1+nv]*U[MX1+nv]/U[RHO] + U[BX1+nv]*U[BX1+nv]*wt_B);
    }
    U[ENG] = FLD_EXT(model,f,PRS,k,j,i)/(g_gamma-1)/UNIT_ENG + k2_b2;

    //load Vc
    Vc = sm_data.Vc;
    
    Vc[RHO][k][j][i] = U[RHO];
    CPN_LOOP(nv){
      Vc[VX1+nv][k][j][i] = U[MX1+nv]/U[RHO];
      Vc[BX1+nv][k][j][i] = U[BX1+nv]*wt_B;
    }
    Vc[PRS][k][j][i] = FLD_EXT(model,f,PRS,k,j,i)/UNIT_ENG;
  }

  //Solving
  sm_riemann = mhd_HLL_Solver;

  mhd_Integrate(&sm_data, sm_riemann, &sm_dts, grid);

  g_dt_temp = MIN(g_dt_temp, sm_dt);

  //load U of state_machine to var_out
  f = HI, wt_B = 0.0;
  DOM_LOOP(k,j,i){
    U = sm_data.Uc[k][j][i];

    FLD_EXT(&VAR_out,f,NDT,k,j,i) = U[RHO]*UNIT_DENSITY/ATOM_MASS(model,f);
    k2_b2 = 0;
    CPN_LOOP(nv){
      FLD_EXT(&VAR_out,f,MX1+nv,k,j,i) = U[MX1+nv]*UNIT_MMT;
      k2_b2 += 0.5*(U[MX1+nv]*U[MX1+nv]/U[RHO] + U[BX1+nv]*U[BX1+nv]*wt_B);
    }
    FLD_EXT(&VAR_out,f,PRS,k,j,i) = (U[ENG]-k2_b2)*(g_gamma-1)*UNIT_ENG;
  }



  //Adv for HII

  //load model(c,g,s) to Vc and Uc in sm_data(dimentionless in pluto)
  f = HII; wt_B = 1.0;

  TOT_LOOP(k,j,i){
    //load Uc
    U = sm_data.Uc[k][j][i];

    U[RHO] = FLD_EXT(model,f,NDT,k,j,i)*ATOM_MASS(model,f)/UNIT_DENSITY;
    k2_b2 = 0;
    CPN_LOOP(nv){
      U[MX1+nv] = FLD_EXT(model,f,MX1+nv,k,j,i)/UNIT_MMT;
      U[BX1+nv] = MAG_EXT(model,nv,k,j,i)/UNIT_B;
      k2_b2 += 0.5*(U[MX1+nv]*U[MX1+nv]/U[RHO] + U[BX1+nv]*U[BX1+nv]);
    }
    U[ENG] = FLD_EXT(model,f,PRS,k,j,i)/(g_gamma-1)/UNIT_ENG + k2_b2;

    //load Vc
    Vc = sm_data.Vc;
    
    Vc[RHO][k][j][i] = U[RHO];
    CPN_LOOP(nv){
      Vc[VX1+nv][k][j][i] = U[MX1+nv]/U[RHO];
      Vc[BX1+nv][k][j][i] = U[BX1+nv];
    }
    Vc[PRS][k][j][i] = FLD_EXT(model,f,PRS,k,j,i)/UNIT_ENG;
  }
  
  //Solving
  sm_riemann = mhd_HLL_Solver;

  mhd_Integrate(&sm_data, sm_riemann, &sm_dts, grid);

  g_dt_temp = MIN(g_dt_temp, sm_dt);


  //load U of state_machine to var_out
  DOM_LOOP(k,j,i){
    U = sm_data.Uc[k][j][i];

    FLD_EXT(&VAR_out,f,NDT,k,j,i) = U[RHO]*UNIT_DENSITY/ATOM_MASS(model,f);
    k2_b2 = 0;
    CPN_LOOP(nv){
      FLD_EXT(&VAR_out,f,MX1+nv,k,j,i) = U[MX1+nv]*UNIT_MMT;
      k2_b2 += 0.5*(U[MX1+nv]*U[MX1+nv]/U[RHO] + U[BX1+nv]*U[BX1+nv]*wt_B);
    }
    FLD_EXT(&VAR_out,f,PRS,k,j,i) = (U[ENG]-k2_b2)*(g_gamma-1)*UNIT_ENG;
  }

  DOM_LOOP(k,j,i){ //B field
    U = sm_data.Uc[k][j][i];

    CPN_LOOP(nv){
      MAG_EXT(&VAR_out,nv,k,j,i) = U[BX1+nv]*UNIT_B;
    }
  }

  //Advection END

//Get RK_Lf
  RK_VAR_LOOP(nv)DOM_LOOP(k,j,i)RK_lf[nv][k][j][i] = 0;

  //Lf of baryon fluid
  FLD_LOOP(f)HD_NVAR_LOOP(nv)DOM_LOOP(k,j,i) 
    RK_EXT(RK_lf, HD_NVAR*f+nv, k,j,i) = FLD_EXT(&VAR_out,f,nv,k,j,i) - FLD_EXT(model,f,nv,k,j,i);

  //Lf of electron fluid
  #if ELC_ON
  #endif

  //Lf of B field
  #if MAG_ON
   #ifdef STAGGERED_MHD
    CPN_LOOP(nv)TOT_LOOP(k,j,i) RK_EXT(RK_lf,RK_MAG+nv,k,j,i) = sm_data.Vs[BX1][k][j][i] - MGS_EXT(model,nv,k,j,i);
   #else
    CPN_LOOP(nv)TOT_LOOP(k,j,i) RK_EXT(RK_lf,RK_MAG+nv,k,j,i) = sm_data.Uc[k][j][i][BX1+nv] - MAG_EXT(model,nv,k,j,i);
   #endif
  #endif

  //Check Part

  first_call = 0;
  return;
}




















