#include "pluto.h"
//note... fluid 0:HI, fluid 1:HE, fluid 2:HII, fluid 3: HE+, fluid 4: Electric
//solver aviable hd_LF_Solver, hd_Roe_Solver, hd_HLL_Solver, hd_HLLC_Solver;
//sm: state machine

static void fluid_to_prim(Fluid*, Fluid*, int nf);
static void prim_to_fluid(Fluid*, Fluid*, Model*);


void mmd_get_lf(Model* model, Grid* grid, double ****RK_lf)
{
  int f,nv,k,j,i;
  double ****Vc, ****Vc1, ****Vc2;
  double t_HII, t_HEII, t_TOT, scrh, scrh1, dt_next;

  static Fluid *fluid_prim, *fluid_new;
  static double ***wt_HII, ***wt_HEII;
  static int first_call = 1; 

  if(first_call){
    wt_HII  = ARRAY_3D(NX3_TOT,NX2_TOT,NX1_TOT,double);
    wt_HEII = ARRAY_3D(NX3_TOT,NX2_TOT,NX1_TOT,double);

    fluid_prim = ARRAY_1D(model->nfluid, Fluid);
    fluid_new  = ARRAY_1D(model->nfluid, Fluid);

    for(f = 0;f<model->nfluid;f++){
      fluid_prim[f].Vc = ARRAY_4D(MHD_NVAR,NX3_TOT,NX2_TOT,NX1_TOT,double);
      fluid_new[f].Vc = ARRAY_4D(MHD_NVAR,NX3_TOT,NX2_TOT,NX1_TOT,double);

      fluid_prim[f].f  = model->fluid[f].f;
      fluid_prim[f].ptc.charge = model->fluid[f].ptc.charge;
      fluid_prim[f].ptc.mass   = model->fluid[f].ptc.mass;
      fluid_prim[f].ptc.q_m    = model->fluid[f].ptc.q_m;
    } 
  fluid_to_prim(model->fluid, fluid_new, model->nfluid);//new_intialize
  }

  //NDT V TMP (c.g.s) -> RHO V PRS (code unit)
  
  fluid_to_prim(model->fluid, fluid_prim, model->nfluid);

  //get ion fluids charge weight this value is between 0 and 1
  TOT_LOOP(k,j,i){
   t_HII  = MODEL_EXT(model,HII,NDT,k,j,i);
   t_HEII = MODEL_EXT(model,HEII,NDT,k,j,i);
   t_TOT  = t_HII + t_HEII;

   wt_HII[k][j][i]  = t_HII/t_TOT;
   wt_HEII[k][j][i] = t_HEII/t_TOT; 
  }

  dt_next = 100;

  //Advection BEG
  //HI
  Vc = fluid_prim[HI].Vc;

  HD_NVAR_LOOP(nv)TOT_LOOP(k,j,i){
    sm_data.Vc[nv][k][j][i] = Vc[nv][k][j][i];
  }

  sm_riemann = hd_HLL_Solver; 

  hd_Integrate(&sm_data, sm_riemann, &sm_dts, grid);

  dt_next = MIN(dt_next, sm_dt);

  HD_NVAR_LOOP(nv)DOM_LOOP(k,j,i){
    FLUID_EXT(fluid_new,HI,nv,k,j,i) = sm_du[nv][k][j][i];
  }

///*

  //HEI
  Vc = fluid_prim[HEI].Vc;

  HD_NVAR_LOOP(nv)TOT_LOOP(k,j,i){
    sm_data.Vc[nv][k][j][i] = Vc[nv][k][j][i];
  }

  sm_riemann = hd_HLL_Solver; 

  hd_Integrate(&sm_data, sm_riemann, &sm_dts, grid);

  dt_next = MIN(dt_next, sm_dt);

  HD_NVAR_LOOP(nv)DOM_LOOP(k,j,i){
    FLUID_EXT(fluid_new,HEI,nv,k,j,i) = sm_du[nv][k][j][i];
  }


  //HII
  Vc = fluid_prim[HII].Vc, Vc1 = fluid_prim[ELCI].Vc;

    //partI: solve rho,vx1 
  TOT_LOOP(k,j,i){
   
    sm_data.Vc[RHO][k][j][i] = Vc[RHO][k][j][i];
    sm_data.Vc[VX1][k][j][i] = Vc[VX1][k][j][i];
    sm_data.Vc[PRS][k][j][i] = Vc[PRS][k][j][i] + wt_HII[k][j][i]*Vc1[PRS][k][j][i];
  }

  sm_riemann = hd_HLL_Solver; 

  hd_Integrate(&sm_data, sm_riemann, &sm_dts, grid);

  dt_next = MIN(dt_next, sm_dt);

  DOM_LOOP(k,j,i){
    FLUID_EXT(fluid_new,HII,RHO,k,j,i) = sm_du[RHO][k][j][i];
    FLUID_EXT(fluid_new,HII,VX1,k,j,i) = sm_du[VX1][k][j][i];
  }

    //partII: solve prs
  TOT_LOOP(k,j,i){
    scrh = Vc[PRS][k][j][i]/(Vc[PRS][k][j][i] + wt_HII[k][j][i]*Vc1[PRS][k][j][i]);

    sm_data.Vc[RHO][k][j][i] = Vc[RHO][k][j][i] * scrh;
    sm_data.Vc[VX1][k][j][i] = Vc[VX1][k][j][i];
    sm_data.Vc[PRS][k][j][i] = Vc[PRS][k][j][i];
  }

  sm_riemann = hd_HLL_Solver; 

  hd_Integrate(&sm_data, sm_riemann, &sm_dts, grid);

  dt_next = MIN(dt_next, sm_dt);

  DOM_LOOP(k,j,i){
    FLUID_EXT(fluid_new,HII,PRS,k,j,i) = sm_du[PRS][k][j][i];
  }



  //HEII
  Vc = fluid_prim[HEII].Vc, Vc1 = fluid_prim[ELCI].Vc;

  //partI: solve rho,vx1 
  TOT_LOOP(k,j,i){
   
    sm_data.Vc[RHO][k][j][i] = Vc[RHO][k][j][i];
    sm_data.Vc[VX1][k][j][i] = Vc[VX1][k][j][i];
    sm_data.Vc[PRS][k][j][i] = Vc[PRS][k][j][i] + wt_HEII[k][j][i]*Vc1[PRS][k][j][i];
  }

  sm_riemann = hd_Roe_Solver; 

  hd_Integrate(&sm_data, sm_riemann, &sm_dts, grid);

  dt_next = MIN(dt_next, sm_dt);

  DOM_LOOP(k,j,i){
    FLUID_EXT(fluid_new,HEII,RHO,k,j,i) = sm_du[RHO][k][j][i];
    FLUID_EXT(fluid_new,HEII,VX1,k,j,i) = sm_du[VX1][k][j][i];
  }

  //partII: solve prs
  TOT_LOOP(k,j,i){
    scrh = Vc[PRS][k][j][i]/(Vc[PRS][k][j][i] + wt_HEII[k][j][i]*Vc1[PRS][k][j][i]);

    sm_data.Vc[RHO][k][j][i] = Vc[RHO][k][j][i] * scrh;
    sm_data.Vc[VX1][k][j][i] = Vc[VX1][k][j][i];
    sm_data.Vc[PRS][k][j][i] = Vc[PRS][k][j][i];
  }

  sm_riemann = hd_HLL_Solver; 

  hd_Integrate(&sm_data, sm_riemann, &sm_dts, grid);

  dt_next = MIN(dt_next, sm_dt);

  DOM_LOOP(k,j,i){
    FLUID_EXT(fluid_new,HEII,PRS,k,j,i) = sm_du[PRS][k][j][i];
  }

 
  //ELECTRON fluid only solve pressure
  Vc = fluid_prim[ELCI].Vc, Vc1 = fluid_prim[HII].Vc, Vc2 = fluid_prim[HEII].Vc;

  TOT_LOOP(k,j,i){
    sm_data.Vc[RHO][k][j][i] = Vc1[RHO][k][j][i]*wt_HII[k][j][i]*Vc[PRS][k][j][i]/(Vc1[PRS][k][j][i]+wt_HII[k][j][i]*Vc[PRS][k][j][i])
                              +Vc2[RHO][k][j][i]*wt_HEII[k][j][i]*Vc[PRS][k][j][i]/(Vc2[PRS][k][j][i]+wt_HEII[k][j][i]*Vc[PRS][k][j][i]);
    sm_data.Vc[VX1][k][j][i] = Vc1[VX1][k][j][i]*wt_HII[k][j][i] + Vc2[VX1][k][j][i]*wt_HEII[k][j][i];
    sm_data.Vc[PRS][k][j][i] = Vc[PRS][k][j][i];
  }

  sm_riemann = hd_HLL_Solver; 

  hd_Integrate(&sm_data, sm_riemann, &sm_dts, grid);

  dt_next = MIN(dt_next, sm_dt);

  DOM_LOOP(k,j,i){
    FLUID_EXT(fluid_new,ELCI,PRS,k,j,i) = sm_du[PRS][k][j][i];
  }

  //Advection END
 // */

  //Get RK_Lf
  prim_to_fluid(fluid_new, fluid_prim, model);

  for(nv=0 ;nv<TOT_VAR; nv++)DOM_LOOP(k,j,i)RK_lf[nv][k][j][i] = 0;

  TOT_LOOP(k,j,i){
    HD_NVAR_LOOP(nv)for(f=0;f<(model->nfluid-1);f++)
      RK_EXT(RK_lf,f,nv,k,j,i) = FLUID_EXT(fluid_prim,f,nv,k,j,i) - MODEL_EXT(model,f,nv,k,j,i);

    RK_EXT(RK_lf,ELCI,0,k,j,i) = FLUID_EXT(fluid_prim,ELCI,TMP,k,j,i) - MODEL_EXT(model,ELCI,TMP,k,j,i);
  }

  

  g_dt = dt_next;
  first_call = 0;
  return;
}

static void fluid_to_prim(Fluid *fluid, Fluid *fluid_prim, int nf)
{
  int f,nv,k,j,i;
  const double UNIT_PRS = UNIT_DENSITY*UNIT_VELOCITY*UNIT_VELOCITY;

  for(f=0;f<nf;f++)TOT_LOOP(k,j,i){
    fluid_prim[f].Vc[PRS][k][j][i] = fluid[f].Vc[PRS][k][j][i]/UNIT_PRS;
    fluid_prim[f].Vc[VX1][k][j][i] = fluid[f].Vc[VX1][k][j][i]/UNIT_VELOCITY;
    fluid_prim[f].Vc[RHO][k][j][i] = fluid[f].Vc[NDT][k][j][i]*fluid[f].ptc.mass*CONST_amu/UNIT_DENSITY;
  }
  return;
}

static void prim_to_fluid(Fluid *fluid_prim, Fluid *fluid, Model* model)
{
  int f,nv,k,j,i;
  double scrh;
  const double UNIT_PRS = UNIT_DENSITY*UNIT_VELOCITY*UNIT_VELOCITY;

  DOM_LOOP(k,j,i){
    for(f=0;f<((model->nfluid)-1);f++){
      fluid[f].Vc[NDT][k][j][i] = fluid_prim[f].Vc[RHO][k][j][i]*UNIT_DENSITY / (fluid[f].ptc.mass*CONST_amu);
      fluid[f].Vc[VX1][k][j][i] = fluid_prim[f].Vc[VX1][k][j][i]*UNIT_VELOCITY;
      fluid[f].Vc[PRS][k][j][i] = fluid_prim[f].Vc[PRS][k][j][i]*UNIT_PRS;
    }
    scrh = FLUID_EXT(fluid,HII,NDT,k,j,i)+FLUID_EXT(fluid,HEII,NDT,k,j,i);
    FLUID_EXT(fluid,ELCI,PRS,k,j,i) = fluid_prim[ELCI].Vc[PRS][k][j][i]*UNIT_PRS;
  }

  return;
}


// printf("%d %12.4e%12.4e%12.4e\n", i,Vc[RHO][k][j][i] ,Vc[VX1][k][j][i] ,Vc[PRS][k][j][i]);//for debug
























