#include "pluto.h"
#define RK_SELFADD(rk,a,rk0) RK_multiplyAdd(rk,1.0,rk,a,rk0);
#define RK_EQUAL(rk,a,rk0) RK_multiplyAdd(rk,a,rk0,0.0,NULL);

//formate set 1:data 2:riemann  
//riemann solvers: HD: hd_LF_Solver, hd_Roe_Solver, hd_HLL_Solver, hd_HLLC_Solver;
//                 MHD: mhd...

void get_du(Model *, Grid *, double ****);

static void get_drift(double* , double* );
static void cell_load(double*, double****, int,int,int);
static void cell_pump(double****, double*, int,int,int);

const double rk_alpha[3] = {0.0, 0.75, 1.0/3.0};
const double rk_beta [3] = {1.0, 0.25, 2.0/3.0};
 
int mmd_AdvanceStep (Model *model, Grid *grid)
{
  double ****Vc;
  double grid_rate, scrh;
  int k,j,i,nv,rk_step,f,d;

  static double ****RK_u, ****RK_u0, ****RK_du, ****RK_out;
  static int first_call = 1;
 
  if(first_call){
    RK_u   = ARRAY_4D(RK_TOT, NX3_TOT, NX2_TOT, NX1_TOT, double);
    RK_u0  = ARRAY_4D(RK_TOT, NX3_TOT, NX2_TOT, NX1_TOT, double);
    RK_du  = ARRAY_4D(RK_TOT, NX3_TOT, NX2_TOT, NX1_TOT, double);
    RK_out = ARRAY_4D(RK_TOT, NX3_TOT, NX2_TOT, NX1_TOT, double);
  }

  //RK_beg

  RK_eqto_Model(RK_u, model);
  
  RK_EQUAL(RK_u0, 1.0, RK_u);

  for(rk_step=0; rk_step<3; rk_step++){
    #if ELC_ON
    quasi_neutrality(model);
    #endif

    mmd_model_BC(model, grid);

    get_du(model, grid, RK_du);  

    RK_SELFADD(RK_du, 1.0, RK_u);

    RK_multiplyAdd(RK_u, rk_alpha[rk_step],RK_u0, rk_beta[rk_step], RK_du);

    Model_eqto_RK(model, RK_u);
  }

  g_dt = g_dt_temp;

  first_call = 0;
  return 0;
}



void get_du(Model *model, Grid *grid, double ****RK_du)
{

  int f,nv,k,j,i,d;
  static double *sigma_HI, *sigma_HEI;
  static double ****RK_lf, ****RK_lc, ****RK_u;
  static int first_call = 1;
  //implicit cell tools
  static double *cell_u0, *cell_drift;
  static double *cell_lf, *cell_lc0,  *cell_lc, *cell_du, **cell_mtx; 
    
  enum{LEFT=-1, RIGHT=1};

  
  if(first_call){
    RK_u = ARRAY_1D(RK_TOT, double***);
    RK_lf = ARRAY_4D(RK_TOT, NX3_TOT, NX2_TOT, NX1_TOT, double);
    RK_lc = ARRAY_4D(RK_TOT, NX3_TOT, NX2_TOT, NX1_TOT, double); 
    //mtx initialize
    g_mtx_size = g_mtx_nx1 = g_mtx_nx2 = IM_TOT;

    //cell intialize
    cell_u0 = ARRAY_1D(RK_TOT, double);
    cell_drift = ARRAY_1D(RK_TOT, double);

    cell_lf = ARRAY_1D(RK_TOT, double);
    cell_lc = ARRAY_1D(RK_TOT, double);
    cell_lc0 = ARRAY_1D(RK_TOT,double);
    cell_du = ARRAY_1D(RK_TOT, double);
    cell_mtx = ARRAY_2D(RK_TOT, RK_TOT, double);
  }

  mmd_get_lf(model, grid, RK_lf);


  get_xuv_flux(model, grid);


//Implicit upgrade
  RK_lockTo_Model(RK_u, model);
  DOM_LOOP(k,j,i){
    mmd_get_lc(model, grid, RK_lc, RK_lf, k,j,i);    
    //if((j==125)&(i==125))printf("Check %.2e, \n", RK_EXT(RK_lc, HII*HD_NVAR+VX1, k,j,i));

    cell_load(cell_lc0,RK_lc,k,j,i);

    cell_load(cell_lf, RK_lf,k,j,i);
    
    cell_load(cell_u0, RK_u ,k,j,i);
    
    //get du
    COL_CAL(cell_du,=,cell_lc0);
    COL_CAL(cell_du,+=,cell_lf);

    //get lc matrix: I-dlc/du
    get_drift(cell_u0, cell_drift); //get drift according to u0
    mtx_clean(cell_mtx);

    IM_VAR_LOOP(nv){
      RK_u[nv][k][j][i] += cell_drift[nv];//RK_u has locked to Model

      mmd_get_lc(model, grid, RK_lc, RK_lf, k,j,i);

      cell_load(cell_lc,RK_lc,k,j,i);
      //get dlc/du
      COL_CAL(cell_lc, -=, cell_lc0);
      COL_PARA(cell_lc, /=, cell_drift[nv]); 

      MTX_COL_CAL(cell_mtx,nv,-=,cell_lc);

      cell_pump(RK_u, cell_u0,k,j,i);//RK_u reset to cell_u0
    }


    lu_solver(cell_mtx, cell_du);

    cell_pump(RK_du,cell_du,k,j,i);
  }

  //Get dB of RK_du (in both staggered or centered magnetic-field)
  CPN_LOOP(nv)TOT_LOOP(k,j,i){
    RK_EXT(RK_du, RK_MAG+nv,k,j,i) = RK_EXT(RK_lf, RK_MAG+nv,k,j,i);
  }


  first_call = 0;
  return;
}



//When fluid structs changed here must be changed too
void RK_eqto_Model(double****rk, Model *model)
{

  int nv,k,j,i,f,d;

  FLD_LOOP(f)HD_NVAR_LOOP(nv)TOT_LOOP(k,j,i){
    rk[f*HD_NVAR+nv][k][j][i] = FLD_EXT(model,f,nv,k,j,i);
  }

  #if ELC_ON
  TOT_LOOP(k,j,i) rk[RK_ELC][k][j][i] = ELC_EXT(model,PRS,k,j,i);
  #endif
  
  #if MAG_ON
   #ifdef Staggered_MHD
    CPN_LOOP(nv)TOT_LOOP(k,j,i) 
      rk[RK_MAG+nv][k][j][i] = MGS_EXT(model,nv,k,j,i);
   
   #else
    CPN_LOOP(nv)TOT_LOOP(k,j,i) 
      rk[RK_MAG+nv][k][j][i] = MAG_EXT(model,nv,k,j,i);
   
   #endif
  #endif

  return;
}


void Model_eqto_RK(Model *model, double****rk)
{

  int nv,k,j,i,f,d;

  FLD_LOOP(f)HD_NVAR_LOOP(nv)TOT_LOOP(k,j,i){
    FLD_EXT(model,f,nv,k,j,i) = rk[f*HD_NVAR+nv][k][j][i];
  }

  //+Quasi-neutrality

  #if ELC_ON
  TOT_LOOP(k,j,i) model->E_fluid[PRS][k][j][i] = rk[RK_ELC][k][j][i];
  #endif

  #if MAG_ON
   #ifdef Staggered_MHD
    CPN_LOOP(nv)TOT_LOOP(k,j,i) MGS_EXT(model,nv,k,j,i) = rk[RK_MAG+nv][k][j][i];
   #else
    CPN_LOOP(nv)TOT_LOOP(k,j,i) MAG_EXT(model,nv,k,j,i) = rk[RK_MAG+nv][k][j][i];
   #endif
  #endif

  return;
}


void RK_lockTo_Model(double****rk, Model *model)
{

  int nv,k,j,i,f;
  //Is there a better way to map Model and rk? 
  FLD_LOOP(f)HD_NVAR_LOOP(nv){
    rk[f*HD_NVAR+nv] = model->fluid[f].Vc[nv];
  }

  #if ELC_ON
  rk[RK_ELC] = model->E_fluid[PRS];
  #endif


  #if MAG_ON
   #ifdef Staggered_MHD
    CPN_LOOP(nv) rk[RK_MAG+nv] = model->Bs_field[nv];
   #else
    CPN_LOOP(nv) rk[RK_MAG+nv] = model->B_field[nv];
   #endif
  #endif
  return;
}
//end


void RK_multiplyAdd(double ****rk, double a0, double ****rk0, double a1, double ****rk1)
{
  int nv,k,j,i;
  if(rk1!=NULL){
    for(nv=0;nv<RK_TOT;nv++)TOT_LOOP(k,j,i){
      rk[nv][k][j][i] = a0*rk0[nv][k][j][i]+a1*rk1[nv][k][j][i];
    }
  }else{
    for(nv=0;nv<RK_TOT;nv++)TOT_LOOP(k,j,i){
      rk[nv][k][j][i] = a0*rk0[nv][k][j][i];
    }
  }

  return;
}


static void get_drift(double* cell_u, double* cell_drift)
{
  const double factor = 1.e-5;
  int nv, f;

  IM_VAR_LOOP(nv){
    cell_drift[nv] = MAX(cell_u[nv]*factor, 1.e-17);
  }  

  FLD_LOOP(f){//Adjust for Velocities  
    CPN_LOOP(nv)
      cell_drift[f*HD_NVAR+VX1+nv] = MAX(1,fabs(cell_u[f*HD_NVAR+VX1+nv]))*factor;
  }

  return;
}


static void cell_load(double* cell_var, double****rk, int k, int j, int i)
{
  
  int nv;
  IM_VAR_LOOP(nv){
    cell_var[nv] = rk[nv][k][j][i];
  }
  return;
}

static void cell_pump(double****rk, double* cell_var, int k, int j, int i)
{
  int nv;
  IM_VAR_LOOP(nv){
    rk[nv][k][j][i] = cell_var[nv];
  }
  return;
}


/* ********* RK4 Version *********** */
/* 
  const double rk_alpha[5] = {256, 1.0/6.0, 2.0/6.0, 2.0/6.0, 1.0/6.0};
  const double rk_beta [5] = {256, 0.0, 0.5, 0.5, 1.0};

//RK BEG
  RK_lockTo_Model(RK_u, model);
  RK_EQUAL(RK_u0, 1.0, RK_u);
  RK_EQUAL(RK_out,1.0, RK_u);

  get_du(model, grid, RK_du);
  for(rk_step=1; rk_step<=4; rk_step++){

    RK_multiplyAdd(RK_u, 1.0,RK_u0, rk_beta[rk_step],RK_du);

    quasi_neutrality(model);
    get_du(model, grid, RK_du);

    RK_SELFADD(RK_out, rk_alpha[rk_step],RK_du);
  }
  RK_EQUAL(RK_u, 1.0, RK_out);
//RK END
*/















