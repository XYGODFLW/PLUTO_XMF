#include "pluto.h"
#define RK_SELFADD(rk,a,rk0) RK_multiplyAdd(rk,1.0,rk,a,rk0);
#define RK_EQUAL(rk,a,rk0) RK_multiplyAdd(rk,a,rk0,0.0,NULL);

//formate set 1:data 2:riemann  
//riemann solvers: HD: hd_LF_Solver, hd_Roe_Solver, hd_HLL_Solver, hd_HLLC_Solver;
//                 MHD: mhd...

void get_du(Model *, Grid *, double ****);

static void RK_multiplyAdd(double ****, double , double ****, double a1, double ****);
static void RK_eqto_Model(double****, Model*);
static void Model_eqto_RK(Model*, double****);
static void RK_lockTo_Model(double ****, Model *);

static void get_drift(double* , double* );
static void uDrift(double ****, double *, int, double, int,int,int);
static void cell_load(double*, double****, int,int,int);
static void cell_pump(double****, double*, int,int,int);

const double rk_alpha[3] = {0.0, 0.75, 1.0/3.0};
const double rk_beta [3] = {1.0, 0.25, 2.0/3.0};
 
int mmd_AdvanceStep (Model *model, Grid *grid)
{
  double ****Vc;
  double grid_rate, scrh;
  int k,j,i,nv,rk_step,f;

  static double ****RK_u, ****RK_u0, ****RK_du, ****RK_out;
  static int first_call = 1;
 
  if(first_call){
    RK_u   = ARRAY_4D(TOT_VAR, NX3_TOT, NX2_TOT, NX1_TOT, double);
    RK_u0  = ARRAY_4D(TOT_VAR, NX3_TOT, NX2_TOT, NX1_TOT, double);
    RK_du  = ARRAY_4D(TOT_VAR, NX3_TOT, NX2_TOT, NX1_TOT, double);
    RK_out = ARRAY_4D(TOT_VAR, NX3_TOT, NX2_TOT, NX1_TOT, double);
  }


  //RK_beg
  RK_eqto_Model(RK_u, model);

  RK_EQUAL(RK_u0, 1.0, RK_u);

  for(rk_step=0; rk_step<3; rk_step++){
    quasi_neutrality(model);

    get_du(model, grid, RK_du);  

    //RK_EQUAL(RK_du, 0.5, RK_du);    

    RK_SELFADD(RK_du, 1.0, RK_u);
    
    RK_multiplyAdd(RK_u, rk_alpha[rk_step],RK_u0, rk_beta[rk_step],RK_du);
    
    Model_eqto_RK(model, RK_u);
  }

  //diffuse term (just for stable)
/*
  grid_rate = grid[IDIR].dx[IBEG+1]/grid[IDIR].dx[IBEG];
  scrh = 1.e-2;
  DOM_LOOP(k,j,i){
    for(f=0;f<model->nfluid;f++){
      
      MODEL_EXT(model,f,VX1,k,j,i) += scrh*(MODEL_EXT(model,f,VX1,k,j,i-1)*grid_rate+MODEL_EXT(model,f,VX1,k,j,i+1) 
                                            -(1+grid_rate)*MODEL_EXT(model,f,VX1,k,j,i) ); 

      MODEL_EXT(model,f,TMP,k,j,i) += scrh*( MODEL_EXT(model,f,TMP,k,j,i-1)*grid_rate+MODEL_EXT(model,f,TMP,k,j,i+1)
                                            -(1+grid_rate)*MODEL_EXT(model,f,TMP,k,j,i) );
    }
    //MODEL_EXT(model,HEII,NDT,k,j,i) += 0.1*( MODEL_EXT(model,HEII,NDT,k,j,i-1)*grid_rate+MODEL_EXT(model,HEII,NDT,k,j,i+1)
                                          //  -(1+grid_rate)*MODEL_EXT(model,HEII,NDT,k,j,i) );
  }
  //RK_beg
*/
  first_call = 0;
  return 0;
}



void get_du(Model *model, Grid *grid, double ****RK_du)
{
  int f,nv,k,j,i;
  static double *sigma_HI, *sigma_HEI;
  static double ****RK_lf, ****RK_lc, ****RK_u;
  static int first_call = 1;
  //implicit cell tools
  static double *cell_u0, *cell_drift;
  static double *cell_lf, *cell_lc0,  *cell_lc, *cell_du, **cell_mtx; 
    
  enum{LEFT=-1, RIGHT=1};

  if(first_call){
    RK_u = ARRAY_1D(TOT_VAR, double***);
    RK_lf = ARRAY_4D(TOT_VAR, NX3_TOT, NX2_TOT, NX1_TOT, double);
    RK_lc = ARRAY_4D(TOT_VAR, NX3_TOT, NX2_TOT, NX1_TOT, double); 
    //mtx initialize
    g_mtx_size = g_mtx_nx1 = g_mtx_nx2 = TOT_VAR;

    //cell intialize
    cell_u0 = ARRAY_1D(TOT_VAR, double);
    cell_drift = ARRAY_1D(TOT_VAR, double);

    cell_lf = ARRAY_1D(TOT_VAR, double);
    cell_lc = ARRAY_1D(TOT_VAR, double);
    cell_lc0 = ARRAY_1D(TOT_VAR,double);
    cell_du = ARRAY_1D(TOT_VAR, double);
    cell_mtx = ARRAY_2D(TOT_VAR, TOT_VAR, double);
  }
  g_mtx_size = g_mtx_nx1 = g_mtx_nx2 = TOT_VAR;

  mmd_get_lf(model, grid, RK_lf);
  get_xuv_flux(model, grid);    

  RK_lockTo_Model(RK_u, model);

  DOM_LOOP(k,j,i){

    mmd_get_lc(model, grid, RK_lc, RK_lf, k,j,i);
     
    cell_load(cell_lc0,RK_lc,k,j,i);
    cell_load(cell_lf, RK_lf,k,j,i);
    cell_load(cell_u0, RK_u, k,j,i);
    //get du
    COL_CAL(cell_du,=,cell_lc0);
    COL_CAL(cell_du,+=,cell_lf);
    
    //get lc matrix: I-dlc/du
    get_drift(cell_u0, cell_drift);
    mtx_clean(cell_mtx);
    
    for(nv=0;nv<TOT_VAR;nv++){
      //get lc(u+du)
      uDrift(RK_u, cell_u0, nv, cell_drift[nv], k,j,i);
      mmd_get_lc(model, grid, RK_lc, RK_lf, k,j,i);
      cell_load(cell_lc,RK_lc,k,j,i);
      
      //get dlc/du
      COL_CAL(cell_lc, -=, cell_lc0);
      COL_PARA(cell_lc, /=, cell_drift[nv]); 

      MTX_COL_CAL(cell_mtx,nv,-=,cell_lc);
    }
    //if((i==50) && !(g_stepNumber%10000));show_matrix(cell_mtx, "%12.4e");//for debug
    lu_solver(cell_mtx, cell_du);
    cell_pump(RK_du,cell_du,k,j,i);
  }
  
  first_call = 0;
  return;
}


//When fluid structs changed here must be changed to
static void RK_eqto_Model(double****rk, Model *model)
{
  int nv,k,j,i;

  HD_NVAR_LOOP(nv)TOT_LOOP(k,j,i){
    rk[HI*3+nv][k][j][i]   = MODEL_EXT(model,HI,nv,k,j,i);
    rk[HEI*3+nv][k][j][i]  = MODEL_EXT(model,HEI,nv,k,j,i);
    rk[HII*3+nv][k][j][i]  = MODEL_EXT(model,HII,nv,k,j,i);
    rk[HEII*3+nv][k][j][i] = MODEL_EXT(model,HEII,nv,k,j,i);
  }
  TOT_LOOP(k,j,i)rk[ELCI*3+0][k][j][i] = MODEL_EXT(model,ELCI,TMP,k,j,i);
 
  return;
}

static void Model_eqto_RK(Model *model, double****rk)
{
  int nv,k,j,i;

  HD_NVAR_LOOP(nv)TOT_LOOP(k,j,i){
    MODEL_EXT(model,HI,nv,k,j,i) = rk[HI*3+nv][k][j][i];
    MODEL_EXT(model,HEI,nv,k,j,i) = rk[HEI*3+nv][k][j][i];
    MODEL_EXT(model,HII,nv,k,j,i) = rk[HII*3+nv][k][j][i];
    MODEL_EXT(model,HEII,nv,k,j,i) = rk[HEII*3+nv][k][j][i];
  }

  
  TOT_LOOP(k,j,i)MODEL_EXT(model,ELCI,TMP,k,j,i) = rk[ELCI*3+0][k][j][i];

  return;
}

static void RK_lockTo_Model(double****rk, Model *model)
{
  int nv,k,j,i;

//Is There a good way to map Model and rk? 
  HD_NVAR_LOOP(nv){
    rk[HI*3+nv]   = model->fluid[HI].Vc[nv];
    rk[HEI*3+nv]  = model->fluid[HEI].Vc[nv];
    rk[HII*3+nv]  = model->fluid[HII].Vc[nv];
    rk[HEII*3+nv] = model->fluid[HEII].Vc[nv];
  }
  rk[ELCI*3+0] = model->fluid[ELCI].Vc[TMP];

  return;
}
//end

static void RK_multiplyAdd(double ****rk, double a0, double ****rk0, double a1, double ****rk1)
{
  int nv,k,j,i;
  if(rk1!=NULL){
    for(nv=0;nv<TOT_VAR;nv++)DOM_LOOP(k,j,i){
      rk[nv][k][j][i] = a0*rk0[nv][k][j][i]+a1*rk1[nv][k][j][i];
    }
  }else{
    for(nv=0;nv<TOT_VAR;nv++)DOM_LOOP(k,j,i){
      rk[nv][k][j][i] = a0*rk0[nv][k][j][i];
    }
  }

  return;
}

static void get_drift(double* cell_u, double* cell_drift)
{
  const double factor = 1.e-7;
  int nv;
  for(nv=0;nv<TOT_VAR;nv++){
    cell_drift[nv] = MAX(cell_u[nv]*factor,1.e-17);
  }  

  cell_drift[1] = MAX(1,fabs(cell_u[1]))*factor;
  cell_drift[4] = MAX(1,fabs(cell_u[4]))*factor;
  cell_drift[7] = MAX(1,fabs(cell_u[7]))*factor;
  cell_drift[10] = MAX(1,fabs(cell_u[10]))*factor;

  return;
}

static void uDrift(double ****rku, double *u0, int var, double drift, int k,int j, int i)
/* **********************************
 *  out\ u dirft:
 *  int\ rku0, var, drift
 *
 *  #discription
 *  u from rku0, in var, 
 *  drift dirction: left or right 
 * **********************************/
{
  cell_pump(rku,u0,k,j,i);
  
  rku[var][k][j][i] += drift;

  return;
}

static void cell_load(double* cell_var, double****rk, int k, int j, int i)
{
  int nv;
  for(nv=0;nv<TOT_VAR;nv++){
    cell_var[nv] = rk[nv][k][j][i];
  }
  return;
}

static void cell_pump(double****rk, double* cell_var, int k, int j, int i)
{
  int nv;
  for(nv=0;nv<TOT_VAR;nv++){
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















