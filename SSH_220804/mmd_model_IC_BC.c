/* *********************************
  NUMBER_DENSITY_BEG      
  TEMPTURE_BEG            
  TREND_H                 
  TREND_HE                
 ********************************** */

#include "pluto.h"
#define IBEG_BOX(k,j,i) for(k=0,j=0,i=IBEG-1;i>=0;i--)
#define IEND_BOX(k,j,i) for(k=0,j=0,i=IEND+1;i<NX1_TOT;i++)

#define PRESSURE(k,j,i) CONST_kB*Vc[NDT][k][j][i]*Vc[TMP][k][j][i]
#define MOMENTON(k,j,i) Vc[NDT][k][j][i]*Vc[VX1][k][j][i]
#define KINETIC_ENERGY(k,j,i) Vc[NDT][k][j][i]*Vc[VX1][k][j][i]*Vc[VX1][k][j][i]
#define DENSITY_FLOW(k,j,i) (Vc[NDT][k][j][i]*Vc[VX1][k][j][i]*x1[i]*x1[i])
#define V_PROFILE(am,trend,k,j,i) 1.e6/(x1[i]*x1[i]*IC_decline(x1[i],trend,am))* \
                            (5*5*IC_decline(5,trend,am));

#define IC_SET(f,trend,am,NDT_COEF) \
    Vc = model->fluid[f].Vc; \
    DOM_LOOP(k,j,i){  \
      Vc[NDT][k][j][i] = (NDT_COEF)*(NDT_IBEG)*IC_decline(x1[i], trend, am); \
      Vc[VX1][k][j][i] = 0; \
      Vc[PRS][k][j][i] = (NDT_COEF)*(NDT_IBEG)*CONST_kB*TMP_IBEG*IC_decline(x1[i], trend, am); \
    }
//V_PROFILE(am,trend,k,j,i);

#define VX1_IBEG(am) 0
#define ID_SIZE (IEND-IBEG+1)
#define OUTPUT_NVAR (3*(model->nfluid)+2)
//V_PROFILE(am,k,j,i); 


//double TMP_IBEG, NDT_IBEG;

//NDT_IBEG = g_inputParam[NUMBER_DENSITY_BEG];
//TMP_IBEG = g_inputParam[TEMPTURE_BEG];

static void read_dbl2D(FILE*, char*, int, int, double**);

static void outflow_BC(double ****, Grid *, int);
static double IC_decline(double, double, double);

void mmd_model_IC(Model *model, Grid *grid, Runtime *ini)//runtime *ini
{
  int k,j,i,b,nv,f;
  double *x1,*x2,*x3, *dx1;
  double ****Vc;

  double NDT_IBEG = g_inputParam[NUMBER_DENSITY_BEG], 
         TMP_IBEG = g_inputParam[TEMPTURE_BEG], 
         ION_RATE = g_inputParam[BASE_ION_RATE];

  static double *VX1_profile_NEU, *VX1_profile_ION;
  //load pointer
  x3 = grid[KDIR].x, x2 = grid[JDIR].x, x1 = grid[IDIR].x;

  //set globals and consts
  g_gamma = g_inputParam[GAMMA];
  

  if(VX1_profile_NEU==NULL)VX1_profile_NEU = ARRAY_1D(NX1_TOT, double);
  if(VX1_profile_ION==NULL)VX1_profile_ION = ARRAY_1D(NX1_TOT, double);

if(g_inputParam[IC_FROM_FILE]){
  
  #ifdef PARALLEL
  int mynode, totalnodes, node; 
  FILE *fp;
  static double **mtx;
  MPI_Status status;
  
  MPI_Comm_size(MPI_COMM_WORLD, &totalnodes);
  MPI_Comm_rank(MPI_COMM_WORLD, &mynode);

  //allocate memory for every core
  if(mynode==0){
    // read file
    if((fp = fopen(ini->IC_file, "r"))==NULL){printf("ERROR: !IC_FILE:%s do not exist, please check your settings\n",ini->IC_file); QUIT_PLUTO(1);}
    mtx = ARRAY_2D(ID_SIZE*totalnodes, OUTPUT_NVAR, double);
  }else{
    mtx = ARRAY_2D(ID_SIZE, OUTPUT_NVAR, double);
  }

  //sed data to other processors
  if(mynode==0){
    read_dbl2D(fp, "TMP_ELC", ID_SIZE*totalnodes, OUTPUT_NVAR, mtx);
    for(node=1;node<totalnodes;node++){
      MPI_Send(mtx[node*ID_SIZE], ID_SIZE*OUTPUT_NVAR, MPI_DOUBLE, node, 1, MPI_COMM_WORLD);
    }
  }else{
    MPI_Recv(mtx[0], ID_SIZE*OUTPUT_NVAR, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, &status);
  }

  DOM_LOOP(k,j,i)for(f=0;f<model->nfluid;f++)for(nv=0;nv<HD_NVAR;nv++){
    MODEL_EXT(model,f,nv,k,j,i) = mtx[i-IBEG][2+3*f+nv];
  }
  DOM_LOOP(k,j,i)for(f=0;f<model->nfluid;f++){
    MODEL_EXT(model,f,PRS,k,j,i) = CONST_kB*MODEL_EXT(model,f,NDT,k,j,i)*MODEL_EXT(model,f,TMP,k,j,i);
  }
  
  FreeArray2D(mtx);
  #endif
}else{
  //IC for f0 HI
  IC_SET(HI,g_inputParam[TREND_H], 1.0, 1-g_inputParam[HE_RATIO]);
  
  //IC for f1 HE
  IC_SET(HEI,g_inputParam[TREND_HE], g_inputParam[HE_SH_COEF], g_inputParam[HE_RATIO]);
  
  //IC for f2 HII
  IC_SET(HII,g_inputParam[TREND_ION], 1.0, (1-g_inputParam[HE_RATIO])*ION_RATE);
  
  //IC for f3 HE+
  IC_SET(HEII,g_inputParam[TREND_ION],g_inputParam[HE_SH_COEF],g_inputParam[HE_RATIO]*ION_RATE);
  
  //IC for f4 ELC
  Vc = model->fluid[ELCI].Vc;
  DOM_LOOP(k,j,i){
    Vc[NDT][k][j][i] = MODEL_EXT(model,HII,NDT,k,j,i) + MODEL_EXT(model,HEII,NDT,k,j,i);
    Vc[VX1][k][j][i] = 0;
    Vc[TMP][k][j][i] = (MODEL_EXT(model,HII,NDT,k,j,i) + MODEL_EXT(model,HEII,NDT,k,j,i)) * CONST_kB*TMP_IBEG;  
  }
}
  
  
  //IC for xuv flux
  DOM_LOOP(k,j,i){
    for(b=0;b<model->xuv.nbeam;b++){
      model->xuv.flux[b][k][j][i]=0;
    }
  }

  return;
}






void mmd_model_BC(Model *model, Grid *grid)
{
  int k,j,i,nv,nf,f,b;
 
  double *x1,*x2,*x3,*dx1,scrh;
  double ****Vc;
  double NDT_IBEG = g_inputParam[NUMBER_DENSITY_BEG], 
         TMP_IBEG = g_inputParam[TEMPTURE_BEG], 
         ION_RATE = g_inputParam[BASE_ION_RATE];
  double grid_rate = grid[IDIR].x[IBEG] / grid[IDIR].x[IBEG+1];

  static double G_IBEG;
  static int par_dim[3] = {0, 0, 0};
  static int first_call=1;
  if(first_call){  
    grid_rate = grid[IDIR].x[IBEG] / grid[IDIR].x[IBEG+1];
    G_IBEG =  -CONST_G*g_inputParam[M_Planet]*CONST_Mjup/UNIT_LENGTH/UNIT_LENGTH;
    D_EXPAND(par_dim[0] = grid[IDIR].nproc > 1;  ,
             par_dim[1] = grid[JDIR].nproc > 1;  ,
             par_dim[2] = grid[KDIR].nproc > 1;)
  }  
  x1 = grid[IDIR].x, dx1 = grid[IDIR].dx;

  //Inner boundary condition..just for calculation stable
  
/* stop here */
   DOM_LOOP(k,j,i){


    if(g_time<20){
      MODEL_EXT(model,HI,VX1,k,j,i) = MAX(0,MODEL_EXT(model,HI,VX1,k,j,i));
      MODEL_EXT(model,HEI,VX1,k,j,i) = MAX(0,MODEL_EXT(model,HEI,VX1,k,j,i));
      MODEL_EXT(model,HII,VX1,k,j,i) = MAX(0,MODEL_EXT(model,HII,VX1,k,j,i));
      MODEL_EXT(model,HEII,VX1,k,j,i) = MAX(0,MODEL_EXT(model,HEII,VX1,k,j,i));
    }

    scrh = g_inputParam[TEMPTURE_BEG];
    if(grid[IDIR].x[i]<2){
      MODEL_EXT(model,HI,PRS,k,j,i) = MAX(MODEL_EXT(model,HI,PRS,k,j,i), scrh*MODEL_EXT(model,HI,NDT,k,j,i)*CONST_kB);
      MODEL_EXT(model,HEI,PRS,k,j,i) = MAX(MODEL_EXT(model,HEI,PRS,k,j,i), scrh*MODEL_EXT(model,HEI,NDT,k,j,i)*CONST_kB);
      MODEL_EXT(model,HII,PRS,k,j,i) = MAX(MODEL_EXT(model,HII,PRS,k,j,i), scrh*MODEL_EXT(model,HII,NDT,k,j,i)*CONST_kB);
      MODEL_EXT(model,HEII,PRS,k,j,i) = MAX(MODEL_EXT(model,HEII,PRS,k,j,i), scrh*MODEL_EXT(model,HEII,NDT,k,j,i)*CONST_kB);
      MODEL_EXT(model,ELCI,PRS,k,j,i) = MAX(MODEL_EXT(model,ELCI,PRS,k,j,i), scrh*MODEL_EXT(model,ELCI,NDT,k,j,i)*CONST_kB);
    }

    if(grid[IDIR].x[i]<1.3){
      MODEL_EXT(model,HEII,NDT,k,j,i) = MAX(10000, MODEL_EXT(model,HEII,NDT,k,j,i));
      MODEL_EXT(model,HEII,PRS,k,j,i) = MAX(MODEL_EXT(model,HEII,PRS,k,j,i), scrh*MODEL_EXT(model,HEII,NDT,k,j,i)*CONST_kB);
    }
/*
    if(grid[IDIR].x[i]<2.0){
      MODEL_EXT(model,HEI,NDT,k,j,i) = MAX(MODEL_EXT(model,HEI,NDT,k,j,i),1.e-5*MODEL_EXT(model,HI,NDT,k,j,i));
      MODEL_EXT(model,HEI,PRS,k,j,i) = MAX(MODEL_EXT(model,HEI,PRS,k,j,i),1.e-5*MODEL_EXT(model,HI,PRS,k,j,i));
    }

    if( (scrh=MODEL_EXT(model,HI,NDT,k,j,i)+MODEL_EXT(model,HEI,NDT,k,j,i))>1.e12 ){
      MODEL_EXT(model,HI,NDT,k,j,i) = scrh*0.9;
      MODEL_EXT(model,HEI,NDT,k,j,i) = scrh*0.1;
    }
    //MODEL_EXT(model,HEI,VX1,k,j,i) = MAX(0,MODEL_EXT(model,HEI,VX1,k,j,i));
*/
  }

  //for 2D or 3D case... loading box is nessesary
  //BC for f0 HI
  Vc = model->fluid[0].Vc; 
  IBEG_BOX(k,j,i){
    Vc[NDT][k][j][i] = NDT_IBEG*(1-g_inputParam[HE_RATIO]);
    Vc[VX1][k][j][i] = Vc[VX1][k][j][i+1];
    Vc[TMP][k][j][i] = Vc[NDT][k][j][i]*CONST_kB*TMP_IBEG;
  }

  outflow_BC(Vc,grid,0);

  //BC for f1 HE
  Vc = model->fluid[1].Vc; 
  IBEG_BOX(k,j,i){
    Vc[NDT][k][j][i] = NDT_IBEG*g_inputParam[HE_RATIO];
    Vc[VX1][k][j][i] = Vc[VX1][k][j][i+1];
    Vc[TMP][k][j][i] = Vc[NDT][k][j][i]*CONST_kB*TMP_IBEG;
  }

  outflow_BC(Vc,grid,0);

  //BC for f2 HII
  Vc = model->fluid[2].Vc; 
  IBEG_BOX(k,j,i){
    Vc[NDT][k][j][i] = MAX(1.e5,MODEL_EXT(model,HII,NDT,k,j,i+1));
    Vc[VX1][k][j][i] = Vc[VX1][k][j][i+1];
    Vc[TMP][k][j][i] = Vc[NDT][k][j][i]*CONST_kB*TMP_IBEG;
  }

  outflow_BC(Vc,grid,0);

  //BC for f3 HE+
  Vc = model->fluid[3].Vc; 
  IBEG_BOX(k,j,i){
    Vc[NDT][k][j][i] = MAX(1.e4,MODEL_EXT(model,HEII,NDT,k,j,i+1));
    Vc[VX1][k][j][i] = Vc[VX1][k][j][i+1];
    Vc[TMP][k][j][i] = Vc[NDT][k][j][i]*CONST_kB*TMP_IBEG;
  }

  outflow_BC(Vc,grid,0);

  //BC for f4 ELC
  Vc = model->fluid[4].Vc; 
  IBEG_BOX(k,j,i){
    Vc[NDT][k][j][i] = MODEL_EXT(model,HII,NDT,k,j,i) + MODEL_EXT(model,HEII,NDT,k,j,i);
    Vc[VX1][k][j][i] = Vc[VX1][k][j][i+1];
    Vc[TMP][k][j][i] = Vc[NDT][k][j][i]*CONST_kB*TMP_IBEG;
  }

  outflow_BC(Vc,grid,0);
  first_call = 0;
  
  //BC for xuv
  if(g_time>g_inputParam[FXUV_TIME])
    for(b=0;b<model->xuv.nbeam;b++){
      model->xuv.flux[b][0][0][IEND+1] = model->xuv.sed[b];
    }
  else{
    for(b=0;b<model->xuv.nbeam;b++){
      model->xuv.flux[b][0][0][IEND+1] = model->xuv.sed[b]*g_inputParam[FXUV_INI];
    }
  }
 
  //exchange data in boundary #prototype from PLUTO 4.2 
  

  #ifdef PARALLEL
   MPI_Barrier (MPI_COMM_WORLD);

   for (f = 0;f<model->nfluid;f++){
     for (nv = 0; nv < HD_NVAR; nv++){
       AL_Exchange_dim ((char *)model->fluid[f].Vc[nv][0][0], par_dim, SZ);
   }}

   for(b=0;b<model->xuv.nbeam;b++){
     AL_Exchange_dim((char *) model->xuv.flux[b][0][0], par_dim, SZ);
   }
   //B will be here
   //MPI_Barrier (MPI_COMM_WORLD);
  #endif

  return;
}


static void outflow_BC(double ****Vc, Grid *grid, int side)
{
  int k,j,i;
  
  double *x1;
  static int first_call = 1;
  static double grid_rate; 

  if(first_call){
    grid_rate = grid[IDIR].x[IEND] / grid[IDIR].x[IEND-1];
  }  
  x1 = grid[IDIR].x;

  IEND_BOX(k,j,i){
    Vc[NDT][k][j][i] = pow(Vc[NDT][k][j][i-1],1+grid_rate) / pow(Vc[NDT][k][j][i-2],grid_rate);
    Vc[VX1][k][j][i] = MAX(0,(1+grid_rate)*Vc[VX1][k][j][i-1] - grid_rate*Vc[VX1][k][j][i-2]);
    Vc[PRS][k][j][i] = pow(Vc[PRS][k][j][i-1],1+grid_rate) / pow(Vc[PRS][k][j][i-2],grid_rate);
  }

  first_call = 0;
  return;
}

static double IC_decline(double x1, double incline, double A_mu)
{
  double stop_point = g_inputParam[STOP_POINT];
  double scale_H;

  scale_H = CONST_kB*g_inputParam[TEMPTURE_BEG]/(CONST_G*g_inputParam[M_Planet]*CONST_Mjup/UNIT_LENGTH/UNIT_LENGTH*CONST_amu*A_mu);

  if(x1<stop_point){
    return exp(-(x1-1)*UNIT_LENGTH/scale_H);
  }else{
    return exp(-(stop_point-1)*UNIT_LENGTH/scale_H)* pow(x1-0.99,incline)*pow(100,stop_point*incline)/(pow(stop_point-0.999,incline)*pow(1000,stop_point*incline));
  }
}

static void read_dbl2D(FILE* fp, char* term, int nj, int ni, double** output)
{
  char s[128];  
  int j,i;

  rewind(fp);
  while(1){
    fscanf(fp,"%s",s);
    if(!strcmp(s,term)){
      for(j=0;j<nj;j++)for(i=0;i<ni;i++) fscanf(fp, "%lf", output[j]+i);
      return;
    }
    if(feof(fp)){printf("> error: term %s is NG\n", term); return;}
  }
}

























