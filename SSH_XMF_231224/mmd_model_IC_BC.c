/* *********************************
  NUMBER_DENSITY_BEG      
  TEMPTURE_BEG            
  TREND_H                 
  TREND_HE                
 ********************************** */

#include "pluto.h"
#define IBEG_BOX(k,j,i) KTOT_LOOP(k)JTOT_LOOP(j)for(i=IBEG-1;i>=0;i--)
#define IEND_BOX(k,j,i) KTOT_LOOP(k)JTOT_LOOP(j)for(i=IEND+1;i<NX1_TOT;i++)
#define JBEG_BOX(k,j,i) KTOT_LOOP(k)ITOT_LOOP(i)for(j=JBEG-1;j>=0;j--)
#define JEND_BOX(k,j,i) KTOT_LOOP(k)ITOT_LOOP(i)for(j=JEND+1;j<NX2_TOT;j++)
#define KBEG_BOX(k,j,i) JTOT_LOOP(j)ITOT_LOOP(i)for(k=KBEG-1;k>=0;k--)
#define KEND_BOX(k,j,i) JTOT_LOOP(j)ITOT_LOOP(i)for(k=KEND+1;k<NX3_TOT;k++)


#define IBEG_SBOX(k,j,i) KTOT_LOOP(k)JTOT_LOOP(j)for(i=IBEG-2;i>=0;i--)
#define IEND_SBOX(k,j,i) KTOT_LOOP(k)JTOT_LOOP(j)for(i=IEND+1;i<NX1_TOT;i++)
#define JBEG_SBOX(k,j,i) KTOT_LOOP(k)ITOT_LOOP(i)for(j=JBEG-2;j>=0;j--)
#define JEND_SBOX(k,j,i) KTOT_LOOP(k)ITOT_LOOP(i)for(j=JEND+1;j<NX2_TOT;j++)
#define KBEG_SBOX(k,j,i) JTOT_LOOP(j)ITOT_LOOP(i)for(k=KBEG-2;k>=0;k--)
#define KEND_SBOX(k,j,i) JTOT_LOOP(j)ITOT_LOOP(i)for(k=KEND+1;k<NX3_TOT;k++)


static void read_dbl2D(FILE*, char*, int, int, double**);

static void outflow_BC(double ****, Grid *, int);

void mmd_model_IC(Model *model, Grid *grid, Runtime *ini)//runtime *ini
{
  int k,j,i,nv,f;
  double *x3, *x2,*x1, z,y,x, scrh, tmp;
  int mynode, totalnodes, node;
  x3 = grid[KDIR].x, x2 = grid[JDIR].x, x1 = grid[IDIR].x;
  g_gamma = g_inputParam[GAMMA];
  tmp = g_inputParam[LB_TMP];

  static int first_call = 1;
  static CSV d0;
  
  if(first_call){
    csv_read("para_files/ionrate_struct.csv", &d0);
    first_call = 0;
  }

  TOT_LOOP(k,j,i){
    x = x1[i], y = x2[j]; 
    //scrh = MAX(0.20, line_ITP(d0.data, d0.nlog, x));
    scrh = g_inputParam[H_ION_RATE];
    //IC for HI fluid
    f=HI;
    FLD_EXT(model,f,MX1,k,j,i) = 1.e-10;
    FLD_EXT(model,f,MX2,k,j,i) = 0;

    FLD_EXT(model,f,NDT,k,j,i) = (1-scrh)*(g_inputParam[LB_RHO]*UNIT_DENSITY)/ATOM_MASS(model,f)
                                 *pow((x-0.9)*10, -1.5);
    
    FLD_EXT(model,f,PRS,k,j,i) = FLD_EXT(model,f,NDT,k,j,i)*CONST_kB*tmp;


    //IC for HII fluid
    f=HII;
    FLD_EXT(model,f,MX1,k,j,i) = 1.e-10;
    FLD_EXT(model,f,MX2,k,j,i) = 0;

    FLD_EXT(model,f,NDT,k,j,i) = scrh*(g_inputParam[LB_RHO]*UNIT_DENSITY)/ATOM_MASS(model,f)
                                 *pow((x-0.9)*10, -1.5);
    FLD_EXT(model,f,PRS,k,j,i) = FLD_EXT(model,f,NDT,k,j,i)*CONST_kB*tmp;
    

    //IC for B(centered)
    MAG_EXT(model,BM1,k,j,i) = g_inputParam[LB_B]/CONST_s4PI*pow(x, -3)*cos(y);
    MAG_EXT(model,BM2,k,j,i) = g_inputParam[LB_B]/CONST_s4PI*pow(x, -3)*sin(y)*0.5;
  }


  return;
}



void mmd_model_BC(Model *model, Grid *grid)
{
  int k,j,i,nv,nf,f,b;
 
  double *x1,*x2,*x3,*dx1,*dx2,*dx3;
  double scrh, ndt_tot, u_HI[3], u_HII[3], prs_tot, x, y, z, tmp;
  double ****Vc, ****Bs, ****Bc;
  double NDT_IBEG = g_inputParam[NDT_BEG], 
         TMP_IBEG = g_inputParam[TMP_BEG], 
         ION_RATE = g_inputParam[BASE_ION_RATE];

  static double G_IBEG, grid_rate;
  static int par_dim[3] = {0, 0, 0};
  static int first_call=1;
  static CSV d0;

  if(first_call){  
    grid_rate = grid[IDIR].dx[IEND] / grid[IDIR].dx[IEND-1];
    G_IBEG =  -CONST_G*g_inputParam[M_Planet]*CONST_Mjup/UNIT_LENGTH/UNIT_LENGTH;
    D_EXPAND(par_dim[0] = grid[IDIR].nproc > 1;  ,
             par_dim[1] = grid[JDIR].nproc > 1;  ,
             par_dim[2] = grid[KDIR].nproc > 1;)

    csv_read("para_files/ionrate_struct.csv", &d0);
  }

  x1 = grid[IDIR].x, dx1 = grid[IDIR].dx;
  x2 = grid[JDIR].x, dx2 = grid[JDIR].dx;


  //inner BC
  DOM_LOOP(k,j,i){
/*
    x = x1[i], y = x2[j]; 
    //scrh = MAX(0.20, line_ITP(d0.data, d0.nlog, x));
    scrh = g_inputParam[H_ION_RATE];

    ndt_tot = FLD_EXT(model,HI,NDT,k,j,i) + FLD_EXT(model,HII,NDT,k,j,i);
    CPN_LOOP(nv){
      u_HI[nv]  = FLD_EXT(model,HI,MX1+nv,k,j,i)/FLD_EXT(model,HI,NDT,k,j,i);
      u_HII[nv] = FLD_EXT(model,HII,MX1+nv,k,j,i)/FLD_EXT(model,HII,NDT,k,j,i);
    }
    prs_tot = FLD_EXT(model,HI,PRS,k,j,i) + FLD_EXT(model,HII,PRS,k,j,i);


    FLD_EXT(model,HI,NDT,k,j,i)  = ndt_tot*(1-scrh);
    FLD_EXT(model,HII,NDT,k,j,i) = ndt_tot*scrh;
    FLD_EXT(model,HI,PRS,k,j,i)  = prs_tot*(1-scrh);
    FLD_EXT(model,HII,PRS,k,j,i) = prs_tot*scrh;


    FLD_EXT(model,HI,MX1,k,j,i)  = MAX(0,FLD_EXT(model,HI,NDT,k,j,i)*u_HI[0]);
    FLD_EXT(model,HI,MX2,k,j,i)  = FLD_EXT(model,HI,NDT,k,j,i)*u_HI[1];
    FLD_EXT(model,HII,MX1,k,j,i) = MAX(0,FLD_EXT(model,HII,NDT,k,j,i)*u_HII[0]);
    FLD_EXT(model,HII,MX2,k,j,i) = FLD_EXT(model,HII,NDT,k,j,i)*u_HII[1];
*/
    FLD_EXT(model,HI,MX1,k,j,i)  = MAX(0,FLD_EXT(model,HI,MX1,k,j,i));
    FLD_EXT(model,HII,MX1,k,j,i) = MAX(0,FLD_EXT(model,HII,MX1,k,j,i));

  }


  //Global load
  scrh = g_inputParam[H_ION_RATE];
  tmp = g_inputParam[LB_TMP];
  
  //BC for HI
  f = HI, Vc = model->fluid[f].Vc; 
  IBEG_BOX(k,j,i){
    Vc[NDT][k][j][i] = (1-scrh)*(g_inputParam[LB_RHO]*UNIT_DENSITY)/ATOM_MASS(model,f);
    EXPAND(
    Vc[MX1][k][j][i] = 0;,
    Vc[MX2][k][j][i] = 0;,
    Vc[MX3][k][j][i] = 0;
    )
    Vc[PRS][k][j][i] = FLD_EXT(model,f,NDT,k,j,i)*CONST_kB*tmp;
  }

  IEND_BOX(k,j,i){
    Vc[NDT][k][j][i] = pow(Vc[NDT][k][j][i-1],1+grid_rate) / pow(Vc[NDT][k][j][i-2],grid_rate);

    EXPAND(
    Vc[MX1][k][j][i] = MAX(0, (Vc[MX1][k][j][i-1])/Vc[NDT][k][j][i-1] * Vc[NDT][k][j][i]);,

    Vc[MX2][k][j][i] = (1+grid_rate)*Vc[MX2][k][j][IEND-1] - grid_rate*Vc[MX2][k][j][IEND-2];,
                       //( (1+grid_rate)*Vc[MX1][k][j][i-1]/Vc[NDT][k][j][i-1] 
                       //-grid_rate*Vc[MX1][k][j][i-2]/Vc[NDT][k][j][i-2])*Vc[NDT][k][j][i];,
                       
    Vc[MX3][k][j][i] = (1+grid_rate)*Vc[MX3][k][j][i-1] - grid_rate*Vc[MX3][k][j][i-2]; 
    )

    Vc[PRS][k][j][i] = pow(Vc[PRS][k][j][i-1],1+grid_rate) / pow(Vc[PRS][k][j][i-2],grid_rate);
  }

  JBEG_BOX(k,j,i){
    Vc[NDT][k][j][i] = Vc[NDT][k][JBEG][i];
    EXPAND(
    Vc[MX1][k][j][i] = Vc[MX1][k][JBEG][i];,
    Vc[MX2][k][j][i] = -Vc[MX2][k][JBEG][i];,
    Vc[MX3][k][j][i] = Vc[MX3][k][JBEG][i];
    )
    Vc[PRS][k][j][i] = Vc[PRS][k][JBEG][i];
  }

  JEND_BOX(k,j,i){
    Vc[NDT][k][j][i] = Vc[NDT][k][JEND][i];
    EXPAND(
    Vc[MX1][k][j][i] = Vc[MX1][k][JEND][i];,
    Vc[MX2][k][j][i] = -Vc[MX2][k][JEND][i];,
    Vc[MX3][k][j][i] = Vc[MX3][k][JEND][i];
    )
    Vc[PRS][k][j][i] = Vc[PRS][k][JEND][i];
  }




  //BC for HII
  f = HII, Vc = model->fluid[f].Vc; 
  IBEG_BOX(k,j,i){
    Vc[NDT][k][j][i] = scrh*(g_inputParam[LB_RHO]*UNIT_DENSITY)/ATOM_MASS(model,f);
    EXPAND(
    Vc[MX1][k][j][i] = 0;,
    Vc[MX2][k][j][i] = 0;,
    Vc[MX3][k][j][i] = 0;
    )
    Vc[PRS][k][j][i] = FLD_EXT(model,f,NDT,k,j,i)*CONST_kB*tmp;
  }

  IEND_BOX(k,j,i){
    Vc[NDT][k][j][i] = pow(Vc[NDT][k][j][i-1],1+grid_rate) / pow(Vc[NDT][k][j][i-2],grid_rate);
    EXPAND(
    Vc[MX1][k][j][i] = MAX(0, (Vc[MX1][k][j][i-1])/Vc[NDT][k][j][i-1] * Vc[NDT][k][j][i]);,

    Vc[MX2][k][j][i] = (1+grid_rate)*Vc[MX2][k][j][IEND-1] - grid_rate*Vc[MX2][k][j][IEND-2];,

    Vc[MX3][k][j][i] = (1+grid_rate)*Vc[MX3][k][j][i-1] - grid_rate*Vc[MX3][k][j][i-2];
    )
    Vc[PRS][k][j][i] = pow(Vc[PRS][k][j][i-1],1+grid_rate) / pow(Vc[PRS][k][j][i-2],grid_rate);
  }

  JBEG_BOX(k,j,i){
    Vc[NDT][k][j][i] = Vc[NDT][k][JBEG][i];
    EXPAND(
    Vc[MX1][k][j][i] = Vc[MX1][k][JBEG][i];,
    Vc[MX2][k][j][i] = -Vc[MX2][k][JBEG][i];,
    Vc[MX3][k][j][i] = Vc[MX3][k][JBEG][i];
    )
    Vc[PRS][k][j][i] = Vc[PRS][k][JBEG][i];
  }

  JEND_BOX(k,j,i){
    Vc[NDT][k][j][i] = Vc[NDT][k][JEND][i];
    EXPAND(
    Vc[MX1][k][j][i] = Vc[MX1][k][JEND][i];,
    Vc[MX2][k][j][i] = -Vc[MX2][k][JEND][i];,
    Vc[MX3][k][j][i] = Vc[MX3][k][JEND][i];
    )
    Vc[PRS][k][j][i] = Vc[PRS][k][JEND][i];
  }

  //BC for E fluid
#if ELC_ON
  Vc = model->E_fluid;
  IBEG_BOX(k,j,i){
    Vc[PRS][k][j][i] = Vc[PRS][k][j][i+1];
  }

  IEND_BOX(k,j,i){
    Vc[PRS][k][j][i] = Vc[PRS][k][j][i-1];
  }

#endif

  //BC for B field
#if MAG_ON == 1
 #ifdef STAGGERED_MHD
/*
  Bs = model->Bs; 
  IBEG_SBOX(k,j,i){
    EXPAND(
    Bs[BX1s][k][j][i] = Bs[BX1s][k][j][i+1];,
    Bs[BX2s][k][j][i] = Bs[BX2s][k][j][i+1];,
    Bs[BX3s][k][j][i] = Bs[BX3s][k][j][i+1];
    )
  }

  IEND_SBOX(k,j,i){
    EXPAND(
    Bs[BX1s][k][j][i] = Bs[BX1s][k][j][i-1];,
    Bs[BX2s][k][j][i] = Bs[BX2s][k][j][i-1];,
    Bs[BX3s][k][j][i] = Bs[BX3s][k][j][i-1];
    )
  }
*/
 #else
  Bc = model->B_field;
  IBEG_BOX(k,j,i){
    EXPAND(
    Bc[BM1][k][j][i] = g_inputParam[LB_B]/CONST_s4PI*pow(x1[i], -3)*cos(x2[j]);,
    Bc[BM2][k][j][i] = g_inputParam[LB_B]/CONST_s4PI*pow(x1[i], -3)*sin(x2[j])*0.5;,
    Bc[BM3][k][j][i] = 0;
    )
  }

  IEND_BOX(k,j,i){
    EXPAND(
    Bc[BM1][k][j][i] = (1+grid_rate)*Bc[BM1][k][j][IEND] - grid_rate*Bc[BM1][k][j][IEND-1];,
    Bc[BM2][k][j][i] = (1+grid_rate)*Bc[BM2][k][j][IEND] - grid_rate*Bc[BM2][k][j][IEND-1];,
    Bc[BM3][k][j][i] = 0;
    )
  }

  JBEG_BOX(k,j,i){
    EXPAND(
    Bc[BM1][k][j][i] = Bc[BM1][k][JBEG][i];,
    Bc[BM2][k][j][i] = -Bc[BM2][k][JBEG][i];,
    Bc[BM3][k][j][i] = 0;
    )
  }

  JEND_BOX(k,j,i){
    EXPAND(
    Bc[BM1][k][j][i] = Bc[BM1][k][JEND][i];,
    Bc[BM2][k][j][i] = -Bc[BM2][k][JEND][i];,
    Bc[BM3][k][j][i] = 0;
    )
  }

 #endif
#endif

  //BC for xuv
  for(b=0;b<model->xuv.nbeam;b++){
    model->xuv.flux[b][0][0][IEND+1] = model->xuv.sed[b];
  }



  //Data exchange between porcessors. Prototypes from PLUTO 4.2 
  #ifdef PARALLEL
   MPI_Barrier (MPI_COMM_WORLD);
   //Baryon fluids
   FLD_LOOP(f)HD_NVAR_LOOP(nv){
     AL_Exchange_dim ((char *)model->fluid[f].Vc[nv][0][0], par_dim, SZ);
   }

   #if MAG_ON == 1
   //B_field
    CPN_LOOP(nv){
      AL_Exchange_dim ((char *)model->B_field[nv][0][0], par_dim, SZ);
    }

   //Bs_field
    #ifdef STAGGERED_MHD 
     EXPAND(
      AL_Exchange_dim ((char *)(model->Bs[BX1s][0][0] - 1), par_dim, SZ_stagx);  ,
      AL_Exchange_dim ((char *)model->Bs[BX2s][0][-1]     , par_dim, SZ_stagy);  ,
      AL_Exchange_dim ((char *)model->Bs[BX3s][-1][0]     , par_dim, SZ_stagz);
     )
    #endif
   #endif

   //Electron fluid
/*
   for(b=0;b<model->xuv.nbeam;b++){
     AL_Exchange_dim((char *) model->xuv.flux[b][0][0], par_dim, SZ);
   }
*/

  MPI_Barrier (MPI_COMM_WORLD);
  #endif
 
  first_call = 0;
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
    Vc[MX1][k][j][i] = MAX(0,(1+grid_rate)*Vc[MX1][k][j][i-1] - grid_rate*Vc[MX1][k][j][i-2]);
    Vc[PRS][k][j][i] = pow(Vc[PRS][k][j][i-1],1+grid_rate) / pow(Vc[PRS][k][j][i-2],grid_rate);
  }

  first_call = 0;
  return;
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

























