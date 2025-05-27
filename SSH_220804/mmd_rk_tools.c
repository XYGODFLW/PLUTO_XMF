#include "pluto.h"
//allcate memory to U0
//Vmx is on planning

void RK_AllocateMemory(int physics, int nptc, RK_var *var)
{
  var->physics = physics;
  var->nptc    = nptc;
  
  if(var->Vc != NULL)return;

  var->Vc  = ARRAY_4D(MHD_NVAR, NX3_TOT, NX2_TOT, NX1_TOT, double);
  //var->Vmx = ARRAY_4D(nptc, NX3_TOT, NX2_TOT, NX1_TOT, double);

  #ifdef STAGGERED_MHD
  if(physics == MHD)
    var->Vs = ARRAY_4D(DIMENSIONS, NX3_TOT, NX2_TOT, NX1_TOT, double);
  #endif

  return;
}

void RK_Copy(RK_var *var0, RK_var *var1)
//copy var0 to var1
{
  int k,j,i,nv;

  if(var1->Vc==NULL) RK_AllocateMemory(var0->physics,var0->nptc,var1);

  TOT_LOOP(k,j,i){
    if(var0->physics == MHD) MHD_NVAR_LOOP(nv)var1->Vc[nv][k][j][i] = var0->Vc[nv][k][j][i];
    if(var0->physics == HD)   HD_NVAR_LOOP(nv)var1->Vc[nv][k][j][i] = var0->Vc[nv][k][j][i];

    #ifdef STAGGERED_MHD 
    if(var1->physics == MHD)
      DIM_LOOP(nv)var1->Vs[nv][k][j][i] = var0->Vs[nv][k][j][i];   
    #endif
  }

  return;
}
 
void RK_to_F(RK_var *var0, Fluid *fluid)
//copy Vs,Vc,Vmx
{
  int k,j,i,nv;

  TOT_LOOP(k,j,i){
    if(var0->physics == MHD) MHD_NVAR_LOOP(nv)fluid->d.Vc[nv][k][j][i] = var0->Vc[nv][k][j][i];
    if(var0->physics == HD)   HD_NVAR_LOOP(nv)fluid->d.Vc[nv][k][j][i] = var0->Vc[nv][k][j][i];
    
    #ifdef STAGGERED_MHD
    if(fluid->physics == MHD)
      DIM_LOOP(nv)fluid->d.Vs[nv][k][j][i] = var0->Vs[nv][k][j][i];  
    #endif
  }
  return;
}


void F_to_RK(Fluid *fluid, RK_var *var0)
{
  int k,j,i,nv;

  if(var0->Vc==NULL) RK_AllocateMemory(fluid->physics, fluid->d.nptc,var0);

  TOT_LOOP(k,j,i){
    if(var0->physics == MHD) MHD_NVAR_LOOP(nv) var0->Vc[nv][k][j][i] = fluid->d.Vc[nv][k][j][i];
    if(var0->physics == HD)   HD_NVAR_LOOP(nv) var0->Vc[nv][k][j][i] = fluid->d.Vc[nv][k][j][i];  

    #ifdef STAGGERED_MHD
    if(var0->physics == MHD)
      DIM_LOOP(nv) var0->Vs[nv][k][j][i] = fluid->d.Vs[nv][k][j][i];  
    #endif
  }

  return;
}

void RK_MultiplyAdd(RK_var *var, RK_var *var0, double m, RK_var *var1)
//add var0 to var1 and renew var1
{
  int k,j,i,nv;
  //checking

  TOT_LOOP(k,j,i){
    if(var0->physics == MHD) MHD_NVAR_LOOP(nv) var->Vc[nv][k][j][i] = var0->Vc[nv][k][j][i]+m*var1->Vc[nv][k][j][i];
    if(var0->physics == HD)   HD_NVAR_LOOP(nv) var->Vc[nv][k][j][i] = var0->Vc[nv][k][j][i]+m*var1->Vc[nv][k][j][i];

    #ifdef STAGGERED_MHD    
    if(var->physics == MHD) DIM_LOOP(nv) var->Vs[nv][k][j][i] = var0->Vs[nv][k][j][i]+m*var1->Vs[nv][k][j][i];
    #endif
  }  
  return;  
}

void RK_sMultiplyAdd(RK_var *var0, double m, RK_var *var1)
//add var0 to var1 and renew var1
{
  int k,j,i,nv;
  //checking

  TOT_LOOP(k,j,i){
    if(var0->physics == MHD) MHD_NVAR_LOOP(nv) var0->Vc[nv][k][j][i] += m*var1->Vc[nv][k][j][i];
    if(var0->physics == HD)   HD_NVAR_LOOP(nv) var0->Vc[nv][k][j][i] += m*var1->Vc[nv][k][j][i];

    #ifdef STAGGERED_MHD    
    if(var0->physics == MHD)
      DIM_LOOP(nv) var0->Vs[nv][k][j][i] += m*var1->Vs[nv][k][j][i];
    #endif
  }  
  return;  
}


void RK_PointerLock(Fluid *fluid, RK_var *var0)
{
  int k,j,i,nv;
  var0->physics = fluid->physics;
  var0->nptc = fluid->d.nptc;  

  var0->Vc = fluid->d.Vc;
  var0->Vs = fluid->d.Vs;
  
  return;
}

double dbl_min(double a, double b)
{
  if(a<b)return a;else return b;
}

double dbl_max(double a, double b)
{
  if(a>b)return a;else return b;
}














