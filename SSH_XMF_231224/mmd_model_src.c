/* *********************************
  Source terms for N,V,P of neutral, 
  ion and electron fluids which needs 
  implicit treatment.
   
  2023-2-23 by Lei Xing (PLUTO-XMF)
 *********************************** */

#define MMT_TRANS(f1,f2, VAR) coef[f1][f2]*\
                              (FLD_EXT(model,f2,MX1+(VAR),k,j,i)/ATOM_MASS_AM(model,f2)*FLD_EXT(model,f1,NDT,k,j,i)-\
                               FLD_EXT(model,f1,MX1+(VAR),k,j,i)/ATOM_MASS_AM(model,f1)*FLD_EXT(model,f2,NDT,k,j,i) );

#define GET_TMP(f,k,j,i) (FLD_EXT(model,f,PRS,k,j,i)/FLD_EXT(model,f,NDT,k,j,i)/CONST_kB)

#include "pluto.h"

void mmd_get_lc(Model *model, Grid *grid, double ****RK_lc, double****RK_lf, int k, int j, int i)
{

  int f,nv,b;

  static double G_COST, TIDAL_COST, tot_ndt;
  static double ion_ndt[32], rec_ndt[32], coef[32][32],  heating, cooling, scrh, dv_G, source[32], tid_f;
  //static double tmp_reduce[32][32], mass_reduce[32][32],
  static double *sigma_HI, *sigma_HEI;
  static int first_call = 1;
  static double *x1, *x2, *x3;

  x1 = grid[IDIR].x;
  x2 = grid[JDIR].x;
  tid_f = 1.5;//non rotating fram is 1.5, rotating fram is 2.5

  RK_VAR_LOOP(nv)RK_lc[nv][k][j][i] = 0;


  if(first_call){
    //gravity and tidal
    G_COST = CONST_G*g_inputParam[M_Planet]*CONST_Mjup/UNIT_LENGTH/UNIT_LENGTH;// G*M/r2
    TIDAL_COST = CONST_G*(g_inputParam[M_Star]*CONST_Msun + g_inputParam[M_Planet]*CONST_Mjup) 
                 *pow(CONST_au*g_inputParam[R_Orbit],-3)*UNIT_LENGTH;
  }

  dv_G = -G_COST/grid[IDIR].x[i]/grid[IDIR].x[i]*G_DT;
         //*(g_inputParam[G_INI]+(g_inputParam[G_FINAL]-g_inputParam[G_INI])
         //*MIN(1, g_time/g_inputParam[G_TIME]))*G_DT;
  f = HI;
  RK_EXT(RK_lc, f*HD_NVAR+MX1, k,j,i)   += (dv_G + TIDAL_COST*x1[i]*(tid_f*pow(sin(x2[j]),2)-1))
                                            * FLD_EXT(model,f,NDT,k,j,i) * ATOM_MASS(model,f);

  RK_EXT(RK_lc, f*HD_NVAR+MX2, k,j,i)   += 0.5*TIDAL_COST*x1[i]*tid_f*sin(2*x2[j])
                                            * FLD_EXT(model,f,NDT,k,j,i) * ATOM_MASS(model,f);

  f = HII;
  RK_EXT(RK_lc, f*HD_NVAR+MX1, k,j,i)  += (dv_G + TIDAL_COST*x1[i]*(tid_f*pow(sin(x2[j]),2)-1))
                                            * FLD_EXT(model,f,NDT,k,j,i) * ATOM_MASS(model,f);

  RK_EXT(RK_lc, f*HD_NVAR+MX2, k,j,i)  += 0.5*TIDAL_COST*x1[i]*tid_f*sin(2*x2[j])
                                            * FLD_EXT(model,f,NDT,k,j,i) * ATOM_MASS(model,f);

  //collision between HI and HII
  scrh = 0.5*(GET_TMP(HI,k,j,i)+GET_TMP(HII,k,j,i));
  coef[HI][HII] = coef[HII][HI] = 2.67e-10*pow(scrh, 0.5)*pow(1-0.083*log10(scrh) ,2) * G_DT;

  CPN_LOOP(nv){
    RK_EXT(RK_lc, HD_NVAR*HI+MX1+nv, k, j, i) += MMT_TRANS(HI, HII, nv);
    RK_EXT(RK_lc, HD_NVAR*HII+MX1+nv, k, j, i) += MMT_TRANS(HII, HI, nv);
  }

  first_call = 0;
  return;
}



  /*
  #define CHE_PRS_TRANS(f1,f2,nloss,nget) (-nloss/MODEL_EXT(model,f1,NDT,k,j,i) * MODEL_EXT(model,f1,PRS,k,j,i) \
                                          + nget/MODEL_EXT(model,f2,NDT,k,j,i) * MODEL_EXT(model,f2,PRS,k,j,i)\
                                        +1./3.*nget*ATOM_MASS(model,f1)*pow(MODEL_EXT(model,f2,VX1,k,j,i)-MODEL_EXT(model,f1,VX1,k,j,i),2) )

  #define COL_PRS_TRANS(f1,f2) (2*(MODEL_EXT(model,f1,NDT,k,j,i)*MODEL_EXT(model,f2,PRS,k,j,i) \
                               -MODEL_EXT(model,f2,NDT,k,j,i)*MODEL_EXT(model,f1,PRS,k,j,i) )  \
                             +2./3.*ATOM_MASS(model,f2)*MODEL_EXT(model,f1,NDT,k,j,i)*MODEL_EXT(model,f2,NDT,k,j,i)  \
                               *pow(MODEL_EXT(model,f1,VX1,k,j,i)-MODEL_EXT(model,f2,VX1,k,j,i) ,2) )

  #define CAM(f) MODEL_EXT(model,f,NDT,k,j,i)/(MODEL_EXT(model,f,NDT,k,j,i)+RK_EXT(RK_lc,f,NDT,k,j,i)+RK_EXT(RK_lf,f,NDT,k,j,i))


  #define T_TRANS(f1,f2) 0

  int f,nv,b;

  static double G_COST, TIDAL_COST, tot_ndt;
  static double ion_ndt[32], rec_ndt[32], coef[32][32],  heating, cooling, scrh, dv_G, source[32];//tmp_reduce[32][32], mass_reduce[32][32],
  static double *sigma_HI, *sigma_HEI;
  static int first_call = 1;
 
  if(first_call){
    //gravity
    G_COST = -CONST_G*g_inputParam[M_Planet]*CONST_Mjup/UNIT_LENGTH/UNIT_LENGTH;// G*M/r2*dt
    TIDAL_COST = 3*CONST_G*g_inputParam[M_Star]*CONST_Msun* pow(CONST_au*g_inputParam[R_Orbit],-3)*UNIT_LENGTH;
  }
  //initialize
  quasi_neutrality(model);
  sigma_HI = model->xuv.sigma[HI], sigma_HEI = model->xuv.sigma[HEI];
  
  //RK_lc = 0

  for(nv=0;nv<RK_TOT;nv++) RK_lc[nv][k][j][i] = 0;
 

  //printf("G_DT = %12.4e\n",G_DT);

  //gravity
  dv_G = G_COST/grid[IDIR].x[i]/grid[IDIR].x[i]*G_DT
         *(g_inputParam[G_INI]+(g_inputParam[G_FINAL]-g_inputParam[G_INI])*MIN(1, g_time/g_inputParam[G_TIME]));
  
  dv_G += TIDAL_COST*grid[IDIR].x[i] * G_DT;


  RK_EXT(RK_lc,HI,VX1,k,j,i)   += dv_G;
  RK_EXT(RK_lc,HEI,VX1,k,j,i)  += dv_G;
  RK_EXT(RK_lc,HII,VX1,k,j,i)  += dv_G;
  RK_EXT(RK_lc,HEII,VX1,k,j,i) += dv_G; 


  //get every term
  ion_ndt[HI] = 0, ion_ndt[HEI] = 0, heating = 0;
  rec_ndt[HII] = 0, rec_ndt[HEII] = 0;

  //photo_ion_source and heating
  BEAM_LOOP(b,model->xuv.nbeam){
    scrh = sigma_HI[b] * FLUX_EXT(model,b,k,j,i)*MODEL_EXT(model,HI,NDT,k,j,i)*G_DT; 
    ion_ndt[HI]  += scrh/model->xuv.hmu[b];
    heating += scrh;

    scrh = sigma_HEI[b] * FLUX_EXT(model,b,k,j,i)*MODEL_EXT(model,HEI,NDT,k,j,i)*G_DT; 
    ion_ndt[HEI]  += scrh/model->xuv.hmu[b];
    heating += scrh;
  }


  cooling = -7.5e-19*MODEL_EXT(model,HI,NDT,k,j,i)*MODEL_EXT(model,ELCI,NDT,k,j,i)*exp(-118348/GET_TMP(ELC,k,j,i))*G_DT;

  //charge exchange in HE/H+ and H/HE+
    //set 1
  scrh = MIN(1.0,exp(g_time-20));

  ion_ndt[HI]  += scrh*1.25e-15*pow(300/TMP_REDUCE(model,HI,HEII,k,j,i), -0.25)
                  *MODEL_EXT(model,HI,NDT,k,j,i)*MODEL_EXT(model,HEII,NDT,k,j,i) *G_DT;

  rec_ndt[HEII] += scrh*1.25e-15*pow(300/TMP_REDUCE(model,HI,HEII,k,j,i), -0.25)
                   *MODEL_EXT(model,HI,NDT,k,j,i)*MODEL_EXT(model,HEII,NDT,k,j,i) *G_DT;
  
  ion_ndt[HEI] += scrh*1.75e-11*pow(300/TMP_REDUCE(model,HII,HEI,k,j,i), 0.75)*exp(-128000/TMP_REDUCE(model,HII,HEI,k,j,i))
                  *MODEL_EXT(model,HII,NDT,k,j,i)*MODEL_EXT(model,HEI,NDT,k,j,i) *G_DT;

  rec_ndt[HII]  += scrh*1.75e-11*pow(300/TMP_REDUCE(model,HII,HEI,k,j,i), 0.75)*exp(-128000/TMP_REDUCE(model,HII,HEI,k,j,i))
                   *MODEL_EXT(model,HII,NDT,k,j,i)*MODEL_EXT(model,HEI,NDT,k,j,i) *G_DT; 

    //set 2
    //...

  //electron collision ionization
  scrh = 13.6*CONST_eV/(1.5*CONST_kB*GET_TMP(ELC,k,j,i));
  ion_ndt[HI]  += 2.91e-8*pow(scrh,0.39)/(0.232+scrh)*exp(-scrh)*MODEL_EXT(model,HI,NDT,k,j,i)*MODEL_EXT(model,ELCI,NDT,k,j,i)*G_DT;
  scrh = 24.6*CONST_eV/(1.5*CONST_kB*GET_TMP(ELC,k,j,i));
  ion_ndt[HEI] += 1.75e-8*pow(scrh,0.35)/(0.180+scrh)*exp(-scrh)*MODEL_EXT(model,HEI,NDT,k,j,i)*MODEL_EXT(model,ELCI,NDT,k,j,i)*G_DT;

  //reconbinations
  rec_ndt[HII] += 4.0e-12*pow(300/GET_TMP(ELC,k,j,i),0.64) *MODEL_EXT(model,HII,NDT,k,j,i)*MODEL_EXT(model,ELCI,NDT,k,j,i) *G_DT; 
  rec_ndt[HEII] += 4.6e-12*pow(300/GET_TMP(ELC,k,j,i),0.64) *MODEL_EXT(model,HEII,NDT,k,j,i)*MODEL_EXT(model,ELCI,NDT,k,j,i) *G_DT; 

  for(f=0,tot_ndt=0; f<model->nfluid; f++) tot_ndt += MODEL_EXT(model,f,NDT,k,j,i);    

  //basic chemistry.. 

  RK_EXT(RK_lc,HI,NDT,k,j,i) += -ion_ndt[HI] + rec_ndt[HII];
 
  RK_EXT(RK_lc,HI,VX1,k,j,i) +=  rec_ndt[HII]/MODEL_EXT(model,HI,NDT,k,j,i) * DIFF_VX1(model,HII,HI,k,j,i);

  RK_EXT(RK_lc,HI,PRS,k,j,i) +=  CHE_PRS_TRANS(HI,HII,ion_ndt[HI],rec_ndt[HII]);
                                
                                 


  RK_EXT(RK_lc,HEI,NDT,k,j,i) += -ion_ndt[HEI] + rec_ndt[HEII];

  RK_EXT(RK_lc,HEI,VX1,k,j,i) += rec_ndt[HEII]/(MODEL_EXT(model,HEI,NDT,k,j,i)) * DIFF_VX1(model,HEII,HEI,k,j,i);

  RK_EXT(RK_lc,HEI,PRS,k,j,i) += CHE_PRS_TRANS(HEI,HEII,ion_ndt[HEI],rec_ndt[HEII]);
                                 
                                  


  RK_EXT(RK_lc,HII,NDT,k,j,i) += ion_ndt[HI] - rec_ndt[HII];

  RK_EXT(RK_lc,HII,VX1,k,j,i) += ion_ndt[HI]/(MODEL_EXT(model,HII,NDT,k,j,i)) * DIFF_VX1(model,HI,HII,k,j,i);
    
  RK_EXT(RK_lc,HII,PRS,k,j,i) += CHE_PRS_TRANS(HII,HI,rec_ndt[HII],ion_ndt[HI]);
                                 



  RK_EXT(RK_lc,HEII,NDT,k,j,i) += ion_ndt[HEI] - rec_ndt[HEII];

  RK_EXT(RK_lc,HEII,VX1,k,j,i) += ion_ndt[HEI]/MODEL_EXT(model,HEII,NDT,k,j,i) * DIFF_VX1(model,HEI,HEII,k,j,i);

  RK_EXT(RK_lc,HEII,PRS,k,j,i) += CHE_PRS_TRANS(HEII,HEI,rec_ndt[HEII],ion_ndt[HEI]);
                                 
    

  RK_EXT(RK_lc,ELCI,0,k,j,i) += (g_gamma-1)*(g_inputParam[HEATING_EFFICIENCY]*heating+cooling);
  

  //charge exchange and collision
  scrh = TMP_REDUCE(model,HI,HEI,k,j,i);
  coef[HI][HEI]   = CONST_kB*scrh/(1.04e18*pow(scrh, 0.732) * CONST_amu) * G_DT;
             
  scrh = 0.5*(GET_TMP(HI,k,j,i)+GET_TMP(HII,k,j,i));
  coef[HI][HII]   = 2.67e-10*pow(scrh, 0.5) * G_DT*pow(1-0.083*log10(scrh) ,2);
  
  scrh = 0.5*(GET_TMP(HEI,k,j,i)+GET_TMP(HEII,k,j,i));
  coef[HEI][HEII] = 3.50e-10*pow(scrh, 0.5) * G_DT*pow(1-0.093*log10(scrh) ,2);
  
  scrh = TMP_REDUCE(model,HII,HEII,k,j,i);  
  coef[HII][HEII] = 1.15*pow(scrh,-1.5)*G_DT;
  
  coef[HII][ELCI] = 0.0299*pow(GET_TMP(ELC,k,j,i),-1.5)*G_DT;
 
  coef[HEII][ELCI] = 0.0299*pow(GET_TMP(ELC,k,j,i),-1.5)*G_DT;


  RK_EXT(RK_lc,HI,NDT,k,j,i) += 0;
 
  RK_EXT(RK_lc,HI,VX1,k,j,i) += coef[HI][HII]/ATOM_MASS_AM(model,HI) * MODEL_EXT(model,HII,NDT,k,j,i) * DIFF_VX1(model,HII,HI,k,j,i)+
                                coef[HI][HEI]/ATOM_MASS_AM(model,HI) * MODEL_EXT(model,HEI,NDT,k,j,i) * DIFF_VX1(model,HEI,HI,k,j,i);
                         
  RK_EXT(RK_lc,HI,PRS,k,j,i) += coef[HI][HII]/(ATOM_MASS_AM(model,HI)+ATOM_MASS_AM(model,HII)) * COL_PRS_TRANS(HI,HII) +
                                coef[HI][HEI]/(ATOM_MASS_AM(model,HI)+ATOM_MASS_AM(model,HEI)) * COL_PRS_TRANS(HI,HEI);

  RK_EXT(RK_lc,HEI,NDT,k,j,i) += 0;

  RK_EXT(RK_lc,HEI,VX1,k,j,i) += coef[HEI][HEII]/ATOM_MASS_AM(model,HEI) *MODEL_EXT(model,HEII,NDT,k,j,i)*DIFF_VX1(model,HEII,HEI,k,j,i)+
                                 coef[HI][HEI]  /ATOM_MASS_AM(model,HEI) *MODEL_EXT(model,HI,NDT,k,j,i) * DIFF_VX1(model,HI,HEI,k,j,i);

  RK_EXT(RK_lc,HEI,PRS,k,j,i) += coef[HEI][HEII]/(ATOM_MASS_AM(model,HEI)+ATOM_MASS_AM(model,HEII)) * COL_PRS_TRANS(HEI,HEII) +
                                 coef[HI][HEI]/(ATOM_MASS_AM(model,HEI)+ATOM_MASS_AM(model,HI)) * COL_PRS_TRANS(HEI,HI);

  RK_EXT(RK_lc,HII,NDT,k,j,i) += 0.0;

  RK_EXT(RK_lc,HII,VX1,k,j,i) += coef[HI][HII]  /ATOM_MASS_AM(model,HII) *MODEL_EXT(model,HI,NDT,k,j,i)* DIFF_VX1(model,HI,HII,k,j,i)+
                                 coef[HII][HEII]/ATOM_MASS_AM(model,HII) *MODEL_EXT(model,HEII,NDT,k,j,i) * DIFF_VX1(model,HEII,HII,k,j,i);//
    
  RK_EXT(RK_lc,HII,PRS,k,j,i) += coef[HI][HII]/(ATOM_MASS_AM(model,HII)+ATOM_MASS_AM(model,HI)) * COL_PRS_TRANS(HII,HI) +
                                 coef[HII][HEII]/(ATOM_MASS_AM(model,HII)+ATOM_MASS_AM(model,HEII))* COL_PRS_TRANS(HII,HEII)+
                                 coef[HII][ELCI]/(ATOM_MASS_AM(model,HII)+ATOM_MASS_AM(model,ELCI))* COL_PRS_TRANS(HII,ELCI)
                                 ;


  RK_EXT(RK_lc,HEII,NDT,k,j,i) += 0.0;

  RK_EXT(RK_lc,HEII,VX1,k,j,i) += coef[HEI][HEII]/ATOM_MASS_AM(model,HEII)* MODEL_EXT(model,HEI,NDT,k,j,i)*DIFF_VX1(model,HEI,HEII,k,j,i)+
                                  coef[HII][HEII]/ATOM_MASS_AM(model,HEII)* MODEL_EXT(model,HII,NDT,k,j,i)*DIFF_VX1(model,HII,HEII,k,j,i);
                                  

  RK_EXT(RK_lc,HEII,PRS,k,j,i) += coef[HEI][HEII]/(ATOM_MASS_AM(model,HEII)+ATOM_MASS_AM(model,HEI)) * COL_PRS_TRANS(HEII,HEI)+
                                  coef[HII][HEII]/(ATOM_MASS_AM(model,HEII)+ATOM_MASS_AM(model,HII))  * COL_PRS_TRANS(HEII,HII)+
                                  coef[HEII][ELCI]/(ATOM_MASS_AM(model,HEII)+ATOM_MASS_AM(model,ELCI)) * COL_PRS_TRANS(HEII,ELCI)
                                  ;

  RK_EXT(RK_lc,ELCI,0,k,j,i) += coef[HII][ELCI]/(ATOM_MASS_AM(model,ELCI)+ATOM_MASS_AM(model,HII))  * COL_PRS_TRANS(ELCI,HII)+
                                coef[HEII][ELCI]/(ATOM_MASS_AM(model,ELCI)+ATOM_MASS_AM(model,HEII)) * COL_PRS_TRANS(ELCI,HEII)
                                ;



//first_call = 0;
*/


/* ***** source term template ****** */
/*
     RK_EXT(RK_lc,HI,NDT,k,j,i) += 0;
 
    RK_EXT(RK_lc,HI,VX1,k,j,i) += 0;

    RK_EXT(RK_lc,HI,TMP,k,j,i) += 0;
*/



void get_xuv_flux(Model *model, Grid *grid)
{
/*
  int k,j,i,b;
  int nbeam;
  double *sed, *dx1;
  double **sigma;
  double ****flux;  
  static double SECTION_COST, scrh;
  static int first_call = 1;

  //initialize
  if(first_call){
    SECTION_COST = UNIT_DENSITY*UNIT_LENGTH/CONST_amu;
  }

  dx1 = grid[IDIR].dx;
  flux = model->xuv.flux;
  nbeam = model->xuv.nbeam;
  sigma = model->xuv.sigma;
  
  //get flux in every cell
  for(b=0;b<nbeam;b++)
    for(i = IEND;i>=IBEG;i--){
      scrh =   MODEL_EXT(model,HI,RHO,0,0,i)* dx1[i]*UNIT_LENGTH *sigma[HI][b];
      scrh +=  MODEL_EXT(model,HEI,RHO,0,0,i)* dx1[i]*UNIT_LENGTH *sigma[HEI][b];

      flux[b][0][0][i] = flux[b][0][0][i+1]*exp(-scrh);
    }
  
  first_call = 0;
  return;
*/
}

void quasi_neutrality(Model *model)
{
/*
  int k,j,i;
  Fluid *fluid;
  fluid = model->fluid;

  DOM_LOOP(k,j,i){
    FLUID_EXT(fluid,ELCI,NDT,k,j,i) =  FLUID_EXT(fluid,HII,NDT,k,j,i) + FLUID_EXT(fluid,HEII,NDT,k,j,i);//fully is NDT*charge
    FLUID_EXT(fluid,ELCI,VX1,k,j,i) = (FLUID_EXT(fluid,HII,NDT,k,j,i)*FLUID_EXT(fluid,HII,VX1,k,j,i)+
                                       FLUID_EXT(fluid,HEII,NDT,k,j,i)*FLUID_EXT(fluid,HEII,VX1,k,j,i))/
                                       FLUID_EXT(fluid,ELCI,NDT,k,j,i);
  }
  return;
*/
}




















