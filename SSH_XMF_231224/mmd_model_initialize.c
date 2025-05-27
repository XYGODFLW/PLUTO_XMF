#include "pluto.h"
static void ptc_ini(Partical *ptc, char*, int, int, int, double);
static void show_configure(Model *model);

void mmd_model_initialize(Model *model, Runtime *ini){
  int k,j,i,nv,f;
  int nf;
  Fluid *fluid;
//pointer lock
  model->nfluid = nf = NFLD;
  model->fluid = ARRAY_1D(nf, Fluid);

//reading file get settings...(set by handy right now) 

  fluid = model->fluid;

  //allocate memoray for every fluid
  for(f=0;f<nf;f++){
    fluid[f].f = f;
    fluid[f].Vc = ARRAY_4D(HD_NVAR, NX3_TOT, NX2_TOT, NX1_TOT, double);
  } 
  #if ELC_ON
  model->E_fluid = ARRAY_4D(HD_NVAR, NX3_TOT, NX2_TOT, NX1_TOT, double);
  #endif

  #if MAG_ON //if STAGGRED_MHD ON, this variable ONLY for OUTPUT
  model->B_field = ARRAY_4D(COMPONENTS, NX3_TOT, NX2_TOT, NX1_TOT, double);
  #endif  


  #ifdef STAGGERED_MHD
  model->Bs = ARRAY_1D(COMPONENTS, double ***);
   D_EXPAND(
     model->Bs[BX1s] = ArrayBox( 0, NX3_TOT-1, 0, NX2_TOT-1,-1, NX1_TOT-1); ,
     model->Bs[BX2s] = ArrayBox( 0, NX3_TOT-1,-1, NX2_TOT-1, 0, NX1_TOT-1); ,
     model->Bs[BX3s] = ArrayBox(-1, NX3_TOT-1, 0, NX2_TOT-1, 0, NX1_TOT-1);)
  #endif
  
  //setting 
  //sequence: label element charge type mass q_m
  ptc_ini(&((fluid+HI)->ptc),		"HI",	H,	0,	NEU,	1.00794);//
  ptc_ini(&((fluid+HII)->ptc),	"HII",	H,	1,	ION,	0.503695);//0.503695if electron pressure is not considered,
                                                                               //the mass of ion will be shared with electron

  //XUV initialize 
  xuv_initialize(&(model->xuv), ini->SED_file, "SED_lists/photoIon_sigma.dat");

  //show_settings()
  show_configure(model);
  return;
}

static void ptc_ini(Partical *ptc, char *label, int element, int charge, int type, double mass)
{
  sprintf(ptc->label, label);
  ptc->element = element;
  ptc->charge = charge;
  ptc->type = type;
  ptc->mass = mass;
  ptc->q_m = charge/mass;

  return;
}

static void show_configure(Model *model)
{
  int f,k,j,i,b,mynode;
  int nf;
  Partical *ptc;

  nf = model->nfluid;
  
  #ifdef PARALLEL
  MPI_Comm_rank(MPI_COMM_WORLD, &mynode);
  if(mynode==0){
  #endif

  //show fluid setting 
  print1("\n> MODEL SETTING LIST\n\n");
  print1("  NUMBER OF FLUIDS = %d\n",nf);

  for(f=0;f<nf;f++){
    ptc = &(model->fluid[f].ptc);
    print1("  fluid: %-12dlabel: %-12scharge: %-12dMass: %-12.6fQ/M: %-12.6f\n",
            f,ptc->label,ptc->charge,ptc->mass,ptc->q_m);   
  }  
  print1("\n");

  //show xuv setting
  print1("\n> XUV SETTING LIST\n");
  print1("> SED_FILE: %s\n", model->xuv.SED_file);
  print1("%12s%12s%12s%12s%12s\n", "wavelength","hmu","SED","sigma_H","sigma_He");
  //for(b=0;b<model->xuv.nbeam;b++){
    //print1("%12.4f%12.4e%12.4e%12.4e%12.4e\n", 
  //         model->xuv.wavelength[b], model->xuv.hmu[b], model->xuv.sed[b], model->xuv.sigma[HI][b], model->xuv.sigma[HEI][b]);  
  //}
  print1("\n> SETTING END \n\n");

  #ifdef PARALLEL
  }
  #endif
  
  return;
}













































