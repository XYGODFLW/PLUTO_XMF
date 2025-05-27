#include "pluto.h"
static void ptc_ini(Partical *ptc, char*, int, int, int, double);
static void show_configure(Model *model);

void mmd_model_initialize(Model *model, Runtime *ini){
  int k,j,i,nv,f;
  int nf;
  Fluid *fluid;
//pointer lock
  model->nfluid = nf = 5;
  model->fluid = ARRAY_1D(nf, Fluid);

//reading file get settings...(set by handy right now) 

  fluid = model->fluid;

  //allocate memoray for every fluid
  for(f=0;f<nf;f++){
    fluid[f].f = f;
    fluid[f].Vc = ARRAY_4D(3, NX3_TOT, NX2_TOT, NX1_TOT, double);
  } 

  //setting 
  //sequence: label element charge type mass q_m
  ptc_ini(&((fluid+HI)->ptc),	"H",	H,	0,	NEU,	1.00794);//
  ptc_ini(&((fluid+HEI)->ptc),	"HE",	HE,	0,	NEU,	4.002602);//
  ptc_ini(&((fluid+HII)->ptc),	"H+",	H,	1,	ION,	1.00739);//
  ptc_ini(&((fluid+HEII)->ptc),	"HE+",	HE,	1,	ION,	4.002053);//
  ptc_ini(&((fluid+ELCI)->ptc), "ELC",	ELECTRON,-1,	ELC,	5.4858e-4);//  

  //XUV initialize 
  xuv_initialize(&(model->xuv), ini->SED_file, "SED_lists/photoIon_sigma.dat");


  //electron collision ionization  


  //collision load

  //show_setting()
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
  printf("\n> MODEL SETTING LIST\n\n");
  printf("  NUMBER OF FLUIDS = %d\n",nf);

  for(f=0;f<nf;f++){
    ptc = &(model->fluid[f].ptc);
    printf("  fluid: %-12dlabel: %-12scharge: %-12dMass: %-12.6fQ/M: %-12.6f\n",
            f,ptc->label,ptc->charge,ptc->mass,ptc->q_m);   
  }  
  printf("\n");

  //show xuv setting
  printf("\n> XUV SETTING LIST\n");
  printf("> SED_FILE: %s\n", model->xuv.SED_file);
  printf("%12s%12s%12s%12s%12s\n", "wavelength","hmu","SED","sigma_H","sigma_He");
  for(b=0;b<model->xuv.nbeam;b++){
    printf("%12.4f%12.4e%12.4e%12.4e%12.4e\n", 
           model->xuv.wavelength[b], model->xuv.hmu[b], model->xuv.sed[b], model->xuv.sigma[HI][b], model->xuv.sigma[HEI][b]);  
  }
  printf("\n> SETTING END \n\n");

  #ifdef PARALLEL
  }
  #endif
  
  return;
}













































