//mmd 
typedef struct INTRODUCTION{
  int nline;
  char* lines[];
}Intro;

/* Fluid related structs */
typedef struct PARTICAL{
  char label[128];  
  int element;
  int charge;
  int type;
  double mass;////UNIT: relative atomic mass
  double q_m;//charge-to-mass-ratio

}Partical;


typedef 
struct RK_VAR{//MHD,HD will not be seperated
  double ****Vc;
}RK_var;



typedef struct FLUID{
  int f;//number of this fluid
  Partical ptc;
  double ****Vc;

}Fluid;



typedef struct XUV{
  Intro information;
  char SED_file[128];

  int nbeam;
  double *wavelength;//in unit A
  double *hmu;
  double *sed;
  double **sigma;

  double flux_tot;
  double spectral_index_504A;

  double ****flux;//flux[nbeam][k][j][i]
}Xuv;

typedef struct MODEL{
  Intro information;
  
  int nfluid;
  //int IF_ION,IF_ELC;
  //int fNEU[128],fION[128],fELC[128];

  Fluid *fluid;    

  Xuv xuv; 

}Model;
























