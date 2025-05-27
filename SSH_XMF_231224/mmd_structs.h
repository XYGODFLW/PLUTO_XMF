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

  Fluid *fluid;
  double ****E_fluid;   
  double ****B_field;
  double ****Bs;

  Xuv xuv;

}Model;

#ifdef PARALLEL
typedef struct NODEDATA{ //collect data from all processors. Caution: this struct only used in core0
  double ****data;       // RK_data model
  double *ck,*cj,*ci;    // center position of grid in x,y,z director
  int bk,bj,bi, ek,ej,ei, lek,lej,lei; // Global beg and end of this node
  int ltk,ltj,lti;
}NodeData;

typedef struct GLOBALDATA{ //reformed data. Caution: this struct only used in core0
  double ****data;    // RK_data
  double *ck, *cj, *ci; // center position of grid in x,y,z director
}GlobalData;


#endif



























