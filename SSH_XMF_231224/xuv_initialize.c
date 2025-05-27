#include "pluto.h"
static char next_line(FILE*);
static void read_int0D(FILE*, char*, int*);
static void read_dbl0D(FILE*, char*, double*);
static void read_dbl1D(FILE*, char*, int, double*);
static void read_dbl2D(FILE*, char*, int, int, double**);

void get_sigma(double **, int , char *);


void xuv_initialize(Xuv *xuv, char* sed_file, char* sigma_file)
{
/*
  int k,j,i,b,mynode;
  int nbeam;
  double flux_504A = 0;
  double **mtx;
  FILE *fp;

  if((fp = fopen(sed_file, "r"))==NULL){printf("ERROR: !SED_FILE:%s do not exist, please check your setting\n",sed_file); QUIT_PLUTO(1);}

  read_int0D(fp, "nbeam", &(xuv->nbeam));
  nbeam = xuv->nbeam;

  sprintf(xuv->SED_file,"%s",sed_file);  

  //xuv struct initialize
  mtx = ARRAY_2D(nbeam, 2, double);
  xuv->wavelength = ARRAY_1D(nbeam, double);
  xuv->sed = ARRAY_1D(nbeam, double);  
  xuv->hmu = ARRAY_1D(nbeam, double);
  xuv->sigma = ARRAY_2D(2, nbeam, double);//2: equal to number of neutral fluid

  xuv->flux = ARRAY_4D(nbeam, NX3_TOT, NX2_TOT, NX1_TOT, double);
  xuv->flux_tot = 0;

  //load file
  read_dbl2D(fp, "flux", nbeam, 2, mtx);
  BEAM_LOOP(b,nbeam){
    xuv->wavelength[b] = mtx[b][0];
    xuv->sed[b] = mtx[b][1];
    xuv->hmu[b] = CONST_h*CONST_c/xuv->wavelength[b]*1.e8;
    xuv->flux_tot += xuv->sed[b];
    if(xuv->wavelength[b]<504)flux_504A += xuv->sed[b];
  }
  xuv->spectral_index_504A = flux_504A/xuv->flux_tot;
  get_sigma(xuv->sigma, nbeam, sigma_file);

  #ifdef PARALLEL
  MPI_Comm_rank(MPI_COMM_WORLD, &mynode);
  if(mynode==0){
    printf("> XUV initialize\n");
    printf("  nbeam = %d\n", nbeam);
    printf("  flux_tot(1AU) = %.2f[erg/cm2/s]\n", xuv->flux_tot);
    printf("  spectral_index_504 = %.2f\n", xuv->spectral_index_504A);
  }  
  #endif

if(g_inputParam[XUV_RESIZE]){
  //resize SED based on file "pluto.ini"
  BEAM_LOOP(b,nbeam){
    if(xuv->wavelength[b]<504){
      xuv->sed[b] *= g_inputParam[XUV_TOT_FLUX]/xuv->flux_tot * g_inputParam[SPECTRAL_INDEX_504A]/xuv->spectral_index_504A;
    }else{
      xuv->sed[b] *= g_inputParam[XUV_TOT_FLUX]/xuv->flux_tot * (1-g_inputParam[SPECTRAL_INDEX_504A])/(1-xuv->spectral_index_504A);
    }
  }
}
  //resize SED in planet position

  xuv->flux_tot = 0;
  flux_504A = 0;
  BEAM_LOOP(b,nbeam){
    xuv->sed[b] *= 1/g_inputParam[R_Orbit]/g_inputParam[R_Orbit]*g_inputParam[XUV_FACTOR];
    xuv->flux_tot += xuv->sed[b];
    if(xuv->wavelength[b]<504) flux_504A += xuv->sed[b];
  }


  #ifdef PARALLEL
  MPI_Comm_rank(MPI_COMM_WORLD, &mynode);
  if(mynode==0){
    printf("> XUV in planet position\n");
    printf("  nbeam = %d\n", nbeam);
    printf("  flux_tot(1AU) = %.2f[erg/cm2/s]\n", 
           xuv->flux_tot/g_inputParam[XUV_FACTOR]*pow(g_inputParam[R_Orbit],2));
    printf("  flux_tot(%.4fAU) = %.2f[erg/cm2/s]\n",
           g_inputParam[R_Orbit], xuv->flux_tot/g_inputParam[XUV_FACTOR]);
    printf("  spectral_index_504 = %.2f\n", flux_504A/xuv->flux_tot);
  }  
  #endif



  FreeArray2D(mtx);
  fclose(fp);
  return;
*/
}


void get_sigma(double **sigma, int nbeam, char *file_name)
{
/*
  FILE *fp; 
 
  if((fp = fopen(file_name,"r")) == NULL)printf("readfile %s is fail\n", file_name);

  read_dbl1D(fp, "sigma_H",  nbeam, sigma[HI]);

  read_dbl1D(fp, "sigma_He", nbeam, sigma[HEI]);

  fclose(fp);
*/
  return;
}





static void read_int0D(FILE* fp, char* term, int* output)
{
  char s[128];  

  rewind(fp);
  while(1){
    fscanf(fp,"%s",s);
    if(!strcmp(s,term)){
      fscanf(fp, "%d", output);
      return;
    }
    if(next_line(fp)){printf("> error: term %s is NG\n", term); return;}
  }
}

static void read_dbl1D(FILE* fp, char* term, int n, double* output)
{
  static char s[128];  
  int i;

  rewind(fp);
  while(1){
    fscanf(fp,"%s",s);
    if(!strcmp(s,term)){
      for(i=0;i<n;i++) fscanf(fp, "%lf,", output+i);
      return;
    }else if(feof(fp)){
      printf("> error: term %s is NG\n", term);exit(1); 
      return;
    }
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



static char next_line(FILE* fp)
{
  while( (fgetc(fp)!='\n') && (!feof(fp)) );
  return feof(fp);    
}





























