#include"pluto.h"


void mmd_write_save(RK_var* u, char* file_name)
{
  int f,k,j,i,nv;
  FILE *fp = fopen(file_name, "w");

    FMHD_LOOP(f){
      MHD_NVAR_LOOP(nv)KDOM_LOOP(k)JDOM_LOOP(j) {
        fwrite(u[f].Vc[nv][k][j]+IBEG, sizeof(double), (IEND-IBEG+1), fp);
      }
      #ifdef STAGGERED_MHD
      DIM_LOOP(nv)KDOM_LOOP(k)JDOM_LOOP(j) {
        fwrite(u[f].Vs[nv][k][j]+IBEG, sizeof(double), (IEND-IBEG+1), fp);
      }
      #endif
    }
 
    FHD_LOOP(f){
      HD_NVAR_LOOP(nv)KDOM_LOOP(k)JDOM_LOOP(j) {
        fwrite(u[f].Vc[nv][k][j]+IBEG, sizeof(double), (IEND-IBEG+1), fp);
      }
    }

  return;
}

void mmd_load_save(char *file_name, RK_var* u)
{
  int f,k,j,i,nv;
  FILE *fp = fopen(file_name, "r");

    FMHD_LOOP(f){
      MHD_NVAR_LOOP(nv)KDOM_LOOP(k)JDOM_LOOP(j) {
        fread(u[f].Vc[nv][k][j]+IBEG, sizeof(double), (IEND-IBEG+1), fp);
      }

      #ifdef STAGGERED_MHD
      DIM_LOOP(nv)KDOM_LOOP(k)JDOM_LOOP(j) {
        fread(u[f].Vs[nv][k][j]+IBEG, sizeof(double), (IEND-IBEG+1), fp);
      }
      #endif
    }
 
    FHD_LOOP(f){
      HD_NVAR_LOOP(nv)KDOM_LOOP(k)JDOM_LOOP(j) {
        fread(u[f].Vc[nv][k][j]+IBEG, sizeof(double), (IEND-IBEG+1), fp);
      }
    }

  return;
}

























