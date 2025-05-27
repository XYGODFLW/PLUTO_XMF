#include"matrix_tools.h"
#include"state_globals.h"
//creat a state machine for matrix solving

void lu_solver(double **mtx_A, double *b)
/* ************************************************
 * LU Decomposition method solving y for equation
 * 
 *  a * y = b 
 *
 * the result stored in b 
 * meanwhile discard the b's orignal data
 *
 * ************************************************ */
{
    int i, j, k;
    int mtx_size = g_mtx_size;

    static double tiny = 1.e-18;

    MTX_LOOP(i,j){
        if(i>j){  
            ARRAY_LOOP(k){
                if(k<j)
                    mtx_A[i][j] -= mtx_A[i][k]*mtx_A[k][j];//get L matrix
            }
            if(fabs(mtx_A[j][j])<tiny)printf("Over TINY in L_U_ decomposition\n"),mtx_A[j][j] = tiny;
            mtx_A[i][j] /= mtx_A[j][j];
        }
        else{
              ARRAY_LOOP(k)
                if(k<i)
                    mtx_A[i][j] -= mtx_A[i][k]*mtx_A[k][j];//get U matrix
        }
    }   

    //substitution and backsubstitution
    for(i=0;i<mtx_size;i++)       
        for(j=0;j<i;j++)
            b[i] -= mtx_A[i][j]*b[j];

    
    for(i = mtx_size - 1; i>=0 ;i--){
        for(j = mtx_size - 1; j>i; j--){
            b[i] -= mtx_A[i][j]*b[j];
        }
        b[i] = b[i]/mtx_A[i][i];
    }
}

void show_matrix(double **a, char* type)
{
   int k,j,i;
   NX2_LOOP(j){
     NX1_LOOP(i)printf(type ,a[j][i]);
     putchar('\n');
   }
   putchar('\n');
}

void mtx_clean(double **a)
{
  int k,j,i;
  MTX_LOOP(j,i){
    if(j==i)
      a[j][i] = 1;
    else a[j][i] = 0;
  }
}





















