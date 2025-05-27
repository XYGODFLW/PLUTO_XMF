#include "pluto.h"

/* *************************************************************** */
void ComputeUserVar (const Data *d, Grid *grid)
/*
 *
 *  PURPOSE
 *
 *    Define user-defined output variables
 *
 *
 *
 ***************************************************************** */
{
  int k, j, i;  
  double *x1,*x2,*x3;
  double ***r, ***theta;
    
  x1 = grid[IDIR].x;  
  x2 = grid[JDIR].x;
  x3 = grid[KDIR].x;  

  r     = d->Vuser[0];
  theta = d->Vuser[1];

  DOM_LOOP(k,j,i){
    r[k][j][i]     = x1[i];
    theta[k][j][i] = x2[j];
  }

 
  return;
}
/* ************************************************************* */
void ChangeDumpVar ()
/* 
 *
 * 
 *************************************************************** */
{ 
  Image *image;

}





