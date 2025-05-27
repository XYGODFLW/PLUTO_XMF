/*
  Only for number data csv
  csv file format
  (var_name:) ..,..,..,
  (datas:)    ..,..,..,

              ......

              ..,..,..,
  by Xing,L 2023.9 (xinglei@ynao.ac.cn)
*/

#include<stdio.h>
#include<stdlib.h>

double line_ITP(double **data, int lnmax, double xout)
{
  const int LN=0, DT=1;//LN:line number, DT: data
  int i, iL=0, iR=lnmax-1;

  for(i=0;i<lnmax;i++){
    if(xout<data[i][LN]){
      iR=i;break;
    }else{
      iL=i;
    }
  }

  if(iR==iL){
    return data[iR][DT];
  }else{
    return ( (data[iR][LN]-xout)*data[iL][DT] + (xout-data[iL][LN])*data[iR][DT] )
           /(data[iR][LN]-data[iL][LN]);
  }  
}





















