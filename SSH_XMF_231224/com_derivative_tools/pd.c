#include"pluto.h"
#include"pd.h"
#include"pd_globals.h"
double pd_1D(double f(double*), double* x, int var)
/* **********************
   \param [in] f   function
   \param [in] x   var list
   \param [in] var derivation of var
   **********************/
{
  int i,j,k;
  static double tiny = 1.e-6;
  double tiny2 = tiny/2;
  static double x0[32], x1[32],dx;
  if(var >= PD_VAR_END){printf("ERROR: pd var over range\n");return -1;}

  PD_LOOP(i){
    x0[i]=x[i],x1[i]=x[i];
  }

  if(fabs(x[var])<tiny){
    x0[var] -= tiny*tiny;
    x1[var] += tiny*tiny;
    dx = 2*tiny*tiny;
  }else{
    x0[var] *= (1-tiny);
    x1[var] *= (1+tiny);
    dx = x[var]*tiny2;
  }
  return (f(x1)-f(x0)) / dx;
} 

double pd_2D(double f(double**), double** x, int var0, int var1)
{
  int i,j,k;
  static double tiny = 1.e-4;
  double tiny2 = 2*tiny;
  static double **x0, **x1, dx;
  static int first_call = 1;

  if(first_call){
    x0 = ARRAY_2D(32,32,double);
    x1 = ARRAY_2D(32,32,double);
  }

  if(var0 >= PD_VAR0_END && var1 >= PD_VAR1_END){printf("ERROR: pd var over range\n");return -1;}

  PD_TOT_LOOP(j,i){
    x0[j][i]=x[j][i], x1[j][i]=x[j][i];
  }

  if(fabs(x[var0][var1])<tiny*tiny){
    x0[var0][var1] -= tiny*tiny;
    x1[var0][var1] += tiny*tiny;
    dx = 2*tiny*tiny;
  }else{
    x0[var0][var1] *= (1-tiny);
    x1[var0][var1] *= (1+tiny);
    dx = x[var0][var1]*tiny2;
  }

  first_call = 0;
  return (f(x1)-f(x0)) / dx;
}
 
void pd_setting(int var_end){PD_VAR_END=var_end;}
void pd_2D_setting(int var0_end, int var1_end){PD_VAR0_END=var0_end, PD_VAR1_END=var1_end;}

















