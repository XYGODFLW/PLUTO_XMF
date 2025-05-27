#ifndef MATRIX_TOOLS_H
#define MATRIX_TOOLS_H

#include<math.h>
#include<stdio.h>
#include"matrix_prototype.h"

extern int g_mtx_size;
extern int g_mtx_nx1;
extern int g_mtx_nx2;

#include"matrix_macro.h"

#define MTX_COL_CAL(mtx,column,operator,array) for(int nv_=0;nv_<g_mtx_nx2;nv_++)(mtx[nv_][column])operator(array[nv_])
#define MTX_COL_PARA(mtx,column,operator,parameter) for(int nv_=0;nv_<g_mtx_nx2;nv_++) mtx[nv_][column]operator(parameter)
#define COL_CAL(array0,operator,array1) for(int nv_=0;nv_<g_mtx_nx1;nv_++) ((array0)[nv_])operator((array1)[nv_]) 
#define COL_PARA(array0,operator,parameter) for(int nv_=0;nv_<g_mtx_nx1;nv_++) ((array0)[nv_])operator(parameter) 

#endif
