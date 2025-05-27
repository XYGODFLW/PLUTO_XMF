/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Set/retrieve output data attributes.
  
  The function SetOutput() sets, for each output data type (DBL, FLT, 
  VTK etc..) the default attributes of the corresponding ::Output structures.
  These include the variable name, a pointer to the actual 
  3D array, the centering of the variable (center/staggered), a 
  conditional inclusion flag (telling if the corresponding variable has
  to be written in the specified format), and so on.
  
  The function SetDumpVar() can be used to include or exclude a given 
  variable to be written using a particular output format.
  
  The function GetUserVar() returns the memory address to a 
  user-defined 3D array.

  \note Starting with PLUTO 4.1 velocity and magnetic field components 
        will be saved as scalars when writing VTK output. 
        If this is not what you want and prefer to save them as vector 
        fields (VTK VECTOR attribute), set VTK_VECTOR_DUMP to YES
        in your definitions.h.        
  
  \authors A. Mignone (mignone@ph.unito.it)
  \date    Aug 24, 2015
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

#ifndef VTK_VECTOR_DUMP  
 #define VTK_VECTOR_DUMP NO
#endif

static Output *all_outputs;
/* ********************************************************************* */
void hd_SetOutput (Data *d, Runtime *runtime)
/*!
 *  Set default attributes (variable names, pointers to data structures, 
 *  filename extensions, etc...) of the output structures.
 *
 * \param [in] d        pointer to Data structure
 * \param [in] runtime  pointer to Runtime structure
 *
 *********************************************************************** */
{
  int nv, i, k;
  Output *output;
  
  if (runtime->user_var > 0)
    d->Vuser = ARRAY_4D(runtime->user_var, NX3_TOT, NX2_TOT, NX1_TOT, double);
  else
    d->Vuser = NULL;

  //printf(">>>>>>>>>>> %d\n", d->Vuser);//debug
  all_outputs = runtime->output;

/* ---------------------------------------------
          Loop on output types 
   --------------------------------------------- */

  for (k = 0; k < MAX_OUTPUT_TYPES; k++){ 
    output = runtime->output + k;
    output->var_name = ARRAY_2D(64,128,char);
    output->stag_var = ARRAY_1D(64, int);
    output->dump_var = ARRAY_1D(64, int);
    strcpy(output->dir, runtime->output_dir); /* output directory is the same     */
                                            /* for all outputs (easy to change) */
    output->nfile    = -1;  

  /* -- set variables names -- */

    hd_SetDefaultVarNames(output);

  /* -- Set array pointers -- */

    for (nv = 0; nv < MHD_NVAR; nv++){
      output->V[nv]        = d->Vc[nv];
      output->stag_var[nv] = -1; /* -- means cell centered -- */ 
    }
    nv = HD_NVAR;

    output->nvar = nv;

  /* -- repeat for user defined vars -- */

    for (i = 0; i < runtime->user_var; i++){
      sprintf (output->var_name[i + nv], "%s", runtime->user_var_name[i]);
      output->V[i + nv] = d->Vuser[i];
      output->stag_var[i + nv] = -1; /* -- assume cell-centered -- */
    }

  /* -- add user vars to total number of variables -- */

    output->nvar += runtime->user_var;
 
  /* -- select which variables are going to be dumped to disk  -- */

    for (nv = output->nvar; nv--; ) output->dump_var[nv] = YES;
    #if ENTROPY_SWITCH
     output->dump_var[ENTR] = NO;
    #endif

    switch (output->type){
      case DBL_OUTPUT:   /* -- dump ALL variables -- */
        sprintf (output->ext,"dbl");
        break;
      case FLT_OUTPUT:   /* -- do not dump staggered fields (below)-- */
        sprintf (output->ext,"flt");
        break;
      case DBL_H5_OUTPUT:   /* -- dump ALL variables -- */
        sprintf (output->ext,"dbl.h5");
        break;
      case FLT_H5_OUTPUT:   /* -- do not dump staggered fields (below)-- */
        sprintf (output->ext,"flt.h5");
        break;
      case VTK_OUTPUT:   /* -- do not dump staggered fields (below) -- */
        sprintf (output->ext,"vtk");
        #if VTK_VECTOR_DUMP == YES
         D_EXPAND(output->dump_var[VX1] = VTK_VECTOR;  ,
                  output->dump_var[VX2] = NO;          ,
                  output->dump_var[VX3] = NO;)
        #endif
        break;
      case TAB_OUTPUT:   /* -- do not dump staggered fields -- */
        sprintf (output->ext,"tab");
        break;
      case PPM_OUTPUT:   /* -- dump density only  -- */
        sprintf (output->ext,"ppm");
        for (nv = output->nvar; nv--; ) output->dump_var[nv] = NO;
        break;
      case PNG_OUTPUT:   /* -- dump density only  -- */
        sprintf (output->ext,"png");
        for (nv = output->nvar; nv--; ) output->dump_var[nv] = NO;
        break;
    }
    
  /* ---------------------------------------------------------------
      for divergence cleaning never dump the scalar psi unless
      the output type can be potentially used for restart
     --------------------------------------------------------------- */
   
  }
 
  ChangeDumpVar();
}

