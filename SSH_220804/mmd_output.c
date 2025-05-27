#include"pluto.h"
#define TOT_SIZE HD_NVAR*NX1_TOT
#define XUV_SIZE nbeam*NX1_TOT

#define VAR_INDEX(nv,k,j,i) ((( (nv)*NX3_TOT+(k))*NX2_TOT+(j))*NX1_TOT+(i) ) 
#define N_IDIR_LOOP(node,i) for(i = ibeg_list[node];i<=iend_list[node];i++)
//#define XUV_SIZE 

void mmd_model_output(Model *model, Grid *grid, Runtime *ini, double t_step, char* folder)
{
  int k,j,i,nv,f,b,node;
  Xuv *xuv;
  static int step = -1, nstep;
  double scrh = 0;
  char file_name[128];

  Fluid *fluid;

  static int nfluid, nbeam, first_call=1;
  static int nlabel = 0;
  static char label[128][32];
  FILE *fp;

  nstep = (int)(g_time/t_step);

  if(first_call){
    nfluid = model->nfluid;
    nbeam = model->xuv.nbeam;
    sprintf(label[0],"I"), nlabel++;
    sprintf(label[1],"X1"), nlabel++;

    for(f=0;f<nfluid;f++){
      nlabel += 3;
      sprintf(label[2+3*f], "%s_%s", "NDT", model->fluid[f].ptc.label);
      sprintf(label[3+3*f], "%s_%s", "VX1", model->fluid[f].ptc.label);
      sprintf(label[4+3*f], "%s_%s", "TMP", model->fluid[f].ptc.label);
    }
  }


  #ifdef PARALLEL

  int mynode; 
  int totalnodes;
  static int *ibeg_list, *iend_list, iscrh;
  static double ***fluid_recv, **xuv_recv, **grid_recv[3];
  MPI_Status status;
  
  MPI_Comm_size(MPI_COMM_WORLD, &totalnodes);
  MPI_Comm_rank(MPI_COMM_WORLD, &mynode);

  if(mynode==0){
    if(fluid_recv==NULL) fluid_recv = ARRAY_3D(model->nfluid, totalnodes, TOT_SIZE, double);
    if(grid_recv[IDIR]==NULL) grid_recv[IDIR] = ARRAY_2D(totalnodes, NX1_TOT, double);
  //if(grid_recv[JDIR]==NULL) grid_recv[JDIR] = ARRAY_2D(totalnodes, NX2_TOT, double);
  //if(grid_recv[KDIR]==NULL) grid_recv[KDIR] = ARRAY_2D(totalnodes, NX3_TOT, double);
    if(xuv_recv==NULL) xuv_recv = ARRAY_2D(totalnodes, XUV_SIZE, double); 
    if(ibeg_list==NULL) ibeg_list = ARRAY_1D(totalnodes, int);
    if(iend_list==NULL) iend_list = ARRAY_1D(totalnodes, int);
  }

  //collect data from all processor
  //fluid
  for(f=0; f<model->nfluid; f++){
    MPI_Barrier (MPI_COMM_WORLD);
  
    if(mynode==0){ 
      for(node=1;node<totalnodes;node++){
        MPI_Recv(fluid_recv[f][node], TOT_SIZE, MPI_DOUBLE, node, 1, MPI_COMM_WORLD, &status);
      }
      memcpy(fluid_recv[f][0], model->fluid[f].Vc[0][0][0], TOT_SIZE*sizeof(double));

    }else{
      MPI_Send(model->fluid[f].Vc[0][0][0], TOT_SIZE, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
    }
  }


  //ibeg_list, iend_list  
  MPI_Barrier(MPI_COMM_WORLD);
  if(mynode==0){
    for(node=1;node<totalnodes;node++){
      MPI_Recv(ibeg_list+node, 1, MPI_INT, node, 1, MPI_COMM_WORLD, &status);
      MPI_Recv(iend_list+node, 1, MPI_INT, node, 2, MPI_COMM_WORLD, &status);
    }
 
  iscrh = IBEG, memcpy(ibeg_list, &iscrh, sizeof(int));
  iscrh = IEND, memcpy(iend_list, &iscrh, sizeof(int));
  }else{
    iscrh = IBEG, MPI_Send(&iscrh, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);
    iscrh = IEND, MPI_Send(&iscrh, 1, MPI_INT, 0, 2, MPI_COMM_WORLD);
  }
  

  //xuv
  MPI_Barrier (MPI_COMM_WORLD);
  if(mynode==0){ 
    for(node=1;node<totalnodes;node++){
      MPI_Recv(xuv_recv[node], XUV_SIZE, MPI_DOUBLE, node, 1, MPI_COMM_WORLD, &status);
    }
    memcpy(xuv_recv[0], model->xuv.flux[0][0][0], XUV_SIZE*sizeof(double));

  }else{
    MPI_Send(model->xuv.flux[0][0][0], XUV_SIZE, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
  }

  //grid
  MPI_Barrier (MPI_COMM_WORLD);
  if(mynode==0){ 
    for(node=1;node<totalnodes;node++){
      MPI_Recv(grid_recv[IDIR][node], NX1_TOT, MPI_DOUBLE, node, 1, MPI_COMM_WORLD, &status);
    }
    memcpy(grid_recv[IDIR][0], grid[IDIR].x, NX1_TOT*sizeof(double));

  }else{
    MPI_Send(grid[IDIR].x, NX1_TOT, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
  }



  //write data
  if(mynode==0)if(nstep != step){
    step = nstep;
    //write data file
    if(g_inputParam[OUTO_NAME]) 
      sprintf(file_name,"%s/%s_fxuv%.1f_Index%.3f_He%.4f_Heat%.3f_%.4d_data.dat", 
              ini->output_dir, ini->obj_name, model->xuv.flux_tot/g_inputParam[XUV_FACTOR], 
              model->xuv.spectral_index_504A, g_inputParam[HE_RATIO],
              g_inputParam[HEATING_EFFICIENCY], step);
    else sprintf(file_name,"%s/%.4d_data.dat", folder, step);   

    printf("> %s out\n", file_name);
    
    fp = fopen(file_name, "w+");
    for(i=0;i<nlabel;i++) fprintf(fp, "%-12s", label[i]); 
    fprintf(fp, "%s", "\n");

    for(node=0; node<totalnodes; node++) N_IDIR_LOOP(node,i){
      fprintf(fp, "%-12d", node*(IEND-IBEG+1)+i);
      fprintf(fp, "%-12.4e", grid_recv[IDIR][node][i]*UNIT_LENGTH);

      for(f=0;f<nfluid;f++){
        fprintf(fp, "%-12.4e", fluid_recv[f][node][NDT*NX1_TOT+i]);//VAR_INDEX(NDT,k,j,i)
        fprintf(fp, "%-12.4e", fluid_recv[f][node][VX1*NX1_TOT+i]);
        fprintf(fp, "%-12.4e", fluid_recv[f][node][PRS*NX1_TOT+i]/fluid_recv[f][node][NDT*NX1_TOT+i]/CONST_kB);
      }
      fprintf(fp, "%s", "\n");
    }



    //write xuv file
    xuv = &(model->xuv);
       
    sprintf(file_name, "xuv.dat");
    
    printf("> %s out\n", file_name);

    fp = freopen(file_name, "w", fp);
   
    fprintf(fp, "%-12s", label[0]);
    fprintf(fp, "%-12s", label[1]);
    for(b=0;b<xuv->nbeam;b++){
      fprintf(fp, "%-12.2f", xuv->wavelength[b]);
    }
    fprintf(fp, "%-12.2s", "FLUX_TOT");
    fprintf(fp, "%s", "\n");

    for(node=0; node<totalnodes; node++) N_IDIR_LOOP(node,i){
      fprintf(fp, "%-12d", node*(IEND-IBEG+1)+i);
      fprintf(fp, "%-12.4e", grid_recv[IDIR][node][i]*UNIT_LENGTH);
        
      scrh = 0;
      for(b=0;b<xuv->nbeam;b++){
        fprintf(fp, "%-12.4e",  xuv_recv[node][VAR_INDEX(b,0,0,i)]);
        scrh += xuv_recv[node][VAR_INDEX(b,0,0,i)];
      }
      fprintf(fp, "%-12.4e", scrh);
      fprintf(fp, "%s", "\n");
    }

    fclose(fp);
/*
*/
  }
  #else
  
  if(nstep != step){
    step = nstep;
    //write data file
    sprintf(file_name,"%s/%.4d_data.dat", folder, step);//tmp way
    
    printf("> %s out\n", file_name);
    
    fp = fopen(file_name, "w+");
    for(i=0;i<nlabel;i++) fprintf(fp, "%-12s", label[i]); 
    fprintf(fp, "%s", "\n");

    DOM_LOOP(k,j,i){
      fprintf(fp, "%-12d", i);
      fprintf(fp, "%-12.4e", grid[IDIR].x[i]*UNIT_LENGTH);

      for(f=0;f<nfluid;f++){
        fprintf(fp, "%-12.4e", MODEL_EXT(model,f,NDT,k,j,i));
        fprintf(fp, "%-12.4e", MODEL_EXT(model,f,VX1,k,j,i));
        fprintf(fp, "%-12.4e", MODEL_EXT(model,f,TMP,k,j,i));
      }
      fprintf(fp, "%s", "\n");
    }

    //write xuv file
    xuv = &(model->xuv);
       

    sprintf(file_name, "%s/%.4d_xuv.dat",folder, step);
    
    printf("> %s out\n", file_name);

    fp = freopen(file_name, "w", fp);

    for(b=0;b<xuv->nbeam;b++){
      fprintf(fp, "%-12.2f", xuv->wavelength[b]);
    }
    fprintf(fp, "%-12.2s", "FLUX_TOT");
    fprintf(fp, "%s", "\n");

    DOM_LOOP(k,j,i){
      scrh = 0;
      for(b=0;b<xuv->nbeam;b++){
        fprintf(fp, "%-12.4e",  xuv->flux[b][k][j][i]);
        scrh += xuv->flux[b][k][j][i];
      }
      fprintf(fp, "%-12.4e", scrh);
      fprintf(fp, "%s", "\n");
    }

    fclose(fp);
  }
#endif
  
  first_call = 0;
  return;
}











































