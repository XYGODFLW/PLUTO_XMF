#include"pluto.h"
#define TOT_SIZE NX3_TOT*NX2_TOT*NX1_TOT
#define XUV_SIZE nbeam  *NX3_TOT*NX2_TOT*NX1_TOT

#define VAR_INDEX(nv,k,j,i) ((( (nv)*NX3_TOT+(k))*NX2_TOT+(j))*NX1_TOT+(i) )


//#define XUV_SIZE 
static void collect_data();

void mmd_model_output(Model *model, Grid *grid, Runtime *ini, double t_step, char* folder)
{

  int k,j,i,nv,f,b;
  Xuv *xuv;
  static int step = -1, nstep;
  double scrh = 0;
  char file_name[128];

  Fluid *fluid;

  static int nfluid, nbeam, first_call=1;
  static int nlabel = 0;
  static char labels[256][32];
  FILE *fp;

  nstep = (int)(g_time/t_step);
 

  //Assigning labels for every variable
  if(first_call){
    nfluid = model->nfluid;
    nbeam = model->xuv.nbeam;
    RD_EXPAND(
      sprintf(labels[nlabel++],"K");,
      sprintf(labels[nlabel++],"J");,
      sprintf(labels[nlabel++],"I");)

    RD_EXPAND(
      sprintf(labels[nlabel++],"X3");,
      sprintf(labels[nlabel++],"X2");,
      sprintf(labels[nlabel++],"X1");)//follow the state structrue

    for(f=0;f<nfluid;f++){
      sprintf(labels[nlabel++], "%s_%s", "NDT", model->fluid[f].ptc.label);
      EXPAND(
        sprintf(labels[nlabel++], "%s_%s", "VX1", model->fluid[f].ptc.label);,
        sprintf(labels[nlabel++], "%s_%s", "VX2", model->fluid[f].ptc.label);,
        sprintf(labels[nlabel++], "%s_%s", "VX3", model->fluid[f].ptc.label);)

      sprintf(labels[nlabel++], "%s_%s", "PRS", model->fluid[f].ptc.label);
    }

    #if ELC_ON
      sprintf(labels[nlabel++], "ELC_NDT");
      EXPAND(
        sprintf(labels[nlabel++], "ELC_VX1");,
        sprintf(labels[nlabel++], "ELC_VX2");,
        sprintf(labels[nlabel++], "ELC_VX3");)

      sprintf(labels[nlabel++], "ELC_PRS", model->fluid[f].ptc.label);
    #endif

    #if MAG_ON
      EXPAND(
        sprintf(labels[nlabel++], "BX1");,
        sprintf(labels[nlabel++], "BX2");,
        sprintf(labels[nlabel++], "BX3");)
    #endif
  }


  #ifdef PARALLEL
  int mynode, totalnodes, nd; 
  int li_tot, lj_tot, lk_tot;
  static int *ibeg_list, *iend_list, iscrh;
  static double ****RK_reform;
  static NodeData *node_data;
  static GlobalData g_data;
  MPI_Status status;
  
  MPI_Comm_size(MPI_COMM_WORLD, &totalnodes);
  MPI_Comm_rank(MPI_COMM_WORLD, &mynode);  
  
  if(RK_reform==NULL){
    RK_reform = ARRAY_4D(RK_TOT, NX3_TOT, NX2_TOT, NX1_TOT, double);
    g_data.data = ARRAY_4D(RK_TOT, grid[KDIR].np_int_glob, grid[JDIR].np_int_glob, grid[IDIR].np_int_glob, double);
    g_data.ck = ARRAY_1D(grid[KDIR].np_int_glob, double); 
    g_data.cj = ARRAY_1D(grid[JDIR].np_int_glob, double); 
    g_data.ci = ARRAY_1D(grid[IDIR].np_int_glob, double); 
  }


  if(first_call){//initialize node_data and GData
  if(mynode==0){
    node_data = ARRAY_1D(totalnodes, NodeData);

    //get g_beg and g_end of all nodes
    node_data[0].bi = grid[IDIR].beg, node_data[0].bj = grid[JDIR].beg, node_data[0].bk = grid[KDIR].beg;
    node_data[0].ei = grid[IDIR].end, node_data[0].ej = grid[JDIR].end, node_data[0].ek = grid[KDIR].end;

    for(nd=1;nd<totalnodes;nd++){
      MPI_Recv(&(node_data[nd].bi), 1, MPI_INT, nd, 1, MPI_COMM_WORLD, &status);
      MPI_Recv(&(node_data[nd].bj), 1, MPI_INT, nd, 1, MPI_COMM_WORLD, &status);
      MPI_Recv(&(node_data[nd].bk), 1, MPI_INT, nd, 1, MPI_COMM_WORLD, &status);
      MPI_Recv(&(node_data[nd].ei), 1, MPI_INT, nd, 1, MPI_COMM_WORLD, &status);
      MPI_Recv(&(node_data[nd].ej), 1, MPI_INT, nd, 1, MPI_COMM_WORLD, &status);
      MPI_Recv(&(node_data[nd].ek), 1, MPI_INT, nd, 1, MPI_COMM_WORLD, &status);
    }
  }else{
      MPI_Send(&(grid[IDIR].beg), 1, MPI_INT, 0, 1, MPI_COMM_WORLD);
      MPI_Send(&(grid[JDIR].beg), 1, MPI_INT, 0, 1, MPI_COMM_WORLD);
      MPI_Send(&(grid[KDIR].beg), 1, MPI_INT, 0, 1, MPI_COMM_WORLD);
      MPI_Send(&(grid[IDIR].end), 1, MPI_INT, 0, 1, MPI_COMM_WORLD);
      MPI_Send(&(grid[JDIR].end), 1, MPI_INT, 0, 1, MPI_COMM_WORLD);
      MPI_Send(&(grid[KDIR].end), 1, MPI_INT, 0, 1, MPI_COMM_WORLD);
  }



  //allocate memory
  if(mynode==0)for(nd = 0;nd<totalnodes;nd++){
    node_data[nd].lek = KBEG + (node_data[nd].ek-node_data[nd].bk);
    node_data[nd].lej = JBEG + (node_data[nd].ej-node_data[nd].bj);
    node_data[nd].lei = IBEG + (node_data[nd].ei-node_data[nd].bi);

    node_data[nd].ltk = 2*grid[KDIR].nghost + 1 + (node_data[nd].ek-node_data[nd].bk);
    node_data[nd].ltj = 2*grid[JDIR].nghost + 1 + (node_data[nd].ej-node_data[nd].bj);
    node_data[nd].lti = 2*grid[IDIR].nghost + 1 + (node_data[nd].ei-node_data[nd].bi);

    node_data[nd].data = ARRAY_4D(RK_TOT, node_data[nd].ltk, node_data[nd].ltj, node_data[nd].lti, double);
    node_data[nd].ck = ARRAY_1D(node_data[nd].ltk, double);
    node_data[nd].cj = ARRAY_1D(node_data[nd].ltj, double);
    node_data[nd].ci = ARRAY_1D(node_data[nd].lti, double);
  }


  
  //collect grid infromation from all nodes
  if(mynode==0){ 
    for(nd=1;nd<totalnodes;nd++){
      MPI_Recv(node_data[nd].ck, node_data[nd].ltk, MPI_DOUBLE, nd, 1, MPI_COMM_WORLD, &status);
      MPI_Recv(node_data[nd].cj, node_data[nd].ltj, MPI_DOUBLE, nd, 1, MPI_COMM_WORLD, &status);
      MPI_Recv(node_data[nd].ci, node_data[nd].lti, MPI_DOUBLE, nd, 1, MPI_COMM_WORLD, &status);
    }
      memcpy(node_data[0].ck, grid[KDIR].x, NX3_TOT*sizeof(double));
      memcpy(node_data[0].cj, grid[JDIR].x, NX2_TOT*sizeof(double));
      memcpy(node_data[0].ci, grid[IDIR].x, NX1_TOT*sizeof(double));  
  }else{
      MPI_Send(grid[KDIR].x, NX3_TOT, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
      MPI_Send(grid[JDIR].x, NX2_TOT, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
      MPI_Send(grid[IDIR].x, NX1_TOT, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
  }

  }//ini node_data end

if(nstep != step){//check step
   step = nstep;

  //collect data from all processors
  if(mynode==0){ 
    for(nd=1;nd<totalnodes;nd++){//collect atmosphere
      MPI_Recv(node_data[nd].data[0][0][0], RK_TOT* node_data[nd].ltk * node_data[nd].ltj
              *node_data[nd].lti, MPI_DOUBLE, nd, 1, MPI_COMM_WORLD, &status);
    }
      RK_eqto_Model(RK_reform, model);
      memcpy(node_data[0].data[0][0][0], RK_reform[0][0][0], RK_TOT*NX3_TOT*NX2_TOT*NX1_TOT*sizeof(double));

  }else{
      RK_eqto_Model(RK_reform, model);
      MPI_Send(RK_reform[0][0][0], RK_TOT*NX3_TOT*NX2_TOT*NX1_TOT, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
  }
  //collect xuv data

  
  //reform data from node_data to g_data
  if(mynode==0){
    for(nd=0;nd<totalnodes;nd++){//reform grid value from node_data to g_data
      if(first_call){
        for(k=KBEG; k<= node_data[nd].lek ;k++) g_data.ck[node_data[nd].bk+k-2*KBEG] = node_data[nd].ck[k];
        for(j=JBEG; j<= node_data[nd].lej ;j++) g_data.cj[node_data[nd].bj+j-2*JBEG] = node_data[nd].cj[j];
        for(i=IBEG; i<= node_data[nd].lei ;i++) g_data.ci[node_data[nd].bi+i-2*IBEG] = node_data[nd].ci[i];
      }
    }

    for(nd=0;nd<totalnodes;nd++)
    RK_VAR_LOOP(nv)  
    for(k=KBEG; k<= node_data[nd].lek ;k++)
    for(j=JBEG; j<= node_data[nd].lej ;j++)
    for(i=IBEG; i<= node_data[nd].lei ;i++)
      g_data.data[nv][node_data[nd].bk+k-2*KBEG][node_data[nd].bj+j-2*JBEG][node_data[nd].bi+i-2*IBEG] 
        = node_data[nd].data[nv][k][j][i];
  }
    
  //write begin
  if(mynode==0){
    //write tab file
    if(g_inputParam[OUTO_NAME]) 
      sprintf(file_name,"%s/%s_fxuv%.1f_Index%.3f_He%.4f_Heat%.3f_%.4d_data.dat", 
              ini->output_dir, ini->obj_name, model->xuv.flux_tot/g_inputParam[XUV_FACTOR], 
              model->xuv.spectral_index_504A, g_inputParam[HE_RATIO],
              g_inputParam[HEATING_EFFICIENCY], step);
    else sprintf(file_name,"%s/%.4d_data.dat", folder, step);  

    print1("> %s out\n", file_name);
  
    fp = fopen(file_name, "w");

    for(i=0;i<nlabel;i++) fprintf(fp, "%-15s", labels[i]);
    fprintf(fp, "\n");

    G_DOM_LOOP(k,j,i){
      RD_EXPAND(
        fprintf(fp, "%-15d", k);,
        fprintf(fp, "%-15d", j);, 
        fprintf(fp, "%-15d", i);      
      )

      RD_EXPAND(
        fprintf(fp, "%-15.4e", g_data.ck[k]);,
        fprintf(fp, "%-15.4e", g_data.cj[j]);,
        fprintf(fp, "%-15.4e", g_data.ci[i]);
      )
      
      for(f=0;f<nfluid;f++){
        scrh = g_data.data[HD_NVAR*f+NDT][k][j][i];
        fprintf(fp, "%-15.4e", scrh);

        CPN_LOOP(nv)
          fprintf(fp, "%-15.4e", g_data.data[HD_NVAR*f+MX1+nv][k][j][i]/(scrh*ATOM_MASS(model,f)) );

        fprintf(fp, "%-15.4e", g_data.data[HD_NVAR*f+PRS][k][j][i]);
      }

    #if ELC_ON
      fprintf(fp, "%-15.4e", g_data.data[RK_ELC][k][j][i]);
    #endif

    #if MAG_ON
      CPN_LOOP(nv)
        fprintf(fp, "%-15.4e", g_data.data[RK_MAG+nv][k][j][i]*CONST_s4PI);
    #endif

      fprintf(fp, "\n");
    } 
    fclose(fp);//write tab end




    //write vtk file
    #if DIMENSIONS>=2

    static double ***d_store, ****vect_store, vect_name[128];  
    if(d_store==NULL){
      d_store = ARRAY_3D(KNP_GLOB, JNP_GLOB, INP_GLOB, double);
      vect_store = ARRAY_4D(3, KNP_GLOB, JNP_GLOB, INP_GLOB, double);
    }

    if(g_inputParam[OUTO_NAME]) 
      sprintf(file_name,"%s/%s_fxuv%.1f_Index%.3f_He%.4f_Heat%.3f_%.4d_data.vtk", 
              ini->output_dir, ini->obj_name, model->xuv.flux_tot/g_inputParam[XUV_FACTOR], 
              model->xuv.spectral_index_504A, g_inputParam[HE_RATIO],
              g_inputParam[HEATING_EFFICIENCY], step);
    else sprintf(file_name,"%s/%.4d_data.vtk", folder, step);
    
    print1("> %s out\n", file_name);

    fp = OpenBinaryFile(file_name, SZ_Float_Vect, "w");

    WriteVTK_Header(fp, grid);

    //WriteVTK_Scalar
    FLD_LOOP(f){
      G_DOM_LOOP(k,j,i) d_store[k][j][i] = g_data.data[HD_NVAR*f+NDT][k][j][i];
      WriteVTK_Scalar (fp, d_store, 1, labels[2*DIMENSIONS+HD_NVAR*f+NDT], grid);

      CPN_LOOP(nv){
        G_DOM_LOOP(k,j,i) d_store[k][j][i] = vect_store[nv][k][j][i] = 
                          g_data.data[HD_NVAR*f+MX1+nv][k][j][i]/(g_data.data[HD_NVAR*f+NDT][k][j][i]*ATOM_MASS(model,f));

        WriteVTK_Scalar (fp, d_store, 1, labels[2*DIMENSIONS+HD_NVAR*f+MX1+nv], grid);
      }

      sprintf(vect_name, "%s_velocity", model->fluid[f].ptc.label);
      WriteVTK_Vector (fp, vect_store, 1, vect_name, grid);

      G_DOM_LOOP(k,j,i) d_store[k][j][i] = g_data.data[HD_NVAR*f+PRS][k][j][i];
      WriteVTK_Scalar (fp, d_store, 1, labels[2*DIMENSIONS+HD_NVAR*f+PRS], grid);
    }
    //WriteVTK_Vector
    


    CPN_LOOP(nv){
      G_DOM_LOOP(k,j,i) d_store[k][j][i] = vect_store[nv][k][j][i] = g_data.data[RK_MAG+nv][k][j][i]*CONST_s4PI;
      WriteVTK_Scalar (fp, d_store, 1, labels[2*DIMENSIONS+RK_MAG+nv], grid);
    }

    WriteVTK_Vector (fp, vect_store, 1, "B", grid);


    
    



    CloseBinaryFile(fp, SZ_float);
    #endif

  }
  #else
  //single core output
  #endif
}


  first_call = 0;
  return;
}






/*
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

  //write xuv file
  xuv = &(model->xuv);
       
  sprintf(file_name, "%s/%.4d_xuv.dat",folder, step);
    
  printf("> %s out\n", file_name);

  fp = freopen(file_name, "w", fp);
   
    fprintf(fp, "%-12s", label[0]);
    fprintf(fp, "%-12s", label[1]);
    for(b=0;b<xuv->nbeam;b++){
      fprintf(fp, "%-12.2f", xuv->wavelength[b]);
    }
    fprintf(fp, "%-12.2s", "FLUX_TOT");
    fprintf(fp, "%s", "\n");


    //write xuv file
    xuv = &(model->xuv);
       
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
  }
 */









































