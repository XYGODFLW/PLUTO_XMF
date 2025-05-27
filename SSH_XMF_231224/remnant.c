/*
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
    printf("origin data, node:%d, IBEG = %d, IEND = %d\n", mynode, IBEG, IEND);
  }
  */
  /*
  MPI_Barrier(MPI_COMM_WORLD);
  if(first_call)if(mynode==0)for(node=0;node<totalnodes;node++){
    printf("node:%d, ibeg:%d, iend:%d \n", node, ibeg_list[node], iend_list[node]);
  }
  */
