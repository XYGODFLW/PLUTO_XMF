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
#include"aly.h"
int read_log(FILE* , char log[256][128], int *end);



int read_log(FILE* fp, char log[256][128], int *log_len)
{
  int k,j,i,end;
  char t_c;
  
  //get var_num;
  i = 0, j = 0, end = 0, log[0][0]=='\0';
  while(1){
    t_c = fgetc(fp);

    while(t_c==','||t_c==' '){
      if(t_c==','){log[j++][i] = '\0',i=0;}
      t_c = fgetc(fp);

      if(t_c=='\n'||t_c==EOF){if(t_c==EOF) end=1; j--;break;}
    }

    if(t_c == '\n'||t_c==EOF){if(t_c==EOF) end=1; break;}
    log[j][i++] = t_c;
  }
  
  (*log_len) = j+1;
  if(!j&&(log[0][0]=='\0')) (*log_len) = 0;

  return end||(!(*log_len));
}


int csv_read(char *file_name, CSV *d0)
{
  FILE *fp;
  char t_c, log[256][128];
  int k,j,i, nvar;

  //Check CSV file, find data_size   
  fp = fopen(file_name, "r");

  //initialize for d0
  j=0, d0->head=1;

  if(d0->data==NULL){
    while(!read_log(fp, log, &nvar)){j++, d0->nvar = nvar;}
    d0->nlog = j-1;//If no head, d0->nlog=j
    d0->data = (double**)malloc(d0->nlog*sizeof(double*));
    for(j=0;j<d0->nlog;j++)d0->data[j] = (double*)malloc(d0->nvar*sizeof(double));
  }
  //load data
  rewind(fp);

  read_log(fp, log, &nvar);
  for(i=0;i<d0->nvar;i++){ 
    sscanf(log[i], "%s", d0->var_names[i]);
  }
  
  for(j=0;j<d0->nlog;j++){
    read_log(fp, log, &nvar);
    for(i=0;i<d0->nvar;i++){ 
      sscanf(log[i], "%lf", d0->data[j]+i);
    }
  }

  fclose(fp);
  return 0;
}



int csv_write(char *file_name ,CSV *d)
{}
























