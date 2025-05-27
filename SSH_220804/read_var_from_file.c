#include"pluto.h"

static char streqr(char*, char*);
static char next_line(FILE*);

void read_file_int(char* file_name, char* term, int* output)
{
  char s[128];  
  FILE *fp;
    
  if((fp = fopen(file_name, "r")) == NULL){
    printf("read_line_int: error in reading file %s \n", file_name);return;
  }

  while(1){

    fscanf(fp,"%s",s);

    if(!strcmp(s,term)){
      fscanf(fp, "%d", output);
      fclose(fp);
      return;
    }
    if(next_line(fp)){printf("> error: term %s is NG\n", term); fclose(fp); return;}

  }

}

void read_file_str(char* file_name, char* term, char* output)
{
  char s[128];  
  FILE *fp;
    
  if((fp = fopen(file_name, "r")) == NULL){
    printf("read_line_int: error in reading file %s \n", file_name);return;
  }

  while(1){

    fscanf(fp,"%s",s);

    if(!strcmp(s,term)){
      fscanf(fp, "%s", output);
      fclose(fp);
      return;
    }
    if(next_line(fp)){printf("> error: term %s is NG\n", term); fclose(fp); return;}

  }
}

void read_file_dbl(char* file_name, char* term, double* output)
{
  char s[128];  
  FILE *fp;
    
  if((fp = fopen(file_name, "r")) == NULL){
    printf("read_line_int: error in reading file %s \n", file_name);return;
  }

  while(1){

    fscanf(fp,"%s",s);

    if(!strcmp(s,term)){
      fscanf(fp, "%lf", output);
      fclose(fp);
      return;
    }
    if(next_line(fp)){printf("> error: term %s is NG\n", term); fclose(fp); return;}

  }
}
//can be rewritten to macro

static char next_line(FILE* fp)
{
  while( (fgetc(fp)!='\n') && (!feof(fp)) );
  return feof(fp);    
}























