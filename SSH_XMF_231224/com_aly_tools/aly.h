/*
  Only for number data csv
  csv file format
  (var_name:) ..,..,..,
  (datas:)    ..,..,..,

              ......

              ..,..,..,
  by Xing,L 2023.9 (xinglei@ynao.ac.cn)
*/

typedef struct{
  char data_name[128];
  char var_names[256][128];
  int nvar, nlog, head;
  double **data;
} CSV;

double line_ITP(double **data, int n, double pst);
int csv_read(char*, CSV *);
int csv_write(char*, CSV *);























