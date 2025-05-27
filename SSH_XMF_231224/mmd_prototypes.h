#include "pluto.h"

/* *******************
 (struct reference)
 CMD_LINE
 DATA
 STATE_1D
 GRID
 TIME_STEP
 OUTPUT
 RUNTIME
 RESTART
 RGB
 IMAGE
 FLOAT_VECT
 INDEX
 RBOX
 *********************/
void mmd_initialize_public(int, char **, Runtime *, Grid *, Cmd_Line *);
void mhd_initialize(Fluid *, Runtime *, Grid *, Cmd_Line *);
void hd_initialize(Fluid *, Runtime *, Grid *, Cmd_Line *);

void mmd_time_step_fluids_ini(Runtime *, Time_Step *);
int  mmd_AdvanceStep (Model *, Grid *);

void AssignFluidsAttributes(Fluid *, Grid *);

//mmd restart
void mmd_write_save(double****, char*);

void mmd_load_save(char*, double****);

//mmd model
void mmd_model_initialize(Model*, Runtime*);
void mmd_model_IC(Model*, Grid*, Runtime*);
void mmd_model_BC(Model*, Grid*);
void mmd_model_output(Model*, Grid*, Runtime*ini, double, char*);
void mmd_get_lf(Model*, Grid*, double ****);
void mmd_get_lc(Model*, Grid*, double ****, double ****, int, int, int);
void quasi_neutrality(Model *);

 //xuv and sed...
void xuv_initialize(Xuv *, char* ,char*);
void get_xuv_flux(Model *, Grid *);

//RK_tools for rk and 
void RK_multiplyAdd(double ****, double , double ****, double a1, double ****);
void RK_eqto_Model(double****, Model*);
void Model_eqto_RK(Model*, double****);
void RK_lockTo_Model(double ****, Model *);
























