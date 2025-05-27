//prototype
void mmd_sm_initialize(Runtime* ini);

//ports declarition
extern Data sm_data;
extern double sm_dt, ****sm_du;    
extern Time_Step sm_dts;
extern char ***g_offpoint;
extern Riemann_Solver *sm_riemann;




