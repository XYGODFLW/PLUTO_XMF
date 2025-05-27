//mhd_convection_set prototypes and some definitions

int mhd_AdvanceStep(const Data *, Riemann_Solver *, Time_Step *, Grid *);
void mhd_UpdateStage(const Data *, Data_Arr, double **, Riemann_Solver *,
                 double, Time_Step *, Grid *);

int mhd_Integrate (Data*, Riemann_Solver *, Time_Step *, Grid *);
double mhd_NextTimeStep (Time_Step *, Grid *);

void mhd_PrimToCons3D(Data_Arr, Data_Arr, RBox *);
void mhd_ConsToPrim3D(Data_Arr, Data_Arr, unsigned char ***, RBox *);

void mhd_SetIndexes (Index *indx, Grid *grid);

void mhd_SetOutput (Data *d, Runtime *input);

void mhd_Startup (Data *, Grid *);

//commen functions
int  SetDumpVar (char *, int, int);
double  ***GetUserVar (char *);




