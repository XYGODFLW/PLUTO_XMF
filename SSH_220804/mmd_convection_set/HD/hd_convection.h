//mhd_convection_set prototypes and some definitions

int hd_AdvanceStep(const Data *, Riemann_Solver *, Time_Step *, Grid *);
void hd_UpdateStage(const Data *, Data_Arr, double **, Riemann_Solver *,
                 double, Time_Step *, Grid *);

int hd_Integrate (Data *, Riemann_Solver *, Time_Step *, Grid *);
double hd_NextTimeStep (Time_Step *, Grid *);

void hd_PrimToCons3D(Data_Arr, Data_Arr, RBox *);
void hd_ConsToPrim3D(Data_Arr, Data_Arr, unsigned char ***, RBox *);

void hd_SetIndexes (Index *indx, Grid *grid);

void hd_Boundary   (const Data *, int, Grid *);

void hd_Startup (Data *, Grid *);
void hd_SetOutput (Data *d, Runtime *input);

void hd_SetDefaultVarNames(Output *);





