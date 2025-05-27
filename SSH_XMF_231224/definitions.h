#define  PHYSICS                 MHD
#define  DIMENSIONS              2
#define  COMPONENTS              2
#define  GEOMETRY                SPHERICAL
#define  BODY_FORCE              NO
#define  COOLING                 NO
#define  RECONSTRUCTION          WENO3
#define  TIME_STEPPING           RK3
#define  DIMENSIONAL_SPLITTING   NO
#define  NTRACER                 0
#define  USER_DEF_PARAMETERS     32

/* -- physics dependent declarations -- */

#define  EOS                     IDEAL
#define  ENTROPY_SWITCH          NO
#define  DIVB_CONTROL            EIGHT_WAVES //CONSTRAINED_TRANSPORT //#1D model is not compatible with CT.
#define  BACKGROUND_FIELD        NO
#define  RESISTIVITY             NO
#define  THERMAL_CONDUCTION      NO
#define  VISCOSITY               NO
#define  ROTATING_FRAME          NO

/* -- user-defined parameters (labels) -- */

#define  R_Planet                0
#define  M_Planet                1
#define  R_Orbit                 2
#define  M_Star                  3
#define  Y_Star                  4
#define  LB_RHO                  5
#define  LB_TMP                  6
#define  LB_B                    7
#define  GAMMA                   8
#define  H_ION_RATE              9
#define  XUV_FACTOR              10
#define  BASE_ION_RATE           11
#define  NDT_BEG                 12
#define  TREND_NEU               13
#define  TREND_ION               14
#define  TMP_BEG                 15
#define  HE_RATIO                16
#define  XUV_TOT_FLUX            17
#define  SPECTRAL_INDEX_504A     18
#define  HEATING_EFFICIENCY      19
#define  IC_FROM_FILE            20
#define  OUTO_NAME               21
#define  XUV_RESIZE              22
#define  TREND_HE                23
#define  TREND_H                 24
#define  OUTPUT_INTERVAL         25
#define  G_INI                   26
#define  G_FINAL                 27
#define  G_TIME                  28
#define  CONTINUE                29
#define  PARA_14                 30
#define  PARA_15                 31
#define  PARA_16                 32
#define  PARA_17                 33


/* [Beg] user-defined constants (do not change this line) */

#define  UNIT_DENSITY            1.e-12
#define  UNIT_LENGTH             (CONST_Rjup*g_inputParam[R_Planet])
#define  UNIT_VELOCITY           1.e6
#define  UNIT_TIME               (UNIT_LENGTH/UNIT_VELOCITY)
#define  UNIT_MMT                (UNIT_DENSITY*UNIT_VELOCITY)
#define  UNIT_ENG                (UNIT_MMT*UNIT_VELOCITY)
#define  UNIT_B                  (pow(UNIT_ENG, 0.5)) 

#define  CONST_Rjup              7.1492e9             /**<  Jupiter Radius.   */
#define  CONST_Mjup              1.8982e30            /**<  Jupiter Mass      */
#define  CONST_s4PI              3.544907701811032054 /**<  sqrt(4*PI)        */                 


/* [End] user-defined constants (do not change this line) */


/* -- supplementary constants (user editable) -- */ 

#define  INITIAL_SMOOTHING         NO
#define  WARNING_MESSAGES          NO
#define  PRINT_TO_FILE             YES
#define  INTERNAL_BOUNDARY         NO
#define  SHOCK_FLATTENING          NO
#define  CHAR_LIMITING             NO
#define  LIMITER                   DEFAULT
#define  CT_EMF_AVERAGE            UCT_HLL
#define  CT_EN_CORRECTION          NO
#define  ASSIGN_VECTOR_POTENTIAL   NO
#define  UPDATE_VECTOR_POTENTIAL   NO
