/* ------------------------------------- 
 * Aim to build a comprehensive radiative HD(MHD) code
 *
 * Struct tree:
 * Model:{Fluid, XUV, B(future), E(future), ...}
 *   Fluid:{status, var}
 *   XUV  :{status, var}
 * 
 * 2021-3 Xing
 * ------------------------------------- */




/* Fluid related structs */
typedef struct PARTICAL{
  char* ptc_label;  
  int ptc_element;
  int ptc_charge;
  int ptc_type;
  double ptc_m;//in_code_unit

  double E_weight; //  FE_weight = charge/ptc_m, charge density. for electron this value is 1
                   //  Pe impact on this ion fluid : Pe/RHOC_e * rho * E_weight
                   //  RHOC_e = SIGMA(rho[i]*FE_weight[i])
                   //  (PB may have similar form?)
  //int energy_level;

}Partical;


typedef struct FLUID{
  int f;
  Partical ptc;
  double ****Vc; 

}Fluid;


/* XUV related structs */
typedef struct X_STATUS{
  int nbeam;
  int wavelength;
}X_status;



typedef struct XUV{
  Intro infro;
  X_status status;

  X_var sed;
  double ****flux;  
}Xuv;


/* MODEL struct */
typedef struct M_STATUS{
  int nfluids;
  
}M_status;


typedef struct MODEL{
  Intro infro;
  M_status;

  Fluid *fluid;
  Xuv xuv;
  
}Model;




























