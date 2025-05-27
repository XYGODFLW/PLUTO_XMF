# Makefile for the Magnetohydrodynamics (MHD) Module 

VPATH        +=  mmd_riemanns/MHD
INCLUDE_DIRS +=  -Immd_riemanns/MHD

HEADERS += mhd_mod_defs.h

OBJ   +=  mhd_bckgrnd_field.o mhd_eigenv.o mhd_fluxes.o mhd_mappers.o mhd_prim_eqn.o \
          mhd_hll.o mhd_hllc.o mhd_hlld.o mhd_hll_speed.o mhd_roe.o mhd_tvdlf.o \
          mhd_set_solver.o mhd_source.o    mhd_rhs.o


# !!! NEW !!!

OBJ += mhd_get_area_flux.o   mhd_rhs_source.o


