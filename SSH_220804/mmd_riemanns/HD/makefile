# Makefile for the HydroDynamics (HD) Module
#
# Note that we add Src/MHD to VPATH since HD and MHD
# share some files (see below)

VPATH += mmd_riemanns/HD

INCLUDE_DIRS += -Immd_riemanns/HD

HEAD += hd_mod_defs.h

OBJ +=  hd_eigenv.o hd_fluxes.o hd_mappers.o  hd_prim_eqn.o \
        hd_hll_speed.o hd_set_solver.o \
        hd_hll.o  hd_hllc.o hd_tvdlf.o hd_roe.o   

OBJ +=  hd_get_area_flux.o hd_rhs.o hd_rhs_source.o

