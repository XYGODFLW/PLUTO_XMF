# *********************************************************
#                                                     
#                  PLUTO 4.2  Makefile  
#                                                     
# *********************************************************

pluto:                              # Default target

ARCH         = Linux.mpicc.defs
PLUTO_DIR    = /home/guojh/SSH_220804/PLUTO
SRC          = $(PLUTO_DIR)/Src
INCLUDE_DIRS = -I. -I$(SRC) 
VPATH        = ./:$(SRC)/New:$(SRC):$(SRC)/Time_Stepping

include $(PLUTO_DIR)/Config/$(ARCH)

# ---------------------------------------------------------
#         Set headers and object files 
# ---------------------------------------------------------
HEADERS = pluto.h prototypes.h structs.h definitions.h macros.h plm_coeffs.h \
          mmd.h mmd_prototypes.h mmd_structs.h \
          matrix_tools.h
      

OBJ = adv_flux.o arrays.o check_states.o  \
      cmd_line_opt.o entropy_switch.o  \
      flag_shock.o get_nghost.o   \
      int_bound_reset.o input_data.o \
      parse_file.o plm_coeffs.o rbox.o \
      set_geometry.o tools.o \

OBJ += bin_io.o colortable.o jet_domain.o \
       main.o restart.o runtime_setup.o show_config.o  \
       set_image.o set_grid.o split_source.o \
       userdef_output.o write_data.o write_tab.o \
       write_img.o write_vtk.o 

OBJ += com_initialize.o \
       mmd_model_initialize.o mmd_model_IC_BC.o \
       mmd_output.o xuv_initialize.o

OBJ += mmd_integrate.o mmd_model_adv.o mmd_model_src.o

include $(SRC)/Math_Tools/makefile

# ---------------------------------------------------------
#  Define macros by adding -D<name> where <name> has been
#  set to TRUE in the system configuration file (.defs) 
# ---------------------------------------------------------

ifeq ($(strip $(PARALLEL)), TRUE)
 CFLAGS += -I$(SRC)/Parallel -DPARALLEL
 include $(SRC)/Parallel/makefile
endif

ifeq ($(strip $(USE_HDF5)), TRUE)
 CFLAGS += -DUSE_HDF5
 OBJ    += hdf5_io.o
endif
      
ifeq ($($strip $(USE_PNG)), TRUE)
 CFLAGS += -DUSE_PNG
endif

-include local_make

# ---------------------------------------------------------
#   Additional_header_files_here   ! dont change this line
# ---------------------------------------------------------


# ---------------------------------------------------------
#   Additional_object_files_here   ! dont change this line
# ---------------------------------------------------------

OBJ += vec_pot_diff.o
OBJ += vec_pot_update.o

include mmd_convection_set/MHD/makefile
include mmd_convection_set/HD/makefile

include mmd_riemanns/MHD/makefile
include mmd_riemanns/MHD/CT/makefile

include mmd_riemanns/HD/makefile


include mmd_states/MHD/makefile
include mmd_states/HD/makefile


include $(SRC)/EOS/Ideal/makefile

#common tools

include com_derivative_tools/makefile
include com_matrix_tools/makefile

#state machine 
include mmd_state_machine/makefile

# ---------------------------------------------------------
#    PLUTO target rule
# ---------------------------------------------------------

pluto: $(OBJ) 
	$(CC) $(OBJ) $(LDFLAGS) -o $@

# ---------------------------------------------------------
#                    Suffix rule
# ---------------------------------------------------------

.c.o:
	$(CC) $(CFLAGS) $(INCLUDE_DIRS) $<

clean:
	@rm -f	*.o
	@echo make clean: done

# ---------------------------------------------------------
#          Dependencies for object files
# ---------------------------------------------------------

$(OBJ):  $(HEADERS)

