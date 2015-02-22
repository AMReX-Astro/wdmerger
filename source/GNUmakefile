SHELL=/bin/bash

# Default inputs and probin names

inputs=inputs
probin=probin

# Bring in the CASTRO inputs

include $(WDMERGER_HOME)/source/inputs

# Bring in the probin parameters

include $(WDMERGER_HOME)/source/probin

# Makefile options. We'll use the ?= method so that
# any including makefile can override these for custom tests
# (e.g. changing the number of dimensions or the EOS).

PRECISION ?= DOUBLE
PROFILE ?= FALSE

DEBUG ?= FALSE

USE_MPI ?= TRUE
USE_OMP ?= FALSE

DIM ?= 3

USE_GRAV ?= TRUE
USE_REACT ?= FALSE
USE_MODELPARSER ?= TRUE
USE_ROTATION ?= TRUE

# Override default compiler for certain machines.

UNAMEN := $(shell uname -n)

ifeq ($(findstring titan, $(UNAMEN)), titan)
	# OLCF Titan
	COMP ?= Cray
	FCOMP ?= Cray
else ifeq ($(findstring h2o, $(UNAMEN)), h2o)
	# NCSA Blue Waters
	COMP ?= Cray
	FCOMP ?= Cray
else
	COMP ?= g++
	FCOMP ?= gfortran
endif

# This sets the EOS directory in $(CASTRO_DIR)/EOS
EOS_dir     := helmeos

# This sets the network directory in $(NETWORK_HOME)
Network_dir := general_null
GENERAL_NET_INPUTS := $(CASTRO_DIR)/Networks/general_null/triple_alpha_plus_o.net 

# Determine where our source files are coming from
# based on the relevant file in the source directory.
# For building the full program, this will be irrelevant
# and we will simply compile everything in the 
# WDMERGER_HOME source directory.

diag_locs = $(WDMERGER_HOME)/source/diagnostics/
diag_pack = $(WDMERGER_HOME)/source/diagnostics/Make.package

root_locs = $(WDMERGER_HOME)/source/
root_pack = $(WDMERGER_HOME)/source/Make.package

prob_locs = ../source/
prob_pack = ../source/Make.package

ifneq ($(wildcard ../source/full_source),)
	Blocs := $(prob_locs) $(diag_locs)
	Bpack := $(prob_pack) $(diag_pack)
else
	Blocs := $(root_locs) $(diag_locs)
	Bpack := $(root_pack) $(diag_pack)
endif

BUILD_GIT_NAME := wdmerger
BUILD_GIT_DIR := $(WDMERGER_HOME)

# Make sure that we copy the inputs and probin files,
# and then build the CASTRO executable
# (this target is 'all' in the CASTRO makefile).

build: make_inputs make_probin all

make_inputs:

	if [ ! -e $(inputs) ] ; \
	then \
		if [ ! -e "../source/inputs" ] ; \
		then \
			touch $(inputs) ; \
			echo "amr.probin_file = $(amr.probin_file)" >> $(inputs) ; \
			echo "max_step = $(max_step)" >> $(inputs) ; \
			echo "stop_time = $(stop_time)" >> $(inputs) ; \
			echo "geometry.is_periodic = $(geometry.is_periodic)" >> $(inputs) ; \
			echo "geometry.coord_sys = $(geometry.coord_sys)" >> $(inputs) ; \
			echo "geometry.prob_lo = $(geometry.prob_lo)" >> $(inputs) ; \
			echo "geometry.prob_hi = $(geometry.prob_hi)" >> $(inputs) ; \
			echo "castro.center = $(castro.center)" >> $(inputs) ; \
			echo "castro.dt_cutoff = $(castro.dt_cutoff)" >> $(inputs) ; \
			echo "castro.cfl = $(castro.cfl)" >> $(inputs) ; \
			echo "castro.init_shrink = $(castro.init_shrink)" >> $(inputs) ; \
			echo "castro.change_max = $(castro.change_max)" >> $(inputs) ; \
			echo "castro.sum_interval = $(castro.sum_interval)" >> $(inputs) ; \
			echo "castro.lo_bc = $(castro.lo_bc)" >> $(inputs) ; \
			echo "castro.hi_bc = $(castro.hi_bc)" >> $(inputs) ; \
			echo "amr.n_cell = $(amr.n_cell)" >> $(inputs) ; \
			echo "amr.max_level = $(amr.max_level)" >> $(inputs) ; \
			echo "amr.ref_ratio = $(amr.ref_ratio)" >> $(inputs) ; \
			echo "amr.regrid_int = $(amr.regrid_int)" >> $(inputs) ; \
			echo "amr.n_error_buf = $(amr.n_error_buf)" >> $(inputs) ; \
			echo "amr.max_grid_size = $(amr.max_grid_size)" >> $(inputs) ; \
			echo "amr.blocking_factor = $(amr.blocking_factor)" >> $(inputs) ; \
			echo "amr.grid_eff = $(amr.grid_eff)" >> $(inputs) ; \
			echo "castro.allow_untagging = $(castro.allow_untagging)" >> $(inputs) ; \
			echo "castro.do_hydro = $(castro.do_hydro)" >> $(inputs) ; \
			echo "castro.do_grav = $(castro.do_grav)" >> $(inputs) ; \
			echo "castro.do_react = $(castro.do_react)" >> $(inputs) ; \
			echo "castro.do_sponge = $(castro.do_sponge)" >> $(inputs) ; \
			echo "castro.add_ext_src = $(castro.add_ext_src)" >> $(inputs) ; \
			echo "castro.do_rotation = $(castro.do_rotation)" >> $(inputs) ; \
			echo "castro.rotational_period = $(castro.rotational_period)" >> $(inputs) ; \
			echo "castro.ppm_type = $(castro.ppm_type)" >> $(inputs) ; \
			echo "castro.ppm_reference = $(castro.ppm_reference)" >> $(inputs) ; \
			echo "castro.ppm_reference_eigenvectors = $(castro.ppm_reference_eigenvectors)" >> $(inputs) ; \
			echo "castro.ppm_reference_edge_limit = $(castro.ppm_reference_edge_limit)" >> $(inputs) ; \
			echo "castro.use_colglaz = $(castro.use_colglaz)" >> $(inputs) ; \
			echo "castro.ppm_trace_grav = $(castro.ppm_trace_grav)" >> $(inputs) ; \
			echo "castro.ppm_trace_rot = $(castro.ppm_trace_rot)" >> $(inputs) ; \
			echo "castro.ppm_temp_fix = $(castro.ppm_temp_fix)" >> $(inputs) ; \
			echo "castro.grav_source_type = $(castro.grav_source_type)" >> $(inputs) ; \
			echo "castro.rot_source_type = $(castro.rot_source_type)" >> $(inputs) ; \
			echo "castro.small_temp = $(castro.small_temp)" >> $(inputs) ; \
			echo "castro.small_dens = $(castro.small_dens)" >> $(inputs) ; \
			echo "castro.allow_negative_energy = $(castro.allow_negative_energy)" >> $(inputs) ; \
			echo "castro.dual_energy_update_E_from_e = $(castro.dual_energy_update_E_from_e)" >> $(inputs) ; \
			echo "castro.dual_energy_eta1 = $(castro.dual_energy_eta1)" >> $(inputs) ; \
			echo "castro.dual_energy_eta2 = $(castro.dual_energy_eta2)" >> $(inputs) ; \
			echo "gravity.gravity_type = $(gravity.gravity_type)" >> $(inputs) ; \
			echo "gravity.max_multipole_order = $(gravity.max_multipole_order)" >> $(inputs) ; \
			echo "gravity.direct_sum_bcs = $(gravity.direct_sum_bcs)" >> $(inputs) ; \
			echo "gravity.ml_tol = $(gravity.ml_tol)" >> $(inputs) ; \
			echo "gravity.sl_tol = $(gravity.sl_tol)" >> $(inputs) ; \
			echo "castro.show_center_of_mass = $(castro.show_center_of_mass)" >> $(inputs) ; \
			echo "castro.print_energy_diagnostics = $(castro.print_energy_diagnostics)" >> $(inputs) ; \
			echo "castro.sum_interval = $(castro.sum_interval)" >> $(inputs) ; \
			echo "amr.plot_files_output = $(amr.plot_files_output)" >> $(inputs) ; \
			echo "amr.checkpoint_files_output = $(amr.checkpoint_files_output)" >> $(inputs) ; \
			echo "amr.check_file = $(amr.check_file)" >> $(inputs) ; \
			echo "amr.check_per = $(amr.check_per)" >> $(inputs) ; \
			echo "amr.plot_file = $(amr.plot_file)" >> $(inputs) ; \
			echo "amr.plot_per = $(amr.plot_per)" >> $(inputs) ; \
			echo "amr.v = $(amr.v)" >> $(inputs) ; \
			echo "castro.v = $(castro.v)" >> $(inputs) ; \
			echo "gravity.v = $(gravity.v)" >> $(inputs) ; \
			echo "mg.v = $(mg.v)" >> $(inputs) ; \
			echo "amr.derive_plot_vars = $(amr.derive_plot_vars)" >> $(inputs) ; \
			echo "castro.plot_phiGrav = $(castro.plot_phiGrav)">> $(inputs) ; \
			echo "amr.data_log = $(amr.data_log)" >> $(inputs) ; \
		else \
			cp ../source/inputs . ; \
		fi \
	fi

make_probin:

	if [ ! -e $(probin) ]; \
	then \
		if [ ! -e "../source/probin" ] ; \
		then \
			touch $(probin) ; \
			echo "&fortin" >> $(probin) ; \
			echo "" >> $(probin) ; \
			echo "  mass_P = $(mass_P)" >> $(probin) ; \
			echo "  mass_S = $(mass_S)" >> $(probin) ; \
			echo "  nsub = $(nsub)" >> $(probin) ; \
			echo "  orbital_kick = $(orbital_kick)" >> $(probin) ; \
			echo "  damping = $(damping)" >> $(probin) ; \
			echo "  damping_alpha = $(damping_alpha)" >> $(probin) ; \
			echo "  interp_temp = $(interp_temp)" >> $(probin) ; \
			echo "  bulk_velx = $(bulk_velx)" >> $(probin) ; \
			echo "  bulk_vely = $(bulk_vely)" >> $(probin) ; \
			echo "  bulk_velz = $(bulk_velz)" >> $(probin) ; \
			echo "  maxTaggingRadius = $(maxTaggingRadius)" >> $(probin) ; \
			echo "" >> $(probin) ; \
			echo "/" >> $(probin) ; \
			echo "" >> $(probin) ; \
			echo "&tagging" >> $(probin) ; \
			echo "" >> $(probin) ; \
			echo "  max_denerr_lev = $(max_denerr_lev)" >> $(probin) ; \
			echo "  denerr = $(denerr)" >> $(probin) ; \
			echo "  max_temperr_lev = $(max_temperr_lev)" >> $(probin) ; \
			echo "  temperr = $(temperr)" >> $(probin) ; \
			echo "" >> $(probin) ; \
			echo "/" >> $(probin) ; \
			echo "" >> $(probin) ; \
			echo "&extern" >> $(probin) ; \
			echo "" >> $(probin) ; \
			echo "  use_eos_coulomb = $(use_eos_coulomb)" >> $(probin) ; \
			echo "  eos_input_is_constant = $(eos_input_is_constant)" >> $(probin) ; \
			echo "" >> $(probin) ; \
			echo "/" >> $(probin) ; \
			echo "" >> $(probin) ; \
		else \
			cp ../source/probin . ; \
		fi \
	fi

# Clean executable up with make clean;
# clean up all build files including objects with make realclean;
# clean up everything, including all results and analysis, with make cleanall.

clean_targets = *.mod helm_table.dat $(inputs) $(probin) 
cleanall_targets = $(clean_targets) ../results/ ../plots/

# This is a hack that we need to tell BoxLib
# not to build the dependency list before we execute the rule.
# MAKECMDGOALS is the name of the rule we gave to make
# at the command line. We probably shouldn't be messing 
# with it, but this works.

ifeq ($(MAKECMDGOALS),cleanall)
	MAKECMDGOALS = realclean
endif

.PHONY cleanall: realclean	
	rm -rf $(cleanall_targets)

clean::
	rm -rf $(clean_targets)

include $(CASTRO_DIR)/Exec/Make.Castro
