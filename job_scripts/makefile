SHELL=/bin/bash

merger_makefile=$(WDMERGER_HOME)/source/GNUmakefile
local_makefile=castro_makefile
inputs=inputs
probin=probin
run_utils=$(WDMERGER_HOME)/job_scripts/run_utils.sh

# Bring in the CASTRO inputs

include $(WDMERGER_HOME)/job_scripts/inputs

# Bring in the probin parameters

include $(WDMERGER_HOME)/job_scripts/probin

# Makefile options. We'll use the ?= method so that
# the including makefile can override these for custom tests
# (e.g. changing the number of dimensions or the EOS).

USE_MPI ?= TRUE
USE_OMP ?= FALSE
DIM ?= 3
EOS_dir ?= helmeos
USE_GRAV ?= TRUE
USE_ROTATION ?= TRUE
DEBUG ?= FALSE

include $(BOXLIB_HOME)/Tools/C_mk/Make.defs

ifeq ($(WHICHLINUX), BLUE_WATERS)
	COMP = Cray
	FCOMP = Cray
else ifeq ($(WHICHLINUX), TITAN)
	COMP = Cray
	FCOMP = Cray
else
	COMP = g++
	FCOMP = gfortran
endif

# Determine where our source files are coming from
# based on the relevant file in the source directory.

root_locs = $(WDMERGER_HOME)/source/
root_pack = $(WDMERGER_HOME)/source/Make.package

mods_locs = $(WDMERGER_HOME)/source/modules/
mods_pack = $(WDMERGER_HOME)/source/modules/Make.package

diag_locs = $(WDMERGER_HOME)/source/diagnostics/
diag_pack = $(WDMERGER_HOME)/source/diagnostics/Make.package

curr_locs = ../source/
curr_pack = ../source/Make.package

ifneq ($(wildcard ../source/full_source),)
	Blocs := $(curr_locs) $(diag_locs)
	Bpack := $(curr_pack) $(diag_pack)
else ifneq ($(wildcard ../source/partial_source),)
	Blocs := $(curr_locs) $(root_locs) $(mods_locs) $(diag_locs)
	Bpack := $(curr_pack) $(root_pack) $(mods_pack) $(diag_pack)
else ifneq ($(wildcard ../source/modules_only),)
	Blocs := $(curr_locs) $(mods_locs) $(diag_locs)
	Bpack := $(curr_pack) $(mods_pack) $(diag_pack)
else
	Blocs := $(root_locs) $(mods_locs) $(diag_locs)
	Bpack := $(root_pack) $(mods_pack) $(diag_pack)
endif

# Sets up the CASTRO makefile in the local directory
# and inserts our preferred make options like compiler
# and whether to use MPI and/or OpenMP. Should be called
# by the main make target in any run directory.

make_castro: make_inputs make_probin

	# Check if makefile exists; if not, copy it from the source directory,
	# then edit the makefile variables that we want to use, including
	# compiler options. If we want to use a custom setup for a particular
	# problem, instead of using the main source files, check for a file named
	# full_problem in the source directory.
	# The < and > characters in sed mark off word boundaries, and are useful
	# so that COMP does not overwrite FCOMP, which also matches.

	if [ ! -e $(local_makefile) ] ; \
	then \
		cp $(merger_makefile) $(local_makefile) ; \
		sed -i "/\<COMP\>/c COMP = $(COMP)" $(local_makefile) ; \
		sed -i "/\<FCOMP\>/c FCOMP = $(FCOMP)" $(local_makefile) ; \
		sed -i "/DIM/c DIM = $(DIM)" $(local_makefile) ; \
		sed -i "/EOS_dir/c EOS_dir := $(EOS_dir)" $(local_makefile) ; \
		sed -i "/USE_GRAV/c USE_GRAV = $(USE_GRAV)" $(local_makefile) ; \
		sed -i "/USE_ROTATION/c USE_ROTATION = $(USE_ROTATION)" $(local_makefile) ; \
		sed -i "/USE_MPI/c USE_MPI = $(USE_MPI)" $(local_makefile) ; \
		sed -i "/USE_OMP/c USE_OMP = $(USE_OMP)" $(local_makefile) ; \
		sed -i "/DEBUG/c DEBUG = $(DEBUG)" $(local_makefile) ; \
		sed -i "/Blocs/c Blocs := $(Blocs)" $(local_makefile) ; \
		sed -i "/Bpack/c Bpack := $(Bpack)" $(local_makefile) ; \
	fi

	make -f $(local_makefile)

make_inputs:

	if [ ! -e $(inputs) ] ; \
	then \
		touch $(inputs) ; \
		echo "amr.probin_file = $(amr.probin_file)" >> $(inputs) ; \
		echo "max_step = $(max_step)" >> $(inputs) ; \
		echo "stop_time = $(stop_time)" >> $(inputs) ; \
		echo "geometry.is_periodic = $(geometry.is_periodic)" >> $(inputs) ; \
		echo "geometry.coord_sys = $(geometry.coord_sys)" >> $(inputs) ; \
		echo "geometry.prob_lo = $(geometry.prob_lo)" >> $(inputs) ; \
		echo "geometry.prob_hi = $(geometry.prob_hi)" >> $(inputs) ; \
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
			echo "  inertial = $(inertial)" >> $(probin) ; \
			echo "  damping = $(damping)" >> $(probin) ; \
			echo "  damping_alpha = $(damping_alpha)" >> $(probin) ; \
			echo "  interp_temp = $(interp_temp)" >> $(probin) ; \
			echo "  starBuffer = $(starBuffer)" >> $(probin) ; \
			echo "  boundaryBuffer = $(boundaryBuffer)" >> $(probin) ; \
			echo "" >> $(probin) ; \
			echo "/" >> $(probin) ; \
			echo "" >> $(probin) ; \
			echo "&tagging" >> $(probin) ; \
			echo "" >> $(probin) ; \
			echo "  max_denerr_lev = $(max_denerr_lev)" >> $(probin) ; \
			echo "  denerr = $(denerr)" >> $(probin) ; \
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

# Clean everything up with make clean

clean_targets = *.mod *.run helm_table.dat $(inputs) $(probin) 
realclean_targets = ../results/ ../plots/

clean:
	if [ -e $(local_makefile) ]; then \
		make realclean -f $(local_makefile); \
		rm -f $(local_makefile); \
	fi

	rm -rf $(clean_targets)

realclean: clean
	rm -rf $(realclean_targets)

# Prints a variable VAR when makefile is called as
# make print-VAR

print-%: ; @echo $* is $($*)
