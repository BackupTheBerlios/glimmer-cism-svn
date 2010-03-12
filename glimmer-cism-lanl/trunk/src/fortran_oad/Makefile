F90_FILES= \
 glimmer_vers.F90 \
 glimmer_global.F90 \
 glimmer_sparse_type.F90 \
 glimmer_log.F90 \
 glimmer_sparse_slap.F90 \
 glimmer_sparse_umfpack.F90 \
 glimmer_sparse_pardiso.F90 \
 glimmer_filenames.F90 \
 glimmer_sparse.F90 \
 glimmer_physcon.F90 \
 glimmer_paramets.F90 \
 glimmer_ncdf.F90 \
 isostasy_types.F90 \
 profile.F90 \
 glimmer_coordinates.F90 \
 glimmer_map_types.F90 \
 glide_grids.F90 \
 xls.F90 \
 remap_glamutils.F90 \
 glide_types.F90 \
 glide_deriv.F90 \
 glide_nonlin.F90 \
 glide_mask.F90 \
 ice3d_lib.F90 \
 glam_strs2.F90 \
 glide_vertint.F90 \
 glide_thckmask.F90 \
 glide_velo_higher.F90 \
 fo_upwind_advect.F90 \
 glimmer_config.F90 \
 glimmer_map_init.F90 \
 glimmer_map_CFproj.F90 \
 glimmer_ncio.F90 \
 glide_io.F90 \
 glide_stop.F90 \
 glide_nc_custom.F90 \
 glide_setup.F90 \
 glimmer_utils.F90 \
 glide_lithot1d.F90 \
 glide_lithot3d.F90 \
 glide_lithot.F90 \
 glide_profile.F90 \
 isostasy_setup.F90 \
 kelvin.F90 \
 isostasy_el.F90 \
 isostasy.F90 \
 glimmer_ncparams.F90 \
 glide_velo.F90 \
 glide_thck.F90 \
 glide_bwater.F90 \
 glide_temp.F90 \
 glide_ground.F90 \
 remap_advection.F90 \
 glam.F90 \
 glide_stress.F90 \
 glimmer_restart_common.F90 \
 glimmer_daily_pdd.F90 \
 glimmer_integrate.F90 \
 glimmer_pdd.F90 \
 glimmer_ts.F90 \
 glimmer_restart.F90 \
 glide.F90 \
 glide_diagnostics.F90 \
 glissade_constants.F90 \
 glissade_remap.F90 \
 glissade.F90 \
 glissade_velo.F90 \
 glam_thck_ppm.F90 \
 glimmer_commandline.F90 \
 glimmer_routing.F90 \
 glint_global_grid.F90 \
 glint_global_interp.F90 \
 ncdf_utils.F90 \
 glimmer_map_proj4.F90 \
 glimmer_map_trans.F90 \
 glimmer_restart_statscal.F90 \
 glimmer_restart_statarr.F90 \
 glimmer_restart_pointarr.F90 \
 glimmer_anomcouple.F90 \
 glimmer_writestats.F90 \
 glint_mpinterp.F90 \
 glint_interp.F90 \
 glint_smb.F90 \
 glint_constants.F90 \
 glint_mbal.F90 \
 glint_mbal_coupling.F90 \
 glint_type.F90 \
 glint_precip_param.F90 \
 glint_climate.F90 \
 glint_io.F90 \
 glint_mbal_io.F90 \
 glint_initialise.F90 \
 glint_timestep.F90 \
 glint_main.F90 \
 glint_example_clim.F90 \
 glex_ebm_clim.F90 \
 glint_commandline.F90 \
 eismint3_types.F90 \
 eismint3_io.F90 \
 eismint3_forcing.F90 \
 testsfg.F90 \
 simple_forcing.F90 \
 glimmer_searchcircle.F90 \
 eis_cony.F90 \
 eis_ela.F90 \
 eis_temp.F90 \
 eis_slc.F90 \
 eis_types.F90 \
 eis_io.F90 \
 eis_forcing.F90 \
 spin_temp.F90 \
 spin_mb.F90 \
 spin_slc.F90 \
 spin_types.F90 \
 spin_io.F90 \
 spin_forcing.F90 \
 simple_glide.F90 \
 eis_glide.F90 \
 glint_example.F90 \
 glex_ebm.F90 \
 eismint3_glide.F90 \
 nc2config.F90 \
 spin_glide.F90 \
 thermal_equib.F90 

INC_FILES=config.inc glide_mask.inc  glide_nan.inc  glimmer_rst_body.inc  glimmer_rst_head.inc

default: sources numCore.pre.xb.x2w.w2f.post.f90

numCore.pre.xb.x2w.w2f.post.f90 : toolChain

toolChain : numCore.f90
	openad -c -m f $<

numCore.f90:  iso_c_binding.f90 $(F90_FILES:.F90=.f90)
	cat $^ > $@

sources: $(INC_FILES) $(F90_FILES:.F90=.f90)


config.inc:
	ln -sf ../../$@ $@

%.inc: 
	ln -sf ../fortran/$@ $@

%.F90: 
	ln -sf ../fortran/$@ $@

%.f90: %.F90 
	gfortran -o $@ -I ./ -E -Xpreprocessor -P $^

clean: 
	rm -f $(F90_FILES:.F90=.f90) $(F90_FILES) *.inc 
	rm -f numCore.f90 *.mod
