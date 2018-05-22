&sbetr_driver
  simulator_name = 'standalone'
  run_type='sbgc'
  param_file='../tools/sbgc.ecacnp_pars.04232018.nc'
  case_id='exp.noadv'
  continue_run=.false.
/

&betr_parameters
  reaction_method = 'simic'
  advection_on = .true.
  diffusion_on = .true.
  reaction_on  = .true.
  ebullition_on= .true.
  input_only   = .false.
/

&betr_time
  delta_time=1800.
  stop_n = 60
  hist_freq=30
  stop_option='nyears'
/

&forcing_inparm
  forcing_filename = '../input_data/single_pt_example.bgc.halfhour.forc.nc'
  use_rootsoit=.false.
/

&betr_grid
  grid_data_filename = '../input_data/single_pt_example.grid.cdl.nc'
/

&regression_test
/
