&sbetr_driver
  simulator_name = 'standalone'
  run_type='sbgc'
  case_id='exp.noadv'
  continue_run=.false.
  is_nitrogen_active=.false.
  is_phosphorus_active =.false.
/

&betr_parameters
  reaction_method = 'ecacnp'
  advection_on = .true.
  diffusion_on = .true.
  reaction_on  = .true.
  ebullition_on= .true.
  input_only   = .false.
  bgc_param_file='../tools/sbgc.ecacnp_pars.06292018.nc'
/

&betr_time
  delta_time=1800.
  stop_n = 5
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
