&jar_driver
  jarmodel_name = 'ecacnp'
  phosphorus_stress=.false.
  nitrogen_stress=.true.
  case_id='exp1'
  hist_freq='day'
/

&betr_time
  delta_time=1800.
  stop_n = 30
  stop_option='nyears'
/

&forcing_inparm
  forcing_filename = '../input_data/jarmodel_example.halfhour.forcing.nc'
/
