set(BETR_SCLM_SOURCES
  CNCarbonFluxType.F90
  CNCarbonStateType.F90
  CNDecompCascadeConType.F90
  CNNitrogenFluxType.F90
  CNNitrogenStateType.F90
  CNStateType.F90
  CanopyStateType.F90
  ChemStateType.F90
  ColumnType.F90
  EcophysConType.F90
  GridcellType.F90
  LandunitType.F90
  PatchType.F90
  PhosphorusFluxType.F90
  PhosphorusStateType.F90
  SoilHydrologyType.F90
  SoilStateType.F90
  TemperatureType.F90
  WaterFluxType.F90
  WaterStateType.F90
  abortutils.F90
  accumulMod.F90
  atm2lndType.F90
  clm_initializeMod.F90
  clm_instMod.F90
  clm_time_manager.F90
  clm_varcon.F90
  clm_varctl.F90
  clm_varpar.F90
  clmgridMod.F90
  column_varcon.F90
  debuginfoMod.F90
  decompMod.F90
  fileutils.F90
  getdatetime.F90
  landunit_varcon.F90
  ncdio_pio.F90
  pftvarcon.F90
  restUtilMod.F90
  soilorder_varcon.F90
  spmdMod.F90
  subgridAveMod.F90
  )

include_directories(${CMAKE_BINARY_DIR}/src/shr)
include_directories(${CMAKE_BINARY_DIR}/src/esmf_wrf_timemgr)

include(add_betr_library)
add_betr_library(betr_sclm ${BETR_SCLM_SOURCES})

set(BETR_LIBRARIES betr_sclm;${BETR_LIBRARIES})

#X#add_subdirectory(tests)

if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS betr_math DESTINATION lib)
  file(GLOB HEADERS *.h)
  install(FILES ${HEADERS} DESTINATION include/sclm/main)
endif()
