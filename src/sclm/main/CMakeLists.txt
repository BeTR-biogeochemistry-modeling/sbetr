set(BETR_SCLM_SOURCES
  ColumnType.F90
  EcophysConType.F90
  GridcellType.F90
  LandunitType.F90
  PatchType.F90
  atm2lndType.F90
  clm_initializeMod.F90
  clm_instMod.F90
  clm_varcon.F90
  clm_varctl.F90
  clm_varpar.F90
  clmgridMod.F90
  column_varcon.F90
  landunit_varcon.F90
  pftvarcon.F90
  restUtilMod.F90
  soilorder_varcon.F90
  subgridAveMod.F90
  )

include_directories(${CMAKE_BINARY_DIR}/src/shr)

include(add_betr_library)
add_betr_library(betr_sclm ${BETR_SCLM_SOURCES})

set(BETR_LIBRARIES betr_sclm;${BETR_LIBRARIES})

#X#add_subdirectory(tests)

if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS betr_math DESTINATION lib)
  file(GLOB HEADERS *.h)
  install(FILES ${HEADERS} DESTINATION include/sclm/main)
endif()
