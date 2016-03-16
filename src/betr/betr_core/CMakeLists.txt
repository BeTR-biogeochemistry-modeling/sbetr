set(BETR_CORE_SOURCES
  BeTR_CNStateType.F90
  BeTR_CarbonFluxType.F90
  BeTR_aerocondType.F90
  BeTR_WaterstateType.F90
  BeTR_WaterfluxType.F90
  BeTR_TemperatureType.F90
  BeTR_SoilHydrologyType.F90
  BeTR_atm2lndType.F90
  BeTR_CanopyStateType.F90
  BeTR_ChemStateType.F90
  BeTR_SoilStateType.F90
  BeTRTracerType.F90
  KineticsMod.F90
  TracerBoundaryCondType.F90
  TracerCoeffType.F90
  TracerFluxType.F90
  TracerParamsMod.F90
  TracerStateType.F90
  Tracer_varcon.F90
  TransportMod.F90
  )

include_directories(${CMAKE_BINARY_DIR}/src/stub_clm)
include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_util)
include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_math)
include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_grid)

include(add_betr_library)
add_betr_library(betr_core ${BETR_CORE_SOURCES})

set(BETR_LIBRARIES betr_core;${BETR_LIBRARIES} PARENT_SCOPE)
set(BETR_LIBRARIES betr_core;${BETR_LIBRARIES})

add_subdirectory(tests)

if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS betr_core DESTINATION lib)
  file(GLOB HEADERS *.h)
  install(FILES ${HEADERS} DESTINATION include/betr/core)
endif()
