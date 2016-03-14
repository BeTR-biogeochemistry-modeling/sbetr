set(BETR_DTYPE_SOURCES
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
  )

include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_util)


include(add_betr_library)
add_betr_library(betr_dtype ${BETR_DTYPE_SOURCES})

set(BETR_LIBRARIES betr_dtype;${BETR_LIBRARIES} PARENT_SCOPE)
set(BETR_LIBRARIES betr_dtype;${BETR_LIBRARIES})

#X#add_subdirectory(tests)

if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS betr_dtype DESTINATION lib)
  file(GLOB HEADERS *.h)
  install(FILES ${HEADERS} DESTINATION include/betr/dtype)
endif()
