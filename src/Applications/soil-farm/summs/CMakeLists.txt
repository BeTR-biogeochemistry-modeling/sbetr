set(SUMMS_SOURCES
  BgcSummsMath.F90
  BgcSummsDebType.F90
  BgcSummsIndexType.F90
  BgcSummsForcType.F90
  BgcSummsCompetType.F90
  BgcSummsDecompType.F90
  BgcSummsNitDenType.F90
  DebGrowMod.F90
  BgcSummsSOMType.F90
  BgcSummsType.F90
  SOMStateVarUpdateSumms.F90
  BgcReactionsSummsType.F90
  PlantSoilBgcSummsType.F90
)

if (BETR_BGC EQUAL 1)
  include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_util)
  include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_math)
  include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_bgc)
  include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_grid)
  include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_dtype)
  include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_core)
  include_directories(${CMAKE_BINARY_DIR}/src/Applications/soil-farm/bgcfarm_util)

  include(add_betr_library)
  add_betr_library(summs ${SUMMS_SOURCES})

  set(BETR_LIBRARIES summs;${BETR_LIBRARIES} PARENT_SCOPE)
  set(BETR_LIBRARIES summs;${BETR_LIBRARIES})

#X#add_subdirectory(tests)

  if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
    install(TARGETS summs DESTINATION lib)
    file(GLOB HEADERS *.h)
    install(FILES ${HEADERS} DESTINATION include/soil-farm/summs)
  endif()
else()
  include_directories(${CMAKE_BINARY_DIR}/src/betr_share/betr_util)
  include_directories(${CMAKE_BINARY_DIR}/src/betr_share/betr_math)
  include_directories(${CMAKE_BINARY_DIR}/src/betr_share/betr_grid)
  include_directories(${CMAKE_BINARY_DIR}/src/betr_share/betr_dtype)
  include_directories(${CMAKE_BINARY_DIR}/src/betr_share/betr_core)
  include_directories(${CMAKE_BINARY_DIR}/src/bgcfarm/bgcfarm_util)

  include(add_bgcfarm_library)
  add_bgcfarm_library(summs ${SUMMS_SOURCES})

  set(BGCFARM_LIBRARIES summs;${BGCFARM_LIBRARIES} PARENT_SCOPE)
  set(BGCFARM_LIBRARIES summs;${BGCFARM_LIBRARIES})

#add_subdirectory(tests)
  if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
    install(TARGETS summs DESTINATION lib)
    file(GLOB HEADERS *.h)
    install(FILES ${HEADERS} DESTINATION include/bgcfarm/summs)
  endif()
endif()
