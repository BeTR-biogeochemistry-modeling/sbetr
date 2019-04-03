set(SUMMS1LAYER_SOURCES
  BgcSummsIndexType.F90
  BgcSummsCompetType.F90
  BgcSummsDecompType.F90
  BgcSummsNitDenType.F90
  BgcSummsSOMType.F90
  BgcSummsType.F90
  BgcSummsMath.F90
  BgcSummsDebType.F90
  DebGrowMod.F90
)


  include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_util)
  include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_math)
  include_directories(${CMAKE_BINARY_DIR}/src/Applications/soil-farm/bgcfarm_util)
  include_directories(${CMAKE_BINARY_DIR}/src/Applications/soil-farm/SUMMS/summsPara)
  include(add_betr_library)
  add_betr_library(summs1layer ${SUMMS1LAYER_SOURCES})

  set(BETR_LIBRARIES summs1layer;${BETR_LIBRARIES} PARENT_SCOPE)
  set(BETR_LIBRARIES summs1layer;${BETR_LIBRARIES})

#X#add_subdirectory(tests)

  if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
    install(TARGETS summs1layer DESTINATION lib)
    file(GLOB HEADERS *.h)
    install(FILES ${HEADERS} DESTINATION include/soil-farm/SUMMS/summs1layer)
  endif()
||||||| merged common ancestors
=======
set(SUMMS1LAYER_SOURCES
  BgcSummsIndexType.F90
  BgcSummsCompetType.F90
  BgcSummsDecompType.F90
  BgcSummsNitDenType.F90
  BgcSummsSOMType.F90
  BgcSummsType.F90
)


  include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_util)
  include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_math)
  include_directories(${CMAKE_BINARY_DIR}/src/Applications/soil-farm/bgcfarm_util)
  include_directories(${CMAKE_BINARY_DIR}/src/Applications/soil-farm/SUMMS/summsPara)
  include(add_betr_library)
  add_betr_library(summs1layer ${SUMMS1LAYER_SOURCES})

  set(BETR_LIBRARIES summs1layer;${BETR_LIBRARIES} PARENT_SCOPE)
  set(BETR_LIBRARIES summs1layer;${BETR_LIBRARIES})

#X#add_subdirectory(tests)

  if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
    install(TARGETS summs1layer DESTINATION lib)
    file(GLOB HEADERS *.h)
    install(FILES ${HEADERS} DESTINATION include/soil-farm/SUMMS/summs1layer)
  endif()
>>>>>>> rzacplsbetr_cmupdated
