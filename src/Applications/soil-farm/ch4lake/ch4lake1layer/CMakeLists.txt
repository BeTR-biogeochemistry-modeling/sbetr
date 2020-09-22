set(CH4LAKE1LAYER_SOURCES
  ch4lakeBGCIndexType.F90
  ch4lakeBGCCompetType.F90
  ch4lakeBGCDecompType.F90
  ch4lakeBGCNitDenType.F90
  ch4lakeBGCSOMType.F90
  ch4lakeBGCType.F90
)


  include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_util)
  include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_math)
  include_directories(${CMAKE_BINARY_DIR}/src/Applications/soil-farm/bgcfarm_util)
  include_directories(${CMAKE_BINARY_DIR}/src/Applications/soil-farm/ch4lake/ch4lakePara)
  include(add_betr_library)
  add_betr_library(ch4lake1layer ${CH4LAKE1LAYER_SOURCES})

  set(BETR_LIBRARIES ch4lake1layer;${BETR_LIBRARIES} PARENT_SCOPE)
  set(BETR_LIBRARIES ch4lake1layer;${BETR_LIBRARIES})

#X#add_subdirectory(tests)

  if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
    install(TARGETS ch4lake1layer DESTINATION lib)
    file(GLOB HEADERS *.h)
    install(FILES ${HEADERS} DESTINATION include/soil-farm/ch4lake/ch4lake1layer)
  endif()
