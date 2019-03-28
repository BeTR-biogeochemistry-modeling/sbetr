set(CH4SOIL1LAYER_SOURCES
  ch4soilBGCIndexType.F90
  ch4soilBGCCompetType.F90
  ch4soilBGCDecompType.F90
  ch4soilBGCNitDenType.F90
  ch4soilBGCSOMType.F90
  ch4soilBGCType.F90
)


  include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_util)
  include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_math)
  include_directories(${CMAKE_BINARY_DIR}/src/Applications/soil-farm/bgcfarm_util)
  include_directories(${CMAKE_BINARY_DIR}/src/Applications/soil-farm/ch4soil/ch4soilPara)
  include(add_betr_library)
  add_betr_library(ch4soil1layer ${CH4SOIL1LAYER_SOURCES})

  set(BETR_LIBRARIES ch4soil1layer;${BETR_LIBRARIES} PARENT_SCOPE)
  set(BETR_LIBRARIES ch4soil1layer;${BETR_LIBRARIES})

#X#add_subdirectory(tests)

  if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
    install(TARGETS ch4soil1layer DESTINATION lib)
    file(GLOB HEADERS *.h)
    install(FILES ${HEADERS} DESTINATION include/soil-farm/ch4soil/ch4soil1layer)
  endif()
