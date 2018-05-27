set(CDOM1LAYER_SOURCES
  cdomBGCIndexType.F90
  cdomBGCCompetType.F90
  cdomBGCDecompType.F90
  cdomBGCNitDenType.F90
  cdomBGCSOMType.F90
  cdomBGCType.F90
)


  include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_util)
  include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_math)
  include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_para)
  include_directories(${CMAKE_BINARY_DIR}/src/Applications/soil-farm/bgcfarm_util)
  include_directories(${CMAKE_BINARY_DIR}/src/Applications/soil-farm/cdom/cdomPara)
  include(add_betr_library)
  add_betr_library(cdom1layer ${CDOM1LAYER_SOURCES})

  set(BETR_LIBRARIES cdom1layer;${BETR_LIBRARIES} PARENT_SCOPE)
  set(BETR_LIBRARIES cdom1layer;${BETR_LIBRARIES})

#X#add_subdirectory(tests)

  if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
    install(TARGETS cdom1layer DESTINATION lib)
    file(GLOB HEADERS *.h)
    install(FILES ${HEADERS} DESTINATION include/soil-farm/cdom/cdom1layer)
  endif()
