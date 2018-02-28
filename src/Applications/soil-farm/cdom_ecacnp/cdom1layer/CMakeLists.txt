set(CDOM1LAYER_SOURCES
  BgccdomCnpIndexType.F90
  BgccdomCnpCompetType.F90
  BgccdomCnpDecompType.F90
  BgccdomCnpNitDenType.F90
  BgccdomCnpSOMType.F90
  BgccdomCnpType.F90
)


  include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_util)
  include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_math)
  include_directories(${CMAKE_BINARY_DIR}/src/Applications/soil-farm/bgcfarm_util)
  include_directories(${CMAKE_BINARY_DIR}/src/Applications/soil-farm/cdom_ecacnp/cdomPara)
  include(add_betr_library)
  add_betr_library(cdom1layer ${CDOM1LAYER_SOURCES})

  set(BETR_LIBRARIES cdom1layer;${BETR_LIBRARIES} PARENT_SCOPE)
  set(BETR_LIBRARIES cdom1layer;${BETR_LIBRARIES})

#X#add_subdirectory(tests)

  if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
    install(TARGETS cdom1layer DESTINATION lib)
    file(GLOB HEADERS *.h)
    install(FILES ${HEADERS} DESTINATION include/soil-farm/cdom_ecacnp/cdom1layer)
  endif()
