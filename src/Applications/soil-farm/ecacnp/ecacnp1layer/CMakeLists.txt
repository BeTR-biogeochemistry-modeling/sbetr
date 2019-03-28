set(ECACNP1LAYER_SOURCES
  ecacnpBGCIndexType.F90
  ecacnpBGCCompetType.F90
  ecacnpBGCDecompType.F90
  ecacnpBGCNitDenType.F90
  ecacnpBGCSOMType.F90
  ecacnpBGCType.F90
)


  include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_util)
  include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_math)
  include_directories(${CMAKE_BINARY_DIR}/src/Applications/soil-farm/bgcfarm_util)
  include_directories(${CMAKE_BINARY_DIR}/src/Applications/soil-farm/ecacnp/ecacnpPara)
  include(add_betr_library)
  add_betr_library(ecacnp1layer ${ECACNP1LAYER_SOURCES})

  set(BETR_LIBRARIES ecacnp1layer;${BETR_LIBRARIES} PARENT_SCOPE)
  set(BETR_LIBRARIES ecacnp1layer;${BETR_LIBRARIES})

#X#add_subdirectory(tests)

  if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
    install(TARGETS ecacnp1layer DESTINATION lib)
    file(GLOB HEADERS *.h)
    install(FILES ${HEADERS} DESTINATION include/soil-farm/ecacnp/ecacnp1layer)
  endif()
