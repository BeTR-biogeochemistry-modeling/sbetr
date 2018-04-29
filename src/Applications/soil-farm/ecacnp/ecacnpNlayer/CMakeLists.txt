set(ECACNPNLAYER_SOURCES
  ecacnpBGCReactionsType.F90
  ecacnpPlantSoilBGCType.F90
)

  include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_util)
  include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_math)
  include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_bgc)
  include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_grid)
  include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_dtype)
  include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_core)
  include_directories(${CMAKE_BINARY_DIR}/src/Applications/soil-farm/bgcfarm_util)
  include_directories(${CMAKE_BINARY_DIR}/src/Applications/soil-farm/ecacnp/ecacnpPara)
  include_directories(${CMAKE_BINARY_DIR}/src/Applications/soil-farm/ecacnp/ecacnp1layer)
  include(add_betr_library)
  add_betr_library(ecacnpNlayer ${ECACNPNLAYER_SOURCES})

  set(BETR_LIBRARIES ecacnpNlayer;${BETR_LIBRARIES} PARENT_SCOPE)
  set(BETR_LIBRARIES ecacnpNlayer;${BETR_LIBRARIES})

#X#add_subdirectory(tests)

  if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
    install(TARGETS ecacnpNlayer DESTINATION lib)
    file(GLOB HEADERS *.h)
    install(FILES ${HEADERS} DESTINATION include/soil-farm/ecacnp/ecacnpNlayer)
  endif()
