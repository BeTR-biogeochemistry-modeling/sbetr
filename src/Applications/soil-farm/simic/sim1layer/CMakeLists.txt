set(SIM1LAYER_SOURCES
  SimicBGCType.F90
)


  include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_util)
  include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_math)
  include_directories(${CMAKE_BINARY_DIR}/src/Applications/soil-farm/bgcfarm_util)
  include_directories(${CMAKE_BINARY_DIR}/src/Applications/soil-farm/simic/simpara)
  include(add_betr_library)
  add_betr_library(sim1layer ${SIM1LAYER_SOURCES})

  set(BETR_LIBRARIES sim1layer;${BETR_LIBRARIES} PARENT_SCOPE)
  set(BETR_LIBRARIES sim1layer;${BETR_LIBRARIES})

#X#add_subdirectory(tests)

  if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
    install(TARGETS sim1layer DESTINATION lib)
    file(GLOB HEADERS *.h)
    install(FILES ${HEADERS} DESTINATION include/soil-farm/simic/sim1layer)
  endif()
