set(SIMICNLAYER_SOURCES
  simicBGCReactionsType.F90
  simicPlantSoilBGCType.F90
)

  include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_util)
  include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_math)
  include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_bgc)
  include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_grid)
  include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_dtype)
  include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_core)
  include_directories(${CMAKE_BINARY_DIR}/src/Applications/soil-farm/bgcfarm_util)
  include_directories(${CMAKE_BINARY_DIR}/src/Applications/soil-farm/simic/simicPara)
  include_directories(${CMAKE_BINARY_DIR}/src/Applications/soil-farm/simic/simic1layer)
  include(add_betr_library)
  add_betr_library(simicNlayer ${SIMICNLAYER_SOURCES})

  set(BETR_LIBRARIES simicNlayer;${BETR_LIBRARIES} PARENT_SCOPE)
  set(BETR_LIBRARIES simicNlayer;${BETR_LIBRARIES})

#X#add_subdirectory(tests)

  if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
    install(TARGETS simicNlayer DESTINATION lib)
    file(GLOB HEADERS *.h)
    install(FILES ${HEADERS} DESTINATION include/soil-farm/simic/simicNlayer)
  endif()
