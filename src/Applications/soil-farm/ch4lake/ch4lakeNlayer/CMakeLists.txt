set(CH4LAKENLAYER_SOURCES
  ch4lakeBGCReactionsType.F90
  ch4lakePlantSoilBGCType.F90
)

  include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_util)
  include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_math)
  include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_bgc)
  include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_grid)
  include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_dtype)
  include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_core)
  include_directories(${CMAKE_BINARY_DIR}/src/Applications/soil-farm/bgcfarm_util)
  include_directories(${CMAKE_BINARY_DIR}/src/Applications/soil-farm/ch4lake/ch4lakePara)
  include_directories(${CMAKE_BINARY_DIR}/src/Applications/soil-farm/ch4lake/ch4lake1layer)
  include(add_betr_library)
  add_betr_library(ch4lakeNlayer ${CH4LAKENLAYER_SOURCES})

  set(BETR_LIBRARIES ch4lakeNlayer;${BETR_LIBRARIES} PARENT_SCOPE)
  set(BETR_LIBRARIES ch4lakeNlayer;${BETR_LIBRARIES})

#X#add_subdirectory(tests)

  if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
    install(TARGETS ch4lakeNlayer DESTINATION lib)
    file(GLOB HEADERS *.h)
    install(FILES ${HEADERS} DESTINATION include/soil-farm/ch4lake/ch4lakeNlayer)
  endif()
