set(V1ECANLAYER_SOURCES
  v1ecaBGCReactionsType.F90
  v1ecaPlantSoilBGCType.F90
)

  include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_util)
  include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_math)
  include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_bgc)
  include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_grid)
  include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_dtype)
  include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_core)
  include_directories(${CMAKE_BINARY_DIR}/src/Applications/soil-farm/bgcfarm_util)
  include_directories(${CMAKE_BINARY_DIR}/src/Applications/soil-farm/v1eca/v1ecaPara)
  include_directories(${CMAKE_BINARY_DIR}/src/Applications/soil-farm/v1eca/v1eca1layer)
  include(add_betr_library)
  add_betr_library(v1ecaNlayer ${V1ECANLAYER_SOURCES})

  set(BETR_LIBRARIES v1ecaNlayer;${BETR_LIBRARIES} PARENT_SCOPE)
  set(BETR_LIBRARIES v1ecaNlayer;${BETR_LIBRARIES})

#X#add_subdirectory(tests)

  if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
    install(TARGETS v1ecaNlayer DESTINATION lib)
    file(GLOB HEADERS *.h)
    install(FILES ${HEADERS} DESTINATION include/soil-farm/v1eca/v1ecaNlayer)
  endif()
