set(CH4SOILNLAYER_SOURCES
  ch4soilBGCReactionsType.F90
  ch4soilPlantSoilBGCType.F90
)

  include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_util)
  include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_math)
  include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_bgc)
  include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_grid)
  include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_dtype)
  include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_core)
  include_directories(${CMAKE_BINARY_DIR}/src/Applications/soil-farm/bgcfarm_util)
  include_directories(${CMAKE_BINARY_DIR}/src/Applications/soil-farm/ch4soil/ch4soilPara)
  include_directories(${CMAKE_BINARY_DIR}/src/Applications/soil-farm/ch4soil/ch4soil1layer)
  include(add_betr_library)
  add_betr_library(ch4soilNlayer ${CH4SOILNLAYER_SOURCES})

  set(BETR_LIBRARIES ch4soilNlayer;${BETR_LIBRARIES} PARENT_SCOPE)
  set(BETR_LIBRARIES ch4soilNlayer;${BETR_LIBRARIES})

#X#add_subdirectory(tests)

  if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
    install(TARGETS ch4soilNlayer DESTINATION lib)
    file(GLOB HEADERS *.h)
    install(FILES ${HEADERS} DESTINATION include/soil-farm/ch4soil/ch4soilNlayer)
  endif()
