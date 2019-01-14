set(SUMMSNLAYER_SOURCES
  BgcReactionsSummsType.F90
  PlantSoilBgcSummsType.F90
)

  include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_util)
  include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_math)
  include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_bgc)
  include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_grid)
  include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_dtype)
  include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_core)
  include_directories(${CMAKE_BINARY_DIR}/src/Applications/soil-farm/bgcfarm_util)
  include_directories(${CMAKE_BINARY_DIR}/src/Applications/soil-farm/SUMMS/summsPara)
  include_directories(${CMAKE_BINARY_DIR}/src/Applications/soil-farm/SUMMS/summs1layer)
  include(add_betr_library)
  add_betr_library(summsNlayer ${SUMMSNLAYER_SOURCES})

  set(BETR_LIBRARIES summsNlayer;${BETR_LIBRARIES} PARENT_SCOPE)
  set(BETR_LIBRARIES summsNlayer;${BETR_LIBRARIES})

#X#add_subdirectory(tests)

  if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
    install(TARGETS summsNlayer DESTINATION lib)
    file(GLOB HEADERS *.h)
    install(FILES ${HEADERS} DESTINATION include/soil-farm/SUMMS/summsNlayer)
  endif()
