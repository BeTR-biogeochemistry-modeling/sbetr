set(CENT1LAYER_SOURCES
  BgcCentCnpIndexType.F90
  BgcCentCnpCompetType.F90
  BgcCentCnpDecompType.F90
  BgcCentCnpNitDenType.F90
  BgcCentCnpSOMType.F90
  BgcCentCnpType.F90
)


  include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_util)
  include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_math)
  include_directories(${CMAKE_BINARY_DIR}/src/Applications/soil-farm/bgcfarm_util)
  include_directories(${CMAKE_BINARY_DIR}/src/Applications/soil-farm/CENT_ECACNP/centPara)
  include(add_betr_library)
  add_betr_library(cent1layer ${CENT1LAYER_SOURCES})

  set(BETR_LIBRARIES cent1layer;${BETR_LIBRARIES} PARENT_SCOPE)
  set(BETR_LIBRARIES cent1layer;${BETR_LIBRARIES})

#X#add_subdirectory(tests)

  if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
    install(TARGETS cent1layer DESTINATION lib)
    file(GLOB HEADERS *.h)
    install(FILES ${HEADERS} DESTINATION include/soil-farm/CENT_ECACNP/cent1layer)
  endif()
