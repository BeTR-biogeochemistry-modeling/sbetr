set(V1ECA1LAYER_SOURCES
  v1ecaBGCIndexType.F90
  v1ecaBGCCompetType.F90
  v1ecaBGCDecompType.F90
  v1ecaBGCNitDenType.F90
  v1ecaBGCSOMType.F90
  v1ecaBGCType.F90
)


  include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_util)
  include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_math)
  include_directories(${CMAKE_BINARY_DIR}/src/Applications/soil-farm/bgcfarm_util)
  include_directories(${CMAKE_BINARY_DIR}/src/Applications/soil-farm/v1eca/v1ecaPara)
  include(add_betr_library)
  add_betr_library(v1eca1layer ${V1ECA1LAYER_SOURCES})

  set(BETR_LIBRARIES v1eca1layer;${BETR_LIBRARIES} PARENT_SCOPE)
  set(BETR_LIBRARIES v1eca1layer;${BETR_LIBRARIES})

#X#add_subdirectory(tests)

  if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
    install(TARGETS v1eca1layer DESTINATION lib)
    file(GLOB HEADERS *.h)
    install(FILES ${HEADERS} DESTINATION include/soil-farm/v1eca/v1eca1layer)
  endif()
