set(KECA1LAYER_SOURCES
  kecaBGCIndexType.F90
  kecaBGCCompetType.F90
  kecaBGCDecompType.F90
  kecaBGCNitDenType.F90
  kecaBGCSOMType.F90
  kecaBGCType.F90
)


  include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_util)
  include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_math)
  include_directories(${CMAKE_BINARY_DIR}/src/Applications/soil-farm/bgcfarm_util)
  include_directories(${CMAKE_BINARY_DIR}/src/Applications/soil-farm/keca/kecaPara)
  include(add_betr_library)
  add_betr_library(keca1layer ${KECA1LAYER_SOURCES})

  set(BETR_LIBRARIES keca1layer;${BETR_LIBRARIES} PARENT_SCOPE)
  set(BETR_LIBRARIES keca1layer;${BETR_LIBRARIES})

#X#add_subdirectory(tests)

  if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
    install(TARGETS keca1layer DESTINATION lib)
    file(GLOB HEADERS *.h)
    install(FILES ${HEADERS} DESTINATION include/soil-farm/keca/keca1layer)
  endif()
