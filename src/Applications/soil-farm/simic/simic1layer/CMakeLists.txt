set(SIMIC1LAYER_SOURCES
  simicBGCIndexType.F90
  simicBGCType.F90
)


  include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_util)
  include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_math)
  include_directories(${CMAKE_BINARY_DIR}/src/Applications/soil-farm/bgcfarm_util)
  include_directories(${CMAKE_BINARY_DIR}/src/Applications/soil-farm/simic/simicPara)
  include(add_betr_library)
  add_betr_library(simic1layer ${SIMIC1LAYER_SOURCES})

  set(BETR_LIBRARIES simic1layer;${BETR_LIBRARIES} PARENT_SCOPE)
  set(BETR_LIBRARIES simic1layer;${BETR_LIBRARIES})

#X#add_subdirectory(tests)

  if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
    install(TARGETS simic1layer DESTINATION lib)
    file(GLOB HEADERS *.h)
    install(FILES ${HEADERS} DESTINATION include/soil-farm/simic/simic1layer)
  endif()
