set(BETR_STANDALONE_COUPLER_SOURCES
  BGCReactionsFactoryMod.F90
  betr_standalone_cpl.F90
  )

include_directories(${CMAKE_BINARY_DIR}/src/shr)
include_directories(${CMAKE_BINARY_DIR}/src/stub_clm)
include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_math)
include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_core)
include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_grid)
include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_main)


include(add_betr_library)
add_betr_library(betr_standalone_coupler ${BETR_STANDALONE_COUPLER_SOURCES})

set(BETR_LIBRARIES betr_standalone_coupler;${BETR_LIBRARIES})

#X#add_subdirectory(tests)

if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS betr_standalone_coupler DESTINATION lib)
  file(GLOB HEADERS *.h)
  install(FILES ${HEADERS} DESTINATION include/betr/standalone_coupler)
endif()
