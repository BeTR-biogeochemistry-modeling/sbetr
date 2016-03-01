set(BETR_MAIN_SOURCES
  BGCReactionsMod.F90
  BetrBGCMod.F90
  PlantSoilBGCMod.F90
  TracerBalanceMod.F90
  )

include_directories(${CMAKE_BINARY_DIR}/src/shr)
include_directories(${CMAKE_BINARY_DIR}/src/stub_clm)
include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_math)
include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_core)


include(add_betr_library)
add_betr_library(betr_main ${BETR_MAIN_SOURCES})

set(BETR_LIBRARIES betr_main;${BETR_LIBRARIES} PARENT_SCOPE)
set(BETR_LIBRARIES betr_main;${BETR_LIBRARIES})

#X#add_subdirectory(tests)

if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS betr_main DESTINATION lib)
  file(GLOB HEADERS *.h)
  install(FILES ${HEADERS} DESTINATION include/betr/main)
endif()
