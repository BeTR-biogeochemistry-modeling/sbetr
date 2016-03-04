set(BETR_DRIVER_SOURCES
  BeTRSimulation.F90
  BeTRSimulationFactory.F90
  BeTRSimulationStandalone.F90
  ReactionsFactoryStandalone.F90
  BGCReactionsMockRunType.F90
  PlantSoilBGCMockRunType.F90
  CLMForcType.F90
  sbetrDriverMod.F90
  )

include_directories(${CMAKE_BINARY_DIR}/src/shr)
include_directories(${CMAKE_BINARY_DIR}/src/stub_clm)
include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_util)
include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_math)
include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_core)
include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_main)


include(add_betr_library)
add_betr_library(betr_driver ${BETR_DRIVER_SOURCES})


set(BETR_LIBRARIES betr_driver;${BETR_LIBRARIES} PARENT_SCOPE)
set(BETR_LIBRARIES betr_driver;${BETR_LIBRARIES})

link_libraries(${BETR_LIBRARIES};${BETR_TPLS})
include(add_betr_executable)
add_betr_executable(sbetr sbetr.F90)

#X#add_subdirectory(tests)

if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS betr_driver DESTINATION lib)
  file(GLOB HEADERS *.h)
  install(FILES ${HEADERS} DESTINATION include/betr/driver)
endif()
