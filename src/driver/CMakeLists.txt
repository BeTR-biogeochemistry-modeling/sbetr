set(BETR_DRIVER_SOURCES
  BeTRType.F90
  bhistFileMod.F90
  BeTRStatusSimType.F90
  BeTRSimulation.F90
  BeTRSimulationFactory.F90
  BeTRSimulationStandalone.F90
  BeTRSimulationCLM.F90
  BeTRSimulationALM.F90
  ForcingDataType.F90
  GridMod.F90
  sbetrDriverMod.F90
  )

include_directories(${CMAKE_BINARY_DIR}/src/shr)
include_directories(${CMAKE_BINARY_DIR}/src/stub_clm)
include_directories(${CMAKE_BINARY_DIR}/src/Applications)
include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_util)
include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_math)
include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_dtype)
include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_core)
include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_main)
include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_grid)
include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_rxns)

include(add_betr_library)
add_betr_library(betr_driver ${BETR_DRIVER_SOURCES})

#set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -w")

set(BETR_LIBRARIES betr_driver;${BETR_LIBRARIES} PARENT_SCOPE)
set(BETR_LIBRARIES betr_driver;${BETR_LIBRARIES})

link_libraries(${BETR_LIBRARIES};${BETR_TPLS})
include(add_betr_executable)
add_betr_executable(sbetr sbetr.F90)

#X#add_subdirectory(tests)

if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS sbetr DESTINATION bin)
  file(GLOB HEADERS *.mod)
  install(FILES ${HEADERS} DESTINATION include/betr/driver)
endif()
