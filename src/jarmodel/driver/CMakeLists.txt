set(JARMODEL_DRIVER_SOURCES
  jarmodel.F90
)
include_directories(${CMAKE_BINARY_DIR}/src/shr)
include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_util)
include_directories(${CMAKE_BINARY_DIR}/src/io_util)
include_directories(${CMAKE_BINARY_DIR}/src/jarmodel/forcing)
include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_math)
include_directories(${CMAKE_BINARY_DIR}/src/Applications/soil-farm/bgcfarm_util)
include_directories(${CMAKE_BINARY_DIR}/src/Applications/soil-farm/CENT_ECACNP/centPara)
include_directories(${CMAKE_BINARY_DIR}/src/Applications/soil-farm/CENT_ECACNP/cent1layer)
include_directories(${CMAKE_BINARY_DIR}/src/jarmodel/forcing)

include(add_betr_library)
add_betr_library(jarmodel_driver ${JARMODEL_DRIVER_SOURCES})

set(BETR_LIBRARIES jarmodel_driver;${BETR_LIBRARIES} PARENT_SCOPE)
set(BETR_LIBRARIES jarmodel_driver;${BETR_LIBRARIES})

link_libraries(${BETR_LIBRARIES};${BETR_TPLS})
include(add_betr_executable)
add_betr_executable(jarmodel jarmodel.F90)

if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLLED")
  install(TARGETS jarmodel DESTINATION bin)
  file(GLOB HEADS *.h)
  install(FILES ${HEADS} DESTINATION include/jarmodel/driver)
endif()
