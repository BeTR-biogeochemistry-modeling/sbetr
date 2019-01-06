set(JARMODEL_DRIVER_SOURCES
  JarModelFactory.F90
  jarmodel.F90
)
include_directories(${CMAKE_BINARY_DIR}/src/shr)
include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_util)
include_directories(${CMAKE_BINARY_DIR}/src/io_util)
include_directories(${CMAKE_BINARY_DIR}/src/jarmodel/forcing)
include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_math)
include_directories(${CMAKE_BINARY_DIR}/src/Applications/soil-farm/bgcfarm_util)
#begin_appadd
include_directories(${CMAKE_BINARY_DIR}/src/Applications/soil-farm/ch4soil/ch4soilPara)
include_directories(${CMAKE_BINARY_DIR}/src/Applications/soil-farm/ch4soil/ch4soil1layer)
include_directories(${CMAKE_BINARY_DIR}/src/Applications/soil-farm/ecacnp/ecacnpPara)
include_directories(${CMAKE_BINARY_DIR}/src/Applications/soil-farm/ecacnp/ecacnp1layer)
include_directories(${CMAKE_BINARY_DIR}/src/Applications/soil-farm/cdom/cdomPara)
include_directories(${CMAKE_BINARY_DIR}/src/Applications/soil-farm/cdom/cdom1layer)
include_directories(${CMAKE_BINARY_DIR}/src/Applications/soil-farm/simic/simicPara)
include_directories(${CMAKE_BINARY_DIR}/src/Applications/soil-farm/simic/simic1layer)
include_directories(${CMAKE_BINARY_DIR}/src/Applications/soil-farm/keca/kecaPara)
include_directories(${CMAKE_BINARY_DIR}/src/Applications/soil-farm/keca/keca1layer)
include_directories(${CMAKE_BINARY_DIR}/src/Applications/soil-farm/v1eca/v1ecaPara)
include_directories(${CMAKE_BINARY_DIR}/src/Applications/soil-farm/v1eca/v1eca1layer)
#end_appadd
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
