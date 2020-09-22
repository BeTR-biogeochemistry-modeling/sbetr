set(JARMODEL_TEST_RDPAR_SOURCES
  readparatest.F90
)
include_directories(${CMAKE_BINARY_DIR}/src/shr)
include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_util)
include_directories(${CMAKE_BINARY_DIR}/src/io_util)
include_directories(${CMAKE_BINARY_DIR}/src/jarmodel/forcing)
include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_math)
include_directories(${CMAKE_BINARY_DIR}/src/Applications/soil-farm/bgcfarm_util/)
include_directories(${CMAKE_BINARY_DIR}/src/Applications/soil-farm/ecacnp/ecacnpPara/)
include_directories(${CMAKE_BINARY_DIR}/src/Applications/soil-farm/cdom/cdomPara/)
include_directories(${CMAKE_BINARY_DIR}/src/Applications/soil-farm/keca/kecaPara/)


include(add_betr_library)
add_betr_library(jarmodel_test_rdpar ${JARMODEL_TEST_RDPAR_SOURCES})

set(BETR_LIBRARIES jarmodel_test_rdpar;${BETR_LIBRARIES} PARENT_SCOPE)
set(BETR_LIBRARIES jarmodel_test_rdpar;${BETR_LIBRARIES})

link_libraries(${BETR_LIBRARIES};${BETR_TPLS})
include(add_betr_executable)
add_betr_executable(readparatest readparatest.F90)

if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLLED")
  install(TARGETS readparatest DESTINATION bin)

  file(GLOB HEADS *.h)
  install(FILES ${HEADS} DESTINATION include/jarmodel/test/readpar)
endif()
