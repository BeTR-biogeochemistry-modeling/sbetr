set(CENTPARA_SOURCES
 CentParaType.F90
)

include_directories(${CMAKE_BINARY_DIR}/src/shr)
include_directories(${CMAKE_BINARY_DIR}/src/io_util)
include_directories(${CMAKE_BINARY_DIR}/src/stub_clm)
include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_util)
include_directories(${CMAKE_BINARY_DIR}/src/Applications/soil-farm/bgcfarm_util)
include(add_betr_library)
add_betr_library(centPara ${CENTPARA_SOURCES})

set(BETR_LIBRARIES centPara;${BETR_LIBRARIES} PARENT_SCOPE)
set(BETR_LIBRARIES centPara;${BETR_LIBRARIES})

#X#add_subdirectory(tests)

if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS centPara DESTINATION lib)
  file(GLOB HEADERS *.h)
  install(FILES ${HEADERS} DESTINATION include/soil-farm/CENT_ECACNP/centPara)
endif()
