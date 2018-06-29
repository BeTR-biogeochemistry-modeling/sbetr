set(SIMICPARA_SOURCES
 simicParaType.F90
)

include_directories(${CMAKE_BINARY_DIR}/src/shr)
include_directories(${CMAKE_BINARY_DIR}/src/io_util)
include_directories(${CMAKE_BINARY_DIR}/src/stub_clm)
include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_util)
include_directories(${CMAKE_BINARY_DIR}/src/Applications/soil-farm/bgcfarm_util)
include(add_betr_library)
add_betr_library(simicPara ${SIMICPARA_SOURCES})

set(BETR_LIBRARIES simicPara;${BETR_LIBRARIES} PARENT_SCOPE)
set(BETR_LIBRARIES simicPara;${BETR_LIBRARIES})

#X#add_subdirectory(tests)

if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS simicPara DESTINATION lib)
  file(GLOB HEADERS *.h)
  install(FILES ${HEADERS} DESTINATION include/soil-farm/simic/simicPara)
endif()
