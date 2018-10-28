set(FAKE1LAYER_SOURCES
  fakeBGCIndexType.F90
  fakeBGCType.F90
)

include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_util)
include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_math)
include_directories(${CMAKE_BINARY_DIR}/src/Applications/soil-farm/bgcfarm_util)
include_directories(${CMAKE_BINARY_DIR}/src/Applications/soil-farm/fake/fakePara)
include(add_betr_library)
add_betr_library(fake1layer ${FAKE1LAYER_SOURCES})

set(BETR_LIBRARIES fake1layer;${BETR_LIBRARIES} PARENT_SCOPE)
set(BETR_LIBRARIES fake1layer;${BETR_LIBRARIES})

if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
   install(TARGETS fake1layer DESTINATION lib)
   file(GLOB HEADERS *.h)
   install(FILES ${HEADERS} DESTINATION include/soil-farm/fake/fake1layer)
endif()
