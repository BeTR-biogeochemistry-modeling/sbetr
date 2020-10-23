set(ECOSYS1LAYER_SOURCES
  ecosysBGCIndexType.F90
  ecosysBGCType.F90
)

include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_util)
include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_math)
include_directories(${CMAKE_BINARY_DIR}/src/Applications/soil-farm/bgcfarm_util)
include_directories(${CMAKE_BINARY_DIR}/src/Applications/soil-farm/ecosys/ecosysPara)
include(add_betr_library)
add_betr_library(ecosys1layer ${ECOSYS1LAYER_SOURCES})

set(BETR_LIBRARIES ecosys1layer;${BETR_LIBRARIES} PARENT_SCOPE)
set(BETR_LIBRARIES ecosys1layer;${BETR_LIBRARIES})

if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
   install(TARGETS ecosys1layer DESTINATION lib)
   file(GLOB HEADERS *.h)
   install(FILES ${HEADERS} DESTINATION include/soil-farm/ecosys/ecosys1layer)
endif()
