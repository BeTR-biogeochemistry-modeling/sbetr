set(FAKENLAYER_SOURCES
  fakeBGCReactionsType.F90
  fakePlantSoilBGCType.F90
)

include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_util)
include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_math)
include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_bgc)
include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_grid)
include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_dtype)
include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_core)
include_directories(${CMAKE_BINARY_DIR}/src/Applications/soil-farm/bgcfarm_util)
include_directories(${CMAKE_BINARY_DIR}/src/Applications/soil-farm/fake/fakePara)
include_directories(${CMAKE_BINARY_DIR}/src/Applications/soil-farm/fake/fake1layer)
include(add_betr_library)
add_betr_library(fakeNlayer ${FAKENLAYER_SOURCES})
set(BETR_LIBRARIES fakeNlayer;${BETR_LIBRARIES} PARENT_SCOPE)
set(BETR_LIBRARIES fakeNlayer;${BETR_LIBRARIES})
if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS fakeNlayer DESTINATION lib)
  file(GLOB HEADERS *.h)
  install(FILES ${HEADERS} DESTINATION include/soil-farm/fake/fakeNlayer)
endif()
