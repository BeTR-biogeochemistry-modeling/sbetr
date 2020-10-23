set(ECOSYSNLAYER_SOURCES
  ecosysBGCReactionsType.F90
  ecosysPlantSoilBGCType.F90
)

include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_util)
include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_math)
include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_bgc)
include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_grid)
include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_dtype)
include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_core)
include_directories(${CMAKE_BINARY_DIR}/src/Applications/soil-farm/bgcfarm_util)
include_directories(${CMAKE_BINARY_DIR}/src/Applications/soil-farm/ecosys/ecosysPara)
include_directories(${CMAKE_BINARY_DIR}/src/Applications/soil-farm/ecosys/ecosys1layer)
include(add_betr_library)
add_betr_library(ecosysNlayer ${ECOSYSNLAYER_SOURCES})
set(BETR_LIBRARIES ecosysNlayer;${BETR_LIBRARIES} PARENT_SCOPE)
set(BETR_LIBRARIES ecosysNlayer;${BETR_LIBRARIES})
if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS ecosysNlayer DESTINATION lib)
  file(GLOB HEADERS *.h)
  install(FILES ${HEADERS} DESTINATION include/soil-farm/ecosys/ecosysNlayer)
endif()
