set(BGCFARM_UTIL_SOURCES
 BiogeoConType.F90
 JarBgcForcType.F90
 BeTRJarModel.F90
 SOMStateVarUpdateMod.F90
 gBGCIndexType.F90
 GeoChemAlgorithmMod.F90
 EcosysMicDynParamMod.F90
)

include_directories(${CMAKE_BINARY_DIR}/src/shr)
include_directories(${CMAKE_BINARY_DIR}/src/io_util)
include_directories(${CMAKE_BINARY_DIR}/src/stub_clm)
include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_math)
include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_dtype)

include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_util)
include(add_betr_library)
add_betr_library(bgcfarm_util ${BGCFARM_UTIL_SOURCES})

set(BETR_LIBRARIES bgcfarm_util;${BETR_LIBRARIES} PARENT_SCOPE)
set(BETR_LIBRARIES bgcfarm_util;${BETR_LIBRARIES})

#X#add_subdirectory(tests)

if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS bgcfarm_util DESTINATION lib)
  file(GLOB HEADERS *.h)
  install(FILES ${HEADERS} DESTINATION include/soil-farm/bgcfarm_util)
endif()
