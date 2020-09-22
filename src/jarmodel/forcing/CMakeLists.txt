set(BETR_JARFORC_SOURCES
  analforc.F90
  SoilForcType.F90
  AtmForcType.F90
  OMForcType.F90
  NutForcType.F90
  ForcDataType.F90
  SetJarForcMod.F90
  )

include_directories(${CMAKE_BINARY_DIR}/src/shr)
include_directories(${CMAKE_BINARY_DIR}/src/io_util)
include_directories(${CMAKE_BINARY_DIR}/src/stub_clm)
include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_util)
include_directories(${CMAKE_BINARY_DIR}/src/Applications/soil-farm/bgcfarm_util)
#begin_appadd
include_directories(${CMAKE_BINARY_DIR}/src/Applications/soil-farm/ecacnp/ecacnpPara)
include_directories(${CMAKE_BINARY_DIR}/src/Applications/soil-farm/keca/kecaPara)
include_directories(${CMAKE_BINARY_DIR}/src/Applications/soil-farm/v1eca/v1ecaPara)
#end_appadd
include(add_betr_library)
add_betr_library(betr_jar_forc ${BETR_JARFORC_SOURCES})

# disable all compiler warnings for clm code
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -w")

set(BETR_LIBRARIES betr_jar_forc;${BETR_LIBRARIES} PARENT_SCOPE)
set(BETR_LIBRARIES betr_jar_forc;${BETR_LIBRARIES})

#X#add_subdirectory(tests)

if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS betr_jar_forc DESTINATION lib)
  file(GLOB HEADERS *.h)
  install(FILES ${HEADERS} DESTINATION include/jarmodel/forc)
endif()
