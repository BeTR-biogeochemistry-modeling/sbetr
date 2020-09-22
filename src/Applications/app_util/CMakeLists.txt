set(BETR_APP_UTIL_SOURCES
   ApplicationsFactory.F90
  )


include_directories(${CMAKE_BINARY_DIR}/src/shr)
include_directories(${CMAKE_BINARY_DIR}/src/io_util)
include_directories(${CMAKE_BINARY_DIR}/src/stub_clm)
include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_util)
include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_core)
include_directories(${CMAKE_BINARY_DIR}/src/Applications/soil-farm/bgcfarm_util)
#begin_appadd
include_directories(${CMAKE_BINARY_DIR}/src/Applications/soil-farm/ecacnp/ecacnpPara)
include_directories(${CMAKE_BINARY_DIR}/src/Applications/soil-farm/ecacnp/ecacnp1layer)
include_directories(${CMAKE_BINARY_DIR}/src/Applications/soil-farm/ecacnp/ecacnpNlayer)
include_directories(${CMAKE_BINARY_DIR}/src/Applications/soil-farm/ch4soil/ch4soilPara)
include_directories(${CMAKE_BINARY_DIR}/src/Applications/soil-farm/ch4soil/ch4soil1layer)
include_directories(${CMAKE_BINARY_DIR}/src/Applications/soil-farm/ch4soil/ch4soilNlayer)
include_directories(${CMAKE_BINARY_DIR}/src/Applications/soil-farm/cdom/cdomPara)
include_directories(${CMAKE_BINARY_DIR}/src/Applications/soil-farm/cdom/cdom1layer)
include_directories(${CMAKE_BINARY_DIR}/src/Applications/soil-farm/cdom/cdomNlayer)
include_directories(${CMAKE_BINARY_DIR}/src/Applications/soil-farm/simic/simicPara)
include_directories(${CMAKE_BINARY_DIR}/src/Applications/soil-farm/simic/simic1layer)
include_directories(${CMAKE_BINARY_DIR}/src/Applications/soil-farm/simic/simicNlayer)
include_directories(${CMAKE_BINARY_DIR}/src/Applications/soil-farm/keca/kecaPara)
include_directories(${CMAKE_BINARY_DIR}/src/Applications/soil-farm/keca/keca1layer)
include_directories(${CMAKE_BINARY_DIR}/src/Applications/soil-farm/keca/kecaNlayer)
include_directories(${CMAKE_BINARY_DIR}/src/Applications/soil-farm/v1eca/v1ecaPara)
include_directories(${CMAKE_BINARY_DIR}/src/Applications/soil-farm/v1eca/v1eca1layer)
include_directories(${CMAKE_BINARY_DIR}/src/Applications/soil-farm/v1eca/v1ecaNlayer)
#end_appadd
include(add_betr_library)
add_betr_library(app_util ${BETR_APP_UTIL_SOURCES})


set(BETR_LIBRARIES app_util;${BETR_LIBRARIES} PARENT_SCOPE)
set(BETR_LIBRARIES app_util;${BETR_LIBRARIES})

link_libraries(${BETR_LIBRARIES};${BETR_TPLS})

#X#add_subdirectory(tests)


if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS app_util DESTINATION lib)
  file(GLOB HEADERS *.h)
  install(FILES ${HEADERS} DESTINATION include/betr/Applications/app_util)
endif()
