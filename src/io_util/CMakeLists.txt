set(BETR_IO_UTIL_SOURCES
  ncdio_pio.F90
  histMod.F90
  )

include_directories(${CMAKE_BINARY_DIR}/src/shr)
include_directories(${CMAKE_BINARY_DIR}/src/esmf_wrf_timemgr)
include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_util)
include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_math)

include(add_betr_library)
add_betr_library(betr_io_util ${BETR_IO_UTIL_SOURCES})

# disable all compiler warnings for clm code
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -w")

set(BETR_LIBRARIES betr_io_util;${BETR_LIBRARIES} PARENT_SCOPE)
set(BETR_LIBRARIES betr_io_util;${BETR_LIBRARIES})

#X#add_subdirectory(tests)

if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS betr_io_util DESTINATION lib)
  file(GLOB HEADERS *.h)
  install(FILES ${HEADERS} DESTINATION include/betr/io_util)
endif()
