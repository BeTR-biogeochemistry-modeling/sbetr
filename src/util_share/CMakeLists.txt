set(BETR_UTIL_SHARE_SOURCES
  abortutils.F90
  accumulMod.F90
  clm_time_manager.F90
  debuginfoMod.F90
  decompMod.F90
  fileutils.F90
  getdatetime.F90
  ncdio_pio.F90
  spmdMod.F90
  )

include(add_betr_library)
add_betr_library(betr_util_share ${BETR_UTIL_SHARE_SOURCES})

set(BETR_LIBRARIES betr_util_share;${BETR_LIBRARIES})

#X#add_subdirectory(tests)

if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS betr_util_share DESTINATION lib)
  file(GLOB HEADERS *.h)
  install(FILES ${HEADERS} DESTINATION include/util_share)
endif()
