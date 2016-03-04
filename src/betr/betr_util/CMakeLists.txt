set(BETR_UTIL_SOURCES
  betr_constants.F90
  )

include(add_betr_library)
add_betr_library(betr_util ${BETR_UTIL_SOURCES})

set(BETR_LIBRARIES betr_util;${BETR_LIBRARIES} PARENT_SCOPE)
set(BETR_LIBRARIES betr_util;${BETR_LIBRARIES})

#X#add_subdirectory(tests)

if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS betr_util DESTINATION lib)
  file(GLOB HEADERS *.h)
  install(FILES ${HEADERS} DESTINATION include/betr/util)
endif()
