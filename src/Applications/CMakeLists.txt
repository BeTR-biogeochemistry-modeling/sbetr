set(BETR_APP_SOURCES
   ApplicationsFactory.F90
  )

include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_util)
include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_core)


include(add_betr_library)
add_betr_library(Applications ${BETR_APP_SOURCES})


set(BETR_LIBRARIES Applications;${BETR_LIBRARIES} PARENT_SCOPE)
set(BETR_LIBRARIES Applications;${BETR_LIBRARIES})

link_libraries(${BETR_LIBRARIES};${BETR_TPLS})

#X#add_subdirectory(tests)

if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS Applications DESTINATION lib)
  file(GLOB HEADERS *.h)
  install(FILES ${HEADERS} DESTINATION include/betr/Applications)
endif()
