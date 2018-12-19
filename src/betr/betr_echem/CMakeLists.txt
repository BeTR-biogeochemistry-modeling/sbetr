set(BETR_ECHEM_SOURCES
  eChemConstMod.F90
  eChemMod.F90
  )

include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_util)
include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_math)
include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_dtype)

include(add_betr_library)
add_betr_library(betr_echem ${BETR_ECHEM_SOURCES})

set(BETR_LIBRARIES betr_echem;${BETR_LIBRARIES} PARENT_SCOPE)
set(BETR_LIBRARIES betr_echem;${BETR_LIBRARIES})

#X#add_subdirectory(tests)

if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS betr_echem DESTINATION lib)
  file(GLOB HEADERS *.h)
  install(FILES ${HEADERS} DESTINATION include/betr/echem)
endif()
