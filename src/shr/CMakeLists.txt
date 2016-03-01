set(BETR_SHR_SOURCES shr_const_mod.F90 shr_mpi_mod.F90
  shr_file_mod.F90 shr_nl_mod.F90
  shr_infnan_mod.F90  shr_string_mod.F90
  shr_kind_mod.F90 shr_sys_mod.F90
  shr_assert_mod.F90  shr_log_mod.F90 shr_timer_mod.F90
  )

configure_file(shr_assert.h shr_assert.h COPYONLY)
configure_file(mpif.h mpif.h COPYONLY)
configure_file(dtypes.h dtypes.h COPYONLY)

include(add_betr_library)
add_betr_library(betr_shr ${BETR_SHR_SOURCES})

set(BETR_LIBRARIES betr_shr;${BETR_LIBRARIES} PARENT_SCOPE)
set(BETR_LIBRARIES betr_shr;${BETR_LIBRARIES})

#X#add_subdirectory(tests)

if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS betr_math DESTINATION lib)
  file(GLOB HEADERS *.h)
  install(FILES ${HEADERS} DESTINATION include/betr/shr)
endif()
