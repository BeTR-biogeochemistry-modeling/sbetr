include_directories(${PROJECT_BINARY_DIR}/include)

# Build pfunit, a fortran unit testing library.
if (NOT EXISTS ${PROJECT_BINARY_DIR}/lib/libpfunit.a)
  message("Configuring pfunit...")
  execute_process(COMMAND cmake -E make_directory ${CMAKE_CURRENT_BINARY_DIR}/pfunit
    RESULT_VARIABLE stat)

  if (HAVE_MPI EQUAL 1)
    message(FATAL_ERROR "BeTR pUFnit MPI build is not yet supported.")

    set(PFUNIT_CMAKE_OPTS ${PFUNIT_CMAKE_OPTS} -DMPI=YES -DMPI_MPIFC=${CMAKE_Fortran_COMPILER})
  else()
    set(PFUNIT_CMAKE_OPTS ${PFUNIT_CMAKE_OPTS} -DMPI=NO)
  endif()
  
  set(PFUNIT_CMAKE_OPTS
    -DCMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR}
    -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
    -DCMAKE_Fortran_FLAGS=${CMAKE_Fortran_FLAGS}
    -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
    )
  execute_process(COMMAND cmake
    ${PFUNIT_CMAKE_OPTS}
    ${CMAKE_CURRENT_SOURCE_DIR}/pfunit
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/pfunit
    OUTPUT_VARIABLE pfunit_config.log ERROR_VARIABLE pfunit_config_err
    RESULT_VARIABLE stat)

  file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/pfunit_config.log ${pfunit_config_log})

  if (NOT stat EQUAL 0)
    file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/pfunit_config_errors.log ${pfunit_config_err})
    message(FATAL_ERROR "Configuration of pfunit library failed.")
  endif()

  message("Building and installing pfunit...")
  execute_process(COMMAND make install -j4
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/pfunit
                  OUTPUT_VARIABLE pfunit_build_log ERROR_VARIABLE pfunit_build_err
                  RESULT_VARIABLE stat)
  file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/pfunit_build.log ${pfunit_build_log})
  if (NOT stat EQUAL 0)
    file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/pfunit_build_errors.log ${pfunit_build_err})
    message(FATAL_ERROR "Build of pfunit library failed.")
  endif()

endif()
add_library(pfunit STATIC IMPORTED GLOBAL)
set_target_properties(pfunit PROPERTIES IMPORTED_LOCATION ${PROJECT_BINARY_DIR}/lib/libpfunit.a)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/3rd-party/pfunit/include/" PARENT_SCOPE)
set(ENV{PFUNIT} ${PROJECT_BINARY_DIR} CACHE)



# Add all the libraries to the build system at large.
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${BETR_TP_C_FLAGS}" PARENT_SCOPE)
set(BETR_LIBRARIES ${BETR_TP_LIBS};${BETR_LIBRARIES} PARENT_SCOPE)
set(BETR_INCDIRS ${BETR_INCDIRS};${BETR_TP_INCDIRS} PARENT_SCOPE)

