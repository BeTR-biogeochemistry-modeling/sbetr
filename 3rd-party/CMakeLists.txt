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

# Build zlib for compression.
if (NOT EXISTS ${Z_LIBRARY})
  set(ZLIB_CONFIG_OPTS --prefix=${PROJECT_BINARY_DIR} --static)

  message("Preparing zlib...")
  execute_process(COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/zlib ${CMAKE_CURRENT_BINARY_DIR}/zlib
                  RESULT_VARIABLE stat)
  if (NOT stat EQUAL 0)
    message(FATAL_ERROR "Copying library source failed.")
  endif()

  message("Configuring zlib...")
  execute_process(COMMAND env CC=${CMAKE_C_COMPILER} CFLAGS=${CMAKE_C_FLAGS} ./configure ${ZLIB_CONFIG_OPTS}
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/zlib
                  OUTPUT_VARIABLE zlib_config_log ERROR_VARIABLE zlib_config_err
                  RESULT_VARIABLE stat)
  file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/zlib_config.log ${zlib_config_log})
  if (NOT stat EQUAL 0)
    if (TRAVIS_CI)
      message(WARNING "Configuration of zlib library failed:")
      message(FATAL_ERROR ${zlib_config_err})
    else()
      file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/zlib_config_errors.log ${zlib_config_err})
      message(FATAL_ERROR "Configuration of zlib library failed. See ${CMAKE_CURRENT_BINARY_DIR}/zlib_config_*.log.")
    endif()
  endif()

  message("Building and installing zlib (${NUM_BUILD_THREADS} threads)...")
  execute_process(COMMAND make install -j${NUM_BUILD_THREADS}
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/zlib
                  OUTPUT_VARIABLE zlib_build_log ERROR_VARIABLE zlib_build_err
                  RESULT_VARIABLE stat)
  file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/zlib_build.log ${zlib_build_log})
  if (NOT stat EQUAL 0)
    if (TRAVIS_CI)
      message(WARNING "Build of zlib library failed:")
      message(FATAL_ERROR ${zlib_build_err})
    else()
      file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/zlib_build_errors.log ${zlib_build_err})
      message(FATAL_ERROR "Build of zlib library failed. See ${CMAKE_CURRENT_BINARY_DIR}/zlib_build_*.log.")
    endif()
  endif()

endif()
set(BETR_TP_LIBS z;${BETR_TP_LIBS})

if (NOT EXISTS ${HDF5_LIBRARY})
  message("Configuring hdf5...")
  set(HDF5_CMAKE_OPTS -DCMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR} -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER})
  set(HDF5_CMAKE_OPTS ${HDF5_CMAKE_OPTS} -DBUILD_SHARED_LIBS=${BUILD_SHARED_LIBS} -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}) 
  set(HDF5_CMAKE_OPTS ${HDF5_CMAKE_OPTS} -DZLIB_INCLUDE_DIRS=${PROJECT_BINARY_DIR}/include -DZLIB_LIBARIES=z -DZLIB_USE_EXTERNAL=OFF -DHDF5_ENABLE_PARALLEL=${HAVE_MPI})
  set(HDF5_CMAKE_OPTS ${HDF5_CMAKE_OPTS} -DHDF5_BUILD_CPP_LIB=OFF)
  set(HDF5_CMAKE_OPTS ${HDF5_CMAKE_OPTS} -DHDF5_BUILD_FORTRAN=OFF)
  set(HDF5_CMAKE_OPTS ${HDF5_CMAKE_OPTS} -DBUILD_TESTING=OFF -DHDF5_BUILD_EXAMPLES=OFF)
  if (EXISTS ${CMAKE_CURRENT_BINARY_DIR}/hdf5)
    execute_process(COMMAND ${CMAKE_COMMAND} -E remove ${CMAKE_CURRENT_BINARY_DIR}/hdf5/CMakeCache.txt)
  else()
    execute_process(COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_CURRENT_BINARY_DIR}/hdf5)
  endif()
  execute_process(COMMAND ${CMAKE_COMMAND} ${CMAKE_CURRENT_SOURCE_DIR}/hdf5 ${HDF5_CMAKE_OPTS}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/hdf5
    OUTPUT_VARIABLE hdf5_config_log ERROR_VARIABLE hdf5_config_err
    RESULT_VARIABLE stat)
  file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/hdf5_config.log ${hdf5_config_log})
  if (NOT stat EQUAL 0)
    if (TRAVIS_CI)
      message(WARNING "Configuration of hdf5 library failed:")
      message(FATAL_ERROR ${hdf5_config_err})
    else()
      file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/hdf5_config_errors.log ${hdf5_config_err})
      message(FATAL_ERROR "Configuration of hdf5 library failed. See ${CMAKE_CURRENT_BINARY_DIR}/hdf5_config_*.log.")
    endif()
  endif()

  message("Building and installing hdf5 (${NUM_BUILD_THREADS} threads)...")
  execute_process(COMMAND make install -j${NUM_BUILD_THREADS}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/hdf5
    OUTPUT_VARIABLE hdf5_build_log ERROR_VARIABLE hdf5_build_err
    RESULT_VARIABLE stat)
  file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/hdf5_build.log ${hdf5_build_log})
  if (NOT stat EQUAL 0)
    if (TRAVIS_CI)
      message(WARNING "Build of hdf5 library failed:")
      message(FATAL_ERROR ${hdf5_build_err})
    else()
      file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/hdf5_build_errors.log ${hdf5_build_err})
      message(FATAL_ERROR "Build of hdf5 library failed. See ${CMAKE_CURRENT_BINARY_DIR}/hdf5_build_*.log.")
    endif()
  endif()
  
  # On Macs we need to change the rpath of shared libraries so we can install them elsewhere.
  if (APPLE AND BUILD_SHARED_LIBS)
    execute_process(COMMAND otool -DX ${PROJECT_BINARY_DIR}/lib/lib${HDF5_LIB_NAME}${LIB_SUFFIX} 
      OUTPUT_VARIABLE libhdf5_id)
    string(STRIP ${libhdf5_id} libhdf5_id)
    execute_process(COMMAND install_name_tool -id @rpath/lib${HDF5_LIB_NAME}${LIB_SUFFIX} ${PROJECT_BINARY_DIR}/lib/${HDF5_LIB_NAME}${LIB_SUFFIX})
    execute_process(COMMAND install_name_tool -id @rpath/lib${HDF5_HL_LIB_NAME}${LIB_SUFFIX} ${PROJECT_BINARY_DIR}/lib/lib${HDF5_HL_LIB_NAME}${LIB_SUFFIX})
    execute_process(COMMAND install_name_tool -change ${libhdf5_id} @rpath/libhdf5${LIB_SUFFIX} ${PROJECT_BINARY_DIR}/lib/${HDF5_HL_LIB_NAME}${LIB_SUFFIX})
  endif()
  
endif()
if (HDF5_LIBRARY STREQUAL "${PROJECT_BINARY_DIR}/lib/lib${HDF5_LIB_NAME}${LIB_SUFFIX}")
  set(POLYMEC_HDF5_LIBRARY "${CMAKE_INSTALL_PREFIX}/lib/lib${HDF5_LIB_NAME}${LIB_SUFFIX}" PARENT_SCOPE)
  set(POLYMEC_HDF5_HL_LIBRARY "${CMAKE_INSTALL_PREFIX}/lib/lib${HDF5_HL_LIB_NAME}${LIB_SUFFIX}" PARENT_SCOPE)
  set(POLYMEC_HDF5_INCLUDE_DIR "${CMAKE_INSTALL_PREFIX}/include" PARENT_SCOPE)
else()
  set(POLYMEC_HDF5_LIBRARY ${HDF5_LIBRARY} PARENT_SCOPE)
  set(POLYMEC_HDF5_HL_LIBRARY ${HDF5_HL_LIBRARY} PARENT_SCOPE)
  set(POLYMEC_HDF5_INCLUDE_DIR ${HDF5_INCLUDE_DIR} PARENT_SCOPE)
endif()

# If our HDF5 library is named something weird, copy it to a file with the 
# regular boring name. This prevents libraries like Silo and Exodus (which 
# take a while to catch up to recent developments) from freaking out.
if (NOT HDF5_LIB_NAME STREQUAL "hdf5" AND HDF5_LIBRARY MATCHES "${PROJECT_BINARY_DIR}/lib" AND NOT EXISTS ${PROJECT_BINARY_DIR}/lib/libhdf5${LIB_SUFFIX})
  execute_process(COMMAND ${CMAKE_COMMAND} -E copy ${HDF5_LIBRARY} ${PROJECT_BINARY_DIR}/lib/libhdf5${LIB_SUFFIX})
endif()

set(POLYMEC_TP_LIBS ${HDF5_LIBRARIES};${POLYMEC_TP_LIBS})


# Add all the libraries to the build system at large.
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${BETR_TP_C_FLAGS}" PARENT_SCOPE)
set(BETR_LIBRARIES ${BETR_TP_LIBS};${BETR_LIBRARIES} PARENT_SCOPE)
set(BETR_INCDIRS ${BETR_INCDIRS};${BETR_TP_INCDIRS} PARENT_SCOPE)

