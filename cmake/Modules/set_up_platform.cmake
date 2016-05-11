# This macro identifies compilers and third-party library needs 
# for particular hosts.
macro(set_up_platform)

  # Set library suffix based on whether we're building shared/static.
  # FIXME: We have to hack this together here, since CMAKE_SHARED_LIBRARY_SUFFIX
  # FIXME: isn't available before project() is called, which is when set_up_platform()
  # FIXME: is invoked. Gross.
  if (BUILD_SHARED_LIBS)
    if (APPLE)
      set(LIB_SUFFIX .dylib)
    elseif (WIN32)
      set(LIB_SUFFIX .dll)
    else()
      set(LIB_SUFFIX .so)
    endif()
  else()
    set(LIB_SUFFIX .a)
  endif()

  # Set defaults for the various third-party libraries. These defaults
  # are hardwired because the project can't have been defined before 
  # this macro is executed, and so PROJECT_BINARY_DIR is unavailable.
  set(Z_LIBRARY "${CMAKE_CURRENT_BINARY_DIR}/lib/libz.a")
  set(Z_INCLUDE_DIR "${CMAKE_CURRENT_BINARY_DIR}/include")
  get_filename_component(Z_LIBRARY_DIR ${Z_LIBRARY} DIRECTORY)
  set(HDF5_LIB_NAME hdf5)
  set(HDF5_HL_LIB_NAME hdf5_hl)
  if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(HDF5_LIB_NAME ${HDF5_LIB_NAME}_debug)
    set(HDF5_HL_LIB_NAME ${HDF5_HL_LIB_NAME}_debug)
  endif()
  set(HDF5_LIBRARY "${CMAKE_CURRENT_BINARY_DIR}/lib/lib${HDF5_LIB_NAME}${LIB_SUFFIX}")
  set(HDF5_HL_LIBRARY "${CMAKE_CURRENT_BINARY_DIR}/lib/lib${HDF5_HL_LIB_NAME}${LIB_SUFFIX}")
  set(HDF5_LIBRARIES ${HDF5_HL_LIB_NAME};${HDF5_LIB_NAME})
  set(HDF5_INCLUDE_DIR "${CMAKE_CURRENT_BINARY_DIR}/include")
  get_filename_component(HDF5_LIBRARY_DIR ${Z_LIBRARY} DIRECTORY)
#  set(SILO_LIBRARY "${CMAKE_CURRENT_BINARY_DIR}/lib/libsiloh5.a")
#  set(SILO_LIBRARIES siloh5)
  if (APPLE)
    set(NEED_LAPACK FALSE)
  else()
    set(NEED_LAPACK TRUE)
  endif()

  # Certain tools (e.g. patch) require TMPDIR to be defined. If it is not, 
  # we do so here.
  set(TMPDIR_VAR $ENV{TMPDIR})
  if (NOT TMPDIR_VAR)
    # FIXME: Does this exist everywhere?
    set(ENV{TMPDIR} "/tmp")
  endif()

  # Get the hostname for this machine. 
  site_name(HOSTNAME)

  if (HOSTNAME MATCHES "cori") # NERSC Cori phase1
    #ndk : I had to "module load cray-hdf5 silo" 
    #ndk  make config debug=1 mpi=1 prefix=$SCRATCH/polymec
    # (Intel's compilers don't do C11.).
    set(CMAKE_C_COMPILER $ENV{CC})
    set(CMAKE_CXX_COMPILER $ENV{CXX})
    set(CMAKE_Fortran_COMPILER $ENV{FC})

   #    set(CMAKE_ADJUSTED_C_COMPILER cc -fPIE)
   #    set(CMAKE_ADJUSTED_C_COMPILER_ID Intel)
   #    set(CMAKE_ADJUSTED_CXX_COMPILER CC -fPIE)
   #    set(CMAKE_ADJUSTED_CXX_COMPILER_ID Intel)

    # We are cared for mathematically.
    set(NEED_LAPACK TRUE)

    # We expect the following libraries to be available.
#    set(Z_LIBRARY /usr/lib64/libz.a)
#    set(Z_INCLUDE_DIR /usr/include)
#    get_filename_component(Z_LIBRARY_DIR ${Z_LIBRARY} DIRECTORY)

    # Note that we use the hdf5 module and not cray-hdf5, since the silo 
    # module (below) is linked against hdf5 and not cray-hdf5.
    # FIXME: Use hdf5-parallel for parallel builds.
    # Set up HDF5.
#    if (USE_MPI EQUAL 1)
#      set(HDF5_LOC $ENV{HDF5_DIR})
#      if (NOT HDF5_LOC)
#        message(FATAL_ERROR "HDF5_DIR not found. Please load the hdf5-parallel module.")
#      endif()
#      include_directories(${HDF5_LOC}/include)
#      link_directories(${HDF5_LOC}/lib)
#      set(HDF5_LIBRARY ${HDF5_LOC}/lib/libhdf5_parallel.a)
#      set(HDF5_LIBRARIES hdf5_parallel;hdf5_hl_parallel)
#    else()
#      set(HDF5_LOC $ENV{HDF5_DIR})
#      if (NOT HDF5_LOC)
#        message(FATAL_ERROR "HDF5_DIR not found. Please load the hdf5 module.")
#      endif()
#      include_directories(${HDF5_LOC}/include)
#      list(APPEND EXTRA_LINK_DIRECTORIES ${HDF5_LOC}/lib)
#      link_directories(${HDF5_LOC}/lib)
#      set(HDF5_LIBRARY ${HDF5_LOC}/lib/libhdf5.a)
#      get_filename_component(HDF5_LIBRARY_DIR ${HDF5_LIBRARY} DIRECTORY)
#    endif()

#    set(SILO_LOC $ENV{SILO_DIR})
#    if (NOT SILO_LOC)
#      message(FATAL_ERROR "SILO_DIR not found. Please load the silo module.")
#    endif()

#    if (EXISTS ${SILO_LOC}/lib/libsiloh5.a)
#      include_directories(${SILO_LOC}/include)
#      link_directories(${SILO_LOC}/lib)
#      list(APPEND EXTRA_LINK_DIRECTORIES ${SILO_LOC}/lib)
#      set(SILO_LIBRARY ${SILO_LOC}/lib/libsiloh5.a)
#      set(SILO_LIBRARIES siloh5)
#    endif()

  elseif (HOSTNAME MATCHES "edison") # NERSC Edison
    # Edison likes Intel's compilers
    # (but Intel's compilers don't do C11.).
    set(CMAKE_C_COMPILER $ENV{CC})
    set(CMAKE_CXX_COMPILER $ENV{CXX})
    set(CMAKE_Fortran_COMPILER $ENV{FC})

#    set(CMAKE_ADJUSTED_C_COMPILER cc -fPIE)
#    set(CMAKE_ADJUSTED_C_COMPILER_ID Intel)
#    set(CMAKE_ADJUSTED_CXX_COMPILER CC -fPIE)
#    set(CMAKE_ADJUSTED_CXX_COMPILER_ID Intel)

    # We are cared for mathematically.
    set(NEED_LAPACK TRUE)

    # We expect the following libraries to be available.
    set(Z_LIBRARY /usr/lib64/libz.a)
    set(Z_INCLUDE_DIR /usr/include)
    get_filename_component(Z_LIBRARY_DIR ${Z_LIBRARY} DIRECTORY)

    # Note that we use the hdf5 module and not cray-hdf5, since the silo 
    # module (below) is linked against hdf5 and not cray-hdf5.
    # FIXME: Use hdf5-parallel for parallel builds.
    # Set up HDF5.
    if (USE_MPI EQUAL 1)
      set(HDF5_LOC $ENV{HDF5_DIR})
      if (NOT HDF5_LOC)
        message(FATAL_ERROR "HDF5_DIR not found. Please load the hdf5-parallel module.")
      endif()
      include_directories(${HDF5_LOC}/include)
      link_directories(${HDF5_LOC}/lib)
      set(HDF5_LIBRARY ${HDF5_LOC}/lib/libhdf5_parallel.a)
      set(HDF5_LIBRARIES hdf5_parallel;hdf5_hl_parallel)
    else()
      set(HDF5_LOC $ENV{HDF5_DIR})
      if (NOT HDF5_LOC)
        message(FATAL_ERROR "HDF5_DIR not found. Please load the hdf5 module.")
      endif()
      include_directories(${HDF5_LOC}/include)
      list(APPEND EXTRA_LINK_DIRECTORIES ${HDF5_LOC}/lib)
      link_directories(${HDF5_LOC}/lib)
      set(HDF5_LIBRARY ${HDF5_LOC}/lib/libhdf5.a)
      get_filename_component(HDF5_LIBRARY_DIR ${HDF5_LIBRARY} DIRECTORY)
    endif()

#    set(SILO_LOC $ENV{SILO_DIR})
#    if (NOT SILO_LOC)
#      message(FATAL_ERROR "SILO_DIR not found. Please load the silo module.")
#    endif()

#    if (EXISTS ${SILO_LOC}/lib/libsiloh5.a)
#      include_directories(${SILO_LOC}/include)
#      link_directories(${SILO_LOC}/lib)
#      list(APPEND EXTRA_LINK_DIRECTORIES ${SILO_LOC}/lib)
#      set(SILO_LIBRARY ${SILO_LOC}/lib/libsiloh5.a)
#      set(SILO_LIBRARIES siloh5)
#    endif()

  elseif(HOSTNAME MATCHES "hopper") # NERSC Hopper

    # Hopper is being decommissioned soon (Dec 2015), so we aren't super 
    # concerned with anything aesthetic here.
    set(CMAKE_C_COMPILER cc)
    set(CMAKE_CXX_COMPILER CC)
    set(CMAKE_Fortran_COMPILER ftn)

    set(NEED_LAPACK FALSE)

    # We expect the following libraries to be available.
    set(Z_LIBRARY /usr/lib64/libz.a)

    # Silo on Hopper doesn't use HDF5, so never mind that stuff.

#    set(SILO_LOC $ENV{SILO_DIR})
#    if (NOT SILO_LOC)
#      message(FATAL_ERROR "SILO_DIR not found. Please load the silo module.")
#    endif()
#    include_directories(${SILO_LOC}/include)
#    link_directories(${SILO_LOC}/lib)
#    list(APPEND EXTRA_LINK_DIRECTORIES ${SILO_LOC}/lib)
#    link_directories(/opt/pgi/default/linux86-64/default/lib) # built with PGI
#    list(APPEND EXTRA_LINK_DIRECTORIES /opt/pgi/default/linux86-64/default/lib) # built with PGI
#    set(SILO_LIBRARY ${SILO_LOC}/lib/libsilo.a)
#    set(SILO_LIBRARIES silo;pgc)

  endif()

endmacro()
