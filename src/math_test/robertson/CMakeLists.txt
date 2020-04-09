set(MATH_TEST_ROB_SOURCES
  robertsonMod.F90
  robertson.F90
  )

include_directories(${CMAKE_BINARY_DIR}/src/shr)
include_directories(${CMAKE_BINARY_DIR}/src/io_util)
include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_util)
include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_math)


include(add_betr_library)
add_betr_library(math_test_rob ${MATH_TEST_ROB_SOURCES})

#set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -w")

set(BETR_LIBRARIES math_test_rob;${BETR_LIBRARIES} PARENT_SCOPE)
set(BETR_LIBRARIES math_test_rob;${BETR_LIBRARIES})

link_libraries(${BETR_LIBRARIES};${BETR_TPLS})
include(add_betr_executable)
add_betr_executable(robertson robertson.F90)
#X#add_subdirectory(tests)

if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS robertson DESTINATION bin)
  file(GLOB HEADERS *.mod)
  install(FILES ${HEADERS} DESTINATION include/math_test/robertson)
endif()
