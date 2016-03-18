set(BETR_UTILS_TESTS
  test_betr_utils.pfunit)

configure_file(testSuites.inc testSuites.inc COPYONLY)

set(BETR_UTILS_INCLUDE_DIRS
  src/betr/betr_util
  )

include(add_betr_pfunit_test)
add_betr_pfunit_test(test_betr_utils ${BETR_UTILS_TESTS} ${BETR_UTILS_INCLUDE_DIRS})
