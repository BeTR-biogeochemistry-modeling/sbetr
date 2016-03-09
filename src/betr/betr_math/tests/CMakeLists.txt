set(BETR_MATH_TESTS
  test_betr_mathfunc.pfunit)

configure_file(testSuites.inc testSuites.inc COPYONLY)

set(BETR_MATH_INCLUDE_DIRS
  src/shr
  src/betr/betr_math
  )

include(add_betr_pfunit_test)
add_betr_pfunit_test(test_betr_mathfunc ${BETR_MATH_TESTS} ${BETR_MATH_INCLUDE_DIRS})

