include(PreprocessMacro)
add_pfunit_sources(MATHFUNC_SRCS test_betr_mathfunc.pfunit)

configure_file(testSuites.inc testSuites.inc COPYONLY)

include(add_betr_pfunit_test)
add_betr_pfunit_test(test_betr_mathfunc ${MATHFUNC_SRCS})

