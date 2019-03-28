add_subdirectory(forcing)
add_subdirectory(test)
add_subdirectory(driver)

set(BETR_LIBRARIES forcing test driver;${BETR_LIBRARIES} PARENT_SCOPE)
