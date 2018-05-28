add_subdirectory(summsPara)
add_subdirectory(summs1layer)
add_subdirectory(summsNlayer)

set(BETR_LIBRARIES summsPara summs1layer summsNlayer;${BETR_LIBRARIES} PARENT_SCOPE)
