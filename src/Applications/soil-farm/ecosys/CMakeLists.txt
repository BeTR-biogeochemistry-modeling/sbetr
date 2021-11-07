add_subdirectory(ecosysPara)
add_subdirectory(ecosys1layer)
add_subdirectory(ecosysNlayer)

set(BETR_LIBRARIES ecosysPara ecosys1layer ecosysNlayer;${BETR_LIBRARIES} PARENT_SCOPE)
