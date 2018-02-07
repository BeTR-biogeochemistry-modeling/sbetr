add_subdirectory(centPara)
add_subdirectory(cent1layer)
add_subdirectory(centNlayer)

set(BETR_LIBRARIES centPara cent1layer centNlayer;${BETR_LIBRARIES} PARENT_SCOPE)
