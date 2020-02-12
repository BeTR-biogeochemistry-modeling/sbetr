add_subdirectory(v1ecaPara)
add_subdirectory(v1eca1layer)
add_subdirectory(v1ecaNlayer)

set(BETR_LIBRARIES v1ecaPara v1eca1layer v1ecaNlayer;${BETR_LIBRARIES} PARENT_SCOPE)
