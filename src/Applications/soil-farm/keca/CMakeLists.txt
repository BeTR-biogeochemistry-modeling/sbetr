add_subdirectory(kecaPara)
add_subdirectory(keca1layer)
add_subdirectory(kecaNlayer)

set(BETR_LIBRARIES kecaPara keca1layer kecaNlayer;${BETR_LIBRARIES} PARENT_SCOPE)
