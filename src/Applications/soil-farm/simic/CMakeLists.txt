add_subdirectory(simpara)
add_subdirectory(sim1layer)
add_subdirectory(simNlayer)

set(BETR_LIBRARIES simpara sim1layer simNlayer;${BETR_LIBRARIES} PARENT_SCOPE)
