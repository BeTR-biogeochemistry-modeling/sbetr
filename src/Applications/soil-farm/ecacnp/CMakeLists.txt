add_subdirectory(ecacnpPara)
add_subdirectory(ecacnp1layer)
add_subdirectory(ecacnpNlayer)

set(BETR_LIBRARIES ecacnpPara ecacnp1layer ecacnpNlayer;${BETR_LIBRARIES} PARENT_SCOPE)
