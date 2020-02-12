add_subdirectory(bgcfarm_util)
#begin_appadd
add_subdirectory(ecacnp)
add_subdirectory(cdom)
add_subdirectory(simic)
add_subdirectory(keca)
add_subdirectory(ch4soil)
add_subdirectory(v1eca)
#end_appadd
set(BETR_LIBRARIES bgcfarm_util;${BETR_LIBRARIES} PARENT_SCOPE)
