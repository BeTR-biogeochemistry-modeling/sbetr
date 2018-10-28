add_subdirectory(bgcfarm_util)
#begin_appadd_directories
add_subdirectory(ecacnp)
add_subdirectory(cdom)
add_subdirectory(simic)
add_subdirectory(keca)
#end_appadd_directories
set(BETR_LIBRARIES bgcfarm_util;${BETR_LIBRARIES} PARENT_SCOPE)
