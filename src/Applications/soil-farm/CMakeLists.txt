add_subdirectory(bgcfarm_util)
add_subdirectory(ecacnp)
add_subdirectory(cdom)
add_subdirectory(simic)
# FIXME(bja, 201603) is this manual list necessary with the subdir
# adding to parent scope?

set(BETR_LIBRARIES bgcfarm_util;${BETR_LIBRARIES} PARENT_SCOPE)
