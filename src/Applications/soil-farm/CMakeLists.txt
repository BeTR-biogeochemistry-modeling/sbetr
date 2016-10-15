add_subdirectory(bgcfarm_util)
add_subdirectory(CENT_ECACNP)
# FIXME(bja, 201603) is this manual list necessary with the subdir
# adding to parent scope?

set(BETR_LIBRARIES bgcfarm_util;cent_ecacnp;${BETR_LIBRARIES} PARENT_SCOPE)
