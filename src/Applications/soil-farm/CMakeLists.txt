add_subdirectory(bgcfarm_util)
add_subdirectory(CENT_ECACNP)

if (UGM)
  add_subdirectory(UGM-M)
endif()
#add_subdirectory(CENT_ECADCNP)
# FIXME(bja, 201603) is this manual list necessary with the subdir
# adding to parent scope?

#set(BETR_LIBRARIES bgcfarm_util;cent_ecacnp;cent_ecadcnp;${BETR_LIBRARIES} PARENT_SCOPE)
set(BETR_LIBRARIES bgcfarm_util;cent_ecacnp;UGM-M;${BETR_LIBRARIES} PARENT_SCOPE)
