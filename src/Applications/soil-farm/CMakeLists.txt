add_subdirectory(bgcfarm_util)
<<<<<<< HEAD
#begin_appadd
add_subdirectory(ecacnp)
add_subdirectory(cdom)
add_subdirectory(simic)
add_subdirectory(keca)
add_subdirectory(ch4soil)
<<<<<<< HEAD
add_subdirectory(SUMMS)
=======
add_subdirectory(v1eca)
>>>>>>> jinyun_rr
#end_appadd
||||||| merged common ancestors
add_subdirectory(CENT_ECACNP)
add_subdirectory(cdom_ecacnp)
if (UGM)
  add_subdirectory(UGM-M)
endif()
#add_subdirectory(CENT_ECADCNP)
# FIXME(bja, 201603) is this manual list necessary with the subdir
# adding to parent scope?

#set(BETR_LIBRARIES bgcfarm_util;cent_ecadcnp;${BETR_LIBRARIES} PARENT_SCOPE)
if (UGM)
set(BETR_LIBRARIES bgcfarm_util;UGM-M;${BETR_LIBRARIES} PARENT_SCOPE)
else()
=======
add_subdirectory(CENT_ECACNP)
add_subdirectory(cdom_ecacnp)
add_subdirectory(SUMMS)
if (UGM)
  add_subdirectory(UGM-M)
endif()
#add_subdirectory(CENT_ECADCNP)
# FIXME(bja, 201603) is this manual list necessary with the subdir
# adding to parent scope?

#set(BETR_LIBRARIES bgcfarm_util;cent_ecadcnp;${BETR_LIBRARIES} PARENT_SCOPE)
if (UGM)
set(BETR_LIBRARIES bgcfarm_util;UGM-M;${BETR_LIBRARIES} PARENT_SCOPE)
else()
>>>>>>> rzacplsbetr_cmupdated
set(BETR_LIBRARIES bgcfarm_util;${BETR_LIBRARIES} PARENT_SCOPE)
