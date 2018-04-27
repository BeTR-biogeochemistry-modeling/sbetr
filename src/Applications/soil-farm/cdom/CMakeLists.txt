add_subdirectory(cdomPara)
add_subdirectory(cdom1layer)
add_subdirectory(cdomNlayer)

set(BETR_LIBRARIES cdomPara cdom1layer cdomNlayer;${BETR_LIBRARIES} PARENT_SCOPE)
