add_subdirectory(ch4lakePara)
add_subdirectory(ch4lake1layer)
add_subdirectory(ch4lakeNlayer)

set(BETR_LIBRARIES ch4lakePara ch4lake1layer ch4lakeNlayer;${BETR_LIBRARIES} PARENT_SCOPE)
