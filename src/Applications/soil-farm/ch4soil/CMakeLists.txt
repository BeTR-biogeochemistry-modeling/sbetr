add_subdirectory(ch4soilPara)
add_subdirectory(ch4soil1layer)
add_subdirectory(ch4soilNlayer)

set(BETR_LIBRARIES ch4soilPara ch4soil1layer ch4soilNlayer;${BETR_LIBRARIES} PARENT_SCOPE)
