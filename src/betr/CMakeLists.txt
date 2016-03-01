add_subdirectory(betr_math)
add_subdirectory(betr_core)
add_subdirectory(betr_main)
set(BETR_LIBRARIES betr_math;betr_core;betr_main;${BETR_LIBRARIES} PARENT_SCOPE)

