add_subdirectory(soil-farm)
add_subdirectory(app_util)
set(BETR_LIBRARIES app_util;${BETR_LIBRARIES} PARENT_SCOPE)
