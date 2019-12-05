APPEND_SOURCE_FILES(SOURCES)

add_library(DPHMPhaseField ${SOURCES})
target_link_libraries(DPHMPhaseField PUBLIC ProcessLib)

include(Tests.cmake)
