APPEND_SOURCE_FILES(SOURCES)

add_library(PhaseFieldInSitu ${SOURCES})
target_link_libraries(PhaseFieldInSitu PUBLIC ProcessLib)

include(Tests.cmake)
