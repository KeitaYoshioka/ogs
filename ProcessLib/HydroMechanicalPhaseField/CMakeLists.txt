APPEND_SOURCE_FILES(SOURCES)

add_library(HydroMechanicalPhaseField ${SOURCES})
target_link_libraries(HydroMechanicalPhaseField PUBLIC ProcessLib)

include(Tests.cmake)
