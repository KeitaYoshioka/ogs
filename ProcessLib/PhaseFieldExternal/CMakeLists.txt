APPEND_SOURCE_FILES(SOURCES)

add_library(PhaseFieldExternal ${SOURCES})
if(BUILD_SHARED_LIBS)
    install(TARGETS PhaseFieldExternal LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif()
target_link_libraries(PhaseFieldExternal PUBLIC ProcessLib)

include(Tests.cmake)
