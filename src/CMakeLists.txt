if (GSL_FOUND)
set(SOURCES
INIReader.cpp
eos.cpp
eos_base.cpp
eos_neos.cpp
grid.cpp
main.cpp
parameter.cpp
reader.cpp
ini.cpp  
)
endif()

add_library(${libname} SHARED ${SOURCES})
set(CompileFlags "")
if (APPLE)
    set(CompileFlags "${CompileFlags} -DAPPLE")
endif (APPLE)
if (GSL_FOUND)
    set(CompileFlags "${CompileFlags} -DGSL")
endif()


add_executable (${exename} main.cpp ${SOURCES})

if (GSL_FOUND)
    target_link_libraries (${exename} ${GSL_LIBRARIES})
endif (GSL_FOUND)

set_target_properties (${exename} PROPERTIES COMPILE_FLAGS "${CompileFlags}")
install(TARGETS ${exename} DESTINATION ${CMAKE_HOME_DIRECTORY})



