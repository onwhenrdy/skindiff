# cmake -DCMAKE_BUILD_TYPE=Release
# make

cmake_minimum_required(VERSION 2.8)

project(dskincmd)

# version settings
set(MAJOR_V "1")
set(MINOR_V "5")
set(PATCH_V "0")

# hg version
find_program(HG hg DOC "Mercurial executable file")
if (HG)
    execute_process(COMMAND ${HG} id -i WORKING_DIRECTORY ${dskincmd_SOURCE_DIR} OUTPUT_VARIABLE REVISION OUTPUT_STRIP_TRAILING_WHITESPACE)
    set(BUILD_NR "${REVISION}")
else ()
    set(BUILD_NR "unknown")
endif(HG)

if (UNIX)
    # warnings
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Werror -Wcast-qual -Wconversion-null\
        -Wmissing-declarations -Woverlength-strings\
        -Wpointer-arith -Wunused-local-typedefs -Wunused-result\
        -Wvarargs -Wvla -Wwrite-strings")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-sign-compare ")
    # std and performance
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -fno-rtti -march=native")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -flto")
    # linker
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lz")
    set(CMAKE_EXE_LINKER_FLAGS_RELEASE "${CMAKE_EXE_LINKER_FLAGS_RELEASE} -flto")
endif(UNIX)

if (WINDOWS)
    # std etc
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /std:c++latest /nologo")
    # warnings
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4 /WX /permissive-")
    # performance
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /O2 /GL /GR- /GS-")
    # linker
    set(CMAKE_EXE_LINKER_FLAGS_RELEASE "/LTCG")
endif (WINDOWS)


# defines
add_definitions(-DDSCMDAPPNAME=${PROJECT_NAME})
add_definitions(-DDSCMDMAJVERSION=${MAJOR_V})
add_definitions(-DDSCMDMINVERSION=${MINOR_V})
add_definitions(-DDSCMDPATCHLEVEL=${PATCH_V})
add_definitions(-DDSCMDBUILDNUMBER=${BUILD_NR})

include_directories(SYSTEM includes)
set(project_sources
   cmdlineparser.cpp
   geometry.cpp
   compartment.cpp
   sink.cpp
   system.cpp
   helper.cpp
   tdmatrix.cpp
   algorithms.cpp
   matrixbuilder.cpp
   compartmentlog2d.cpp
   session.cpp
   parameter.cpp
   jsonparser.cpp
   consoleprogressbar.cpp
   compartmentlog3d.cpp
   versioninfo.cpp
   systemcmd.cpp
)

FILE(GLOB_RECURSE libFiles "*.h")
FILE(GLOB_RECURSE incFiles "*.hpp")
add_custom_target(headers SOURCES ${libFiles})
add_custom_target(includes SOURCES ${incFiles})
add_custom_target(sources SOURCES ${project_sources})
add_custom_target(test_files SOURCES ${testFiles})

add_executable(${PROJECT_NAME} main.cpp ${project_sources})
