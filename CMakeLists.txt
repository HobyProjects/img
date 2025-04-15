cmake_minimum_required(VERSION 3.15)
project(image LANGUAGES C CXX)

##############################################################################
# Compiler Configuration
##############################################################################
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_C_STANDARD 17)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS OFF)

##############################################################################
# Version Configuration
##############################################################################
set(IMG_VERSIONS_CONTROL_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/config)
set(IMG_VERSION_CONTROL_FILE ${IMG_VERSIONS_CONTROL_INCLUDE_DIR}/version.txt)

file(READ ${IMG_VERSION_CONTROL_FILE} VERSION_STR)
string(STRIP "${VERSION_STR}" VERSION_STR)  # remove newline

# Split into major.minor.patch
string(REPLACE "." ";" VERSION_LIST ${VERSION_STR})
list(GET VERSION_LIST 0 IMG_VERSION_MAJOR)
list(GET VERSION_LIST 1 IMG_VERSION_MINOR)
list(GET VERSION_LIST 2 IMG_VERSION_PATCH)

# Convert to number and increment patch
math(EXPR IMG_VERSION_PATCH "${IMG_VERSION_PATCH} + 1")

# Handle overflow
if(IMG_VERSION_PATCH GREATER 999)
    set(IMG_VERSION_PATCH 0)
    math(EXPR IMG_VERSION_MINOR "${IMG_VERSION_MINOR} + 1")
endif()

if(IMG_VERSION_MINOR GREATER 999)
    set(IMG_VERSION_MINOR 0)
    math(EXPR IMG_VERSION_MAJOR "${IMG_VERSION_MAJOR} + 1")
endif()

if(IMG_VERSION_MAJOR GREATER 999)
    set(IMG_VERSION_MAJOR 0)
    message(WARNING "Version major number overflowed, That's mean it's time to say goodbye to this project.")
endif()

# Update version.txt with the new version
file(WRITE ${IMG_VERSION_CONTROL_FILE}
    "${IMG_VERSION_MAJOR}.${IMG_VERSION_MINOR}.${IMG_VERSION_PATCH}\n")

set(IMG_VERION_STRING "${IMG_VERSION_MAJOR}.${IMG_VERSION_MINOR}.${IMG_VERSION_PATCH}")
message(STATUS "${PROJECT_NAME} Version: ${IMG_VERION_STRING}")

set(IMG_VERSION_IN ${CMAKE_CURRENT_SOURCE_DIR}/config/img_version.hpp.in)
set(IMG_VERSION_OUT ${CMAKE_BINARY_DIR}/config/img_version.hpp)
set(IMG_VERSION_OUT_DIR ${CMAKE_BINARY_DIR}/config)

# Generate the header file
configure_file(
    ${IMG_VERSION_IN}
    ${IMG_VERSION_OUT}
    @ONLY
)

###############################################################################
# Platform Specific Configuration
###############################################################################
if(CMAKE_SYSTEM_NAME STREQUAL "Windows")
    add_compile_definitions(IMG_PLATFORM_WINDOWS)
elseif(CMAKE_SYSTEM_NAME STREQUAL "Linux")
    add_compile_definitions(IMG_PLATFORM_LINUX)
endif()

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    add_compile_definitions(IMG_DEBUG)
else()
    add_compile_definitions(IMG_RELEASE)
endif()


###############################################################################
# Generate the library
###############################################################################
add_library(${PROJECT_NAME} STATIC img.hpp img.cpp)
target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${IMG_VERSION_OUT_DIR})
target_precompile_headers(${PROJECT_NAME} PUBLIC img_pch.hpp)