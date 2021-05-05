CMAKE_MINIMUM_REQUIRED(VERSION 3.9)
#CMAKE_POLICY(VERSION 3.12)

PROJECT(shapelydrugs LANGUAGES CXX VERSION 0.0.1.0)

# Stop GCC builds failing on MacOS
IF(APPLE AND (CMAKE_CXX_COMPILER_ID STREQUAL "GNU"))
SET(CMAKE_OSX_SYSROOT "/")
SET(CMAKE_OSX_DEPLOYMENT_TARGET "")
ENDIF()

# We require c++-17 compiler
SET(CMAKE_CXX_STANDARD_REQUIRED ON)
SET(CMAKE_CXX_STANDARD 17)
SET(CMAKE_CXX_EXTENSIONS OFF)
SET(CMAKE_CXX_VISIBILITY_PRESET hidden)

# Export the compiler commands
SET(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Default to release mode build. Does Debug actually do anything?
IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE Release)
  MESSAGE(STATUS "Defaulting to Release mode build")
ENDIF()

# Setup some compile flags, for all versions and specific versions
IF("${CMAKE_SYSTEM_NAME}" STREQUAL "Linux")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
    #SET(CMAKE_AR "gcc-ar")
    #SET(CMAKE_CXX_ARCHIVE_CREATE "<CMAKE_AR> qcs <TARGET> <LINK_FLAGS> <OBJECTS>")
    #SET(CMAKE_CXX_ARCHIVE_FINISH   true)
ENDIF()
SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O3 -DDEBUG")
SET(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -fno-fast-math -DNDEBUG -flto ")
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/include/ ${CMAKE_CURRENT_SOURCE_DIR}/external/stdx/include ${CMAKE_CURRENT_SOURCE_DIR}/external/eigen ${CMAKE_CURRENT_SOURCE_DIR}/external/dynamic_bitset/include ${CMAKE_CURRENT_SOURCE_DIR}/external/dynamic_bitset/extlibs/libpopcnt ${CMAKE_CURRENT_SOURCE_DIR}/external/cereal/include)

# Setup install directories
INCLUDE(GNUInstallDirs)

# List all the files to compile into the physdock library
SET(SHAPELY_SRCS
  src/coordinates.cpp
  test/main_dev.cpp
)

# Build the executable
FIND_PACKAGE(OpenMP)
ADD_EXECUTABLE(shapely_tester ${SHAPELY_SRCS})
TARGET_LINK_LIBRARIES(shapely_tester PUBLIC stdc++fs OpenMP::OpenMP_CXX)
SET_TARGET_PROPERTIES(shapely_tester PROPERTIES COMPILE_FLAGS "-fopenmp")

