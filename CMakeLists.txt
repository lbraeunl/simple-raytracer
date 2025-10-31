cmake_minimum_required(VERSION 3.28)
project(CMakeSFMLProject LANGUAGES CXX)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

include(FetchContent)

# SFML
FetchContent_Declare(SFML
    GIT_REPOSITORY https://github.com/SFML/SFML.git
    GIT_TAG 3.0.2
    GIT_SHALLOW ON
    EXCLUDE_FROM_ALL
    SYSTEM)
FetchContent_MakeAvailable(SFML)

# tinyobjloader
add_subdirectory(external/tinyobjloader)

# Executable
add_executable(main src/main.cpp)
target_compile_features(main PRIVATE cxx_std_17)

# GLM
target_include_directories(main PRIVATE external/glm)

# Link libraries
target_link_libraries(main PRIVATE SFML::Graphics tinyobjloader)
