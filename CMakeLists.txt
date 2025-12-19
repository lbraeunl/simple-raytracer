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

add_subdirectory(external/tinyobjloader)
add_subdirectory(external/yaml-cpp)

# Executable
add_executable(main src/main.cpp src/Geometry.cpp src/Model.cpp src/Model.hpp src/Scene.cpp src/Scene.hpp src/BVH.hpp src/BVH.cpp src/Renderer.hpp src/Renderer.cpp src/Postprocessing.hpp src/Postprocessing.cpp)
target_compile_features(main PRIVATE cxx_std_17)

target_include_directories(main PRIVATE external/glm external/tinycolormap/include external/stb)

# Link libraries
target_link_libraries(main PRIVATE SFML::Graphics tinyobjloader yaml-cpp)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

