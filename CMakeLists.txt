# Minimum CMake version required
cmake_minimum_required(VERSION 3.10)

# Project name and version 
project(GotoSimplePathTracer VERSION 1.0 LANGUAGES CXX C)

# Specify C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Add the 'src' directory
add_subdirectory(src)
