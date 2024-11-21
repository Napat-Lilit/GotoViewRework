
# Set the source files, including both C++ and C files
set(SOURCES
    main.cpp
    
    error.cpp
    vec3.cpp
    transform.cpp
)

# Specify where to place the executable & Add the executable with the listed source files
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR})
add_executable(${PROJECT_NAME} ${SOURCES})

# Add include directories
target_include_directories(${PROJECT_NAME} PRIVATE
    ${PROJECT_SOURCE_DIR}/src/include
)

# Add specific compile options/flags
target_compile_options(${PROJECT_NAME} PRIVATE 
    -fopenmp
)
