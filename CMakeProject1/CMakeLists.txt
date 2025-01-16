# CMake minimum version
cmake_minimum_required(VERSION 3.8)

# Project name
project(CMakeProject1)

# Set the C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Include RDKit directories
set(RDKIT_INCLUDE_DIR "C:/Users/Abdullah/anaconda3/envs/rdkit_env/include")
set(RDKIT_LIBRARY_DIR "C:/Users/Abdullah/anaconda3/envs/rdkit_env/lib")

include_directories(${RDKIT_INCLUDE_DIR})
link_directories(${RDKIT_LIBRARY_DIR})

# Add source files
add_executable(CMakeProject1 "CMakeProject1.cpp" "CMakeProject1.h")

# Link RDKit libraries
target_link_libraries(CMakeProject1 
    RDKitSmilesParse 
    RDKitFileParsers 
    RDKitGraphMol 
    RDKitRDGeometryLib
)

# Optional: Add debug information or other customizations
if (CMAKE_VERSION VERSION_GREATER 3.12)
  set_property(TARGET CMakeProject1 PROPERTY CXX_STANDARD 20)
endif()
