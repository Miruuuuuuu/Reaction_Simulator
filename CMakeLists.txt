# CMakeList.txt : Top-level CMake project file, do global configuration
# and include sub-projects here.
#
cmake_minimum_required (VERSION 3.8)

# Enable Hot Reload for MSVC compilers if supported.
if (POLICY CMP0141)
  cmake_policy(SET CMP0141 NEW)
  set(CMAKE_MSVC_DEBUG_INFORMATION_FORMAT "$<IF:$<AND:$<C_COMPILER_ID:MSVC>,$<CXX_COMPILER_ID:MSVC>>,$<$<CONFIG:Debug,RelWithDebInfo>:EditAndContinue>,$<$<CONFIG:Debug,RelWithDebInfo>:ProgramDatabase>>")
endif()

project ("CMakeProject1")

# Include sub-projects.
add_subdirectory ("CMakeProject1")


# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Add the RDKit libraries
# Update these paths according to your Conda environment
set(RDKIT_INCLUDE_DIR "C:/Users/Abdullah/anaconda3/envs/rdkit_env/include")
set(RDKIT_LIBRARY_DIR "C:/Users/Abdullah/anaconda3/envs/rdkit_env/lib")

# Include RDKit headers
include_directories(${RDKIT_INCLUDE_DIR})

# Link RDKit libraries
link_directories(${RDKIT_LIBRARY_DIR})

# Add your executable
add_executable(CMakeProject1 main.cpp) # Replace main.cpp with your actual source file(s)

# Link RDKit libraries to your project
target_link_libraries(CMakeProject1 RDKitSmilesParse RDKitFileParsers RDKitGraphMol RDKitRDGeometryLib)