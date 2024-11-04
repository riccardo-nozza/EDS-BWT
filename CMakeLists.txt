cmake_minimum_required(VERSION 3.13.4)
project(EDSBWT)

set(CMAKE_CXX_STANDARD 20)
set(CC g++)
set(OMP 0)
set(FASTQ 0)
set(RECOVERBW 1)

if(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  add_compile_options(-Wall -Wextra -Wpedantic -O3 -Ofast)
endif()

add_definitions(-DRECOVERBW=${RECOVERBW} -DOMP=${OMP} -DFASTQ=${FASTQ})

include_directories(/home/include 
Move-r/include/move_r
#Move-r/external/malloc_count
#${CMAKE_SOURCE_DIR}
)
link_directories(/home/lib
#submodules/Move-r
)

add_executable(MOVE_EDSBWTSearch 
mainMove_EDSBWT.cpp
MOVE_EDSBWTSearch.cpp
Sorting.cpp
#malloc_count/malloc_count.c
)
add_executable(build_MLF
build_MLF.cpp
)

add_subdirectory(Move-r)

target_link_libraries(MOVE_EDSBWTSearch PRIVATE sdsl divsufsort divsufsort64 dl move_r)  # Link required libraries
target_link_libraries(build_MLF PRIVATE sdsl divsufsort divsufsort64 dl move_r)  # Link required libraries