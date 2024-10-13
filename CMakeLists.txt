cmake_minimum_required(VERSION 3.13.4)
project(EDSBWT)

set(CMAKE_CXX_STANDARD 20)
set(CC g++)
set(OMP 0)
set(FASTQ 0)
set(RECOVERBW 1)

if(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  add_compile_options(-Wall -Wextra -Wpedantic)
endif()

add_definitions(-DRECOVERBW=${RECOVERBW} -DOMP=${OMP} -DFASTQ=${FASTQ})

include_directories(/home/riccardo/include 
#submodules/Move-r/include
)
link_directories(/home/riccardo/lib
#submodules/Move-r
)


add_executable(da_to_everything da_to_everything.cpp)
add_executable(M_LF_build M_LF_build.cpp malloc_count/malloc_count.c)
#add_executable(eds_to_fasta eds_to_fasta.cpp)
add_executable(stringCheck stringCheck.cpp)
add_executable(EDSBWTSearch
  mainEDS-BWT.cpp
  EDSBWTsearch.cpp
  Sorting.cpp
  malloc_count/malloc_count.c
)


add_subdirectory(submodules/Move-r)

target_link_libraries(da_to_everything PRIVATE sdsl divsufsort divsufsort64 dl)  # Link required libraries
target_link_libraries(M_LF_build PRIVATE sdsl divsufsort divsufsort64 dl move_r)  # Link required libraries
#target_link_libraries(eds_to_fasta PRIVATE sdsl divsufsort divsufsort64 dl)  # Link required libraries
target_link_libraries(stringCheck PRIVATE sdsl divsufsort divsufsort64 dl)  # Link required libraries
target_link_libraries(EDSBWTSearch PRIVATE sdsl divsufsort divsufsort64 dl)  # Link required libraries
