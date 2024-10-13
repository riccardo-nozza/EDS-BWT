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

include_directories(/home/include 
submodules/Move-r/include/move_r
)
link_directories(/home/lib
#submodules/Move-r
)


add_executable(M_LF_build M_LF_build.cpp)


add_subdirectory(submodules/Move-r)

target_link_libraries(M_LF_build PRIVATE sdsl divsufsort divsufsort64 dl move_r)  # Link required libraries
