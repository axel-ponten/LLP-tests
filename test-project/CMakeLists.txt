cmake_minimum_required(VERSION 2.6)
set (CMAKE_CXX_STANDARD 11)

add_executable(example src/example.cxx)

find_library(PROPOSAL_LIBRARIES REQUIRED NAMES PROPOSAL)

if(PROPOSAL_LIBRARIES)
    target_link_libraries(example ${PROPOSAL_LIBRARIES})
endif()

