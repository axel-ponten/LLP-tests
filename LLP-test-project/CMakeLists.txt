cmake_minimum_required(VERSION 2.6)
set (CMAKE_CXX_STANDARD 11)

add_executable(main src/main.cxx)

find_library(PROPOSAL_LIBRARIES REQUIRED NAMES PROPOSAL)

if(PROPOSAL_LIBRARIES)
    target_link_libraries(main ${PROPOSAL_LIBRARIES})
endif()

