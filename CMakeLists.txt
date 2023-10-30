cmake_minimum_required(VERSION 3.12)
project(scope)

set(CMAKE_CXX_STANDARD 14)
IF(NOT APPLE)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
ENDIF()
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address -O1 -fno-omit-frame-pointer")

include_directories(
        ${PROJECT_SOURCE_DIR}/config
        ${PROJECT_SOURCE_DIR}/utility
        ${PROJECT_SOURCE_DIR}/utility/automorphism
        ${PROJECT_SOURCE_DIR}/graph
        ${PROJECT_SOURCE_DIR}/plan
        ${PROJECT_SOURCE_DIR}/counting
)

link_directories(
        ${PROJECT_SOURCE_DIR}/config
        ${PROJECT_SOURCE_DIR}/utility
        ${PROJECT_SOURCE_DIR}/utility/automorphism
        ${PROJECT_SOURCE_DIR}/graph
        ${PROJECT_SOURCE_DIR}/plan
        ${PROJECT_SOURCE_DIR}/counting
)

# link to GLPK
IF(NOT APPLE)
    include_directories($ENV{HOME}/.local/include)
    link_directories($ENV{HOME}/.local/lib)
ELSE()
    include_directories(/opt/homebrew/Cellar/glpk/5.0/include)
    link_directories(/opt/homebrew/Cellar/glpk/5.0/lib)
ENDIF()

add_subdirectory(config)
add_subdirectory(graph)
add_subdirectory(utility)
add_subdirectory(executable)
add_subdirectory(plan)
add_subdirectory(counting)

add_executable(test.out test.cpp)
target_link_libraries(test.out graph utility plan counting)
