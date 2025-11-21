set(CMAKE_CXX_COMPILER "/opt/homebrew/opt/llvm/bin/clang++")
set(CMAKE_C_COMPILER   "/opt/homebrew/opt/llvm/bin/clang")

cmake_minimum_required(VERSION 3.30)
project(lalib LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_library(lalib INTERFACE)

target_include_directories(lalib INTERFACE ${CMAKE_SOURCE_DIR}/include)

include(CTest)

add_subdirectory(tests)
