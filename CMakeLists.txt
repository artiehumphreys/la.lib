set(CMAKE_CXX_COMPILER "/opt/homebrew/opt/llvm/bin/clang++")
set(CMAKE_C_COMPILER   "/opt/homebrew/opt/llvm/bin/clang")

cmake_minimum_required(VERSION 3.30)
project(lalib LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include(CTest)

add_subdirectory(tests)
