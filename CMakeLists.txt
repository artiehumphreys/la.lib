set(CMAKE_CXX_COMPILER "/opt/homebrew/opt/llvm/bin/clang++")
set(CMAKE_C_COMPILER "/opt/homebrew/opt/llvm/bin/clang")

cmake_minimum_required(VERSION 3.30)
project(lalib LANGUAGES CXX)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_library(lalib INTERFACE)

target_include_directories(lalib INTERFACE ${CMAKE_SOURCE_DIR}/include)

include(CTest)

add_subdirectory(tests)

if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    target_compile_options(lalib INTERFACE -ftime-trace)

    target_compile_options(lalib INTERFACE
        -Rpass=loop-vectorize
        -Rpass-missed=loop-vectorize
        -Rpass-analysis=loop-vectorize
    )

    target_compile_options(lalib INTERFACE
        -Rpass=inline
        -Rpass-missed=inline
        -Rpass-analysis=inline
    )
endif()
