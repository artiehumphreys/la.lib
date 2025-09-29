add_executable(test_scalar_math test_scalar_math.cpp)
target_include_directories(test_scalar_math PRIVATE ${PROJECT_SOURCE_DIR}/include)

add_test(NAME scalar_math COMMAND $<TARGET_FILE:test_scalar_math>)
