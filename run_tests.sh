rm -rf build
cmake -B build -S . -DCMAKE_BUILD_TYPE=Debug
cmake --build build --target test_scalar_math
ctest --test-dir build --output-on-failure
