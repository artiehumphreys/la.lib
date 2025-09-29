rm -rf build
cmake -S . -B build
cmake --build build --target test_scalar_math
ctest --test-dir build --output-on-failure
