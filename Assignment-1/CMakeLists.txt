find_package(Eigen3 CONFIG REQUIRED)

add_executable(test-assignment-1 Assignment-1-main.cpp)

target_link_libraries(test-assignment-1 PRIVATE Eigen3::Eigen)