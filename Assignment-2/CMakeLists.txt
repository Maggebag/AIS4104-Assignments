find_package(Eigen3 CONFIG REQUIRED)

add_executable(test-assignment-2 Assignment-2-main.cpp)

target_link_libraries(test-assignment-2 PRIVATE Eigen3::Eigen)