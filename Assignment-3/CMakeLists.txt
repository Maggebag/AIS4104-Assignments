find_package(Eigen3 CONFIG REQUIRED)

add_library(assignment3_lib
        Assignment-3.cpp
)

target_include_directories(assignment3_lib PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(assignment3_lib PUBLIC assignment1_lib assignment2_lib Eigen3::Eigen)

add_executable(test-assignment-3
        main-3.cpp
)

target_link_libraries(test-assignment-3 PRIVATE assignment3_lib)