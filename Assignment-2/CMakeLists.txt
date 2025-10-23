find_package(Eigen3 CONFIG REQUIRED)

add_library(assignment2_lib
        Assignment-2.cpp
)

target_include_directories(assignment2_lib PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(assignment2_lib PUBLIC assignment1_lib Eigen3::Eigen)

add_executable(test-assignment-2
        main-2.cpp
)

target_link_libraries(test-assignment-2 PRIVATE assignment2_lib)