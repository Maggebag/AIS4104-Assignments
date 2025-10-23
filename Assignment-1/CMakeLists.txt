find_package(Eigen3 CONFIG REQUIRED)

add_library(assignment1_lib
                Assignment-1.cpp
)

target_include_directories(assignment1_lib PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(assignment1_lib PUBLIC Eigen3::Eigen)

add_executable(test-assignment-1
        main-1.cpp
)

target_link_libraries(test-assignment-1 PRIVATE assignment1_lib)