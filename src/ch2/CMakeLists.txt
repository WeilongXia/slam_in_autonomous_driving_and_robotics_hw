# find_package(Eigen3 REQUIRED)
# find_package(OpenCV REQUIRED)
# include_directories(${EIGEN3_INCLUDE_DIRS})
# include_directories(${OpenCV_INCLUDE_DIRS})
add_executable(motion motion.cc)
target_link_libraries(motion
        ${PROJECT_NAME}.common
        ${PROJECT_NAME}.tools
        )
