cmake_minimum_required(VERSION 3.5)

project(StateEstimationFilter VERSION 1.0 LANGUAGES CXX)

set(CMAKE_INCLUDE_CURRENT_DIR ON)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Eigen3 REQUIRED)
#find_package(OpenCV COMPONENTS core highgui imgproc REQUIRED)
#find_package(Qt5 COMPONENTS Widgets Core Gui OpenGL REQUIRED)
#find_package(OpenGL REQUIRED)

include_directories(
  include
  ${EIGEN3_INCLUDE_DIR}
  #${OpenCV_INCLUDE_DIRS}
  #${QT_QTOPENGL_INCLUDE_DIR}
  #${OPENGL_INCLUDE_DIR}
)

#message(STATUS "FOUND OpenCV with version: ${OpenCV_VERSION}")

link_directories(
  #${OPENCV_LIBRARY_DIRS}
  #C:\\OpenCV\\opencv-4.5.4\\build\\install\\x64\\vc16\\lib
)

set(PROJECT_SOURCES
  FilterParameters.hpp
  FilterUtils.hpp
  FilterUtils.cpp
  RingBuffer.hpp
  Measurement.hpp
  State.hpp
  State2d.hpp
  State3d.hpp
  Gps.hpp
  GpsConversion.hpp
  GpsConversion.cpp
  Ekf.hpp
  Ekf.cpp
  Ekf2d.hpp
  Ekf2d.cpp
  Ekf3d.hpp
  Ekf3d.cpp
)

add_executable(unit_test unit_test.cpp State.hpp State3d.hpp State2d.hpp FilterUtils.hpp Measurement.hpp)
target_link_libraries(unit_test PRIVATE)

add_executable(test_ekf test_ekf.cpp ${PROJECT_SOURCES})
target_link_libraries(test_ekf PRIVATE)

#set_target_properties(state_est_ekf PROPERTIES
#    MACOSX_BUNDLE_GUI_IDENTIFIER my.example.com
#    MACOSX_BUNDLE_BUNDLE_VERSION ${PROJECT_VERSION}
#    MACOSX_BUNDLE_SHORT_VERSION_STRING ${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}
#)

