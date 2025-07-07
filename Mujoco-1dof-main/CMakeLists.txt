cmake_minimum_required(VERSION 3.10)
project(mujoco_1dof)

add_compile_options(-std=c++17)
set(CMAKE_CXX_FLAGS_RELEASE_INIT "-Wall")
enable_language(C)
enable_language(CXX)
set(CMAKE_BUILD_TYPE Release)

set(MUJOCO_PATH "${PROJECT_SOURCE_DIR}/include/mujoco-2.3.1")
set(SIMULATE_PATH "${MUJOCO_PATH}/simulate")

#for visualstudio cmake configuration
set(pinocchio_DIR "/opt/openrobots/lib/cmake/pinocchio/")
set(hpp-fcl_DIR "/opt/openrobots/lib/cmake/hpp-fcl/")
set(example-robot-data_DIR "/opt/openrobots/lib/cmake/example-robot-data/")

find_package(pinocchio)
find_package(Curses REQUIRED)
find_package(glfw3 REQUIRED)
find_package(Threads REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(OpenGL REQUIRED COMPONENTS OpenGL EGL)
find_package(Boost REQUIRED COMPONENTS system)
find_package(pinocchio REQUIRED)


########## for ros ################
find_package(catkin REQUIRED COMPONENTS
  geometry_msgs
  roscpp
  rospy
  sensor_msgs
  std_msgs
  trajectory_msgs
  actionlib_msgs 
  actionlib
  # message_generation
)

catkin_package(
INCLUDE_DIRS include
# #  LIBRARIES state_listener
CATKIN_DEPENDS geometry_msgs roscpp rospy std_msgs trajectory_msgs actionlib_msgs actionlib
 DEPENDS system_lib
)

include_directories(
    # include
    ${catkin_INCLUDE_DIRS}
)
########## for ros END ##############

find_library(MUJOCO
  NAMES libmujoco.so.2.3.1
  PATHS ${MUJOCO_PATH}/lib
  REQUIRED
)

add_library(lodepng STATIC ${SIMULATE_PATH}/lodepng.h ${SIMULATE_PATH}/lodepng.cpp)
add_library(libsimulate STATIC)
target_sources(
  libsimulate
  PUBLIC  ${SIMULATE_PATH}/simulate.h
  PRIVATE ${SIMULATE_PATH}/simulate.cc
          ${SIMULATE_PATH}/array_safety.h
          ${SIMULATE_PATH}/uitools.h
          ${SIMULATE_PATH}/uitools.cc
          ${SIMULATE_PATH}/glfw_dispatch.h
          ${SIMULATE_PATH}/glfw_dispatch.cc
)
target_link_libraries(libsimulate PUBLIC lodepng ${MUJOCO} ${catkin_LIBRARIES} )
target_link_options(libsimulate PRIVATE ${MUJOCO_SIMULATE_LINK_OPTIONS})
target_compile_options(libsimulate PUBLIC ${MUJOCO_SIMULATE_COMPILE_OPTIONS})
target_include_directories(libsimulate 
  PUBLIC 
    ${SIMULATE_PATH}
    ${MUJOCO_PATH}/include
    ${MUJOCO_PATH}/simulate
    ${MUJOCO_PATH}/sample
    ${catkin_INCLUDE_DIRS}
)
target_compile_features(libsimulate PRIVATE cxx_std_17)

###########
## Build ##
###########

## Specify additional locations of header files
## Your package locations should be listed before other locations
include_directories(
  include
  ${EIGEN3_INCLUDE_DIRS}
  ${catkin_INCLUDE_DIRS}
  ${MUJOCO_PATH}/include/mujoco
  ${MUJOCO_PATH}/include
  ${MUJOCO_PATH}/simulate
  ${MUJOCO_PATH}/sample
  ${GLFW_INCLUDE_DIR}
  ${OPENGL_INCLUDE_DIR}
  ${CURSES_INCLUDE_DIR}
)


## Declare a C++ library
# add_library(${PROJECT_NAME}
#   src/simulate_flexible_panda.cc
#   ${UITOOLS_PATH}/uitools.c
#   ${LUCA_DYN}
#   ${TINYXML}
#   ${LUCA_PATH}/include/tinyxml/tinyxml.h
# )

## Add cmake target dependencies of the library
## as an example, code may need to be generated before libraries
## either from message generation or dynamic reconfigure
# add_dependencies(${PROJECT_NAME} ${${PROJECT_NAME}_EXPORTED_TARGETS} ${catkin_EXPORTED_TARGETS})

## Declare a C++ executable
## With catkin_make all packages are built within a single CMake context
## The recommended prefix ensures that target names across packages don't collide


link_libraries(
  ${catkin_LIBRARIES}
  lodepng
  Eigen3::Eigen
  glfw 
  ${GLFW3_LIBRARY} 
  libsimulate
  ${CURSES_LIBRARY}
  -lpthread  
  ${MUJOCO}
  pinocchio::pinocchio
  
)

add_library(MujocoController STATIC src/Controller/ControllerNode.cpp  
                                    src/Controller/Controller.cpp
                                    src/Dynamics/Dynamics.cpp)

add_executable(simulate_1dof
  src/main.cpp
  ${SIMULATE_PATH}/array_safety.h
  ${SIMULATE_PATH}/glfw_dispatch.cc
  ${SIMULATE_PATH}/glfw_dispatch.h
)

target_link_libraries(simulate_1dof
MujocoController)

