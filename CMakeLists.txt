cmake_minimum_required(VERSION 3.1)
project(ImplicitNonlinearComplementarity)
set(CMAKE_CXX_STANDARD 17)
set(OpenGL_GL_PREFERENCE "GLVND")

#CHAI-3D
FIND_PACKAGE(CHAI3D REQUIRED CONFIG)
ADD_DEFINITIONS(${CHAI3D_DEFINITIONS})
if (CMAKE_SYSTEM_NAME MATCHES Darwin)
    ADD_DEFINITIONS(-D_MACOSX -DMACOSX)
elseif (CMAKE_SYSTEM_NAME MATCHES Linux)
    ADD_DEFINITIONS(-D_LINUX -DLINUX)
endif () #NOTE: Windows is currently unsupported
LINK_DIRECTORIES(${CHAI3D_LIBRARY_DIRS})
INCLUDE_DIRECTORIES(${CHAI3D_INCLUDE_DIRS})

#EIGEN3
find_package(Eigen3 3.3 REQUIRED NO_MODULE)

# GLFW
if (EXISTS ${PROJECT_SOURCE_DIR}/extras/GLFW)
    add_subdirectory (${PROJECT_SOURCE_DIR}/extras/GLFW)
endif ()

find_package(drake CONFIG REQUIRED PATHS /opt/drake)

get_filename_component(DRAKE_PYTHONPATH "${drake_DIR}" DIRECTORY)
get_filename_component(DRAKE_PYTHONPATH "${DRAKE_PYTHONPATH}" DIRECTORY)
set(DRAKE_PYTHONPATH
        "${DRAKE_PYTHONPATH}/python${PYTHON_VERSION_MAJOR_MINOR}/site-packages"
        )

add_subdirectory(external/collisiondetection)
include_directories(extras/GLFW/include)
include_directories(external/collisiondetection/include)
include_directories(external/tetgen)
include_directories(include/imgui)
include_directories(include/imgui/backends)
include_directories(src)
include_directories(external/libigl/include)

file(GLOB_RECURSE IMPLICIT_LCP_SOURCES  src/*.cpp src/*.hpp external/tetgen/*.cpp external/tetgen/*.h)
add_library(ImplicitLCP ${IMPLICIT_LCP_SOURCES})
add_definitions(-DTETLIBRARY)
find_package(OpenMP)

target_link_libraries(ImplicitLCP PUBLIC
        ${GNUPLOT_LIBRARIES}
        collisiondetection
        ${CHAI3D_LIBRARIES}
        ${GLFW_LIBRARIES}
        Eigen3::Eigen
        drake::drake
        OpenMP::OpenMP_CXX)

# Add the matplotlib-cpp library as an external project
add_subdirectory(external/matplotlib-cpp)
include_directories(external/matplotlib-cpp)
include_directories("/usr/include/python3.8")
#/home/aldo/anaconda3/lib/python3.9/site-packages/numpy/core/include


add_executable(PegInHole simulations/PegInHole/main.cpp simulations/PegInHole/PegInHole.cpp)
#add_executable(RigidBodyDemo simulations/RigidBodyDemo/main.cpp simulations/RigidBodyDemo/RigidBodyDemo.cpp)
#add_executable(DeformableBodyDemo simulations/DeformableBodyDemo/main.cpp simulations/DeformableBodyDemo/DeformableBodyDemo.cpp)
#add_executable(HapticDemo simulations/HapticDemo/main.cpp simulations/HapticDemo/HapticDemo.cpp)
add_executable(PenaltyDemo simulations/PenaltyDemo/main.cpp simulations/PenaltyDemo/PenaltyDemo.cpp)

target_link_libraries(PegInHole PUBLIC ImplicitLCP)
#target_link_libraries(RigidBodyDemo PUBLIC ImplicitLCP matplotlib_cpp)
#target_link_libraries(DeformableBodyDemo PUBLIC ImplicitLCP)
#target_link_libraries(HapticDemo PUBLIC ImplicitLCP)
target_link_libraries(PenaltyDemo PUBLIC ImplicitLCP)



