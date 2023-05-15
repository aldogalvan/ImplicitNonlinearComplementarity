cmake_minimum_required(VERSION 3.25)
project(ImplicitNonlinearComplementarity)

set(CMAKE_CXX_STANDARD 17)


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
set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -D__CUDACC__")

# GLFW
if (EXISTS ${PROJECT_SOURCE_DIR}/extras/GLFW)
    add_subdirectory (${PROJECT_SOURCE_DIR}/extras/GLFW)
endif ()

add_executable(ImplicitNonlinearComplementarity main.cpp)

target_link_libraries(ImplicitNonlinearComplementarity PUBLIC
        ${CHAI3D_LIBRARIES}
        ${GLFW_LIBRARIES}
        Eigen3::Eigen)

