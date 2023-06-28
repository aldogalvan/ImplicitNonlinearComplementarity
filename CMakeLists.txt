cmake_minimum_required(VERSION 3.25)
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
set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -D__CUDACC__")

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
include_directories(external/collisiondetection/include)

#add_executable(ImplicitNonlinearComplementarity main1.cpp implicit_lcp.cpp implicit_lcp.hpp lcp.cpp lcp.hpp collision.cpp collision.hpp aabb.cpp aabb.hpp objects.hpp helper.hpp)
#add_executable(ImplicitNonlinearComplementarity main2.cpp implicit_lcp.cpp implicit_lcp.hpp lcp.cpp lcp.hpp collision.cpp collision.hpp aabb.cpp aabb.hpp objects.hpp helper.hpp)
#add_executable(ImplicitNonlinearComplementarity main3.cpp tetgen/*.cpp tetgen/*.h
#        implicit_lcp.cpp implicit_lcp.hpp lcp.cpp lcp.hpp collision.cpp collision.hpp
#        aabb.cpp aabb.hpp objects.hpp helper.hpp)

add_library(ImplicitLCP src/lcp.cpp src/lcp.hpp src/implicit_lcp.cpp src/implicit_lcp.hpp
        src/collision.cpp src/collision.hpp src/aabb.cpp src/aabb.hpp
        src/objects.cpp src/objects.hpp src/helper.hpp src/solvers.cpp
        src/pgs.hpp src/pgs.cpp)

target_link_libraries(ImplicitLCP PUBLIC
        collisiondetection
        ${CHAI3D_LIBRARIES}
        ${GLFW_LIBRARIES}
        Eigen3::Eigen
        drake::drake)

add_executable(RigidBodiesDemo rigid_bodies_demo.cpp)
add_executable(HapticsDemoQuasistatic haptics_demo_quasistatic.cpp)
add_executable(HapticsDemo haptics_demo.cpp)

add_definitions(-DTETLIBRARY)
target_link_libraries(RigidBodiesDemo PUBLIC
        ImplicitLCP)
target_link_libraries(HapticsDemoQuasistatic PUBLIC
        ImplicitLCP)
target_link_libraries(HapticsDemo PUBLIC
        ImplicitLCP)


