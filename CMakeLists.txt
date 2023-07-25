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

file(GLOB_RECURSE IMPLICIT_LCP_SOURCES  src/*.cpp src/*.hpp external/tetgen/*.cpp external/tetgen/*.h)
add_library(ImplicitLCP ${IMPLICIT_LCP_SOURCES})
add_definitions(-DTETLIBRARY)

target_link_libraries(ImplicitLCP PUBLIC
        collisiondetection
        ${CHAI3D_LIBRARIES}
        ${GLFW_LIBRARIES}
        Eigen3::Eigen
        drake::drake)


add_executable(PegInHole simulations/PegInHole/main.cpp simulations/PegInHole/PegInHole.cpp)
add_executable(RigidBodyDemo simulations/RigidBodyDemo/main.cpp simulations/RigidBodyDemo/RigidBodyDemo.cpp)
add_executable(DeformableBodyDemo simulations/DeformableBodyDemo/main.cpp simulations/DeformableBodyDemo/DeformableBodyDemo.cpp)

#target_link_libraries(RigidBodiesDemo PUBLIC
#        ImplicitLCP)
#target_link_libraries(HapticsDemoQuasistatic PUBLIC ImplicitLCP)
#target_link_libraries(HapticsDemo PUBLIC ImplicitLCP)
#target_link_libraries(HapticsDemoRigid PUBLIC ImplicitLCP)
#target_link_libraries(HapticsDemoDeformable PUBLIC ImplicitLCP)
target_link_libraries(PegInHole PUBLIC ImplicitLCP)
target_link_libraries(RigidBodyDemo PUBLIC ImplicitLCP)
target_link_libraries(DeformableBodyDemo PUBLIC ImplicitLCP)


add_subdirectory(simulations/PegInHoleQt)



