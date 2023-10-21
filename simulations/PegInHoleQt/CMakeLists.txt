cmake_minimum_required(VERSION 3.1)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)
set(CMAKE_AUTOUIC ON)


find_package(Qt5 COMPONENTS
        Core
        Gui
        Widgets
        OpenGL
        REQUIRED)


# Add a compiler flag
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")

# Add the QT executables
add_executable(PegInHoleQt main.cpp PegInHole.cpp Application.cpp Application.h Application.qrc Interface.cpp Interface.h Interface.ui)

# Add the Qt5 Widgets for linking
target_link_libraries(PegInHoleQt
        Qt5::OpenGL
        Qt5::Core
        Qt5::Gui
        Qt5::Widgets
        ImplicitLCP)