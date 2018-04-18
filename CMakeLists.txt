cmake_minimum_required(VERSION 3.7.0)
project(OSVRClientTest)

set (CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR})

add_subdirectory(LibAvWrapper)

set(BOOST_ROOT "C:/Program Files (x86)/boost_1_57_0")
set(BOOST_INCLUDE_DIR "C:/Program Files (x86)/boost_1_57_0")
set(BOOST_LIBRARYDIR "C:/Program Files (x86)/boost_1_57_0/stage/win32/lib")
find_package( Boost 1.57 )

find_package(osvr REQUIRED)

set(SDL2_LIBRARIES "C:/Program Files (x86)/SDL2-2.0.7/lib/x86/SDL2.lib")
set(SDL2_INCLUDE_DIR "C:/Program Files (x86)/SDL2-2.0.7/include")
find_package(SDL2 REQUIRED)

find_package(OpenGL REQUIRED)

#set(OPENGL_gl_LIBRARY "C:/Program Files (x86)/Windows Kits/10/Lib/10.0.16299.0/um/x86/OpenGL32.Lib")
#set(OPENGL_glu_LIBRARY "C:/Program Files (x86)/Windows Kits/10/Lib/10.0.16299.0/um/x86/GlU32.Lib")
#set(OPENGL_INCLUDE_DIRS "C:/Program Files (x86)/Windows Kits/10/Include/10.0.16299.0/um")

set(GLEW_INCLUDE_DIR "C:/Program Files (x86)/glew-2.1.0/include")
set(GLEW_LIBRARY "C:/Program Files (x86)/glew-2.1.0/lib/Release/Win32/glew32.lib")
find_package(GLEW REQUIRED)

find_package(osvrRenderManager REQUIRED)

set(ZMQ_INCLUDE_DIRS "C:/Program Files (x86)/libzmq/include")
set(ZMQ_LIBRARIES "C:/Program Files (x86)/libzmq/bin/Win32/Release/v120/dynamic/libzmq.lib")
find_package(ZMQ REQUIRED)
#osvrrm_find_SDL2()
#osvrrm_find_GLEW()

find_package(Eigen3 REQUIRED)


set(Boost_LIBRARY_DIR "C:/Program Files (x86)/boost_1_57_0/stage/win32/lib")
set(OSVR_LIBRARY_DIRS "C:/Program Files (x86)/osvrcore/lib")
set(OSVR_INCLUDE_DIR "C:/Program Files (x86)/osvrcore/include")
#set(OPENGL_INCLUDE_DIRS "C:/Program Files (x86)/glew-2.1.0/include" )
set(SDL2_INCLUDE_DIRS "C:/Program Files (x86)/SDL2-2.0.7/include" )
set(SOIL_LIBRARIES "C:/Users/ShunHuaiYao/Documents/GitHub/Simple-OpenGL-Image-Library/VS2013_build/Release/SOIL.lib" )
set(SOIL_INCLUDE_DIRS "C:/Users/ShunHuaiYao/Documents/GitHub/Simple-OpenGL-Image-Library/src" )


LINK_DIRECTORIES(${OSVR_LIBRARY_DIRS} ${EIGEN3_LIBRARY_DIRS} ${Boost_LIBRARY_DIR})

FILE(GLOB OSVRClientTrackerSrc src/*.cpp)

add_executable(OSVRClientTest
  ${OSVRClientTrackerSrc}
)

MESSAGE( STATUS "BOOST_LIBRARIES:         " ${BOOST_LIBRARIES} )
MESSAGE( STATUS "Boost_LIBRARIES:         " ${Boost_LIBRARIES} )
MESSAGE( STATUS "Boost_LIBRARY_DIR:         " ${Boost_LIBRARY_DIR} )
MESSAGE( STATUS "Boost_INCLUDE_DIR:         " ${Boost_INCLUDE_DIR} )
MESSAGE( STATUS "OSVR_LIBRARIES:         " ${OSVR_LIBRARIES} )
MESSAGE( STATUS "OSVR_LIBRARY_DIRS:         " ${OSVR_LIBRARY_DIRS} )
MESSAGE( STATUS "OSVR_INCLUDE_DIR:         " ${OSVR_INCLUDE_DIR} )
MESSAGE( STATUS "EIGEN3_LIBRARIES:         " ${EIGEN3_LIBRARIES} )
MESSAGE( STATUS "ENGEN3_LIBRARY_DIRS:         " ${EIGEN3_LIBRARY_DIRS} )
MESSAGE( STATUS "EIGEN3_INCLUDE_DIR:         " ${EIGEN3_INCLUDE_DIR} )
MESSAGE( STATUS "OPENGL_LIBRARIES:         " ${OPENGL_LIBRARIES} )
MESSAGE( STATUS "OPENGL_gl_LIBRARY:         " ${OPENGL_gl_LIBRARY} )
MESSAGE( STATUS "OPENGL_INCLUDE_DIRS:         " ${OPENGL_INCLUDE_DIRS} )
MESSAGE( STATUS "SDL2_LIBRARIES:         " ${SDL2_LIBRARIES} )
MESSAGE( STATUS "SDL2_INCLUDE_DIRS:         " ${SDL2_INCLUDE_DIRS} )
MESSAGE( STATUS "SOIL_LIBRARIES:         " ${SOIL_LIBRARIES} )
MESSAGE( STATUS "SOIL_INCLUDE_DIRS:         " ${SOIL_INCLUDE_DIRS} )
MESSAGE( STATUS "ZMQ_LIBRARIES:         " ${ZMQ_LIBRARIES} )
MESSAGE( STATUS "ZMQ_INCLUDE_DIRS:         " ${ZMQ_INCLUDE_DIRS} )

target_include_directories(OSVRClientTest PUBLIC
LibAvWrapper/inc
${Boost_INCLUDE_DIR}
${OSVR_INCLUDE_DIR}
${EIGEN3_INCLUDE_DIR}
${OPENGL_INCLUDE_DIRS}
${SDL2_INCLUDE_DIRS}
${SOIL_INCLUDE_DIRS}
${ZMQ_INCLUDE_DIRS}
)

target_link_libraries(OSVRClientTest
    osvr::osvrClientKitCpp
    osvrRenderManager::osvrRenderManagerCpp
    LibAvWrapper
    C:/Program Files (x86)/boost_1_57_0/stage/win32/lib/libboost_program_options-vc120-mt-gd-1_57.lib
    ${Boost_LIBRARIES}
    ${EIGEN3_LIBRARIES}
    ${GLEW_LIBRARIES}
    ${SOIL_LIBRARIES}
    ${SDL2_LIBRARIES}
    ${OPENGL_LIBRARIES}
    ${OSVR_LIBRARIES}
    ${ZMQ_LIBRARIES}
)
