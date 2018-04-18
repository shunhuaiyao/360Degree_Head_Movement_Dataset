set (CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

LIST(APPEND CMAKE_MODULE_PATH "." "..")

find_package( FFMPEG REQUIRED )
MESSAGE( STATUS "FFMPEG found:         " ${FFMPEG_FOUND} )

set(GLEW_INCLUDE_DIR "C:/Program Files (x86)/glew-2.1.0/include")
set(GLEW_LIBRARY "C:/Program Files (x86)/glew-2.1.0/lib/Release/Win32/glew32.lib")
#set(OPENGL_gl_LIBRARY "C:/Program Files (x86)/Windows Kits/10/Lib/10.0.16299.0/um/x86/OpenGL32.Lib")
#set(OPENGL_glu_LIBRARY "C:/Program Files (x86)/Windows Kits/10/Lib/10.0.16299.0/um/x86/GlU32.Lib")
#set(OPENGL_INCLUDE_DIRS "C:/Program Files (x86)/Windows Kits/10/Include/10.0.16299.0/um")

find_package(GLEW)
find_package(OpenGL)
# if (${OpenCV_FOUND} AND ${GLEW_FOUND})
#   message(STATUS "Set USE_OPENGl")
#   add_definitions(-DUSE_OPENGl)
# endif()

set(SDL2_LIBRARIES "C:/Program Files (x86)/SDL2-2.0.7/lib/x86/SDL2.lib")
set(SDL2_INCLUDE_DIR "C:/Program Files (x86)/SDL2-2.0.7/include")
find_package(SDL2)

MESSAGE( STATUS "OPENGL_LIBRARIES:         " ${OPENGL_LIBRARIES} )
MESSAGE( STATUS "OPENGL_gl_LIBRARY:         " ${OPENGL_gl_LIBRARY} )
MESSAGE( STATUS "OPENGL_glu_LIBRARY:         " ${OPENGL_glu_LIBRARY} )
MESSAGE( STATUS "OPENGL_INCLUDE_DIRS:         " ${OPENGL_INCLUDE_DIRS} )
MESSAGE( STATUS "OpenCV_FOUND:         " ${OpenCV_FOUND} )

FILE(GLOB LibAvWrapperSrc src/*.cpp)

add_library(LibAvWrapper ${LibAvWrapperSrc})
target_include_directories(LibAvWrapper PUBLIC inc ${FFMPEG_INCLUDE_DIR} ${OPENGL_INCLUDE_DIRS} ${GLEW_INCLUDE_DIRS} ${SDL2_INCLUDE_DIR})
#target_compile_features(LibAvWrapper PRIVATE cxx_range_for cxx_nullptr cxx_auto_type)
target_link_libraries( LibAvWrapper -pthread  ${FFMPEG_LIBRARIES} avresample ${OPENGL_LIBRARIES} ${GLEW_LIBRARIES})

MESSAGE( STATUS "FFMPEG_LIBRARIES:         " ${FFMPEG_LIBRARIES} )
MESSAGE( STATUS "OPENGL_LIBRARIES:         " ${OPENGL_LIBRARIES} )
MESSAGE( STATUS "GLEW_LIBRARIES:         " ${GLEW_LIBRARIES} )
