CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

set (CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

IF(NOT CMAKE_BUILD_TYPE)
  # SET(CMAKE_BUILD_TYPE "DEBUG")
  SET(CMAKE_BUILD_TYPE "RELEASE")
  #SET(CMAKE_BUILD_TYPE "RELWITHDEBINFO")
  #SET(CMAKE_BUILD_TYPE "MINSIZEREL")
ENDIF()

find_package(OpenMP QUIET)
if(OPENMP_FOUND)
   MESSAGE( STATUS "xc: OpenMP found" )
   set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
   set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

# Find necessary packages
FIND_PACKAGE( PythonLibs 3.4 REQUIRED )

FIND_PACKAGE(Boost 1.57.0)
IF(Boost_FOUND)

  set(Boost_INCLUDE_DIRS "C:/Program Files (x86)/boost_1_57_0")
  set(PYTHON_INCLUDE_DIRS "C:/Users/ShunHuaiYao/AppData/Local/Programs/Python/Python34/include")

  INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS} ${PYTHON_INCLUDE_DIRS})
  SET(Boost_USE_STATIC_LIBS OFF)
  SET(Boost_USE_MULTITHREADED ON)
  SET(Boost_USE_STATIC_RUNTIME OFF)
  FIND_PACKAGE(Boost 1.57.0 COMPONENTS python3)

  set(Boost_LIBRARIES "C:/Program Files (x86)/boost_1_57_0/stage/win32/lib/boost_python-vc140-mt-1_57.lib")
  set(PYTHON_LIBRARIES "C:/Users/ShunHuaiYao/AppData/Local/Programs/Python/Python34/libs/python34.lib")
  
  MESSAGE( STATUS "Boost_LIBRARIES:         " ${Boost_LIBRARIES} )
  MESSAGE( STATUS "Boost_INCLUDE_DIRS:      " ${Boost_INCLUDE_DIRS} )
  MESSAGE( STATUS "PYTHON_LIBRARIES:        " ${PYTHON_LIBRARIES} )
  MESSAGE( STATUS "PYTHON_INCLUDE_DIRS:     " ${PYTHON_INCLUDE_DIRS} )
  
  ADD_LIBRARY(CQuaternion SHARED CQuaternion.cpp)
  #TARGET_LINK_LIBRARIES(CQuaternion ${Boost_LIBRARIES} "boost_python3" ${PYTHON_LIBRARIES})
  TARGET_LINK_LIBRARIES(CQuaternion ${Boost_LIBRARIES} ${PYTHON_LIBRARIES})
  set_target_properties( CQuaternion PROPERTIES PREFIX "" )
  #add_custom_command(TARGET CQuaternion POST_BUILD COMMAND cp CQuaternion.dll ..)
ELSEIF(NOT Boost_FOUND)
  MESSAGE(FATAL_ERROR "Unable to find correct Boost version. Did you set BOOST_ROOT?")
ENDIF()

#IF(CMAKE_COMPILER_IS_GNUCXX)
#  ADD_DEFINITIONS("-Wall")
#ELSE()
#  MESSAGE(FATAL_ERROR "CMakeLists.txt has not been tested/written for your compiler.")
#ENDIF()
