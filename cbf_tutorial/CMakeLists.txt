cmake_minimum_required(VERSION 3.10)
project(cbf_tutorial)

set(CMAKE_CXX_STANDARD 11)

find_package(osqp REQUIRED CONFIG)

if(UNIX)
    find_package(Eigen3 REQUIRED)
    if(${CMAKE_VERSION} VERSION_LESS "3.12.0") 
         find_package(PythonInterp)
         find_package(PythonLibs)
         include_directories(include
                             ${EIGEN3_INCLUDE_DIR}
                             ${PYTHON_INCLUDE_DIRS})
    else()
         find_package(Python COMPONENTS Interpreter Development)
         include_directories(include
                             ${EIGEN3_INCLUDE_DIR}
                             ${Python_INCLUDE_DIRS})
    endif()
elseif(WIN32)
    list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")
    find_package(Eigen3 REQUIRED)

    set(Python_FIND_REGISTRY "LAST")
    find_package(Python COMPONENTS Interpreter Development NumPy)
  
    include_directories(include
                        ${EIGEN3_INCLUDE_DIR}
                        ${Python_INCLUDE_DIRS}
                        ${Python_NumPy_INCLUDE_DIRS})
endif (UNIX)

add_executable(double_integrator src/double_integrator.cpp src/cbf.cpp src/qpwrapper_abstract.cpp src/qpwrapper_osqp.cpp)
add_executable(single_integrator src/single_integrator.cpp src/cbf.cpp src/qpwrapper_abstract.cpp src/qpwrapper_osqp.cpp)
add_executable(single_integrator_2dof src/single_integrator_2dof.cpp src/cbf.cpp src/qpwrapper_abstract.cpp src/qpwrapper_osqp.cpp)
add_executable(manipulator src/manipulator.cpp src/cbf.cpp src/qpwrapper_abstract.cpp src/qpwrapper_osqp.cpp)

if(UNIX)
    if(${CMAKE_VERSION} VERSION_LESS "3.12.0") 
         target_link_libraries(double_integrator PUBLIC osqp::osqpstatic
                               ${PYTHON_LIBRARIES})
         target_link_libraries(single_integrator PUBLIC osqp::osqpstatic
                               ${PYTHON_LIBRARIES})
         target_link_libraries(single_integrator_2dof PUBLIC osqp::osqpstatic
                               ${PYTHON_LIBRARIES})
         target_link_libraries(manipulator PUBLIC osqp::osqpstatic
                               ${PYTHON_LIBRARIES})
    else()
         target_link_libraries(double_integrator PUBLIC osqp::osqpstatic
                               ${Python_LIBRARIES})
         target_link_libraries(single_integrator PUBLIC osqp::osqpstatic
                               ${Python_LIBRARIES})
         target_link_libraries(single_integrator_2dof PUBLIC osqp::osqpstatic
                               ${Python_LIBRARIES})
         target_link_libraries(manipulator PUBLIC osqp::osqpstatic
                               ${Python_LIBRARIES})
    endif()
elseif(WIN32)
    target_link_libraries(double_integrator PUBLIC osqp::osqpstatic
            ${Python_LIBRARIES})
    target_link_libraries(single_integrator PUBLIC osqp::osqpstatic
            ${Python_LIBRARIES})
    target_link_libraries(single_integrator_2dof PUBLIC osqp::osqpstatic
            ${Python_LIBRARIES})
    target_link_libraries(manipulator PUBLIC osqp::osqpstatic
            ${Python_LIBRARIES})
endif (UNIX)
