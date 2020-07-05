PROJECT(PPPLib_v1.0)
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
CMAKE_POLICY(SET CMP0015 NEW)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_CXX_FLAGS "-std=c++11 -fPIC")
set(CMAKE_CXX_COMPILE "g++")
set(CMAKE_C_COMPILE "gcc")
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)
set(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)

include_directories(${PROJECT_SOURCE_DIR}/include
                    ${PROJECT_SOURCE_DIR}/include/Core
                    ${PROJECT_SOURCE_DIR}/include/Core/Logger
                    ${PROJECT_SOURCE_DIR}/include/Core/Config
                    ${PROJECT_SOURCE_DIR}/include/Core/CmnFunc
                    ${PROJECT_SOURCE_DIR}/include/Core/ReadFiles
                    ${PROJECT_SOURCE_DIR}/include/Core/ParSetting
                    ${PROJECT_SOURCE_DIR}/include/Core/Solver
                    ${PROJECT_SOURCE_DIR}/include/Core/AdjModel
                    ${PROJECT_SOURCE_DIR}/include/Core/QuaCtrl
                    ${PROJECT_SOURCE_DIR}/include/Core/OutSolution
                    ${PROJECT_SOURCE_DIR}/include/GNSS
                    ${PROJECT_SOURCE_DIR}/include/GNSS/AmbModel
                    ${PROJECT_SOURCE_DIR}/include/GNSS/AntModel
                    ${PROJECT_SOURCE_DIR}/include/GNSS/EphModel
                    ${PROJECT_SOURCE_DIR}/include/GNSS/IonModel
                    ${PROJECT_SOURCE_DIR}/include/GNSS/BiaModel
                    ${PROJECT_SOURCE_DIR}/include/GNSS/TidModel
                    ${PROJECT_SOURCE_DIR}/include/GNSS/TrpModel
                    ${PROJECT_SOURCE_DIR}/include/GNSS/AttModel)


add_subdirectory(${PROJECT_SOURCE_DIR}/src/Core)
add_subdirectory(${PROJECT_SOURCE_DIR}/src/GNSS)
add_subdirectory(${PROJECT_SOURCE_DIR}/exec)




