PROJECT(exec)
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
CMAKE_POLICY(SET CMP0015 NEW)

add_executable(PPPMain PPPMain.cc)
target_link_libraries(PPPMain CoreLib)

add_executable(test_unit test_unit.cc)
target_link_libraries(test_unit CoreLib)
