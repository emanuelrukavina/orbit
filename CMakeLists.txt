cmake_minimum_required(VERSION 3.16)
project(liborbit VERSION 1.0.0 LANGUAGES CXX)

add_library(orbit STATIC
    src/liborbit.cpp
    src/math.cpp
    src/fit.cpp
    src/statistics.cpp
)
add_library(liborbit::orbit ALIAS orbit)
target_compile_features(orbit PUBLIC cxx_std_17)
target_include_directories(orbit
    PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
           $<INSTALL_INTERFACE:include>
    PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)
set_target_properties(orbit PROPERTIES OUTPUT_NAME orbit)

add_executable(orbit_example examples/main.cpp)
target_link_libraries(orbit_example PRIVATE liborbit::orbit)

include(CTest)
if(BUILD_TESTING)
    add_executable(orbit_tests tests/test_orbit.cpp)
    target_link_libraries(orbit_tests PRIVATE liborbit::orbit)
    add_test(NAME orbit_tests COMMAND orbit_tests)
endif()

install(TARGETS orbit ARCHIVE DESTINATION lib)
install(FILES include/liborbit.hpp DESTINATION include)
