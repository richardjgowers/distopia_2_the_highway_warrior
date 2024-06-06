
# option to link to external HIGHWAY library
if (DISTOPIA_USE_EXTERNAL_HWY)
    message(STATUS "Using external HWY library")
    find_package(HWY 1.2.0 REQUIRED)

else()
    message(STATUS "Using internal HWY library")
    set(HWY_ENABLE_TESTS OFF CACHE BOOL "Enable HWY tests" FORCE)
    add_subdirectory("highway")
endif()

set(BENCHMARK_ENABLE_TESTING OFF CACHE BOOL "Enable googlebench tests" FORCE)
add_subdirectory("googlebench")

set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)


add_library(distopia)
target_sources(distopia PRIVATE "src/distopia.cpp")
target_link_libraries(distopia PRIVATE hwy)
target_include_directories(distopia PUBLIC
        "highway"
        ${CMAKE_CURRENT_LIST_DIR}
        "${CMAKE_CURRENT_LIST_DIR}/include")

if (DISTOPIA_BUILD_BENCHMARKS)
    add_executable(bench)
    target_sources(bench PRIVATE "test/bench.cpp")
    target_link_libraries(bench PUBLIC benchmark::benchmark)
    target_link_libraries(bench PUBLIC distopia)
    target_include_directories(bench PUBLIC  ${CMAKE_CURRENT_LIST_DIR})
endif()

if (DISTOPIA_BUILD_TESTS)
    Include(GoogleTest)
    add_subdirectory("googletest")
    enable_testing()
    include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})
    
    add_executable(test)
    target_sources(test PRIVATE "test/test.cpp")
    target_link_libraries(test PUBLIC gtest gtest_main)
    target_link_libraries(test PUBLIC distopia)
    target_include_directories(test PUBLIC  ${CMAKE_SOURCE_DIR})
    
    add_executable(test_mda_match)
    target_sources(test_mda_match PRIVATE "test/test_mda_match.cpp")
    target_link_libraries(test_mda_match PUBLIC gtest gtest_main)
    target_link_libraries(test_mda_match PUBLIC distopia)
    target_include_directories(test_mda_match PUBLIC  ${CMAKE_SOURCE_DIR})
endif()