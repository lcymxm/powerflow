cmake_minimum_required(VERSION 3.16)

project(PowerFlow
    VERSION 1.0.0
    DESCRIPTION "AC power-flow solver: Newton-Raphson, CSV/JSON export, SVG diagrams"
    LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

option(POWERFLOW_BUILD_TESTS "Build the unit test suite" ON)
option(POWERFLOW_WARNINGS_AS_ERRORS "Treat compiler warnings as errors" OFF)

# ---------------------------------------------------------------- warnings ----
function(powerflow_enable_warnings target)
    if(MSVC)
        target_compile_options(${target} PRIVATE /W4 /permissive-)
        if(POWERFLOW_WARNINGS_AS_ERRORS)
            target_compile_options(${target} PRIVATE /WX)
        endif()
    else()
        target_compile_options(${target} PRIVATE -Wall -Wextra -Wpedantic -Wshadow)
        if(POWERFLOW_WARNINGS_AS_ERRORS)
            target_compile_options(${target} PRIVATE -Werror)
        endif()
    endif()
endfunction()

# ------------------------------------------------------------------- core ----
add_library(powerflow_core
    src/case_file.cpp
    src/contingency.cpp
    src/dc_opf.cpp
    src/jacobian.cpp
    src/network.cpp
    src/power_flow.cpp
    src/report.cpp
    src/sparse_lu.cpp
    src/svg_renderer.cpp
    src/synthetic.cpp)

target_include_directories(powerflow_core PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>)

powerflow_enable_warnings(powerflow_core)

# --------------------------------------------------------------- executable ----
add_executable(powerflow src/main.cpp)
target_link_libraries(powerflow PRIVATE powerflow_core)
powerflow_enable_warnings(powerflow)

# -------------------------------------------------------------------- tests ----
include(CTest)
if(POWERFLOW_BUILD_TESTS AND BUILD_TESTING)
    add_executable(powerflow_tests
        tests/test_main.cpp
        tests/test_case_file.cpp
        tests/test_contingency.cpp
        tests/test_dc_opf.cpp
        tests/test_dense.cpp
        tests/test_network.cpp
        tests/test_power_flow.cpp
        tests/test_report.cpp
        tests/test_sparse.cpp)

    target_link_libraries(powerflow_tests PRIVATE powerflow_core)
    target_include_directories(powerflow_tests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/tests)
    target_compile_definitions(powerflow_tests PRIVATE
        POWERFLOW_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/data")
    powerflow_enable_warnings(powerflow_tests)

    add_test(NAME powerflow.unit COMMAND powerflow_tests)

    set(case14 ${CMAKE_CURRENT_SOURCE_DIR}/data/case14.m)
    set(case30 ${CMAKE_CURRENT_SOURCE_DIR}/data/case30.m)
    set(case57 ${CMAKE_CURRENT_SOURCE_DIR}/data/case57.m)

    add_test(NAME powerflow.cli.help COMMAND powerflow --help)
    add_test(NAME powerflow.cli.version COMMAND powerflow --version)
    add_test(NAME powerflow.cli.ieee14 COMMAND powerflow ${case14} --no-tables)
    add_test(NAME powerflow.cli.ieee30_json COMMAND powerflow ${case30} --no-tables
             --json ${CMAKE_CURRENT_BINARY_DIR}/case30.json)
    add_test(NAME powerflow.cli.ieee57_csv COMMAND powerflow ${case57} --no-tables
             --csv ${CMAKE_CURRENT_BINARY_DIR}/csv_out)
    add_test(NAME powerflow.cli.ieee14_svg COMMAND powerflow ${case14} --no-tables
             --svg ${CMAKE_CURRENT_BINARY_DIR}/case14.svg)
    add_test(NAME powerflow.cli.compare COMMAND powerflow ${case14} --method both --no-tables)
    add_test(NAME powerflow.cli.dense COMMAND powerflow ${case14} --no-tables --linear dense)
    add_test(NAME powerflow.cli.n1 COMMAND powerflow ${case14} --n1 --n1-worst 3)
    # The scan deliberately exits non-zero when it finds a contingency that breaks
    # a limit (case14 does), so this smoke test only checks that it produced a
    # report rather than judging the network's security.
    set_tests_properties(powerflow.cli.n1 PROPERTIES PASS_REGULAR_EXPRESSION
        "N-1 contingency scan")
    add_test(NAME powerflow.cli.opf COMMAND powerflow ${case30} --opf)
    add_test(NAME powerflow.cli.benchmark COMMAND powerflow --benchmark 120)
endif()

# ------------------------------------------------------------------ install ----
install(TARGETS powerflow RUNTIME DESTINATION bin)
install(DIRECTORY include/powerflow DESTINATION include)
install(DIRECTORY data/ DESTINATION share/powerflow/data)
