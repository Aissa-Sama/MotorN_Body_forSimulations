cmake_minimum_required(VERSION 3.15)
project(nbody_sim)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(MSVC)
    add_compile_options(/fp:precise)
    add_compile_options(/fp:except-)
else()
    add_compile_options(-ffp-contract=off)
    add_compile_options(-fno-unsafe-math-optimizations)
    add_compile_options(-fno-associative-math)
endif()

include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/core
    ${CMAKE_CURRENT_SOURCE_DIR}/analysis
    ${CMAKE_CURRENT_SOURCE_DIR}/integrators
    ${CMAKE_CURRENT_SOURCE_DIR}/io
    ${CMAKE_CURRENT_SOURCE_DIR}/regularization
    ${CMAKE_CURRENT_SOURCE_DIR}/regularization/ks
    ${CMAKE_CURRENT_SOURCE_DIR}/regularization/chain
    ${CMAKE_CURRENT_SOURCE_DIR}/regularization/hierarchy
)

file(GLOB LIB_SOURCES
    "analysis/*.cpp"
    "integrators/*.cpp"
    "io/*.cpp"
    "regularization/ks/*.cpp"
    "regularization/chain/*.cpp"
    "initial_conditions.cpp"
    "core/units.cpp"
    "core/nbody_system.cpp"
    "regularization/hierarchy/*.cpp"
    "regularization/chain/chain3_bs_integrator.cpp"
)

add_library(nbody_lib STATIC ${LIB_SOURCES})

add_executable(nbody_sim main.cpp)
target_link_libraries(nbody_sim PRIVATE nbody_lib)

# ── Tests ────────────────────────────────────────────────────────────────────
foreach(_t
    test_ks
    test_chain_ks
    test_chain3_state
    test_chain3_integration
    test_chain3_diagnostic
    test_hybrid_chain3
    test_hierarchical
    test_archain3
    test_archain3_bs
    test_chain3_bs
    test_bs_vs_rk4
    test_hierarchical_pn
    test_archain_n_ks
    test_gw_observables
    test_physical_benchmarks
    test_validation_benchmarks
    test_archain_n
    test_archain_n_pn
    test_block_timestep
)
    if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/tests/${_t}.cpp)
        add_executable(${_t} tests/${_t}.cpp)
        target_link_libraries(${_t} PRIVATE nbody_lib)
    endif()
endforeach()

# ── Python binding ────────────────────────────────────────────────────────────
if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/extern/pybind11/CMakeLists.txt
   AND EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/python/nbody_bindings.cpp)

    add_subdirectory(extern/pybind11)
    pybind11_add_module(nbody_core python/nbody_bindings.cpp)
    target_include_directories(nbody_core PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}
        ${CMAKE_CURRENT_SOURCE_DIR}/core
        ${CMAKE_CURRENT_SOURCE_DIR}/analysis
        ${CMAKE_CURRENT_SOURCE_DIR}/integrators
        ${CMAKE_CURRENT_SOURCE_DIR}/io
        ${CMAKE_CURRENT_SOURCE_DIR}/regularization
        ${CMAKE_CURRENT_SOURCE_DIR}/regularization/ks
        ${CMAKE_CURRENT_SOURCE_DIR}/regularization/chain
        ${CMAKE_CURRENT_SOURCE_DIR}/regularization/hierarchy
    )
    target_link_libraries(nbody_core PRIVATE nbody_lib)
    if(MSVC)
        target_compile_options(nbody_core PRIVATE /W4 /utf-8 /wd4100 /wd4505)
    endif()
    message(STATUS "Target nbody_core.pyd configurado (pybind11 local)")
else()
    message(STATUS "extern/pybind11 o python/nbody_bindings.cpp no encontrado — binding omitido")
endif()

# ── Warnings ─────────────────────────────────────────────────────────────────
set(ALL_TEST_TARGETS
    test_ks test_chain_ks test_chain3_state test_chain3_integration
    test_chain3_diagnostic test_hybrid_chain3 test_hierarchical
    test_archain3 test_archain3_bs test_chain3_bs test_bs_vs_rk4
    test_validation_benchmarks test_archain_n test_archain_n_pn
    test_hierarchical_pn test_archain_n_ks test_gw_observables
    test_physical_benchmarks test_block_timestep
)

if(MSVC)
    target_compile_options(nbody_lib PRIVATE /W4)
    target_compile_options(nbody_sim PRIVATE /W4 /utf-8)
    foreach(_t ${ALL_TEST_TARGETS})
        if(TARGET ${_t})
            target_compile_options(${_t} PRIVATE /W4 /utf-8)
        endif()
    endforeach()
else()
    target_compile_options(nbody_lib PRIVATE -Wall -Wextra)
    target_compile_options(nbody_sim PRIVATE -Wall -Wextra)
    foreach(_t ${ALL_TEST_TARGETS})
        if(TARGET ${_t})
            target_compile_options(${_t} PRIVATE -Wall -Wextra)
        endif()
    endforeach()
endif()
