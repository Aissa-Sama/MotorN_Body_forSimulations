cmake_minimum_required(VERSION 3.15)
project(nbody_sim)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# ============================================================================
# 1. CREAR UNA LIBRERÍA CON TODO EL CÓDIGO (excepto main)
# ============================================================================

# Directorios de include
include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/core
    ${CMAKE_CURRENT_SOURCE_DIR}/analysis
    ${CMAKE_CURRENT_SOURCE_DIR}/integrators
    ${CMAKE_CURRENT_SOURCE_DIR}/io
    ${CMAKE_CURRENT_SOURCE_DIR}/regularization
    ${CMAKE_CURRENT_SOURCE_DIR}/regularization/ks
)

# Recopilar todos los archivos fuente (EXCLUYENDO main.cpp y tests)
file(GLOB LIB_SOURCES 
    "analysis/*.cpp"
    "integrators/*.cpp"
    "io/*.cpp"
    "regularization/ks/*.cpp"
    "initial_conditions.cpp"
    "core/units.cpp"
)

# Crear una librería ESTÁTICA con todo el código
add_library(nbody_lib STATIC ${LIB_SOURCES})

# ============================================================================
# 2. CREAR EL EJECUTABLE PRINCIPAL
# ============================================================================
add_executable(nbody_sim main.cpp)
target_link_libraries(nbody_sim PRIVATE nbody_lib)

# ============================================================================
# 3. CREAR EL EJECUTABLE DE TEST (si existe)
# ============================================================================
if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_ks.cpp)
    add_executable(test_ks tests/test_ks.cpp)
    target_link_libraries(test_ks PRIVATE nbody_lib)
endif()

# Opciones de compilación
if(MSVC)
    target_compile_options(nbody_lib PRIVATE /W4)
    target_compile_options(nbody_sim PRIVATE /W4)
    if(TARGET test_ks)
        target_compile_options(test_ks PRIVATE /W4)
    endif()
else()
    target_compile_options(nbody_lib PRIVATE -Wall -Wextra)
    target_compile_options(nbody_sim PRIVATE -Wall -Wextra)
    if(TARGET test_ks)
        target_compile_options(test_ks PRIVATE -Wall -Wextra)
    endif()
endif()