cmake_minimum_required(VERSION 3.15)
project(nbody_sim)

# ============================================================================
# CONFIGURACIÓN BÁSICA DEL LENGUAJE
# ============================================================================
set(CMAKE_CXX_STANDARD 17)           # Usar C++17
set(CMAKE_CXX_STANDARD_REQUIRED ON)  # Es obligatorio

# ============================================================================
# 1. CREAR UNA LIBRERÍA CON TODO EL CÓDIGO (excepto main)
# ============================================================================

# Directorios de include (para que los #include funcionen)
include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}                # Raíz del proyecto
    ${CMAKE_CURRENT_SOURCE_DIR}/core           # Vec3, Vec4, Body, NBodySystem
    ${CMAKE_CURRENT_SOURCE_DIR}/analysis       # BinaryDetector, BinaryState
    ${CMAKE_CURRENT_SOURCE_DIR}/integrators    # HybridIntegrator, RK4, etc.
    ${CMAKE_CURRENT_SOURCE_DIR}/io             # Loggers, Snapshots
    ${CMAKE_CURRENT_SOURCE_DIR}/regularization # Base de regularización
    ${CMAKE_CURRENT_SOURCE_DIR}/regularization/ks  # KS regularización
    ${CMAKE_CURRENT_SOURCE_DIR}/regularization/chain # Chain regularization
    ${CMAKE_CURRENT_SOURCE_DIR}/regularization/hierarchy # HierarchyBuilder, HierarchyNode
)

# Recopilar todos los archivos fuente (.cpp) EXCLUYENDO main.cpp y tests
# GLOB significa "global" - busca todos los archivos que coincidan
file(GLOB LIB_SOURCES 
    "analysis/*.cpp"                # Todos los .cpp en analysis/
    "integrators/*.cpp"              # Todos los .cpp en integrators/
    "io/*.cpp"                       # Todos los .cpp en io/
    "regularization/ks/*.cpp"        # Todos los .cpp en regularization/ks/
    "regularization/chain/*.cpp"     # Todos los .cpp en regularization/chain/
    "initial_conditions.cpp"         # Archivo específico
    "core/units.cpp"                  # Archivo específico
    "core/nbody_system.cpp"
    "regularization/hierarchy/*.cpp"   # HierarchyBuilder
    "regularization/chain/chain3_bs_integrator.cpp"
)

# Crear una librería ESTÁTICA (un archivo .lib) con todo el código fuente
# Esto permite que múltiples ejecutables (tests) compartan el mismo código
add_library(nbody_lib STATIC ${LIB_SOURCES})

# ============================================================================
# 2. CREAR EL EJECUTABLE PRINCIPAL
# ============================================================================
add_executable(nbody_sim main.cpp)                    # El simulador principal
target_link_libraries(nbody_sim PRIVATE nbody_lib)    # Enlaza con la librería

# ============================================================================
# 3. CREAR LOS EJECUTABLES DE TEST
# ============================================================================
# Cada test es un ejecutable independiente que enlaza con la misma librería
# El IF(EXISTS ...) verifica que el archivo existe antes de intentar compilarlo

# Test de KS original (binarias)
if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_ks.cpp)
    add_executable(test_ks tests/test_ks.cpp)
    target_link_libraries(test_ks PRIVATE nbody_lib)
endif()

# Test de transformaciones KS para cadena (funciones básicas)
if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_chain_ks.cpp)
    add_executable(test_chain_ks tests/test_chain_ks.cpp)
    target_link_libraries(test_chain_ks PRIVATE nbody_lib)
endif()

# Test de transformaciones cadena ↔ físico (initialize/write_back)
if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_chain3_state.cpp)
    add_executable(test_chain3_state tests/test_chain3_state.cpp)
    target_link_libraries(test_chain3_state PRIVATE nbody_lib)
endif()

# Test de integración de cadena (figura-8 completa)
if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_chain3_integration.cpp)
    add_executable(test_chain3_integration tests/test_chain3_integration.cpp)
    target_link_libraries(test_chain3_integration PRIVATE nbody_lib)
endif()

# 🔴 NUEVO: Test de diagnóstico para chain regularization (A, B, C)
if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_chain3_diagnostic.cpp)
    add_executable(test_chain3_diagnostic tests/test_chain3_diagnostic.cpp)
    target_link_libraries(test_chain3_diagnostic PRIVATE nbody_lib)
endif()

if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_hybrid_chain3.cpp)
    add_executable(test_hybrid_chain3 tests/test_hybrid_chain3.cpp)
    target_link_libraries(test_hybrid_chain3 PRIVATE nbody_lib)
    if(MSVC)
        target_compile_options(test_hybrid_chain3 PRIVATE /W4)
    else()
        target_compile_options(test_hybrid_chain3 PRIVATE -Wall -Wextra)
    endif()
endif()

if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_hierarchical.cpp)
    add_executable(test_hierarchical tests/test_hierarchical.cpp)
    target_link_libraries(test_hierarchical PRIVATE nbody_lib)
    if(MSVC)
        target_compile_options(test_hierarchical PRIVATE /W4)
    else()
        target_compile_options(test_hierarchical PRIVATE -Wall -Wextra)
    endif()
endif()

if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_chain3_bs.cpp)
    add_executable(test_chain3_bs tests/test_chain3_bs.cpp)
    target_link_libraries(test_chain3_bs PRIVATE nbody_lib)
    if(MSVC)
        target_compile_options(test_chain3_bs PRIVATE /W4)
    else()
        target_compile_options(test_chain3_bs PRIVATE -Wall -Wextra)
    endif()

if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_bs_vs_rk4.cpp)
    add_executable(test_bs_vs_rk4 tests/test_bs_vs_rk4.cpp)
    target_link_libraries(test_bs_vs_rk4 PRIVATE nbody_lib)
    if(MSVC)
        target_compile_options(test_bs_vs_rk4 PRIVATE /utf-8)
    endif()
endif()
endif()

# ============================================================================
# 4. OPCIONES DE COMPILACIÓN (SEPARADAS POR COMPILADOR)
# ============================================================================
if(MSVC)
    # ========================================================================
    # OPCIONES PARA MICROSOFT VISUAL C++ (MSVC)
    # ========================================================================
    # /W4 = advertencias nivel 4 (las más estrictas pero razonables)
    
    # Opciones para la librería
    target_compile_options(nbody_lib PRIVATE /W4)
    
    # Opciones para el ejecutable principal
    target_compile_options(nbody_sim PRIVATE /W4)
    
    # Opciones para cada test (solo si el objetivo TARGET existe)
    if(TARGET test_ks)
        target_compile_options(test_ks PRIVATE /W4)
    endif()
    
    if(TARGET test_chain_ks)
        target_compile_options(test_chain_ks PRIVATE /W4)
    endif()
    
    if(TARGET test_chain3_state)
        target_compile_options(test_chain3_state PRIVATE /W4)
    endif()
    
    if(TARGET test_chain3_integration)
        target_compile_options(test_chain3_integration PRIVATE /W4)
    endif()
    
    # 🔴 NUEVO: Opciones para el test de diagnóstico
    if(TARGET test_chain3_diagnostic)
        target_compile_options(test_chain3_diagnostic PRIVATE /W4)
    endif()

else()
    # ========================================================================
    # OPCIONES PARA COMPILADORES GCC/CLANG (Linux/macOS)
    # ========================================================================
    # -Wall = todas las advertencias
    # -Wextra = advertencias extra
    
    target_compile_options(nbody_lib PRIVATE -Wall -Wextra)
    target_compile_options(nbody_sim PRIVATE -Wall -Wextra)
    
    if(TARGET test_ks)
        target_compile_options(test_ks PRIVATE -Wall -Wextra)
    endif()
    
    if(TARGET test_chain_ks)
        target_compile_options(test_chain_ks PRIVATE -Wall -Wextra)
    endif()
    
    if(TARGET test_chain3_state)
        target_compile_options(test_chain3_state PRIVATE -Wall -Wextra)
    endif()
    
    if(TARGET test_chain3_integration)
        target_compile_options(test_chain3_integration PRIVATE -Wall -Wextra)
    endif()
    
    # 🔴 NUEVO: Opciones para GCC/Clang
    if(TARGET test_chain3_diagnostic)
        target_compile_options(test_chain3_diagnostic PRIVATE -Wall -Wextra)
    endif()
endif()

# ============================================================================
# NOTA: Los tests se compilarán automáticamente cuando hagas:
#   cmake --build .
#
# Para ejecutar un test específico:
#   .\test_chain3_diagnostic.exe
#   .\test_chain3_integration.exe
#   etc.
# ============================================================================