cmake_minimum_required(VERSION 3.25)
project(mldsa-decompose LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# nP is a template parameter, so each party count needs its own binary. The
# Makefile appends whatever `n=` asks for; these are just the handy defaults.
set(MLDSA_PARTIES "2;3;5" CACHE STRING "Party counts to build main<n> for")
set(MLDSA_PARAM   "MLDSA_65" CACHE STRING "Parameter set: MLDSA_44 | MLDSA_65 | MLDSA_87")
set(MLDSA_SIGMA   "128" CACHE STRING "Ring authentication slack sigma (shares/MACs live in Z_2^(nu+sigma))")
option(MLDSA_TEST "Open secret values and run PrepSign self-checks" OFF)
option(MLDSA_SLOT "Use the slot-restricted boundary test (OFF = two-round baseline)" ON)
option(MLDSA_TAMPER_C "Compile the adversarial c-share corruption test" OFF)
list(REMOVE_DUPLICATES MLDSA_PARTIES)

# Pulls emp-tool and emp-ot in as dependencies automatically.
find_package(emp-ag REQUIRED)

# PrepSign step-4 downstream (open c -> A2B -> Decompose) under n-party
# malicious security. One binary per party count (nP is a template parameter).
foreach(_np ${MLDSA_PARTIES})
    add_executable(main${_np} main.cpp)
    target_include_directories(main${_np} PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src/infra
        ${CMAKE_CURRENT_SOURCE_DIR}/src/circuit
        ${CMAKE_CURRENT_SOURCE_DIR}/src/protocol)
    target_compile_definitions(main${_np}
        PRIVATE MLDSA_NP=${_np} MLDSA_PARAM=${MLDSA_PARAM}
                MLDSA_SIGMA=${MLDSA_SIGMA}
                MLDSA_SLOT=$<BOOL:${MLDSA_SLOT}>)
    if(MLDSA_TEST)
        target_compile_definitions(main${_np} PRIVATE TEST)
    endif()
    if(MLDSA_TAMPER_C)
        target_compile_definitions(main${_np} PRIVATE TAMPER_C)
    endif()
    target_link_libraries(main${_np} PRIVATE emp-ag::emp-ag)
endforeach()
