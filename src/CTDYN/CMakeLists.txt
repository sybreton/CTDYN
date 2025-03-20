find_package(fortran_stdlib REQUIRED)

# Default search for external BLAS/LAPACK library
# (e.g. openBLAS)
if ($ENV{USE_DEFAULT_EXTERNAL} EQUAL "1")
    find_package(BLAS 
                 REQUIRED
                 )
    find_package(LAPACK 
                 REQUIRED
                 )
    add_compile_definitions(STDLIB_EXTERNAL_BLAS STDLIB_EXTERNAL_LAPACK)
# Attempt at using LAPACK/BLAS version from the MESA SDK if available and
# USE_MESASDK is set to 1
elseif ($ENV{USE_MESASDK} EQUAL "1" AND EXISTS "$ENV{MESASDK_ROOT}")
    set(BLAS_DIR "$ENV{MESASDK_ROOT}/math-slots/default/lib")
    set(LAPACK_DIR "$ENV{MESASDK_ROOT}/math-slots/default/lib")
    find_library(NAMES BLAS 
                PATHS ${BLAS_DIR}
                REQUIRED
                NO_DEFAULT_PATH 
                NO_CMAKE_FIND_ROOT_PATH
                )
    find_library(NAMES LAPACK
                PATHS ${LAPACK_DIR}
                REQUIRED
                NO_DEFAULT_PATH
                NO_CMAKE_FIND_ROOT_PATH
                )
    add_compile_definitions(STDLIB_EXTERNAL_BLAS STDLIB_EXTERNAL_LAPACK)
endif()

add_executable(ctdyn "")

if ($ENV{USE_DEFAULT_EXTERNAL} OR $ENV{USE_MESASDK} EQUAL "1")
    target_link_libraries(
        ctdyn  
        PRIVATE
        fortran_stdlib::fortran_stdlib
    )
    target_link_libraries(
        ctdyn  
        PRIVATE
        BLAS::BLAS -lblas
    )
    target_link_libraries(
        ctdyn  
        PRIVATE
        LAPACK::LAPACK -llapack
    )
else ()
    target_link_libraries(
        ctdyn  
        PUBLIC
        fortran_stdlib::fortran_stdlib
    )
endif()

target_sources(
    ctdyn 
    PRIVATE
    kind_parameters.f90
    mathfunc.f90
    util.f90
    cio.f90
    stellar_profiles.f90
    write_outputs.f90
    dyna.f90
    zbrent.f90
    main.f90
)
