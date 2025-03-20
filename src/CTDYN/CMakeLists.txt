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
endif()

add_executable(ctdyn "")

if ($ENV{USE_DEFAULT_EXTERNAL} EQUAL "1")
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
        -llapack
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
