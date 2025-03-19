find_package(fortran_stdlib REQUIRED)

add_executable(ctdyn "")

target_link_libraries(
    ctdyn  
    PUBLIC
    fortran_stdlib::fortran_stdlib
)

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
