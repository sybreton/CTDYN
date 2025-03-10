find_package(LAPACK)
find_package(BLAS)

add_library(kind_parameters kind_parameters.f90)
add_library(util util.f90)
add_library(cio cio.f90)
add_library(b2func b2func.f90)
add_library(b3func b3func.f90)
add_library(b4func b4func.f90)
add_library(b5func b5func.f90)
add_library(c1s2func c1s2func.f90)
add_library(c3s2func c3s2func.f90)
add_library(plg plg.f90)
add_library(bessel bessel.f90)
add_library(cdfunc cdfunc.f90)
add_library(stellar_profiles stellar_profiles.f90)
add_library(func_flow func_flow.f90)
add_library(write_outputs write_outputs.f90)
add_library(dyna dyna.f90)
add_library(zbrent zbrent.f90)

add_executable(ctdyn "")

target_link_libraries(ctdyn LAPACK::LAPACK -llapack)
target_link_libraries(ctdyn BLAS::BLAS -lblas)


target_sources(ctdyn PRIVATE
               kind_parameters.f90
               util.f90
               cio.f90
               b2func.f90
               b3func.f90
               b4func.f90
               b5func.f90
               c1s2func.f90
               c3s2func.f90
               plg.f90
               bessel.f90
               cdfunc.f90
               stellar_profiles.f90
               func_flow.f90
               write_outputs.f90
               dyna.f90
               zbrent.f90
               main.f90)
