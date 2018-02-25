# rm *.o *.mod *.exe
rm *.o *.exe
# rm *.exe

gfortran -c icefront_run.f90
gfortran -c icefront_ini_fixed.f90
gfortran -c mod_kinds.f90
gfortran -c mod_param.f90

gfortran *.o -o icefront_run.exe
