rm *.o *.mod *.exe
# rm *.o *.exe
# rm *.exe

gfortran -c mod_param_iceplume.f90
gfortran -c iceplume_run.f90
gfortran -c opkda1.F opkda2.F opkdmain.F

gfortran *.o -o iceplume_test.exe
