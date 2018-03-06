rm *.o *.mod *.exe

gfortran -c mod_param_iceplume.f90
gfortran -c iceplume_meltrate.f90
gfortran -c iceplume_plume_model.f90
gfortran -c iceplume_calc.f90
gfortran -c iceplume_write.f90
gfortran -c opkda1.F opkda2.F opkdmain.F
gfortran -c iceplume_run.f90

gfortran *.o -o iceplume_test.exe

rm *.o *.mod
