rm *.o *.mod *.exe

gfortran -c mod_kinds.f90
gfortran -c mod_iceplume.f90
gfortran -c opkd.f90
gfortran -c iceplume_write.f90
gfortran -c iceplume.f90

gfortran *.o -o iceplume_test.exe

rm *.o *.mod
