# rm *.o *.mod *.exe
rm *.o *.exe
# rm *.exe

gfortran -c read_test.f90
gfortran read_test.o -o read_test.exe
