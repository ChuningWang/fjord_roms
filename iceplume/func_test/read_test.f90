program test_read

   implicit none
   integer, parameter :: r8 = SELECTED_REAL_KIND(12,300)  ! 64-bit
   integer, parameter :: VLI_K = selected_int_kind (18)
   integer, parameter :: DR_K = selected_real_kind (14)
   real(r8), parameter :: PI = 4.0d0*atan(1.0d0)

   integer(VLI_K) :: i
   real(DR_K) :: a, b, c

   open(unit=15, file='iceplume_test_input.txt', action='read')

   100 format(f12.6, f12.6, f12.6)
   200 format(f100.90)
   read (15, 100)  a, b, c
   write (*, 100)  a, b, c
   write (*, 200)  PI

end program test_read
