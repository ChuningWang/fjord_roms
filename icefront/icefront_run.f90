! Test the icefront model.

PROGRAM icefront
    USE mod_kinds
    USE mod_param

	implicit none

	print *, "IceFront."

    CALL icefront_init_fixed()

END PROGRAM icefront


!     SUBROUTINE icefront_init_fixed()
! 
! ! =========================================================
! !
! ! Initialize ICEFRONT parameters and variables.
! !
! ! =========================================================
! !         USE mod_kinds
! !         USE mod_param
! 
!         implicit none
!         
!         integer :: iter
! 
!     END SUBROUTINE icefront_init_fixed
