! Test the icefront model.

SUBROUTINE iceplume()

    USE mod_grid
    USE mod_ocean
    USE mod_stepping

    ! ==================================================================
    ! read in profiles at the grid cell

    CALL iceplume_calc
    ! CALL iceplume_write(fileName)

END
