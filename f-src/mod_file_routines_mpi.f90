!------------------------------------------------------------------------------
! MODULE: file routines 
!------------------------------------------------------------------------------
!> @author Johannes Gebert - HLRS - NUM - gebert@hlrs.de
!
! DESCRIPTION: 
!> Module containing file I/O
!------------------------------------------------------------------------------
MODULE file_routines_mpi

USE MPI
USE ISO_FORTRAN_ENV
USE global_std
USE strings

IMPLICIT NONE

CONTAINS

!------------------------------------------------------------------------------
! SUBROUTINE: mpi_err
!------------------------------------------------------------------------------  
!> @author Ralf Schneider - HLRS - NUM - schneider@hlrs.de
!
!> @brief
!> Evaluates allocation errors
!
!> @param[in] ierr Errorcode 
!> @param[out] text Message to print
!------------------------------------------------------------------------------  
subroutine mpi_err(ierr, text)

   !-- Dummy parameters
   integer(kind=mik), intent(in) :: ierr
   character(len=*), intent(in) :: text
   
   if (ierr /= MPI_SUCCESS) THEN
      write(*, "(100('!'))")
      write(*, '(A,I0,A)') 'MPI ERROR :', ierr,";"
      write(*, '(A)') trim(text)
      write(*, "(100('!'))")
      write(*, *) 'Exit ...'
      stop
    end if
   
end subroutine mpi_err

!------------------------------------------------------------------------------
! SUBROUTINE: write_histo_csv
!------------------------------------------------------------------------------  
!> @author Johannes Gebert - HLRS - NUM - gebert@hlrs.de
!
!> @brief
!> Write the csv file of the histogram
!
!> @description
!> !!! Routine needs extensive rework due to meta file format
!
!> @param[in] fh File handle
!> @param[in] filename File name
!> @param[in] hbnds Histogram boundaries
!> @param[in] mov_avg_width Width of the moving average
!> @param[in] histogram Actual histogram data
!------------------------------------------------------------------------------  
 SUBROUTINE write_histo_csv (fh, filename, hbnds, mov_avg_width, histogram)
   ! Arg_divider acts as a parameter defining a moving average (!)
   ! It has an immediate effect like a filtered graph.
   INTEGER  (KIND=ik), INTENT(IN) :: fh
   CHARACTER(len=*)  , INTENT(IN) :: filename
   INTEGER  (KIND=ik), DIMENSION(3), INTENT(IN) :: hbnds    ! histogram lower/upper bounds
   INTEGER  (KIND=ik)              , INTENT(IN) :: mov_avg_width
   INTEGER  (KIND=ik), DIMENSION(:), ALLOCATABLE, INTENT(IN) :: histogram

   INTEGER  (KIND=ik) :: ii, avg, span, step
   
   span = mov_avg_width/2 ! int division
   IF (mov_avg_width .EQ. 0_ik) step = 1_ik

   OPEN(UNIT = fh, FILE=TRIM(filename), ACTION="WRITE", STATUS="new")
   WRITE(fh,'(A)') "scaledHU, Voxels"
   DO ii=hbnds(1)+span, hbnds(2)-span, step
      avg = SUM(histogram(ii-span:ii+span))/(mov_avg_width+1)
           IF ( histogram(ii) .GT. 0_ik ) THEN 
                   WRITE(fh,'(I18,A,I18)') ii," , ",avg
           END IF
   END DO
   CLOSE(fh)
 END SUBROUTINE write_histo_csv

END MODULE file_routines_mpi
