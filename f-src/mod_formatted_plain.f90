!------------------------------------------------------------------------------
! MODULE: formatted_plain
!------------------------------------------------------------------------------
!> @author Johannes Gebert - HLRS - NUM - gebert@hlrs.de
!
! @description: 
!> "Formatted I/O for plain ascii files that are not related to PureDat, Meta, 
!> and log/mon messages. 
!------------------------------------------------------------------------------
MODULE formatted_plain

USE ISO_FORTRAN_ENV
USE global_std
USE math

IMPLICIT NONE

INTERFACE write_matrix
  MODULE PROCEDURE write_matrix_int
  MODULE PROCEDURE write_matrix_real 
END INTERFACE write_matrix

CONTAINS

!------------------------------------------------------------------------------
! SUBROUTINE: write_matrix_real
!------------------------------------------------------------------------------  
!> @author Ralf Schneider, schneider@hlrs.de, HLRS/NUM
!> @author Johannes Gebert,   gebert@hlrs.de, HLRS/NUM
!
!> @brief
!> Subroutine to print regular tensors respectively matrices.
!
!> @Description
!> Please provide mat_real OR mat_in :-)
!> Automatically writes "sym" if hide_zeros .EQV. .TRUE. and left triangle = 0
!> Hide zeroes is set as default.
!> Accepted formats: 'std'/'standard' for scientific formatting and
!> 'spl'/'simple' for traditional formatting
!
!> @param[in] fh Handle of file to print to
!> @param[in] dim1 Object to print
!> @param[in] dim2 Object to print
!> @param[in] name Name of the object to print
!> @param[in] mat_real Dimensions of the 2nd rank tensor, double precision
!> @param[in] mat_int  Dimensions of the 2nd rank tensor, integer kind = 4
!> @param[in] fmt Formatting of the data
!> @param[in] unit Physical unit of the information to print
!> @param[in] hide_zeros Whether to suppress zeros for printing matrices
!------------------------------------------------------------------------------
SUBROUTINE write_matrix_real (fh, name, fmt, unit, mat)

INTEGER(KIND=ik), INTENT(IN) :: fh   
CHARACTER(LEN=*), INTENT(IN) :: name 
REAL   (KIND=rk), DIMENSION(:, :), INTENT(IN) :: mat    
CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: fmt 
CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: unit 

! Internal variables 
INTEGER(KIND=ik)   :: prec , fw, nm_fmt_lngth, ii, jj, kk, dim1, dim2
CHARACTER(LEN=mcl) :: fmt_a, sep, nm_fmt
CHARACTER(LEN=mcl) :: text
CHARACTER(LEN=mcl) :: fmt_u
REAL(KIND=rk) :: sym_out
LOGICAL :: sym_u

!------------------------------------------------------------------------------
! Initialize and check for presence of the variables
!------------------------------------------------------------------------------
dim1 = SIZE(mat, 1)
dim2 = SIZE(mat, 2)
fmt_u = 'standard'
mssg='' 
text = ''

prec = PRECISION(mat)
fw = prec+8
sym_u = .FALSE.

IF (PRESENT(unit)) THEN
    IF (unit /= '') text = " Unit: ("//TRIM(unit)//")"
END IF

IF (dim1 == dim2) THEN
    CALL check_sym(fh, mat, name, sym_out=sym_out)
    sym_u = .TRUE.
END IF

!------------------------------------------------------------------------------
! Generate formats
!------------------------------------------------------------------------------
IF(PRESENT(fmt)) fmt_u = fmt

SELECT CASE (TRIM(fmt_u))
   CASE ('std', 'standard')
        WRITE(fmt_a, "(3(A,I0),A)") "(",dim2,"(E",fw,".",prec,"E2))"
        WRITE(sep  , "(A,I0,A)")    "(",fw*dim2,"('-'))"

        ! Calculate text and unit length. If name to long - overflow formaming to the right
        nm_fmt_lngth  = fw*dim2-4-LEN_TRIM(name)-LEN_TRIM(text)

   CASE ('spl', 'simple')
        WRITE(fmt_a,  "(3(A,I0),A)") "(",dim2,"(F10.3))"
        WRITE(sep  ,  "(A,I0,A)")    "(",dim2*10,"('-'))"        

        ! Calculate text and unit length. If name to long - overflow formaming to the right
        nm_fmt_lngth  = dim2*10-4-2-LEN_TRIM(name)-LEN_TRIM(text) 

   CASE('wxm', 'wxmaxima')
       WRITE(fmt_a, '(5(A,I0),A)')  "(' [',",dim2-1,"(E",fw,".",prec,"E2,','),E",fw,".",prec,"E2,'],' )"

       WRITE(fh,"(A,A)")TRIM(name),": matrix("

       DO kk = 1, dim1 - 1
          WRITE(fh, fmt_a) mat(kk,:)
       END DO

       WRITE(fmt_a,'(5(A,I0),A)')  "(' [',",dim2-1,"(E",fw,".",prec,"E2,','),E",fw,".",prec,"E2,']);' )"

       WRITE(fh, fmt_a) mat(dim1, :)
END SELECT

IF (nm_fmt_lngth .LT. 1_ik) nm_fmt_lngth = 1_ik
WRITE(nm_fmt, "(A,I0,A)")  "(2('-') ,3A,", nm_fmt_lngth ,"('-'), A)"    



!------------------------------------------------------------------------------
! Write output
!------------------------------------------------------------------------------
WRITE(fh, '(A)')
WRITE(fh, sep)                                    ! Separator
WRITE(fh, nm_fmt) ' ',TRIM(name), ' ', TRIM(text) ! Named separator

DO ii=1, dim1
DO jj=1, dim2

    IF ((sym_u) .AND. (ii==dim1) .AND. (jj==1)) THEN
        SELECT CASE(fmt_u)
        CASE('spl', 'simple')
            WRITE(fh, '(A)', ADVANCE='NO') "symmetric "
        CASE('std', 'standard')
            WRITE(fh, '(A)', ADVANCE='NO') "   symmetric           "
        END SELECT

    ELSE IF ((sym_u) .AND. (ii==dim1) .AND. (jj==2)) THEN
        IF (ABS(sym_out) <=  10E-08) sym_out = 0._rk
        WRITE(fh, fmt_a, ADVANCE='NO') sym_out          
    ELSE
        IF ((ABS(mat(ii,jj)) >=  10E-08) .AND. ((.NOT. sym_u) .OR. ((sym_u) .AND. (jj .GE. ii)))) THEN 
            WRITE(fh, fmt_a, ADVANCE='NO') mat (ii,jj)
        ELSE
            SELECT CASE(fmt_u)
            CASE('spl', 'simple')
                WRITE(fh, '(A)', ADVANCE='NO') "      .   "
            CASE('std', 'standard')
                WRITE(fh, '(A)', ADVANCE='NO') "   .                   "
            END SELECT
        END IF
    END IF

END DO
WRITE(fh,'(A)') ''
END DO        

WRITE(fh, '(A)') ''                               ! Newline & Carriage return
End Subroutine write_matrix_real


!------------------------------------------------------------------------------
! SUBROUTINE: write_matrix_int
!------------------------------------------------------------------------------  
!> @author Ralf Schneider, schneider@hlrs.de, HLRS/NUM
!> @author Johannes Gebert,   gebert@hlrs.de, HLRS/NUM
!
!> @brief
!> Subroutine to print regular tensors respectively matrices.
!
!> @Description
!> Please provide mat_real OR mat_in :-)
!> Accepted formats: 'std'/'standard' for scientific formatting and
!> 'spl'/'simple' for traditional formatting
!
!> @param[in] fh Handle of file to print to
!> @param[in] name Name of the object to print
!> @param[in] mat_real Dimensions of the 2nd rank tensor, double precision
!> @param[in] mat_int  Dimensions of the 2nd rank tensor, integer kind = 4
!> @param[in] fmt Formatting of the data
!> @param[in] unit Physical unit of the information to print
!> @param[in] hide_zeros Whether to suppress zeros for printing matrices
!------------------------------------------------------------------------------
SUBROUTINE write_matrix_int (fh, name, fmt, unit, mat)

INTEGER(KIND=ik), INTENT(IN) :: fh   
CHARACTER(LEN=*), INTENT(IN) :: name 
INTEGER(KIND=ik), DIMENSION(:, :), INTENT(IN) :: mat    
CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: fmt 
CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: unit 

! Internal variables 
INTEGER(KIND=ik)   :: nm_fmt_lngth, ii, jj, dim1, dim2
CHARACTER(LEN=mcl) :: fmt_a, sep, nm_fmt
CHARACTER(LEN=mcl) :: text
CHARACTER(LEN=mcl) :: fmt_u
REAL(KIND=rk) :: sym_out
LOGICAL :: sym_u

!------------------------------------------------------------------------------
! Initialize and check for presence of the variables
!------------------------------------------------------------------------------
dim1 = SIZE(mat, 1)
dim2 = SIZE(mat, 2)
fmt_u = 'standard'
sym_u = .FALSE.
mssg='' 
text = ''

IF (PRESENT(unit)) THEN
    IF (unit /= '') text = " Unit: ("//TRIM(unit)//")"
END IF

IF (dim1 == dim2) THEN
    CALL check_sym(fh, REAL(mat, KIND=rk), name, sym_out=sym_out)
    sym_u = .TRUE.
END IF

!------------------------------------------------------------------------------
! Generate formats
!------------------------------------------------------------------------------
IF(PRESENT(fmt)) fmt_u = fmt

SELECT CASE (TRIM(fmt_u))
   CASE ('std', 'standard')

        WRITE(sep, "(A,I0,A)") "(",dim2,"('-'))"

        ! Calculate text and unit length. If name to long - overflow formatting to the right
        nm_fmt_lngth = dim2-4-2-LEN_TRIM(name)-LEN_TRIM(text)

   CASE ('spl', 'simple')

        WRITE(sep, "(A,I0,A)") "(",dim2*10,"('-'))"        

        ! Calculate text and unit length. If name to long - overflow formaming to the right
        nm_fmt_lngth  = dim2*10-4-LEN_TRIM(name)-LEN_TRIM(text) 
END SELECT

IF (nm_fmt_lngth .LT. 1_ik) nm_fmt_lngth = 1_ik
WRITE(nm_fmt, "(A,I0,A)")  "(2('-') ,3A,", nm_fmt_lngth ,"('-'), A)"    

WRITE(fmt_a, "(3(A,I0),A)") "(",dim2,"(I10))"

!------------------------------------------------------------------------------
! Write output
!------------------------------------------------------------------------------
WRITE(fh, sep)                                    ! Separator
WRITE(fh, nm_fmt) ' ',TRIM(name), ' ', TRIM(text) ! Named separator

DO ii=1, dim1
DO jj=1, dim2

    IF ((sym_u) .AND. (ii==dim1) .AND. (jj==1)) THEN
        WRITE(fh, '(A)', ADVANCE='NO') " symmetric"
    ELSE IF ((sym_u) .AND. (ii==dim1) .AND. (jj==2)) THEN
        IF (ABS(sym_out) <=  10E-08) sym_out = 0._rk
        WRITE(fh, fmt_a, ADVANCE='NO') sym_out          
    ELSE
        IF ((mat(ii,jj) /= 0) .AND. ((.NOT. sym_u) .OR. ((sym_u) .AND. (jj .GE. ii)))) THEN 
            WRITE(fh, fmt_a, ADVANCE='NO') mat (ii,jj)
        ELSE
            WRITE(fh, '(A)', ADVANCE='NO') "         ."
        END IF
    END IF

END DO
WRITE(fh,'(A)') ''
END DO        

WRITE(fh, '(A)') ''
End Subroutine write_matrix_int

!------------------------------------------------------------------------------
! SUBROUTINE: extract_histogram_scalar_array
!------------------------------------------------------------------------------  
!> @author Johannes Gebert - HLRS - NUM - gebert@hlrs.de
!
!> @brief
!> Extracts a histogram from a 3-dimensional array
!> !!! Scaling may need extrensive Rework !!! 
!
!> @param[in] array Actual image
!> @param[in] hbnds Boundary values
!> @param[out] histogram Returns the histogram
!------------------------------------------------------------------------------  
SUBROUTINE extract_histogram_scalar_array (array, hbnds, histogram)
! This is an inherently unflexible subroutine. It delivers exactly this kind of histogram and nothing else...
INTEGER(KIND=ik), DIMENSION(:,:,:) :: array
INTEGER(KIND=ik), DIMENSION(3)             , INTENT(IN)  :: hbnds    ! histogram lower/upper bounds
INTEGER(KIND=ik), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: histogram

! Internal variables
INTEGER(KIND=ik) :: ii, jj, kk
INTEGER(KIND=ik), DIMENSION(3) :: shp

ALLOCATE(histogram(hbnds(1):hbnds(2)))

histogram(:) = 0_ik

shp = SHAPE(array)

! Take care of sign of hmin!! Not that intuitive
DO kk=1, shp(3)
DO jj=1, shp(2)
DO ii=1, shp(1)
  histogram( array(ii, jj, kk) ) = histogram(  array(ii, jj, kk) ) + 1_ik
END DO
END DO
END DO

END SUBROUTINE extract_histogram_scalar_array

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
!> @param[in] hdr_str String of the histograms header
!> @param[in] hbnds Histogram boundaries
!> @param[in] mov_avg_width Width of the moving average
!> @param[in] histogram Actual histogram data
!------------------------------------------------------------------------------  
 SUBROUTINE write_histo_csv (fh, filename, hdr_str, hbnds, mov_avg_width, histogram)
   ! Arg_divider acts as a parameter defining a moving average (!)
   ! It has an immediate effect like a filtered graph.
   INTEGER(KIND=ik), INTENT(IN) :: fh
   CHARACTER(len=*), INTENT(IN) :: filename, hdr_str
   INTEGER(KIND=ik), DIMENSION(3), INTENT(IN) :: hbnds    ! histogram lower/upper bounds
   INTEGER(KIND=ik)              , INTENT(IN) :: mov_avg_width
   INTEGER(KIND=ik), DIMENSION(:), INTENT(IN) :: histogram

   INTEGER  (KIND=ik) :: ii, avg, span, step
   
   span = mov_avg_width / 2 ! int division

   IF (mov_avg_width == 0_ik) step = 1_ik
  
   OPEN(UNIT=fh, FILE=TRIM(filename), ACTION="WRITE", STATUS="NEW")
   
   WRITE(fh,'(A)') TRIM(ADJUSTL(hdr_str)) ! "scaledHU, Voxels"
   
   DO ii = hbnds(1)+span, hbnds(2)-span, step
      avg = SUM(histogram(ii-span:ii+span))/(mov_avg_width+1)

      IF (histogram(ii) >= 0_ik ) WRITE(fh,'(I18,A,I18)') ii," , ",avg
   END DO

   CLOSE(fh)

 END SUBROUTINE write_histo_csv


!------------------------------------------------------------------------------
! SUBROUTINE: underscore_to_blank
!------------------------------------------------------------------------------  
!> @author Johannes Gebert - HLRS - NUM - gebert@hlrs.de
!
!> @brief
!> Replaces underscores with blanks.
!
!> @param[in] instring Input string
!> @param[out] outstring Output string
!------------------------------------------------------------------------------  
SUBROUTINE underscore_to_blank (instring, outstring)
  ! This whole subroutine is a workaround :-)
  CHARACTER(LEN=*) :: instring
  CHARACTER(LEN=*) :: outstring
  INTEGER(KIND=ik) :: ii

  outstring=instring
  DO ii=1, LEN_TRIM(instring)
    IF (instring(ii:ii) == '_')  outstring(ii:ii) = ' '
  END DO

  outstring=ADJUSTL(TRIM(outstring))
END SUBROUTINE underscore_to_blank

END MODULE formatted_plain
