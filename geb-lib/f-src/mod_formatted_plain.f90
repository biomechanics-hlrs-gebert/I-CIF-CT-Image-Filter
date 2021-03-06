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
USE mechanical
USE strings
USE user_interaction
USE math

IMPLICIT NONE

CONTAINS

!------------------------------------------------------------------------------
! SUBROUTINE: write_matrix
!------------------------------------------------------------------------------  
!> @author Johannes Gebert - HLRS - NUM - gebert@hlrs.de
!> @author Ralf Schneider - HLRS - NUM - schneider@hlrs.de
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
!> @param[in] name Name of the object to print
!> @param[in] mat Actual matrix
!> @param[in] fmti Formatting of the data
!> @param[in] unit Physical unit of the information to print
!------------------------------------------------------------------------------
SUBROUTINE write_matrix(fh, name, mat, fmti, unit)

INTEGER(INT64), INTENT(IN) :: fh   
REAL(rk), INTENT(IN), DIMENSION(:, :) :: mat    
CHARACTER(*), INTENT(IN) :: name 
CHARACTER(*), INTENT(IN), OPTIONAL :: fmti, unit 

! Internal variables 
CHARACTER(mcl) :: fmt_a, sep, nm_fmt, text, fmt_u
INTEGER(ik) :: prec , fw, nm_fmt_lngth, ii, jj, kk, dim1, dim2
REAL(rk) :: sym_out
LOGICAL :: sym

!------------------------------------------------------------------------------
! Initialize and check for presence of the variables
!------------------------------------------------------------------------------
dim1 = SIZE(mat, 1)
dim2 = SIZE(mat, 2)
fmt_u = 'standard'
mssg = '' 
text = ''

prec = PRECISION(mat)
fw = prec+8
sym = .FALSE.


IF (PRESENT(unit)) THEN
    IF (unit /= '') text = " Unit: ("//TRIM(unit)//")"
END IF

!------------------------------------------------------------------------------
! Check symmetry
!------------------------------------------------------------------------------
IF (dim1 == dim2) THEN
    sym = .TRUE.
    CALL check_sym(mat, sym_out)
END IF     
!------------------------------------------------------------------------------
! Generate formats
!------------------------------------------------------------------------------
IF(PRESENT(fmti)) fmt_u = fmti

SELECT CASE (TRIM(fmt_u))
   CASE ('std', 'standard')
        WRITE(fmt_a, "(3(A,I0),A)") "(",dim2,"(E",fw,".",prec,"E2))"
        WRITE(sep  , "(A,I0,A)")    "(",fw*dim2,"('-'))"

        ! Calculate text and unit length. If name to long - overflow formaming to the right
        nm_fmt_lngth  = fw*dim2-4-LEN_TRIM(name)-LEN_TRIM(text)

   CASE ('spl', 'simple')
        WRITE(fmt_a,  "(3(A,I0),A)") "(",dim2,"(F15.4))"
        WRITE(sep  ,  "(A,I0,A)")    "(",dim2*15,"('-'))"        

        ! Calculate text and unit length. If name to long - overflow formaming to the right
        nm_fmt_lngth  = dim2*10-4-LEN_TRIM(name)-LEN_TRIM(text) 

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
WRITE(fh, sep)                                    ! Separator
WRITE(fh, nm_fmt) ' ',TRIM(name), ' ', TRIM(text) ! Named separator

DO ii=1, dim1
    DO jj=1, dim2

        IF ((sym) .AND. (ii==dim1) .AND. (jj==1)) THEN
            SELECT CASE(fmt_u)
            CASE('spl', 'simple')
                WRITE(fh, '(A)', ADVANCE='NO') "symmetric "
            CASE('std', 'standard')
                WRITE(fh, '(A)', ADVANCE='NO') "   symmetric           "
            END SELECT

        ELSE IF ((sym) .AND. (ii==dim1) .AND. (jj==2)) THEN
            !------------------------------------------------------------------------------
            ! Numbers are returned = 0._rk, but formatting is special here.
            !------------------------------------------------------------------------------
            IF (ABS(sym_out) <=  num_zero) THEN
                sym_out = 0._rk
            END IF

            WRITE(fh, fmt_a, ADVANCE='NO') sym_out          
        ELSE
            IF ((ABS(mat(ii,jj)) >=  num_zero) .AND. ((.NOT. sym) .OR. ((sym) .AND. (jj .GE. ii)))) THEN 
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

WRITE(fh, '(A)') '' ! Newline & Carriage return
END SUBROUTINE write_matrix

!------------------------------------------------------------------------------
! SUBROUTINE: write_matrix_int
!------------------------------------------------------------------------------  
!> @author Johannes Gebert - HLRS - NUM - gebert@hlrs.de
!> @author Ralf Schneider - HLRS - NUM - schneider@hlrs.de
!
!> @brief
!> Subroutine to print matrices of integers. Lean functionality. 
!> Choose reals for other functionality.
!
!> @Description
!> Please provide mat_real OR mat_in :-)
!> Automatically writes "sym" if hide_zeros .EQV. .TRUE. and left triangle = 0
!> Hide zeroes is set as default.
!> Accepted formats: 'std'/'standard' for scientific formatting and
!> 'spl'/'simple' for traditional formatting
!
!> @param[in] fh Handle of file to print to
!> @param[in] name Name of the object to print
!> @param[in] mat Actual matrix
!> @param[in] fmti Formatting of the data
!> @param[in] unit Physical unit of the information to print
!------------------------------------------------------------------------------
SUBROUTINE write_matrix_int(fh, name, mat, unit)

INTEGER(ik), INTENT(IN) :: fh, mat(:,:)
CHARACTER(*), INTENT(IN) :: name, unit

CHARACTER(mcl) :: fmt_a, sep, nm_fmt, text, fmt_u
INTEGER(ik) :: prec , fw, nm_fmt_lngth, ii, jj, kk, dim2
REAL(rk) :: sym_out
LOGICAL :: sym

!------------------------------------------------------------------------------
! Initialize and check for presence of the variables
!------------------------------------------------------------------------------
text = ''
dim2 = SIZE(mat, 2)
fw=8_ik
IF (unit /= '') text = " Unit: ("//TRIM(unit)//")"

!------------------------------------------------------------------------------
! Generate formats
!------------------------------------------------------------------------------
WRITE(fmt_a, "(2(A,I0),A)") "(",dim2,"(I", fw, "))s"
WRITE(sep,   "(A,I0,A)")    "(",dim2*fw,"('-'))"

! Calculate text and unit length. If name to long - overflow formaming to the right
nm_fmt_lngth  = fw*dim2-4-LEN_TRIM(name)-LEN_TRIM(text)

IF (nm_fmt_lngth <= 1_ik) nm_fmt_lngth = 1_ik
WRITE(nm_fmt, "(A,I0,A)")  "(2('-') ,3A,", nm_fmt_lngth ,"('-'), A)"    

!------------------------------------------------------------------------------
! Write output
!------------------------------------------------------------------------------
WRITE(fh, sep)                                    ! Separator
WRITE(fh, nm_fmt) ' ',TRIM(name), ' ', TRIM(text) ! Named separator

DO ii=1, SIZE(mat, 1)
    DO jj=1, SIZE(mat, 2)
        WRITE(fh, fmt_a, ADVANCE='NO') mat (ii,jj)
    END DO

    WRITE(fh,'(A)') ''
END DO        

WRITE(fh, '(A)') '' ! Newline & Carriage return
END SUBROUTINE write_matrix_int

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
  CHARACTER(*) :: instring, outstring
  INTEGER(ik) :: ii

  outstring=instring
  DO ii=1, LEN_TRIM(instring)
    IF (instring(ii:ii) == '_')  outstring(ii:ii) = ' '
  END DO

  outstring=ADJUSTL(TRIM(outstring))
END SUBROUTINE underscore_to_blank

!------------------------------------------------------------------------------
! SUBROUTINE: check_tensor_2nd_rank_R66_header
!------------------------------------------------------------------------------  
!> @author Johannes Gebert, gebert@hlrs.de, HLRS/NUM
!
!> @brief
!> Validate the header of a tensor_2nd_rank_R66 file
!
!> @param[in] header Header of the tensor_2nd_rank_R66 file
!> @param[out] abort Feedback whether to stop the program
!------------------------------------------------------------------------------  
SUBROUTINE check_tensor_2nd_rank_R66_header(header, abort)

CHARACTER(*), INTENT(IN) :: header
LOGICAL, INTENT(OUT) :: abort

CHARACTER(mcl) :: tokens(100)
INTEGER(ik) :: ntokens

abort = .FALSE.

CALL parse(TRIM(ADJUSTL(header)), ",", args=tokens, nargs=ntokens)

IF(ntokens /= 45_ik) abort = .TRUE.

!------------------------------------------------------------------------------
! Implementation might look silly. But calculating the indices during runtime, 
! to translate them into string and compare them with an if-function or 
! parsing the strings to numbers for comparison will consume way more time.
!------------------------------------------------------------------------------
IF(TRIM(ADJUSTL(tokens( 1))) /= "Domain")     abort = .TRUE.
IF(TRIM(ADJUSTL(tokens( 2))) /= "Criteria")   abort = .TRUE.
IF(TRIM(ADJUSTL(tokens( 3))) /= "Density")    abort = .TRUE.
IF(TRIM(ADJUSTL(tokens( 4))) /= "DoA_Zener")  abort = .TRUE.
IF(TRIM(ADJUSTL(tokens( 5))) /= "DoA_Gebert") abort = .TRUE.
IF(TRIM(ADJUSTL(tokens( 6))) /= "Sym_dev")    abort = .TRUE.
IF(TRIM(ADJUSTL(tokens( 7))) /= "pos_alpha")  abort = .TRUE.
IF(TRIM(ADJUSTL(tokens( 8))) /= "pos_eta")    abort = .TRUE.
IF(TRIM(ADJUSTL(tokens( 9))) /= "pos_phi")    abort = .TRUE.
IF(TRIM(ADJUSTL(tokens(10))) /= "S11") abort = .TRUE.
IF(TRIM(ADJUSTL(tokens(11))) /= "S21") abort = .TRUE.
IF(TRIM(ADJUSTL(tokens(12))) /= "S31") abort = .TRUE.
IF(TRIM(ADJUSTL(tokens(13))) /= "S41") abort = .TRUE.
IF(TRIM(ADJUSTL(tokens(14))) /= "S51") abort = .TRUE.
IF(TRIM(ADJUSTL(tokens(15))) /= "S61") abort = .TRUE.
IF(TRIM(ADJUSTL(tokens(16))) /= "S12") abort = .TRUE.
IF(TRIM(ADJUSTL(tokens(17))) /= "S22") abort = .TRUE.
IF(TRIM(ADJUSTL(tokens(18))) /= "S32") abort = .TRUE.
IF(TRIM(ADJUSTL(tokens(19))) /= "S42") abort = .TRUE.
IF(TRIM(ADJUSTL(tokens(20))) /= "S52") abort = .TRUE.
IF(TRIM(ADJUSTL(tokens(21))) /= "S62") abort = .TRUE.
IF(TRIM(ADJUSTL(tokens(22))) /= "S13") abort = .TRUE.
IF(TRIM(ADJUSTL(tokens(23))) /= "S23") abort = .TRUE.
IF(TRIM(ADJUSTL(tokens(24))) /= "S33") abort = .TRUE.
IF(TRIM(ADJUSTL(tokens(25))) /= "S43") abort = .TRUE.
IF(TRIM(ADJUSTL(tokens(26))) /= "S53") abort = .TRUE.
IF(TRIM(ADJUSTL(tokens(27))) /= "S63") abort = .TRUE.
IF(TRIM(ADJUSTL(tokens(28))) /= "S14") abort = .TRUE.
IF(TRIM(ADJUSTL(tokens(29))) /= "S24") abort = .TRUE.
IF(TRIM(ADJUSTL(tokens(30))) /= "S34") abort = .TRUE.
IF(TRIM(ADJUSTL(tokens(31))) /= "S44") abort = .TRUE.
IF(TRIM(ADJUSTL(tokens(32))) /= "S54") abort = .TRUE.
IF(TRIM(ADJUSTL(tokens(33))) /= "S64") abort = .TRUE.
IF(TRIM(ADJUSTL(tokens(34))) /= "S15") abort = .TRUE.
IF(TRIM(ADJUSTL(tokens(35))) /= "S25") abort = .TRUE.
IF(TRIM(ADJUSTL(tokens(36))) /= "S35") abort = .TRUE.
IF(TRIM(ADJUSTL(tokens(37))) /= "S45") abort = .TRUE.
IF(TRIM(ADJUSTL(tokens(38))) /= "S55") abort = .TRUE.
IF(TRIM(ADJUSTL(tokens(39))) /= "S65") abort = .TRUE.
IF(TRIM(ADJUSTL(tokens(40))) /= "S16") abort = .TRUE.
IF(TRIM(ADJUSTL(tokens(41))) /= "S26") abort = .TRUE.
IF(TRIM(ADJUSTL(tokens(42))) /= "S36") abort = .TRUE.
IF(TRIM(ADJUSTL(tokens(43))) /= "S46") abort = .TRUE.
IF(TRIM(ADJUSTL(tokens(44))) /= "S56") abort = .TRUE.
IF(TRIM(ADJUSTL(tokens(45))) /= "S66") abort = .TRUE.

END SUBROUTINE check_tensor_2nd_rank_R66_header


!------------------------------------------------------------------------------
! SUBROUTINE: write_tensor_2nd_rank_R66_header
!------------------------------------------------------------------------------  
!> @author Johannes Gebert, gebert@hlrs.de, HLRS/NUM
!
!> @brief
!> Write the header of a tensor_2nd_rank_R66 file.
!
!> @param[in] fh File handle, the header will be written to
!------------------------------------------------------------------------------  
SUBROUTINE write_tensor_2nd_rank_R66_header(fh)

INTEGER(ik), INTENT(IN) :: fh

!------------------------------------------------------------------------------
! Implementation might look silly. But calculating the indices during runtime, 
! to translate them into string and compare them with an if-function or 
! parsing the strings to numbers for comparison will consume way more time.
!------------------------------------------------------------------------------
WRITE(fh, '(A)', ADVANCE='NO')  "Domain, "
WRITE(fh, '(A)', ADVANCE='NO')  "Criteria, "
WRITE(fh, '(A)', ADVANCE='NO')  "Density, "
WRITE(fh, '(A)', ADVANCE='NO')  "DoA_Zener, "
WRITE(fh, '(A)', ADVANCE='NO')  "DoA_Gebert, "
WRITE(fh, '(A)', ADVANCE='NO')  "Sym_dev, "
WRITE(fh, '(A)', ADVANCE='NO')  "pos_alpha, "
WRITE(fh, '(A)', ADVANCE='NO')  "pos_eta, "
WRITE(fh, '(A)', ADVANCE='NO')  "pos_phi, "
WRITE(fh, '(A)', ADVANCE='NO')  "S11, "
WRITE(fh, '(A)', ADVANCE='NO')  "S21, "
WRITE(fh, '(A)', ADVANCE='NO')  "S31, "
WRITE(fh, '(A)', ADVANCE='NO')  "S41, "
WRITE(fh, '(A)', ADVANCE='NO')  "S51, "
WRITE(fh, '(A)', ADVANCE='NO')  "S61, "
WRITE(fh, '(A)', ADVANCE='NO')  "S12, "
WRITE(fh, '(A)', ADVANCE='NO')  "S22, "
WRITE(fh, '(A)', ADVANCE='NO')  "S32, "
WRITE(fh, '(A)', ADVANCE='NO')  "S42, "
WRITE(fh, '(A)', ADVANCE='NO')  "S52, "
WRITE(fh, '(A)', ADVANCE='NO')  "S62, "
WRITE(fh, '(A)', ADVANCE='NO')  "S13, "
WRITE(fh, '(A)', ADVANCE='NO')  "S23, "
WRITE(fh, '(A)', ADVANCE='NO')  "S33, "
WRITE(fh, '(A)', ADVANCE='NO')  "S43, "
WRITE(fh, '(A)', ADVANCE='NO')  "S53, "
WRITE(fh, '(A)', ADVANCE='NO')  "S63, "
WRITE(fh, '(A)', ADVANCE='NO')  "S14, "
WRITE(fh, '(A)', ADVANCE='NO')  "S24, "
WRITE(fh, '(A)', ADVANCE='NO')  "S34, "
WRITE(fh, '(A)', ADVANCE='NO')  "S44, "
WRITE(fh, '(A)', ADVANCE='NO')  "S54, "
WRITE(fh, '(A)', ADVANCE='NO')  "S64, "
WRITE(fh, '(A)', ADVANCE='NO')  "S15, "
WRITE(fh, '(A)', ADVANCE='NO')  "S25, "
WRITE(fh, '(A)', ADVANCE='NO')  "S35, "
WRITE(fh, '(A)', ADVANCE='NO')  "S45, "
WRITE(fh, '(A)', ADVANCE='NO')  "S55, "
WRITE(fh, '(A)', ADVANCE='NO')  "S65, "
WRITE(fh, '(A)', ADVANCE='NO')  "S16, "
WRITE(fh, '(A)', ADVANCE='NO')  "S26, "
WRITE(fh, '(A)', ADVANCE='NO')  "S36, "
WRITE(fh, '(A)', ADVANCE='NO')  "S46, "
WRITE(fh, '(A)', ADVANCE='NO')  "S56, "
WRITE(fh, '(A)'              )  "S66"

END SUBROUTINE write_tensor_2nd_rank_R66_header

!------------------------------------------------------------------------------
! SUBROUTINE: parse_tensor_2nd_rank_R66_row
!------------------------------------------------------------------------------  
!> @author Johannes Gebert, gebert@hlrs.de, HLRS/NUM
!
!> @brief
!> Parse the header of a tensor_2nd_rank_R66 file
!
!> @description
!> Function returns a negative first column (domain numer) if the data 
!> within a row does not meet the only valid amount of entries.
!
!> @param[in] row An unparsed row of the tensor_2nd_rank_R66 file
!> @param[out] tensor_of_row The resulting TYPE tensor_2nd_rank_R66 file
!> @param[out] invalid Whether the data is corrupt
!------------------------------------------------------------------------------  
SUBROUTINE parse_tensor_2nd_rank_R66_row(row, tensor_of_row, invalid)

CHARACTER(*), INTENT(IN) :: row
TYPE(tensor_2nd_rank_R66), INTENT(OUT) :: tensor_of_row
LOGICAL, INTENT(OUT) :: invalid

CHARACTER(mcl) :: tokens(100)
INTEGER(ik) :: ntokens, ii, jj, tkn

invalid = .FALSE.

CALL parse(TRIM(ADJUSTL(row)), ",", args=tokens, nargs=ntokens)

IF(ntokens /= 45_ik) invalid = .TRUE.

!------------------------------------------------------------------------------
! Required, even if row is invalid. Does not check whether tokens(1) actually
! is an integer/domain number.
!------------------------------------------------------------------------------
READ(tokens(1), '(I20)') tensor_of_row%dmn 

IF(invalid) THEN   
    WRITE(std_out, FMT_WRN_xAI0) &
        "Invalid 'domain' "//TRIM(tokens(1))//" amount of tokens: ", ntokens

    tensor_of_row%dmn        = -tensor_of_row%dmn
    tensor_of_row%crit       = ""
    tensor_of_row%density    = 0._rk
    tensor_of_row%doa_zener  = 0._rk
    tensor_of_row%doa_gebert = 0._rk
    tensor_of_row%pos        = 0._rk
    tensor_of_row%mat        = 0._rk    
ELSE   
    READ(tokens(2), '(A)')      tensor_of_row%crit
    READ(tokens(3), '(F39.10)') tensor_of_row%density
    READ(tokens(4), '(F39.10)') tensor_of_row%doa_zener
    READ(tokens(5), '(F39.10)') tensor_of_row%doa_gebert
    READ(tokens(6), '(F39.10)') tensor_of_row%sym

    DO ii=1, 3
        READ(tokens(6_ik + ii), '(F39.10)') tensor_of_row%pos(ii)
    END DO

    tkn = 10_ik
    DO jj=1, 6
        DO ii=1, 6
            READ(tokens(tkn), '(F39.10)') tensor_of_row%mat(ii, jj)
            tkn = tkn + 1_ik
        END DO
    END DO
END IF

END SUBROUTINE parse_tensor_2nd_rank_R66_row

!------------------------------------------------------------------------------
! SUBROUTINE: write_tensor_2nd_rank_R66_row
!------------------------------------------------------------------------------  
!> @author Johannes Gebert, gebert@hlrs.de, HLRS/NUM
!
!> @brief
!> Write a tensor_2nd_rank_R66 into a row of a tensor_2nd_rank_R66 file.
!
!> @param[in] fh File handle, the header will be written to
!> @param[in] tensor_of_row The resulting TYPE tensor_2nd_rank_R66 file
!------------------------------------------------------------------------------  
SUBROUTINE write_tensor_2nd_rank_R66_row(fh, tensor_of_row)

INTEGER(ik), INTENT(IN) :: fh
TYPE(tensor_2nd_rank_R66), INTENT(IN) :: tensor_of_row

INTEGER(ik) :: ii, jj

!------------------------------------------------------------------------------  
! The formats are independent of other plain text formats to extract a higher
! number of digits, which maintains accuracy better.
!------------------------------------------------------------------------------  
CHARACTER(Len=*), PARAMETER :: FINT = "(I0, A)"
CHARACTER(Len=*), PARAMETER :: FREAL = "(F0.10, A)"
CHARACTER(Len=*), PARAMETER :: FSCI  = "(E17.10, A)"
CHARACTER(Len=*), PARAMETER :: FCHAR = "(2A)"

WRITE(fh, FINT , ADVANCE='NO') tensor_of_row%dmn, ", "
WRITE(fh, FCHAR, ADVANCE='NO') TRIM(tensor_of_row%crit), ", "
WRITE(fh, FREAL, ADVANCE='NO') tensor_of_row%density, ", "
WRITE(fh, FREAL, ADVANCE='NO') tensor_of_row%doa_zener, ", "
WRITE(fh, FREAL, ADVANCE='NO') tensor_of_row%doa_gebert, ", "
WRITE(fh, FSCI , ADVANCE='NO') tensor_of_row%sym, ", "

DO ii=1, 3
    WRITE(fh, FREAL, ADVANCE='NO') tensor_of_row%pos(ii), ", "
END DO

DO jj=1, 6
    DO ii=1, 6
        WRITE(fh, FREAL, ADVANCE='NO') tensor_of_row%mat(ii,jj), ", "
    END DO
END DO

!------------------------------------------------------------------------------
! Newline, since mat(6,6) gets printed with ADVANCE='NO'
!------------------------------------------------------------------------------
WRITE(fh, '(A)') 

END SUBROUTINE write_tensor_2nd_rank_R66_row

!------------------------------------------------------------------------------
! SUBROUTINE: parse_tensor_2nd_rank_R66
!------------------------------------------------------------------------------  
!> @author Johannes Gebert, gebert@hlrs.de, HLRS/NUM
!
!> @brief
!> Parse a "csv" file with the structure of a tensor_2nd_rank_R66 data type.
!
!> @param[in] fh Input file handle (file existane was checked before).
!> @param[in] filename Filename of the input file
!> @param[in] amnt_lines Number of lines contained in the file
!> @param[out] tensors_in Array of all input tensors
!> @param[out] invalid_entries Amount of invalid entries
!------------------------------------------------------------------------------  
SUBROUTINE parse_tensor_2nd_rank_R66(fh, filename, amnt_lines, tensors_in, invalid_entries)

INTEGER(ik), INTENT(IN) :: fh, amnt_lines
CHARACTER(*), INTENT(IN) :: filename
TYPE(tensor_2nd_rank_R66), DIMENSION(:), INTENT(OUT) :: tensors_in
INTEGER(ik), INTENT(OUT) :: invalid_entries
!------------------------------------------------------------------------------
! Line of csv may be pretty long (>40 floats)
!------------------------------------------------------------------------------
CHARACTER(10_ik*mcl) :: header, line
INTEGER(ik) :: ii
LOGICAL :: abort, invalid

invalid_entries = 0_ik
invalid = .FALSE. 

!------------------------------------------------------------------------------
! Parse the header of the "csv" data
! Error if the header does not indicate 43 entries (determined by TYPE)
!------------------------------------------------------------------------------
READ(fh, '(A)') header
CALL check_tensor_2nd_rank_R66_header(header, abort)

IF(abort) THEN
    mssg = "The input file "//TRIM(ADJUSTL(filename))//&
      &" is not a valid 'tensor_2nd_rank_R66' file."
    CALL print_err_stop(fh, mssg, 1)
END IF

DO ii = 2_ik, amnt_lines
    READ(fh, '(A)') line

    CALL parse_tensor_2nd_rank_R66_row(line, tensors_in(ii-1_ik), invalid)
    IF(invalid) invalid_entries = invalid_entries + 1_ik 
END DO

END SUBROUTINE parse_tensor_2nd_rank_R66

!------------------------------------------------------------------------------
! SUBROUTINE: write_tensor_2nd_rank_R66
!------------------------------------------------------------------------------  
!> @author Johannes Gebert, gebert@hlrs.de, HLRS/NUM
!
!> @brief
!> Write a "csv" file with the structure of a tensor_2nd_rank_R66 data type.
!
!> @param[in] fh Input file handle (file existane was checked before).
!> @param[in] amnt_lines Number of lines contained in the file
!> @param[out] tensors_in Array of all input tensors
!------------------------------------------------------------------------------  
SUBROUTINE write_tensor_2nd_rank_R66(fh, amnt_lines, tensors_in)

INTEGER(ik), INTENT(IN) :: fh, amnt_lines
TYPE(tensor_2nd_rank_R66), DIMENSION(:), INTENT(IN) :: tensors_in

INTEGER(ik) :: ii

CALL write_tensor_2nd_rank_R66_header(fh)

!------------------------------------------------------------------------------  
! Written this way (2...amnt_lines and ii-1 as indice) to remember how to  
! deal with the header. Must be kept in mind.
!------------------------------------------------------------------------------  
DO ii = 2_ik, amnt_lines
    CALL write_tensor_2nd_rank_R66_row(fh, tensors_in(ii-1_ik))
END DO

END SUBROUTINE write_tensor_2nd_rank_R66


!------------------------------------------------------------------------------
! SUBROUTINE: std_stop
!------------------------------------------------------------------------------  
!> @author Johannes Gebert - HLRS - NUM - gebert@hlrs.de
!
!> @brief
!> Stop program if a subroutine signals an error.
!
!> @param[in] stat Status integer
!> @param[out] abrt Whether to abort the program.
!------------------------------------------------------------------------------
SUBROUTINE std_stop(stat, abrt)

CHARACTER(*), INTENT(IN) :: stat
LOGICAL, INTENT(INOUT) :: abrt

INTEGER(mik) :: ierr, my_rank


CALL MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierr)

!------------------------------------------------------------------------------
! Set my_rank = 0 to enable error printing subseqently.
!------------------------------------------------------------------------------
IF(ierr /= 0_ik) my_rank = 0_ik

IF((stat /= "") .AND. (my_rank==0)) THEN
    WRITE(std_out, FMT_ERR) "Error in keyword '"//TRIM(stat)//"'."
    abrt = .TRUE.
END IF 

END SUBROUTINE std_stop
    

END MODULE formatted_plain
