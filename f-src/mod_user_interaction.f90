!------------------------------------------------------------------------------
! MODULE: mod_user_interaction
!------------------------------------------------------------------------------
!> @author Johannes Gebert - HLRS - NUM - gebert@hlrs.de
!
! @Description:
!> Module containing the formatting of all (error) messages
!------------------------------------------------------------------------------
MODULE user_interaction

    USE ISO_FORTRAN_ENV
    USE MPI
    USE global_std
    USE strings

IMPLICIT NONE

!------------------------------------------------------------------------------
! Parameters of (error) messages
!------------------------------------------------------------------------------
INTEGER, PARAMETER :: mw = 90  ! Message width, including leading and trailing descriptors

CHARACTER(Len=*), PARAMETER :: FMT_ERR      = "('EE ', A)"
CHARACTER(Len=*), PARAMETER :: FMT_ERR_STOP = "('EE PROGRAM STOPPED.')"
CHARACTER(Len=*), PARAMETER :: FMT_ERR_SEP  = "('EE ', 76('='))"

CHARACTER(Len=*), PARAMETER :: FMT_ERR_AI0  = "('EE ', *(A,I0))"  

!------------------------------------------------------------------------------
! Text formats
!------------------------------------------------------------------------------
CHARACTER(Len=*), PARAMETER :: FMT_TXT      = "('-- ',10A)"
CHARACTER(Len=*), PARAMETER :: FMT_TXT_SEP  = "(80('-'))"

CHARACTER(Len=*), PARAMETER :: FMT_TXT_AF0A  = "('-- ',A,1X,F0.6,1x,A)"
CHARACTER(Len=*), PARAMETER :: FMT_TXT_AF15A = "('-- ',A,1X,F15.6,1x,A)"
CHARACTER(Len=*), PARAMETER :: FMT_TXT_A3I0  = "('-- ',A,3(1x,I0))"

!------------------------------------------------------------------------------
! Message/debug formats
!------------------------------------------------------------------------------
CHARACTER(Len=*), PARAMETER :: FMT_MSG      = "('MM ',10A)"
CHARACTER(Len=*), PARAMETER :: FMT_MSG_SEP  = "(80('-'))"
!
CHARACTER(Len=*), PARAMETER :: FMT_MSG_AxI0 = "('MM ',A,*(1x, I0))"
CHARACTER(Len=*), PARAMETER :: FMT_MSG_AI0  = "('MM ',*(A,1X,I0,1X))"
CHARACTER(Len=*), PARAMETER :: FMT_MSG_AI0A = "('MM ',A,1X,I0,1X,A)"
CHARACTER(Len=*), PARAMETER :: FMT_MSG_2AI0 = "('MM ',2(A,1X,I0,1X))"
CHARACTER(Len=*), PARAMETER :: FMT_MSG_A3I0 = "('MM ',A,3(1x,I0))"
!
CHARACTER(Len=*), PARAMETER :: FMT_MSG_AF0  = "('MM ',A,1X,F0.6)"
CHARACTER(Len=*), PARAMETER :: FMT_MSG_AF0A = "('MM ',A,1X,F0.6,1x,A)"
CHARACTER(Len=*), PARAMETER :: FMT_MSG_A2F0 = "('MM ',A,2(1X,F0.6))"
CHARACTER(Len=*), PARAMETER :: FMT_MSG_A3F0 = "('MM ',A,3(1X,F0.6))"
!
CHARACTER(Len=*), PARAMETER :: FMT_MSG_AL   = "('MM ',A,1x,L1)"

!------------------------------------------------------------------------------
! Warning formats
!------------------------------------------------------------------------------
CHARACTER(Len=*), PARAMETER :: FMT_WRN      = "('WW ',A)"
CHARACTER(Len=*), PARAMETER :: FMT_WRN_SEP  = "(80('-'))"
!
CHARACTER(Len=*), PARAMETER :: FMT_WRN_AI0  = "('WW ',A,1X,I0)"
CHARACTER(Len=*), PARAMETER :: FMT_WRN_AI0A = "('WW ',A,1X,I0,1X,A)"
CHARACTER(Len=*), PARAMETER :: FMT_WRN_AF0  = "('WW ',A,1X,F0.6)"

!------------------------------------------------------------------------------
! Debug formats
!------------------------------------------------------------------------------
CHARACTER(Len=*), PARAMETER :: FMT_DBG_SEP = "('#DBG#',75('='))"

!------------------------------------------------------------------------------
! Provide colors on std_out (!) 
! Needs to compile with -fbackslash 
! Use of a requires resetting it.
!------------------------------------------------------------------------------
CHARACTER(LEN=*), PARAMETER ::  FMT_Blck    = "\x1B[30m"
CHARACTER(LEN=*), PARAMETER ::  FMT_Red     = "\x1B[31m"
CHARACTER(LEN=*), PARAMETER ::  FMT_Green   = "\x1B[32m"
CHARACTER(LEN=*), PARAMETER ::  FMT_Orange  = "\x1B[33m"
CHARACTER(LEN=*), PARAMETER ::  FMT_Blue    = "\x1B[34m"
CHARACTER(LEN=*), PARAMETER ::  FMT_Purple  = "\x1B[35m"
CHARACTER(LEN=*), PARAMETER ::  FMT_Cyan    = "\x1B[36m"
CHARACTER(LEN=*), PARAMETER ::  FMT_Gray    = "\x1B[37m"
CHARACTER(LEN=*), PARAMETER ::  FMT_nocolor = "\x1B[0m"

CONTAINS

!------------------------------------------------------------------------------
! SUBROUTINE: get_cmd_args
!------------------------------------------------------------------------------  
!> @author Johannes Gebert - HLRS - NUM - gebert@hlrs.de
!
!> @brief
!> Parse cmd args of the process chain. Standard approach.
!
!> @param[out] binary Name of the program
!> @param[out] infile Name of the input file (not meta; XTOM...)
!> @param[out] stp Stop the program?
!> @param[out] restart Check whether a restart is required
!> @param[out] cmd_arg_history Ravision of the program
!------------------------------------------------------------------------------
SUBROUTINE get_cmd_args(binary, infile, stp, restart, cmd_arg_history)

CHARACTER(LEN=*), INTENT(OUT) :: binary, infile
CHARACTER(LEN=*), INTENT(OUT), OPTIONAL :: restart, cmd_arg_history
LOGICAL :: stp

CHARACTER(LEN=mcl) :: cmd_arg
INTEGER(KIND=ik) :: ii

stp = .FALSE.

CALL GET_COMMAND_ARGUMENT(0, binary)

IF (command_argument_count() == 0) THEN 

    CALL usage(binary)

    mssg = "No command argument given."
    stp = .TRUE. 
ELSE
    DO ii=1, 15 ! Read up to 15 command arguments.
    
        CALL GET_COMMAND_ARGUMENT(ii, cmd_arg)

        IF (cmd_arg == '') EXIT

        infile = TRIM(cmd_arg)
        
        cmd_arg_history = TRIM(cmd_arg_history)//' '//TRIM(cmd_arg)

        IF (cmd_arg(1:1) == '-') THEN
            SELECT CASE( cmd_arg(2:LEN_TRIM(cmd_arg)) )
            CASE('-restart', '-Restart')
                restart = 'Y'
            CASE('v', '-Version', '-version')
                CALL show_title()
                stp = .TRUE. 
            CASE('h', '-Help', '-help')
                CALL usage(binary)
                stp = .TRUE. 
            END SELECT
            !
            SELECT CASE( cmd_arg(3:4) )
            CASE('NO', 'no', 'No', 'nO');  restart = 'N'
            END SELECT
        END IF
    END DO

    IF(TRIM(infile) == '') THEN
        mssg = "No input file given via command argument."
        stp = .TRUE. 
    END IF 
END IF

END SUBROUTINE get_cmd_args

!------------------------------------------------------------------------------
! SUBROUTINE: show_title
!------------------------------------------------------------------------------  
!> @author Johannes Gebert - HLRS - NUM - gebert@hlrs.de
!
!> @brief
!> Show brief information about the program. Variables from global_std module!
!------------------------------------------------------------------------------
SUBROUTINE show_title()
WRITE(std_out, FMT_TXT_SEP)
WRITE(std_out, FMT_TXT) 'High-Performance Computing Center | Stuttgart (HLRS)'
WRITE(std_out, FMT_TXT) ''
WRITE(std_out, FMT_TXT) TRIM(ADJUSTL(longname))//' '//TRIM(ADJUSTL(revision))
WRITE(std_out, FMT_TXT) ''     
WRITE(std_out, FMT_TXT) 'Developer & maintainer: Johannes Gebert, M.Sc. (HLRS, NUM)'
WRITE(std_out, FMT_TXT_SEP)
END SUBROUTINE show_title

!------------------------------------------------------------------------------
! FUNCITON: usage
!------------------------------------------------------------------------------  
!> @author Johannes Gebert - HLRS - NUM - gebert@hlrs.de
!
!> @brief
!> Print program usage. 
!
!> @param[in] this_binary Name of the binary - for example parsed from cmd args
!------------------------------------------------------------------------------  
SUBROUTINE usage(this_binary)

CHARACTER(LEN=*), INTENT(IN) :: this_binary

WRITE(std_out, FMT_TXT) 'Usage:'
WRITE(std_out, FMT_TXT) './'//TRIM(ADJUSTL(this_binary))//' '//'<flags> <basename.meta>'
WRITE(std_out, FMT_TXT) ''
WRITE(std_out, FMT_TXT) '-h/ --help      This message.'
WRITE(std_out, FMT_TXT) '-v/ --version   Version of the program'
WRITE(std_out, FMT_TXT) '--restart       Overwrite restart keyword'
WRITE(std_out, FMT_TXT) '--no-restart    Overwrite restart keyword'
WRITE(std_out, FMT_TXT) ''
WRITE(std_out, FMT_TXT) 'The meta-file must be the last command argument.'
WRITE(std_out, FMT_TXT_SEP)

END SUBROUTINE usage


!------------------------------------------------------------------------------
! FUNCITON: determine_stout
!------------------------------------------------------------------------------  
!> @author Johannes Gebert - HLRS - NUM - gebert@hlrs.de
!
!> @brief
!> Check whether the program can access a std_out. Usually given by the
!> environment settings. Fall back always results in writing to a file.
!
!> @return fh_std_out File handle of the »real« std_out
!------------------------------------------------------------------------------  
FUNCTION determine_stout() RESULT(fh_std_out)

INTEGER(KIND=ik) :: fh_std_out

INTEGER(KIND=ik) :: stat
CHARACTER(LEN=scl) :: use_std_out

CALL GET_ENVIRONMENT_VARIABLE(NAME='USE_STD_OUT', VALUE=use_std_out, STATUS=stat)

IF ((stat == 0) .AND. (use_std_out == 'YES')) THEN
    fh_std_out = 6 ! Standard std_out
ELSE
    fh_std_out = give_new_unit()
END IF
 
END FUNCTION determine_stout



!------------------------------------------------------------------------------
! SUBROUTINE: give_new_unit
!------------------------------------------------------------------------------  
!> @author Ralf Schneider - HLRS - NUM - schneider@hlrs.de
!> @author Johannes Gebert - HLRS - NUM - gebert@hlrs.de
!
!> @brief
!> Function which returns a new free unit
!
!> @return new_unit A new file unit
!------------------------------------------------------------------------------
function give_new_unit() result(new_unit)

Integer(KIND=ik) :: new_unit
Integer(KIND=ik) :: ii
Logical :: unit_is_open

Do ii = 300, huge(new_unit)-1

    inquire(unit=ii, opened=unit_is_open)

    if( .not. unit_is_open ) then
        new_unit = ii
        Exit
    end if

End Do

if ( unit_is_open ) then
    mssg = 'Something bad and unexpected happened during search for free unit: &
    &Could not find a new unit between 100 and huge(Int(kind=4))'
    CALL print_err_stop(std_out, mssg, 1_ik)
END IF

End function give_new_unit

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
! SUBROUTINE: print_err_stop
!------------------------------------------------------------------------------  
!> @author Johannes Gebert - HLRS - NUM - gebert@hlrs.de
!
!> @brief
!> Aborts if err > 0. err is required to call this routine after another one 
!> with a status feedback.
!
!> @param[in] fh Handle of file to print to
!> @param[in] txt Error message to print
!> @param[in] err Errorcode / status of the message
!------------------------------------------------------------------------------  
SUBROUTINE print_err_stop(fh, txt, err) ! , pro_path, pro_name

INTEGER(KIND=ik), INTENT(IN) :: fh , err
CHARACTER(LEN=*), INTENT(IN) :: txt

IF (err /= 0) THEN
   WRITE(fh, FMT_ERR) TRIM(txt)
   WRITE(fh, FMT_ERR_STOP)
   WRITE(*,FMT_ERR) "Can't stop gracefully."
   STOP 
END IF

END SUBROUTINE print_err_stop

END MODULE user_interaction