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
! Formats
!------------------------------------------------------------------------------
CHARACTER(Len=*), PARAMETER :: FMT_SEP = "(80('-'))"
CHARACTER(Len=*), PARAMETER :: TAB_WDTH = "4"
CHARACTER(Len=*), PARAMETER :: FMT_INT = "I0"
CHARACTER(Len=*), PARAMETER :: FMT_SCI = "E0.6"
CHARACTER(Len=*), PARAMETER :: FMT_REAL = "F0.6"
!
CHARACTER(Len=*), PARAMETER :: TXT = "('-- ',"
CHARACTER(Len=*), PARAMETER :: DBG = "('DD ',"
CHARACTER(Len=*), PARAMETER :: MSG = "('MM ',"
CHARACTER(Len=*), PARAMETER :: WRN = "('WW ',"
CHARACTER(Len=*), PARAMETER :: ERR = "('EE ',"
!
CHARACTER(LEN=*), PARAMETER :: FMT     = "*(A))"
!
CHARACTER(LEN=*), PARAMETER :: AI0xAI0 = "(A,"//FMT_INT//",1x,*(A,1x,"//FMT_INT//")))"
CHARACTER(LEN=*), PARAMETER :: AI0AxI0 = "(A,"//FMT_INT//",1x,A,*(1x,"//FMT_INT//")))"
CHARACTER(LEN=*), PARAMETER :: xAI0    = "*(A,1x,"//FMT_INT//",1x))"
CHARACTER(LEN=*), PARAMETER :: AxI0    = "*(A,1x,*("//FMT_INT//",1x)))"
!
CHARACTER(LEN=*), PARAMETER :: AI0xAF0 = "(A,"//FMT_INT//",1x,*(A,1x,"//FMT_REAL//")))"
CHARACTER(LEN=*), PARAMETER :: AI0AxF0 = "(A,"//FMT_INT//",1x,A,*(1x,"//FMT_REAL//")))"
CHARACTER(LEN=*), PARAMETER :: xAF0    = "*(A,1x,"//FMT_REAL//",1x))"
CHARACTER(LEN=*), PARAMETER :: AxF0    = "*(A,1x,*("//FMT_REAL//",1x)))"
!
CHARACTER(LEN=*), PARAMETER :: AI0xAE0 = "(A,"//FMT_INT//",1x,*(A,1x,"//FMT_SCI//")))"
CHARACTER(LEN=*), PARAMETER :: AI0AxE0 = "(A,"//FMT_INT//",1x,A,*(1x,"//FMT_SCI//")))"
CHARACTER(LEN=*), PARAMETER :: xAE0    = "*(A,1x,"//FMT_SCI//",1x))"
CHARACTER(LEN=*), PARAMETER :: AxE0    = "*(A,1x,*("//FMT_SCI//",1x)))"
!
CHARACTER(LEN=*), PARAMETER :: xAL     = "*(A,1x,L1,1x))"
!
!------------------------------------------------------------------------------
! The following formats are wrappers to mask a direct use of the format 
! string concatenation. 
!
! For example
! WRITE(*, FMT_ERR_xAI0) "Lorem Ipsum" and
! WRITE(*, ERR//xAI0) "Lorem Ipsum" result in the same output.
!------------------------------------------------------------------------------
! Error formats
!------------------------------------------------------------------------------
CHARACTER(Len=*), PARAMETER :: FMT_ERR_STOP = "('EE PROGRAM STOPPED.')"
!
CHARACTER(Len=*), PARAMETER :: FMT_ERR         = ERR//FMT
CHARACTER(Len=*), PARAMETER :: FMT_ERR_SEP     = FMT_SEP ! "('EE ',80('='))"
!
CHARACTER(Len=*), PARAMETER :: FMT_ERR_AI0xAI0 = ERR//AI0xAI0
CHARACTER(Len=*), PARAMETER :: FMT_ERR_AI0AxI0 = ERR//AI0AxI0
CHARACTER(Len=*), PARAMETER :: FMT_ERR_xAI0    = ERR//xAI0
CHARACTER(Len=*), PARAMETER :: FMT_ERR_AxI0    = ERR//AxI0
!
CHARACTER(Len=*), PARAMETER :: FMT_ERR_AI0xAF0 = ERR//AI0xAF0
CHARACTER(Len=*), PARAMETER :: FMT_ERR_AI0AxF0 = ERR//AI0AxF0
CHARACTER(Len=*), PARAMETER :: FMT_ERR_xAF0    = ERR//xAF0
CHARACTER(Len=*), PARAMETER :: FMT_ERR_AxF0    = ERR//AxF0
!
CHARACTER(Len=*), PARAMETER :: FMT_ERR_AI0xAE0 = ERR//AI0xAE0
CHARACTER(Len=*), PARAMETER :: FMT_ERR_AI0AxE0 = ERR//AI0AxE0
CHARACTER(Len=*), PARAMETER :: FMT_ERR_xAE0    = ERR//xAE0
CHARACTER(Len=*), PARAMETER :: FMT_ERR_AxE0    = ERR//AxE0
!
CHARACTER(Len=*), PARAMETER :: FMT_ERR_xAL     = ERR//xAL

!------------------------------------------------------------------------------
! Text formats
!------------------------------------------------------------------------------
CHARACTER(Len=*), PARAMETER :: FMT_TXT         = TXT//FMT
CHARACTER(Len=*), PARAMETER :: FMT_TXT_SEP     = FMT_SEP ! "('-- ',80('-'))"
!
CHARACTER(Len=*), PARAMETER :: FMT_TXT_AI0xAI0 = TXT//AI0xAI0
CHARACTER(Len=*), PARAMETER :: FMT_TXT_AI0AxI0 = TXT//AI0AxI0
CHARACTER(Len=*), PARAMETER :: FMT_TXT_xAI0    = TXT//xAI0
CHARACTER(Len=*), PARAMETER :: FMT_TXT_AxI0    = TXT//AxI0
!
CHARACTER(Len=*), PARAMETER :: FMT_TXT_AI0xAF0 = TXT//AI0xAF0
CHARACTER(Len=*), PARAMETER :: FMT_TXT_AI0AxF0 = TXT//AI0AxF0
CHARACTER(Len=*), PARAMETER :: FMT_TXT_xAF0    = TXT//xAF0
CHARACTER(Len=*), PARAMETER :: FMT_TXT_AxF0    = TXT//AxF0
!
CHARACTER(Len=*), PARAMETER :: FMT_TXT_AI0xAE0 = TXT//AI0xAE0
CHARACTER(Len=*), PARAMETER :: FMT_TXT_AI0AxE0 = TXT//AI0AxE0
CHARACTER(Len=*), PARAMETER :: FMT_TXT_xAE0    = TXT//xAE0
CHARACTER(Len=*), PARAMETER :: FMT_TXT_AxE0    = TXT//AxE0
!
CHARACTER(Len=*), PARAMETER :: FMT_TXT_xAL     = TXT//xAL

!------------------------------------------------------------------------------
! Message/debug formats
!------------------------------------------------------------------------------
CHARACTER(Len=*), PARAMETER :: FMT_MSG         = MSG//FMT
CHARACTER(Len=*), PARAMETER :: FMT_MSG_SEP     = FMT_SEP ! "('MM ',80('-'))"
!
CHARACTER(Len=*), PARAMETER :: FMT_MSG_AI0xAI0 = MSG//AI0xAI0
CHARACTER(Len=*), PARAMETER :: FMT_MSG_AI0AxI0 = MSG//AI0AxI0
CHARACTER(Len=*), PARAMETER :: FMT_MSG_xAI0    = MSG//xAI0
CHARACTER(Len=*), PARAMETER :: FMT_MSG_AxI0    = MSG//AxI0
!
CHARACTER(Len=*), PARAMETER :: FMT_MSG_AI0xAF0 = MSG//AI0xAF0
CHARACTER(Len=*), PARAMETER :: FMT_MSG_AI0AxF0 = MSG//AI0AxF0
CHARACTER(Len=*), PARAMETER :: FMT_MSG_xAF0    = MSG//xAF0
CHARACTER(Len=*), PARAMETER :: FMT_MSG_AxF0    = MSG//AxF0
!
CHARACTER(Len=*), PARAMETER :: FMT_MSG_AI0xAE0 = MSG//AI0xAE0
CHARACTER(Len=*), PARAMETER :: FMT_MSG_AI0AxE0 = MSG//AI0AxE0
CHARACTER(Len=*), PARAMETER :: FMT_MSG_xAE0    = MSG//xAE0
CHARACTER(Len=*), PARAMETER :: FMT_MSG_AxE0    = MSG//AxE0
!
CHARACTER(Len=*), PARAMETER :: FMT_MSG_xAL     = MSG//xAL

!------------------------------------------------------------------------------
! Warning formats
!------------------------------------------------------------------------------
CHARACTER(Len=*), PARAMETER :: FMT_WRN         = WRN//FMT
CHARACTER(Len=*), PARAMETER :: FMT_WRN_SEP     = FMT_SEP ! "('WW ',80('-'))"
!
CHARACTER(Len=*), PARAMETER :: FMT_WRN_AI0xAI0 = WRN//AI0xAI0
CHARACTER(Len=*), PARAMETER :: FMT_WRN_AI0AxI0 = WRN//AI0AxI0
CHARACTER(Len=*), PARAMETER :: FMT_WRN_xAI0    = WRN//xAI0
CHARACTER(Len=*), PARAMETER :: FMT_WRN_AxI0    = WRN//AxI0
!
CHARACTER(Len=*), PARAMETER :: FMT_WRN_AI0xAF0 = WRN//AI0xAF0
CHARACTER(Len=*), PARAMETER :: FMT_WRN_AI0AxF0 = WRN//AI0AxF0
CHARACTER(Len=*), PARAMETER :: FMT_WRN_xAF0    = WRN//xAF0
CHARACTER(Len=*), PARAMETER :: FMT_WRN_AxF0    = WRN//AxF0
!
CHARACTER(Len=*), PARAMETER :: FMT_WRN_AI0xAE0 = WRN//AI0xAE0
CHARACTER(Len=*), PARAMETER :: FMT_WRN_AI0AxE0 = WRN//AI0AxE0
CHARACTER(Len=*), PARAMETER :: FMT_WRN_xAE0    = WRN//xAE0
CHARACTER(Len=*), PARAMETER :: FMT_WRN_AxE0    = WRN//AxE0
!
CHARACTER(Len=*), PARAMETER :: FMT_WRN_xAL     = WRN//xAL

!------------------------------------------------------------------------------
! Debug formats
!------------------------------------------------------------------------------
CHARACTER(Len=*), PARAMETER :: FMT_DBG         = DBG//FMT
CHARACTER(Len=*), PARAMETER :: FMT_DBG_SEP     = FMT_SEP ! "('DD ',80('-'))"
!
CHARACTER(Len=*), PARAMETER :: FMT_DBG_AI0xAI0 = DBG//AI0xAI0
CHARACTER(Len=*), PARAMETER :: FMT_DBG_AI0AxI0 = DBG//AI0AxI0
CHARACTER(Len=*), PARAMETER :: FMT_DBG_xAI0    = DBG//xAI0
CHARACTER(Len=*), PARAMETER :: FMT_DBG_AxI0    = DBG//AxI0
!
CHARACTER(Len=*), PARAMETER :: FMT_DBG_AI0xAF0 = DBG//AI0xAF0
CHARACTER(Len=*), PARAMETER :: FMT_DBG_AI0AxF0 = DBG//AI0AxF0
CHARACTER(Len=*), PARAMETER :: FMT_DBG_xAF0    = DBG//xAF0
CHARACTER(Len=*), PARAMETER :: FMT_DBG_AxF0    = DBG//AxF0
!
CHARACTER(Len=*), PARAMETER :: FMT_DBG_AI0xAE0 = DBG//AI0xAE0
CHARACTER(Len=*), PARAMETER :: FMT_DBG_AI0AxE0 = DBG//AI0AxE0
CHARACTER(Len=*), PARAMETER :: FMT_DBG_xAE0    = DBG//xAE0
CHARACTER(Len=*), PARAMETER :: FMT_DBG_AxE0    = DBG//AxE0
!
CHARACTER(Len=*), PARAMETER :: FMT_DBG_xAL     = DBG//xAL

!------------------------------------------------------------------------------
! Provide colors on std_out (!) 
! Needs to compile with -fbackslash 
! Use of a requires resetting it.
! Will interfere with exporting "\", especially in the context of *.tex
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
CHARACTER(LEN=*), INTENT(OUT) :: restart, cmd_arg_history
LOGICAL :: stp

CHARACTER(LEN=mcl) :: cmd_arg
INTEGER(KIND=ik) :: ii

stp = .FALSE.
restart=''

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
WRITE(std_out, FMT_TXT) TRIM(ADJUSTL(longname))
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
!> Stop a program.
!
!> @description
!> Aborts non-gracefull with a stop on one processor if err > 0. 
!> Err /= 0 is required to call this routine after another one 
!> with a status feedback.
!
!> @param[in] fh Handle of file to print to
!> @param[in] text Error message to print
!> @param[in] error Errorcode / status of the message
!------------------------------------------------------------------------------  
SUBROUTINE print_err_stop(fh, text, error) ! , pro_path, pro_name

INTEGER(KIND=ik), INTENT(IN) :: fh , error
CHARACTER(LEN=*), INTENT(IN) :: text

IF (error > 0) THEN
   WRITE(fh, FMT_ERR) TRIM(text)
   WRITE(fh, FMT_ERR_STOP)
   WRITE(fh,FMT_ERR) "Can't stop gracefully."
   STOP 
END IF

END SUBROUTINE print_err_stop
  

!------------------------------------------------------------------------------
! SUBROUTINE: estimated_time_of_arrival
!------------------------------------------------------------------------------  
!> @author Johannes Gebert - HLRS - NUM - gebert@hlrs.de
!
!> @brief
!> Calculates the estimated time of arrival. Best suited for Julius or other 
!> interactive notebooks/desktops/...
!
!> @param[in] sec Handle of file to print to
!> @param[out] string Output string of the ETA. Printed to std_out if not given
!------------------------------------------------------------------------------  
SUBROUTINE estimated_time_of_arrival(sec, string)

    INTEGER(KIND=ik), INTENT(IN) :: sec
    CHARACTER(LEN=scl), INTENT(OUT), OPTIONAL :: string

    INTEGER(KIND=ik) :: mins, hours, secs, seconds, mins_s, remainder
    CHARACTER(LEN=scl) :: str

    mins_s    = MODULO(sec,  3600_ik)
    seconds   = MODULO(mins_s, 60_ik)
    remainder = MODULO(seconds, 1_ik)

    hours = (sec-mins_s) / 3600_ik
    mins  = (mins_s-seconds) / 60_ik
    secs  = (seconds-remainder)

    write(string,'(I3,2(A,I2),A)') hours,":",mins,":",secs," hhh:mm:ss"

    IF(PRESENT(string)) THEN
        string = str
    ELSE
        WRITE(std_out, FMT_TXT) "ETA: "//TRIM(ADJUSTL(str))
    END IF

END SUBROUTINE estimated_time_of_arrival

END MODULE user_interaction