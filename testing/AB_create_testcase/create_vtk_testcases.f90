  !---------------------------------------------------------------------------------------------------
  !> \mainpage HLRS - vtk STURCUTRED_POINTS data manipulation tool
  !>
  !> <hr>
  !> \section desc Description
  !>
  !> This program creates testcases to check whether specific parts work properly.
  !> It's more like a script which needs to be modified according to your specific needs.
  !>
  !> <hr>
  !> \section developers Developers
  !>
  !> Johannes Gebert
  !>
  !>  \section modified Last modified:
  !>  by: Johannes Gebert \n
  !>  on: 20.01.2021
  !---------------------------------------------------------------------------------------------------


PROGRAM create_scalar_vtk_cases

  USE iso_fortran_env

  IMPLICIT NONE

  INTEGER      (KIND=INT64) , PARAMETER                       :: ik=8, rk=8
  INTEGER      (KIND=INT16) , DIMENSION (:,:,:), ALLOCATABLE  :: array
  INTEGER      (KIND=INT64) , DIMENSION (3)                   :: dims, flh
  REAL         (KIND=REAL64), DIMENSION (3)                   :: spcng
  INTEGER      (KIND=INT64)                                   :: ii, xx, yy, zz
  INTEGER      (KIND=INT64)                                   :: xxx, yyy, zzz
  INTEGER      (KIND=INT64)                                   :: steps   , step_width_sf
  INTEGER      (KIND=INT64)                                   ::           step_width_sv
  INTEGER      (KIND=INT64)                                   :: remainder_sv, remainder_sf
  INTEGER      (KIND=INT64) , PARAMETER                       :: fl_in_un=22, fl_out_un=23
  CHARACTER    (LEN=512)                                      :: fl_in, fl_out, struc
  INTEGER      (KIND=INT64)                                   :: baseline_sv, max_sv, mov_val


  !-- Manual parametrization -----------------------------------------------------------------------
  fl_out       = "/home/GEB/AAB_programs/cCT_mfCT_Spatial_Registration_Tool/testing/AB_create_testcase&
       &/SR_TC1_staircase_x_0000_1000_2000.vtk"
  !-- Structural details
  dims      = (/ 800 , 800 , 800 /)                         ! vtk header
  spcng     = (/ 0.0346000, 0.0741000, 0.0425400 /)         ! vtk header
  !-- Manual parametrization -----------------------------------------------------------------------
  baseline_sv =    0_ik
  max_sv      = 2000_ik
  steps       =  800_ik

  !-- Target structure
! struc = "UNIFORM"          ! Uniform scalar field
  struc = "STAIRCASE_X"      ! Value increasing from x=0 to x=dim(1)
                             ! Make sure, steps fit into dimensions without fractions :-)
! struc = "CHEQUERED"        ! Voxel by Voxel - lo - hi - lo - hi - lo - ...
! struc = "ZEBRA"            ! Changing baseline_sv/max_sv plane by plane

  !-- Programs task
  ALLOCATE(array(dims(1),dims(2),dims(3)))

  !-- UNIFORM
  IF(struc == "UNIFORM") THEN
     array(:,:,:) = baseline_sv
  END IF

  !-- STAIRCASE_X
  IF(struc == "STAIRCASE_X") THEN
     remainder_sv  = MODULO(max_sv-baseline_sv,steps)
     step_width_sv = ((max_sv-baseline_sv)-remainder_sv)/steps

     remainder_sf  = MODULO(dims(1),steps)
     step_width_sf = (dims(1)-remainder_sf)/steps
     DO ii=1, steps
        array(1_ik+ii*step_width_sf-step_width_sf:ii*step_width_sf,:,:) = (ii-1_ik)*step_width_sv+baseline_sv
     END DO

     IF(remainder_sf/=0_ik) THEN
        array(dims(1)-remainder_sf:dims(1),:,:) = max_sv
     END IF
  END IF

  WRITE(*,'(A,I6)') "baseline_sv: ",baseline_sv
  WRITE(*,'(A,I6)') "     max_sv: ",max_sv

  !-- CHEQUERED
  IF(struc == "CHEQUERED") THEN
     array (:,:,:) = max_sv
     DO xxx=1,dims(1)-1_ik
        IF (MODULO(xxx,2_ik) == 1_ik) THEN
           xx=xxx
        ELSE
           xx=xxx+1_ik
        END IF
        DO yyy=1,dims(2)-1_ik
           IF (MODULO(yyy,2_ik) == 1_ik) THEN
              yy=yyy
           ELSE
              yy=yyy+1_ik
           END IF
           DO zzz=1,dims(3)-1_ik,2
              array(xx,yy,zzz) = baseline_sv
           END DO
        END DO
     END DO
  END IF

  !-- ZEBRA
  IF(struc=="ZEBRA")THEN
     array (:,:,:) = max_sv
     DO xx=1,dims(1),2
        array(xx,:,:)=baseline_sv
     END DO
  END IF







  !-- Write new vtk File
  OPEN(UNIT=fl_out_un, FILE=TRIM(fl_out), ACTION="write", STATUS="new")
  WRITE(fl_out_un,'(A)')            "# vtk DataFile Version 5.1"
  WRITE(fl_out_un,'(A)')            "vtk output"
  WRITE(fl_out_un,'(A)')            "BINARY"
  WRITE(fl_out_un,'(A)')            "DATASET STRUCTURED_POINTS"
  WRITE(fl_out_un,'(A,3(I5,A))')    "DIMENSIONS", dims(1)," ",dims(2)," ",dims(3),""
  WRITE(fl_out_un,'(A,3(F11.6,A))') "SPACING ", spcng(1)," ",spcng(2)," ",spcng(3)
  WRITE(fl_out_un,'(A)')            "ORIGIN 0 0 0"
  WRITE(fl_out_un,'(A, I11)')       "POINT_DATA", dims(1)*dims(2)*dims(3)
  WRITE(fl_out_un,'(A)')            "SCALARS DICOMImage short"            ! tricky part of this stuff!
  WRITE(fl_out_un,'(A)')            "LOOKUP_TABLE default"
  CLOSE(UNIT=fl_out_un)

  OPEN(UNIT=fl_out_un, FILE=TRIM(fl_out), CONVERT='big_endian', &
       ACCESS="stream", FORM="unformatted", STATUS="old", POSITION="append")
  WRITE(UNIT=fl_out_un) array(:,:,:)
  CLOSE(UNIT=fl_out_un)

  OPEN(UNIT=fl_out_un, FILE=TRIM(fl_out), ACTION="write", STATUS="old", POSITION="append")
  WRITE(fl_out_un ,'(A)')          "METADATA"
  WRITE(fl_out_un ,'(A)')          "INFORMATION 0"
  WRITE(fl_out_un ,'(A)')
  CLOSE(UNIT=fl_out_un)

END PROGRAM create_scalar_vtk_cases
