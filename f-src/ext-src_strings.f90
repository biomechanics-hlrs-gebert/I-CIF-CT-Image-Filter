!------------------------------------------------------------------------------
! HLRS - Numerical Methods and Libraries, Nobelstraße 19, 70569 Stumgart
!------------------------------------------------------------------------------
!
! MODULE: strings
!
!> @author George Benthien
!> @author Johannes Gebert, gebert@hlrs.de, HLRS/NUM
!
!> @Descrition: 
!> Several functions to get some string work done. Mainly developed by 
!> George Benthien
!> http://numerical.recipes/forum/showthread.php?t=1659
!> http://www.gbenthien.net/strings/index.html
!> slight modifications by Johannes Gebert on 06.03.2021
!------------------------------------------------------------------------------
MODULE strings

USE global_std

CONTAINS

subroutine parse(str,delims,args,nargs)

! Parses the string 'str' into arguments args(1), ..., args(nargs) based on
! the delimiters contained in the string 'delims'. Preceding a delimiter in
! 'str' by a backslash (\) makes this particular instance not a delimiter.
! The integer output variable nargs contains the number of arguments found.

character(len=*) :: str,delims
character(len=len_trim(str)) :: strsav
character(len=mcl),dimension(:) :: args
integer(kind=ik) :: nargs

strsav=str
call compact(str)
na=size(args)
do i=1,na
  args(i)=' '
end do
nargs=0
lenstr=len_trim(str)
if(lenstr==0) return
k=0

do
   if(len_trim(str) == 0) exit
   nargs=nargs+1
   call split(str,delims,args(nargs))
   call removebksl(args(nargs))
end do
str=strsav

end subroutine parse

!----------------------------------------------------------------------------

subroutine compact(str)

! Converts multiple spaces and tabs to single spaces; deletes control characters;
! removes initial spaces.

character(len=*):: str
character(len=1):: ch
character(len=len_trim(str)):: outstr

str=adjustl(str)
lenstr=len_trim(str)
outstr=' '
isp=0
k=0

do i=1,lenstr
  ch=str(i:i)
  ich=iachar(ch)

  select case(ich)

    case(9,32)     ! space or tab character
      if(isp==0) then
        k=k+1
        outstr(k:k)=' '
      end if
      isp=1

    case(33:)      ! not a space, quote, or control character
      k=k+1
      outstr(k:k)=ch
      isp=0

  end select

end do

str=adjustl(outstr)

end subroutine compact

!----------------------------------------------------------------------------

subroutine removesp(str)

! Removes spaces, tabs, and control characters in string str

character(len=*):: str
character(len=1):: ch
character(len=len_trim(str))::outstr

str=adjustl(str)
lenstr=len_trim(str)
outstr=' '
k=0

do i=1,lenstr
  ch=str(i:i)
  ich=iachar(ch)
  select case(ich)
    case(0:32)  ! space, tab, or control character
         cycle
    case(33:)
      k=k+1
      outstr(k:k)=ch
  end select
end do

str=adjustl(outstr)

end subroutine removesp

!----------------------------------------------------------------------------

subroutine value_dr(str,rnum,ios)

! Converts number string to a double precision real number

character(len=*)::str
real(rk)::rnum
integer :: ios

ilen=len_trim(str)
ipos=scan(str,'Ee')
if(.not.is_digit(str(ilen:ilen)) .and. ipos/=0) then
   ios=3
   return
end if
read(str,*,iostat=ios) rnum

end subroutine value_dr

!----------------------------------------------------------------------------

subroutine value_sr(str,rnum,ios)

! Converts number string to a single precision real number

character(len=*)::str
real(rk) :: rnum
real(rk) :: rnumd 

call value_dr(str,rnumd,ios)
if( abs(rnumd) > huge(rnum) ) then
  ios=15
  return
end if
if( abs(rnumd) < tiny(rnum) ) rnum=0.0_rk
rnum=rnumd

end subroutine value_sr

!----------------------------------------------------------------------------

subroutine value_di(str,inum,ios)

! Converts number string to a double precision integer value

character(len=*)::str
integer(ik) :: inum
real(rk) :: rnum

call value_dr(str,rnum,ios)
if(abs(rnum)>huge(inum)) then
  ios=15
  return
end if
inum=nint(rnum,ik)

end subroutine value_di

!----------------------------------------------------------------------------

subroutine value_si(str,inum,ios)

! Converts number string to a single precision integer value

character(len=*)::str
integer(ik) :: inum
real(rk) :: rnum

call value_dr(str,rnum,ios)
if(abs(rnum)>huge(inum)) then
  ios=15
  return
end if
inum=nint(rnum,ik)

end subroutine value_si

!----------------------------------------------------------------------------

subroutine shiftstr(str,n)

! Shifts characters in in the string 'str' n positions (positive values
! denote a right shift and negative values denote a left shift). Characters
! that are shifted off the end are lost. Positions opened up by the shift 
! are replaced by spaces.

character(len=*):: str

lenstr=len(str)
nabs=iabs(n)
if(nabs>=lenstr) then
  str=repeat(' ',lenstr)
  return
end if
if(n<0) str=str(nabs+1:)//repeat(' ',nabs)  ! shift left
if(n>0) str=repeat(' ',nabs)//str(:lenstr-nabs)  ! shift right 
return

end subroutine shiftstr

! !----------------------------------------------------------------------------

subroutine insertstr(str,strins,loc)

! Inserts the string 'strins' into the string 'str' at position 'loc'. 
! Characters in 'str' starting at position 'loc' are shifted right to
! make room for the inserted string. Trailing spaces of 'strins' are 
! removed prior to insertion

character(len=*):: str,strins
character(len=len(str))::tempstr

lenstrins=len_trim(strins)
tempstr=str(loc:)
call shiftstr(tempstr,lenstrins)
tempstr(1:lenstrins)=strins(1:lenstrins)
str(loc:)=tempstr
return

end subroutine insertstr

!----------------------------------------------------------------------------*

subroutine write_dr(rnum,str,fmt)

! Writes double precision real number rnum to string str using format fmt

real(rk) :: rnum
character(len=*) :: str,fmt
character(len=80) :: formt

formt='('//trim(fmt)//')'
write(str,formt) rnum
str=adjustl(str)

end subroutine write_dr

!----------------------------------------------------------------------------*

subroutine write_sr(rnum,str,fmt)

! Writes single precision real number rnum to string str using format fmt

real(rk) :: rnum
character(len=*) :: str,fmt
character(len=80) :: formt

formt='('//trim(fmt)//')'
write(str,formt) rnum
str=adjustl(str)

end subroutine write_sr

!----------------------------------------------------------------------------*

subroutine write_di(inum,str,fmt)

! Writes double precision integer inum to string str using format fmt

integer(ik) :: inum
character(len=*) :: str,fmt
character(len=80) :: formt

formt='('//trim(fmt)//')'
write(str,formt) inum
str=adjustl(str)

end subroutine write_di

!----------------------------------------------------------------------------*

subroutine write_si(inum,str,fmt)

! Writes single precision integer inum to string str using format fmt

integer(ik) :: inum
character(len=*) :: str,fmt
character(len=80) :: formt

formt='('//trim(fmt)//')'
write(str,formt) inum
str=adjustl(str)

end subroutine write_si

!----------------------------------------------------------------------------*

subroutine trimzero(str)

! Deletes nonsignificant trailing zeroes from number string str. If number
! string ends in a decimal point, one trailing zero is added.

character(len=*) :: str
character :: ch
character(len=10) :: exp

ipos=scan(str,'eE')
if(ipos>0) then
   exp=str(ipos:)
   str=str(1:ipos-1)
endif
lstr=len_trim(str)
do i=lstr,1,-1
   ch=str(i:i)
   if(ch=='0') cycle          
   if(ch=='.') then
      str=str(1:i)//'0'
      if(ipos>0) str=trim(str)//trim(exp)
      exit
   endif
   str=str(1:i)
   exit
end do
if(ipos>0) str=trim(str)//trim(exp)

end subroutine trimzero

!----------------------------------------------------------------------------

function is_digit(ch) result(res)

! Returns .true. if ch is a digit (0,1,...,9) and .false. otherwise

character :: ch
logical :: res

select case(ch)
case('0':'9')
  res=.true.
case default
  res=.false.
end select
return

end function is_digit

!----------------------------------------------------------------------------

subroutine split(str,delims,before,sep)

! Routine finds the first instance of a character from 'delims' in the
! the string 'str'. The characters before the found delimiter are
! output in 'before'. The characters after the found delimiter are
! output in 'str'. The optional output character 'sep' contains the 
! found delimiter. A delimiter in 'str' is treated like an ordinary 
! character if it is preceded by a backslash (\). If the backslash 
! character is desired in 'str', then precede it with another backslash.

character(len=*) :: str,delims
character,optional :: sep
logical :: pres
character :: ch,cha
character(len=mcl) :: before

pres=present(sep)
str=adjustl(str)
call compact(str)
lenstr=len_trim(str)
if(lenstr == 0) return        ! string str is empty
k=0
ibsl=0                        ! backslash initially inactive
before=' '
do i=1,lenstr
   ch=str(i:i)
   if(ibsl == 1) then          ! backslash active
      k=k+1
      before(k:k)=ch
      ibsl=0
      cycle
   end if
   if(ch == '\\' ) then          ! backslash with backslash inactive
      k=k+1
      before(k:k)=ch
      ibsl=1
      cycle
   end if
   ipos=index(delims,ch)
   if(ipos == 0) then          ! character is not a delimiter
      k=k+1
      before(k:k)=ch
      cycle
   end if
   if(ch /= ' ') then          ! character is a delimiter that is not a space
      str=str(i+1:)
      if(pres) sep=ch
      exit
   end if
   cha=str(i+1:i+1)            ! character is a space delimiter
   iposa=index(delims,cha)
   if(iposa > 0) then          ! next character is a delimiter
      str=str(i+2:)
      if(pres) sep=cha
      exit
   else
      str=str(i+1:)
      if(pres) sep=ch
      exit
   end if
end do
if(i >= lenstr) str=''
str=adjustl(str)              ! remove initial spaces
return

end subroutine split

!----------------------------------------------------------------------------

subroutine removebksl(str)

! Removes backslash (\) characters. Double backslashes (\\) are replaced
! by a single backslash.

character(len=*):: str
character(len=1):: ch
character(len=len_trim(str))::outstr

str=adjustl(str)
lenstr=len_trim(str)
outstr=' '
k=0
ibsl=0                        ! backslash initially inactive

do i=1,lenstr
  ch=str(i:i)
  if(ibsl == 1) then          ! backslash active
   k=k+1
   outstr(k:k)=ch
   ibsl=0
   cycle
  end if
  if(ch == '\\') then          ! backslash with backslash inactive
   ibsl=1
   cycle
  end if
  k=k+1
  outstr(k:k)=ch              ! non-backslash with backslash inactive
end do

str=adjustl(outstr)

end subroutine removebksl

!----------------------------------------------------------------------------

END MODULE strings  


