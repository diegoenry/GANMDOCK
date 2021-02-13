!@author
!> Diego Enry B. Gomes
!@brief
!>Error Handling
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Error(n,arg1,arg2)
use colors
integer,intent(in) :: n
character(len=30)    :: arg1,arg2
write(*,*) "---------------------------------------------------------------"
write(*,*) "Program stopped with the following error:"

select case(n)
  case(1)
    write(*,*) "No arguments."

  case(2)
    write(*,*) "No arguments for ",red,arg1

  case(3)
    write(*,*) "ERROR: Improper value for variable '",green,trim(arg1),normal,"' => '",red,trim(arg2),normal,"'"

  case(4)
    write(*,*) "ERROR: -nmrec and -nmlig mandatory" !by now
end select

write(*,*) normal
write(*,*) "run "//red//"./gamini"//green//" -h"//normal//" for help"

STOP
end subroutine

