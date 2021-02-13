! This is a simple fortran 90 program to read any PDB
! and split the chains by writing "TER" between them.
!
! Written by
! Diego Enry B. Gomes (1,2) | diego@biof.ufrj.br
! Elodie Japaud (2,3) | 
! David Perahia (2)   |
! 
! 1) Universidade Federal do Rio de Janeiro, Rio de Janeiro - Brazil
! 2) Universite Paris Diderot,   Paris    - France
! 3) Ecole Normale Superieure de Cachan,     Cachan   - France
!
! Wed Apr 20 18:20:09 CEST 2011

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program fixpdb
character(len=80)  :: line
character(len=80)  :: inputfile
character(len=1)   :: chain

! gets the argument from the command line
CALL GETARG(1 , inputfile)

! uses "grep" to remove all unwanted lines
!call system("grep ATOM "//inputfile//"> "//trim(inputfile)//".tmp.pdb")

open(1,file=trim(inputfile))

! reads the first chain identifier
read(1,"(a80)") line
chain=line(22:22) 
backspace(1)

! does the magic
do
  read(1,"(a80)",err=100,end=100) line 

    if ( line(22:22) /= chain ) then

     write(*,"(a3)") "TER" 
   chain=line(22:22)

    endif

write(*,"(a80)") line 

enddo
100 continue
    write(*,"(a3)") "END"

end

