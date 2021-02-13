!######################################################################
SUBROUTINE Dock_ChoosePair(current,&
                           receptor,&
                           ligand,&
                           nPersonReceptor,&
                           nPersonLigand,&
                           nPersonDock)
use dock_vectors
integer,intent(in)   :: current        !current value in loop that calls choose pair.
!integer,intent(in)   :: nPerson
integer,intent(in)   :: nPersonReceptor
integer,intent(in)   :: nPersonLigand
integer,intent(out)  :: receptor
integer,intent(out)  :: ligand

10 continue
! select the structures to dock
receptor  =   int(rand(0)*(nPersonReceptor))+1
ligand    =   int(rand(0)*(nPersonLigand))+1

! This will garantee we're not repeating individuals
! However it will be an infinite loop if population size bigger than all avaible
! combinations ( n_receptor * n_ligand ).

ndupes=0

do i=1,current-1
  if (pop(i)%receptor==receptor.and.pop(i)%ligand==ligand) then
    write(*,*) "Found a duplicate",current 
    ndupes=ndupes+1
    if (ndupes>=nPersonDock) then 
      write(*,*) "ndupes=nPerson"
      STOP
    endif
    goto 10
  endif
enddo

END subroutine
!######################################################################
