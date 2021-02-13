!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine GADock
use vectors
use dock_vectors

open(66,file="docklogfile.dat")

crossover_method=1
EnergyTableLast=0

if ( (nPersonReceptor*nPersonLigand) < nPersonDock ) then
 nPersonDock=(nPersonReceptor*nPersonLigand)
endif

allocate(pop(nPersonDock))
allocate(pop_tmp(nPersonDock*2))
allocate(EnergyTable(nPersonDock*6)) 
allocate(exclude(nPersonDock*2))
allocate(duplicate(2,nPersonDock*2))

!reset EnergyTable
do i=1,nPersonDock*6
  EnergyTable(i)%receptor  = 0
  EnergyTable(i)%ligand    = 0
  EnergyTable(i)%fitness   = 0.0
  do j=1,10
  EnergyTable(i)%energy(j) = 0.0
  EnergyTable(i)%rms(j)    = 0.0
  enddo
enddo

if (EnergyTableFile=="") then
  EnergyTableFile="etable.dat"
  EnergyTableLast=0
endif
!Open and read energy table
open(55,file=EnergyTableFile)
!count number of lines
  do 100
    read(55,*,end=200,err=200) i
    EnergyTableLast=EnergyTableLast+1
100 continue
200 rewind(55)

  do i=1,EnergyTableLast
    read(55,*) EnergyTable(i)
  enddo

open(2,file="gadock.dat")
do i=0,ngenerations

  write(*,*) "Generation: ",i

  call Dock_GeneratePopulation(i)

  ! call StoreEnergytable(i) !must come before SortFitness

  call WriteLog(i)

enddo

close(2)

end subroutine GADock


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine CheckDuplicate
use vectors, only : nPersonDock
use dock_vectors

do i=1,nPersonDock*2
  duplicate(1,i) = 0
  duplicate(2,i) = 0
enddo

do i=nPersonDock+2,nPersonDock*2 !nao compara o 1o elemento neh
  do j=nPersonDock+1,i-1
    if (pop_tmp(i)%receptor==pop_tmp(j)%receptor.and.pop_tmp(i)%ligand==pop_tmp(j)%ligand) then
      write(*,*) "Found a duplicate in pop_tmp, position =",j," pair:",pop_tmp(j)%receptor,pop_tmp(j)%ligand
      duplicate(1,i)=i  !
      duplicate(2,i)=j  !
      exclude(i)=1      !exclude calculation
      goto 10
    endif
  enddo
10 continue
enddo

end subroutine CheckDuplicate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine GetDuplicateResult
use vectors, only : nPersonDock
use dock_vectors

do i=nPersonDock+2,nPersonDock*2
  if (duplicate(1,i)/=0) then
    pop_tmp(i) = pop_tmp(duplicate(2,i) ) !transfers whole object
  endif
enddo

end subroutine GetDuplicateResult
