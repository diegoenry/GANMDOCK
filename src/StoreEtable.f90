subroutine StoreEnergytable(generation)
use vectors!, only : nPerson,EnergyTableLast
use dock_vectors
implicit none
integer             :: i
integer,intent(in)  :: generation

if (generation==0) then
  do i=1,nPersonDock
    if (exclude(i)==0) then
      EnergyTableLast=EnergyTableLast+1
      EnergyTable(EnergyTableLast)=pop(i) ! transfers the whole object
      write(55,*) EnergyTable(EnergyTableLast)
      write(*,*) "Etable",EnergyTableLast,"rec",EnergyTable(EnergyTableLast)%receptor,"lig",EnergyTable(EnergyTableLast)%ligand
   endif
  enddo

else

  do i=nPersonDock+1,nPersonDock*2
    if (exclude(i)==0) then
      EnergyTableLast=EnergyTableLast+1
      EnergyTable(EnergyTableLast)=pop_tmp(i) ! transfers the whole object
      write(55,*) EnergyTable(EnergyTableLast)
      write(*,*) "Etable",EnergyTableLast,"rec",EnergyTable(EnergyTableLast)%receptor,"lig",EnergyTable(EnergyTableLast)%ligand
    endif
  enddo

endif

end subroutine StoreEnergytable