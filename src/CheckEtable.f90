!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine CheckEnergyTable(generation)
use vectors, only : nPerson,EnergyTableLast,nPersonLigand,nPersonReceptor,nPersonDock
use dock_vectors
integer,intent(in) :: generation
integer  :: receptor,ligand

!reset excluded list
exclude(:)=0 

if (generation==0) then
  do i=1,nPersonDock
    receptor = pop(i)%receptor
    ligand   = pop(i)%ligand
    do j=1,EnergyTableLast
      if (EnergyTable(j)%receptor==receptor.and.EnergyTable(j)%ligand==ligand) then
        pop(i)=EnergyTable(j) ! transfers the whole object
        exclude(i)=1
        write(*,"(i5,i5,A)") receptor,ligand, " found in the Energy Table"
        exit
      endif
    enddo
  enddo

else

  do i=nPersonDock+1,nPersonDock*2
    receptor = pop_tmp(i)%receptor
    ligand   = pop_tmp(i)%ligand
    do j=1,EnergyTableLast
      if (EnergyTable(j)%receptor==receptor.and.EnergyTable(j)%ligand==ligand) then
        pop_tmp(i)=EnergyTable(j) ! transfers the whole object
        exclude(i)=1
        write(*,"(i5,i5,A)") receptor,ligand, " found in the Energy Table"
        exit
      endif
    enddo
  enddo

endif

end subroutine CheckEnergyTable