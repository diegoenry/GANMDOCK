!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine crossover
use dock_vectors
use vectors, only : nPersonDock
integer ::  crossover_1
integer ::  crossover_2
integer ::  posicao

write(*,"(A)") "crossover"

do i=1,nPersonDock,2

  call roulette(crossover_1,crossover_2)

  !crossover
  pop_tmp(i+nPersonDock)%receptor    =   pop(crossover_1)%receptor
  pop_tmp(i+nPersonDock)%ligand      =   pop(crossover_2)%ligand

  pop_tmp(i+nPersonDock+1)%receptor  =   pop(crossover_2)%receptor
  pop_tmp(i+nPersonDock+1)%ligand    =   pop(crossover_1)%ligand

enddo

write(*,*) "Parents"
do i=1,nPersonDock
  write(*,*) pop_tmp(i)%receptor,pop_tmp(i)%ligand
enddo

write(*,*) "offspring"
do i=1,nPersonDock
 write(*,*) pop_tmp(i+nPersonDock)%receptor,pop_tmp(i+nPersonDock)%ligand
enddo

END subroutine crossover
