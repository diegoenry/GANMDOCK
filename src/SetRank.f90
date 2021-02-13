! Rank structures by fitness but does not sort
SUBROUTINE SetRank(generation)
use vectors, only : nPersonDock
use dock_vectors
real  ::  max_energy  ! it is the worst energy
real  ::  fitness     ! it is the worst energy

if (generation==0) then
  max_energy     =   maxval( pop(:)%energy(1) )
  do i=1,nPersonDock
    fitness = max_energy - pop(i)%energy(1)
    pop(i)%fitness = abs (fitness)
  enddo
  fitness_sum = sum(pop(:)%fitness)
endif

if (generation>0) then
  max_energy     =   maxval( pop_tmp(:)%energy(1) )
  do i=1,nPersonDock*2
    fitness = max_energy - pop_tmp(i)%energy(1)
    pop_tmp(i)%fitness = abs (fitness)
  enddo
  fitness_sum = sum(pop_tmp(:)%fitness)
endif
end subroutine SetRank
