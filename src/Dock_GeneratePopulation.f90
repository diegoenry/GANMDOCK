! Generates the initial population
!########################################################################
SUBROUTINE Dock_GeneratePopulation(generation)
use dock_vectors
use vectors
integer,intent(in)   :: generation

write(*,*) "Dock_GeneratePopulation",generation

open(99,file="runhex")
if (pbs==0) write(99,1000)
if (pbs==1) write(99,2000) !nproc

if(generation==0) then

  i=1 !bias to include the Right result.
  pop(i)%receptor=1
  pop(i)%ligand=1

  do i=2,nPersonDock
    call Dock_ChoosePair(i,pop(i)%receptor ,pop(i)%ligand   ,&
                           nPersonReceptor ,nPersonLigand   ,&
                           nPersonDock)
  enddo
  
  call CheckEnergyTable(generation)

  do i=1,nPersonDock
    if (exclude(i)==0) then
      write(*,*) pop(i)%receptor,pop(i)%ligand,exclude(i),"running"
      write(99,"(A,1x,3(i5))") "./dock_and_rerank.bash ", i, pop(i)%receptor,pop(i)%ligand
    else
      write(*,*) pop(i)%receptor,pop(i)%ligand,exclude(i),"Not running"
    endif
  enddo
  write(99,"(A3)") "EOF"
  close(99)

!Actually run the dockings
call system("cat runhex")
call system("bash runhex")

call Zrank_Get_RESULTS(generation,nPersonDock)
call SetRank(generation)
call StoreEnergytable(generation) !must come before SortFitness
call SortFitness(generation)

!transfer best results from pop_tmp to pop
  do i=1,nPersonDock
    pop_tmp(i)%receptor  = pop(i)%receptor
    pop_tmp(i)%ligand    = pop(i)%ligand
    pop_tmp(i)%fitness   = pop(i)%fitness
    do j=1,10
    pop_tmp(i)%energy(j) = pop(i)%energy(j)
    pop_tmp(i)%rms(j)    = pop(i)%rms(j)
    enddo
  enddo

endif


if (generation>0) then
  write(*,*) "Dock_GeneratePopulation",generation

!Scheme of crossover when there is Elitism
  call crossover
  call CheckEnergyTable(generation)
  call CheckDuplicate

  do i=nPersonDock+1,nPersonDock*2
    if (exclude(i)==0) then
      write(*,*) pop_tmp(i)%receptor,pop_tmp(i)%ligand,exclude(i),"running"
      write(99,"(A,1x,3(i5))") "./dock_and_rerank.bash ", i-nPersonDock, pop_tmp(i)%receptor,pop_tmp(i)%ligand
    else
      write(*,*) pop_tmp(i)%receptor,pop_tmp(i)%ligand,exclude(i),"NOT running"
    endif
  enddo
  write(99,"(A3)") "EOF"
  close(99)


!Actually run the dockings
  call system("cat runhex")
  call system("bash runhex")
  call Zrank_Get_RESULTS(generation,nPersonDock)
  call GetDuplicateResult
  call SetRank(generation)
  call StoreEnergytable(generation) !must come before SortFitness
  call SortFitness(generation)

!transfer best results from pop_tmp to pop
  if (elitism==0) then
    do i=1,nPersonDock
      pop(i)%receptor  = pop_tmp(i)%receptor
      pop(i)%ligand    = pop_tmp(i)%ligand
      pop(i)%fitness   = pop_tmp(i)%fitness
      do j=1,10
      pop(i)%energy(j) = pop_tmp(i)%energy(j) 
      pop(i)%rms(j)    = pop_tmp(i)%rms(j)
      enddo
    enddo
  else 
!with elitism
    do i=1,int(nPersonDock*elitism_rate)
      pop(i)%receptor  = pop_tmp(i)%receptor
      pop(i)%ligand    = pop_tmp(i)%ligand
      pop(i)%fitness   = pop_tmp(i)%fitness

      do j=1,10
        pop(i)%energy(j) = pop_tmp(i)%energy(j)
        pop(i)%rms(j)    = pop_tmp(i)%rms(j)
      enddo

    enddo

!WHERE IS THE COMPETITION
    call RouletteElitism(nPersonDock)

  endif

endif

1000  format("parallel -j1 <<EOF")
2000  format("parallel -j1 --sshloginfile $PBS_NODEFILE --wd $PWD <<EOF")
!2000  format("parallel -j",i5," --sshloginfile $PBS_NODEFILE --wd $PWD <<EOF") ! buggy
! running multiple instances of HEX on the same node is unstable
! I guess they're from too many attempts to access the same file
! Anyway, I'm sick of the errors.


END subroutine Dock_GeneratePopulation
