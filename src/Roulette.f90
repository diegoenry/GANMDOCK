!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Roulette(crossover_1,crossover_2)
use vectors, only : nPersonDock
use dock_vectors
real        ::  number_1
real        ::  number_2
real        ::  tmp_sum
integer     ::  counter=0
integer,intent(out) :: crossover_1
integer,intent(out) :: crossover_2

write(*,"(A)") "Roulette"
 
  number_1  =   (rand(0)*Fitness_sum)+1
  tmp_sum=0.0d0

do i=1,nPersonDock
  crossover_1=i
  tmp_sum = tmp_sum + pop(i)%Fitness
  if( tmp_sum > number_1 ) exit
enddo

100   number_2  =   (rand(0)*Fitness_sum)+1.0
 tmp_sum=0.0d0

do i=1,nPersonDock
  crossover_2=i
  tmp_sum = tmp_sum + pop(i)%Fitness
  if( tmp_sum > number_2 ) exit
enddo
  if( crossover_1 == crossover_2 ) then
  counter = counter + 1
  if (counter > nPersonDock) then
    write(*,*) "[Roulette routine] Trapped in an infinite loop, will accepting any crossover."
    write(*,*) "[Roulette routine] Don’t worry, population will converge anyway."
    goto 200
  endif
  goto 100
  endif

200 continue
END subroutine Roulette