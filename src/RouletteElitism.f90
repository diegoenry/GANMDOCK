Subroutine RouletteElitism(nPersonDock) ! Competicao
use dock_vectors
implicit none
integer  :: i,j,k
real     :: number_1
real     :: FitnessSumElitism=0.0
real     :: tmp_sum
integer,intent(in)   :: nPersonDock

!computa a soma do fitness para o resto da populacao.
do i=(int(nPersonDock*ElitismRate)+1),nPersonDock*2
  FitnessSumElitism = pop_tmp(i)%fitness + FitnessSumElitism
enddo

!faz para todos elementos do "resto"
do j=(int(nPersonDock*ElitismRate)+1),nPersonDock
number_1 = ( rand(0)*FitnessSumElitism ) + 1
tmp_sum=0.0

!desta forma podem haver duplicatas, vamos ver se vai bem se não mudamos.
  do i=(int(nPersonDock*ElitismRate)+1), nPersonDock*2
    tmp_sum = tmp_sum + pop_tmp(i)%fitness
    
    if(tmp_sum > number_1) then

      pop(j)%receptor  = pop_tmp(i)%receptor
      pop(j)%ligand    = pop_tmp(i)%ligand
      pop(j)%fitness   = pop_tmp(i)%fitness
    
      do k=1,10
        pop(j)%energy(j) = pop_tmp(i)%energy(k)
        pop(j)%rms(j)    = pop_tmp(i)%rms(k)
      enddo
      
    endif

  enddo

enddo

end subroutine RouletteElitism
