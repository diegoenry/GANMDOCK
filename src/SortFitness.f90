!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine SortFitness(generation)
use vectors, only : nPersonDock
use dock_vectors
type(valor) ::  pop_swap

if (generation==0) then
  do i=1,nPersonDock-1
     do j=1,nPersonDock-1
          if( (pop(j)%fitness) < (pop(j+1)%fitness) ) then
              pop_swap    =   pop(j)      ! transfers the whole object
              pop(j)      =   pop(j+1)    ! transfers the whole object
              pop(j+1)    =   pop_swap    ! transfers the whole object
          endif
    enddo
 enddo
endif

if (generation>0) then
 do i=1,nPersonDock*2-1
    do j=1,nPersonDock*2-1
          if( (pop_tmp(j)%fitness) < (pop_tmp(j+1)%fitness) ) then
              pop_swap      =   pop_tmp(j)      ! transfers the whole object
              pop_tmp(j)    =   pop_tmp(j+1)    ! transfers the whole object
              pop_tmp(j+1)  =   pop_swap        ! transfers the whole object
          endif
    enddo
 enddo

endif

end subroutine SortFitness
