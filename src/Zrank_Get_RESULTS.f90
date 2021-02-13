subroutine Zrank_Get_RESULTS(generation,nPersonDock)
use dock_vectors
integer,intent(in)   :: nPersonDock
integer,intent(in)   :: generation
character(len=5)     :: receptor
character(len=5)     :: ligand
character(len=40)    :: trash
character(len=256)   :: arquivo

if (generation==0) then
do i=1,nPersonDock
  if (exclude(i)==0) then
    write(receptor,"(i5)") (pop(i)%receptor)
    write(ligand,"(i5)")   (pop(i)%ligand)
    arquivo="dockout/"//trim(adjustl(receptor))//"_"//trim(adjustl(ligand))//".list.zr.reranked"
    open(99,file=trim(arquivo))
    pop(i)%energy(:)=0.0 !acochambracao p caso ter menos de 10
    do j=1,10 !10best
      read(99,*,end=11,err=11) trash, pop(i)%energy(j)
    enddo ; goto 11
11  continue
    close(99)
  endif
enddo

else

do i=nPersonDock+1,nPersonDock*2
  if(exclude(i)==0) then
    write(receptor,"(i5)") (pop_tmp(i)%receptor)
    write(ligand,"(i5)")   (pop_tmp(i)%ligand)
    arquivo="dockout/"//trim(adjustl(receptor))//"_"//trim(adjustl(ligand))//".list.zr.reranked"
    open(99,file=trim(arquivo))
    pop_tmp(i)%energy(:)=0.0 !acochambracao p caso ter menos de 10
    do j=1,10 !10best
      read(99,*,end=13,err=13) trash, pop_tmp(i)%energy(j)
    enddo ; goto 13
13  continue
    close(99)
  endif
enddo
endif

end