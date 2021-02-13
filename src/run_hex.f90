SUBROUTINE run_hex(i,generation,receptor,ligand)
use vectors , only   :  nPersonDock,pbs,nproc,nnodes
use dock_vectors, only  : exclude
integer,intent(in)   :: i
integer,intent(in)   :: generation
integer,intent(in)   :: receptor
integer,intent(in)   :: ligand

!if(generation==0) then
  if (i==1) then
    open(99,file="runhex")
    if (pbs==0) write(99,1000)
    if (pbs==1) write(99,2000)
    if (pbs==2) write(99,3000)
  endif

   if(exclude(i)==0)  write(99,"(A,1x,3(i5))") "./dock_and_rerank.bash ", i, receptor, ligand

  if(i==nPersonDock) then
!    if(pbs/=0) write(99,"(A3)") "EOF"
    write(99,"(A3)") "EOF"
    close(99)
  endif
!endif

1000  format("parallel -j1 <<EOF")
2000  format("parallel -j1 --sshloginfile $PBS_NODEFILE --wd $PWD <<EOF")
3000  format("parallel -j1 --sshloginfile $PBS_NODEFILE --wd $PWD <<EOF")
end subroutine run_hex

subroutine WriteHex
write(*,*)
end subroutine WriteHex