program SuperRunISscore
character(len=1000) :: str
character(len=128)  :: filename
character(len=128)  :: command
character :: lixo
real  :: hex,zrank,rms
real  :: isscore,isscore_rms,isscore_irms
character(len=5) c_rec,c_lig

open(1,file="docklogfile.dat")
open(3,file="finalresult.dat")
10 continue

read(1,"(a1000)",end=99,err=99) str
if(index(str,"Gen")/=0) then
  write(3,*) trim(str)
  read(1,"(a1000)") str
endif

if(index(str,"#  ")/=0) then
  read(str,*) lixo,nrec,nlig
  write(3,*) "#  ", nrec,nlig
  write(*,"(7(A12))") "# File","Hex","Zrank","RMS","IS-Score","IS-RMS","IS-iRMS"
endif

!int to char
write(c_rec,"(i5)") nrec
write(c_lig,"(i5)") nlig

!Write finalresult header.
write(3,*) "@ filename | hex | zrank | rms | isscore | isscore_rms | isscore_irms "
do i=1,100
  read(1,*) filename,hex,zrank,rms
  command="bash runisscore.bash "//trim(filename)//" "//c_rec//" "//c_lig
  call system(command)
  
  open(2,file="isscore.out")
  read(2,*) isscore, isscore_rms, isscore_irms
  close(2)
  write(*,"(A20,6f12.5)") trim(filename),hex,zrank,rms,isscore,isscore_rms,isscore_irms
  write(3,"(A20,6f12.5)") trim(filename),hex,zrank,rms,isscore,isscore_rms,isscore_irms
enddo

goto 10

99 continue
! # Generation            0
! #          14          13
! 14_13_0001.pdb  -390.37000       26.055401       21.750000
end
