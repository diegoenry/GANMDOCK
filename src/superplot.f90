program SuperPlot
character(len=1000) :: str
character(len=128) :: filename
character(len=128) :: command
character :: lixo
real  :: hex,zrank,rms
real  :: isscore, isscore_rms, isscore_irms
integer :: gen
character(len=5) :: genstr
character(len=10) :: outputfile

open(1,file="finalresult.teste")

gen=-1

10 continue
read(1,"(a1000)",end=99,err=99) str
if(index(str,"Gen")/=0) then
gen=gen+1
!  write(2,*) trim(str)
  read(1,"(a1000)") str
  write(unit=genstr, fmt='(I5)') gen
  write(unit=outputfile, fmt='(A10)') "gen."//trim(adjustl(genstr))
  open(2,file=trim(adjustl(outputfile)))
endif

if(index(str,"#  ")/=0) then
  read(str,*) lixo,nrec,nlig
endif

do i=1,100
  read(1,*) filename,hex,zrank,rms,isscore, isscore_rms,isscore_irms
  write(2,*) hex,zrank,rms,isscore, isscore_rms,isscore_irms
enddo

goto 10
close(2)

99 continue
end
