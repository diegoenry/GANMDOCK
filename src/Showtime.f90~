!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine credits

write(*,*) "###############################################################"
write(*,*) "#                   GANMDOCK - version 1.0rc3                 #"
write(*,*) "#               Mon Nov 19 16:46:08 BRST 2012                 #"
write(*,*) "#                                                             #"
write(*,*) "# Diego E B Gomes(1,2), Pedro G Pascutti(2),Paulo M Bisch(2), #"
write(*,*) "# David Perahia(3), Luis P B Scott(1)  | diego@biof.ufrj.br   #"
write(*,*) "#                                                             #"
write(*,*) "# 1 - Universidade Federal do ABC                             #"
write(*,*) "# 2 - Univerisdade Federal do Rio de Janeiro                  #"
write(*,*) "# 3 - Ecole Normale Superieure de Cachan                      #"
write(*,*) "###############################################################"

end


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine ShowDisplacements
use vectors
  write(*,*) "Receptor"
  do i=1,nModesReceptor
    write(*,*) "mode:",NormalModeReceptor(i)%ind
    write(*,*) "disp:",NormalModeReceptor(i)%list
    write(*,*)
    enddo
  write(*,*) "Ligand"
  do i=1,nModesLigand
    write(*,*) "mode:",NormalModeLigand(i)%ind
    write(*,*)"disp:",NormalModeLigand(i)%list
    write(*,*)
  enddo
end subroutine ShowDisplacements



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine ShowPopulation
use vectors

  write(*,*) "Receptors"
  
  do i=1,nPersonReceptor
    write(*,*) "receptor:",i,(Receptor(i)%Disp(j),j=1,nModesReceptor)
  enddo

  write(*,*) "Ligands"
  do i=1,nPersonLigand
    write(*,*) "ligand:",i,(Ligand(i)%Disp(j),j=1,nModesLigand)
  enddo

end subroutine ShowPopulation


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine ShowDockingPopulation
use vectors
use dock_vectors
  do i=1,nPersonDock
    write(*,*) i,pop(i)%receptor,pop(i)%ligand
  enddo
end subroutine ShowDockingPopulation



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine ShowEnergyTable
use vectors
use dock_vectors
write(*,*) "Etable has ", EnergyTableLast, "elements"
  do j=1,EnergyTableLast
    write(*,*)  i,&
                EnergyTable(j)%receptor,&
                EnergyTable(j)%ligand,&
                EnergyTable(j)%energy(1),&
                EnergyTable(j)%fitness
  enddo
end subroutine ShowEnergyTable



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine WriteLog(g)
use dock_vectors
use vectors
character(len=128)   :: filename
character(len=128)   :: filename2
character(len=128)   :: filename3
character(len=128)   :: command
integer,intent(in)  :: g
integer  :: p !population size ( 0 or >0)
character(len=5) :: rec
character(len=5) :: lig
real  :: rms(1000)
real  :: energy(1000)
real  :: zrank(1000)
character(len=128)  :: str(1000)
real  :: isscore(1000)
real  :: isscoreRMS(1000)
real  :: isscoreIRMS(1000)

if (g==0) p=nPersonDock
if (g/=0)  p=nPersonDock*2


write(66,*) "# Generation ",g

DO j=1,p

write(66,*) "#",pop_tmp(j)%receptor,pop_tmp(j)%ligand
write(unit=rec, fmt='(I5)') pop_tmp(j)%receptor
write(unit=lig, fmt='(I5)') pop_tmp(j)%ligand

filename="dockout/"//trim(adjustl(rec))//"_"//trim(adjustl(lig))//".hex"
filename2="dockout/"//trim(adjustl(rec))//"_"//trim(adjustl(lig))//".list.zr.out"

open(1,file=trim(filename))

do 
read(1,*) str(1)
  if (str(1)=="1") exit
enddo
backspace(1)

n=0
do i=1,100
  read(1,*,end=10,err=10) ind,nlixo,rms(i),energy(i)
  n=n+1
enddo
10 continue
close(1)

open(1,file=trim(filename2))
do i=1,n
  read(1,*) str(i),zrank(i)
enddo
close(1)

goto 99 ! This is not working

filename3="dockout/"//trim(adjustl(rec))//"_"//trim(adjustl(lig))//".isscore.rms"
open(1,file=trim(filename3))
do i=1,n
  read(1,*) isscoreRMS(i)
enddo
close(1)

filename3="dockout/"//trim(adjustl(rec))//"_"//trim(adjustl(lig))//".isscore.irms"
open(1,file=trim(filename3))
do i=1,n
  read(1,*) isscoreIRMS(i)
enddo
close(1)

!filename3="dockout/"//trim(adjustl(rec))//"_"//trim(adjustl(lig))//".isscore"
!open(1,file=trim(filename3))
!do i=1,n
!  read(1,*) isscore(i)
!enddo
!close(1)

99 continue 
!str = docking result
!
do i=1,n
!  write(66,*) str(i),energy(i),zrank(i),isscore(i),rms(i),isscoreRMS(i),isscoreIRMS(i)
  write(66,*) trim(adjustl(str(i))),energy(i),zrank(i),rms(i)!,isscoreRMS(i),isscoreIRMS(i)
enddo

ENDDO

end subroutine WriteLog
