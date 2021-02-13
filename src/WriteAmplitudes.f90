subroutine WriteAmplitudes
use vectors

write(*,*) "Avaible values for each mode for: Receptor "
do i=1,nmodes
  write(*,*) "MODE ",i+6, NormalModeReceptor(i)%list(:)
enddo

write(*,*)
write(*,*) "Chosed values for each Receptor"
write(*,*) "Receptor  | ","  MODE 7 | ","  MODE 8 |","  MODE 9"
do i=1,nPerson
  write(*,*) i, Receptor(i)%Disp(:)
enddo

write(*,*)
write(*,*) "Avaible values for each mode for: Ligand "
do i=1,nmodes
  write(*,*) "MODE ",i+6, NormalModeLigand(i)%list(:)
enddo

write(*,*)
write(*,*) "Chosed values for each Ligand"
write(*,*) "Ligand  | ","  MODE 7  | ","  MODE 8  |","  MODE 9"
do i=1,nPerson
  write(*,*) i, Ligand(i)%Disp(:)
enddo

end subroutine WriteAmplitudes
