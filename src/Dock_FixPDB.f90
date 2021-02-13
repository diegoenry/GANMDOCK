!Use convpdb.pl to remove hydrogens and fix the pdb for docking
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Dock_FixPDB(nPersonReceptor,nPersonLigand)
character(len=5)    :: string
character(len=512)  :: command

do i=1,nPersonReceptor

  write(*,"(2(A,i5,2x),'\r',$)") "fixing pdb files: Receptor:",i
  write(unit=string, fmt='(I5)') i

  command="convpdb.pl -setchain R -out charmm22noh output/r_"//trim(adjustl(string))//".pdb > dockinput/fixed.r_"//trim(adjustl(string))//".pdb"
  call system(command)

enddo

do i=1,nPersonLigand

  write(*,"(2(A,i5,2x),'\r',$)") "fixing pdb files: Ligand:",i
  write(unit=string, fmt='(I5)') i

  command="convpdb.pl -setchain L -out charmm22noh output/l_"//trim(adjustl(string))//".pdb > dockinput/fixed.l_"//trim(adjustl(string))//".pdb"
  call system(command)

enddo

end subroutine