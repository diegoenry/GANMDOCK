convpdb.pl -setchain R -out charmm22noh ~/PROGRAMAS/GANMDOCK.1.0.RELEASE/teste/1AY7/1AY7_r_b.pdb > rb.pdb 
convpdb.pl -setchain L -out charmm22noh ~/PROGRAMAS/GANMDOCK.1.0.RELEASE/teste/1AY7/1AY7_l_b.pdb > lb.pdb
echo TER > tmp
cat rb.pdb TER lb.pdb > complex.pdb

#ganmdock -nodock -nmrec nmrec -nmlig nmlig -nproc 2
ganmdock -onlydock -nmrec nmrec -nmlig nmlig -nproc 2

ganmdock -nmrec nmrec -nmlig nmlig -nproc 2 -onlydock -etable etable.dat

ganmdock -onlydock -nmrec nmrec -nmlig nmlig -nproc 64 -pbs 1 -etable etable.dat -psize 24 -pdock 64
