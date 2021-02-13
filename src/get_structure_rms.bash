nPerson=40

for ((i=1;i<=nPerson;i++)) do
rms[$i]=`rms.pl -out ca -fit -fitsel ca -align fasta ligand.pdb output/l_${i}.pdb`
done


rm -rf rms.ligand.dat
for ((i=1;i<=nPerson;i++)) do
   echo $i ${rms[$i]} >> rms.ligand.dat
done

for ((i=1;i<=nPerson;i++)) do
   rms[$i]=`rms.pl -out ca -fit -fitsel ca -align fasta receptor.pdb output/r_${i}.pdb`
done


rm -rf rms.receptor.dat
for ((i=1;i<=nPerson;i++)) do
      echo $i ${rms[$i]} >> rms.receptor.dat
   done


