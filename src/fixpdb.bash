#!/bin/bash
export root=$PWD
nperson=$1

cd ${root}/output/
for ((pdb=1; pdb<=${nperson}; pdb++)) do
   echo -ne "Fixing receptor $pdb \r"
   convpdb.pl  -setchain R -out generic r_${pdb}.pdb > fixed.r_${pdb}.pdb
done

cd ${root}/output/
for ((pdb=1; pdb<=${nperson}; pdb++)) do
    echo -ne "Fixing ligand $pdb \r"
    convpdb.pl  -setchain L -out generic l_${pdb}.pdb > fixed.l_${pdb}.pdb
done


