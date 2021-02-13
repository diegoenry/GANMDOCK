echo a | isscore.pl -w scratch ./dockout/${2}_${3}/${1} RL complex.pdb RL > isscore/${1}.out 

isscore_rms=`awk '/RMSD / {print $6}' isscore.out`
isscore_irms=`awk '/RMSD=/ {print $9}' isscore.out`
isscore=`awk '/IS-score    =/ {print $3}' isscore.out`

echo "IS-Score |  RMS  |  iRMS  "
echo "$isscore  $isscore_rms $isscore_irms" > isscore.out 
