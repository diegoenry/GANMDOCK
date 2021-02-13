nPerson=40
rm -rf disp.ligand.dat
rm -rf disp.receptor.dat

for ((i=1;i<=nPerson;i++)) do

   disp1=`awk '/RESTRAINT   1 MOD/ {print $5}' output/l_$i-final.out`
   disp2=`awk '/RESTRAINT   2 MOD/ {print $5}' output/l_$i-final.out`
   disp3=`awk '/RESTRAINT   3 MOD/ {print $5}' output/l_$i-final.out`
   disp4=`awk '/RESTRAINT   4 MOD/ {print $5}' output/l_$i-final.out`
   disp5=`awk '/RESTRAINT   5 MOD/ {print $5}' output/l_$i-final.out`

   echo $i $disp1 $disp2 $disp3 $disp4 $disp5 >> disp.ligand.dat

done

for ((i=1;i<=nPerson;i++)) do
      
      disp1=`awk '/RESTRAINT   1 MOD/ {print $5}' output/r_$i-final.out`
      disp2=`awk '/RESTRAINT   2 MOD/ {print $5}' output/r_$i-final.out`
      disp3=`awk '/RESTRAINT   3 MOD/ {print $5}' output/r_$i-final.out`
      disp4=`awk '/RESTRAINT   4 MOD/ {print $5}' output/r_$i-final.out`
      disp5=`awk '/RESTRAINT   5 MOD/ {print $5}' output/r_$i-final.out`

      echo $i $disp1 $disp2 $disp3 $disp4 $disp5 >> disp.receptor.dat
  
done

