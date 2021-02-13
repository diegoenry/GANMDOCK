export FLAGS="-ffree-form -fno-backslash -ffree-line-length-none -fbackslash"
rm -rf teste.f90

for target in \
  Vectors.f90 \
  Dock_Vectors.f90 \
  CheckEtable.f90 \
  Dock_ChoosePair.f90 \
  Dock_FixPDB.f90 \
  Dock_GeneratePopulation.f90 \
  GADock.f90 \
  GeneratePopulation.f90 \
  GenerateStructures.f90 \
  InitAmplitudes.f90 \
  Roulette.f90 \
  SetDefaults.f90 \
  SetRank.f90 \
  SortFitness.f90 \
  StoreEtable.f90 \
  Zrank_Get_RESULTS.f90 \
  crossover.f90 \
  run_hex.f90 \
  Error.f90 \
  ReadCommandLine.f90 \
  RouletteElitism.f90 \
  WriteAmplitudes.f90 \
  Help.f90 \
  Showtime.f90 \
  main.f90 ; do
echo "Compiling ${target}"
gfortran -O2 -c ${target} $FLAGS
#gfortran-mp-4.4 -O2 -c ${target} $FLAGS

cat ${target} >> teste.f90

done

#link object
gfortran  -O2 *.o -o ../bin/ganmdock
#gfortran-mp-4.4 -O2 *.o -o ../bin/ganmdock

#FixPDB
gfortran -O2 fixpdb.f90 -o ../bin/fixpdb

#clean up
rm -rf *.o

