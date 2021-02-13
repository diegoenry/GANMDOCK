export root=$PWD
cp complex.pdb /dev/shm/

#go to dockout folder
cd ${root}/dockout/

#list available docks
list=`ls -d */ | sort -n`

for dir in $list ; do
  cd ${root}/dockout/$dir
  mkdir -p /dev/shm/$dir/isscore
  echo $PWD
  ls *.pdb | parallel "echo n | isscore.pl -w /dev/shm/scratch {} RL /dev/shm/complex.pdb RL > /dev/shm/$dir/isscore/{}.out"
  cat /dev/shm/$dir/isscore/*.out > all.isscore.out
  trap "echo Exited!; exit;" SIGINT SIGTERM
done
