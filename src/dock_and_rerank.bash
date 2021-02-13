#!/bin/bash
# Written by
# Diego Enry B. Gomes (1,2)  | diego@biof.ufrj.br
# Elodie Japaud (2,3)        |
# Luis Paulo Scott (4)       |
# Pedro Pascutti (1)         |
# Paulo Mascarello Bisch (1) |
# David Perahia (2)          |
#
# 1) Universidade Federal do Rio de Janeiro, Rio de Janeiro - Brazil
# 2) Universite Paris Diderot, Paris - France
# 3) Ecole Normale Superieure de Cachan, Cachan - France
# 4) Universidade Federal do ABC - Santo André - Brasil


echo "Running on $HOSTNAME"
sleep 2

#### HEX ####
export HEX_ROOT=/home/dgomes/software64/hex/ # required
export HEX_PDB=.
export HEX_MACROS=.
export HEX_COLOUR=.
export HEX_DATA=$HEX_ROOT/data/
export HEX_CACHE=$HEX_ROOT/hex_cache
#export HEX_STARTUP=/var/lmdm/dgomes/hex/examples/startup_v5.mac
export HEX_CPUS=8
export HEX_GPUS=1
export HEX_MESSAGES=yes
#export HEX_STEREO=yes
export PATH=$PATH:$HEX_ROOT/bin    # required
#### HEX ####i
export PATH=${PATH}:"/home/dgomes/software64/parallel/bin/"
export PATH=${PATH}:"/home/dgomes/GANMDOCK.1.0.RELEASE/bin/"


export root=$PWD
person=$1
receptor=$2
ligand=$3


write_hex_input() {
echo "
close_all
open_receptor ./dockinput/fixed.r_${receptor}.pdb
open_ligand   ./dockinput/fixed.l_${ligand}.pdb
open_complex  ./complex.pdb

DOCKING_RECEPTOR_SAMPLES  492
DOCKING_LIGAND_SAMPLES    492
DOCKING_ALPHA_SAMPLES     128

RECEPTOR_RANGE_ANGLE       30
LIGAND_RANGE_ANGLE         30
TWIST_RANGE_ANGLE          90


R12_RANGE                  21
R12_STEP                 0.75
GRID_SIZE                0.6

DOCKING_MAIN_SCAN          16
DOCKING_MAIN_SEARCH        25
max_docking_solutions      1000

#pede a GPU
docking_fft_device 1

# Actually do the docking
activate_docking

# Save best results
save_range 1 100   /dev/shm/ ${receptor}_${ligand}_ pdb
save_summary      ./dockout/${receptor}_${ligand}.sum
save_transform    ./dockout/${receptor}_${ligand}.hex

EXIT
" > /dev/shm/hex.${person}.macro
}

exec_hex() {
echo "running ${receptor} vs ${ligand}"
$HEX_ROOT/bin/hex -noexec -nogui -e /dev/shm/hex.${person}.macro -l /dev/shm/${receptor}"_"${ligand}".log" >/dev/null
}

fix_pdb() {
# HEX is stupid, it merges the PDB files but discarts the TER and END records
echo "Fixing PDB $person"
#cd ${root}/dockout/
cd /dev/shm/
for pdbfile in ${receptor}"_"${ligand}"_"*.pdb ; do
  grep ATOM $pdbfile > tmp.$pdbfile
  fixpdb tmp.$pdbfile > $pdbfile
  rm -rf tmp.$pdbfile
done
}

re_score() {
echo "ReScoring $person"
#cd ${root}/dockout/
cd /dev/shm/
# run ZRANK to re-rank the results
list=${receptor}"_"${ligand}".list"
ls ${receptor}"_"${ligand}*pdb > $list
zrank $list
}

re_rank() {
echo "ReRanking $person"
#cd ${root}/dockout/
cd /dev/shm/
sort -gk 2 ${list}.zr.out > $list.zr.reranked
head -1 $list.zr.reranked > ${list}.best.energy
}

move_files() {
cd ${root}/dockout/
rm -rf ${receptor}"_"${ligand}/
cd /dev/shm/
mkdir ${receptor}"_"${ligand}
mv ${receptor}"_"${ligand}_*pdb ${receptor}"_"${ligand}/
tar cfj ${root}/dockout/${receptor}"_"${ligand}.tar.bz2 ${receptor}"_"${ligand}/
rm -rf ${receptor}"_"${ligand}/
mv /dev/shm/${list}                        ${root}/dockout/
mv /dev/shm/${receptor}"_"${ligand}".log"  ${root}/dockout/
mv /dev/shm/${list}.zr.out                 ${root}/dockout/
mv /dev/shm/${list}.zr.reranked            ${root}/dockout/
mv /dev/shm/${list}.best.energy            ${root}/dockout/
}

write_hex_input ; sleep 1
exec_hex        ; sleep 1
fix_pdb         ; sleep 1
re_score        ; sleep 1
re_rank         ; sleep 1
move_files      ; 

