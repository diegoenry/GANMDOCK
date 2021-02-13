subroutine SetDefaults
use vectors

pbs=0
nPerson=10
ngenerations=5
gridspacing=0.5
nproc=1
nnodes=1
nmodes=1

!docking
ElitismRate=0.1

!extra
nodock=0
onlydock=0
rmscutoff=0.5
nPersonDock=20

EnergyTableFile=""

end subroutine SetDefaults