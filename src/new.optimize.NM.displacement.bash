#!/bin/bash
###########################################################
# programa capes/cofecub 2009-2011
# 
# diego e. barreto gomes 1,2  | diego@biof.ufrj.br
# luis p. b. scott       3    |
# pedro g. pascutti      1    | pascutti@biof.ufrj.br
# paulo m. bisch         1    | pmbisch@biof.ufrj.br
# david perahia          2    | david.perahia@ens-cachan.fr
# 
###########################################################

# this is a script to generate a charmm input file optimizing one
# normal mode and one associated amplitude and generates one pdb file
# containing a minimized structure.
#
# optimization of the structure starts from the crystal structure and 
# follows a protocol of molecular dynamics simulation at 30 kelvin, 
# then a final steepest descent minimization.
#
# fri dec  9 11:30:28 cet 2011

# Qui Jun 21 17:25:59 BRT 2012
# Modifiquei o script para nao fazer a 3a MD ja que o Maucosta mostrou que a estrutura ficaria melhor.




echo "* # this is a script to generate a charmm input file optimizing one
* # normal mode and one associated amplitude and generates one pdb file
* # containing a minimized structure.
* #
* # optimization of the structure starts from the crystal structure and 
* # follows a protocol of molecular dynamics simulation at 30 kelvin, 
* # then a final steepest descent minimization.
* #

bomlev -2

! faster routines are used
faster on

!==========
!directory name where the coordinates and trajectory will be written
!set dirnam output
!set file @filename

!default parameters
set temperature 30.0
set kmod 10000
!set nummod !!  will come from disp.dat
!set modnu  !!  will come from command line
!set disp   !!  not used
!set istr   !!  will come from command line  
!normal mode file name
set modnam @filename.mod

open read file unit 40 name @modnam

! initial structure
set minic @filename.pdb

!read topology and parameter files ! protein topology and parameter
open read card unit 10 name $CHARMMDATA/top_all22_prot.rtf
read  rtf card unit 10

open read card unit 20 name $CHARMMDATA/par_all22_prot.prm
read para card unit 20 flex

!read the psf
open unit 4 read card name @filename.psf
read psf card unit 4
close unit 4

! read the initial coordinates (energy minimized or x-ray one)
! for which the modes were computed
open read card unit 12 name @minic
read coor pdb resi unit 12
close unit 12
coor copy comp

! define energy model
! parameters from paulo ricardo batista 
! define force field 2
set par1 8.0
set par2 10.0
set par3 12.0
set par4 2.0
updat  inbfrq -1 -
atom rdiel switch vswitch -
ctonnb @par1 ctofnb @par2 cutnb @par3 eps @par4 e14fac 1.0 wmin 1.5

open read card unit 12 name @minic
read coor pdb resi unit 12
close unit 12
coor copy comp

"
nummod=0
while read line
do
   let nummod=$nummod+1
   q[$nummod]=$line
done<disp.dat

echo "
! set number of modes that will be combined
set nummod $nummod

!open a trajectory file for the normal mode coordinates
open write card unit 17 - 
   name @dirnam/@istr.qcr
"

for ((i=1; i<=nummod ;i++)) do
   echo "set q${i}  ${q[$i]}"
done

echo "
!------------
!minimization
!------------
! set up mmfp (umbrella) potential
! mxmd = maximum number of modes (read only once)
! kmdn = force constant for all the modes
! qn  =constrained coordinate value for mode imdn
! umdn = unit number for the modes
! imdn = mode number
! uout = unit for storing the normal coordinates
! krot = force constant for the overall rotation
! kcgr=force constant for the overall translation
! nsvq = normal coordinates saving frequency

mmfp
vmod init mxmd @nummod krot 0.000001 kcgr 1000 umdn 40 uout 17 nsvq 100 
"
for ((i=1; i<=nummod ;i++)) do
   let j=$i+6 #from mode 7
   echo "vmod add  qn @q${i} imdn ${j} kmdn 1000"
done

echo "
end

!quick minimization 
mini sd nstep 1000 tolgrd 0.01


#goto fim 


!##### md1 just relax a little #####
open write file unit 14 name @dirnam/@istr.md1.dcd
open write card unit 17 name @dirnam/@istr.md1.qcr
open write card unit 19 name @dirnam/@istr.md1.rst

dyna start nstep 1000  timestep 0.001  ntrfrq 10  -
    iseed $RANDOM -
    iprfrq 1000  ihtfrq 0 ieqfrq 500 imgfrq -1 inbfrq -1 -
    iunwri 19 iunrea -16  iuncrd 14 -
    nprint 100 nsavc 1000 echeck 100000 -
    ihbfrq 0 isvfrq 1000 -
    iasors 0 iasvel 1 iscvel 0 ichecw 0 -
    firstt 30.0 finalt @temperature teminc 0.0 tstruc @temperature -
    tconstant tcouple 0.1 treference @temperature

set i 2000

!conferir close unit 14
!conferir close unit 17
close unit 19

!##### md2 increasing kmod #####

label kmdn_loop

open read  card unit 16 name @dirnam/@istr.md1.rst
open write file unit 14 name @dirnam/@istr.md2.dcd
open write card unit 17 name @dirnam/@istr.md2.qcr
open write card unit 19 name @dirnam/@istr.md2.rst

mmfp
"
for ((i=1; i<=nummod ;i++)) do
  let j=$i+6 #from mode 7
  echo "vmod change nrestraint  ${i}  qn @q${i}  kmdn @i"
done
  
echo "
end

dyna cpt restart nstep  2000  timestep 0.001  ntrfrq 10  -
    iprfrq 1000  ihtfrq 0 ieqfrq 500 imgfrq -1 inbfrq -1 -
    iunwri 19 iunrea 16  iuncrd 14 -
    nprint 100 nsavc 1000 echeck 100000 -
    ihbfrq 0 isvfrq 1000 -
    iasors 0 iasvel 0 iscvel 0 ichecw 0 -
    firstt @temperature finalt @temperature teminc 0.0 tstruc @temperature -
    tconstant tcouple 0.1 treference @temperature

increment i by 1000

!we need to restart from a restart file :(
envi md1 @dirnam/@istr.md1.rst
envi md2 @dirnam/@istr.md2.rst
system \"MD1=\`echo \$MD1 |tr [A-Z] [a-z]\` ; MD2=\`echo \$MD2 |tr [A-Z] [a-z]\` ; cp \$MD2 \$MD1\"

if i .lt. 10000 goto kmdn_loop

close unit 14
close unit 16
close unit 17
close unit 19

!##### MINIMIZATION #####
!# here we increase (put 10000) the kmod back to the original value and run the simulation
set i @kmod
mmfp
"
for ((i=1; i<=nummod ;i++)) do
   let j=$i+6 #from mode 7
   echo "vmod change nrestraint  ${i}  qn @q${i}  kmdn @i"
done

echo "
end

open read  card unit 16 name @dirnam/@istr.md2.rst
open write file unit 14 name @dirnam/@istr.mini.dcd
open write card unit 17 name @dirnam/@istr.mini.qcr
open write card unit 19 name @dirnam/@istr.mini.rst

!wrap up with a final minimization
mini conj nstep 1000 tolgrd 0.01


close unit 14
close unit 16
close unit 17
close unit 19

label fim

open write card unit 19 name -
   @dirnam/@istr-final.out

outu 19

! compute the energy of the final structure
ener

! compute the rmsd from the initial structure
coor rms

mmfp
vmod print
end

outu 6

! write the final minimized coordinates
open write card unit 12 name - 
   @dirnam/@istr.pdb
write coor pdb unit 12
* system name
* minimized coordinates
* total_e= ?ener rmsd= ?rms
*

time diff 
! close the normal mode file
close unit 40

close unit 12
close unit 17
close unit 19

outu 6
stop
"
