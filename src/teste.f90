module vectors

!Genetic algorithm stuff
!integer  ::  nPerson  !> @var Number of individuals in the population 
!integer  ::  nmodes   !> @var Number of normal modes to optimize
!integer  ::  ngenerations !> @var Number of generations
!real     ::  MutationRate !> @var Mutation Rate
!real     ::  ElitismRate  	!> @var Elitism Rate
!integer  ::  crossover_1  	!> @var First individual to crossover
!integer  ::  crossover_2  	!> @var Second individual to crossover
!real     ::  Fitness_sum  	!> @var Sum of the fitness

!!!!!!!!!
integer  :: pbs         	!> @var check if we will run locally or using a cluster
integer  :: nproc       	!> @var number of processors per job
integer  :: nPerson     	!> @var Population size
integer  :: nmodes      	!> @var Number of normal modes to optimize
real     :: gridspacing 	!> @var grid spacing
integer  :: ngenerations 	!> @var number of generations
integer  :: nnodes      	!> @var useless
integer  :: modnum      	!> @var index to reference a particular mode
integer  :: nPersonReceptor
integer  :: nPersonLigand
integer  :: nPersonDock !> Size of the Docking population


type type_NormalMode
  integer           :: ind      !> index
  real              :: maxdisp 	!> Maximum displacement
  real              :: mindisp  !> min
  integer           :: ndisp    !> number of available displacements (related to gridspacing)
  real,allocatable  :: list(:)  !> list   of available displacements
  integer			:: nmodes	!> number of normal modes
end type type_NormalMode

type(type_NormalMode),allocatable :: NormalMode(:)
type(type_NormalMode),allocatable :: NormalModeReceptor(:)
type(type_NormalMode),allocatable :: NormalModeLigand(:)


type chromossome
  integer            :: ind ! index
  real,allocatable   :: Disp(:)
  real               :: Energy
  real               :: Fitness
end type chromossome

type(chromossome),allocatable :: Receptor(:)
type(chromossome),allocatable :: Ligand(:)

!type(chromossome),allocatable :: Person(:)
!type(chromossome),allocatable :: pop_tmp(:)
!type(chromossome),allocatable :: BestEnergy(:) !19Dez2011
!type(chromossome)             :: crossover(2)
!type(chromossome)             :: swap


!Keep all the sampled displacements in a list
real, allocatable :: ReceptorTable(:,:) !nmodes,ndisp
real, allocatable :: LigandTable(:,:)   !nmodes,ndisp

integer :: nodock, onlydock
real    :: rmscutoff


!Energy Table Stuff
integer     ::  nEnergyTable      !> @var Number of values in the Energy Table
integer     ::  EnergyTableLast   !> @var last value inserted in the energy table (beta)
integer     ::  nBest             !> @var number of best structures to write out
character   ::  EnergyTableFile*128


!Files with the vectors to optimize
character(len=128)  ::  NormalModeReceptorFile !> @var File containing the NM vectors and min-max amplitudes for the receptor
character(len=128)  ::  NormalModeLigandFile   !> @var File containing the NM vectors and min-max amplitudes for the ligand
integer             ::  nModesReceptor !unused so far
integer             ::  nModesLigand   !unused so far
end module vectors


module colors
implicit none
character       ::  red*5
character       ::  green*5
character       ::  blue*5
character       ::  normal*5
character       ::  underline*4
contains
    subroutine setcolors
        red=char(27)//"[31m"
        green=char(27)//"[32m"
        blue=char(27)//"[34m"
        normal=char(27)//"[0m"
        underline=char(27)//"[4m"
    end subroutine
end module colors!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module dock_vectors
type    ::  valor
  integer     ::  receptor
  integer     ::  ligand
  real        ::  energy(10)
  real        ::  fitness
  real        ::  rms(10)
!  real        ::  fraction(10) !> Fraction of known contacts
!  look for residues either in receptor or ligand (or a pair) MUST be in the contact region
end type valor

type(valor),allocatable ::  pop(:)
type(valor),allocatable ::  pop_tmp(:)
type(valor),allocatable ::  EnergyTable(:)

integer  :: crossover_method
real     :: Fitness_sum
real     :: ElitismRate
integer,allocatable  :: exclude(:)

integer,allocatable  :: duplicate(:,:)

end module dock_vectors!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine CheckEnergyTable(generation)
use vectors, only : nPerson,EnergyTableLast,nPersonLigand,nPersonReceptor,nPersonDock
use dock_vectors
integer,intent(in) :: generation
integer  :: receptor,ligand

!reset excluded list
exclude(:)=0 

if (generation==0) then
  do i=1,nPersonDock
    receptor = pop(i)%receptor
    ligand   = pop(i)%ligand
    do j=1,EnergyTableLast
      if (EnergyTable(j)%receptor==receptor.and.EnergyTable(j)%ligand==ligand) then
        pop(i)=EnergyTable(j) ! transfers the whole object
        exclude(i)=1
        write(*,"(i5,i5,A)") receptor,ligand, " found in the Energy Table"
        exit
      endif
    enddo
  enddo

else

  do i=nPersonDock+1,nPersonDock*2
    receptor = pop_tmp(i)%receptor
    ligand   = pop_tmp(i)%ligand
    do j=1,EnergyTableLast
      if (EnergyTable(j)%receptor==receptor.and.EnergyTable(j)%ligand==ligand) then
        pop_tmp(i)=EnergyTable(j) ! transfers the whole object
        exclude(i)=1
        write(*,"(i5,i5,A)") receptor,ligand, " found in the Energy Table"
        exit
      endif
    enddo
  enddo

endif

end subroutine CheckEnergyTable
!######################################################################
SUBROUTINE Dock_ChoosePair(current,&
                           receptor,&
                           ligand,&
                           nPersonReceptor,&
                           nPersonLigand,&
                           nPersonDock)
use dock_vectors
integer,intent(in)   :: current        !current value in loop that calls choose pair.
!integer,intent(in)   :: nPerson
integer,intent(in)   :: nPersonReceptor
integer,intent(in)   :: nPersonLigand
integer,intent(out)  :: receptor
integer,intent(out)  :: ligand

10 continue
! select the structures to dock
receptor  =   int(rand(0)*(nPersonReceptor))+1
ligand    =   int(rand(0)*(nPersonLigand))+1

! This will garantee we're not repeating individuals
! However it will be an infinite loop if population size bigger than all avaible
! combinations ( n_receptor * n_ligand ).

ndupes=0

do i=1,current-1
  if (pop(i)%receptor==receptor.and.pop(i)%ligand==ligand) then
    write(*,*) "Found a duplicate",current 
    ndupes=ndupes+1
    if (ndupes>=nPersonDock) then 
      write(*,*) "ndupes=nPerson"
      STOP
    endif
    goto 10
  endif
enddo

END subroutine
!######################################################################
!Use convpdb.pl to remove hydrogens and fix the pdb for docking
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Dock_FixPDB(nPersonReceptor,nPersonLigand)
character(len=5)    :: string
character(len=512)  :: command

do i=1,nPersonReceptor

  write(*,"(2(A,i5,2x),'\r',$)") "fixing pdb files: Receptor:",i
  write(unit=string, fmt='(I5)') i

  command="convpdb.pl -setchain R -out charmm22noh output/r_"//trim(adjustl(string))//".pdb > dockinput/fixed.r_"//trim(adjustl(string))//".pdb"
  call system(command)

enddo

do i=1,nPersonLigand

  write(*,"(2(A,i5,2x),'\r',$)") "fixing pdb files: Ligand:",i
  write(unit=string, fmt='(I5)') i

  command="convpdb.pl -setchain L -out charmm22noh output/l_"//trim(adjustl(string))//".pdb > dockinput/fixed.l_"//trim(adjustl(string))//".pdb"
  call system(command)

enddo

end subroutine! Generates the initial population
!########################################################################
SUBROUTINE Dock_GeneratePopulation(generation)
use dock_vectors
use vectors
integer,intent(in)   :: generation

write(*,*) "Dock_GeneratePopulation",generation

open(99,file="runhex")
if (pbs==0) write(99,1000)
if (pbs==1) write(99,2000) !nproc

if(generation==0) then

  i=1 !bias to include the Right result.
  pop(i)%receptor=1
  pop(i)%ligand=1

  do i=2,nPersonDock
    call Dock_ChoosePair(i,pop(i)%receptor ,pop(i)%ligand   ,&
                           nPersonReceptor ,nPersonLigand   ,&
                           nPersonDock)
  enddo
  
  call CheckEnergyTable(generation)

  do i=1,nPersonDock
    if (exclude(i)==0) then
      write(*,*) pop(i)%receptor,pop(i)%ligand,exclude(i),"running"
      write(99,"(A,1x,3(i5))") "./dock_and_rerank.bash ", i, pop(i)%receptor,pop(i)%ligand
    else
      write(*,*) pop(i)%receptor,pop(i)%ligand,exclude(i),"Not running"
    endif
  enddo
  write(99,"(A3)") "EOF"
  close(99)

!Actually run the dockings
call system("cat runhex")
call system("bash runhex")

call Zrank_Get_RESULTS(generation,nPersonDock)
call SetRank(generation)
call StoreEnergytable(generation) !must come before SortFitness
call SortFitness(generation)

!transfer best results from pop_tmp to pop
  do i=1,nPersonDock
    pop_tmp(i)%receptor  = pop(i)%receptor
    pop_tmp(i)%ligand    = pop(i)%ligand
    pop_tmp(i)%fitness   = pop(i)%fitness
    do j=1,10
    pop_tmp(i)%energy(j) = pop(i)%energy(j)
    pop_tmp(i)%rms(j)    = pop(i)%rms(j)
    enddo
  enddo

endif


if (generation>0) then
  write(*,*) "Dock_GeneratePopulation",generation

!Scheme of crossover when there is Elitism
  call crossover
  call CheckEnergyTable(generation)
  call CheckDuplicate

  do i=nPersonDock+1,nPersonDock*2
    if (exclude(i)==0) then
      write(*,*) pop_tmp(i)%receptor,pop_tmp(i)%ligand,exclude(i),"running"
      write(99,"(A,1x,3(i5))") "./dock_and_rerank.bash ", i-nPersonDock, pop_tmp(i)%receptor,pop_tmp(i)%ligand
    else
      write(*,*) pop_tmp(i)%receptor,pop_tmp(i)%ligand,exclude(i),"NOT running"
    endif
  enddo
  write(99,"(A3)") "EOF"
  close(99)


!Actually run the dockings
  call system("cat runhex")
  call system("bash runhex")
  call Zrank_Get_RESULTS(generation,nPersonDock)
  call GetDuplicateResult
  call SetRank(generation)
  call StoreEnergytable(generation) !must come before SortFitness
  call SortFitness(generation)

!transfer best results from pop_tmp to pop
  if (elitism==0) then
    do i=1,nPersonDock
      pop(i)%receptor  = pop_tmp(i)%receptor
      pop(i)%ligand    = pop_tmp(i)%ligand
      pop(i)%fitness   = pop_tmp(i)%fitness
      do j=1,10
      pop(i)%energy(j) = pop_tmp(i)%energy(j) 
      pop(i)%rms(j)    = pop_tmp(i)%rms(j)
      enddo
    enddo
  else 
!with elitism
    do i=1,int(nPersonDock*elitism_rate)
      pop(i)%receptor  = pop_tmp(i)%receptor
      pop(i)%ligand    = pop_tmp(i)%ligand
      pop(i)%fitness   = pop_tmp(i)%fitness

      do j=1,10
        pop(i)%energy(j) = pop_tmp(i)%energy(j)
        pop(i)%rms(j)    = pop_tmp(i)%rms(j)
      enddo

    enddo

!WHERE IS THE COMPETITION
    call RouletteElitism(nPersonDock)

  endif

endif

1000  format("parallel -j1 <<EOF")
2000  format("parallel -j1 --sshloginfile $PBS_NODEFILE --wd $PWD <<EOF")
!2000  format("parallel -j",i5," --sshloginfile $PBS_NODEFILE --wd $PWD <<EOF") ! buggy
! running multiple instances of HEX on the same node is unstable
! I guess they're from too many attempts to access the same file
! Anyway, I'm sick of the errors.


END subroutine Dock_GeneratePopulation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine GADock
use vectors
use dock_vectors

open(66,file="docklogfile.dat")

crossover_method=1
EnergyTableLast=0

if ( (nPersonReceptor*nPersonLigand) < nPersonDock ) then
 nPersonDock=(nPersonReceptor*nPersonLigand)
endif

allocate(pop(nPersonDock))
allocate(pop_tmp(nPersonDock*2))
allocate(EnergyTable(nPersonDock*6)) 
allocate(exclude(nPersonDock*2))
allocate(duplicate(2,nPersonDock*2))

!reset EnergyTable
do i=1,nPersonDock*6
  EnergyTable(i)%receptor  = 0
  EnergyTable(i)%ligand    = 0
  EnergyTable(i)%fitness   = 0.0
  do j=1,10
  EnergyTable(i)%energy(j) = 0.0
  EnergyTable(i)%rms(j)    = 0.0
  enddo
enddo

if (EnergyTableFile=="") then
  EnergyTableFile="etable.dat"
  EnergyTableLast=0
endif
!Open and read energy table
open(55,file=EnergyTableFile)
!count number of lines
  do 100
    read(55,*,end=200,err=200) i
    EnergyTableLast=EnergyTableLast+1
100 continue
200 rewind(55)

  do i=1,EnergyTableLast
    read(55,*) EnergyTable(i)
  enddo

open(2,file="gadock.dat")
do i=0,ngenerations

  write(*,*) "Generation: ",i

  call Dock_GeneratePopulation(i)

  ! call StoreEnergytable(i) !must come before SortFitness

  call WriteLog(i)

enddo

close(2)

end subroutine GADock


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine CheckDuplicate
use vectors, only : nPersonDock
use dock_vectors

do i=1,nPersonDock*2
  duplicate(1,i) = 0
  duplicate(2,i) = 0
enddo

do i=nPersonDock+2,nPersonDock*2 !nao compara o 1o elemento neh
  do j=nPersonDock+1,i-1
    if (pop_tmp(i)%receptor==pop_tmp(j)%receptor.and.pop_tmp(i)%ligand==pop_tmp(j)%ligand) then
      write(*,*) "Found a duplicate in pop_tmp, position =",j," pair:",pop_tmp(j)%receptor,pop_tmp(j)%ligand
      duplicate(1,i)=i  !
      duplicate(2,i)=j  !
      exclude(i)=1      !exclude calculation
      goto 10
    endif
  enddo
10 continue
enddo

end subroutine CheckDuplicate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine GetDuplicateResult
use vectors, only : nPersonDock
use dock_vectors

do i=nPersonDock+2,nPersonDock*2
  if (duplicate(1,i)/=0) then
    pop_tmp(i) = pop_tmp(duplicate(2,i) ) !transfers whole object
  endif
enddo

end subroutine GetDuplicateResult

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!It also resizes nPersonLigand and nPersonReceptor
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine GeneratePopulation(gen)
use vectors
integer :: position
integer,intent(in)  :: gen
integer :: distcheck
real    :: sum

1 allocate(Receptor(nPersonReceptor))
do i=1,nPersonReceptor
  allocate(Receptor(i)%Disp(nModesReceptor))
enddo

if (gen>0) goto 100 !or else this is the initial population
!Assign random amplitudes to each mode.

!force no displacements AT ALL
i=1
    Receptor(i)%Disp(:)=0.0

!Receptor
do i=2,nPersonReceptor!nPerson 
  Receptor(i)%ind=i
  nfail=0
10 continue
  do j=1,nModesReceptor 
    position=int(rand(0)*NormalModeReceptor(j)%ndisp+1)
    Receptor(i)%Disp(j)=NormalModeReceptor(j)%list(position)
  enddo
    do j=1,i-1
        sum=0.0
        do nq=1,nModesReceptor
            sum= sum + ( Receptor(i)%Disp(nq) - Receptor(j)%Disp(nq) )**2
        enddo
        if ( sqrt(sum).lt.rmscutoff) then !acho que tem que dividir pelo nmodes
          nfail=nfail+1
          if(nfail>=1000) then
            nPersonReceptor=i-1 !adjust population size
            write(*,*) "Receptor population reduced to ",nPersonReceptor
            deallocate(Receptor)
            goto 1
          endif
          goto 10
        endif
    enddo
enddo

2 allocate(Ligand(nPersonLigand))
do i=1,nPersonLigand
  allocate(Ligand(i)%Disp(nmodes))
enddo

!force no displacements AT ALL
i=1
    Receptor(i)%Disp(:)=0.0

!Ligand
do i=2,nPersonLigand
  Ligand(i)%ind=i
  nfail=0
20 continue
  do j=1,nModesLigand 
    position=int(rand(0)*NormalModeLigand(j)%ndisp+1)
    Ligand(i)%Disp(j)=NormalModeLigand(j)%list(position)
  enddo
    do j=1,i-1
        sum=0.0
        do nq=1,nModesLigand
            sum= sum + ( Ligand(i)%Disp(nq) - Ligand(j)%Disp(nq) )**2
        enddo
        if ( sqrt(sum).lt.rmscutoff) then !acho que tem que dividir pelo nmodes
          nfail=nfail+1
          if(nfail>=1000) then
            nPersonLigand=i-1 !adjust population size
            write(*,*) "Ligand population reduced to ",nPersonLigand
            deallocate(Ligand)
            goto 2
          endif
          goto 20
        endif
    enddo
enddo

100 continue
return !Done with initial population
99  write(*,*) "FAIL!"
    write(*,*) "RMS too small", sqrt(sum),rmscutoff
    write(*,*) "Population can not be bigger than ",i-1
    STOP

    end subroutine GeneratePopulation!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine GenerateStructures
use vectors, only : nPerson,nmodes,Receptor,Ligand,nproc,nnodes,pbs
use vectors, only : nModesReceptor,nModesLigand
use vectors, only : nPersonReceptor,nPersonLigand

open(99,file="runcharmm")

1000  format("parallel -j",i5," <<EOF")
2000  format("parallel -j",i5," --sshloginfile $PBS_NODEFILE --wd $PWD <<EOF")
!3000  format("parallel -j1 --sshloginfile $PBS_NODEFILE -W$PWD <<EOF")

!if (pbs==1) write(99,*) "export CHARMMEXEC='mpirun -n 8 $CHARMMEXEC'"

if (pbs==0) write(99,1000) nproc
if (pbs==1) write(99,2000) nproc
if (pbs==2) write(99,2000) nnodes
!if (pbs==3) write(99,3000) !use MPI and nproc

do i=1,nPersonReceptor
  call runcharmm(nModesReceptor,i,Receptor(i)%disp(:),pbs,nproc,0) ! "0" for Receptor
enddo

do i=1,nPersonLigand
  call runcharmm(nModesLigand,i,Ligand(i)%disp(:),pbs,nproc,1)   ! "1" for Ligand
enddo
  
write(99,"(a3)") "EOF" 

close(99)

!Actually run charmm
call system("bash runcharmm")


end subroutine GenerateStructures




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine runcharmm(nmodes,pos,disp,pbs,nproc,rl)
!script is written in WRITE UNIT "99"
integer,intent(in)   :: pbs          ! Are we using a cluster ? ( 0 = no ; 1 = yes, grid ; 2 = yes, single job )
integer,intent(in)   :: nproc
integer,intent(in)   :: nmodes       ! Number of normal modes to consider
integer,intent(in)   :: pos          ! Position Receptor or Ligand
real,intent(in)      :: disp(nmodes) ! Amplitude for each normal mode
integer,intent(in)   :: rl           ! "0" for Receptor, "1" for Ligand
character            :: string*5
character            :: command*256
character            :: command2*256

  open(1,file="disp.dat")
  write(1,"(f8.3)") disp(:)
  close(1)
  write(unit=string, fmt='(I4)') pos !computes energy for "pos"

  if (rl==0) then
    command="./optimize.NM.displacement.bash > input/r_"//trim(adjustl(string))//".inp"
    command2="$CHARMMEXEC -i input/r_"//trim(adjustl(string))//".inp istr=r_"//trim(adjustl(string))//" -o output/r_"//trim(adjustl(string))//".log filename=Receptor dirnam=output"
  else
    command="./optimize.NM.displacement.bash > input/l_"//trim(adjustl(string))//".inp"
    command2="$CHARMMEXEC -i input/l_"//trim(adjustl(string))//".inp istr=l_"//trim(adjustl(string))//" -o output/l_"//trim(adjustl(string))//".log filename=Ligand dirnam=output"
  endif

!  if (pbs/=3) write(99,*) trim(command), " ; ", trim(command2)
!  if (pbs==3) write(99,"(A,A,A,i3,2x,A)") trim(command), " ; ", "mpirun -n ",nproc, trim(command2)

  call system(trim(command))

  if (pbs/=3) write(99,*) trim(command2)
  if (pbs==3) write(99,"(A,A,A,i3,2x,A)") "mpirun -n ",nproc, trim(command2)

end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine InitAmplitudes
use vectors
character   :: string
integer     :: usemaxminReceptor(nmodes) !> If (max-min) < gridpsacing use max min as values.
integer     :: usemaxminLigand(nmodes) !> If (max-min) < gridpsacing use max min as values.
real        :: MaxMinDiff !> difference between Maxdisp and Mindisp
real        :: MaxMinDiffOverGrid

usemaxminReceptor=0
usemaxminLigand=0

!The following should replace , in case no file is provided and only to compute the number of modes requested by the user.
!so far we will use the VERY SAME number for both Receptor and Ligand
!Receptor
if (trim(NormalModeReceptorFile)=="") then
  nModesReceptor=0
  write(*,*) "Receptor will not be optimized"
else
  call CountlinesNormalModeFile(nModesReceptor,NormalModeReceptorFile)
  
  allocate(NormalModeReceptor(nModesReceptor))

  call ReadNormalModeFile(nModesReceptor,&
                          NormalModeReceptorFile,&
                          NormalModeReceptor%ind,&
                          NormalModeReceptor%mindisp,&
                          NormalModeReceptor%maxdisp)
endif


do i=1,nModesReceptor

  MaxMinDiff= (NormalModeReceptor(i)%maxdisp - NormalModeReceptor(i)%mindisp )

  if  (MaxMinDiff==0.0) then
    NormalModeReceptor(i)%ndisp=1

    else
    call grid(NormalModeReceptor(i)%mindisp,&
              NormalModeReceptor(i)%maxdisp,&
              gridspacing,&
              NormalModeReceptor(i)%ndisp)
  endif

!Allocate LIST and generate a discrete amplitude distribution from 'min' to 'max'.
  allocate(NormalModeReceptor(i)%list  ( NormalModeReceptor(i)%ndisp   ) )

  call GenerateNormalModeList(gridspacing,&
                              NormalModeReceptor(i)%list,&
                              NormalModeReceptor(i)%ndisp,&
                              NormalModeReceptor(i)%mindisp,&
                              NormalModeReceptor(i)%maxdisp)
enddo

!same as above for Ligand
if (trim(NormalModeLigandFile)=="") then
  nModesLigand=0
  write(*,*) "Ligand will not be optimized"
else
  call CountlinesNormalModeFile(nModesLigand,NormalModeLigandFile)

!Allocate vector
  allocate(NormalModeLigand(nModesLigand))

  call ReadNormalModeFile(nModesLigand,&
                          NormalModeLigandFile,&
                          NormalModeLigand%ind,&
                          NormalModeLigand%mindisp,&
                          NormalModeLigand%maxdisp)
endif


do i=1,nModesLigand

  MaxMinDiff= (NormalModeLigand(i)%maxdisp - NormalModeLigand(i)%mindisp )

  if  (MaxMinDiff==0.0) then
    NormalModeLigand(i)%ndisp=1

    else
    call grid(NormalModeLigand(i)%mindisp,&
              NormalModeLigand(i)%maxdisp,&
              gridspacing,&
              NormalModeLigand(i)%ndisp)
  endif

!Allocate LIST and generate a discrete amplitude distribution from 'min' to 'max'.
  allocate(NormalModeLigand(i)%list  ( NormalModeLigand(i)%ndisp   ) )

  call GenerateNormalModeList(gridspacing,&
                              NormalModeLigand(i)%list,&
                              NormalModeLigand(i)%ndisp,&
                              NormalModeLigand(i)%mindisp,&
                              NormalModeLigand(i)%maxdisp)
enddo




end subroutine InitAmplitudes



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Computes ndisp
subroutine grid(x,y,z,i)
integer,intent(out)     :: i
real,intent(in)     :: x,y,z
real    ::  d
real    ::  a
d=y-x

!write(*,*) x,y,d

i=2 !ja inicia com o max e min
a=z
do while (a<d)
   a=a+z
   i=i+1
enddo

return
end subroutine grid


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine CountLinesNormalModeFile(i,filename)
integer, intent(out)    :: i
character(len=128),intent(in)    :: filename
character   ::  string

i=0
  open(1,file=trim(filename))
     do
       read(1,*,end=10,err=11) string
       if (string=="#") then
         cycle
       else
          i=i+1
       endif
     enddo
10   close(1) ; return
11   write(*,*) "CountLinesNormalModeFile ended with error"
STOP
end subroutine CountLinesNormalModeFile


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine ReadNormalModeFile(nm,&
                              filename,&
                              ind,&
                              mindisp,&
                              maxdisp)
integer, intent(in)    :: nm
character(len=128),intent(in)    :: filename
character   ::  string*1
integer :: j
integer,intent(out) :: ind(nm)      !> index
real    :: maxdisp(nm)  !> Maximum displacement
real    :: mindisp(nm)  !> min
j=0
ind(:)=0
maxdisp(:)=0.0
mindisp(:)=0.0

write(*,*) "Reading File: ",trim(filename)
open(1,file=trim(filename))

do while (j<nm)
  read(1,*) string
  if (string=="#") cycle
  backspace(1)
  j=j+1 !modnum
  read(1,*) ind(j), mindisp(j), maxdisp(j)
  write(*,*)ind(j), mindisp(j), maxdisp(j)
enddo
close(1)
return

end subroutine ReadNormalModeFile


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine GenerateNormalModeList(gridspacing,&
                                  list,&
                                  ndisp,&
                                  mindisp,&
                                  maxdisp)
real,intent(in)    :: gridspacing
real,intent(out)   :: list(ndisp)
integer,intent(in) :: ndisp
real,intent(in)    :: mindisp
real,intent(in)    :: maxdisp

list=0.0

  if( ndisp==1) then
    list(1) = mindisp
  else
    list(1) = mindisp
    do j=2,ndisp-1
      list(j) = list(j-1)+gridspacing
    enddo
    list(j) = maxdisp
  endif
end subroutine GenerateNormalModeList
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Roulette(crossover_1,crossover_2)
use vectors, only : nPersonDock
use dock_vectors
real        ::  number_1
real        ::  number_2
real        ::  tmp_sum
integer     ::  counter=0
integer,intent(out) :: crossover_1
integer,intent(out) :: crossover_2

write(*,"(A)") "Roulette"
 
  number_1  =   (rand(0)*Fitness_sum)+1
  tmp_sum=0.0d0

do i=1,nPersonDock
  crossover_1=i
  tmp_sum = tmp_sum + pop(i)%Fitness
  if( tmp_sum > number_1 ) exit
enddo

100   number_2  =   (rand(0)*Fitness_sum)+1.0
 tmp_sum=0.0d0

do i=1,nPersonDock
  crossover_2=i
  tmp_sum = tmp_sum + pop(i)%Fitness
  if( tmp_sum > number_2 ) exit
enddo
  if( crossover_1 == crossover_2 ) then
  counter = counter + 1
  if (counter > nPersonDock) then
    write(*,*) "[Roulette routine] Trapped in an infinite loop, will accepting any crossover."
    write(*,*) "[Roulette routine] Don’t worry, population will converge anyway."
    goto 200
  endif
  goto 100
  endif

200 continue
END subroutine Roulettesubroutine SetDefaults
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

end subroutine SetDefaults! Rank structures by fitness but does not sort
SUBROUTINE SetRank(generation)
use vectors, only : nPersonDock
use dock_vectors
real  ::  max_energy  ! it is the worst energy
real  ::  fitness     ! it is the worst energy

if (generation==0) then
  max_energy     =   maxval( pop(:)%energy(1) )
  do i=1,nPersonDock
    fitness = max_energy - pop(i)%energy(1)
    pop(i)%fitness = abs (fitness)
  enddo
  fitness_sum = sum(pop(:)%fitness)
endif

if (generation>0) then
  max_energy     =   maxval( pop_tmp(:)%energy(1) )
  do i=1,nPersonDock*2
    fitness = max_energy - pop_tmp(i)%energy(1)
    pop_tmp(i)%fitness = abs (fitness)
  enddo
  fitness_sum = sum(pop_tmp(:)%fitness)
endif
end subroutine SetRank
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine SortFitness(generation)
use vectors, only : nPersonDock
use dock_vectors
type(valor) ::  pop_swap

if (generation==0) then
  do i=1,nPersonDock-1
     do j=1,nPersonDock-1
          if( (pop(j)%fitness) < (pop(j+1)%fitness) ) then
              pop_swap    =   pop(j)      ! transfers the whole object
              pop(j)      =   pop(j+1)    ! transfers the whole object
              pop(j+1)    =   pop_swap    ! transfers the whole object
          endif
    enddo
 enddo
endif

if (generation>0) then
 do i=1,nPersonDock*2-1
    do j=1,nPersonDock*2-1
          if( (pop_tmp(j)%fitness) < (pop_tmp(j+1)%fitness) ) then
              pop_swap      =   pop_tmp(j)      ! transfers the whole object
              pop_tmp(j)    =   pop_tmp(j+1)    ! transfers the whole object
              pop_tmp(j+1)  =   pop_swap        ! transfers the whole object
          endif
    enddo
 enddo

endif

end subroutine SortFitness
subroutine StoreEnergytable(generation)
use vectors!, only : nPerson,EnergyTableLast
use dock_vectors
implicit none
integer             :: i
integer,intent(in)  :: generation

if (generation==0) then
  do i=1,nPersonDock
    if (exclude(i)==0) then
      EnergyTableLast=EnergyTableLast+1
      EnergyTable(EnergyTableLast)=pop(i) ! transfers the whole object
      write(55,*) EnergyTable(EnergyTableLast)
      write(*,*) "Etable",EnergyTableLast,"rec",EnergyTable(EnergyTableLast)%receptor,"lig",EnergyTable(EnergyTableLast)%ligand
   endif
  enddo

else

  do i=nPersonDock+1,nPersonDock*2
    if (exclude(i)==0) then
      EnergyTableLast=EnergyTableLast+1
      EnergyTable(EnergyTableLast)=pop_tmp(i) ! transfers the whole object
      write(55,*) EnergyTable(EnergyTableLast)
      write(*,*) "Etable",EnergyTableLast,"rec",EnergyTable(EnergyTableLast)%receptor,"lig",EnergyTable(EnergyTableLast)%ligand
    endif
  enddo

endif

end subroutine StoreEnergytablesubroutine Zrank_Get_RESULTS(generation,nPersonDock)
use dock_vectors
integer,intent(in)   :: nPersonDock
integer,intent(in)   :: generation
character(len=5)     :: receptor
character(len=5)     :: ligand
character(len=40)    :: trash
character(len=256)   :: arquivo

if (generation==0) then
do i=1,nPersonDock
  if (exclude(i)==0) then
    write(receptor,"(i5)") (pop(i)%receptor)
    write(ligand,"(i5)")   (pop(i)%ligand)
    arquivo="dockout/"//trim(adjustl(receptor))//"_"//trim(adjustl(ligand))//".list.zr.reranked"
    open(99,file=trim(arquivo))
    pop(i)%energy(:)=0.0 !acochambracao p caso ter menos de 10
    do j=1,10 !10best
      read(99,*,end=11,err=11) trash, pop(i)%energy(j)
    enddo ; goto 11
11  continue
    close(99)
  endif
enddo

else

do i=nPersonDock+1,nPersonDock*2
  if(exclude(i)==0) then
    write(receptor,"(i5)") (pop_tmp(i)%receptor)
    write(ligand,"(i5)")   (pop_tmp(i)%ligand)
    arquivo="dockout/"//trim(adjustl(receptor))//"_"//trim(adjustl(ligand))//".list.zr.reranked"
    open(99,file=trim(arquivo))
    pop_tmp(i)%energy(:)=0.0 !acochambracao p caso ter menos de 10
    do j=1,10 !10best
      read(99,*,end=13,err=13) trash, pop_tmp(i)%energy(j)
    enddo ; goto 13
13  continue
    close(99)
  endif
enddo
endif

end!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine crossover
use dock_vectors
use vectors, only : nPersonDock
integer ::  crossover_1
integer ::  crossover_2
integer ::  posicao

write(*,"(A)") "crossover"

do i=1,nPersonDock,2

  call roulette(crossover_1,crossover_2)

  !crossover
  pop_tmp(i+nPersonDock)%receptor    =   pop(crossover_1)%receptor
  pop_tmp(i+nPersonDock)%ligand      =   pop(crossover_2)%ligand

  pop_tmp(i+nPersonDock+1)%receptor  =   pop(crossover_2)%receptor
  pop_tmp(i+nPersonDock+1)%ligand    =   pop(crossover_1)%ligand

enddo

write(*,*) "Parents"
do i=1,nPersonDock
  write(*,*) pop_tmp(i)%receptor,pop_tmp(i)%ligand
enddo

write(*,*) "offspring"
do i=1,nPersonDock
 write(*,*) pop_tmp(i+nPersonDock)%receptor,pop_tmp(i+nPersonDock)%ligand
enddo

END subroutine crossover
SUBROUTINE run_hex(i,generation,receptor,ligand)
use vectors , only   :  nPersonDock,pbs,nproc,nnodes
use dock_vectors, only  : exclude
integer,intent(in)   :: i
integer,intent(in)   :: generation
integer,intent(in)   :: receptor
integer,intent(in)   :: ligand

!if(generation==0) then
  if (i==1) then
    open(99,file="runhex")
    if (pbs==0) write(99,1000)
    if (pbs==1) write(99,2000)
    if (pbs==2) write(99,3000)
  endif

   if(exclude(i)==0)  write(99,"(A,1x,3(i5))") "./dock_and_rerank.bash ", i, receptor, ligand

  if(i==nPersonDock) then
!    if(pbs/=0) write(99,"(A3)") "EOF"
    write(99,"(A3)") "EOF"
    close(99)
  endif
!endif

1000  format("parallel -j1 <<EOF")
2000  format("parallel -j1 --sshloginfile $PBS_NODEFILE --wd $PWD <<EOF")
3000  format("parallel -j1 --sshloginfile $PBS_NODEFILE --wd $PWD <<EOF")
end subroutine run_hex

subroutine WriteHex
write(*,*)
end subroutine WriteHex
!@author
!> Diego Enry B. Gomes
!@brief
!>Error Handling
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Error(n,arg1,arg2)
use colors
integer,intent(in) :: n
character(len=30)    :: arg1,arg2
write(*,*) "---------------------------------------------------------------"
write(*,*) "Program stopped with the following error:"

select case(n)
  case(1)
    write(*,*) "No arguments."

  case(2)
    write(*,*) "No arguments for ",red,arg1

  case(3)
    write(*,*) "ERROR: Improper value for variable '",green,trim(arg1),normal,"' => '",red,trim(arg2),normal,"'"

  case(4)
    write(*,*) "ERROR: -nmrec and -nmlig mandatory" !by now
end select

write(*,*) normal
write(*,*) "run "//red//"./gamini"//green//" -h"//normal//" for help"

STOP
end subroutine

subroutine ReadCommandLine
use vectors
integer         ::  n
integer         ::  iargc
character(len=30) :: argv
!number of provided arguments
n = iargc() 

NormalModeReceptorFile=""
NormalModeLigandFile=""

!Check if Help was requested
do i = 1, n
  call getarg( i, argv )  
  if(trim(argv).eq."-h") call help
enddo

!Get the argument values
do i = 1, n
  call getarg( i, argv )  
  if(trim(argv).eq."-pbs")      call GetVal(i,argv,pbs)
  if(trim(argv).eq."-ngen")     call GetVal(i,argv,ngenerations)
  if(trim(argv).eq."-psize")    call GetVal(i,argv,nPerson)
  if(trim(argv).eq."-nbest")    call GetVal(i,argv,nBest)
  if(trim(argv).eq."-nmodes")   call GetVal(i,argv,nModes)
  if(trim(argv).eq."-grid")     call dGetVal(i,argv,gridspacing)
  if(trim(argv).eq."-rms")      call dGetVal(i,argv,rmscutoff)
  if(trim(argv).eq."-onlydock") onlydock=1
  if(trim(argv).eq."-nodock")   nodock=1
  if(trim(argv).eq."-etable")   call StrGetVal(i,argv,EnergyTableFile)
  if(trim(argv).eq."-nmrec")    call StrGetVal(i,argv,NormalModeReceptorFile)
  if(trim(argv).eq."-nmlig")    call StrGetVal(i,argv,NormalModeLigandFile)
  if(trim(argv).eq."-pdock")    call GetVal(i,argv,nPersonDock)
  if(trim(argv).eq."-nproc")    call GetVal(i,argv,nproc)
  if(trim(argv).eq."-nnodes")    call GetVal(i,argv,nnodes)


!  if(trim(argv).eq."-nmodes") call GetVal(i,argv,nmodes) !argument is real
enddo

!Failsafe procedures
if(nPerson<nBest) then
  write(*,*)
  write(*,*) "WARNING ------------------------------------"
  write(*,*) "nbest can not exceed population size (psize)"
  write(*,*) "Setting nbest to:", nPerson
  nBest=nPerson
  write(*,*)
endif

!!!if(NormalModeReceptorFile=="".or.NormalModeLigandFile=="") call error(4)

!end of Failsafe

!Writeout run parameters
if(n.eq.0) then 
  write(*,*) "Running with defaults:"
else
  write(*,*) "Running these parameters:"
endif
  write(*,*) "pbs    = ",pbs
  write(*,*) "ngen   = ",nGenerations
  write(*,*) "psize  = ",nPerson
  write(*,*) "grid   = ",GridSpacing
  write(*,*) "nnodes = ",nNodes
  write(*,*) "nmodes = ",nModes
  write(*,*) "nbest  = ",nBest
  write(*,*) "pdock  = ",nPersonDock
end


subroutine GetVal(i,inargv,nval)
character(len=30),intent(in)    :: inargv
character(len=30)    :: argv
integer,intent(out)  :: nval
integer,intent(in)   :: i
  call getarg(i+1,argv)
!  if(argv.eq." ".or.argv(1:1).eq."-") call error(2,inargv)
  if(trim(argv).eq."") call error(2,inargv)
  read(argv,*,err=10) nval !string to number
return
10 call error(3,inargv,argv)
end subroutine


subroutine dGetVal(i,inargv,nval)
character(len=30),intent(in)    :: inargv
character(len=30)    :: argv
real,intent(out)     :: nval
integer,intent(in)   :: i
  call getarg(i+1,argv)
!  if(argv.eq." ".or.argv(1:1).eq."-") call error(2,inargv)
  if(trim(argv).eq."") call error(2,inargv)
  read(argv,*,err=10) nval !string to number
return
10 call error(3,inargv,argv)
end subroutine


subroutine strGetVal(i,inargv,strArg)
character(len=30)    :: inargv
character(len=30)    :: argv
character(len=30),intent(out)  :: strArg
integer,intent(in)   :: i
  call getarg(i+1,argv)
!  if(argv.eq." ".or.argv(1:1).eq."-") call error(2,inargv)
  if(trim(argv).eq."") call error(2,inargv)
  strArg=argv
return
end subroutine
Subroutine RouletteElitism(nPersonDock) ! Competicao
use dock_vectors
implicit none
integer  :: i,j,k
real     :: number_1
real     :: FitnessSumElitism=0.0
real     :: tmp_sum
integer,intent(in)   :: nPersonDock

!computa a soma do fitness para o resto da populacao.
do i=(int(nPersonDock*ElitismRate)+1),nPersonDock*2
  FitnessSumElitism = pop_tmp(i)%fitness + FitnessSumElitism
enddo

!faz para todos elementos do "resto"
do j=(int(nPersonDock*ElitismRate)+1),nPersonDock
number_1 = ( rand(0)*FitnessSumElitism ) + 1
tmp_sum=0.0

!desta forma podem haver duplicatas, vamos ver se vai bem se não mudamos.
  do i=(int(nPersonDock*ElitismRate)+1), nPersonDock*2
    tmp_sum = tmp_sum + pop_tmp(i)%fitness
    
    if(tmp_sum > number_1) then

      pop(j)%receptor  = pop_tmp(i)%receptor
      pop(j)%ligand    = pop_tmp(i)%ligand
      pop(j)%fitness   = pop_tmp(i)%fitness
    
      do k=1,10
        pop(j)%energy(j) = pop_tmp(i)%energy(k)
        pop(j)%rms(j)    = pop_tmp(i)%rms(k)
      enddo
      
    endif

  enddo

enddo

end subroutine RouletteElitism
subroutine WriteAmplitudes
use vectors

write(*,*) "Avaible values for each mode for: Receptor "
do i=1,nmodes
  write(*,*) "MODE ",i+6, NormalModeReceptor(i)%list(:)
enddo

write(*,*)
write(*,*) "Chosed values for each Receptor"
write(*,*) "Receptor  | ","  MODE 7 | ","  MODE 8 |","  MODE 9"
do i=1,nPerson
  write(*,*) i, Receptor(i)%Disp(:)
enddo

write(*,*)
write(*,*) "Avaible values for each mode for: Ligand "
do i=1,nmodes
  write(*,*) "MODE ",i+6, NormalModeLigand(i)%list(:)
enddo

write(*,*)
write(*,*) "Chosed values for each Ligand"
write(*,*) "Ligand  | ","  MODE 7  | ","  MODE 8  |","  MODE 9"
do i=1,nPerson
  write(*,*) i, Ligand(i)%Disp(:)
enddo

end subroutine WriteAmplitudes
subroutine Help
use vectors

write(*,*) "Usage options:"
write(*,*)
write(*,*) "Options for the Genetic Algorithm"
write(*,*) "-ngen   = Number of Generations"
write(*,*) "-psize  = Population Size"
write(*,*)
write(*,*) "Options for the Normal Mode Search"
write(*,*) "-grid   = Gridspacing between two displacements"
write(*,*) "-nmodes = Number of Normal Modes to optimize"
write(*,*) "-nbest  = Number of Best Structures to keep"
write(*,*) "-nmrec  = File containing Normal Modes to optimize for the Receptor"
write(*,*) "-nmlig  = File containing Normal Modet to optimize for the Ligand"
write(*,*) 
write(*,*) "-pdock  = nPersonDock"
write(*,*) 
write(*,*) "Options for the running on clusters"
write(*,*) "-pbs    = Run locally (0, default) or using a PBS cluster (1)"
write(*,*) "-nnodes = Number of compute nodes to run (default=1)"
write(*,*) "-nproc  = Number of processors at each node (default=1)"
write(*,*)
write(*,*) "Advanced Options"
write(*,*) "-onlydock = Only dock"
write(*,*) "-nodock   = Do not dock"
write(*,*) "-etable   = Provide an energy table file"
write(*,*) 
write(*,*) "The input file for the receptor or ligand normal mode (-nmrec, -nmlig)"
write(*,*) "should be as follows, and '#' serve as comment"
write(*,*) "# mode | min | max"
write(*,*) " 7  -0.1   0.1"
write(*,*) " 8   0.0   0.1"
write(*,*) "# 9   0.0   0.0 (will ignore this mode)"
write(*,*) "10  -0.5   0.9"

STOP
end subroutine Help
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine credits

write(*,*) "###############################################################"
write(*,*) "#                   GANMDOCK - version 1.0rc3                 #"
write(*,*) "#               Mon Nov 19 16:46:08 BRST 2012                 #"
write(*,*) "#                                                             #"
write(*,*) "# Diego E B Gomes(1,2), Pedro G Pascutti(2),Paulo M Bisch(2), #"
write(*,*) "# David Perahia(3), Luis P B Scott(1)  | diego@biof.ufrj.br   #"
write(*,*) "#                                                             #"
write(*,*) "# 1 - Universidade Federal do ABC                             #"
write(*,*) "# 2 - Univerisdade Federal do Rio de Janeiro                  #"
write(*,*) "# 3 - Ecole Normale Superieure de Cachan                      #"
write(*,*) "###############################################################"

end


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine ShowDisplacements
use vectors
  write(*,*) "Receptor"
  do i=1,nModesReceptor
    write(*,*) "mode:",NormalModeReceptor(i)%ind
    write(*,*) "disp:",NormalModeReceptor(i)%list
    write(*,*)
    enddo
  write(*,*) "Ligand"
  do i=1,nModesLigand
    write(*,*) "mode:",NormalModeLigand(i)%ind
    write(*,*)"disp:",NormalModeLigand(i)%list
    write(*,*)
  enddo
end subroutine ShowDisplacements



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine ShowPopulation
use vectors

  write(*,*) "Receptors"
  
  do i=1,nPersonReceptor
    write(*,*) "receptor:",i,(Receptor(i)%Disp(j),j=1,nModesReceptor)
  enddo

  write(*,*) "Ligands"
  do i=1,nPersonLigand
    write(*,*) "ligand:",i,(Ligand(i)%Disp(j),j=1,nModesLigand)
  enddo

end subroutine ShowPopulation


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine ShowDockingPopulation
use vectors
use dock_vectors
  do i=1,nPersonDock
    write(*,*) i,pop(i)%receptor,pop(i)%ligand
  enddo
end subroutine ShowDockingPopulation



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine ShowEnergyTable
use vectors
use dock_vectors
write(*,*) "Etable has ", EnergyTableLast, "elements"
  do j=1,EnergyTableLast
    write(*,*)  i,&
                EnergyTable(j)%receptor,&
                EnergyTable(j)%ligand,&
                EnergyTable(j)%energy(1),&
                EnergyTable(j)%fitness
  enddo
end subroutine ShowEnergyTable



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine WriteLog(g)
use dock_vectors
use vectors
character(len=128)   :: filename
character(len=128)   :: filename2
character(len=128)   :: filename3
character(len=128)   :: command
integer,intent(in)  :: g
integer  :: p !population size ( 0 or >0)
character(len=5) :: rec
character(len=5) :: lig
real  :: rms(1000)
real  :: energy(1000)
real  :: zrank(1000)
character(len=128)  :: str(1000)
real  :: isscore(1000)
real  :: isscoreRMS(1000)
real  :: isscoreIRMS(1000)

if (g==0) p=nPersonDock
if (g/=0)  p=nPersonDock*2


write(66,*) "# Generation ",g

DO j=1,p

write(66,*) "#",pop_tmp(j)%receptor,pop_tmp(j)%ligand
write(unit=rec, fmt='(I5)') pop_tmp(j)%receptor
write(unit=lig, fmt='(I5)') pop_tmp(j)%ligand

filename="dockout/"//trim(adjustl(rec))//"_"//trim(adjustl(lig))//".hex"
filename2="dockout/"//trim(adjustl(rec))//"_"//trim(adjustl(lig))//".list.zr.out"

open(1,file=trim(filename))

do 
read(1,*) str(1)
  if (str(1)=="1") exit
enddo
backspace(1)

n=0
do i=1,100
  read(1,*,end=10,err=10) ind,nlixo,rms(i),energy(i)
  n=n+1
enddo
10 continue
close(1)

open(1,file=trim(filename2))
do i=1,n
  read(1,*) str(i),zrank(i)
enddo
close(1)

goto 99 ! This is not working

filename3="dockout/"//trim(adjustl(rec))//"_"//trim(adjustl(lig))//".isscore.rms"
open(1,file=trim(filename3))
do i=1,n
  read(1,*) isscoreRMS(i)
enddo
close(1)

filename3="dockout/"//trim(adjustl(rec))//"_"//trim(adjustl(lig))//".isscore.irms"
open(1,file=trim(filename3))
do i=1,n
  read(1,*) isscoreIRMS(i)
enddo
close(1)

!filename3="dockout/"//trim(adjustl(rec))//"_"//trim(adjustl(lig))//".isscore"
!open(1,file=trim(filename3))
!do i=1,n
!  read(1,*) isscore(i)
!enddo
!close(1)

99 continue 
!str = docking result
!
do i=1,n
!  write(66,*) str(i),energy(i),zrank(i),isscore(i),rms(i),isscoreRMS(i),isscoreIRMS(i)
  write(66,*) trim(adjustl(str(i))),energy(i),zrank(i),rms(i)!,isscoreRMS(i),isscoreIRMS(i)
enddo

ENDDO

end subroutine WriteLog
program GANMDOCK
use vectors
  
  call credits

  call SetDefaults

  call ReadCommandLine

  call InitAmplitudes

  !calcula o tamanho da nPerson para caber essas diferencas.
  call CalcNPerson

  call GeneratePopulation(0)
 
!Compute the structures.

  if(onlydock==0) call GenerateStructures(0)
  if(onlydock==0) call Dock_FixPDB(nPersonReceptor,nPersonLigand)

!Dock the structures.
  if(nodock==0) call GADock

end



subroutine CalcNPerson
use vectors
integer     ::  rtemp

nPersonReceptor=1
do i=1,nModesReceptor
  call fac(NormalModeReceptor(i)%ndisp,rtemp)
  nPersonReceptor=nPersonReceptor*rtemp
  if ( nPersonReceptor >= nPerson ) then
    nPersonReceptor=nPerson
    exit
  endif
enddo

nPersonLigand=1
do i=1,nModesLigand
  call fac(NormalModeLigand(i)%ndisp,rtemp)
  nPersonLigand=nPersonLigand*rtemp
  if ( nPersonLigand >= nPerson ) then
    nPersonLigand=nPerson
    exit
  endif
enddo

  write(*,*) "Population Sizes"
  write(*,*) "Receptor:",nPersonReceptor
  write(*,*) "Ligand:",nPersonLigand

end subroutine CalcNPerson

subroutine fac(np,r)
integer,intent(in) :: np
integer :: n
integer,intent(out) :: r
n=np
r=1
if (n==1) goto 10
if (n>=8) then ! se for fatorial mt grande 
  r=100000
  return
endif
do while (n>=1)
 r=r*n
 n=n-1
enddo

10 return
end
