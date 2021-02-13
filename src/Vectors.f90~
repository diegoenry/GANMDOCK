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
end module colors