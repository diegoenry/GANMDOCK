!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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

end module dock_vectors