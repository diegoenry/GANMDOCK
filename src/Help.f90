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
