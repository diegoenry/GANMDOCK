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