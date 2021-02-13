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

    end subroutine GeneratePopulation