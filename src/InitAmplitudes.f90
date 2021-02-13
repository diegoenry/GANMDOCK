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
