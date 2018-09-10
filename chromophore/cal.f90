program main
implicit none

real*8 :: XX(3,3),XY(3),beta(3)
real*8,allocatable :: X(:,:),Yt(:),Y(:)
integer :: i,j,nconfig,nmol,tconfig

namelist /para/ nconfig,nmol,tconfig

open(unit=10,file='input.dat')
read(10,NML=para)
close(10)
allocate(X(nconfig,3),Yt(nconfig),Y(nconfig))

open(unit=20,file='finverse.dat')
open(unit=30,file='XY.dat')
open(unit=40,file='beta.dat')
open(unit=50,file='X.dat')
open(unit=60,file='Yt.dat')
open(unit=70,file='Y.dat')

do i=1,3
   read(20,*) XX(i,:)
end do

read(30,*) XY

do i=1,nconfig
   read(50,*) X(i,:)
end do

do i=1,nconfig
   read(70,*) Y(i)
end do

do i=1,3
   beta(i)=0.0
end do

beta=matmul(XY,XX)

do i=1,3
   write(40,*) beta(i)
end do


Yt=matmul(beta,transpose(X))
do i=1,nconfig
   write(60,*) Yt(i),Y(i)
end do

end program

