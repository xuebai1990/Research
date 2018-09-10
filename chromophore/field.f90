program main
implicit none

character(len=4) :: filename
integer :: nconfig,nmol,i,j,k,num,tconfig
real*8 :: hdis,dis1,dis2
real*8 :: chromx(21),chromy(21),chromz(21),wx(500,4),wy(500,4),wz(500,4),r,f
real*8 :: a,b,c,d,e,g,inv(3,3),XY(3)
real*8,allocatable :: field(:,:),exe(:)
real*8,parameter :: cO=-1.04d0,cH=0.52d0
real*8 :: gc(21),ec(21)

namelist /para/ nconfig,nmol,tconfig

open(unit=10,file='input.dat')
read(10,NML=para)
close(10)

open(unit=30,file='solfield.dat')
open(unit=40,file='nosort1600.dat')
open(unit=50,file='XY.dat')
open(unit=60,file='X.dat')
open(unit=70,file='Y.dat')

allocate(field(nconfig,3),exe(nconfig))

gc(1)=-0.132865;ec(1)=-0.066282
gc(2)=0.113169;ec(2)=0.098383
gc(3)=-0.010488;ec(3)=-0.01147
gc(4)=-0.15142;ec(4)=-0.049799
gc(5)=-0.093404;ec(5)=-0.095423
gc(6)=0.097812;ec(6)=0.1144
gc(7)=-0.046238;ec(7)=-0.14477
gc(8)=0.088212;ec(8)=0.115383
gc(9)=0.102122;ec(9)=0.10706
gc(10)=-0.102322;ec(10)=-0.113737
gc(11)=0.065593;ec(11)=0.015504
gc(12)=0.119949;ec(12)=0.085561
gc(13)=0.115173;ec(13)=0.105958
gc(14)=0.129747;ec(14)=0.109494
gc(15)=-0.194628;ec(15)=-0.260056
gc(16)=-0.112117;ec(16)=-0.12383
gc(17)=0.122929;ec(17)=0.10815
gc(18)=0.123079;ec(18)=0.108673
gc(19)=0.104684;ec(19)=0.103484
gc(20)=0.096809;ec(20)=0.123296
gc(21)=-0.435728;ec(21)=-0.329935

do j=1,nconfig
   do k=1,3
      field(j,k)=0.0
   end do
end do

do i=1,nconfig
   read(40,*) num,exe(i)
   write(70,*) exe(i)
   write(filename,"(I4)") num
   open(unit=20,file='QB_WATER'//trim(adjustl(filename))//'.xyz')
   do j=1,21
      read(20,*) chromx(j),chromy(j),chromz(j)
   end do
   do j=1,nmol
      do k=1,4
         read(20,*) wx(j,k),wy(j,k),wz(j,k)
      end do
   end do
   close(20)

   do j=1,21
      do k=1,nmol
         d=chromx(j)-wx(k,4)
         e=chromy(j)-wy(k,4)
         g=chromz(j)-wz(k,4)
         r=sqrt(d**2+e**2+g**2)
         a=cO*gc(j)/r
         field(i,1)=field(i,1)+a
         d=chromx(j)-wx(k,2)
         e=chromy(j)-wy(k,2)
         g=chromz(j)-wz(k,2)
         r=sqrt(d**2+e**2+g**2)
         a=cH*gc(j)/r
         field(i,1)=field(i,1)+a
         d=chromx(j)-wx(k,3)
         e=chromy(j)-wy(k,3)
         g=chromz(j)-wz(k,3)
         r=sqrt(d**2+e**2+g**2)
         a=cH*gc(j)/r
         field(i,1)=field(i,1)+a
         d=chromx(j)-wx(k,4)
         e=chromy(j)-wy(k,4)
         g=chromz(j)-wz(k,4)
         r=sqrt(d**2+e**2+g**2)
         a=cO*ec(j)/r
         field(i,2)=field(i,2)+a
         d=chromx(j)-wx(k,2)
         e=chromy(j)-wy(k,2)
         g=chromz(j)-wz(k,2)
         r=sqrt(d**2+e**2+g**2)
         a=cH*ec(j)/r
         field(i,2)=field(i,2)+a
         d=chromx(j)-wx(k,3)
         e=chromy(j)-wy(k,3)
         g=chromz(j)-wz(k,3)
         r=sqrt(d**2+e**2+g**2)
         a=cH*ec(j)/r
         field(i,2)=field(i,2)+a
      end do
   end do
   field(i,3)=1.0d0   
end do

do i=1,3
   do j=1,3
      inv(i,j)=0.0
   end do
end do

!do i=1,nconfig
!   do j=1,64
!      do k=1,64
!         inv(j,k)=inv(j,k)+field(i,j)*field(i,k)
!      end do
!   end do
!end do
inv=matmul(transpose(field),field)
!field=field/nconfig

do i=1,3
   XY(j)=0.0
end do

!do i=1,nconfig
!   do j=1,64
!      XY(j)=XY(j)+field(i,j)*exe(i)
!   end do
!end do
XY=matmul(exe,field)
!XY=XY/nconfig


do i=1,3
   write(30,"(3F20.10)") inv(i,:)
end do

write(50,"(3F20.10)") XY

do i=1,nconfig
   write(60,"(3F20.10)") field(i,:)
end do

close(30)
close(40)
close(50)
close(60)
close(70)

end program
   
