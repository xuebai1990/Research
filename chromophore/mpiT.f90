program main
include 'mpif.h'

interface
   function inverse(A) result(Ainv)
      real*8, dimension(:,:), intent(in) :: A
      real*8, dimension(size(A,1),size(A,2)) :: Ainv
   end function
end interface

character(len=4) :: filename
integer :: nconfig,nmol,i,j,k,tconfig,l,ic,m,count
real*8 :: chromx(21),chromy(21),chromz(21),wx(500,4),wy(500,4),wz(500,4),rx,ry,rz,r,charge,cri
real*8 :: indp(4500),elec0(4500),sij,vij,elecQB(4500)
real*8 :: inv(5,5),XY(5)
real*8 :: mA(4500,4500),invmA(4500,4500)
real*8,allocatable :: field(:,:),exe(:),lfield(:,:)
integer,allocatable :: num(:)
real*8 :: cw(4),poQB(21),pow(4)
real*8 :: gc(21),ec(21)
integer :: numtasks,taskid,len,ierr,per
integer status(MPI_STATUS_SIZE)
parameter (MASTER=0)

namelist /para/ nconfig,nmol,tconfig

open(unit=10,file='input.dat')
read(10,NML=para)
close(10)

open(unit=30,file='solfield.dat')
open(unit=40,file='nosort1600.dat')
open(unit=50,file='XY.dat')
open(unit=60,file='X.dat')
open(unit=70,file='Y.dat')

allocate(field(5,nconfig),exe(nconfig))
allocate(num(nconfig))

cw(1)=0.00;cw(2)=0.52;cw(3)=0.52;cw(4)=-1.04
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

poQB(1)=1.5064
poQB(2)=1.5064
poQB(3)=1.5064
poQB(4)=1.5064
poQB(5)=1.5064
poQB(6)=0.5186
poQB(7)=1.5064
poQB(8)=0.5186
poQB(9)=0.5186
poQB(10)=1.5064
poQB(11)=1.5064
poQB(12)=0.5186
poQB(13)=0.5186
poQB(14)=0.5186
poQB(15)=1.1258
poQB(16)=1.5064
poQB(17)=0.5186
poQB(18)=0.5186
poQB(19)=0.5186
poQB(20)=1.5064
poQB(21)=0.9465
pow(1)=0.9465
pow(2)=0.5186
pow(3)=0.5186
pow(4)=0.0000

do j=1,5
   do k=1,nconfig
      field(j,k)=0.0
   end do
end do

call MPI_INIT(ierr)
call MPI_COMM_SIZE(MPI_COMM_WORLD,numtasks,ierr)
call MPI_COMM_RANK(MPI_COMM_WORLD,taskid,ierr)

do i=1,nconfig
   read(40,*) num(i),exe(i)
   if(taskid.eq.MASTER) then
      write(70,*) exe(i)
   end if
end do

per=nconfig/numtasks
allocate(lfield(5,per))

do i=1,5
   call MPI_SCATTER(field(i,:),per,MPI_DOUBLE,lfield(i,:),per,MPI_DOUBLE,MASTER,MPI_COMM_WORLD,ierr)
end do

do i=1,per
   write(*,*) per*taskid+i
   write(filename,"(I4)") num(per*taskid+i)
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

!Contribution from electrostatic interaction
   do j=1,21
      do k=1,nmol
         do l=2,4
            rx=chromx(j)-wx(k,l)
            ry=chromy(j)-wy(k,l)
            rz=chromz(j)-wz(k,l)
            r=sqrt(rx**2+ry**2+rz**2)
            lfield(1,i)=lfield(1,i)+cw(l)*gc(j)/r
            lfield(2,i)=lfield(2,i)+cw(l)*ec(j)/r
         end do
      end do
   end do

!Calculate polarization energy
   do ic=1,2
   do j=1,nmol*9
      elec0(j)=0.0
   end do
   do j=1,nmol*9
      do k=1,nmol*9
         mA(j,k)=0.0
      end do
   end do
   do j=1,nmol
      do k=1,21
         do l=1,3
            rx=chromx(k)-wx(j,l)
            ry=chromy(k)-wy(j,l)
            rz=chromz(k)-wz(j,l)
            r=sqrt(rx**2+ry**2+rz**2)
            if(ic==1) then
               charge=gc(k)
            else
               charge=ec(k)
            end if
            elec0(3*(3*(j-1)+l-1)+1)=elec0(3*(3*(j-1)+l-1)+1)+charge*rx/r**3
            elec0(3*(3*(j-1)+l-1)+2)=elec0(3*(3*(j-1)+l-1)+2)+charge*ry/r**3
            elec0(3*(3*(j-1)+l-1)+3)=elec0(3*(3*(j-1)+l-1)+3)+charge*rz/r**3
         end do
      end do
      elecQB=elec0
      do k=1,nmol
         do m=1,3
            do l=1,3
               if(k==j.and.m==l) cycle
               rx=wx(k,m)-wx(j,l)
               ry=wy(k,m)-wy(j,l)
               rz=wz(k,m)-wz(j,l)
               r=sqrt(rx**2+ry**2+rz**2)
               elec0(3*(3*(j-1)+l-1)+1)=elec0(3*(3*(j-1)+l-1)+1)+cw(m)*rx/r**3
               elec0(3*(3*(j-1)+l-1)+2)=elec0(3*(3*(j-1)+l-1)+2)+cw(m)*ry/r**3
               elec0(3*(3*(j-1)+l-1)+3)=elec0(3*(3*(j-1)+l-1)+3)+cw(m)*rz/r**3  
            end do
         end do
      end do
   end do

   do j=1,nmol
      do k=1,nmol
         do l=1,3
            do m=1,3
               if(j==k.and.l==m) then
                  mA(3*(3*(j-1)+l-1)+1,3*(3*(k-1)+m-1)+1)=1.0/pow(l)
                  mA(3*(3*(j-1)+l-1)+2,3*(3*(k-1)+m-1)+2)=1.0/pow(l)
                  mA(3*(3*(j-1)+l-1)+3,3*(3*(k-1)+m-1)+3)=1.0/pow(l)
               else
                  rx=wx(k,m)-wx(j,l)
                  ry=wy(k,m)-wy(j,l)
                  rz=wz(k,m)-wz(j,l)
                  r=sqrt(rx**2+ry**2+rz**2)
!                  sij=1.662*(pow(l)*pow(m))**(1.0/6)
!                  if(r<sij) then
!                     vij=r/sij
!                  else
                     vij=1
!                  end if
                  mA(3*(3*(j-1)+l-1)+1,3*(3*(k-1)+m-1)+1)=(4*vij**3-3*vij**4)/r**3-3.0*vij**4*rx*rx/r**5
                  mA(3*(3*(j-1)+l-1)+1,3*(3*(k-1)+m-1)+2)=-3.0*vij**4*rx*ry/r**5
                  mA(3*(3*(j-1)+l-1)+1,3*(3*(k-1)+m-1)+3)=-3.0*vij**4*rx*rz/r**5
                  mA(3*(3*(j-1)+l-1)+2,3*(3*(k-1)+m-1)+1)=-3.0*vij**4*ry*rx/r**5
                  mA(3*(3*(j-1)+l-1)+2,3*(3*(k-1)+m-1)+2)=(4*vij**3-3*vij**4)/r**3-3.0*vij**4*ry*ry/r**5
                  mA(3*(3*(j-1)+l-1)+2,3*(3*(k-1)+m-1)+3)=-3.0*vij**4*ry*rz/r**5
                  mA(3*(3*(j-1)+l-1)+3,3*(3*(k-1)+m-1)+1)=-3.0*vij**4*rz*rx/r**5
                  mA(3*(3*(j-1)+l-1)+3,3*(3*(k-1)+m-1)+2)=-3.0*vij**4*rz*ry/r**5
                  mA(3*(3*(j-1)+l-1)+3,3*(3*(k-1)+m-1)+3)=(4*vij**3-3*vij**4)/r**3-3.0*vij**4*rz*rz/r**5
               end if 
            end do
         end do
      end do
   end do

   invmA=inverse(mA)
   indp=matmul(elec0,transpose(invmA))

   lfield(2+ic,i)=lfield(2+ic,i)-0.5*dot_product(indp,elecQB)
   end do

   lfield(5,i)=1.0d0   
end do
 

do i=1,5
   call MPI_GATHER(lfield(i,:),per,MPI_DOUBLE,field(i,:),per,MPI_DOUBLE,MASTER,MPI_COMM_WORLD,ierr)
end do

!call MPI_FINALIZE(ierr)

do i=1,5
   do j=1,5
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
inv=matmul(field,transpose(field))
!field=field/nconfig

do i=1,5
   XY(j)=0.0
end do

!do i=1,nconfig
!   do j=1,64
!      XY(j)=XY(j)+field(i,j)*exe(i)
!   end do
!end do
XY=matmul(field,exe)
!XY=XY/nconfig

if(taskid.eq.MASTER) then
do i=1,5
   write(30,"(5F20.10)") inv(i,:)
end do

write(50,"(5F20.10)") XY

do i=1,nconfig
   write(60,"(5F20.10)") field(:,i)
end do
end if

call MPI_FINALIZE(ierr)

close(30)
close(40)
close(50)
close(60)
close(70)

end program
   
function inverse(A) result(Ainv)
  real*8, dimension(:,:), intent(in) :: A
  real*8, dimension(size(A,1),size(A,2)) :: Ainv

  real*8, dimension(size(A,1)) :: work  ! work array for LAPACK
  integer, dimension(size(A,1)) :: ipiv   ! pivot indices
  integer :: n, info

  ! External procedures defined in LAPACK
  external DGETRF
  external DGETRI

  ! Store A in Ainv to prevent it from being overwritten by LAPACK
  Ainv = A
  n = size(A,1)

  ! DGETRF computes an LU factorization of a general M-by-N matrix A
  ! using partial pivoting with row interchanges.
  call DGETRF(n, n, Ainv, n, ipiv, info)

  if (info /= 0) then
     stop 'Matrix is numerically singular!'
  end if

  ! DGETRI computes the inverse of a matrix using the LU factorization
  ! computed by DGETRF.
  call DGETRI(n, Ainv, n, ipiv, work, n, info)

  if (info /= 0) then
     stop 'Matrix inversion failed!'
  end if
end function inverse

