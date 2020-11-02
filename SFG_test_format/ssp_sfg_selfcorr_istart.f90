program ssp_sfg_selfcorr_istart

! Evry 13/06/2017


implicit none

integer*8 :: t, Nstep,NstepTOT, NO, Natom,istart, i, j, k, m
integer*8 :: ntime, ntt0, nt0, Neff, Nskip,Nskip1,Nskip2
integer*8 :: Nsolid, Nions
real*8, dimension(3) :: dMdRz
real*8, dimension(3,3) :: dAdRz
real*8, dimension(:), allocatable :: avMZ, avAXX, flag, corr, anormal
real*8, dimension(:,:), allocatable :: MZ, AXX
real*8, dimension(:,:,:,:,:), allocatable :: D
real*8, dimension(:,:,:), allocatable :: vz
real*8 :: dx, dy, dz, dx1, dy1, dz1, dx2, dy2, dz2, dr, scalar
real*8 :: xo,  yo,  zo, xh1, yh1, zh1, xh2, yh2, zh2, a, b, dt
real*8, parameter :: conv= 21.8769125400 ! convert vel: bohr/a.u.TIME --> Ang/fs
character(3) :: atom

integer narg

character(len=256) input 
character(len=256) output 
character(len=256) txt 

!---------------------------------------

      
narg=command_argument_count()

if (narg .eq. 2) then
   call get_command_argument(1,input)
   call get_command_argument(2,output)

else
   write (*,*) 
   write (*,*) 'Usage:./ssp_sfg_selfcorr input output'
   write (*,*) '----------------------------------------'
   write (*,*) 
   stop
endif

!-------------------------------------------------------------

open(11,file='BOXDATA',form='formatted')
!     Inizialize the Flag
      txt = ' ' 
!     Loop untill the Flag
      do while (trim(txt).NE.'$BOX-DIMENTIONS')
         read(11,*) txt
      enddo
      print*, txt
      read(11,*) a
      read(11,*) b
close(11)
open(11,file='BOXDATA',form='formatted')
!     Inizialize the Flag
      txt = ' ' 
!     Loop untill the Flag
      do while (trim(txt).NE.'$NO')
         read(11,*) txt
      enddo
      print*, txt
      read(11,*) NO
close(11)
open(11,file='BOXDATA',form='formatted')
!     Inizialize the Flag
      txt = ' ' 
!     Loop untill the Flag
      do while (trim(txt).NE.'$NSTEP')
         read(11,*) txt
      enddo
      print*, txt
      read(11,*) NstepTOT
close(11)
open(11,file='BOXDATA',form='formatted')
!     Inizialize the Flag
      txt = ' '
!     Loop untill the Flag
      do while (trim(txt).NE.'$ISTART')
         read(11,*) txt
      enddo
      print*, txt
      read(11,*) Istart
close(11)

open(11,file='BOXDATA',form='formatted')
!     Inizialize the Flag
      txt = ' ' 
!     Loop untill the Flag
      do while (trim(txt).NE.'$NSOLID')
         read(11,*) txt
      enddo
      print*, txt
      read(11,*) Nsolid
close(11)

open(11,file='BOXDATA',form='formatted')
!     Inizialize the Flag
      txt = ' ' 
!     Loop untill the Flag
      do while (trim(txt).NE.'$NIONS')
         read(11,*) txt
      enddo
      print*, txt
      read(11,*) Nions
close(11)

open(11,file='BOXDATA',form='formatted')
!     Inizialize the Flag
      txt = ' ' 
!     Loop untill the Flag
      do while (trim(txt).NE.'$DT')
         read(11,*) txt
      enddo
      print*, txt
      read(11,*) dt
close(11)
Nskip=Nsolid+Nions  
Nstep=NstepTOT-Istart+1

write(*,*) '-------------END BOXDATA READING--------'





! Input
open (1, FILE="pos_rebuilt.xyz")
open (2, FILE="vel.xyz")
open (10,file=input)

!Output
open (8, FILE=output)
open (99,file="OH_ref")
open (98,file="z_vel_values")
open (89,file="tot_dMZ_dt")
open (88,file="tot_dAXX_dt")

! Skip the Istart steps
!-------------------------

!  pos_rebuilt.xyz
   Nskip1=(Istart-1)*(NO*3+2)
   do t=1, Nskip1
      read(1,*)
   enddo

!  vel_water.xyz
   Nskip2=(Istart-1)*(NO*3+Nions+Nsolid+2)
   do t=1, Nskip2
      read(2,*)
   enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! values from Remi article for : dM(x,y,z)/dRz and dA(x,y,z)/dRz
! where x,y,z are coordinates in the oh_ref
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
dMdRz(1)   =-0.1500000     ! = dM(x)/dRz
dMdRz(2)   =-0.d0          ! = dM(y)/dRz
dMdRz(3)   = 2.1000000     ! = dM(z)/dRz
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
dAdRz(1,1) = 0.4000000     ! = dA(x,x)/dRz
dAdRz(2,2) = 0.5300000     ! = dA(y,y)/dRz
dAdRz(3,3) = 1.5600000     ! = dA(z,z)/dRz

dAdRz(1,2) = 0.d0          ! = dA(x,y)/dRz
dAdRz(2,1) = dAdRz(1,2)    ! = dA(y,x)/dRz

dAdRz(1,3) = 0.0200000     ! = dA(x,z)/dRz
dAdRz(3,1) = dAdRz(1,3)    ! = dA(z,x)/dRz

dAdRz(2,3) = 0.d0          ! = dA(y,z)/dRz
dAdRz(3,2) = dAdRz(2,3)    ! = dA(z,y)/dRz
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! reading the binder file (to select the mol of the interface)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
allocate(flag(NO))
do i=1,NO
 read(10,*) j, flag(i)
end do
close(10)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! reading atom pos and calcule D matrix for each atom at each step
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!  oh_ref  !!!!!!
! it is for each oh_bond of each mol
! z = oh_bond direction
! y = perpendicular to z out of the mol plane
! x = perpendicular to z in the mol plane
!!!!!! D matrix !!!!!!
! D(mol(NO), OHbond(2), step(Nstep), oh_ref(3), cartesian_ref(3) 
allocate(D(NO,2,Nstep,3,3))
!!!!!!!!!!!!!!!!!!!!!!
print*, "evaluating D matrix"
do t=1,Nstep
 print*, "t",t
 read(1,*) Natom
 read(1,*)
 do i=1,NO
  read(1,*) atom, xo,  yo,  zo
  read(1,*) atom, xh1, yh1, zh1
  read(1,*) atom, xh2, yh2, zh2
  !!!!!!!!!!! z(1) !!!!!!!!!!
  dx1=xh1-xo
  dy1=yh1-yo
  dz1=zh1-zo
  dr=dsqrt(dx1**2+dy1**2+dz1**2)
  D(i,1,t,3,1)=dx1/dr
  D(i,1,t,3,2)=dy1/dr
  D(i,1,t,3,3)=dz1/dr
  !!!!!!!!!!! z(2) !!!!!!!!!!
  dx2=xh2-xo
  dy2=yh2-yo
  dz2=zh2-zo 
  dr=dsqrt(dx2**2+dy2**2+dz2**2)
  D(i,2,t,3,1)=dx2/dr
  D(i,2,t,3,2)=dy2/dr
  D(i,2,t,3,3)=dz2/dr
  !!!!!!!!!!! x(1) !!!!!!!!!!
  scalar=dx2*D(i,1,t,3,1)+dy2*D(i,1,t,3,2)+dz2*D(i,1,t,3,3)
  dx=scalar*D(i,1,t,3,1)-dx2
  dy=scalar*D(i,1,t,3,2)-dy2
  dz=scalar*D(i,1,t,3,3)-dz2
  dr=dsqrt(dx**2+dy**2+dz**2)
  D(i,1,t,1,1)=dx/dr
  D(i,1,t,1,2)=dy/dr
  D(i,1,t,1,3)=dz/dr
  !!!!!!!!!!! x(2) !!!!!!!!!!
  scalar=dx1*D(i,2,t,3,1)+dy1*D(i,2,t,3,2)+dz1*D(i,2,t,3,3)
  dx=scalar*D(i,2,t,3,1)-dx1
  dy=scalar*D(i,2,t,3,2)-dy1
  dz=scalar*D(i,2,t,3,3)-dz1
  dr=dsqrt(dx**2+dy**2+dz**2)
  D(i,2,t,1,1)=dx/dr
  D(i,2,t,1,2)=dy/dr
  D(i,2,t,1,3)=dz/dr
  !!!!!!!!!!! y(1) !!!!!!!!!!
  D(i,1,t,2,1)=D(i,1,t,3,2)*D(i,1,t,1,3)-D(i,1,t,3,3)*D(i,1,t,1,2)
  D(i,1,t,2,2)=D(i,1,t,3,3)*D(i,1,t,1,1)-D(i,1,t,3,1)*D(i,1,t,1,3)
  D(i,1,t,2,3)=D(i,1,t,3,1)*D(i,1,t,1,2)-D(i,1,t,3,2)*D(i,1,t,1,1)
  !!!!!!!!!!! y(2) !!!!!!!!!!
  D(i,2,t,2,1)=D(i,2,t,3,2)*D(i,2,t,1,3)-D(i,2,t,3,3)*D(i,2,t,1,2)
  D(i,2,t,2,2)=D(i,2,t,3,3)*D(i,2,t,1,1)-D(i,2,t,3,1)*D(i,2,t,1,3)
  D(i,2,t,2,3)=D(i,2,t,3,1)*D(i,2,t,1,2)-D(i,2,t,3,2)*D(i,2,t,1,1)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!
  write(99,*) i
  write(99,*) D(i,1,t,1,1), D(i,1,t,1,2), D(i,1,t,1,3)
  write(99,*) D(i,1,t,2,1), D(i,1,t,2,2), D(i,1,t,2,3)
  write(99,*) D(i,1,t,3,1), D(i,1,t,3,2), D(i,1,t,3,3)
  write(99,*) "**************************************"
  write(99,*) D(i,2,t,1,1), D(i,2,t,1,2), D(i,2,t,1,3)
  write(99,*) D(i,2,t,2,1), D(i,2,t,2,2), D(i,2,t,2,3)
  write(99,*) D(i,2,t,3,1), D(i,2,t,3,2), D(i,2,t,3,3)
  write(99,*) "**************************************"
  !!!!!!!!!!!!!!!!!!!!!!!!!!!
 end do
end do
close(1)
close(99)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! calculating v(z) (with z of the oh_ref) for each mol
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
print*, "calculating v(z) for each mol"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! v(z) is the projection of v(h)-v(o) on the direction of the oh_bond
allocate(vz(NO,2,Nstep))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
do t=1,Nstep
 print*, "t",t
 read(2,*) Natom
 read(2,*)
 do i=1,NO
  read(2,*) atom, xo,  yo,  zo
  read(2,*) atom, xh1, yh1, zh1
  read(2,*) atom, xh2, yh2, zh2
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  dx=xh1-xo
  dy=yh1-yo
  dz=zh1-zo
  vz(i,1,t)=dx*D(i,1,t,3,1)+dy*D(i,1,t,3,2)+dz*D(i,1,t,3,3)
  vz(i,1,t)= vz(i,1,t)*conv
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  dx=xh2-xo
  dy=yh2-yo
  dz=zh2-zo
  vz(i,2,t)=dx*D(i,2,t,3,1)+dy*D(i,2,t,3,2)+dz*D(i,2,t,3,3)
  vz(i,2,t)= vz(i,2,t)*conv
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  write(98,*) i
  write(98,*) vz(i,1,t)
  write(98,*) vz(i,1,t)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 end do
 do i=1,Nskip
  read(2,*)
 end do
end do
close(2)
close(98)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       ! calculate M(Z) (total dip) (z=cartesian coordinate) 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
allocate(MZ(NO,Nstep))
allocate(avMz(NO))
MZ=0.d0
avMz=0.d0
do t=1,Nstep 
 do m=1,NO
  if (flag(m)==1) then
   do i=1,2     ! do on N° oh_bond (=2) X mol
    do j=1,3     ! do on x,y,z (oh_ref)
     MZ(m,t)=MZ(m,t)+( D(m,i,t,j,3)*dMdRz(j)*vz(m,i,t) )
    end do
   end do
  end if
 end do
end do

do m=1,NO
 do t=1,Nstep
  avMz(m)=avMz(m)+MZ(m,t)
 end do
 avMz(m)=avMz(m)/dble(Nstep)
 write(89,*) avMz(m)
end do
close(89)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       ! calculate A(ZZ) (total pol) (x,y=cartesian coordinate)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
allocate(AXX(NO,Nstep))
allocate(avAXX(Nstep))
AXX=0.d0
avAXX=0.d0
do t=1,Nstep
 do m=1,NO
  if (flag(m)==1) then
   do i=1,2     ! do on N° oh_bond (=2) X mol
    do j=1,3     ! do on x,y,z (oh_ref)
     do k=1,3     ! do on x,y,z (oh_ref)
      AXX(m,t)=AXX(m,t)+( D(m,i,t,j,1)*D(m,i,t,k,1)*dAdRz(j,k)*vz(m,i,t) )
     end do       ! do on k
    end do       ! do on j
   end do
  end if
 end do
end do

do m=1,NO
 do t=1,Nstep
  avAXX(m)=avAXX(m)+AXX(m,t)
 end do
 avAXX(m)=avAXX(m)/dble(Nstep)
 write(88,*) avAXX(m)
end do
close(88)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

deallocate(D)
deallocate(vz)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! correlation function M(0)(Z)*A(t)(XX or YY)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
allocate(corr(0:Nstep))
allocate(anormal(0:Nstep))
corr=0.d0
anormal=0.d0
Neff=0
do m=1,NO
 print*, "mol", m
 if (flag(m)==1) then
  Neff=Neff+1
  do nt0=1,Nstep
   do ntt0=nt0,Nstep
    ntime = ntt0 - nt0     !time t for correlation function
    corr(ntime)= corr(ntime)+ (AXX(m,ntt0)-avAXX(m))*(MZ(m,nt0)-avMZ(m))
    anormal(ntime) = anormal(ntime) + 1   
   end do
  end do
 end if
end do

do i=0,Nstep-1
 corr(i)=corr(i)/dble(anormal(i))
 corr(i)=corr(i)*(0.208*1.6*dt*Neff)/(a*b*10)
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! the result is times 10 at the minus 23
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 write(8,*)i,corr(i)
end do
deallocate(corr)
deallocate(anormal)
deallocate(AXX)
deallocate(avAXX)
deallocate(MZ)
deallocate(avMZ)
close(8)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

stop
end program
