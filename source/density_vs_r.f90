program density_vs_r2

implicit none

real*8, parameter ::  tol=0.7, bin=5 !bin=how many points in 1 Angstrom
real*8 a, b, c
integer*8 xgrid, ygrid, Npoint, t, NstepTOT, Istart, Nskip1, Nskip2, stepp, NstepReal
integer*8, parameter ::  max=12, min=-2
integer*8 NO, io, io2
integer*8 :: Natom, Nstep, i, j, k, step, m, nop
real*8 :: xo, yo, zo, u, xdiff, ydiff, zdiff, urmin, drmin, uzdiff, dzdiff, r, pippo
character(3) :: atom
real*8, dimension (:), allocatable :: x, y, z, dens, dist, ddens, udens
real*8, dimension (:,:), allocatable :: ulx, uly, ulz, dlx, dly, dlz
real*8 fake1,fake2

!working flag for getting info
character(256)  txt


!-------------------------------------

open(11,file='BOXDATA',form='formatted')
!     Inizialize the Flag
      txt = ' ' 
!     Loop untill the Flag
      do while (trim(txt).NE.'$BOX-DIMENTIONS')
         read(11,*) txt
      enddo
      !print*, txt
      read(11,*) a
      read(11,*) b
      read(11,*) c
close(11)
open(11,file='BOXDATA',form='formatted')
!     Inizialize the Flag
      txt = ' ' 
!     Loop untill the Flag
      do while (trim(txt).NE.'$NO')
         read(11,*) txt
      enddo
      !print*, txt
      read(11,*) NO
close(11)
open(11,file='BOXDATA',form='formatted')
!     Inizialize the Flag
      txt = ' ' 
!     Loop untill the Flag
      do while (trim(txt).NE.'$NSTEP')
         read(11,*) txt
      enddo
      !print*, txt
      read(11,*) NstepTOT
close(11)
open(11,file='BOXDATA',form='formatted')
      txt = ' '
      do while (trim(txt).NE.'$ISTART')
          read(11,*) txt
      enddo
      !print*, txt
      read(11,*) Istart
close(11)

Nstep=NstepTOT-Istart+1


open(12,file='grid_interface',form='formatted')
read(12,*) !Skipe line
read(12,*) !Skipe line
read(12,*) xgrid, fake1, fake2
read(12,*) ygrid, fake1, fake2
close(12)


open (5, FILE="interface.xyz")
  read(5,*) Npoint
close(5)
!-------------------------------


nop=(max-min)*bin+1

allocate (x(NO))
allocate (y(NO))
allocate (z(NO))
allocate (dens(nop))
allocate (udens(nop))
allocate (ddens(nop))
allocate (dist(nop))
allocate (ulx(xgrid,ygrid))
allocate (uly(xgrid,ygrid))
allocate (ulz(xgrid,ygrid))
allocate (dlx(xgrid,ygrid))
allocate (dly(xgrid,ygrid))
allocate (dlz(xgrid,ygrid))

dens=0
udens=0
ddens=0

open (1, FILE="pos_rebuilt.xyz")
open (2, FILE="interface.xyz")
open (3, FILE="density_vs_r")
open (4, FILE="up_density_vs_r")
open (5, FILE="down_density_vs_r")

! skip pos_rebuilt.xyz
Nskip1=(Istart-1)*(NO*3+2)
do t = 1, Nskip1
    read(1,*)
enddo

Nskip2=(Istart-1)*(Npoint+2)
do t=1, Nskip2/10
    read(2,*)
enddo

do step = 1, Nstep/10
 !print*, step
                                     ! reading atomic coordinate for each step
 read(2,*,iostat=io2)
    if (io2 /= 0) then 
        exit 
    endif
 read(2,*)
 do i=1,xgrid
  do j=1,ygrid
   read(2,*) atom, ulx(i,j), uly(i,j), ulz(i,j)
   read(2,*) atom, dlx(i,j), dly(i,j), dlz(i,j)
  end do
 end do
 !write(*,*) step

 do stepp = 1, 10
 j=1
 read(1,*,iostat=io) Natom
 if (io /= 0) then
     exit 
 endif
 read(1,*)
 do while (j<= NO)
  read(1,*) atom, xo, yo, zo
  if (atom == "O") then
   x(j)=xo
   y(j)=yo
   z(j)=zo
   read(1,*)
   read(1,*)
   j = j+1
  end if
 end do
! do i=1, 6930
!  read(1,*)
! end do

                                     ! reading grid points
! do i=1, 28818
!  read(2,*)
! end do
                                     ! searching the nearest grid point for each molecule
 do m=1,NO
  urmin=100
  drmin=100
  do i=1,xgrid
   do j=1,ygrid
                                              ! up
    xdiff=ulx(i,j)-x(m)
    ydiff=uly(i,j)-y(m)
    zdiff=ulz(i,j)-z(m)
                                     ! pbc
    if ( xdiff > a/2) then
     xdiff=xdiff-a
    else if ( xdiff < -a/2) then
     xdiff=xdiff+a
    end if

    if ( ydiff > b/2) then
     ydiff=ydiff-b
    else if ( ydiff < -b/2) then
     ydiff=ydiff+b
    end if

                                     ! end pbc
    r=SQRT(xdiff**2+ydiff**2+zdiff**2)
    if (r<urmin) then
     urmin=r
     uzdiff=zdiff
    end if
                                               ! down
    xdiff=dlx(i,j)-x(m)
    ydiff=dly(i,j)-y(m)
    zdiff=dlz(i,j)-z(m)
                                     ! pbc
    if ( xdiff > a/2) then
     xdiff=xdiff-a
    else if ( xdiff < -a/2) then
     xdiff=xdiff+a
    end if

    if ( ydiff > b/2) then
     ydiff=ydiff-b
    else if ( ydiff < -b/2) then
     ydiff=ydiff+b
    end if
                                     ! end pbc
    r=SQRT(xdiff**2+ydiff**2+zdiff**2)
    if (r<drmin) then
     drmin=r
     dzdiff=zdiff
    end if
   end do
  end do
  if (uzdiff<0) then
   urmin=-urmin
  end if
  if (dzdiff>0) then
   drmin=-drmin
  end if
!  print*, urmin, drmin
                                     ! evaluating density as function rmin
  do i=1,nop
   u=i-1
   dist(i)=u/bin+min
   if (urmin>=dist(i)-tol .and. urmin<dist(i)+tol) then
    dens(i) = dens(i)+1
    udens(i) = udens(i)+1
   end if
   if (drmin>=dist(i)-tol .and. drmin<dist(i)+tol) then
    dens(i) = dens(i)+1
    ddens(i) = ddens(i)+1
   end if
  end do


 end do
 end do ! end stepp (10)
end do



deallocate (x)
deallocate (y)
deallocate (z)
deallocate (ulx)
deallocate (uly)
deallocate (ulz)
deallocate (dlx)
deallocate (dly)
deallocate (dlz)
NstepReal = (step-2)*10
do i=1,nop
 dens(i)=dens(i)/(NstepReal*tol*a*b*4*0.035)
 udens(i)=udens(i)/(NstepReal*tol*a*b*2*0.035)
 ddens(i)=ddens(i)/(NstepReal*tol*a*b*2*0.035)
 write(3,*) dist(i)          ,          dens(i)
 write(4,*) dist(i)          ,          udens(i)
 write(5,*) dist(i)          ,          ddens(i)
end do

deallocate (dens)
deallocate (dist)
deallocate (udens)
deallocate (ddens)

stop 
end program
