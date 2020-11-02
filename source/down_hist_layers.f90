program  down_hist_layers

! 01/03/2017

implicit none

real*8, parameter :: bin=5
real*8  intL0, intL1, intL2, intL3
real*8, parameter :: ulimit=0.33, dlimit=-0.33
integer*8 NO, xgrid, ygrid
integer*8 :: Natom, Nstep, i, j, k, step, m
real*8 :: xo, yo, zo, xdiff, ydiff, zdiff, urmin, drmin, uzdiff, dzdiff
real*8 :: xbis, ybis, zbis, bis, xoh1, xoh2, yoh1, yoh2, zoh1, zoh2, coseno, r, u
real*8 :: uL0=0, uL1=0, uL2=0, uL3=0, dL0=0, dL1=0, dL2=0, dL3=0, pL0=0, pL1=0, pL2=0, pL3=0 
real*8 :: nL0=0, nL1=0, nL2=0, nL3=0, a, b, c
character(3) :: atom
real*8, dimension (:), allocatable :: x, y, z, xh1, yh1, zh1, xh2, yh2, zh2
real*8, dimension (:,:), allocatable :: ulx, uly, ulz, dlx, dly, dlz
real*8 fake1, fake2

!working flag for getting info
character(256)  txt

!###########################################
! STARTS OF THE PROGRAM
!##########################################


open(11,file='BOXDATA',form='formatted')
!     Inizialize the Flag
      txt = ' '
!     
      do while (trim(txt).NE.'$BOX-DIMENTIONS')
         read(11,*) txt
      enddo
      print*, txt
      read(11,*) a
      read(11,*) b
      read(11,*) c
close(11)

open(11,file='BOXDATA',form='formatted')
!     Inizialize the Flag
      txt = ' '
!     Loop untill the Flag
      do while (trim(txt).NE.'$NSTEP')
         read(11,*) txt
      enddo
      print*, txt
      read(11,*) Nstep
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
      do while (trim(txt).NE.'$LAYERS-LIMITS')
         read(11,*) txt
      enddo
      read(11,*) intL0
      read(11,*) intL1
      read(11,*) intL2
      read(11,*) intL3
close(11)



write(*,*) '-------------END BOXDATA READING--------'

open(12,file='grid_interface',form='formatted')
read(12,*) !Skipe line
read(12,*) !Skipe line
read(12,*) xgrid, fake1, fake2
read(12,*) ygrid, fake1, fake2
close(12)






open (1, FILE="pos_rebuilt.xyz")
open (2, FILE="interface.xyz")
open (3, FILE="histogram_L0")
open (4, FILE="histogram_L1")
open (5, FILE="histogram_L2")
open (7, FILE="histogram_L3")

!---------------------------
! Allocation
!-----------------------------


allocate (x(NO))
allocate (y(NO))
allocate (z(NO))
allocate (xh1(NO))
allocate (yh1(NO))
allocate (zh1(NO))
allocate (xh2(NO))
allocate (yh2(NO))
allocate (zh2(NO))
allocate (ulx(xgrid,ygrid))
allocate (uly(xgrid,ygrid))
allocate (ulz(xgrid,ygrid))
allocate (dlx(xgrid,ygrid))
allocate (dly(xgrid,ygrid))
allocate (dlz(xgrid,ygrid))

Nstep = Nstep/10
do step = 1, Nstep
 print*, "step", step
                                    ! reading atomic coordinate for each step
 j=1
 read(1,*) Natom
 read(1,*)
 do while (j<= NO)
  read(1,*) atom, xo, yo, zo
  if (atom == "O") then
   x(j)=xo
   y(j)=yo
   z(j)=zo
  read(1,*) atom, xh1(j), yh1(j), zh1(j)
  read(1,*) atom, xh2(j), yh2(j), zh2(j)
  j = j+1
  end if
 end do
                                     ! reading grid points for each step
 read(2,*)
 read(2,*)
 do i=1,xgrid
  do j=1,ygrid
   read(2,*) atom, ulx(i,j), uly(i,j), ulz(i,j)
   read(2,*) atom, dlx(i,j), dly(i,j), dlz(i,j)
  end do
 end do
                               ! search the nearest grid point for each molecule + evaluate r
 do m=1,NO
  urmin=100
  drmin=100
  do i=1,xgrid
   do j=1,ygrid
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
                                                ! calculate orientation
    xoh1 = xh1(m) - x(m)
    yoh1 = yh1(m) - y(m)
    zoh1 = zh1(m) - z(m)
    xoh2 = xh2(m) - x(m)
    yoh2 = yh2(m) - y(m)
    zoh2 = zh2(m) - z(m)

    xbis=xoh1+xoh2
    ybis=yoh1+yoh2
    zbis=zoh1+zoh2
    bis=SQRT(xbis**2+ybis**2+zbis**2)

                                              ! down
                           ! select only the molecules of interest (between min & max)
                                                ! calculate orientation
    coseno=-zbis/bis

                                                ! evaluating orientation proprieties
    if (drmin <= intL0 ) then
     nL0=nL0+1
     if (coseno>ulimit) then
      uL0=uL0+1
     else if (coseno<dlimit) then
      dL0=dL0+1
     else
      pL0=pL0+1
     end if

    else if (drmin > intL0 .and. drmin < intL1 ) then
     nL1=nL1+1
     if (coseno>ulimit) then
      uL1=uL1+1
     else if (coseno<dlimit) then
      dL1=dL1+1
     else
      pL1=pL1+1
     end if

    else if (drmin >= intL1 .and. drmin <= intL2 ) then
     nL2=nL2+1
     if (coseno>ulimit) then
      uL2=uL2+1
     else if (coseno<dlimit) then
      dL2=dL2+1
     else
      pL2=pL2+1
     end if

    else if (drmin >= intL2 .and. drmin < intL3 ) then
     nL3=nL3+1
     if (coseno>ulimit) then
      uL3=uL3+1
     else if (coseno<dlimit) then
      dL3=dL3+1
     else
      pL3=pL3+1
     end if
    end if

 end do
end do

deallocate (x)
deallocate (y)
deallocate (z)
deallocate (xh1)
deallocate (yh1)
deallocate (zh1)
deallocate (xh2)
deallocate (yh2)
deallocate (zh2)
deallocate (ulx)
deallocate (uly)
deallocate (ulz)
deallocate (dlx)
deallocate (dly)
deallocate (dlz)

uL0=uL0/nL0
uL1=uL1/nL1
uL2=uL2/nL2
uL3=uL3/nL3

pL0=pL0/nL0
pL1=pL1/nL1
pL2=pL2/nL2
pL3=pL3/nL3

dL0=dL0/nL0
dL1=dL1/nL1
dL2=dL2/nL2
dL3=dL3/nL3

write (3,*) "0"    ,    uL0
write (3,*) "1"    ,    pL0
write (3,*) "2"    ,    dL0

write (4,*) "0"    ,    uL1
write (4,*) "1"    ,    pL1
write (4,*) "2"    ,    dL1

write (5,*) "0"    ,    uL2
write (5,*) "1"    ,    pL2
write (5,*) "2"    ,    dL2

write (7,*) "0"    ,    uL3
write (7,*) "1"    ,    pL3
write (7,*) "2"    ,    dL3

stop
end program
