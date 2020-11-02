program interface

implicit none

real*8, parameter :: pi=3.14159265359, E=2.4, int=1, ext=2.5, tollerance=0.004
real*8, dimension (:), allocatable :: x, y, z, lx, ly, lz
real*8, dimension (:), allocatable :: xh1, yh1, zh1, xh2, yh2, zh2
real*8 :: l1, l2, l3, a, b, c, xo, yo, zo, xdiff, ydiff, zdiff, r, p, lzdiff, pdiff
real*8 :: lznear, rmin, l3max, l3min, dmax, dmin, l2max, l2min, l1max, l1min, targup, targdown
integer*8 :: step, Nstep, d, Natom, i, j, k, NO, m, s, st, u, v, w, Nskip
character(3) :: atom
integer*8 :: nlx,nly,nlz,io
real*8 :: dlx, dly, dlz, fx, fy, fz

!working flag for getting info
character(256)  txt

!------------------------------------------
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
      read(11,*) c
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
      read(11,*) Nstep
close(11)


print*, '*********END BOXDATA READING******'





open (1, FILE="pos_rebuilt.xyz")
open (2, FILE="interface.xyz")
open (11, FILE="grid_interface")

allocate (x(NO))
allocate (y(NO))
allocate (z(NO))
allocate (xh1(NO))
allocate (yh1(NO))
allocate (zh1(NO))
allocate (xh2(NO))
allocate (yh2(NO))
allocate (zh2(NO))
allocate (lx(9200))
allocate (ly(9200))
allocate (lz(9200))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!   GRID DEFINITION     !!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

dlx=0.5d0
dly=0.5d0
dlz=0.25d0

nlx=NINT(a/dlx)
nly=NINT(b/dly)
nlz=NINT(c/dlz)

fx=-((a/2)-(dlx/2))
fy=-((b/2)-(dly/2))
fz=-((c/2)-(dlz/2))


write(11,*) "========================="
write(11,*) " Npoint | first | last"
write(11,*) nlx, fx, (nlx-1)*dlx+fx
write(11,*) nly, fy, (nly-1)*dly+fy
write(11,*) nlz, fz, (nlz-1)*dlz+fz
write(11,*) "========================="
close(11)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



do step = 1, Nstep/10
! print*, step
 d=0
 u=0
 v=0
                                      ! reading input file
 j=0
 read(1,*,iostat=io) Natom
 if (io /= 0 ) then
     stop
 endif
 read(1,*)
 do i=1,Natom
  read(1,*) atom, xo, yo, zo
  if (atom == "O") then
   j=j+1
   x(j)=xo
   y(j)=yo
   z(j)=zo
  end if
 end do
 if (j/=NO) then
  print*, "***** ERROR READING INPUT FILE *****"
 end if
 Nskip=9*(Natom+2)
 do i=1, Nskip
  read(1,*)
 end do
 
 print*, step

                                      ! evaluating p
 do i=1,nlx
  l1=dlx*(i-1)+fx
  do j=1,nly
   l2=dly*(j-1)+fy
   l3max = -100
   l3min = 100
   dmax=-100
   dmin=100
   targup=100        ! targup= reference distance from target p value (0.016) for interface up point
   targdown=100      ! targdown= reference distance from target p value (0.016) for interface down point
   do k=1,nlz
    l3=dlz*(k-1)+fz
    p=0  
    do m=1, NO
     xdiff=x(m)-l1
     ydiff=y(m)-l2
     zdiff=z(m)-l3    
                                               ! correction for pbc
     if (xdiff < 0) then
      xdiff = -xdiff
     else 
     end if 
     if (xdiff > a/2) then
      xdiff = a - xdiff
     else 
     end if
     if (ydiff < 0) then
      ydiff = -ydiff
     else
     end if 
     if (ydiff > b/2) then
      ydiff = b - ydiff
     else 
     end if
     if (zdiff < 0) then
      zdiff = -zdiff
     else
     end if 
     if (zdiff > c/2) then
      zdiff = c - zdiff
     else 
     end if
                                               ! end correction for pbc
     r = SQRT(xdiff**2+ydiff**2+zdiff**2)
     if (r<=3*E) then
      p = p + exp(-r**2/(2*E**2))/((2*pi*E**2)**1.5)
     end if    
    end do
                                       ! evaluating if each grid point is part of the interface
    pdiff = 0.016 - p
    if (pdiff < 0) then
     pdiff = -pdiff
    end if    
    if ( pdiff < tollerance) then
     if (l3 > 0  .and. pdiff < targup) then
      dmax = 1
      l3max = l3
      l2max = l2
      l1max = l1
      targup=pdiff
     else if (l3 < 0 .and. pdiff < targdown) then
      dmin = 1
      l3min = l3
      l2min = l2
      l1min = l1
      targdown=pdiff
     end if
    end if
   end do
   if ( dmax > 0) then
    d = d+1
    lx(d) = l1max
    ly(d) = l2max
    lz(d) = l3max
   end if
   if ( dmin > 0) then
    d = d+1
    lx(d) = l1min
    ly(d) = l2min
    lz(d) = l3min
   end if
  end do
 end do

 !do i = 1, 10
  !s=10*(step-1)+i
  write(2,*) d
  write(2,*) "step =", step
  do j = 1, d
   write(2,FMT="(A , F14.6 , F14.6 , F14.6)") "P", lx(j), ly(j), lz(j)
  end do
 !end do
end do


deallocate (x)
deallocate (y)
deallocate (z)
deallocate (lx)
deallocate (ly)
deallocate (lz)
deallocate (xh1)
deallocate (yh1)
deallocate (zh1)
deallocate (xh2)
deallocate (yh2)
deallocate (zh2)

stop
end program
