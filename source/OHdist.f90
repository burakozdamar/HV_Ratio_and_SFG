program OHdist

implicit none
integer*8 :: i, j, step, Natom, numL3, numL2, numL1, numL0, Nstep
real*8 :: xo, yo, zo, dist, xdiff, ydiff, zdiff, coord=0, xdiffintra, ydiffintra, zdiffintra, L3intercoord=0, L2intercoord=0
real*8 :: x1, y1, z1, x2, y2, z2, rohinter, rohintra, coseno, L3L0coord=0, L3L1coord=0, L3L2coord=0, L3L3coord=0
real*8 :: L0u, L1u, L2u, L3u, L2L0coord=0, L2L1coord=0, L2L2coord=0, L2L3coord=0, L1intercoord=0, L0intercoord=0
real*8 :: L1L0coord=0, L1L1coord=0, L1L2coord=0, L1L3coord=0, L0L0coord=0, L0L1coord=0, L0L2coord=0, L0L3coord=0
real*8 :: aL0u, aL1u, aL2u, aL3u, aL2L0coord=0, aL2L1coord=0, aL2L2coord=0, aL2L3coord=0, aL1intercoord=0, aL0intercoord=0
real*8 :: aL1L0coord=0, aL1L1coord=0, aL1L2coord=0, aL1L3coord=0, aL0L0coord=0, aL0L1coord=0, aL0L2coord=0, aL0L3coord=0
real*8 :: dL0u, dL1u, dL2u, dL3u, dL2L0coord=0, dL2L1coord=0, dL2L2coord=0, dL2L3coord=0, dL1intercoord=0, dL0intercoord=0
real*8 :: dL1L0coord=0, dL1L1coord=0, dL1L2coord=0, dL1L3coord=0, dL0L0coord=0, dL0L1coord=0, dL0L2coord=0, dL0L3coord=0
real*8 :: aL3intercoord=0, aL2intercoord=0, aL3L0coord=0, aL3L1coord=0, aL3L2coord=0, aL3L3coord=0
real*8 :: dL3intercoord=0, dL2intercoord=0, dL3L0coord=0, dL3L1coord=0, dL3L2coord=0, dL3L3coord=0
real*8 :: nL1L1=0, nL3L3=0, rL1L1, avrL1L1=0, rL3L3, avrL3L3=0
real*8 :: totL3=0, totL2=0, totL1=0, totL0=0, a, b, c
real*8, parameter :: cut=3.2
real*8, dimension(:) , allocatable :: xL3, yL3, zL3, xL2, yL2, zL2, xL1, yL1, zL1, xL0, yL0, zL0 
real*8, dimension(:) , allocatable :: xh1L3, yh1L3, zh1L3, xh1L2, yh1L2, zh1L2, xh1L1, yh1L1, zh1L1, xh1L0, yh1L0, zh1L0
real*8, dimension(:) , allocatable :: xh2L3, yh2L3, zh2L3, xh2L2, yh2L2, zh2L2, xh2L1, yh2L1, zh2L1, xh2L0, yh2L0, zh2L0
character(5) :: atom, hyd
integer*8 NO

!working flag for getting info
character(256)  txt

!###########################################
! STARTS OF THE PROGRAM
!##########################################


!----------------------------------------
! GET INFO FROM BOXDATA
!-----------------------------------------
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


open (1, FILE="up_down_water_layers.xyz")
open (13, FILE="L3_OHdist.dat")
open (12, FILE="L2_OHdist.dat")
open (11, FILE="L1_OHdist.dat")
open (10, FILE="L0_OHdist.dat")


do step = 1, Nstep
print*, "step", step


                                    ! reading atomic coordinate
 read(1,*) Natom
 read(1,*)
 do i=1, NO
  read(1,*) atom, xo, yo, zo
  read(1,*) hyd, x1, y1, z1
  read(1,*) hyd, x2, y2, z2
  if (atom == "dOL3" .or. atom == "uOL3") then
   xdiff=x1-xo
   ydiff=y1-yo
   zdiff=z1-zo
   dist=SQRT((xdiff)**2+(ydiff)**2+(zdiff)**2)
   write(13,*) dist
   xdiff=x2-xo
   ydiff=y2-yo
   zdiff=z2-zo
   dist=SQRT((xdiff)**2+(ydiff)**2+(zdiff)**2)
   write(13,*) dist

  else if (atom == "dOL2" .or. atom == "uOL2"  ) then
   xdiff=x1-xo
   ydiff=y1-yo
   zdiff=z1-zo
   dist=SQRT((xdiff)**2+(ydiff)**2+(zdiff)**2)
   write(12,*) dist
   xdiff=x2-xo
   ydiff=y2-yo
   zdiff=z2-zo
   dist=SQRT((xdiff)**2+(ydiff)**2+(zdiff)**2)
   write(12,*) dist

  else if (atom == "dOL1" .or. atom == "uOL1"  ) then
   xdiff=x1-xo
   ydiff=y1-yo
   zdiff=z1-zo
   dist=SQRT((xdiff)**2+(ydiff)**2+(zdiff)**2)
   write(11,*) dist
   xdiff=x2-xo
   ydiff=y2-yo
   zdiff=z2-zo
   dist=SQRT((xdiff)**2+(ydiff)**2+(zdiff)**2)
   write(11,*) dist

  else if (atom == "dOL0" .or. atom == "uOL0" ) then
   xdiff=x1-xo
   ydiff=y1-yo
   zdiff=z1-zo
   dist=SQRT((xdiff)**2+(ydiff)**2+(zdiff)**2)
   write(10,*) dist
   xdiff=x2-xo
   ydiff=y2-yo
   zdiff=z2-zo
   dist=SQRT((xdiff)**2+(ydiff)**2+(zdiff)**2)
   write(10,*) dist

  end if
 end do
end do

stop
 end program
