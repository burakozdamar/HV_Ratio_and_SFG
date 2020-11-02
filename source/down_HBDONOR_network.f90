program down_HB_network

implicit none
integer*8 :: i, j, step, Natom, numL3, numL2, numL1, numL0, Nstep, numoL2
real*8 :: xo, yo, zo, dist, xdiff, ydiff, zdiff, coord=0, xdiffintra, ydiffintra, zdiffintra, L3intercoord=0, L2intercoord=0
real*8 :: x1, y1, z1, x2, y2, z2, rohinter, rohintra, coseno, L3L0coord=0, L3L1coord=0, Lcutcoord=0, L3L3coord=0
real*8 :: L0u, L1u, L2u, L3u, L2L0coord=0, L2L1coord=0, L2L2coord=0, L2L3coord=0, L1intercoord=0, L0intercoord=0
real*8 :: L1L0coord=0, L1L1coord=0, L1L2coord=0, L1L3coord=0, L0L0coord=0, L0L1coord=0, L0L2coord=0, L0L3coord=0
real*8 :: aL0u, aL1u, aL2u, aL3u, aL2L0coord=0, aL2L1coord=0, aL2L2coord=0, aL2L3coord=0, aL1intercoord=0, aL0intercoord=0
real*8 :: aL1L0coord=0, aL1L1coord=0, aL1L2coord=0, aL1L3coord=0, aL0L0coord=0, aL0L1coord=0, aL0L2coord=0, aL0L3coord=0
real*8 :: dL0u, dL1u, dL2u, dL3u, dL2L0coord=0, dL2L1coord=0, dL2L2coord=0, dL2L3coord=0, dL1intercoord=0, dL0intercoord=0
real*8 :: dL1L0coord=0, dL1L1coord=0, dL1L2coord=0, dL1L3coord=0, dL0L0coord=0, dL0L1coord=0, dL0L2coord=0, dL0L3coord=0
real*8 :: aL3intercoord=0, aL2intercoord=0, aL3L0coord=0, aL3L1coord=0, aLcutcoord=0, aL3L3coord=0
real*8 :: dL3intercoord=0, dL2intercoord=0, dL3L0coord=0, dL3L1coord=0, dLcutcoord=0, dL3L3coord=0
real*8 :: nL1L1=0, nL3L3=0, rL1L1, avrL1L1=0, rL3L3, avrL3L3=0
real*8 :: totL3=0, totL2=0, totL1=0, totL0=0, a, b, c, zcos
integer*8 NO
real*8, parameter :: cut=3.2
real*8, dimension(:) , allocatable :: xL3, yL3, zL3, xL2, yL2, zL2, xL1, yL1, zL1, xL0, yL0, zL0 
real*8, dimension(:) , allocatable :: xh1L3, yh1L3, zh1L3, xh1L2, yh1L2, zh1L2, xh1L1, yh1L1, zh1L1, xh1L0, yh1L0, zh1L0
real*8, dimension(:) , allocatable :: xh2L3, yh2L3, zh2L3, xh2L2, yh2L2, zh2L2, xh2L1, yh2L1, zh2L1, xh2L0, yh2L0, zh2L0
real*8, dimension(:) , allocatable :: oxL2, oyL2, ozL2
character(5) :: atom, hyd
integer*8 :: io1
integer*8 :: Npoint, t, NstepTOT, Istart 

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
      read(11,*) NO
close(11)
open(11,file='BOXDATA',form='formatted')
!     Inizialize the Flag
      txt = ' ' 
!     Loop untill the Flag
      do while (trim(txt).NE.'$NSTEP')
         read(11,*) txt
      enddo
      read(11,*) NstepTOT
close(11)
open(11,file='BOXDATA',form='formatted')
      txt = ' '
      do while (trim(txt).NE.'$ISTART')
          read(11,*) txt
      enddo
      read(11,*) Istart
close(11)

Nstep=NstepTOT-Istart+1
!-------------------------------
!-------------------------------
! Allocate
!--------------------------------------

allocate (xL0(NO))
allocate (yL0(NO))
allocate (zL0(NO))
allocate (xL1(NO))
allocate (yL1(NO))
allocate (zL1(NO))
allocate (xL2(NO))
allocate (yL2(NO))
allocate (zL2(NO))
allocate (xL3(NO))
allocate (yL3(NO))
allocate (zL3(NO))
allocate (xh1L0(NO))
allocate (yh1L0(NO))
allocate (zh1L0(NO))
allocate (xh2L0(NO))
allocate (yh2L0(NO))
allocate (zh2L0(NO))
allocate (xh1L1(NO))
allocate (yh1L1(NO))
allocate (zh1L1(NO))
allocate (xh2L1(NO))
allocate (yh2L1(NO))
allocate (zh2L1(NO))
allocate (xh1L2(NO))
allocate (yh1L2(NO))
allocate (zh1L2(NO))
allocate (xh2L2(NO))
allocate (yh2L2(NO))
allocate (zh2L2(NO))
allocate (xh1L3(NO))
allocate (yh1L3(NO))
allocate (zh1L3(NO))
allocate (xh2L3(NO))
allocate (yh2L3(NO))
allocate (zh2L3(NO))
allocate (oxL2(NO))
allocate (oyL2(NO))
allocate (ozL2(NO))


!

open (1, FILE="up_down_water_layers.xyz")

open (7, FILE="HBDONOR_network_tot_dL3")
open (8, FILE="HBDONOR_network_inter_dL3")
open (9, FILE="HBDONOR_network_intra_dL3")

open (17, FILE="HBDONOR_network_tot_dL2")
open (18, FILE="HBDONOR_network_inter_dL2")
open (19, FILE="HBDONOR_network_intra_dL2")

open (27, FILE="HBDONOR_network_tot_dL1")
open (28, FILE="HBDONOR_network_inter_dL1")
open (29, FILE="HBDONOR_network_intra_dL1")

open (37, FILE="HBDONOR_network_tot_dL0")
open (38, FILE="HBDONOR_network_inter_dL0")
open (39, FILE="HBDONOR_network_intra_dL0")



do step = 1, Nstep
 numL3=0
 numL2=0
 numL1=0
 numL0=0
 numoL2=0

                                    ! reading atomic coordinate
 read(1,*,iostat=io1) Natom
    if (io1 /= 0) then
        stop
    endif
 read(1,*)
 do i=1, NO
  read(1,*) atom, xo, yo, zo
  read(1,*) hyd, x1, y1, z1
  read(1,*) hyd, x2, y2, z2
  if (atom == "dOL3" .or. atom == "uOL3") then
   numL3=numL3+1
   xL3(numL3)=xo
   yL3(numL3)=yo
   zL3(numL3)=zo
   xh1L3(numL3)=x1
   yh1L3(numL3)=y1
   zh1L3(numL3)=z1
   xh2L3(numL3)=x2
   yh2L3(numL3)=y2
   zh2L3(numL3)=z2

  else if (atom == "dOL2"  ) then
   numL2=numL2+1
   xL2(numL2)=xo
   yL2(numL2)=yo
   zL2(numL2)=zo
   xh1L2(numL2)=x1
   yh1L2(numL2)=y1
   zh1L2(numL2)=z1
   xh2L2(numL2)=x2
   yh2L2(numL2)=y2
   zh2L2(numL2)=z2


  else if (atom == "uOL2"  ) then
   numoL2=numoL2+1
   oxL2(numoL2)=xo
   oyL2(numoL2)=yo
   ozL2(numoL2)=zo




  else if (atom == "dOL1"  ) then
   numL1=numL1+1
   xL1(numL1)=xo
   yL1(numL1)=yo
   zL1(numL1)=zo
   xh1L1(numL1)=x1
   yh1L1(numL1)=y1
   zh1L1(numL1)=z1
   xh2L1(numL1)=x2
   yh2L1(numL1)=y2
   zh2L1(numL1)=z2

  else if (atom == "dOL0"  ) then
   numL0=numL0+1
   xL0(numL0)=xo
   yL0(numL0)=yo
   zL0(numL0)=zo
   xh1L0(numL0)=x1
   yh1L0(numL0)=y1
   zh1L0(numL0)=z1
   xh2L0(numL0)=x2
   yh2L0(numL0)=y2
   zh2L0(numL0)=z2

  end if
 end do
  

                                                                         ! evaluating number of Hbond

                                                                     ! Layer L3
 do i=1, numL3
                                                                 ! with L3
  do j=1,numL3
   rL3L3=0
   xdiff=xL3(i)-xL3(j)
   ydiff=yL3(i)-yL3(j)
   zdiff=zL3(i)-zL3(j)
                                               ! correction for pbc
   if (xdiff > a/2) then
    xdiff = xdiff - a
   else if (xdiff < -a/2) then
    xdiff = xdiff + a
   end if
   if (ydiff > b/2) then
    ydiff = ydiff - b
   else if (ydiff < -b/2) then
    ydiff = ydiff + b
   end if
                                               ! end correction for pbc
   dist=SQRT((xdiff)**2+(ydiff)**2+(zdiff)**2)
   zcos= -zdiff/dist
   if (dist < cut .and. dist>0.1) then
                                                                              ! H1-O
    xdiff=xL3(j)-xh1L3(i)
    ydiff=yL3(j)-yh1L3(i)
    zdiff=zL3(j)-zh1L3(i)
    xdiffintra=xL3(i)-xh1L3(i)
    ydiffintra=yL3(i)-yh1L3(i)
    zdiffintra=zL3(i)-zh1L3(i)
                                               ! correction for pbc
    if (xdiff > a/2) then
     xdiff = xdiff - a
    else if (xdiff < -a/2) then
     xdiff = xdiff + a
    end if
    if (ydiff > b/2) then
     ydiff = ydiff - b
    else if (ydiff < -b/2) then
     ydiff = ydiff + b
    end if
    
    if (xdiffintra > a/2) then
     xdiffintra = xdiffintra - a
    else if (xdiffintra < -a/2) then
     xdiffintra = xdiffintra + a
    end if
    if (ydiffintra > b/2) then
     ydiffintra = ydiffintra - b
    else if (ydiffintra < -b/2) then
     ydiffintra = ydiffintra + b
    end if
                                               ! end correction for pbc
    rohinter=SQRT((xdiff)**2+(ydiff)**2+(zdiff)**2)
    rohintra=SQRT((xdiffintra)**2+(ydiffintra)**2+(zdiffintra)**2)
    coseno=(xdiff*xdiffintra+ydiff*ydiffintra+zdiff*zdiffintra)/(rohinter*rohintra)
    if (coseno < -0.76604) then
     zcos=-zcos
     write(9,*) dist, zcos
     write(7,*) dist, zcos
    end if
                                                                              ! H2-O
    xdiff=xL3(j)-xh2L3(i)
    ydiff=yL3(j)-yh2L3(i)
    zdiff=zL3(j)-zh2L3(i)
    xdiffintra=xL3(i)-xh2L3(i)
    ydiffintra=yL3(i)-yh2L3(i)
    zdiffintra=zL3(i)-zh2L3(i)
                                               ! correction for pbc
    if (xdiff > a/2) then
     xdiff = xdiff - a
    else if (xdiff < -a/2) then
     xdiff = xdiff + a
    end if
    if (ydiff > b/2) then
     ydiff = ydiff - b
    else if (ydiff < -b/2) then
     ydiff = ydiff + b
    end if

    if (xdiffintra > a/2) then
     xdiffintra = xdiffintra - a
    else if (xdiffintra < -a/2) then
     xdiffintra = xdiffintra + a
    end if
    if (ydiffintra > b/2) then
     ydiffintra = ydiffintra - b
    else if (ydiffintra < -b/2) then
     ydiffintra = ydiffintra + b
    end if
                                               ! end correction for pbc
    rohinter=SQRT((xdiff)**2+(ydiff)**2+(zdiff)**2)
    rohintra=SQRT((xdiffintra)**2+(ydiffintra)**2+(zdiffintra)**2)
    coseno=(xdiff*xdiffintra+ydiff*ydiffintra+zdiff*zdiffintra)/(rohinter*rohintra)
    if (coseno < -0.76604) then
     zcos=-zcos
     write(9,*) dist, zcos
     write(7,*) dist, zcos
    end if

   end if
  end do

                                                                 ! with other L2
  do j=1,numoL2
   rL3L3=0
   xdiff=xL3(i)-oxL2(j)
   ydiff=yL3(i)-oyL2(j)
   zdiff=zL3(i)-ozL2(j)
                                               ! correction for pbc
   if (xdiff > a/2) then
    xdiff = xdiff - a
   else if (xdiff < -a/2) then
    xdiff = xdiff + a
   end if
   if (ydiff > b/2) then
    ydiff = ydiff - b
   else if (ydiff < -b/2) then
    ydiff = ydiff + b
   end if
                                               ! end correction for pbc
   dist=SQRT((xdiff)**2+(ydiff)**2+(zdiff)**2)
   zcos= -zdiff/dist
   if (dist < cut .and. dist>0.1) then
                                                                              ! H1-O
    xdiff=oxL2(j)-xh1L3(i)
    ydiff=oyL2(j)-yh1L3(i)
    zdiff=ozL2(j)-zh1L3(i)
    xdiffintra=xL3(i)-xh1L3(i)
    ydiffintra=yL3(i)-yh1L3(i)
    zdiffintra=zL3(i)-zh1L3(i)
                                               ! correction for pbc
    if (xdiff > a/2) then
     xdiff = xdiff - a
    else if (xdiff < -a/2) then
     xdiff = xdiff + a
    end if
    if (ydiff > b/2) then
     ydiff = ydiff - b
    else if (ydiff < -b/2) then
     ydiff = ydiff + b
    end if

    if (xdiffintra > a/2) then
     xdiffintra = xdiffintra - a
    else if (xdiffintra < -a/2) then
     xdiffintra = xdiffintra + a
    end if
    if (ydiffintra > b/2) then
     ydiffintra = ydiffintra - b
    else if (ydiffintra < -b/2) then
     ydiffintra = ydiffintra + b
    end if
                                               ! end correction for pbc
    rohinter=SQRT((xdiff)**2+(ydiff)**2+(zdiff)**2)
    rohintra=SQRT((xdiffintra)**2+(ydiffintra)**2+(zdiffintra)**2)
    coseno=(xdiff*xdiffintra+ydiff*ydiffintra+zdiff*zdiffintra)/(rohinter*rohintra)
    if (coseno < -0.76604) then
     zcos=-zcos
     write(9,*) dist, zcos
     write(7,*) dist, zcos
    end if
                                                                              ! H2-O
    xdiff=oxL2(j)-xh2L3(i)
    ydiff=oyL2(j)-yh2L3(i)
    zdiff=ozL2(j)-zh2L3(i)
    xdiffintra=xL3(i)-xh2L3(i)
    ydiffintra=yL3(i)-yh2L3(i)
    zdiffintra=zL3(i)-zh2L3(i)
                                               ! correction for pbc
    if (xdiff > a/2) then
     xdiff = xdiff - a
    else if (xdiff < -a/2) then
     xdiff = xdiff + a
    end if
    if (ydiff > b/2) then
     ydiff = ydiff - b
    else if (ydiff < -b/2) then
     ydiff = ydiff + b
    end if

    if (xdiffintra > a/2) then
     xdiffintra = xdiffintra - a
    else if (xdiffintra < -a/2) then
     xdiffintra = xdiffintra + a
    end if
    if (ydiffintra > b/2) then
     ydiffintra = ydiffintra - b
    else if (ydiffintra < -b/2) then
     ydiffintra = ydiffintra + b
    end if
                                               ! end correction for pbc
    rohinter=SQRT((xdiff)**2+(ydiff)**2+(zdiff)**2)
    rohintra=SQRT((xdiffintra)**2+(ydiffintra)**2+(zdiffintra)**2)
    coseno=(xdiff*xdiffintra+ydiff*ydiffintra+zdiff*zdiffintra)/(rohinter*rohintra)
    if (coseno < -0.76604) then
     zcos=-zcos
     write(9,*) dist, zcos
     write(7,*) dist, zcos
    end if

   end if
  end do


                                                             ! with L2
  do j=1,numL2
   xdiff=xL3(i)-xL2(j)
   ydiff=yL3(i)-yL2(j)
   zdiff=zL3(i)-zL2(j)
                                               ! correction for pbc
   if (xdiff > a/2) then
    xdiff = xdiff - a
   else if (xdiff < -a/2) then
    xdiff = xdiff + a
   end if
   if (ydiff > b/2) then
    ydiff = ydiff - b
   else if (ydiff < -b/2) then
    ydiff = ydiff + b
   end if
                                               ! end correction for pbc
   dist=SQRT((xdiff)**2+(ydiff)**2+(zdiff)**2)
   zcos= -zdiff/dist
   if (dist < cut .and. dist>0.1) then
                                                                              ! O-H1
    xdiff=xL3(i)-xh1L2(j)
    ydiff=yL3(i)-yh1L2(j)
    zdiff=zL3(i)-zh1L2(j)
    xdiffintra=xL2(j)-xh1L2(j)
    ydiffintra=yL2(j)-yh1L2(j)
    zdiffintra=zL2(j)-zh1L2(j)
                                               ! correction for pbc
    if (xdiff > a/2) then
     xdiff = xdiff - a
    else if (xdiff < -a/2) then
     xdiff = xdiff + a
    end if
    if (ydiff > b/2) then
     ydiff = ydiff - b
    else if (ydiff < -b/2) then
     ydiff = ydiff + b
    end if

    if (xdiffintra > a/2) then
     xdiffintra = xdiffintra - a
    else if (xdiffintra < -a/2) then
     xdiffintra = xdiffintra + a
    end if
    if (ydiffintra > b/2) then
     ydiffintra = ydiffintra - b
    else if (ydiffintra < -b/2) then
     ydiffintra = ydiffintra + b
    end if
                                               ! end correction for pbc
    rohinter=SQRT((xdiff)**2+(ydiff)**2+(zdiff)**2)
    rohintra=SQRT((xdiffintra)**2+(ydiffintra)**2+(zdiffintra)**2)
    coseno=(xdiff*xdiffintra+ydiff*ydiffintra+zdiff*zdiffintra)/(rohinter*rohintra)
    if (coseno < -0.76604) then
     write(18,*) dist, zcos
     write(17,*) dist, zcos      
    end if    
                                                                              ! O-H2
    xdiff=xL3(i)-xh2L2(j)
    ydiff=yL3(i)-yh2L2(j)
    zdiff=zL3(i)-zh2L2(j)
    xdiffintra=xL2(j)-xh2L2(j)
    ydiffintra=yL2(j)-yh2L2(j)
    zdiffintra=zL2(j)-zh2L2(j)
                                               ! correction for pbc
    if (xdiff > a/2) then
     xdiff = xdiff - a
    else if (xdiff < -a/2) then
     xdiff = xdiff + a
    end if
    if (ydiff > b/2) then
     ydiff = ydiff - b
    else if (ydiff < -b/2) then
     ydiff = ydiff + b
    end if

    if (xdiffintra > a/2) then
     xdiffintra = xdiffintra - a
    else if (xdiffintra < -a/2) then
     xdiffintra = xdiffintra + a
    end if
    if (ydiffintra > b/2) then
     ydiffintra = ydiffintra - b
    else if (ydiffintra < -b/2) then
     ydiffintra = ydiffintra + b
    end if
                                               ! end correction for pbc
    rohinter=SQRT((xdiff)**2+(ydiff)**2+(zdiff)**2)
    rohintra=SQRT((xdiffintra)**2+(ydiffintra)**2+(zdiffintra)**2)
    coseno=(xdiff*xdiffintra+ydiff*ydiffintra+zdiff*zdiffintra)/(rohinter*rohintra)
    if (coseno < -0.76604) then
     write(18,*) dist, zcos
     write(17,*) dist, zcos
    end if
                                                                              ! H1-O
    xdiff=xL2(j)-xh1L3(i)
    ydiff=yL2(j)-yh1L3(i)
    zdiff=zL2(j)-zh1L3(i)
    xdiffintra=xL3(i)-xh1L3(i)
    ydiffintra=yL3(i)-yh1L3(i)
    zdiffintra=zL3(i)-zh1L3(i)
                                               ! correction for pbc
    if (xdiff > a/2) then
     xdiff = xdiff - a
    else if (xdiff < -a/2) then
     xdiff = xdiff + a
    end if
    if (ydiff > b/2) then
     ydiff = ydiff - b
    else if (ydiff < -b/2) then
     ydiff = ydiff + b
    end if

    if (xdiffintra > a/2) then
     xdiffintra = xdiffintra - a
    else if (xdiffintra < -a/2) then
     xdiffintra = xdiffintra + a
    end if
    if (ydiffintra > b/2) then
     ydiffintra = ydiffintra - b
    else if (ydiffintra < -b/2) then
     ydiffintra = ydiffintra + b
    end if
                                               ! end correction for pbc
    rohinter=SQRT((xdiff)**2+(ydiff)**2+(zdiff)**2)
    rohintra=SQRT((xdiffintra)**2+(ydiffintra)**2+(zdiffintra)**2)
    coseno=(xdiff*xdiffintra+ydiff*ydiffintra+zdiff*zdiffintra)/(rohinter*rohintra)
    if (coseno < -0.76604) then
     zcos=-zcos
     write(8,*) dist, zcos
     write(7,*) dist, zcos
    end if
                                                                              ! H2-O
    xdiff=xL2(j)-xh2L3(i)
    ydiff=yL2(j)-yh2L3(i)
    zdiff=zL2(j)-zh2L3(i)
    xdiffintra=xL3(i)-xh2L3(i)
    ydiffintra=yL3(i)-yh2L3(i)
    zdiffintra=zL3(i)-zh2L3(i)
                                               ! correction for pbc
    if (xdiff > a/2) then
     xdiff = xdiff - a
    else if (xdiff < -a/2) then
     xdiff = xdiff + a
    end if
    if (ydiff > b/2) then
     ydiff = ydiff - b
    else if (ydiff < -b/2) then
     ydiff = ydiff + b
    end if

    if (xdiffintra > a/2) then
     xdiffintra = xdiffintra - a
    else if (xdiffintra < -a/2) then
     xdiffintra = xdiffintra + a
    end if
    if (ydiffintra > b/2) then
     ydiffintra = ydiffintra - b
    else if (ydiffintra < -b/2) then
     ydiffintra = ydiffintra + b
    end if
                                               ! end correction for pbc
    rohinter=SQRT((xdiff)**2+(ydiff)**2+(zdiff)**2)
    rohintra=SQRT((xdiffintra)**2+(ydiffintra)**2+(zdiffintra)**2)
    coseno=(xdiff*xdiffintra+ydiff*ydiffintra+zdiff*zdiffintra)/(rohinter*rohintra)
    if (coseno < -0.76604) then
     zcos=-zcos
     write(8,*) dist, zcos
     write(7,*) dist, zcos
    end if
   end if
  end do


                                                             ! with L1
  do j=1,numL1
   xdiff=xL3(i)-xL1(j)
   ydiff=yL3(i)-yL1(j)
   zdiff=zL3(i)-zL1(j)
                                               ! correction for pbc
   if (xdiff > a/2) then
    xdiff = xdiff - a
   else if (xdiff < -a/2) then
    xdiff = xdiff + a
   end if
   if (ydiff > b/2) then
    ydiff = ydiff - b
   else if (ydiff < -b/2) then
    ydiff = ydiff + b
   end if
                                               ! end correction for pbc
   dist=SQRT((xdiff)**2+(ydiff)**2+(zdiff)**2)
   zcos= -zdiff/dist
   if (dist < cut .and. dist>0.1) then
                                                                              ! O-H1
    xdiff=xL3(i)-xh1L1(j)
    ydiff=yL3(i)-yh1L1(j)
    zdiff=zL3(i)-zh1L1(j)
    xdiffintra=xL1(j)-xh1L1(j)
    ydiffintra=yL1(j)-yh1L1(j)
    zdiffintra=zL1(j)-zh1L1(j)
                                               ! correction for pbc
    if (xdiff > a/2) then
     xdiff = xdiff - a
    else if (xdiff < -a/2) then
     xdiff = xdiff + a
    end if
    if (ydiff > b/2) then
     ydiff = ydiff - b
    else if (ydiff < -b/2) then
     ydiff = ydiff + b
    end if

    if (xdiffintra > a/2) then
     xdiffintra = xdiffintra - a
    else if (xdiffintra < -a/2) then
     xdiffintra = xdiffintra + a
    end if
    if (ydiffintra > b/2) then
     ydiffintra = ydiffintra - b
    else if (ydiffintra < -b/2) then
     ydiffintra = ydiffintra + b
    end if
                                               ! end correction for pbc
    rohinter=SQRT((xdiff)**2+(ydiff)**2+(zdiff)**2)
    rohintra=SQRT((xdiffintra)**2+(ydiffintra)**2+(zdiffintra)**2)
    coseno=(xdiff*xdiffintra+ydiff*ydiffintra+zdiff*zdiffintra)/(rohinter*rohintra)
    if (coseno < -0.76604) then
     write(28,*) dist, zcos
     write(27,*) dist, zcos
    end if
                                                                              ! O-H2
    xdiff=xL3(i)-xh2L1(j)
    ydiff=yL3(i)-yh2L1(j)
    zdiff=zL3(i)-zh2L1(j)
    xdiffintra=xL1(j)-xh2L1(j)
    ydiffintra=yL1(j)-yh2L1(j)
    zdiffintra=zL1(j)-zh2L1(j)
                                               ! correction for pbc
    if (xdiff > a/2) then
     xdiff = xdiff - a
    else if (xdiff < -a/2) then
     xdiff = xdiff + a
    end if
    if (ydiff > b/2) then
     ydiff = ydiff - b
    else if (ydiff < -b/2) then
     ydiff = ydiff + b
    end if

    if (xdiffintra > a/2) then
     xdiffintra = xdiffintra - a
    else if (xdiffintra < -a/2) then
     xdiffintra = xdiffintra + a
    end if
    if (ydiffintra > b/2) then
     ydiffintra = ydiffintra - b
    else if (ydiffintra < -b/2) then
     ydiffintra = ydiffintra + b
    end if
                                               ! end correction for pbc
    rohinter=SQRT((xdiff)**2+(ydiff)**2+(zdiff)**2)
    rohintra=SQRT((xdiffintra)**2+(ydiffintra)**2+(zdiffintra)**2)
    coseno=(xdiff*xdiffintra+ydiff*ydiffintra+zdiff*zdiffintra)/(rohinter*rohintra)
    if (coseno < -0.76604) then
     write(28,*) dist, zcos
     write(27,*) dist, zcos
    end if
                                                                              ! H1-O
    xdiff=xL1(j)-xh1L3(i)
    ydiff=yL1(j)-yh1L3(i)
    zdiff=zL1(j)-zh1L3(i)
    xdiffintra=xL3(i)-xh1L3(i)
    ydiffintra=yL3(i)-yh1L3(i)
    zdiffintra=zL3(i)-zh1L3(i)
                                               ! correction for pbc
    if (xdiff > a/2) then
     xdiff = xdiff - a
    else if (xdiff < -a/2) then
     xdiff = xdiff + a
    end if
    if (ydiff > b/2) then
     ydiff = ydiff - b
    else if (ydiff < -b/2) then
     ydiff = ydiff + b
    end if

    if (xdiffintra > a/2) then
     xdiffintra = xdiffintra - a
    else if (xdiffintra < -a/2) then
     xdiffintra = xdiffintra + a
    end if
    if (ydiffintra > b/2) then
     ydiffintra = ydiffintra - b
    else if (ydiffintra < -b/2) then
     ydiffintra = ydiffintra + b
    end if
                                               ! end correction for pbc
    rohinter=SQRT((xdiff)**2+(ydiff)**2+(zdiff)**2)
    rohintra=SQRT((xdiffintra)**2+(ydiffintra)**2+(zdiffintra)**2)
    coseno=(xdiff*xdiffintra+ydiff*ydiffintra+zdiff*zdiffintra)/(rohinter*rohintra)
    if (coseno < -0.76604) then
     zcos=-zcos
     write(8,*) dist, zcos
     write(7,*) dist, zcos
    end if
                                                                              ! H2-O
    xdiff=xL1(j)-xh2L3(i)
    ydiff=yL1(j)-yh2L3(i)
    zdiff=zL1(j)-zh2L3(i)
    xdiffintra=xL3(i)-xh2L3(i)
    ydiffintra=yL3(i)-yh2L3(i)
    zdiffintra=zL3(i)-zh2L3(i)
                                               ! correction for pbc
    if (xdiff > a/2) then
     xdiff = xdiff - a
    else if (xdiff < -a/2) then
     xdiff = xdiff + a
    end if
    if (ydiff > b/2) then
     ydiff = ydiff - b
    else if (ydiff < -b/2) then
     ydiff = ydiff + b
    end if

    if (xdiffintra > a/2) then
     xdiffintra = xdiffintra - a
    else if (xdiffintra < -a/2) then
     xdiffintra = xdiffintra + a
    end if
    if (ydiffintra > b/2) then
     ydiffintra = ydiffintra - b
    else if (ydiffintra < -b/2) then
     ydiffintra = ydiffintra + b
    end if
                                               ! end correction for pbc
    rohinter=SQRT((xdiff)**2+(ydiff)**2+(zdiff)**2)
    rohintra=SQRT((xdiffintra)**2+(ydiffintra)**2+(zdiffintra)**2)
    coseno=(xdiff*xdiffintra+ydiff*ydiffintra+zdiff*zdiffintra)/(rohinter*rohintra)
    if (coseno < -0.76604) then
     zcos=-zcos
     write(8,*) dist, zcos
     write(7,*) dist, zcos
    end if
   end if
  end do


                                                             ! with L0
  do j=1,numL0
   xdiff=xL3(i)-xL0(j)
   ydiff=yL3(i)-yL0(j)
   zdiff=zL3(i)-zL0(j)
                                               ! correction for pbc
   if (xdiff > a/2) then
    xdiff = xdiff - a
   else if (xdiff < -a/2) then
    xdiff = xdiff + a
   end if
   if (ydiff > b/2) then
    ydiff = ydiff - b
   else if (ydiff < -b/2) then
    ydiff = ydiff + b
   end if
                                               ! end correction for pbc
   dist=SQRT((xdiff)**2+(ydiff)**2+(zdiff)**2)
   zcos= -zdiff/dist
   if (dist < cut .and. dist>0.1) then
                                                                              ! O-H1
    xdiff=xL3(i)-xh1L0(j)
    ydiff=yL3(i)-yh1L0(j)
    zdiff=zL3(i)-zh1L0(j)
    xdiffintra=xL1(j)-xh1L0(j)
    ydiffintra=yL1(j)-yh1L0(j)
    zdiffintra=zL1(j)-zh1L0(j)
                                               ! correction for pbc
    if (xdiff > a/2) then
     xdiff = xdiff - a
    else if (xdiff < -a/2) then
     xdiff = xdiff + a
    end if
    if (ydiff > b/2) then
     ydiff = ydiff - b
    else if (ydiff < -b/2) then
     ydiff = ydiff + b
    end if

    if (xdiffintra > a/2) then
     xdiffintra = xdiffintra - a
    else if (xdiffintra < -a/2) then
     xdiffintra = xdiffintra + a
    end if
    if (ydiffintra > b/2) then
     ydiffintra = ydiffintra - b
    else if (ydiffintra < -b/2) then
     ydiffintra = ydiffintra + b
    end if
                                               ! end correction for pbc
    rohinter=SQRT((xdiff)**2+(ydiff)**2+(zdiff)**2)
    rohintra=SQRT((xdiffintra)**2+(ydiffintra)**2+(zdiffintra)**2)
    coseno=(xdiff*xdiffintra+ydiff*ydiffintra+zdiff*zdiffintra)/(rohinter*rohintra)
    if (coseno < -0.76604) then
     write(38,*) dist, zcos
     write(37,*) dist, zcos
    end if
                                                                              ! O-H2
    xdiff=xL3(i)-xh2L0(j)
    ydiff=yL3(i)-yh2L0(j)
    zdiff=zL3(i)-zh2L0(j)
    xdiffintra=xL1(j)-xh2L0(j)
    ydiffintra=yL1(j)-yh2L0(j)
    zdiffintra=zL1(j)-zh2L0(j)
                                               ! correction for pbc
    if (xdiff > a/2) then
     xdiff = xdiff - a
    else if (xdiff < -a/2) then
     xdiff = xdiff + a
    end if
    if (ydiff > b/2) then
     ydiff = ydiff - b
    else if (ydiff < -b/2) then
     ydiff = ydiff + b
    end if

    if (xdiffintra > a/2) then
     xdiffintra = xdiffintra - a
    else if (xdiffintra < -a/2) then
     xdiffintra = xdiffintra + a
    end if
    if (ydiffintra > b/2) then
     ydiffintra = ydiffintra - b
    else if (ydiffintra < -b/2) then
     ydiffintra = ydiffintra + b
    end if
                                               ! end correction for pbc
    rohinter=SQRT((xdiff)**2+(ydiff)**2+(zdiff)**2)
    rohintra=SQRT((xdiffintra)**2+(ydiffintra)**2+(zdiffintra)**2)
    coseno=(xdiff*xdiffintra+ydiff*ydiffintra+zdiff*zdiffintra)/(rohinter*rohintra)
    if (coseno < -0.76604) then
     write(38,*) dist, zcos
     write(37,*) dist, zcos
    end if
                                                                              ! H1-O
    xdiff=xL0(j)-xh1L3(i)
    ydiff=yL0(j)-yh1L3(i)
    zdiff=zL0(j)-zh1L3(i)
    xdiffintra=xL3(i)-xh1L3(i)
    ydiffintra=yL3(i)-yh1L3(i)
    zdiffintra=zL3(i)-zh1L3(i)
                                               ! correction for pbc
    if (xdiff > a/2) then
     xdiff = xdiff - a
    else if (xdiff < -a/2) then
     xdiff = xdiff + a
    end if
    if (ydiff > b/2) then
     ydiff = ydiff - b
    else if (ydiff < -b/2) then
     ydiff = ydiff + b
    end if

    if (xdiffintra > a/2) then
     xdiffintra = xdiffintra - a
    else if (xdiffintra < -a/2) then
     xdiffintra = xdiffintra + a
    end if
    if (ydiffintra > b/2) then
     ydiffintra = ydiffintra - b
    else if (ydiffintra < -b/2) then
     ydiffintra = ydiffintra + b
    end if
                                               ! end correction for pbc
    rohinter=SQRT((xdiff)**2+(ydiff)**2+(zdiff)**2)
    rohintra=SQRT((xdiffintra)**2+(ydiffintra)**2+(zdiffintra)**2)
    coseno=(xdiff*xdiffintra+ydiff*ydiffintra+zdiff*zdiffintra)/(rohinter*rohintra)
    if (coseno < -0.76604) then
     zcos=-zcos
     write(8,*) dist, zcos
     write(7,*) dist, zcos
    end if
                                                                              ! H2-O
    xdiff=xL0(j)-xh2L3(i)
    ydiff=yL0(j)-yh2L3(i)
    zdiff=zL0(j)-zh2L3(i)
    xdiffintra=xL3(i)-xh2L3(i)
    ydiffintra=yL3(i)-yh2L3(i)
    zdiffintra=zL3(i)-zh2L3(i)
                                               ! correction for pbc
    if (xdiff > a/2) then
     xdiff = xdiff - a
    else if (xdiff < -a/2) then
     xdiff = xdiff + a
    end if
    if (ydiff > b/2) then
     ydiff = ydiff - b
    else if (ydiff < -b/2) then
     ydiff = ydiff + b
    end if

    if (xdiffintra > a/2) then
     xdiffintra = xdiffintra - a
    else if (xdiffintra < -a/2) then
     xdiffintra = xdiffintra + a
    end if
    if (ydiffintra > b/2) then
     ydiffintra = ydiffintra - b
    else if (ydiffintra < -b/2) then
     ydiffintra = ydiffintra + b
    end if
                                               ! end correction for pbc
    rohinter=SQRT((xdiff)**2+(ydiff)**2+(zdiff)**2)
    rohintra=SQRT((xdiffintra)**2+(ydiffintra)**2+(zdiffintra)**2)
    coseno=(xdiff*xdiffintra+ydiff*ydiffintra+zdiff*zdiffintra)/(rohinter*rohintra)
    if (coseno < -0.76604) then
     zcos=-zcos
     write(8,*) dist, zcos
     write(7,*) dist, zcos
    end if
   end if
  end do
 
 end do



                                                                     ! Layer L2
 do i=1, numL2
                                                                 ! with L2
  do j=1,numL2
   xdiff=xL2(i)-xL2(j)
   ydiff=yL2(i)-yL2(j)
   zdiff=zL2(i)-zL2(j)
                                               ! correction for pbc
   if (xdiff > a/2) then
    xdiff = xdiff - a
   else if (xdiff < -a/2) then
    xdiff = xdiff + a
   end if
   if (ydiff > b/2) then
    ydiff = ydiff - b
   else if (ydiff < -b/2) then
    ydiff = ydiff + b
   end if
                                               ! end correction for pbc
   dist=SQRT((xdiff)**2+(ydiff)**2+(zdiff)**2)
   zcos= -zdiff/dist
   if (dist < cut .and. dist>0.1) then
                                                                              ! H1-O
    xdiff=xL2(j)-xh1L2(i)
    ydiff=yL2(j)-yh1L2(i)
    zdiff=zL2(j)-zh1L2(i)
    xdiffintra=xL2(i)-xh1L2(i)
    ydiffintra=yL2(i)-yh1L2(i)
    zdiffintra=zL2(i)-zh1L2(i)
                                               ! correction for pbc
    if (xdiff > a/2) then
     xdiff = xdiff - a
    else if (xdiff < -a/2) then
     xdiff = xdiff + a
    end if
    if (ydiff > b/2) then
     ydiff = ydiff - b
    else if (ydiff < -b/2) then
     ydiff = ydiff + b
    end if

    if (xdiffintra > a/2) then
     xdiffintra = xdiffintra - a
    else if (xdiffintra < -a/2) then
     xdiffintra = xdiffintra + a
    end if
    if (ydiffintra > b/2) then
     ydiffintra = ydiffintra - b
    else if (ydiffintra < -b/2) then
     ydiffintra = ydiffintra + b
    end if
                                               ! end correction for pbc
    rohinter=SQRT((xdiff)**2+(ydiff)**2+(zdiff)**2)
    rohintra=SQRT((xdiffintra)**2+(ydiffintra)**2+(zdiffintra)**2)
    coseno=(xdiff*xdiffintra+ydiff*ydiffintra+zdiff*zdiffintra)/(rohinter*rohintra)
    if (coseno < -0.76604) then
     zcos=-zcos
     write(19,*) dist, zcos
     write(17,*) dist, zcos
    end if
                                                                              ! H2-O
    xdiff=xL2(j)-xh2L2(i)
    ydiff=yL2(j)-yh2L2(i)
    zdiff=zL2(j)-zh2L2(i)
    xdiffintra=xL2(i)-xh2L2(i)
    ydiffintra=yL2(i)-yh2L2(i)
    zdiffintra=zL2(i)-zh2L2(i)
                                               ! correction for pbc
    if (xdiff > a/2) then
     xdiff = xdiff - a
    else if (xdiff < -a/2) then
     xdiff = xdiff + a
    end if
    if (ydiff > b/2) then
     ydiff = ydiff - b
    else if (ydiff < -b/2) then
     ydiff = ydiff + b
    end if

    if (xdiffintra > a/2) then
     xdiffintra = xdiffintra - a
    else if (xdiffintra < -a/2) then
     xdiffintra = xdiffintra + a
    end if
    if (ydiffintra > b/2) then
     ydiffintra = ydiffintra - b
    else if (ydiffintra < -b/2) then
     ydiffintra = ydiffintra + b
    end if
                                               ! end correction for pbc
    rohinter=SQRT((xdiff)**2+(ydiff)**2+(zdiff)**2)
    rohintra=SQRT((xdiffintra)**2+(ydiffintra)**2+(zdiffintra)**2)
    coseno=(xdiff*xdiffintra+ydiff*ydiffintra+zdiff*zdiffintra)/(rohinter*rohintra)
    if (coseno < -0.76604) then
     zcos=-zcos
     write(19,*) dist, zcos
     write(17,*) dist, zcos
    end if
   end if
  end do


                                                             ! with L1
  do j=1,numL1
   xdiff=xL2(i)-xL1(j)
   ydiff=yL2(i)-yL1(j)
   zdiff=zL2(i)-zL1(j)
                                               ! correction for pbc
   if (xdiff > a/2) then
    xdiff = xdiff - a
   else if (xdiff < -a/2) then
    xdiff = xdiff + a
   end if
   if (ydiff > b/2) then
    ydiff = ydiff - b
   else if (ydiff < -b/2) then
    ydiff = ydiff + b
   end if
                                               ! end correction for pbc
   dist=SQRT((xdiff)**2+(ydiff)**2+(zdiff)**2)
   zcos= -zdiff/dist
   if (dist < cut .and. dist>0.1) then
                                                                              ! O-H1
    xdiff=xL2(i)-xh1L1(j)
    ydiff=yL2(i)-yh1L1(j)
    zdiff=zL2(i)-zh1L1(j)
    xdiffintra=xL1(j)-xh1L1(j)
    ydiffintra=yL1(j)-yh1L1(j)
    zdiffintra=zL1(j)-zh1L1(j)
                                               ! correction for pbc
    if (xdiff > a/2) then
     xdiff = xdiff - a
    else if (xdiff < -a/2) then
     xdiff = xdiff + a
    end if
    if (ydiff > b/2) then
     ydiff = ydiff - b
    else if (ydiff < -b/2) then
     ydiff = ydiff + b
    end if

    if (xdiffintra > a/2) then
     xdiffintra = xdiffintra - a
    else if (xdiffintra < -a/2) then
     xdiffintra = xdiffintra + a
    end if
    if (ydiffintra > b/2) then
     ydiffintra = ydiffintra - b
    else if (ydiffintra < -b/2) then
     ydiffintra = ydiffintra + b
    end if
                                               ! end correction for pbc
    rohinter=SQRT((xdiff)**2+(ydiff)**2+(zdiff)**2)
    rohintra=SQRT((xdiffintra)**2+(ydiffintra)**2+(zdiffintra)**2)
    coseno=(xdiff*xdiffintra+ydiff*ydiffintra+zdiff*zdiffintra)/(rohinter*rohintra)
    if (coseno < -0.76604) then
     write(28,*) dist, zcos
     write(27,*) dist, zcos
    end if
                                                                              ! O-H2
    xdiff=xL2(i)-xh2L1(j)
    ydiff=yL2(i)-yh2L1(j)
    zdiff=zL2(i)-zh2L1(j)
    xdiffintra=xL2(j)-xh2L1(j)
    ydiffintra=yL2(j)-yh2L1(j)
    zdiffintra=zL2(j)-zh2L1(j)
                                               ! correction for pbc
    if (xdiff > a/2) then
     xdiff = xdiff - a
    else if (xdiff < -a/2) then
     xdiff = xdiff + a
    end if
    if (ydiff > b/2) then
     ydiff = ydiff - b
    else if (ydiff < -b/2) then
     ydiff = ydiff + b
    end if

    if (xdiffintra > a/2) then
     xdiffintra = xdiffintra - a
    else if (xdiffintra < -a/2) then
     xdiffintra = xdiffintra + a
    end if
    if (ydiffintra > b/2) then
     ydiffintra = ydiffintra - b
    else if (ydiffintra < -b/2) then
     ydiffintra = ydiffintra + b
    end if
                                               ! end correction for pbc
    rohinter=SQRT((xdiff)**2+(ydiff)**2+(zdiff)**2)
    rohintra=SQRT((xdiffintra)**2+(ydiffintra)**2+(zdiffintra)**2)
    coseno=(xdiff*xdiffintra+ydiff*ydiffintra+zdiff*zdiffintra)/(rohinter*rohintra)
    if (coseno < -0.76604) then
     write(28,*) dist, zcos
     write(27,*) dist, zcos
    end if
                                                                              ! H1-O
    xdiff=xL1(j)-xh1L2(i)
    ydiff=yL1(j)-yh1L2(i)
    zdiff=zL1(j)-zh1L2(i)
    xdiffintra=xL2(i)-xh1L2(i)
    ydiffintra=yL2(i)-yh1L2(i)
    zdiffintra=zL2(i)-zh1L2(i)
                                               ! correction for pbc
    if (xdiff > a/2) then
     xdiff = xdiff - a
    else if (xdiff < -a/2) then
     xdiff = xdiff + a
    end if
    if (ydiff > b/2) then
     ydiff = ydiff - b
    else if (ydiff < -b/2) then
     ydiff = ydiff + b
    end if

    if (xdiffintra > a/2) then
     xdiffintra = xdiffintra - a
    else if (xdiffintra < -a/2) then
     xdiffintra = xdiffintra + a
    end if
    if (ydiffintra > b/2) then
     ydiffintra = ydiffintra - b
    else if (ydiffintra < -b/2) then
     ydiffintra = ydiffintra + b
    end if
                                               ! end correction for pbc
    rohinter=SQRT((xdiff)**2+(ydiff)**2+(zdiff)**2)
    rohintra=SQRT((xdiffintra)**2+(ydiffintra)**2+(zdiffintra)**2)
    coseno=(xdiff*xdiffintra+ydiff*ydiffintra+zdiff*zdiffintra)/(rohinter*rohintra)
    if (coseno < -0.76604) then
     zcos=-zcos
     write(18,*) dist, zcos
     write(17,*) dist, zcos
    end if
                                                                              ! H2-O
    xdiff=xL1(j)-xh2L2(i)
    ydiff=yL1(j)-yh2L2(i)
    zdiff=zL1(j)-zh2L2(i)
    xdiffintra=xL2(i)-xh2L2(i)
    ydiffintra=yL2(i)-yh2L2(i)
    zdiffintra=zL2(i)-zh2L2(i)
                                               ! correction for pbc
    if (xdiff > a/2) then
     xdiff = xdiff - a
    else if (xdiff < -a/2) then
     xdiff = xdiff + a
    end if
    if (ydiff > b/2) then
     ydiff = ydiff - b
    else if (ydiff < -b/2) then
     ydiff = ydiff + b
    end if

    if (xdiffintra > a/2) then
     xdiffintra = xdiffintra - a
    else if (xdiffintra < -a/2) then
     xdiffintra = xdiffintra + a
    end if
    if (ydiffintra > b/2) then
     ydiffintra = ydiffintra - b
    else if (ydiffintra < -b/2) then
     ydiffintra = ydiffintra + b
    end if
                                               ! end correction for pbc
    rohinter=SQRT((xdiff)**2+(ydiff)**2+(zdiff)**2)
    rohintra=SQRT((xdiffintra)**2+(ydiffintra)**2+(zdiffintra)**2)
    coseno=(xdiff*xdiffintra+ydiff*ydiffintra+zdiff*zdiffintra)/(rohinter*rohintra)
    if (coseno < -0.76604) then
     zcos=-zcos
     write(18,*) dist, zcos
     write(17,*) dist, zcos
    end if
   end if
  end do


                                                             ! with L0
  do j=1,numL0
   xdiff=xL2(i)-xL0(j)
   ydiff=yL2(i)-yL0(j)
   zdiff=zL2(i)-zL0(j)
                                               ! correction for pbc
   if (xdiff > a/2) then
    xdiff = xdiff - a
   else if (xdiff < -a/2) then
    xdiff = xdiff + a
   end if
   if (ydiff > b/2) then
    ydiff = ydiff - b
   else if (ydiff < -b/2) then
    ydiff = ydiff + b
   end if
                                               ! end correction for pbc
   dist=SQRT((xdiff)**2+(ydiff)**2+(zdiff)**2)
   zcos= -zdiff/dist
   if (dist < cut .and. dist>0.1) then
                                                                              ! O-H1
    xdiff=xL2(i)-xh1L0(j)
    ydiff=yL2(i)-yh1L0(j)
    zdiff=zL2(i)-zh1L0(j)
    xdiffintra=xL0(j)-xh1L0(j)
    ydiffintra=yL0(j)-yh1L0(j)
    zdiffintra=zL0(j)-zh1L0(j)
                                               ! correction for pbc
    if (xdiff > a/2) then
     xdiff = xdiff - a
    else if (xdiff < -a/2) then
     xdiff = xdiff + a
    end if
    if (ydiff > b/2) then
     ydiff = ydiff - b
    else if (ydiff < -b/2) then
     ydiff = ydiff + b
    end if

    if (xdiffintra > a/2) then
     xdiffintra = xdiffintra - a
    else if (xdiffintra < -a/2) then
     xdiffintra = xdiffintra + a
    end if
    if (ydiffintra > b/2) then
     ydiffintra = ydiffintra - b
    else if (ydiffintra < -b/2) then
     ydiffintra = ydiffintra + b
    end if
                                               ! end correction for pbc
    rohinter=SQRT((xdiff)**2+(ydiff)**2+(zdiff)**2)
    rohintra=SQRT((xdiffintra)**2+(ydiffintra)**2+(zdiffintra)**2)
    coseno=(xdiff*xdiffintra+ydiff*ydiffintra+zdiff*zdiffintra)/(rohinter*rohintra)
    if (coseno < -0.76604) then
     write(38,*) dist, zcos
     write(37,*) dist, zcos
    end if
                                                                              ! O-H2
    xdiff=xL2(i)-xh2L0(j)
    ydiff=yL2(i)-yh2L0(j)
    zdiff=zL2(i)-zh2L0(j)
    xdiffintra=xL0(j)-xh2L0(j)
    ydiffintra=yL0(j)-yh2L0(j)
    zdiffintra=zL0(j)-zh2L0(j)
                                               ! correction for pbc
    if (xdiff > a/2) then
     xdiff = xdiff - a
    else if (xdiff < -a/2) then
     xdiff = xdiff + a
    end if
    if (ydiff > b/2) then
     ydiff = ydiff - b
    else if (ydiff < -b/2) then
     ydiff = ydiff + b
    end if

    if (xdiffintra > a/2) then
     xdiffintra = xdiffintra - a
    else if (xdiffintra < -a/2) then
     xdiffintra = xdiffintra + a
    end if
    if (ydiffintra > b/2) then
     ydiffintra = ydiffintra - b
    else if (ydiffintra < -b/2) then
     ydiffintra = ydiffintra + b
    end if
                                               ! end correction for pbc
    rohinter=SQRT((xdiff)**2+(ydiff)**2+(zdiff)**2)
    rohintra=SQRT((xdiffintra)**2+(ydiffintra)**2+(zdiffintra)**2)
    coseno=(xdiff*xdiffintra+ydiff*ydiffintra+zdiff*zdiffintra)/(rohinter*rohintra)
    if (coseno < -0.76604) then
     write(38,*) dist, zcos
     write(37,*) dist, zcos
    end if
                                                                              ! H1-O                                     
    xdiff=xL0(j)-xh1L2(i)
    ydiff=yL0(j)-yh1L2(i)
    zdiff=zL0(j)-zh1L2(i)
    xdiffintra=xL2(i)-xh1L2(i)
    ydiffintra=yL2(i)-yh1L2(i)
    zdiffintra=zL2(i)-zh1L2(i)
                                               ! correction for pbc
    if (xdiff > a/2) then
     xdiff = xdiff - a
    else if (xdiff < -a/2) then
     xdiff = xdiff + a
    end if
    if (ydiff > b/2) then
     ydiff = ydiff - b
    else if (ydiff < -b/2) then
     ydiff = ydiff + b
    end if

    if (xdiffintra > a/2) then
     xdiffintra = xdiffintra - a
    else if (xdiffintra < -a/2) then
     xdiffintra = xdiffintra + a
    end if
    if (ydiffintra > b/2) then
     ydiffintra = ydiffintra - b
    else if (ydiffintra < -b/2) then
     ydiffintra = ydiffintra + b
    end if
                                               ! end correction for pbc
    rohinter=SQRT((xdiff)**2+(ydiff)**2+(zdiff)**2)
    rohintra=SQRT((xdiffintra)**2+(ydiffintra)**2+(zdiffintra)**2)
    coseno=(xdiff*xdiffintra+ydiff*ydiffintra+zdiff*zdiffintra)/(rohinter*rohintra)
    if (coseno < -0.76604) then
     zcos=-zcos
     write(18,*) dist, zcos
     write(17,*) dist, zcos
    end if
                                                                              ! H2-0
    xdiff=xL0(j)-xh2L2(i)
    ydiff=yL0(j)-yh2L2(i)
    zdiff=zL0(j)-zh2L2(i)
    xdiffintra=xL2(i)-xh2L2(i)
    ydiffintra=yL2(i)-yh2L2(i)
    zdiffintra=zL2(i)-zh2L2(i)
                                               ! correction for pbc
    if (xdiff > a/2) then
     xdiff = xdiff - a
    else if (xdiff < -a/2) then
     xdiff = xdiff + a
    end if
    if (ydiff > b/2) then
     ydiff = ydiff - b
    else if (ydiff < -b/2) then
     ydiff = ydiff + b
    end if

    if (xdiffintra > a/2) then
     xdiffintra = xdiffintra - a
    else if (xdiffintra < -a/2) then
     xdiffintra = xdiffintra + a
    end if
    if (ydiffintra > b/2) then
     ydiffintra = ydiffintra - b
    else if (ydiffintra < -b/2) then
     ydiffintra = ydiffintra + b
    end if
                                               ! end correction for pbc
    rohinter=SQRT((xdiff)**2+(ydiff)**2+(zdiff)**2)
    rohintra=SQRT((xdiffintra)**2+(ydiffintra)**2+(zdiffintra)**2)
    coseno=(xdiff*xdiffintra+ydiff*ydiffintra+zdiff*zdiffintra)/(rohinter*rohintra)
    if (coseno < -0.76604) then
     zcos=-zcos
     write(18,*) dist, zcos
     write(17,*) dist, zcos
    end if
   end if
  end do

 end do



                                                                     ! Layer L1
 do i=1, numL1
                                                                 ! with L1
  do j=1,numL1
   rL1L1=0
   xdiff=xL1(i)-xL1(j)
   ydiff=yL1(i)-yL1(j)
   zdiff=zL1(i)-zL1(j)
                                               ! correction for pbc
   if (xdiff > a/2) then
    xdiff = xdiff - a
   else if (xdiff < -a/2) then
    xdiff = xdiff + a
   end if
   if (ydiff > b/2) then
    ydiff = ydiff - b
   else if (ydiff < -b/2) then
    ydiff = ydiff + b
   end if
                                               ! end correction for pbc
   dist=SQRT((xdiff)**2+(ydiff)**2+(zdiff)**2)
   zcos= -zdiff/dist
   if (dist < cut .and. dist>0.1) then
                                                                              ! H1-O
    xdiff=xL1(j)-xh1L1(i)
    ydiff=yL1(j)-yh1L1(i)
    zdiff=zL1(j)-zh1L1(i)
    xdiffintra=xL1(i)-xh1L1(i)
    ydiffintra=yL1(i)-yh1L1(i)
    zdiffintra=zL1(i)-zh1L1(i)
                                               ! correction for pbc
    if (xdiff > a/2) then
     xdiff = xdiff - a
    else if (xdiff < -a/2) then
     xdiff = xdiff + a
    end if
    if (ydiff > b/2) then
     ydiff = ydiff - b
    else if (ydiff < -b/2) then
     ydiff = ydiff + b
    end if

    if (xdiffintra > a/2) then
     xdiffintra = xdiffintra - a
    else if (xdiffintra < -a/2) then
     xdiffintra = xdiffintra + a
    end if
    if (ydiffintra > b/2) then
     ydiffintra = ydiffintra - b
    else if (ydiffintra < -b/2) then
     ydiffintra = ydiffintra + b
    end if
                                               ! end correction for pbc
    rohinter=SQRT((xdiff)**2+(ydiff)**2+(zdiff)**2)
    rohintra=SQRT((xdiffintra)**2+(ydiffintra)**2+(zdiffintra)**2)
    coseno=(xdiff*xdiffintra+ydiff*ydiffintra+zdiff*zdiffintra)/(rohinter*rohintra)
    if (coseno < -0.76604) then
     zcos=-zcos
     write(29,*) dist, zcos
     write(27,*) dist, zcos
    end if
                                                                              ! H2-O
    xdiff=xL1(j)-xh2L1(i)
    ydiff=yL1(j)-yh2L1(i)
    zdiff=zL1(j)-zh2L1(i)
    xdiffintra=xL1(i)-xh2L1(i)
    ydiffintra=yL1(i)-yh2L1(i)
    zdiffintra=zL1(i)-zh2L1(i)
                                               ! correction for pbc
    if (xdiff > a/2) then
     xdiff = xdiff - a
    else if (xdiff < -a/2) then
     xdiff = xdiff + a
    end if
    if (ydiff > b/2) then
     ydiff = ydiff - b
    else if (ydiff < -b/2) then
     ydiff = ydiff + b
    end if

    if (xdiffintra > a/2) then
     xdiffintra = xdiffintra - a
    else if (xdiffintra < -a/2) then
     xdiffintra = xdiffintra + a
    end if
    if (ydiffintra > b/2) then
     ydiffintra = ydiffintra - b
    else if (ydiffintra < -b/2) then
     ydiffintra = ydiffintra + b
    end if
                                               ! end correction for pbc
    rohinter=SQRT((xdiff)**2+(ydiff)**2+(zdiff)**2)
    rohintra=SQRT((xdiffintra)**2+(ydiffintra)**2+(zdiffintra)**2)
    coseno=(xdiff*xdiffintra+ydiff*ydiffintra+zdiff*zdiffintra)/(rohinter*rohintra)
    if (coseno < -0.76604) then
     zcos=-zcos
     write(29,*) dist, zcos
     write(27,*) dist, zcos
    end if
   
   end if
  end do

                                                                 ! with L0
  do j=1,numL0
   xdiff=xL1(i)-xL0(j)
   ydiff=yL1(i)-yL0(j)
   zdiff=zL1(i)-zL0(j)
                                               ! correction for pbc
   if (xdiff > a/2) then
    xdiff = xdiff - a
   else if (xdiff < -a/2) then
    xdiff = xdiff + a
   end if
   if (ydiff > b/2) then
    ydiff = ydiff - b
   else if (ydiff < -b/2) then
    ydiff = ydiff + b
   end if
                                               ! end correction for pbc
   dist=SQRT((xdiff)**2+(ydiff)**2+(zdiff)**2)
   zcos= -zdiff/dist
   if (dist < cut .and. dist>0.1) then
                                                                              ! O-H1
    xdiff=xL1(i)-xh1L0(j)
    ydiff=yL1(i)-yh1L0(j)
    zdiff=zL1(i)-zh1L0(j)
    xdiffintra=xL0(j)-xh1L0(j)
    ydiffintra=yL0(j)-yh1L0(j)
    zdiffintra=zL0(j)-zh1L0(j)
                                               ! correction for pbc
    if (xdiff > a/2) then
     xdiff = xdiff - a
    else if (xdiff < -a/2) then
     xdiff = xdiff + a
    end if
    if (ydiff > b/2) then
     ydiff = ydiff - b
    else if (ydiff < -b/2) then
     ydiff = ydiff + b
    end if

    if (xdiffintra > a/2) then
     xdiffintra = xdiffintra - a
    else if (xdiffintra < -a/2) then
     xdiffintra = xdiffintra + a
    end if
    if (ydiffintra > b/2) then
     ydiffintra = ydiffintra - b
    else if (ydiffintra < -b/2) then
     ydiffintra = ydiffintra + b
    end if
                                               ! end correction for pbc
    rohinter=SQRT((xdiff)**2+(ydiff)**2+(zdiff)**2)
    rohintra=SQRT((xdiffintra)**2+(ydiffintra)**2+(zdiffintra)**2)
    coseno=(xdiff*xdiffintra+ydiff*ydiffintra+zdiff*zdiffintra)/(rohinter*rohintra)
    if (coseno < -0.76604) then
     write(38,*) dist, zcos
     write(37,*) dist, zcos
    end if
                                                                              ! O-H2
    xdiff=xL1(i)-xh2L0(j)
    ydiff=yL1(i)-yh2L0(j)
    zdiff=zL1(i)-zh2L0(j)
    xdiffintra=xL0(j)-xh2L0(j)
    ydiffintra=yL0(j)-yh2L0(j)
    zdiffintra=zL0(j)-zh2L0(j)
                                               ! correction for pbc
    if (xdiff > a/2) then
     xdiff = xdiff - a
    else if (xdiff < -a/2) then
     xdiff = xdiff + a
    end if
    if (ydiff > b/2) then
     ydiff = ydiff - b
    else if (ydiff < -b/2) then
     ydiff = ydiff + b
    end if

    if (xdiffintra > a/2) then
     xdiffintra = xdiffintra - a
    else if (xdiffintra < -a/2) then
     xdiffintra = xdiffintra + a
    end if
    if (ydiffintra > b/2) then
     ydiffintra = ydiffintra - b
    else if (ydiffintra < -b/2) then
     ydiffintra = ydiffintra + b
    end if
                                               ! end correction for pbc
    rohinter=SQRT((xdiff)**2+(ydiff)**2+(zdiff)**2)
    rohintra=SQRT((xdiffintra)**2+(ydiffintra)**2+(zdiffintra)**2)
    coseno=(xdiff*xdiffintra+ydiff*ydiffintra+zdiff*zdiffintra)/(rohinter*rohintra)
    if (coseno < -0.76604) then
     write(38,*) dist, zcos
     write(37,*) dist, zcos
    end if
                                                                              ! H1-O
    xdiff=xL0(j)-xh1L1(i)
    ydiff=yL0(j)-yh1L1(i)
    zdiff=zL0(j)-zh1L1(i)
    xdiffintra=xL1(i)-xh1L1(i)
    ydiffintra=yL1(i)-yh1L1(i)
    zdiffintra=zL1(i)-zh1L1(i)
                                               ! correction for pbc
    if (xdiff > a/2) then
     xdiff = xdiff - a
    else if (xdiff < -a/2) then
     xdiff = xdiff + a
    end if
    if (ydiff > b/2) then
     ydiff = ydiff - b
    else if (ydiff < -b/2) then
     ydiff = ydiff + b
    end if

    if (xdiffintra > a/2) then
     xdiffintra = xdiffintra - a
    else if (xdiffintra < -a/2) then
     xdiffintra = xdiffintra + a
    end if
    if (ydiffintra > b/2) then
     ydiffintra = ydiffintra - b
    else if (ydiffintra < -b/2) then
     ydiffintra = ydiffintra + b
    end if
                                               ! end correction for pbc
    rohinter=SQRT((xdiff)**2+(ydiff)**2+(zdiff)**2)
    rohintra=SQRT((xdiffintra)**2+(ydiffintra)**2+(zdiffintra)**2)
    coseno=(xdiff*xdiffintra+ydiff*ydiffintra+zdiff*zdiffintra)/(rohinter*rohintra)
    if (coseno < -0.76604) then
     zcos=-zcos
     write(28,*) dist, zcos
     write(27,*) dist, zcos
    end if
                                                                              ! H2-O
    xdiff=xL0(j)-xh2L1(i)
    ydiff=yL0(j)-yh2L1(i)
    zdiff=zL0(j)-zh2L1(i)
    xdiffintra=xL1(i)-xh2L1(i)
    ydiffintra=yL1(i)-yh2L1(i)
    zdiffintra=zL1(i)-zh2L1(i)
                                               ! correction for pbc
    if (xdiff > a/2) then
     xdiff = xdiff - a
    else if (xdiff < -a/2) then
     xdiff = xdiff + a
    end if
    if (ydiff > b/2) then
     ydiff = ydiff - b
    else if (ydiff < -b/2) then
     ydiff = ydiff + b
    end if

    if (xdiffintra > a/2) then
     xdiffintra = xdiffintra - a
    else if (xdiffintra < -a/2) then
     xdiffintra = xdiffintra + a
    end if
    if (ydiffintra > b/2) then
     ydiffintra = ydiffintra - b
    else if (ydiffintra < -b/2) then
     ydiffintra = ydiffintra + b
    end if
                                               ! end correction for pbc
    rohinter=SQRT((xdiff)**2+(ydiff)**2+(zdiff)**2)
    rohintra=SQRT((xdiffintra)**2+(ydiffintra)**2+(zdiffintra)**2)
    coseno=(xdiff*xdiffintra+ydiff*ydiffintra+zdiff*zdiffintra)/(rohinter*rohintra)
    if (coseno < -0.76604) then
     zcos=-zcos
     write(28,*) dist, zcos
     write(27,*) dist, zcos
    end if

   end if
  end do

 end do


                                                                     ! Layer L0
 if (numL0>0) then
 do i=1, numL0
                                                                 ! with L0
  do j=1,numL0
   xdiff=xL0(i)-xL0(j)
   ydiff=yL0(i)-yL0(j)
   zdiff=zL0(i)-zL0(j)
                                               ! correction for pbc
   if (xdiff > a/2) then
    xdiff = xdiff - a
   else if (xdiff < -a/2) then
    xdiff = xdiff + a
   end if
   if (ydiff > b/2) then
    ydiff = ydiff - b
   else if (ydiff < -b/2) then
    ydiff = ydiff + b
   end if
                                               ! end correction for pbc
   dist=SQRT((xdiff)**2+(ydiff)**2+(zdiff)**2)
   zcos= -zdiff/dist
   if (dist < cut .and. dist>0.1) then
                                                                              ! H1-O
    xdiff=xL0(j)-xh1L0(i)
    ydiff=yL0(j)-yh1L0(i)
    zdiff=zL0(j)-zh1L0(i)
    xdiffintra=xL0(i)-xh1L0(i)
    ydiffintra=yL0(i)-yh1L0(i)
    zdiffintra=zL0(i)-zh1L0(i)
                                              ! correction for pbc
    if (xdiff > a/2) then
     xdiff = xdiff - a
    else if (xdiff < -a/2) then
     xdiff = xdiff + a
    end if
    if (ydiff > b/2) then
     ydiff = ydiff - b
    else if (ydiff < -b/2) then
     ydiff = ydiff + b
    end if

    if (xdiffintra > a/2) then
     xdiffintra = xdiffintra - a
    else if (xdiffintra < -a/2) then
     xdiffintra = xdiffintra + a
    end if
    if (ydiffintra > b/2) then
     ydiffintra = ydiffintra - b
    else if (ydiffintra < -b/2) then
     ydiffintra = ydiffintra + b
    end if
                                               ! end correction for pbc
    rohinter=SQRT((xdiff)**2+(ydiff)**2+(zdiff)**2)
    rohintra=SQRT((xdiffintra)**2+(ydiffintra)**2+(zdiffintra)**2)
    coseno=(xdiff*xdiffintra+ydiff*ydiffintra+zdiff*zdiffintra)/(rohinter*rohintra)
    if (coseno < -0.76604) then
     zcos=-zcos
     write(39,*) dist, zcos
     write(37,*) dist, zcos
    end if
                                                                              ! H2-O
    xdiff=xL0(j)-xh2L0(i)
    ydiff=yL0(j)-yh2L0(i)
    zdiff=zL0(j)-zh2L0(i)
    xdiffintra=xL0(i)-xh2L0(i)
    ydiffintra=yL0(i)-yh2L0(i)
    zdiffintra=zL0(i)-zh2L0(i)
                                              ! correction for pbc
    if (xdiff > a/2) then
     xdiff = xdiff - a
    else if (xdiff < -a/2) then
     xdiff = xdiff + a
    end if
    if (ydiff > b/2) then
     ydiff = ydiff - b
    else if (ydiff < -b/2) then
     ydiff = ydiff + b
    end if

    if (xdiffintra > a/2) then
     xdiffintra = xdiffintra - a
    else if (xdiffintra < -a/2) then
     xdiffintra = xdiffintra + a
    end if
    if (ydiffintra > b/2) then
     ydiffintra = ydiffintra - b
    else if (ydiffintra < -b/2) then
     ydiffintra = ydiffintra + b
    end if
                                               ! end correction for pbc
    rohinter=SQRT((xdiff)**2+(ydiff)**2+(zdiff)**2)
    rohintra=SQRT((xdiffintra)**2+(ydiffintra)**2+(zdiffintra)**2)
    coseno=(xdiff*xdiffintra+ydiff*ydiffintra+zdiff*zdiffintra)/(rohinter*rohintra)
    if (coseno < -0.76604) then
     zcos=-zcos
     write(39,*) dist, zcos
     write(37,*) dist, zcos
    end if

   end if
  end do

 end do
 end if
end do

deallocate (xL0)
deallocate (yL0)
deallocate (zL0)
deallocate (xL1)
deallocate (yL1)
deallocate (zL1)
deallocate (xL2)
deallocate (yL2)
deallocate (zL2)
deallocate (xL3)
deallocate (yL3)
deallocate (zL3)
deallocate (xh1L0)
deallocate (yh1L0)
deallocate (zh1L0)
deallocate (xh2L0)
deallocate (yh2L0)
deallocate (zh2L0)
deallocate (xh1L1)
deallocate (yh1L1)
deallocate (zh1L1)
deallocate (xh2L1)
deallocate (yh2L1)
deallocate (zh2L1)
deallocate (xh1L2)
deallocate (yh1L2)
deallocate (zh1L2)
deallocate (xh2L2)
deallocate (yh2L2)
deallocate (zh2L2)
deallocate (xh1L3)
deallocate (yh1L3)
deallocate (zh1L3)
deallocate (xh2L3)
deallocate (yh2L3)
deallocate (zh2L3)
deallocate (oxL2)
deallocate (oyL2)
deallocate (ozL2)

stop
end program
