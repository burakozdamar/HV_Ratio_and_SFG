program down_hbond

! Evry 28/02/2017 

implicit none
integer*8 :: i, j, step, Natom, numL3, numL2, numL1, numL0, Nstep, io1, io2
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
real*8 :: totL3=0, totL2=0, totL1=0, totL0=0, a, b, c
real*8, parameter :: cut=3.2
real*8, dimension(:) , allocatable :: xL3, yL3, zL3, xL2, yL2, zL2, xL1, yL1, zL1, xL0, yL0, zL0 
real*8, dimension(:) , allocatable :: xh1L3, yh1L3, zh1L3, xh1L2, yh1L2, zh1L2, xh1L1, yh1L1, zh1L1, xh1L0, yh1L0, zh1L0
real*8, dimension(:) , allocatable :: xh2L3, yh2L3, zh2L3, xh2L2, yh2L2, zh2L2, xh2L1, yh2L1, zh2L1, xh2L0, yh2L0, zh2L0
character(5) :: atom, hyd
integer*8 :: NO, NstepTOT, istart, NstepReal

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

open (1, FILE="up_down_water_layers.xyz")
open (2, FILE="down_hbond")
open (3, FILE="down_distance_hbond")


do step = 1, Nstep


                                    ! reading atomic coordinate
L0u=0
L1u=0
L2u=0
L3u=0
aL0u=0
aL1u=0
aL2u=0
aL3u=0
dL0u=0
dL1u=0
dL2u=0
dL3u=0
numL3=0
numL2=0
numL1=0
numL0=0
 read(1,*,iostat=io1) Natom
    if (io1 /= 0) then
        exit
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
  
 totL3=totL3+numL3
 totL2=totL2+numL2
 totL1=totL1+numL1
 totL0=totL0+numL0

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
   if (dist < cut .and. dist>0.1) then
                                                                              ! O-H1
    xdiff=xL3(i)-xh1L3(j)
    ydiff=yL3(i)-yh1L3(j)
    zdiff=zL3(i)-zh1L3(j)
    xdiffintra=xL3(j)-xh1L3(j)
    ydiffintra=yL3(j)-yh1L3(j)
    zdiffintra=zL3(j)-zh1L3(j)
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
    L3u=L3u+1
    aL3u=aL3u+1
    nL3L3=nL3L3+1
    rL3L3=dist
    end if
                                                                              ! O-H2
    xdiff=xL3(i)-xh2L3(j)
    ydiff=yL3(i)-yh2L3(j)
    zdiff=zL3(i)-zh2L3(j)
    xdiffintra=xL3(j)-xh2L3(j)
    ydiffintra=yL3(j)-yh2L3(j)
    zdiffintra=zL3(j)-zh2L3(j)
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
     L3u=L3u+1
     aL3u=aL3u+1
     nL3L3=nL3L3+1
     rL3L3=dist
    end if
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
     L3u=L3u+1
     dL3u=dL3u+1
     nL3L3=nL3L3+1
     rL3L3=dist
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
     L3u=L3u+1
     dL3u=dL3u+1
     nL3L3=nL3L3+1
     rL3L3=dist
    end if

   end if
   avrL3L3=avrL3L3+rL3L3
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
     L2u=L2u+1
     aL2u=aL2u+1   
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
     L2u=L2u+1
     aL2u=aL2u+1
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
     L2u=L2u+1
     dL2u=dL2u+1
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
     L2u=L2u+1
     dL2u=dL2u+1
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
     L1u=L1u+1
     aL1u=aL1u+1
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
     L1u=L1u+1
     aL1u=aL1u+1
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
     L1u=L1u+1
     dL1u=dL1u+1
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
     L1u=L1u+1
     dL1u=dL1u+1
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
     L0u=L0u+1
     aL0u=aL0u+1
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
     L0u=L0u+1
     aL0u=aL0u+1
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
     L0u=L0u+1
     dL0u=dL0u+1
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
     L0u=L0u+1
     dL0u=dL0u+1
    end if
   end if
  end do
 
 end do

 L3L0coord=L3L0coord+L0u/numL3
 L3L1coord=L3L1coord+L1u/numL3
 Lcutcoord=Lcutcoord+L2u/numL3
 L3L3coord=L3L3coord+L3u/numL3
 L3intercoord=L3L1coord+Lcutcoord+L3L0coord

 if (numL0>0) then
 L0L3coord=L0L3coord+L0u/numL0
 end if
 L1L3coord=L1L3coord+L1u/numL1
 L2L3coord=L2L3coord+L2u/numL2


 aL3L0coord=aL3L0coord+aL0u/numL3
 aL3L1coord=aL3L1coord+aL1u/numL3
 aLcutcoord=aLcutcoord+aL2u/numL3
 aL3L3coord=aL3L3coord+aL3u/numL3
 aL3intercoord=aL3L1coord+aLcutcoord+aL3L0coord

 if (numL0>0) then
 aL0L3coord=aL0L3coord+dL0u/numL0
 end if
 aL1L3coord=aL1L3coord+dL1u/numL1
 aL2L3coord=aL2L3coord+dL2u/numL2


 dL3L0coord=dL3L0coord+dL0u/numL3
 dL3L1coord=dL3L1coord+dL1u/numL3
 dLcutcoord=dLcutcoord+dL2u/numL3
 dL3L3coord=dL3L3coord+dL3u/numL3
 dL3intercoord=dL3L1coord+dLcutcoord+dL3L0coord

 if (numL0>0) then
 dL0L3coord=dL0L3coord+aL0u/numL0
 end if
 dL1L3coord=dL1L3coord+aL1u/numL1
 dL2L3coord=dL2L3coord+aL2u/numL2


                                                                     ! Layer L2
 L0u=0
 L1u=0
 L2u=0
 L3u=0
 aL0u=0
 aL1u=0
 aL2u=0
 aL3u=0
 dL0u=0
 dL1u=0
 dL2u=0
 dL3u=0
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
   if (dist < cut .and. dist>0.1) then
                                                                              ! O-H1
    xdiff=xL2(i)-xh1L2(j)
    ydiff=yL2(i)-yh1L2(j)
    zdiff=zL2(i)-zh1L2(j)
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
     L2u=L2u+1
     aL2u=aL2u+1
    end if
                                                                              ! O-H2
    xdiff=xL2(i)-xh2L2(j)
    ydiff=yL2(i)-yh2L2(j)
    zdiff=zL2(i)-zh2L2(j)
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
     L2u=L2u+1
     aL2u=aL2u+1
    end if
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
     L2u=L2u+1
     dL2u=dL2u+1
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
     L2u=L2u+1
     dL2u=dL2u+1
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
     L1u=L1u+1
     aL1u=aL1u+1
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
     L1u=L1u+1
     aL1u=aL1u+1
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
     L1u=L1u+1
     dL1u=dL1u+1
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
     L1u=L1u+1
     dL1u=dL1u+1
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
     L0u=L0u+1
     aL0u=aL0u+1
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
     L0u=L0u+1
     aL0u=aL0u+1
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
     L0u=L0u+1
     dL0u=dL0u+1
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
     L0u=L0u+1
     dL0u=dL0u+1
    end if
   end if
  end do

 end do

 L2L0coord=L2L0coord+L0u/numL2
 L2L1coord=L2L1coord+L1u/numL2
 L2L2coord=L2L2coord+L2u/numL2
 L2intercoord=L2L1coord+L2L3coord+L2L0coord

 if (numL0>0) then
 L0L2coord=L0L2coord+L0u/numL0
 end if
 L1L2coord=L1L2coord+L1u/numL1


 aL2L0coord=aL2L0coord+aL0u/numL2
 aL2L1coord=aL2L1coord+aL1u/numL2
 aL2L2coord=aL2L2coord+aL2u/numL2
 aL2intercoord=aL2L1coord+aL2L3coord+aL2L0coord

 if (numL0>0) then
 aL0L2coord=aL0L2coord+dL0u/numL0
 end if
 aL1L2coord=aL1L2coord+dL1u/numL1


 dL2L0coord=dL2L0coord+dL0u/numL2
 dL2L1coord=dL2L1coord+dL1u/numL2
 dL2L2coord=dL2L2coord+dL2u/numL2
 dL2intercoord=dL2L1coord+dL2L3coord+dL2L0coord

 if (numL0>0) then
 dL0L2coord=dL0L2coord+aL0u/numL0
 end if
 dL1L2coord=dL1L2coord+aL1u/numL1


                                                                     ! Layer L1
 L0u=0
 L1u=0
 L2u=0
 L3u=0
 aL0u=0
 aL1u=0
 aL2u=0
 aL3u=0
 dL0u=0
 dL1u=0
 dL2u=0
 dL3u=0
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
   if (dist < cut .and. dist>0.1) then
                                                                              ! O-H1
    xdiff=xL1(i)-xh1L1(j)
    ydiff=yL1(i)-yh1L1(j)
    zdiff=zL1(i)-zh1L1(j)
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
     L1u=L1u+1
     aL1u=aL1u+1
     nL1L1=nL1L1+1
     rL1L1=dist
    end if
                                                                              ! O-H2
    xdiff=xL1(i)-xh2L1(j)
    ydiff=yL1(i)-yh2L1(j)
    zdiff=zL1(i)-zh2L1(j)
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
     L1u=L1u+1
     aL1u=aL1u+1
     nL1L1=nL1L1+1
     rL1L1=dist
    end if
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
     L1u=L1u+1
     dL1u=dL1u+1
     nL1L1=nL1L1+1
     rL1L1=dist
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
     L1u=L1u+1
     dL1u=dL1u+1
     nL1L1=nL1L1+1
     rL1L1=dist
    end if
   
   end if
   avrL1L1=avrL1L1+rL1L1
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
     L0u=L0u+1
     aL0u=aL0u+1
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
     L0u=L0u+1
     aL0u=aL0u+1
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
     L0u=L0u+1
     dL0u=dL0u+1
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
     L0u=L0u+1
     dL0u=dL0u+1
    end if

   end if
  end do

 end do

 L1L0coord=L1L0coord+L0u/numL1
 L1L1coord=L1L1coord+L1u/numL1
 L1intercoord=L1L2coord+L1L3coord+L1L0coord

 if (numL0>0) then
 L0L1coord=L0L1coord+L0u/numL0
 end if

 aL1L0coord=aL1L0coord+aL0u/numL1
 aL1L1coord=aL1L1coord+aL1u/numL1
 aL1intercoord=aL1L2coord+aL1L3coord+aL1L0coord

 if (numL0>0) then
 aL0L1coord=aL0L1coord+dL0u/numL0
 end if

 dL1L0coord=dL1L0coord+dL0u/numL1
 dL1L1coord=dL1L1coord+dL1u/numL1
 dL1intercoord=dL1L2coord+dL1L3coord+dL1L0coord

 if (numL0>0) then
 dL0L1coord=dL0L1coord+aL0u/numL0
 end if

                                                                     ! Layer L0
 L0u=0
 L1u=0
 L2u=0
 L3u=0
 aL0u=0
 aL1u=0
 aL2u=0
 aL3u=0
 dL0u=0
 dL1u=0
 dL2u=0
 dL3u=0

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
   if (dist < cut .and. dist>0.1) then
                                                                              ! O-H1
    xdiff=xL0(i)-xh1L0(j)
    ydiff=yL0(i)-yh1L0(j)
    zdiff=zL0(i)-zh1L0(j)
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
     L0u=L0u+1
     aL0u=aL0u+1

    end if
                                                                              ! O-H2
    xdiff=xL0(i)-xh2L0(j)
    ydiff=yL0(i)-yh2L0(j)
    zdiff=zL0(i)-zh2L0(j)
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
     L0u=L0u+1
     aL0u=aL0u+1
    end if
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
     L0u=L0u+1
     dL0u=dL0u+1
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
     L0u=L0u+1
     dL0u=dL0u+1
    end if

   end if
  end do

 end do

 L0L0coord=L0L0coord+L0u/numL0
 L0intercoord=L0L1coord+L0L3coord+L0L2coord

 aL0L0coord=aL0L0coord+aL0u/numL0
 aL0intercoord=aL0L1coord+aL0L3coord+aL0L2coord

 dL0L0coord=dL0L0coord+dL0u/numL0
 dL0intercoord=dL0L1coord+dL0L3coord+dL0L2coord
 end if
end do
NstepReal = step - 1

L3L0coord=L3L0coord/NstepReal
L3L1coord=L3L1coord/NstepReal
Lcutcoord=Lcutcoord/NstepReal
L3L3coord=L3L3coord/NstepReal
L3intercoord=L3intercoord/NstepReal

L2L0coord=L2L0coord/NstepReal
L2L1coord=L2L1coord/NstepReal
L2L2coord=L2L2coord/NstepReal
L2L3coord=L2L3coord/NstepReal
L2intercoord=L2intercoord/NstepReal

L1L0coord=L1L0coord/NstepReal
L1L1coord=L1L1coord/NstepReal
L1L2coord=L1L2coord/NstepReal
L1L3coord=L1L3coord/NstepReal
L1intercoord=L1intercoord/NstepReal

L0L0coord=L0L0coord/NstepReal
L0L1coord=L0L1coord/NstepReal
L0L2coord=L0L2coord/NstepReal
L0L3coord=L0L3coord/NstepReal
L0intercoord=L0intercoord/NstepReal



aL3L0coord=aL3L0coord/NstepReal
aL3L1coord=aL3L1coord/NstepReal
aLcutcoord=aLcutcoord/NstepReal
aL3L3coord=aL3L3coord/NstepReal
aL3intercoord=aL3intercoord/NstepReal

aL2L0coord=aL2L0coord/NstepReal
aL2L1coord=aL2L1coord/NstepReal
aL2L2coord=aL2L2coord/NstepReal
aL2L3coord=aL2L3coord/NstepReal
aL2intercoord=aL2intercoord/NstepReal

aL1L0coord=aL1L0coord/NstepReal
aL1L1coord=aL1L1coord/NstepReal
aL1L2coord=aL1L2coord/NstepReal
aL1L3coord=aL1L3coord/NstepReal
aL1intercoord=aL1intercoord/NstepReal

aL0L0coord=aL0L0coord/NstepReal
aL0L1coord=aL0L1coord/NstepReal
aL0L2coord=aL0L2coord/NstepReal
aL0L3coord=aL0L3coord/NstepReal
aL0intercoord=aL0intercoord/NstepReal



dL3L0coord=dL3L0coord/NstepReal
dL3L1coord=dL3L1coord/NstepReal
dLcutcoord=dLcutcoord/NstepReal
dL3L3coord=dL3L3coord/NstepReal
dL3intercoord=dL3intercoord/NstepReal

dL2L0coord=dL2L0coord/NstepReal
dL2L1coord=dL2L1coord/NstepReal
dL2L2coord=dL2L2coord/NstepReal
dL2L3coord=dL2L3coord/NstepReal
dL2intercoord=dL2intercoord/NstepReal

dL1L0coord=dL1L0coord/NstepReal
dL1L1coord=dL1L1coord/NstepReal
dL1L2coord=dL1L2coord/NstepReal
dL1L3coord=dL1L3coord/NstepReal
dL1intercoord=dL1intercoord/NstepReal

dL0L0coord=dL0L0coord/NstepReal
dL0L1coord=dL0L1coord/NstepReal
dL0L2coord=dL0L2coord/NstepReal
dL0L3coord=dL0L3coord/NstepReal
dL0intercoord=dL0intercoord/NstepReal

avrL1L1 = avrL1L1/nL1L1
avrL3L3 = avrL3L3/nL3L3

totL3=totL3/NstepReal
totL2=totL2/NstepReal
totL1=totL1/NstepReal
totL0=totL0/NstepReal


write(2,*) "*************************************************"
write(2,*) "Layer L3"  ,  totL3
write(2,*) "*************************************************"
write(2,*) "H-Bond x molecule =",L3L3coord+L3intercoord,"as donor=",dL3L3coord+dL3intercoord,"as acceptor=",aL3L3coord+aL3intercoord
write(2,*) "intralayer H-Bond x molecule =", L3L3coord, "as donor=", dL3L3coord, "as acceptor=", aL3L3coord 
write(2,*) "interlayer H-Bond x molecule =", L3intercoord, "as donor=", dL3intercoord, "as acceptor=", aL3intercoord
write(2,*) "with L0 =", L3L0coord, "as donor=", dL3L0coord, "as acceptor=", aL3L0coord
write(2,*) "with L1 =", L3L1coord, "as donor=", dL3L1coord, "as acceptor=", aL3L1coord
write(2,*) "with L2 =", Lcutcoord, "as donor=", dLcutcoord, "as acceptor=", aLcutcoord
write(2,*) "*************************************************"
write(2,*) "Layer L2"  ,  totL2
write(2,*) "*************************************************"
write(2,*) "H-Bond x molecule =",L2L2coord+L2intercoord,"as donor=",dL2L2coord+dL2intercoord,"as acceptor=",aL2L2coord+aL2intercoord
write(2,*) "intralayer H-Bond x molecule =", L2L2coord, "as donor=",dL2L2coord, "as acceptor=",aL2L2coord
write(2,*) "interlayer H-Bond x molecule =", L2intercoord, "as donor=",dL2intercoord, "as acceptor=",aL2intercoord
write(2,*) "with L0 =", L2L0coord, "as donor=", dL2L0coord, "as acceptor=", aL2L0coord
write(2,*) "with L1 =", L2L1coord, "as donor=", dL2L1coord, "as acceptor=", aL2L1coord
write(2,*) "with L3 =", L2L3coord, "as donor=", dL2L3coord, "as acceptor=", aL2L3coord
write(2,*) "*************************************************"
write(2,*) "Layer L1"  ,  totL1
write(2,*) "*************************************************"
write(2,*) "H-Bond x molecule =",L1L1coord+L1intercoord,"as donor=",dL1L1coord+dL1intercoord,"as acceptor=",aL1L1coord+aL1intercoord
write(2,*) "intralayer H-Bond x molecule =", L1L1coord, "as donor=", dL1L1coord, "as acceptor=", aL1L1coord
write(2,*) "interlayer H-Bond x molecule =", L1intercoord, "as donor=", dL1intercoord, "as acceptor=", aL1intercoord
write(2,*) "with L0 =", L1L0coord, "as donor=", dL1L0coord, "as acceptor=", aL1L0coord
write(2,*) "with L2 =", L1L2coord, "as donor=", dL1L2coord, "as acceptor=", aL1L2coord
write(2,*) "with L3 =", L1L3coord, "as donor=", dL1L3coord, "as acceptor=", aL1L3coord
write(2,*) "*************************************************"
write(2,*) "Layer L0"  ,  totL0
write(2,*) "*************************************************"
write(2,*) "H-Bond x molecule =",L0L0coord+L0intercoord,"as donor=",dL0L0coord+dL0intercoord,"as acceptor=",aL0L0coord+aL0intercoord
write(2,*) "intralayer H-Bond x molecule =", L0L0coord, "as donor=", dL0L0coord, "as acceptor=", aL0L0coord
write(2,*) "interlayer H-Bond x molecule =", L0intercoord, "as donor=", dL0intercoord, "as acceptor=", aL0intercoord
write(2,*) "with L1 =", L0L1coord, "as donor=", dL0L1coord, "as acceptor=", aL0L1coord
write(2,*) "with L2 =", L0L2coord, "as donor=", dL0L2coord, "as acceptor=", aL0L2coord
write(2,*) "with L3 =", L0L3coord, "as donor=", dL0L3coord, "as acceptor=", aL0L3coord
write(2,*) "*************************************************"




write(3,*) "average hbond dist L1-L1 =", avrL1L1
write(3,*) "average hbond dist L3-L3 =", avrL3L3

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
stop
 end program
