program SurfaceHB2L1res

! Evry 28/02/2017 
! Modified 06/09/2020 Massy-Palaiseau ECLA
! Wanlin
implicit none
integer*8 :: i, j, step, Natom, Nsolid, NSO, NWO, io1, io2, numSO
real*8 :: xo, yo, zo, dist, xdiff, ydiff, zdiff, xdiffintra, ydiffintra, zdiffintra
real*8 :: xs, ys, zs
real*8 :: x1, y1, z1, x2, y2, z2, rohinter, rohintra, coseno, Lcutcoord=0
real*8 :: WS, WStot, densityHsurface, surface
real*8 :: totL1=0, totL0=0, totSO=0, a=13.386d0, b=13.286d0, c=85.d0
real*8, parameter :: cut=3.3
real*8, dimension(:) , allocatable :: xWO, yWO, zWO, xSO, ySO, zSO
real*8, dimension(:) , allocatable :: xh1, yh1, zh1 
real*8, dimension(:) , allocatable :: xh2, yh2, zh2
character(5) :: atom, hyd
integer*8 :: numWO = 120
integer*8 :: NstepTOT, Nstep, Nskip, Istart, NstepReal 
!working flag for getting info
character(256)  txt

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
!###########################################
! STARTS OF THE PROGRAM
!##########################################


!----------------------------------------
! GET INFO FROM BOXDATA
!-----------------------------------------
numSO=0
open (12, FILE="pos_rebuilt_solid.xyz")
 read(12,*) Natom
 read(12,*)
 do i=1, numWO*3 
  read(12,*)
 enddo
 Nsolid=Natom-numWO*3
 do i=1, Nsolid
  read(12,*) atom, xs, ys, zs 
   if (atom == "O") then
    numSO=numSO+1 
   endif 
 enddo
close(12)

allocate (xSO(numSO))
allocate (ySO(numSO))
allocate (zSO(numSO))
allocate (xWO(numWO))
allocate (yWO(numWO))
allocate (zWO(numWO))
allocate (xh1(numWO))
allocate (yh1(numWO))
allocate (zh1(numWO))
allocate (xh2(numWO))
allocate (yh2(numWO))
allocate (zh2(numWO))


open (1, FILE="pos_rebuilt_solid.xyz")
open (2, FILE="PEG2Water.dat")

! skip steps

Nskip = (Istart-1) * (Natom+2)

do step = 1, Nskip
    read (1,*)
enddo
WStot = 0
do step = 1, Nstep

 
 NWO = 0
 NSO = 0
 ! reading atomic coordinate
 read(1,*,iostat=io1)
    if (io1 /= 0) then
        exit
    endif
 read(1,*)
 do i=1, numWO 
  read(1,*) atom, xo, yo, zo
  read(1,*) hyd, x1, y1, z1
  read(1,*) hyd, x2, y2, z2
  if (atom == "O"  ) then
   !NWO is number of O of water molecules
   NWO=NWO+1
   xWO(NWO)=xo
   yWO(NWO)=yo
   zWO(NWO)=zo
   xh1(NWO)=x1
   yh1(NWO)=y1
   zh1(NWO)=z1
   xh2(NWO)=x2
   yh2(NWO)=y2
   zh2(NWO)=z2
  end if
 end do
 do i = 1, Nsolid
  read(1,*) atom, xs, ys, zs 
   if (atom == "O") then
    NSO=NSO+1
    xSO(NSO)=xs
    ySO(NSO)=ys
    zSO(NSO)=zs
   endif 
 enddo
   
!HB between water to solid
 WS=0
 do i=1, NWO
                                                                 ! with O on PEG
  do j=1,NSO
   xdiff=xWO(i)-xSO(j)
   ydiff=yWO(i)-ySO(j)
   zdiff=zWO(i)-zSO(j)
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
   if (dist < cut .and. dist>0.1) then         ! criteria of HB distance
                                                                              ! H1 in L1-O on PEG
    xdiff=xSO(j)-xh1(i)
    ydiff=ySO(j)-yh1(i)
    zdiff=zSO(j)-zh1(i)
    xdiffintra=xWO(i)-xh1(i)
    ydiffintra=yWO(i)-yh1(i)
    zdiffintra=zWO(i)-zh1(i)
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
     WS=WS+1
    end if
                                                                              ! H2-O
    xdiff=xSO(j)-xh2(i)
    ydiff=ySO(j)-yh2(i)
    zdiff=zSO(j)-zh2(i)
    xdiffintra=xWO(i)-xh2(i)
    ydiffintra=yWO(i)-yh2(i)
    zdiffintra=zWO(i)-zh2(i)
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
     WS=WS+1
    end if

   end if !end criteria of HB distance
  end do !end j

 end do !end i
 WStot = WStot + WS
end do !end Nstep
NstepReal = step - 1
WStot=WStot/NstepReal

surface=a*b*0.02 ! surface of box
densityHsurface=WStot/surface !number of Hb with pegO /nm2


write(2,*) "*************************************************"
write(2,*) "Nb of tot Hb with PEG =", WStot
write(2,*) "Density of Hb on surface (#Hbs/nm2)=" ,densityHsurface
deallocate (xSO)
deallocate (ySO)
deallocate (zSO)
deallocate (xWO)
deallocate (yWO)
deallocate (zWO)
deallocate (xh1)
deallocate (yh1)
deallocate (zh1)
deallocate (xh2)
deallocate (yh2)
deallocate (zh2)
close(1)
close(2)
 end program
