program create_binder_istart

! Daria Evry 13-06-17


! ************ BE CAREFULL : Z and X are exchanged *********

implicit none

integer*8 NO, Nstep,NstepTOT,Istart,Nskip1,Nskip2,step10, ss
integer*8 Npoint
integer*8 :: step, Natom, i, j, k, m, s, t, uL0=0, uL1=0, uL2=0, uL3=0, dL0=0, dL1=0, dL2=0, dL3=0

!target= % of time spent by a molecule in a layer to be part of it 
real*8 target, targetL0
real*8 intL0, intL1, intL2
real*8 a ,b ,c
real*8 :: xo, yo, zo, xdiff, ydiff, zdiff, r, lzdiff, lznear, rmin
real*8 :: x1, y1, z1, x2, y2, z2
real*8, dimension(:), allocatable ::  L0, L1, L2, L3
integer*8, dimension (:), allocatable :: mL0, mL1, mL2, mL3, box
real*8, dimension(:), allocatable :: lx, ly, lz
real*8, dimension(:,:), allocatable :: x, y, z, xh1, yh1, zh1, xh2, yh2, zh2
character(3) :: atom, point
character(256) txt
integer*8 :: checkup=0, checkdown=0


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    READ INFO FROM BOXDATA
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
      do while (trim(txt).NE.'$LAYERS-LIMITS')
         read(11,*) txt
      enddo
      print*, txt
      read(11,*) intL0
      read(11,*) intL1
      read(11,*) intL2
close(11)
open(11,file='BOXDATA',form='formatted')
!     Inizialize the Flag
      txt = ' '
!     Loop untill the Flag
      do while (trim(txt).NE.'$TARGET')
         read(11,*) txt
      enddo
      print*, txt
      read(11,*) target
close(11)
open(11,file='BOXDATA',form='formatted')
!     Inizialize the Flag
      txt = ' '
!     Loop untill the Flag
      do while (trim(txt).NE.'$TARGETL0')
         read(11,*) txt
      enddo
      print*, txt
      read(11,*) targetL0
close(11)

write(*,*) '-------------END BOXDATA READING--------'


Nstep=NstepTOT-istart+1


! READ Npoint for the interface file
!-------------------------------------
open (2, FILE="interface.xyz")
  read(2,*) Npoint
close(2)

! Allocate the vectors
!-------------------------------------

allocate(mL0(NO))
allocate(mL1(NO))
allocate(mL2(NO))
allocate(mL3(NO))
allocate(box(NO))
allocate(lx(Npoint))
allocate(ly(Npoint))
allocate(lz(Npoint))

allocate(L0(NO))
allocate(L1(NO))
allocate(L2(NO))
allocate(L3(NO))
L0=0
L1=0
L2=0
L3=0 ! daria merda

allocate(x(NO,Nstep))
allocate(y(NO,Nstep))
allocate(z(NO,Nstep))
allocate(xh1(NO,Nstep))
allocate(yh1(NO,Nstep))
allocate(zh1(NO,Nstep))
allocate(xh2(NO,Nstep))
allocate(yh2(NO,Nstep))
allocate(zh2(NO,Nstep))


! Open the files
!-------------------------------------

! Input
open (1, FILE="pos_rebuilt.xyz")
open (2, FILE="interface.xyz")

! Output
open (11, FILE="binder_uL0.xyz")
open (12, FILE="binder_uL1.xyz")
open (13, FILE="binder_uL2.xyz")
open (14, FILE="binder_uL3.xyz")
open (21, FILE="binder_dL0.xyz")
open (22, FILE="binder_dL1.xyz")
open (23, FILE="binder_dL2.xyz")
open (24, FILE="binder_dL3.xyz")
open (99, FILE="layers_population")
open (15, FILE="binder_uL0L1.xyz")
open (25, FILE="binder_dL0L1.xyz")
open (16, FILE="binder_uL0L1L2.xyz")
open (26, FILE="binder_dL0L1L2.xyz")
open (17, FILE="binder_uL0L1L2L3.xyz")
open (27, FILE="binder_dL0L1L2L3.xyz")



! Skip the Istart steps
!-------------------------

!  pos_rebuilt.xyz
   Nskip1=(Istart-1)*(NO*3+2)
   do t=1, Nskip1
      read(1,*)
   enddo

!  Interface.xyz
   Nskip2=(Istart-1)*(Npoint+2) 
   do t=1, Nskip2/10
      read(2,*)
   enddo

do step = 1, Nstep/10
print*, "***** STEP", step, " *****"
                                     ! reading atomic coordinate
 read(2,*)
 read(2,*)
 do i=1, Npoint
  read(2,*) point, lx(i), ly(i), lz(i)
 end do
 
 do step10 = 1, 10
 ss = step10 * step
 j=0
 read(1,*) Natom
 read(1,*)
 do while (j< NO)
  read(1,*) atom, xo, yo, zo
  if (atom == "O") then
   j = j+1
   x(j,ss)=xo
   y(j,ss)=yo
   z(j,ss)=zo
   read(1,*) atom, xh1(j,ss), yh1(j,ss), zh1(j,ss)
   read(1,*) atom, xh2(j,ss), yh2(j,ss), zh2(j,ss)
  end if
 end do
 if (j/=NO) then
  print*, "***** ERROR READING pos_rebuilt.xyz *****"
 end if

                           ! check distance of each  molecules from the interface
 do m=1, NO
  rmin = 100
  do i=1, Npoint
   xdiff = lx(i) - x(m,ss)
   ydiff = ly(i) - y(m,ss)
   zdiff = lz(i) - z(m,ss)
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
   r = SQRT(xdiff**2+ydiff**2+zdiff**2)
   if (r < rmin) then
    rmin = r
    lznear=lz(i)
   end if
  end do
  lzdiff= lznear-z(m,ss)
  if (lznear>=0) then
   box(m)=1
   checkup=checkup+1
   if (lzdiff<0) then
    rmin=-rmin
   end if
  end if
  if (lznear<0) then
   box(m)=2
   checkdown=checkdown+1
   if (lzdiff>0) then
    rmin=-rmin
   end if
  end if
                                                ! assign each water to its layer
  if (rmin <= intL0 ) then
   L0(m)=L0(m)+1
  else if (rmin > intL0 .and. rmin <= intL1 ) then
   L1(m)=L1(m)+1
  else if (rmin > intL1 .and. rmin <= intL2 ) then
   L2(m)=L2(m)+1
  else
   L3(m)=L3(m)+1
  end if

 end do
 end do ! end step10
end do ! end step

do m=1, NO

 L0(m)=L0(m)/Nstep
 L1(m)=L1(m)/Nstep
 L2(m)=L2(m)/Nstep
 L3(m)=L3(m)/Nstep
 
 if (box(m)==1) then
  if (L0(m)>=targetL0) then
    uL0=uL0+1    
    write(11,*) m, "1"
    write(12,*) m, "0"
    write(13,*) m, "0"
    write(14,*) m, "0"
    write(21,*) m, "0"
    write(22,*) m, "0"
    write(23,*) m, "0"
    write(24,*) m, "0"

    write(15,*) m, "1"
    write(25,*) m, "0"
    write(16,*) m, "1"
    write(26,*) m, "0"
    write(17,*) m, "1"
    write(27,*) m, "0"
  else if (L1(m)>=target) then
    uL1=uL1+1
    write(11,*) m, "0"
    write(12,*) m, "1"
    write(13,*) m, "0"
    write(14,*) m, "0"
    write(21,*) m, "0"
    write(22,*) m, "0"
    write(23,*) m, "0"
    write(24,*) m, "0"

    write(15,*) m, "1"
    write(25,*) m, "0"
    write(16,*) m, "1"
    write(26,*) m, "0"
    write(17,*) m, "1"
    write(27,*) m, "0"
  else if (L2(m)>=target) then
    uL2=uL2+1
    write(11,*) m, "0"
    write(12,*) m, "0"
    write(13,*) m, "1"
    write(14,*) m, "0"
    write(21,*) m, "0"
    write(22,*) m, "0"
    write(23,*) m, "0"
    write(24,*) m, "0"

    write(15,*) m, "0"
    write(25,*) m, "0"
    write(16,*) m, "1"
    write(26,*) m, "0"
    write(17,*) m, "1"
    write(27,*) m, "0"
  else if (L3(m)>=target) then
    uL3=uL3+1
    write(11,*) m, "0"
    write(12,*) m, "0"
    write(13,*) m, "0"
    write(14,*) m, "1"
    write(21,*) m, "0"
    write(22,*) m, "0"
    write(23,*) m, "0"
    write(24,*) m, "0"

    write(15,*) m, "0"
    write(25,*) m, "0"
    write(16,*) m, "0"
    write(26,*) m, "0"
    write(17,*) m, "1"
    write(27,*) m, "0"
  else
    write(11,*) m, "0"
    write(12,*) m, "0"
    write(13,*) m, "0"
    write(14,*) m, "0"
    write(21,*) m, "0"
    write(22,*) m, "0"
    write(23,*) m, "0"
    write(24,*) m, "0"

    write(15,*) m, "0"
    write(25,*) m, "0"
    write(16,*) m, "0"
    write(26,*) m, "0"
    write(17,*) m, "0"
    write(27,*) m, "0" 
  end if

 else if (box(m)==2) then
  if (L0(m)>=targetL0) then
    dL0=dL0+1
    write(21,*) m, "1"
    write(22,*) m, "0"
    write(23,*) m, "0"
    write(24,*) m, "0"
    write(11,*) m, "0"
    write(12,*) m, "0"
    write(13,*) m, "0"
    write(14,*) m, "0"

    write(15,*) m, "0"
    write(25,*) m, "1"
    write(16,*) m, "0"
    write(26,*) m, "1"
    write(17,*) m, "0"
    write(27,*) m, "1"
  else if (L1(m)>=target) then
    dL1=dL1+1
    write(21,*) m, "0"
    write(22,*) m, "1"
    write(23,*) m, "0"
    write(24,*) m, "0"
    write(11,*) m, "0"
    write(12,*) m, "0"
    write(13,*) m, "0"
    write(14,*) m, "0"

    write(15,*) m, "0"
    write(25,*) m, "1"
    write(16,*) m, "0"
    write(26,*) m, "1"
    write(17,*) m, "0"
    write(27,*) m, "1"
  else if (L2(m)>=target) then
    dL2=dL2+1
    write(21,*) m, "0"
    write(22,*) m, "0"
    write(23,*) m, "1"
    write(24,*) m, "0"
    write(11,*) m, "0"
    write(12,*) m, "0"
    write(13,*) m, "0"
    write(14,*) m, "0"

    write(15,*) m, "0"
    write(25,*) m, "0"
    write(16,*) m, "0"
    write(26,*) m, "1"
    write(17,*) m, "0"
    write(27,*) m, "1"
  else if (L3(m)>=target) then
    dL3=dL3+1
    write(21,*) m, "0"
    write(22,*) m, "0"
    write(23,*) m, "0"
    write(24,*) m, "1"
    write(11,*) m, "0"
    write(12,*) m, "0"
    write(13,*) m, "0"
    write(14,*) m, "0"

    write(15,*) m, "0"
    write(25,*) m, "0"
    write(16,*) m, "0"
    write(26,*) m, "0"
    write(17,*) m, "0"
    write(27,*) m, "1"
  else
    write(11,*) m, "0"
    write(12,*) m, "0"
    write(13,*) m, "0"
    write(14,*) m, "0"
    write(21,*) m, "0"
    write(22,*) m, "0"
    write(23,*) m, "0"
    write(24,*) m, "0"

    write(15,*) m, "0"
    write(25,*) m, "0"
    write(16,*) m, "0"
    write(26,*) m, "0"
    write(17,*) m, "0"
    write(27,*) m, "0"
  end if
 end if
end do

checkup=checkup/Nstep
checkdown=checkdown/Nstep

write(99,*) "up",checkup , "down",checkdown
write(99,*) "uL0 =",uL0 ,"|", "dL0 =",dL0
write(99,*) "uL1 =",uL1 ,"|", "dL1 =",dL1
write(99,*) "uL2 =",uL2 ,"|", "dL2 =",dL2
write(99,*) "uL3 =",uL3 ,"|", "dL3 =",dL3

deallocate(mL0)
deallocate(mL1)
deallocate(mL2)
deallocate(mL3)
deallocate(box)
deallocate(lx)
deallocate(ly)
deallocate(lz)

deallocate(L0)
deallocate(L1)
deallocate(L2)
deallocate(L3)

deallocate(x)
deallocate(y)
deallocate(z)
deallocate(xh1)
deallocate(yh1)
deallocate(zh1)
deallocate(xh2)
deallocate(yh2)
deallocate(zh2)


stop
end program
