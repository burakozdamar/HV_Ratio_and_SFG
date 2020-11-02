program deconvolve
! only for down L1

implicit none

integer*8 :: t, Nstep,NstepTOT, NO, Natom,istart, i, j, Nskip1
real*8 :: dz1, dz2, zz, a, b, c
real*8 :: xo,  yo,  zo, xh1, yh1, zh1, xh2, yh2, zh2
real*8, dimension(:), allocatable :: totti, flag, dzeko
character(3) :: atom

character(len=256) txt

!---------------------------------------


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

Nstep=NstepTOT-Istart+1

write(*,*) '-------------END BOXDATA READING--------'

open (1, FILE="pos_rebuilt.xyz")
open (10,file="binder_dL1.xyz")
open (2, FILE="binder_up.xyz")
open (3, FILE="binder_down.xyz")
open (4, FILE="BIL_populations")


allocate(totti(NO))
allocate(dzeko(NO))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! reading the binder file (to select the mol of the interface)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
allocate(flag(NO))
do i=1,NO
 read(10,*) j, flag(i)
end do
close(10)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Skip the Istart steps
!-------------------------

!  pos_rebuilt.xyz
   Nskip1=(Istart-1)*(NO*3+2)
   do t=1, Nskip1
      read(1,*)
   enddo

totti=0.d0
dzeko=0.d0
do t=1,Nstep
 print*, "t",t
 read(1,*) Natom
 read(1,*)
 do i=1,NO
  read(1,*) atom, xo,  yo,  zo
  read(1,*) atom, xh1, yh1, zh1
  read(1,*) atom, xh2, yh2, zh2
  dz1=zh1-zo
  dz2=zh2-zo
  zz=dz1+dz2
  if (flag(i)==1) then
   if (zz>=0) then
    totti(i)=totti(i)+1
   else
    dzeko(i)=dzeko(i)+1
   endif
  end if
 end do
end do

a=0
b=0
c=0
do i=1,NO
 totti(i)=totti(i)/Nstep
 dzeko(i)=dzeko(i)/Nstep
 if (flag(i)==1) then
  a=a+1
  if  (totti(i)>0.7) then
   b=b+1
   write(2,*) i, "0"
   write(3,*) i, "1"
  else !if (dzeko(>0.5) then
   c=c+1
   write(2,*) i, "1"
   write(3,*) i, "0"
!  else
!   write(2,*) i, "0"
!   write(3,*) i, "0"
  end if
 else 
  write(2,*) i, "0"
  write(3,*) i, "0"
 end if
end do

write(4,*) "tot", a
write(4,*) "down", b
write(4,*) "up", c

stop
end program
