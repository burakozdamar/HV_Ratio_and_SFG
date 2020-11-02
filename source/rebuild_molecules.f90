program rebuilt_molecules

! if you have ions modify ani or/and cat
! you have to modify the input file (the one it is read by the program --> channel 1)

implicit none

integer*8 :: Nstep, step, i, jo, jh, Natom, Nsolid, NO, hcount,  j

!----------------------------------
real*8, parameter :: hcutoff=1.4
!-----------------------------------------------------

real*8  a, b, c
real*8 :: xfake, yfake, zfake, rmin1, rmin2, xmin1, xmin2, ymin1, ymin2,  zmin1, zmin2
real*8 :: xvfake, yvfake, zvfake, xmin1v, xmin2v, ymin1v, ymin2v,  zmin1v, zmin2v
real*8 :: r, xdiff, ydiff, zdiff, rmin3, rmin4, xmin3, xmin4, ymin3, ymin4,  zmin3, zmin4
character(3) :: atom, atomv
real*8, dimension (:), allocatable ::  xo, yo, zo
real*8, dimension (:), allocatable ::  xov, yov, zov
real*8, dimension (:), allocatable ::  xh, yh, zh
real*8, dimension (:), allocatable ::  xhv, yhv, zhv

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
      do while (trim(txt).NE.'$NATOM')
         read(11,*) txt
      enddo
      print*, txt
      read(11,*) Natom
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
      do while (trim(txt).NE.'$NSTEP')
         read(11,*) txt
      enddo
      print*, txt
      read(11,*) Nstep
close(11)


print*, '*********END BOXDATA READING******'

allocate (xo(NO))
allocate (yo(NO))
allocate (zo(NO))
allocate (xov(NO))
allocate (yov(NO))
allocate (zov(NO))
allocate (xh(NO*2))
allocate (yh(NO*2))
allocate (zh(NO*2))
allocate (xhv(NO*2))
allocate (yhv(NO*2))
allocate (zhv(NO*2))

!Input
!------------------------------
open (1, FILE="pos_recentered.xyz")
open (2, FILE="vel.xyz")

! Output
!-----------------
open (3, FILE="pos_mol_rebuilt.xyz")
open (4, FILE="vel_mol_rebuilt.xyz")



do step = 1, Nstep

!Inizialize the counter
 jo=0
 jh=0

!--Read--

 read(1,*)       !pos
 read(1,*)       !pos

 read(2,*)       !vel
 read(2,*)       !vel

!--Write--


! pos
 write(3,*) Natom
 write(3,*) "step ", step
! vel
 write(4,*) Natom
 write(4,*) "step ", step


!!!!!!!!!!!!!!!!!  SKIP SOLID+IONS !!!!!!!!!!!!!!!!!!!!!
 do i=1,Nsolid
! pos
  read(1,*)   atom, xfake, yfake, zfake
  write(3,*)  atom, xfake, yfake, zfake
! vel
  read(2,*)   atomv, xvfake, yvfake, zvfake
  write(4,*)  atomv, xvfake, yvfake, zvfake
 end do !i
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 do i=1,NO
! pos
  read(1,*) atom, xfake, yfake, zfake
! vel
  read(2,*) atomv, xvfake, yvfake, zvfake
   jo=jo+1 ! update counter
   xo(jo)=xfake
   yo(jo)=yfake
   zo(jo)=zfake
   xov(jo)=xvfake
   yov(jo)=yvfake
   zov(jo)=zvfake
 end do
 do i=1,NO*2 
! pos
  read(1,*) atom, xfake, yfake, zfake
! vel
  read(2,*) atomv, xvfake, yvfake, zvfake
   jh=jh+1 ! update counter
   xh(jh)=xfake
   yh(jh)=yfake
   zh(jh)=zfake
   xhv(jh)=xvfake
   yhv(jh)=yvfake
   zhv(jh)=zvfake
 end do ! i
 if (NO*3 /= (jo+jh)) then
  print*, "********** ERROR READING INPUT FILE **********"
  print*, "********** Problem with H or O  **********"
  stop
 end if

!----------------------------------
! molecules reconstruction
!------------------------------------
 do i=1,jo

  hcount=0

  rmin1=2.0
  rmin2=2.0
  
  do j=1,jh
   xdiff=xo(i)-xh(j)
   ydiff=yo(i)-yh(j)
   zdiff=zo(i)-zh(j)
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

   if ( zdiff > c/2) then
    zdiff=zdiff-c
   else if ( zdiff < -c/2) then
    zdiff=zdiff+c
   end if
                                     ! end pbc

   r=SQRT(xdiff**2+ydiff**2+zdiff**2)
   
   if (r<=hcutoff) then
    hcount=hcount+1
    if (r<rmin1) then
     rmin2=rmin1

!    pos
     xmin2=xmin1
     ymin2=ymin1
     zmin2=zmin1
!    vel
     xmin2v=xmin1v
     ymin2v=ymin1v
     zmin2v=zmin1v

     rmin1=r

!    pos
     xmin1=xh(j)
     ymin1=yh(j)
     zmin1=zh(j)
!    vel
     xmin1v=xhv(j)
     ymin1v=yhv(j)
     zmin1v=zhv(j)

    else if (r<rmin2) then

     rmin2=r

!    pos
     xmin2=xh(j)
     ymin2=yh(j)
     zmin2=zh(j)
!    vel
     xmin2v=xhv(j)
     ymin2v=yhv(j)
     zmin2v=zhv(j)

    end if
   end if
  end do !j
  if (hcount<2) then
   print*, "not enough H ; mol :", i, "step", step
  end if


! pos
  write(3,*) "O", xo(i), yo(i), zo(i)
  write(3,*) "H", xmin1, ymin1, zmin1
  write(3,*) "H", xmin2, ymin2, zmin2

! vel
  write(4,*) "O", xov(i), yov(i), zov(i)
  write(4,*) "H", xmin1v, ymin1v, zmin1v
  write(4,*) "H", xmin2v, ymin2v, zmin2v
 end do !i
end do !t

deallocate (xo)
deallocate (yo)
deallocate (zo)
deallocate (xov)
deallocate (yov)
deallocate (zov)
deallocate (xh)
deallocate (yh)
deallocate (zh)
deallocate (xhv)
deallocate (yhv)
deallocate (zhv)

end program
