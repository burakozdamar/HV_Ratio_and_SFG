program rebuilt_box_KCl

implicit none

integer*8 NO
real*8 a, b, c, zval
integer*8 :: Nstep, step, i, j, Natom, Nsolid, Nion, ii
real*8 :: xo, yo, zo, xoh1, xoh2, yoh1, yoh2, zoh1, zoh2
real*8, dimension (:), allocatable :: x, y, z, xh1, yh1, zh1, xh2, yh2, zh2
character(3) :: atom
real*8 :: xi, yi, zi 

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
open(11,file='BOXDATA',form='formatted')
!     Inizialize the Flag
      txt = ' ' 
!     Loop untill the Flag
      do while (trim(txt).NE.'$NIONS')
         read(11,*) txt
      enddo
      print*, txt
      read(11,*) Nion
close(11)
open(11,file='BOXDATA',form='formatted')
!     Inizialize the Flag
      txt = ' ' 
!     Loop untill the Flag
      do while (trim(txt).NE.'$ZTRASL')
         read(11,*) txt
      enddo
      print*, txt
      read(11,*) zval
close(11)

print*, '********End BOXDATA reading***'

!-------------------------------------------
!-------------------------------------------

open (1, FILE="pos_mol_rebuilt.xyz")
open (2, FILE="pos_rebuilt_KCl.xyz")

allocate (x(NO))
allocate (y(NO))
allocate (z(NO))
allocate (xh1(NO))
allocate (yh1(NO))
allocate (zh1(NO))
allocate (xh2(NO))
allocate (yh2(NO))
allocate (zh2(NO))

Natom=NO*3
do step = 1, Nstep
 print*, step
                                         ! reading input file
 write(2,*) Natom + Nion
 write(2,*) "step = ", step

 j=1
! Skip 2 lines
 read(1,*)
 read(1,*)

! Skip solid lines
 do i=1,Nsolid
  read(1,*)
 end do

! Traslate + Write ions
 do i=1,Nion
   read(1,*) atom, xi, yi, zi
   zi=zi+zval
   if (zi < -c/2) then
     zi=zi+c
   else if (zi > c/2) then
     zi=zi-c
   end if
   write(2,*) atom, xi, yi, zi
 enddo

!-Water
 do while (j<= NO)
  read(1,*) atom, xo, yo, zo
  x(j)=xo
  y(j)=yo
  z(j)=zo
  if (atom == "O") then
   read(1,*) atom, xh1(j), yh1(j), zh1(j)
   read(1,*) atom, xh2(j), yh2(j), zh2(j)
   j = j+1
  end if
 end do


                                         ! box traslation   
 do i=1,NO
 z(i)=z(i)+zval
 zh1(i)=zh1(i)+zval
 zh2(i)=zh2(i)+zval
                                         ! centering water molecules in the box
 !!!!!!!!  pbc  !!!!!!!!!!!!
 do ii=1,2
 if (x(i) < -a/2) then
  x(i)=x(i)+a
 else if (x(i) > a/2) then
  x(i)=x(i)-a
 end if

 if (xh1(i) < -a/2) then
  xh1(i)=xh1(i)+a
 else if (xh1(i) > a/2) then
  xh1(i)=xh1(i)-a
 end if

 if (xh2(i) < -a/2) then
  xh2(i)=xh2(i)+a
 else if (xh2(i) > a/2) then
  xh2(i)=xh2(i)-a
 end if

 if (y(i) < -b/2) then
  y(i)=y(i)+b
 else if (y(i) > b/2) then
  y(i)=y(i)-b
 end if

 if (yh1(i) < -b/2) then
  yh1(i)=yh1(i)+b
 else if (yh1(i) > b/2) then
  yh1(i)=yh1(i)-b
 end if

 if (yh2(i) < -b/2) then
  yh2(i)=yh2(i)+b
 else if (yh2(i) > b/2) then
  yh2(i)=yh2(i)-b
 end if

 if (z(i) < -c/2) then
  z(i)=z(i)+c
 else if (z(i) > c/2) then
  z(i)=z(i)-c
 end if

 if (zh1(i) < -c/2) then
  zh1(i)=zh1(i)+c
 else if (zh1(i) > c/2) then
  zh1(i)=zh1(i)-c
 end if

 if (zh2(i) < -c/2) then
  zh2(i)=zh2(i)+c
 else if (zh2(i) > c/2) then
  zh2(i)=zh2(i)-c
 end if
 end do
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                                               ! rebuilt H2O molecules
  xoh1=xh1(i)-x(i)
  if (xoh1>a/2) then
   xh1(i)=xh1(i)-a
  else if (xoh1<-a/2) then
   xh1(i)=xh1(i)+a
  end if
  yoh1=yh1(i)-y(i)
  if (yoh1>b/2) then
   yh1(i)=yh1(i)-b
  else if (yoh1<-b/2) then
   yh1(i)=yh1(i)+b
  end if
  zoh1=zh1(i)-z(i)
  if (zoh1>c/2) then
   zh1(i)=zh1(i)-c
  else if (zoh1<-c/2) then
   zh1(i)=zh1(i)+c
  end if





  xoh2=xh2(i)-x(i)
  if (xoh2>a/2) then
   xh2(i)=xh2(i)-a
  else if (xoh2<-a/2) then
   xh2(i)=xh2(i)+a
  end if
  yoh2=yh2(i)-y(i)
  if (yoh2>b/2) then
   yh2(i)=yh2(i)-b
  else if (yoh2<-b/2) then
   yh2(i)=yh2(i)+b
  end if
  zoh2=zh2(i)-z(i)
  if (zoh2>c/2) then
   zh2(i)=zh2(i)-c
  else if (zoh2<-c/2) then
   zh2(i)=zh2(i)+c
  end if



  write(2,*) "O", x(i), y(i), z(i)
  write(2,*) "H", xh1(i), yh1(i), zh1(i)
  write(2,*) "H", xh2(i), yh2(i), zh2(i)


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

stop
end program

