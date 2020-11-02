program up_down_water_layers

implicit none

integer*8 npoint, Nlay, io1, io2
real*8, parameter :: pi=3.14159265359
real*8 :: a, b, c, xo, yo, zo, xdiff, ydiff, zdiff, r, p, lzdiff, lznear, rmin
real*8, dimension(:), allocatable ::  intL
real*8, dimension (:), allocatable :: x, y, z, lx, ly, lz
real*8, dimension (:), allocatable :: xh1, yh1, zh1, xh2, yh2, zh2
real*8, dimension (:,:), allocatable ::  uxL, uyL, uzL, uxh1L, uyh1L, uzh1L, uxh2L, uyh2L, uzh2L
real*8, dimension (:,:), allocatable ::  dxL, dyL, dzL, dxh1L, dyh1L, dzh1L, dxh2L, dyh2L, dzh2L
integer*8, dimension(:), allocatable :: uuL, duL
integer*8 :: step, stepp,  Nstep, Natom, i, j, k, NO, m, s
integer*8 :: h1i, h1j, h2i, h2j, Ntot
character(3) :: atom, point
real*8, dimension (:), allocatable :: box   ! 0=up ; 1=down
integer*8 ::t, NstepTOT, Istart, Nskip1, Nskip2

!BOXDATA reading variable
!working flag for getting info
integer narg
character(len=256) arg
character(256)  txt


character(len=1) str
character(len=5) out1,out2,out3
 
!#################################
! STARTING OF THE PROGRAM
!######################################

!-------------------------------------
narg=command_argument_count()

if (narg .eq. 1) then
   call get_command_argument(1,arg)
   read(arg,'(i10)') Nlay
else
   write (*,*) 
   write (*,*) 'Usage: ./up_down_water_layers Nlayers   '
   write (*,*) '----------------------------------------'
   write (*,*) ' Nlayers  = Number of layers            '              
   write (*,*) '----------------------------------------'
   write (*,*) 
   stop
endif

Nlay=Nlay+1

allocate(intL(Nlay))




!-------------------------------------------------------------

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

open(11,file='BOXDATA',form='formatted')
!     Inizialize the Flag
      txt = ' ' 
!     Loop untill the Flag
      do while (trim(txt).NE.'$LAYERS-LIMITS')
         read(11,*) txt
      enddo
      intL(1) = -c/2
      do i=2, Nlay-1
         read(11,*) intL(i)
      enddo
      intL(Nlay) = c
close(11)




!---------------------------------------------------

!read the number of points of the interface
open (2, FILE="interface.xyz")
  read(2,*) Npoint
close(2)

      
 

open (1, FILE="pos_rebuilt.xyz")
open (2, FILE="interface.xyz")
open (4, FILE="up_down_water_layers.xyz")

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
allocate (box(NO))
allocate (x(NO))
allocate (y(NO))
allocate (z(NO))
allocate (xh1(NO))
allocate (yh1(NO))
allocate (zh1(NO))
allocate (xh2(NO))
allocate (yh2(NO))
allocate (zh2(NO))
                        !!! up !!!
allocate (uuL(nLay))
allocate (uxL(NO,Nlay)) 
allocate (uyL(NO,Nlay)) 
allocate (uzL(NO,Nlay)) 
allocate (uxh1L(NO,NLay))
allocate (uyh1L(NO,NLay))
allocate (uzh1L(NO,NLay))
allocate (uxh2L(NO,NLay))
allocate (uyh2L(NO,NLay))
allocate (uzh2L(NO,NLay))

                           !!! down !!!
allocate (duL(nLay))
allocate (dxL(NO,Nlay)) 
allocate (dyL(NO,Nlay)) 
allocate (dzL(NO,Nlay)) 
allocate (dxh1L(NO,NLay))
allocate (dyh1L(NO,NLay))
allocate (dzh1L(NO,NLay))
allocate (dxh2L(NO,NLay))
allocate (dyh2L(NO,NLay))
allocate (dzh2L(NO,NLay))


!!!!!!!!!!!!!!!!!!!!!!!!
allocate (lx(Npoint))
allocate (ly(Npoint))
allocate (lz(Npoint))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! skip pos_rebuilt.xyz
Nskip1=(Istart-1)*(NO*3+2)
do t = 1, Nskip1
    read(1,*)
enddo

Nskip2=(Istart-1)*(Npoint+2)
do t=1, Nskip2/10
    read(2,*)
enddo

do step = 1, Nstep/10
 !s = 1 + (step-1)*10
 read(2,*,iostat=io2)
    if (io2 /= 0) then
        stop
    endif
 read(2,*)
 do i=1,Npoint
  read(2,*) point, lx(i), ly(i), lz(i) 
 end do

 do stepp = 1, 10
    s = (step-1)*10 + stepp
!Inizialize the counters   
 do k=1,Nlay
  uuL(k)=0
  duL(k)=0
 enddo


                                     ! reading atomic coordinate
 j=1
 read(1,*,iostat=io1) Natom
    if (io1 /= 0) then
        stop
    endif
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



! check distance of each  molecules from the interface
 do m=1, NO
  rmin = 100
  do i=1, Npoint
   xdiff = lx(i) - x(m)
   ydiff = ly(i) - y(m)
   zdiff = lz(i) - z(m)
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
  lzdiff= lznear-z(m)
  if (lznear>0) then  
   box(m)=0
   if (lzdiff<0) then
    rmin=-rmin
   end if
  end if
  if (lznear<0) then
   box(m)=1
   if (lzdiff>0) then
    rmin=-rmin
   end if
  end if

                                                ! assign each water to its layer
                                   !!! up !!!
   if (box(m)==0) then
    do k=1,Nlay
     if (rmin > intL(k) .and. rmin <= intL(k+1) ) then
      uuL(k)=uuL(k)+1
      uxL(uuL(k),k) = x(m)
      uyL(uuL(k),k) = y(m)
      uzL(uuL(k),k) = z(m)
      uxh1L(uuL(k),k) = xh1(m)
      uyh1L(uuL(k),k) = yh1(m)
      uzh1L(uuL(k),k) = zh1(m)
      uxh2L(uuL(k),k) = xh2(m)
      uyh2L(uuL(k),k) = yh2(m)
      uzh2L(uuL(k),k) = zh2(m)
     end if
    enddo
   endif
!!! down !!!
   if (box(m)==1) then
    do k=1,Nlay
     if (rmin > intL(k) .and. rmin <= intL(k+1) ) then
      duL(k)=duL(k)+1
      dxL(duL(k),k) = x(m)
      dyL(duL(k),k) = y(m)
      dzL(duL(k),k) = z(m)
      dxh1L(duL(k),k) = xh1(m)
      dyh1L(duL(k),k) = yh1(m)
      dzh1L(duL(k),k) = zh1(m)
      dxh2L(duL(k),k) = xh2(m)
      dyh2L(duL(k),k) = yh2(m)
      dzh2L(duL(k),k) = zh2(m)
     endif
    enddo
   endif
 end do !NO

 write (4,*) Natom
 write (4,*) "step =", s
 
 do k=1,Nlay
   do i=1,uuL(k)
     write(str,'(i1.1)') k-1
     out1='uOL'//str
     out2='uH1L'//str
     out3='uH2L'//str
     write (4,FMT="(A , F14.6 , F14.6 , F14.6)") out1,    uxL(i,k),   uyL(i,k),   uzL(i,k)
     write (4,FMT="(A , F14.6 , F14.6 , F14.6)") out2,  uxh1L(i,k), uyh1L(i,k), uzh1L(i,k)
     write (4,FMT="(A , F14.6 , F14.6 , F14.6)") out3,  uxh2L(i,k), uyh2L(i,k), uzh2L(i,k)
   end do
 enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 do k=1,Nlay
   do i=1,duL(k)
     write(str,'(i1.1)') k-1
     out1='dOL'//str
     out2='dH1L'//str
     out3='dH2L'//str
     write (4,FMT="(A , F14.6 , F14.6 , F14.6)") out1,    dxL(i,k),   dyL(i,k),   dzL(i,k)
     write (4,FMT="(A , F14.6 , F14.6 , F14.6)") out2,  dxh1L(i,k), dyh1L(i,k), dzh1L(i,k)
     write (4,FMT="(A , F14.6 , F14.6 , F14.6)") out3,  dxh2L(i,k), dyh2L(i,k), dzh2L(i,k)
   end do
 enddo
 enddo ! end stepp (10)
end do ! time 

deallocate (box)
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
!!!!!!!!!!!!!!!!!!!!! up !!!!!!!!!!!!!!!!!!!!!!!
deallocate (uuL)
deallocate (uxL) 
deallocate (uyL) 
deallocate (uzL) 
deallocate (uxh1L)
deallocate (uyh1L)
deallocate (uzh1L)
deallocate (uxh2L)
deallocate (uyh2L)
deallocate (uzh2L)

                           !!! down !!!
deallocate (duL)
deallocate (dxL) 
deallocate (dyL) 
deallocate (dzL) 
deallocate (dxh1L)
deallocate (dyh1L)
deallocate (dzh1L)
deallocate (dxh2L)
deallocate (dyh2L)
deallocate (dzh2L)


stop
end program
