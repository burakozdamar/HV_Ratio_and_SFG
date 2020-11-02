program make_2D_graph_for_HB_network

implicit none

! x=dist ; y=cos ; z=num.HB

real*8 :: x, y, bin, tol, xval, yval
real*8 :: xmax, xvmax, xmin, xvmin
real*8 :: ymax, yvmax, ymin, yvmin
integer*8 :: IO, i, j, k, nx, ny
real*8, dimension(:,:), allocatable :: z

integer narg

character(len=256) input 
character(len=256) output 

!!!!!!!!!!!!!!!!!!
bin=0.02
tol=bin/2
!!!!!!!!!!!!!!!!!!
xmax=3.2
xmin=2.3
nx=45
!!!!!!!!!!!!!!!!!!
ymax=1.d0
ymin=-1.d0
ny=100
!!!!!!!!!!!!!!!!!!
allocate(z(nx,ny))
!!!!!!!!!!!!!!!!!!

narg=command_argument_count()

if (narg .eq. 2) then
   call get_command_argument(1,input)
   call get_command_argument(2,output)
else
   write (*,*) 
   write (*,*) 'Usage:./make_2D_graph_for_HB_network input output'
   write (*,*) '----------------------------------------'
   write (*,*) 
   stop
endif

!-------------------------------------------------------------
open (1, FILE=input)
open (2, FILE=output)

k=0
do 
 read(1,*,IOSTAT=IO) x, y
 if (IO /= 0) then
  exit
 end if 
 do i=1,nx
  xvmin=xmin+(i-1)*bin
  xvmax=xmin+(i*bin) 
  do j=1,ny
   yvmin=ymin+(j-1)*bin    
   yvmax=ymin+(j*bin)
   
   if (x>=xvmin .and. x<xvmax .and. y>=yvmin .and. y<yvmax) then
    z(i,j)=z(i,j)+1
   end if

  end do
 end do

 k=k+1
 print*, k
end do 

do i=1,nx
 xval=xmin+(i*bin)-tol
 do j=1,ny
  yval=ymin+(j*bin)-tol
  z(i,j)=z(i,j)/k
  write(2,*) xval, yval, z(i,j)
 end do
end do

deallocate(z)

stop
end program
