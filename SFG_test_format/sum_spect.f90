program sum_up_down

implicit none
integer*8 :: i
real*8 :: spectra, spectrau, spectrad, skip
real*8 :: spe2, spe2u, spe2d
integer narg

character(len=256) input1
character(len=256) input2
character(len=256) output


narg=command_argument_count()

if (narg .eq. 3) then
   call get_command_argument(1,input1)
   call get_command_argument(2,input2)
   call get_command_argument(3,output)
else
 write(*,*) 'Usage: sum input1 input2 output'
 write(*,*) '----------------------------------------'
 write(*,*) 
 stop
endif


open(1,FILE=input1)
open(2,FILE=input2)
open(3,FILE=output)

do i=1,4000
 read(1,*) skip, spectrau, spe2u
 read(2,*) skip, spectrad, spe2d
 spectra=spectrau+spectrad
 spe2= spe2u +spe2d
 write(3,*) i, spectra, spe2
end do

stop
end program
