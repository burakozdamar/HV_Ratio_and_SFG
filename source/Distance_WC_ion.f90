program Distance_ion_WC
!!! 2020.11.10 Evry
!!! Onelin
implicit none

integer*8 :: Natom, NO, Nion, i, j, Nstep, io, io2
integer*8 :: xgrid, ygrid, t, tt
real*8 fake1, fake2 !Dummy variables for the grid reading
real*8 :: a, b, c
real*8 :: x, y, z, xk, yk, zk, xcl, ycl, zcl
real*8 :: xdiffK, ydiffK, zdiffK, xdiffCl, ydiffCl, zdiffCl, xdiffIon, ydiffIon, zdiffIon
real*8 :: dK, dKmin, dCl, dClmin, dKCl, time
character(20) :: atom, out1, txt

real*8, dimension (:,:), allocatable :: ulx, uly, ulz, dlx, dly, dlz



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
      read(11,*) Nstep
close(11)

open(12,file='grid_interface',form='formatted')
read(12,*) !Skipe line
read(12,*) !Skipe line
read(12,*) xgrid, fake1, fake2
read(12,*) ygrid, fake1, fake2
close(12)

allocate (ulx(xgrid,ygrid))
allocate (uly(xgrid,ygrid))
allocate (ulz(xgrid,ygrid))
allocate (dlx(xgrid,ygrid))
allocate (dly(xgrid,ygrid))
allocate (dlz(xgrid,ygrid))

open (1, FILE="pos_rebuilt_solid.xyz")
open (2, FILE="interface.xyz")

out1 = "K-Cl-WC-dis.dat"
open (3, FILE= out1)

write(3,*) "# Dkmin, DClmin, dKCl"

do t = 1, Nstep/10
	write(6,*) t
	read(2,*,iostat=io2)
		if (io2 /= 0) then
			exit
		endif
	read(2,*)
	do i = 1, xgrid
		do j= 1, ygrid
			read(2,*) atom, ulx(i,j), uly(i,j), ulz(i,j)
			read(2,*) atom, dlx(i,j), dly(i,j), dlz(i,j)
		enddo
	enddo

	do tt = 1, 10
		read(1,*,iostat=io) Natom
	        if (io /= 0)  then
	            exit
	        endif
        read(1,*)
        !!!Read water molecules + solid
        do i = 1, (Natom - Nion)
            read(1,*)
        enddo
        !!!Read ions
        do i = 1, Nion
            read(1,*) atom, x, y, z
            if (atom == "K") then
                xk = x
                yk = y
                zk = z
            endif
            if (atom == "Cl") then
                xcl= x
                ycl = y
                zcl= z
            endif
        enddo
	     !!! distance between ions
	    xdiffIon=xcl-xk
	    ydiffIon=ycl-yk
	    zdiffIon=zcl-zk
	    call pbc(xdiffIon, a)
	    call pbc(ydiffIon, b)
	    call pbc(zdiffIon, c)

	    dKCl=(xdiffIon**2+ydiffIon**2+zdiffIon**2)**0.5


	  	dKmin=100
	  	dClmin=100
	  	do i=1,xgrid
		   	do j=1,ygrid
		                                              ! up
			    xdiffK=dlx(i,j)-xk
			    ydiffK=dly(i,j)-yk
			    zdiffK=dlz(i,j)-zk

			    xdiffCl=dlx(i,j)-xcl
			    ydiffCl=dly(i,j)-ycl
			    zdiffCl=dlz(i,j)-zcl

		        call pbc(xdiffK, a)
		        call pbc(ydiffK, b)
		        call pbc(zdiffK, c)

		        call pbc(xdiffCl, a)
		        call pbc(ydiffCl, b)
		        call pbc(zdiffCl, c)

		        dK = sqrt(xdiffK**2+ydiffK**2+zdiffK**2)
		        dCl = sqrt(xdiffCl**2+ydiffCl**2+zdiffCl**2)

		        if (dK <= dKmin) then
		        	dKmin = dK
		        endif
		        if (dCl <= dClmin) then
		        	dClmin = dCl
		        endif
		     enddo !end j
		enddo !end i
		time = ((t-1)*10+tt)*0.4/1000
	    write(3,FMT="(F6.2, F7.3, F7.3, F7.3)") time,  dkmin,  dclmin,  dKCl
    enddo ! end tt 10 steps
enddo ! end t interface step

close(1)
close(2)
close(3)

deallocate (ulx)
deallocate (uly)
deallocate (ulz)
deallocate (dlx)
deallocate (dly)
deallocate (dlz)

end program

subroutine pbc(coord,lattice)
real*8 :: coord, lattice
if (coord > lattice/2) then
    coord = coord - lattice
else if (coord < -lattice/2) then
    coord = coord + lattice
endif
return
end
