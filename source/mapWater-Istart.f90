program mappingWater2WCinterface
!!!2020.07.27
!!!2020.08.27 add gaussian formula
!!!2020.08.31 from C/O(SAM) to water
!!!Wanlin
!!!Palaiseau-ECLA
implicit none
character(10) :: type, p, txt
Integer*8 :: t, tt, i, j, k, h, io,io2, Natom,Npoint, NC, NO, NIdown, NWO
Integer*8 :: nlx, nly,nbgrid
Integer*8 :: NstepTOT,Nstep,Istart,Nskip1,Nskip2,stepp
real*8 :: dlx, dly, ai, bi
real*8 :: x, y, z, time, xI, yI, zI
real*8 :: a=13.386d0, b=13.286d0, c=85.d0, cutoff=3.d0, E=1.8d0, pi=3.14159265359
!!!cutoff is the criteria of the distance between inter-facial C/O(SAM) to the WC interface
real*8 :: xdiffC, ydiffC, zdiffC, xdiffO, ydiffO, zdiffO
real*8 :: rC, rO, rCmin, rOmin, rWO
real*8 :: NCtop, NOtop, NCtop_Mean, NOtop_Mean
real*8, dimension(:), allocatable :: xC, yC, zC, xO, yO, zO,xWO,yWO,zWO
real*8, dimension(:), allocatable :: xIdown, yIdown, zIdown
Real*8, dimension(:), allocatable :: lx1, lx2, ly1, ly2
Real*8, dimension(:,:), allocatable :: NCtopGrid, NOtopGrid, NCtopGrid_Mean, NOtopGrid_Mean


allocate(xIdown(729))
allocate(yIdown(729))
allocate(zIdown(729))
allocate(xC(128), yC(128), zC(128))
allocate(xO(300), yO(300), zO(300))
allocate(xWO(120), yWO(120), zWO(120))
!open(1,file="pos_rebuilt_solid.xyz")
open(3,file="Water3-Gaussian-Map.dat")
!open(4,file="O-3-Gaussian-Map.dat")

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
      txt = ' '
      do while (trim(txt).NE.'$ISTART')
          read(11,*) txt
      enddo
      !print*, txt
      read(11,*) Istart
close(11)
open (5, FILE="interface.xyz")
  read(5,*) Npoint
close(5)
open (6, FILE="pos_rebuilt_solid.xyz")
  read(6,*) Natom
close(6)
Nstep=NstepTOT-Istart+1
nly = 13
nlx = 13
dlx = a/nlx
dly = b/nly
ai = -a/2
bi = b/2




allocate(NCtopGrid(nlx,nly), NOtopGrid(nlx,nly))
allocate(NCtopGrid_Mean(nlx,nly), NOtopGrid_Mean(nlx,nly))
allocate(lx1(nlx), lx2(nlx), ly1(nly), ly2(nly))

do i = 1, nlx
    lx2(i) = ai + dlx * i 
    lx1(i) = lx2(i) -dlx
enddo

do i = 1, nly
    ly2(i) = bi - dly*i
    ly1(i) = ly2(i) + dly
enddo

open(1,file="pos_rebuilt_solid.xyz")
open(2,file="interface.xyz")
NCtop = 0
NOtop = 0
do i = 1, nlx
    do j = 1, nly
        NCtopGrid(i,j) = 0
        NOtopGrid(i,j) = 0
    enddo !end j 
enddo ! end i

Nskip1=(Istart-1)*(Natom+2)
do t =1, Nskip1
    read(1,*)
enddo

Nskip2=(Istart-1)*(Npoint+2)
do t =1, Nskip2/10
    read(2,*)
enddo



do t = 1,Nstep/10
    NIdown = 0
    print*, t
    !!!Read interface files, get downpart WC interface coordinates 
    read(2,*,iostat=io2) Npoint
    if (io2 /= 0) then 
        stop
    endif
    read(2,*)
    do i = 1, Npoint
        read(2,*) p, xI, yI, zI
        if (zI <= 0) then
            NIdown = NIdown  + 1
            xIdown(NIdown) = xI
            yIdown(NIdown) = yI
            zIdown(NIdown) = zI
        endif
    enddo
    !!! pos_rebuild files are 10 times larger than interface files
    do tt = 1, 10 
        NC = 0
        NO = 0
        NWO = 0
        read(1,*,iostat=io) Natom
        if (io /= 0)  then
            stop
        endif
        read(1,*)
        !!!Read water molecules
        do i = 1,360
            read(1,*) type, x, y, z
            if (type == "O") then
            	NWO = NWO + 1
            	xWO(NWO) = x
            	yWO(NWO) = y
            	zWO(NWO) = z
            endif
        enddo
        !!!Read solid molecules and save the coord of C and O
        do i = 361, Natom
            read(1,*)
            !if (type == "C") then
                !NC=NC+1
                !xC(NC) = x
                !yC(NC) = y
                !zC(NC) = z
            !endif
            !if (type == "O") then
            !    NO=NO+1
            !    xO(NO) = x
            !    yO(NO) = y
            !    zO(NO) = z
            !endif
        enddo

    !!!Calculate and find the closest distance of each C and O towards the interface
    do i = 1, NWO
    rCmin = 7.d0
        do j = 1, NIdown
            xdiffC = xWO(i) - xIdown(j)
            ydiffC = yWO(i) - yIdown(j)
            zdiffC = zWO(i) - zIdown(j)
            call pbc(xdiffC, a)
            call pbc(ydiffC, b)
            call pbc(zdiffC, c)
            rWO = sqrt(xdiffC**2+ydiffC**2+zdiffC**2)
            if (rWO <= rCmin) then 
                rCmin = rWO
            endif
        enddo !end j
            if (rCmin <= cutoff) then
                !!!NCtop means nb of carbon atoms within 5Å with respect to WC interface
                NCtop = NCtop + 1
                do k = 1, nlx
                    if (xWO(i) <= lx2(k) .and. xWO(i) >= (lx2(k)-dlx)) then
                        do h = 1, nly
                            if (yWO(i) <= (ly2(h)+dly) .and. yWO(i) >= ly2(h)) then
                                NCtopGrid(k,h) = NCtopGrid(k,h) + exp(-rCmin**2/(2*E**2))/((2*pi*E**2)**1.5)
                                exit
                            endif
                        enddo
                    endif
                enddo
                !xCtop(NCtop) = xC(i)
                !yCtop(NCtop) = yC(i)
                !zCtop(NCtop) = zC(i) 
            endif
    enddo ! end i=NC
    ! do i = 1, NO
    ! rOmin = 11.d0
    !     do j = 1, NIdown
    !         xdiffO = xO(i) - xIdown(j)
    !         ydiffO = yO(i) - yIdown(j)
    !         zdiffO = zO(i) - zIdown(j)
    !         call pbc(xdiffO, a)
    !         call pbc(ydiffO, b)
    !         call pbc(zdiffO, c)
    !         rO = sqrt(xdiffO**2+ydiffO**2+zdiffO**2)
    !         if (rO <= rOmin) then 
    !             rOmin = rO
    !         endif
    !     enddo !end j
    !         if (rOmin < cutoff) then
    !             !!!NOtop means nb of O atoms within 5Å with respect to WC interface
    !             NOtop = NOtop + 1
    !             do k = 1, nlx
    !                 if (xO(i) <= lx2(k) .and. xO(i) >= (lx2(k)-dlx)) then
    !                     do h = 1, nly
    !                         if (yO(i) <= (ly2(h)+dly) .and. yO(i) >= ly2(h)) then
    !                             NOtopGrid(k,h) = NOtopGrid(k,h) + exp(-rOmin**2/(2*E**2))/((2*pi*E**2)**1.5)
    !                             exit
    !                         endif
    !                     enddo
    !                 endif
    !             enddo
    !             !xOtop(NOtop) = xO(i)
    !             !yOtop(NOtop) = yO(i)
    !             !zOtop(NOtop) = zO(i) 
    !         endif
    ! enddo ! end i=NO

    enddo ! end tt
enddo ! end t

     NCtop_Mean = NCtop/Nstep
     !NOtop_Mean = NOtop/Nstep
    print*, NCtop_Mean, NOtop_Mean

    do k = 1, nlx
        do h = 1, nly
            NCtopGrid_Mean(k,h) = NCtopGrid(k,h)/Nstep
            !NOtopGrid_Mean(k,h) = NOtopGrid(k,h)/Nstep
        enddo
    enddo

    do k = 1, nly
        write(3,FMT="(13F9.6)") NCtopGrid_Mean(:,k) 
        !write(4,FMT="(13F9.6)") NOtopGrid_Mean(:,k)
    enddo

close(1)
close(2)
close(3)
!close(4)
deallocate(xIdown)
deallocate(yIdown)
deallocate(zIdown)
deallocate(xC, yC, zC)
deallocate(xO, yO, zO)
deallocate(NCtopGrid, NOtopGrid)
deallocate(NCtopGrid_Mean, NOtopGrid_Mean)
deallocate(lx1, lx2, ly1, ly2)

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
