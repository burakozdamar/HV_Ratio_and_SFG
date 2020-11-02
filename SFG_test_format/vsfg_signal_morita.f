      program vsfg_signal_morita

      implicit real*8(a-h,o-z)
      integer nmax,i, imax,imax2,j
      parameter(nmax=300000)
      real*8 corr(nmax),corr2(nmax)
      real*8 tf(0:nmax),tf2(0:nmax)
      real*8 filter,cel,pi,filter2,filter3
      real*8 coeff,coeff2,cotan,tang
      real*8 kt
      integer iimax
      real*8 imaginary_part,imaginary_part_1
      integer padd
 
      integer narg
      character(len=256) txt

      character(len=256) input
      character(len=256) output,output2

c############################################
c     STARTING OF THE PROGRAM
c############################################


      narg=command_argument_count()
      if (narg .eq. 2) then
        call get_command_argument(1,input)
        call get_command_argument(2,output)
      else
        write (*,*)
        write (*,*) 'Usage:./vsfg_signal_morita.x input output'
        write (*,*) '----------------------------------------'
        write (*,*)
        stop
      endif

!-------------------------------------------------------------



      open(11,file='BOXDATA',form='formatted')
!     Inizialize the Flag
        txt = ' '
!     Loop untill the Flag
      do while (trim(txt).NE.'$DT')
         read(11,*) txt
      enddo
      print*, txt
      read(11,*) dt
      close(11)
      open(11,file='BOXDATA',form='formatted')
!     Inizialize the Flag
        txt = ' '
!     Loop untill the Flag
      do while (trim(txt).NE.'$PADD')
         read(11,*) txt
      enddo
      print*, txt
      read(11,*) padd
      close(11)
      padd=padd+1
      open(11,file='BOXDATA',form='formatted')
!     Inizialize the Flag
        txt = ' '
!     Loop untill the Flag
      do while (trim(txt).NE.'$FREQ-DFREQ')
         read(11,*) txt
      enddo
      print*, txt
      read(11,*) wmax, dw !cm-1
      close(11)

      write(*,*) '-------------END BOXDATA READING--------'


! === Get the classical correlation function between polarisability and dipole
!     (from vsfg_correlation.f)


      open(10,file=input) 
      do i=1,nmax
         read(10,*,end=100)ii,corr(i)
      enddo
 100  continue
      imax=i-1
      write(6,*)'imax=',imax
      close(10)

!#########################################
c   ADD THE ZERO-PADDING
c   add a vector of zero at the end of the
c   correlation function to avoid the discontinuity
c   (and so improve the resolution of the spectra) 
c
c   padd=1 => don't use the ZEro padding 
c   correction

      padd = padd+1

c   Zero Padding 
      if (padd.gt.1) then
        do i = (imax+1), imax*padd
              corr(i) = 0.0D0
        enddo
      endif
      imax=imax*padd



! normalisation pour VDOS
!      do i=2,imax
!         corr(i)=corr(i)/corr(1)
!      enddo
!      corr(1)=1.d0


! ==== Fourier transform the classical correlation function

      do i=1,nmax
         tf(i)=0.d0
         tf2(i)=0.d0
      enddo

      num_tf=wmax/dw
      write(6,*)'num_tf=',num_tf
      if(wmax/dw >= nmax)then
         write(6,*)'!!! dimension table tf !!!'
         stop
      endif

      dt = dt*1.d-15 !s
      cel = 3.d10 !cm/s
      pi=dacos(-1.d0)
      imax2=imax/2

      do j=0,num_tf
         sigma=dble(j)*dw       !cm-1
         omega=2.d0*pi*cel*sigma !s-1
         do i=1,imax
            t=dt*dble(i)
            tf(j)=tf(j)+dcos(omega*t)*corr(i)*filter(i,imax2)
            tf2(j)=tf2(j)+dsin(omega*t)*corr(i)*filter(i,imax2)
         enddo
         tf(j)=tf(j)*2.d0 !because of DFT for symetrical function (see Allen p.336)
         tf2(j)=tf2(j)*2.d0 !because of DFT for symetrical function (see Allen p.336)
      enddo

      open(10,file='fft.dat')
      do j=0,num_tf
        sigma=dble(j)*dw       !cm-1
        write(10,*)sigma,tf(j)-tf(int(wmax/dw)),tf2(j)-tf2(int(wmax/dw))
      enddo
      close(10)

! ======= 

      planck = 1.05457267d-34      !J.s
      kt = 1.3806658d-23 * 300.0d0 !J 
      beta = 1.d0/kt
      
      open(10,file=output)
      output2='square'//trim(output)
      open(11,file=output2)

! Offsets to put the high frequency signal at zero 
      sigma=dble(num_tf)*dw           !cm-1
      omega=2.d0*pi*cel*sigma   !s-1
      real_part_1 = -(1.d0/2.d0)/omega*beta*tf2(num_tf) 
      imaginary_part_1 = (1.d0/2.d0)/omega*beta*tf(num_tf) 
      write(6,*)'Extreme values at num_tf :',real_part_1,
     1    imaginary_part_1


      do j=1,num_tf
         sigma=dble(j)*dw       !cm-1
         omega=2.d0*pi*cel*sigma !s-1

! because of i prefactor : inversion between Real/Imaginary parts of Fourier transform
! 1/2:Fourier-Laplace instead of Fourier 
         real_part = -(1.d0/2.d0)/omega*beta*tf2(j) 
         imaginary_part = (1.d0/2.d0)/omega*beta*tf(j) 

         write(10,*)sigma,real_part,imaginary_part
!         write(11,*)sigma,(real_part)**2+(imaginary_part)**2

         if(sigma.ge.2500)
     1       write(11,*)sigma,
     1       (real_part-real_part_1)**2+
     1       (imaginary_part-imaginary_part_1)**2

      enddo
      close(10)
      close(11)


! ======= 


      STOP
      END

      real*8 function filter(ind,ncor)
      implicit none
      integer ind,ncor
      
      filter =dexp(-0.5d0*100.d0*(dble(ind)/dble(ncor))*
     1     (dble(ind)/dble(ncor)))

      return
      end

      real*8 function filter2(ind,ncor)
      implicit none
      integer ind,ncor
      
      filter2 =1.0d0-dble(ind)/dble(ncor)

      return
      end


      real*8 function filter3(ind,ncor)
      implicit none
      integer ind,ncor
      real*8 pi

      pi=dacos(-1.d0)
      filter3=1.0+0.5*cos(pi*dble(ind)/dble(ncor))

      return
      end
