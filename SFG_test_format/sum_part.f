!ccccccccccccccccccccccccccccccccccccccccccccccccc
!c  
!c  19 June 2017-Evry
!c  Daria 
!cccccccccccccccccccccccccccccccccccccccccccccccccc

      program sum_part

      implicit none


      integer Nf
      integer io

      real*8 ::  Im, Re
      real*8, dimension(:), allocatable ::  ind, Retot, Imtot
      

c     Internal variables
      character(len=2)  str
      character(len=10) inputA
      character(len=8) input
      character(len=15) output


      integer k,i


ccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccc


      integer narg
      character(len=256) arg 
      
      narg=command_argument_count()

      if (narg .eq. 3) then

         call get_command_argument(1,arg)
         read(arg,'(i10)') Nf

         call get_command_argument(2,input)
         call get_command_argument(3,output)

      else
         write (*,*) 
         write (*,*) 'Usage: sum_part Nf input output '
         write (*,*) '----------------------------------------'
         write (*,*) 'Nf =  Number of files that must be sum  '
         write (*,*) ' input = root of the name'
         write (*,*) ' output = root of the name'
         write (*,*) '----------------------------------------'
         write (*,*) 
         stop
      endif

      allocate(Retot(4000))
      allocate(Imtot(4000))
      allocate(ind(4000))



      do k = Ni, Nf
c        Call the k-th file for the field
         write(str,'(i1.1)') k

         inputA = trim(trim(input)//trim(str))
 
         open(10,file=inputA,form='formatted')          
         do i=1,4000
            read(10,*,iostat=io) ind(i),  Re, Im
                !if (io/=0) then
                !       stop
                !endif
            Imtot(i)= imtot(i)+Im
            Retot(i)= Retot(i)+Re
         enddo !i
         close(10)
      enddo !k

c     Normalize
      do i=1, 4000
         Imtot(i)=Imtot(i)/Nf
         Retot(i)=Retot(i)/Nf
      enddo



      open(100,file=output,form='formatted')
      do i=1, 4000
           write(100,*) ind(i), Retot(i), Imtot(i)
      enddo
      close(100)


      deallocate(Retot)
      deallocate(Imtot)
      deallocate(ind)

      Write(77,*) 'Normal termination =)'



      end

