      Program Recenter
      Implicit none

      Integer::narg
      Character(len=60), allocatable, dimension(:):: arg
      Double precision, allocatable, dimension(:,:):: X, Y, Z
      Character(len=4), allocatable, dimension(:,:):: Sb
      Character(len=4) :: atlo="null"
      Character(10) :: txt
      Double precision:: A, B, C
      Integer:: Nat, nmovie
      Integer:: i, j, k
      Double precision, allocatable, dimension(:,:):: CM,axis
      Double precision:: Mass_H,Mass_O,Mass_Si,Mass_Cl,Mass_K,Mass_Al
      Double precision:: CM_SUM_MiRi,CM_SUM_Mi

      !COUNTING THE NUMBER OR ARGUMENTS
      narg=command_argument_count()

      !print*, "Number of arguments =", narg

      If ( narg .lt. 5 ) then
       print*,'-----------------------------------------------------------------&
     &--------------------------------------------------------------------'
       print*,'This program re-build a trajectory-file centering the box in the &
     &origin'
       print*,'  ---------------------------------------------------------------&
     &------------'
       print*,'   Usage: Recenter input_name output_name A B C'
       print*,'  ---------------------------------------------------------------&
     &------------'
       print*,'         where,'
       print*,'                - input_name is the name for the file containing&
     & the XYZ positions.'
       print*,'                - output_name is the name for the file containing&
     & the rebuild trajectory, in XYZ format.'
       print*,'                - A, B and C, parameters of the box.'
       print*,'-----------------------------------------------------------------&
     &---------------------------------------------------------------------'
       stop
      Endif

      Allocate(arg(narg))

      !READING ARGUMENTS
      Do i=1,narg
       call get_command_argument(i,arg(i))
      Enddo

      !READING PARAMETERS OF THE BOX
      Read(arg(3),*) A
      Read(arg(4),*) B
      Read(arg(5),*) C


nmovie=0
      Open(Unit=9,Action="Read",File=arg(1))
      Do
       Read(9,*,end=100) Nat
       Read(9,*,end=100)
       Do i=1,Nat
        Read(9,*,end=100)
       Enddo
       nmovie=nmovie+1
      Enddo
  100 Continue

      Close(9)

      allocate(Sb(nmovie,Nat))
      allocate(X(nmovie,Nat))
      allocate(Y(nmovie,Nat))
      allocate(Z(nmovie,Nat))

      Open(Unit=10,Action="Read",File=arg(1))
      Do i=1,nmovie
      Read(10,*,end=101)
      Read(10,*,end=101)
       Do j=1,Nat
        Read(10,*,end=101) Sb(i,j), X(i,j), Y(i,j), Z(i,j)
        if (i == 1) then
           if ((trim(Sb(i,j)) == "N") .or. (trim(Sb(i,j)) == "B" .or. trim(Sb(i,j)) == "C")) then
           atlo="good"
           endif
       endif
       Enddo
      Enddo
  101 Continue
      Close(10)

      !COMPUTING THE CM FOR EACH STEP
      allocate(CM(nmovie,3))
      allocate(axis(nmovie,Nat))
      !MASS OF EACH TYPE OF ATOM
      !SOURCE CHEMICALELEMENT.COM
      Mass_H=1.00794
      Mass_O=15.9994
      Mass_Si=28.0855
      Mass_Cl=35.4527 
      Mass_K=39.0983
      Mass_Al=26.981539

      ! CALCULATE CM ON EACH AXIS      
      CM_SUM_MiRi=0
      CM_SUM_Mi=0
      
      !INITIALISATION OF AXIS ARRAY
      do i=1,nmovie
       do j=1,Nat
          axis(i,j)=0
       enddo
      enddo
      
      
      do k=1,3
        if (k==1) then
            do i=1,nmovie
               do j=1,Nat
                  axis(i,j)=X(i,j)
               enddo
             enddo           
        elseif (k==2) then
            do i=1,nmovie
               do j=1,Nat
                  axis(i,j)=Y(i,j)
               enddo
             enddo
        elseif (k==3) then
             do i=1,nmovie
               do j=1,Nat
                  axis(i,j)=Z(i,j)
               enddo
             enddo 
        endif
    
        do i=1,nmovie
          !RESET AT EACH STEP
          CM_SUM_MiRi=0
          CM_SUM_Mi=0
   
           do j=1,Nat
               !CASES OF ATOMS TYPE
               if(atlo=="null") then
               if (Sb(i,j)=='K') then
                 CM_SUM_MiRi = CM_SUM_MiRi + Mass_K * axis(i,j)
                 CM_SUM_Mi = CM_SUM_Mi + Mass_K
               endif

               if (Sb(i,j)=='Cl') then
                 CM_SUM_MiRi = CM_SUM_MiRi + Mass_Cl * axis(i,j)
                 CM_SUM_Mi = CM_SUM_Mi + Mass_Cl
               endif
       
               if (Sb(i,j)=='H') then
                 CM_SUM_MiRi = CM_SUM_MiRi + Mass_H * axis(i,j)
                 CM_SUM_Mi = CM_SUM_Mi + Mass_H
               endif

                if (Sb(i,j)=='O') then
                 CM_SUM_MiRi = CM_SUM_MiRi + Mass_O * axis(i,j)
                 CM_SUM_Mi = CM_SUM_Mi + Mass_O
               endif
               endif
              
               if (Sb(i,j)=='Al') then
                 CM_SUM_MiRi = CM_SUM_MiRi + Mass_Al * axis(i,j)
                 CM_SUM_Mi = CM_SUM_Mi + Mass_Al
               endif

               
               if (Sb(i,j)=='Si') then
                 CM_SUM_MiRi = CM_SUM_MiRi + Mass_Si * axis(i,j)
                 CM_SUM_Mi = CM_SUM_Mi + Mass_Si
               endif   
            enddo
            CM(i,k) = CM_SUM_MiRi / CM_SUM_Mi
        enddo 

      enddo
      
      !CENTERING INTO THE BOX WITH CM
      Do i=1,nmovie
           Do j=1,Nat
                X(i,j)=X(i,j)-CM(i,1)
                Y(i,j)=Y(i,j)-CM(i,2)
                Z(i,j)=Z(i,j)-CM(i,3)
           enddo
      enddo

      !CENTERING INTO THE BOX
      Do i=1,nmovie
       Do j=1,Nat
        If (X(i,j).lt.(-A/2.0d0)) then
         X(i,j)=X(i,j)+A
        Endif
        If (X(i,j).gt.(A/2.0d0)) then
         X(i,j)=X(i,j)-A
        Endif
        If (Y(i,j).lt.(-B/2.0d0)) then
         Y(i,j)=Y(i,j)+B
        Endif
        If (Y(i,j).gt.(B/2.0d0)) then
         Y(i,j)=Y(i,j)-B
        Endif
        If (Z(i,j).lt.(-C/2.0d0)) then
         Z(i,j)=Z(i,j)+C
        Endif
        If (Z(i,j).gt.(C/2.0d0)) then
         Z(i,j)=Z(i,j)-C
        Endif
       Enddo
      Enddo

      Open(Unit=15,Action="Write",File=arg(2))
      Do i=1,nmovie
       Write(15,*) Nat
       Write(15,*) 
       Do j=1,Nat
        Write(15,FMT="(A, F14.6, F14.6, F14.6)") Sb(i,j), X(i,j), Y(i,j), Z(i,j)
       Enddo
      Enddo
      Close(15)
      
      End
