        program checknfsymm
!---------------------------------------------------------------
!     Purpose:  Check the summetry of the computed near field
!
!     Author: Fabio Della Sala 2009
!     e-mail: fabio.dellasala@unisalento.it
!
!----------------------------------------------------------------
       implicit none
!----------------------------------------
       character*80  fdip,filename, filename_points, filename_output
       real*8 LAM       
!------------------------
       integer OS
       integer pcx,pcy,pcz
       integer imix,imiy,imiz
       integer imax,imay,imaz
       real*8  rmax,rmay,rmaz
       real*8  rcx,rcy,rcz
       real*8  rmix,rmiy,rmiz
!------------------------------------
       integer*8 mem
!-------------------------------------
        integer wimix,wimiy,wimiz                                    
        integer wimax,wimay,wimaz
       integer ixu,iyu,izu,itot,ierr,i
       real*8   xu,yu,zu,field
       real*8, allocatable ::  tot(:,:,:)
!----------------------------------
        write(*,'(A)') 'CHECK NEAR FIELD SYMMETRY:'
!--------------read input -------------------------
        write(*,*) 'Input filenamedip, DIPY.SAV (or DIPZ.SAV)'
        read(*,*,ERR=902) fdip 
        write(*,*) 'Input  wavelength (in nanometers)' 
        read(*,*,ERR=902) LAM                 ! wave lenght
        write(*,*) 'Input filename_for_list_of_ext_points' 
        read(*,*,ERR=902) filename_points
        write(*,*) 'Input filename_output'        
        read(*,*,ERR=902) filename
!       read(*,*,ERR=902) eincflag
!       read(*,*,ERR=902) iformat  
!------------- -read DIP.SAV----------------------------
        open(UNIT=738,FILE=fdip,STATUS='OLD',FORM='UNFORMATTED')
        read(738)  OS
        read(738)  rmix,rmax,imix,imax
        read(738)  rmiy,rmay,imiy,imay
        read(738)  rmiz,rmaz,imiz,imaz
        read(738)  rcx,rcy,rcz
        read(738)  pcx,pcy,pcz
        write(*,*) 'particle center', pcx/2.d0,pcy/2.d0,pcz/2.d0
        close(738)
!-------------------------------------------------------
        write(filename_output,'(A,".DAT")')
     &  TRIM(FILENAME)
        write(*,*) 'reading',filename_output

        open(UNIT=90,FILE=filename_output,STATUS='OLD')

        wimax=-100000000
        wimay=-100000000
        wimaz=-100000000
        wimix= 100000000
        wimiy= 100000000
        wimiz =100000000
        
        itot=0                
 3332   read(90,*,END=3333)   xu,yu,zu
         ixu=xu
         iyu=yu
         izu=zu
       
         wimax=max(ixu,wimax)
         wimay=max(iyu,wimay)
         wimaz=max(izu,wimaz)
         wimix=min(ixu,wimix)
         wimiy=min(iyu,wimiy)
         wimiz=min(izu,wimiz)
         itot=itot+1
         goto 3332
 3333    continue

         write(*,*) 'points ', itot
         write(*,*) 'regions :' 
         write(*,*) 'x',wimix,wimax
         write(*,*) 'y',wimiy,wimay
         write(*,*) 'z',wimiz,wimaz
         close(90)

         mem=(wimax-wimix+1)*(wimay-wimiy+1)*(wimaz-wimiz+1)
         mem=mem*8
         write(*,*) 'Memory to be allocated',
     &               mem*1.d0/1024.d0/1024.d0,'MB'

         allocate( tot(wimix:wimax,wimiy:wimay,wimiz:wimaz), stat=ierr)
         if (ierr.ne.0) stop 'cannot allocate tot' 
         tot=0.d0

!        reread file and store in tot

         open(UNIT=90,FILE=filename_output,STATUS='OLD')
         do i=1,itot
          read(90,*)   xu,yu,zu, field
          ixu=xu
          iyu=yu
          izu=zu
!         tot  might be overwritten if the external point are on fine grid
!         tot is not fully assigend for interal field check (if not rectualar shape)
!         tot is fully assigned for planes/points not for arbiraty points
          tot(ixu,iyu,izu)=field
         enddo
         close(90)

         call  checkplanes(TOT,wimix,wimax,wimiy,wimay,wimiz,wimaz,
     &                   pcx,pcy,pcz)   
        
         deallocate(TOT)
         write(*,*) 'Bye Bye from CheckNFSymm' 
         stop
 902     write(*,*) ' ERROR in input'
         write(*,*) ' Format of the input file:'
         write(*,*) 
         write(*,*) '  DIPY.SAV (or DIPZ.SAV)'
         write(*,*) '  wavelength (in nanometers)'
         write(*,*) '  filename_for_list_of_external_points'
         write(*,*) '  filename_output'
         end program

