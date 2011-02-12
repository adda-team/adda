      PROGRAM NEARFIELD
!---------------------------------------------------------------
!     Purpose:  Compute the Near Field on 
!                an arbitrary set of points, using  ADDA dipoles
!
!     Author: Fabio Della Sala, Stefania D'Agostino, 2009
!     e-mail: fabio.dellasala@unisalento.it
!
!      S. D'Agostino, P.P. Pompa, R. Chiuri, R. Phaneuf, 
!      D. G. Britti, R. Rinaldi, R. Cingolani   and F. Della Sala
!      "Enhanced fluorescence by metal nanospheres on metal substrates"
!  	Optics Letters Vol. 34, pp.2381 (2009)
!----------------------------------------------------------------
      implicit none
 
      call computenearfield
    
      END PROGRAM

!**************************************************
      subroutine computenearfield
!**************************************************
      IMPLICIT NONE
!...........................................
      real*8 PIGRECO
      parameter(PIGRECO=3.1415926535897932384626433832795d0)

      INCLUDE 'mpif.h'
      integer  numprocs,myid,master

      external normsquared
      real*8   normsquared  
!.................................................

      real*8   xf_xmin,xf_xmax,xf_ymin,xf_ymax
      integer  NGPX,NGPY,what
      real*8   minv,maxv,minvlog,maxvlog,tmpval

      integer OS,NGP

      real*8 dir(3),pol(3)
      REAL*8 AKD,AKD2,vol,gridspa,rrr,iii
      REAL*8 DX(3)


      COMPLEX*16 CXUNIT,CXZERO

      complex*16 :: CXZEM(3),FINALEM(3)

       real*8 eexr,eexi,eeyr,eeyi,eezr,eezi
       complex*16 EE(3)


      complex*16, allocatable :: INCBEAM(:,:)
      complex*16, allocatable :: CXPOL(:,:)
      complex*16, allocatable :: INTFIELD(:,:)
      real*8,     allocatable :: POS(:,:)
       
      REAL*8 xu,yv,zw

      INTEGER ioscommvar,NGP_POINT,i
      integer perc,U
      complex*16 alpha,succe

!.............parameter frominput file..........
      character*80 fdip
      character*80 filename_points 
      character*80 FILENAME 
      REAL*8       LAM           ! in NANOMETERS  
      integer      debug         
      integer      eincflag      ! 0 or 1
      integer      calcinternal
      integer      iformat       ! 0 or 1 or 2
!...................filenames............
       character*80
     & filename_output,
     & filename_output_xf,
     & filename_output_logxf,
     & filename_output_xmgr
    

      integer ierr
      integer T,NT, FINALNT
      COMPLEX*16 G1,
     &     G2,G3,
     &     G4,G5,
     &     G6,G7,
     &     G8,G9
      complex*16 DDX,DDY,DDZ
      complex*16 VVZ,VVY,VVX    
      complex*16 som1,som2,som3
      real*8 som,field,lfie,mazdev

      integer iext  
      integer imax,imay,imaz,pcx,pcy,pcz
      integer imix,imiy,imiz
      real*8 rmax,rmay,rmaz,rcx,rcy,rcz
      real*8 rmix,rmiy,rmiz
      real*8 mindist(3)


!      xfarbe param 
      integer shift,levels,pal
      integer*8 mem,amem

  

! .............start program.............
      master=0
      shift=0
      levels=32
      pal=1

      CALL MPI_INIT(IERR)
      if (ierr.gt.0) stop 'mpierror'
      
      CALL MPI_COMM_RANK(MPI_COMM_WORLD,MYID,IERR)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NUMPROCS,IERR)     
!      write(*,*) 'MYID=',MYID,'  NUMPROCS=',NUMPROCS,' IERR=',IERR
            
          

      if (MYID.EQ.master) then
       write(*,*) 
       write(*,'(A)') 'nearfield :'
       write(*,*)     '   MPI numprocs:',NUMPROCS
       write(*,*) 

       write(*,*) 'Read Input { '
       write(*,*) 'Input filenamedip, DIPY.SAV (or DIPZ.SAV)'
       read(*,*,ERR=902) fdip
       write(*,*) 'Input  wavelength (in nanometers)'
       read(*,*,ERR=902) LAM   ! wave lenght 
       AKD=2.d0*PIGRECO/LAM                             
       write(*,*) 'Input filename_for_list_of_ext_points' 
       read(*,*,ERR=902) filename_points      
       write(*,*) 'Input filename_output'                        
       read(*,*,ERR=902) filename  
       write(*,*) 'Input  EincFlag (0=no/1=yes) '
       read(*,*,ERR=902) eincflag
       write(*,*) 'Input Format_of_the_output ( 0/1/2) '
       read(*,*,ERR=902) iformat
       write(*,*) 'Input DebugFlag ( 0/1/2/3) ' 
       read(*,*,ERR=902) debug        
       write(*,*) 'Input CalcInternalFlag (0=no/1=yes) ' 
       read(*,*,ERR=902) calcinternal
       write(filename_output,'(A,".DAT")')                                          
     &  TRIM(FILENAME)
       write(filename_output_xf,'(A,".xf")')
     &  TRIM(FILENAME)                                                 
       write(*,'(A,A40)') ' FDIP:     ',fdip
       write(*,*)          'AKD',AKD                                                             
                                   
       write(*,'(A,A40)') ' OUTXF:    ', filename_output_xf                    
       write(*,'(A,A40)') ' POINTS:   ', filename_points                                              
       write(*,*)          'DEBUG:    ', debug                                                     
       write(*,*)          'CalcInt:  ', calcinternal 
       write(*,*) 
       write(*,'(A)') ' NearField values will contains'
       if (eincflag.eq.0) 
     & write(*,*) '  Edipolar'                                    
       if (eincflag.eq.1) 
     & write(*,*) '  Etotal=Edipolar+Eincident' 
       write(*,'(A,A40)') ' NearField will be stored on file: ', 
     & filename_output      
       write(*,*)         'using format level ', iformat

!       if (eincflag.eq.0) 
!     & write(*,*) 'Computing EDIP'                                    
!       if (eincflag.eq.1) 
!     & write(*,*) 'Computing ETOT=EDIP+EINC'     
!      note the incident field must be know for all external point

       write(*,*) '} input'

       write(*,*) 
       write(*,'(A,A40)') 'Reading { ', fdip
       open(UNIT=738,FILE=fdip,FORM='UNFORMATTED',ERR=901)
       
       read(738)  OS
       read(738)  rmix,rmax,imix,imax
       read(738)  rmiy,rmay,imiy,imay
       read(738)  rmiz,rmaz,imiz,imaz
       read(738)  rcx,rcy,rcz 
       read(738)  pcx,pcy,pcz
       read(738)  mindist(1),mindist(2),mindist(3)
       write(*,*) 'NumDipoles=',OS
       write(*,'("X",2F12.6,2I6)')  rmix,rmax,imix,imax
       write(*,'("Y",2F12.6,2I6)')  rmiy,rmay,imiy,imay
       write(*,'("Z",2F12.6,2I6)')  rmiz,rmaz,imiz,imaz
       write(*,*) 'domain center', rcx,rcy,rcz 
       write(*,*) 'particle center', pcx/2.d0,pcy/2.d0,pcz/2.d0
       write(*,'("MD", 3F12.6)')  mindist(1),mindist(2),mindist(3)
       DX(1)=MINDIST(1)*1000
       DX(2)=MINDIST(2)*1000
       DX(3)=MINDIST(3)*1000
       write(*,'("DX ",3F12.6)') DX(1),DX(2),DX(3) 

      endif ! master

      call MPI_BCAST(OS, 1, MPI_INTEGER, master,
     &  MPI_COMM_WORLD, ierr)

      call MPI_BCAST(AKD, 1, MPI_DOUBLE_PRECISION, master,
     &  MPI_COMM_WORLD, ierr)

      call MPI_BCAST(DX, 3, MPI_DOUBLE_PRECISION, master,
     &  MPI_COMM_WORLD, ierr)



       CXZERO=(0.d0,0.d0)
 
       CXUNIT=(0.d0,1.d0)

       if (MYID.eq.master) open(UNIT=2145,
     & FILE=filename_output,STATUS='UNKNOWN')

       if (MYID.eq.master) then 
        mem=(OS*8)
        amem=(1+2)*3
        write(*,*) 'Memory to be allocated (all nodes)',
     &               mem*amem*1.d0/1024.d0/1024.d0,'MB'
       endif
 
       allocate ( POS(1:OS,1:3),stat=ierr)          ! dipoles pos
       if (ierr.ne.0) stop 'cannot allocate POS'    

       allocate ( CXPOL(1:OS,1:3),stat=ierr)        ! dipole
       if (ierr.ne.0) stop 'cannot allocate CXPOL'
 
       CXPOL=(0.d0,0.d0)
       POS=0.d0

       if ( MYID.eq.master) then
         read(738) (CXPOL(i,1),i=1,OS)
         read(738) (CXPOL(i,2),i=1,OS)
         read(738) (CXPOL(i,3),i=1,OS)
         read(738) (POS(i,1),i=1,OS)
         read(738) (POS(i,2),i=1,OS)
         read(738) (POS(i,3),i=1,OS)
         read(738) iext
        if (iext.ge.1.and.debug.ge.1) amem=amem+2*3
        if (iext.ge.2.and.debug.ge.2) amem=amem+2*3
         write(*,*) 'Total Memory  allocated (master)',
     &               mem*amem*1.d0/1024.d0/1024.d0,'MB'

         write(*,*) 'Extensions',iext
         if (iext.ge.1.and.debug.ge.1) then
          write(*,*) 'Reading IncBeam'
          allocate (INCBEAM(1:OS,1:3),stat=ierr)
          if (ierr.ne.0) stop 'cannot allocate INCBEAM'
          read(738)  (INCBEAM(i,1),i=1,OS)
          read(738)  (INCBEAM(i,2),i=1,OS)
          read(738)  (INCBEAM(i,3),i=1,OS)
         endif
!         if (einflag.gt.0.and.iext.eq.0) 
!     &    stop ' EincFlag=1, but IncBeam not stored in DIP.SAV'
         if (iext.ge.2.and.debug.ge.2) then
          write(*,*) 'Reading Internal Field'
          allocate (INTFIELD(1:OS,1:3),stat=ierr)
          if (ierr.ne.0) stop 'cannot allocate FIELD'
          read(738)  (INTFIELD(i,1),i=1,OS)
          read(738)  (INTFIELD(i,2),i=1,OS)
          read(738)  (INTFIELD(i,3),i=1,OS)
         endif
         if (debug.ge.3) then 
          write(*,*) 'Output Dipoles'
          do i=1,4
          write(*,'(I8,3F5.1,7D18.10)')  i,POS(i,1),POS(i,2),POS(i,3),
     &     CXPOL(i,1),CXPOL(i,2), CXPOL(i,3),    
     &     normsquared(CXPOL(i,1),CXPOL(i,2), CXPOL(i,3)) 
          enddo
          write(*,*) '....'
          do i=OS-4,OS
          write(*,'(I8,3F5.1,7D18.10)')  i,POS(i,1),POS(i,2),POS(i,3),
     &      CXPOL(i,1),CXPOL(i,2), CXPOL(i,3),
     &       normsquared(CXPOL(i,1),CXPOL(i,2), CXPOL(i,3))    
          enddo
          endif
         write(*,'(A3,A40)') ' } ',fdip  
        endif ! master


       call MPI_BCAST(CXPOL(1,1), OS, MPI_DOUBLE_COMPLEX, master,
     &  MPI_COMM_WORLD, ierr)
       if ( MYID.eq.master) write(*,*) 'bcast 1'        
       call MPI_BCAST(CXPOL(1,2), OS, MPI_DOUBLE_COMPLEX, master,
     &  MPI_COMM_WORLD, ierr)
        if ( MYID.eq.master) write(*,*) 'bcast 2'      
       call MPI_BCAST(CXPOL(1,3), OS, MPI_DOUBLE_COMPLEX, master,
     &  MPI_COMM_WORLD, ierr)
       if ( MYID.eq.master)  write(*,*) 'bcast 3'      
       call MPI_BCAST(POS(1,1), OS, MPI_DOUBLE_PRECISION, master,
     &  MPI_COMM_WORLD, ierr)
        if ( MYID.eq.master) write(*,*) 'bcast 1'      
       call MPI_BCAST(POS(1,2), OS, MPI_DOUBLE_PRECISION, master,
     &  MPI_COMM_WORLD, ierr)
       if ( MYID.eq.master) write(*,*) 'bcast 2'      
       call MPI_BCAST(POS(1,3), OS, MPI_DOUBLE_PRECISION, master,
     &  MPI_COMM_WORLD, ierr)
        if ( MYID.eq.master) write(*,*) 'bcast 3'      


         ioscommvar=-1
         if ( MYID.eq.master) then 

          open(UNIT=432,FILE='commvar',STATUS='OLD',IOSTAT=ioscommvar)
         
          if (ioscommvar.eq.0) then
           write(*,*)
           write(*,*) 'Reading Commvar { '
!          write(*,*) '---------------------'
           read(432,*) AKD2

           read(432,*) dir(1),dir(2),dir(3)
           read(432,*) pol(1),pol(2),pol(3) 
           if (debug.ge.2) then 
            read(432,*) gridspa,vol

            read(432,*) rrr,iii
            alpha=rrr+iii*CXUNIT

            read(432,*) rrr,iii
            succe=(rrr+iii*CXUNIT)/vol

            write(*,*) 'AKD2',AKD2
            if (dabs(AKD2/1000-AKD).gt.1.d-7)  then                
             write(*,*) 'input lambda differs from adda calculations'
             write(*,*) 'adda akd',akd2/1000            
             write(*,*) 'input akd',akd                 
             stop  
            endif 
           endif

           write(*,*) 'DIR', dir(1),dir(2),dir(3)
           write(*,*)  'POL',pol(1),pol(2),pol(3) 
! P=a Eexc = V succe Etot

           close(432)

           write(*,*) '} commvar'    

           if (debug.ge.1.and.iext.ge.1) then
!               use AKD in nm because DX is in nm
            call testbeam(rcx,rcy,rcz,debug,
     &                   POS,OS,INCBEAM,dir,pol,AKD,DX)  
           else 
            if (debug.ge.1.and.iext.eq.0) then
             write(*,*) 'WARNING!: cannot check incident field'
             write(*,*) '          run adda with -store_beam'
             write(*,*) '          and rerun gendip'
            endif
           endif 

           if (debug.ge.2.and.iext.ge.2) then   
            call testinternalfield(INTFIELD,CXPOL,INCBEAM,OS,
     &       AKD2,  ! uses the adda in micromenters
     &       alpha,succe,vol)

             gridspa=gridspa*1000   ! ora in nm
             vol=vol*1.d9           ! ora in nm^3
        
             DX(1)=gridspa  ! in nm
             DX(2)=gridspa  ! in nm
             DX(3)=gridspa  ! in nm
             AKD2=AKD2/1000  ! 2pi * (lam *1000) 
             write(*,*) 'DX',DX(1),DX(2),DX(3)
             write(*,*) 'AKD2',AKD2 
           endif ! debuadda
 
          else ! *commbvar*
           write(*,*) 'cannot find commvar' 
           if (eincflag.eq.1) then
            eincflag=0
            write(*,*) 
     & 'WARNING! : Cannot compute EincField on external points !!'
           endif
          endif        

          if (debug.ge.3) 
     &  call  checksymm(imix,imax,imiy,imay,imiz,imaz,
     &                       pcx,pcy,pcz,OS,CXPOL,POS)   

         endif  ! master

         DO T=1,OS
            CXPOL(T,1)= CXPOL(T,1)*1.d9
            CXPOL(T,2)= CXPOL(T,2)*1.d9
            CXPOL(T,3)= CXPOL(T,3)*1.d9
          enddo

       
       if ( MYID.eq.master) then 
        write(*,*) !'AKD',AKD 
        open(UNIT=2101,   FILE=filename_points,STATUS='OLD')
        read(2101,*) what
        write(*,*) 'point input what',what 
        if (what.eq.2) then 
         read(2101,*) NGPX,xf_xmin,xf_xmax
         read(2101,*) NGPY,xf_ymin,xf_ymax
         read(2101,*) NGP
         write(*,*) 'PLOTTING  A PLANE; NPOINTS=',NGP 
        else if (what.eq.1) then
         read(2101,*) xf_xmin,xf_xmax, NGP
         write(*,*) 'PLOTTING  A LINE; NPOINTS=',NGP 
        else if (what.eq.0) then
         read(2101,*) NGP
         write(*,*) 'PLOTTING ANY; NPOINTS=',NGP 
        endif
        write(*,*) 'Reading File ',TRIM(filename_points)
       endif

       call MPI_BCAST(NGP, 1, MPI_INTEGER, master,
     & MPI_COMM_WORLD, ierr)
               
       SOM=0.d0
       field=0.d0
       NGP_POINT=0

      
       if (MYID.eq.master) then 
         if (what.eq.2) then
        write(filename_output_xf,'(A,".xf")')                                        
     & TRIM(FILENAME)!,numprocs           
        write(filename_output_logxf,'(A,"log.xf")')                                        
     & TRIM(FILENAME)!,numprocs 
    
!           write(*,*) 'prima',MYID 
        call prepxf(2445,filename_output_xf,minv,maxv,shift)
        call prepxf(2446,filename_output_logxf,
     & minvlog,maxvlog,shift)   
!           write(*,*) 'dopo',MYID 
          write(2445,*) NGPY,NGPX
          write(2446,*) NGPY,NGPX
        else if (what.eq.1) then
            write(filename_output_xmgr,'(A,".xmgr")')                                        
     & TRIM(FILENAME)!,numprocs       
         open(UNIT=2445,FILE=filename_output_xmgr,STATUS='UNKNOWN')
        endif

        endif
 
      if (MYID.eq.master) then
       if (eincflag.eq.1) then
        write(*,*) 'Computing Near-Field + E_Inc ....................'
       else
        write(*,*) 'Computing Near-Field ............................'
       endif

      write(*,'(A8,"/",A8,"/",A8,3A10,A16)') 
     &  '  Points',
     &  ' Outside',
     &  '   Total',
     &  '    X   ',
     &  '    Y   ',
     &  '    Z   ',
     &  '      Field     '
      endif

      
      do U=1,NGP ! <<<<<<<<<<<<<<< LOOP EXTERNAL <<<<<<<<<<<<<<<<<<<<<<
! main loop over external points;
! note: skip point INSIDE the target

       FINALEM=CXZERO
       FINALNT=0

!      master read external point
       if ( MYID.eq.master) then 
        READ(2101,*) xu,yv,zw 
       endif

         
!      send external point         
       call MPI_BCAST(xu, 1, MPI_DOUBLE_PRECISION, master,
     &  MPI_COMM_WORLD, ierr)
       call MPI_BCAST(yv, 1, MPI_DOUBLE_PRECISION, master,
     &  MPI_COMM_WORLD, ierr)
       call MPI_BCAST(zw, 1, MPI_DOUBLE_PRECISION, master,
     &  MPI_COMM_WORLD, ierr)


!       task local variables  
        CXZEM(1)=CXZERO
        CXZEM(2)=CXZERO
        CXZEM(3)=CXZERO
        NT=0

                 DO T=1+MYID,OS,NUMPROCS ! <== loop over dipoles
!                 write(*,*) T,OS

                  call  computegreen(xu,yv,zw,      ! external point
     &                   POS(T,1),POS(T,2),POS(T,3), ! dipole point
     &                   AKD,DX,
     &                   g1,g2,g3,g4,g5,g6,g7,g8,g9,NT) 

!               multiply by dipoles
                         som1=-(G1*CXPOL(T,1)+
     &                          G2*CXPOL(T,2)+
     &                          G3*CXPOL(T,3))

                         som2=-(G4*CXPOL(T,1)+
     &                          G5*CXPOL(T,2)+
     &                          G6*CXPOL(T,3))

                         som3=-(G7*CXPOL(T,1)+
     &                          G8*CXPOL(T,2)+
     &                          G9*CXPOL(T,3))
!                finally sum
                         CXZEM(1)=CXZEM(1)+som1
                         CXZEM(2)=CXZEM(2)+som2
                         CXZEM(3)=CXZEM(3)+som3
                 ENDDO  ! <============== loop over dipoles 

!         sum over all tasks
!          final field in FINALEM(1:3)
!          number of dipoles in FINALNT    
         call MPI_REDUCE(CXZEM(1), FINALEM(1), 1, MPI_DOUBLE_COMPLEX,
     &                MPI_SUM, master, MPI_COMM_WORLD, ierr)

         call MPI_REDUCE(CXZEM(2), FINALEM(2), 1, MPI_DOUBLE_COMPLEX,
     &                MPI_SUM, master, MPI_COMM_WORLD, ierr)

         call MPI_REDUCE(CXZEM(3), FINALEM(3), 1, MPI_DOUBLE_COMPLEX,
     &                MPI_SUM, master, MPI_COMM_WORLD, ierr)

         call MPI_REDUCE(NT     , FINALNT,     1, MPI_INTEGER,
     &                MPI_SUM, master, MPI_COMM_WORLD, ierr)


         IF (MYID.eq.master) then
           if (FINALNT.eq.OS) then    
            NGP_POINT=NGP_POINT+1   
           else if (calcinternal.eq.0) then 
!           the external point is close to dipoles => diseregard
            FINALEM(1)=CXZERO
            FINALEM(2)=CXZERO
            FINALEM(3)=CXZERO              
           endif


           if (eincflag.eq.1) then
              call calcbeam(xu,yv,zw,DX,
     &           rcx,rcy,rcz,dir,pol,AKD,EE)
              FINALEM(1)=FINALEM(1)+EE(1)
              FINALEM(2)=FINALEM(2)+EE(2)
              FINALEM(3)=FINALEM(3)+EE(3)
           endif ! eincflag 

 
            field    =     cdabs(FINALEM(1))**2+
     &                     cdabs(FINALEM(2))**2+
     &                     cdabs(FINALEM(3))**2

!           write field to file output
            if (iformat.eq.0) then 
            write(2145,'(3F10.3,1X,E19.10)') 
     &               xu,yv,zw,
     &               field 
            else if (iformat.eq.1) then
             write(2145,'(3F10.3,1X,1E19.10,1X,6E16.8)') 
     &               xu,yv,zw,
     &               field,FINALEM(1),
     &               FINALEM(2),FINALEM(3)        
            else if (iformat.eq.2) then
            write(2145,'(3F10.3,1X,2E19.10,1X,6E19.10)') 
     &               xu,yv,zw,
     &               field,log10(field),FINALEM(1),
     &               FINALEM(2),FINALEM(3) 
            endif

            som=som+field

            if (what.eq.2) then 
             tmpval= field 
             write(2445,'(E16.8)') tmpval
             maxv=max( tmpval,maxv)
             minv=min( tmpval,minv)
             tmpval= log10(field) 
             write(2446,'(E16.8)')  tmpval
             maxvlog=max( tmpval,maxvlog)
             minvlog=min( tmpval,minvlog)
            else if (what.eq.1) then
             tmpval= field 
             write(2445,'(2E16.8)') 
     &        xf_xmin+(U-1)*1.d0/(NGP-1)*(xf_xmax-xf_xmin),
     &        tmpval
            endif

            
         ENDIF ! MYID


         if (MYID.EQ.master) then
                    perc=int(U/NGP)*100
                   if (mod(U,NGP/20).eq.0) then 
                    write(*,'(I8,"/",I8,"/",I8,3F10.3,E16.8)') 
     & U,NGP_POINT,NGP,xu,yv,zw,field !,som
                   endif
                  endif

         enddo    !<<<<<<<<<<<< LOOP EXTERNAL <<<<<<<<<<<<<<

         IF (MYID.eq.master) then 
          write(*,*) 'POINT OUTSIDE TARGET',NGP_POINT
          write(*,*) 'AVERAGE FIELD',som/NGP

          close(2101) ! file of external points
          close(2145) ! file for near field

          if (what.eq.2) then
           call postxf(2445,filename_output_xf,minv,maxv,
     &        shift,pal,levels,
     &      xf_ymin,xf_ymax,xf_xmin,xf_xmax)
           call postxf(2446,filename_output_logxf,minvlog,maxvlog,
     &        shift,pal,levels,
     &      xf_ymin,xf_ymax,xf_xmin,xf_xmax)
          else if (what.eq.1) then
           close(2445)
          endif

!                     
        if (NGP_POINT.eq.0.and.NGP.eq.OS
     &      .and.debug.ge.2.and.iext.ge.2
     &      .and.ioscommvar.eq.0) then
          write(*,*)
          write(*,*) 'Check Near Field with Internal { '

          mazdev=0.d0
!          open(UNIT=2310,FILE='DIP_CHECK',STATUS='UNKNOWN')
          open(UNIT=2145,
     &       FILE=filename_output,STATUS='UNKNOWN') 
          write(*,'(A,A40)') ' reload ',filename_output

         do T=1,NGP 
         read(2145,*)xu,yv,zw,
     &               field,lfie,
     &    eexr,eexi,eeyr,eeyi,eezr,eezi  ! OUT

          DDX=1.d0*(eexr + CXUNIT*eexi) !DIP near
          DDY=1.d0*(eeyr + CXUNIT*eeyi) !DIP near
          DDZ=1.d0*(eezr + CXUNIT*eezi) !DIP near

          if (eincflag.eq.1) then 
! a^-1 P  =  EDIP +EINC
           VVX=   CXPOL(T,1)/alpha/1.d9
           VVY=   CXPOL(T,2)/alpha/1.d9
           VVZ=   CXPOL(T,3)/alpha/1.d9
          else
! -EINC + a^-1 P  =  EDIP 
!          QQX= 1.d0*(INTFIELD(T,1))
!          QQY= 1.d0*(INTFIELD(T,2))
!          QQZ= 1.d0*(INTFIELD(T,3))
           VVX=-1.d0*INCBEAM(T,1) +   CXPOL(T,1)/alpha/1.d9 
           VVY=-1.d0*INCBEAM(T,2) +   CXPOL(T,2)/alpha/1.d9
           VVZ=-1.d0*INCBEAM(T,3) +   CXPOL(T,3)/alpha/1.d9
          endif
          mazdev=max(mazdev,
     &      dabs(cdabs(DDX/VVX)+cdabs(DDY/VVY)+cdabs(DDZ/VVZ)-3.d0))
!         write(2310,'(12D16.8)') DDX/VVX,DDY/VVY,DDZ/VVZ
          enddo
          write(*,*) ' ERROR INF : ',mazdev
          close(2145) ! nearfield
!         close(2310) ! DIP_CHECK
          write(*,*) '} Check Near Field with Internal'
        endif

        endif ! MASTER
!##############################################

 777     continue

!        dellocation
         deallocate(POS)
         deallocate(CXPOL)         
         if (MYID.eq.master) then 
          if (iext.ge.1.and.debug.ge.1) deallocate(INCBEAM)
          if (iext.ge.2.and.debug.ge.2) deallocate(INTFIELD)
         endif

         
         if (MYID.eq.master) then
         write(*,*) 
         write(*,*) 'Bye Bye from nearfield' 
         endif

         call MPI_FINALIZE(ierr)
         return

 901     write(*,*) ' ERROR!: Cannot Open file ',fdip
         stop
 902     write(*,*) ' ERROR in input'
         write(*,*) ' Format of the input file:'
         write(*,*) 
         write(*,*) '  DIPY.SAV (or DIPZ.SAV)'
         write(*,*) '  wavelength (in nanometers)'
         write(*,*) '  filename_for_list_of_external_points'
         write(*,*) '  filename_output'
         write(*,*) '  EincFlag (0=no , 1=yes) '
         write(*,*) '  Format_of_the_output ( 0 or 1 or 2) '
         write(*,*) '  DebugFlag ( 0 or 1 or 2 or 3) '
         write(*,*) '  CalcInternalFlag (0=no 1=yes) '
         stop
         end







!****************************************************************
!****************************************************************

       subroutine computegreen(xu,yv,zw,xi,yj,zk,AKD,DX,
     & g1,g2,g3,g4,g5,g6,g7,g8,g9,NT)
       implicit none
!.............param..............
       REAL*8  xu,yv,zw 
       real*8  xi,yj,zk
       real*8  DX(3)
       REAL*8  X1,X2,X3
       REAL*8 AKD
     
       integer NT ! inout  
! ...............output............
       COMPLEX*16 G1,
     &     G2,G3,
     &     G4,G5,
     &     G6,G7,
     &     G8,G9
! ...............local..........
       complex*16 CXUNIT, CXZERO   
       real*8
     &     R,R2,
     &     R3,a1,
     &     a2,a3,
     &     a4,a5,
     &     a6,a7,
     &     a8,a9
       complex*16 tmp0,tmp1,tmp12,tmp13


                         CXUNIT=(0.d0,1.d0)
                         CXZERO=(0.d0,0.d0)
                         X1=((xu-xi)*DX(1))
                         X2=((yv-yj)*DX(2))
                         X3=((zw-zk)*DX(3))
     
                         R2=X1**2+X2**2+X3**2
                         R=dsqrt(R2)
                         R3=R*R2

                         IF(R/DX(1).GT.0.866d0)THEN
                           NT=NT+1
!              compute matrix  3x3 Rij*Rij
                           a1=X1**2
                           a2=X1*X2
                           a3=X1*X3
                           a5=X2**2
                           a6=X2*X3
                           a9=X3**2
                           a4=a2
                           a7=a3
                           a8=a6
!              compute Green matrix  3x3 Gij
                           tmp0=CXUNIT*AKD*R
                           tmp1=CDEXP(tmp0)/R3
                           tmp12=tmp1*(AKD**2+3*(tmp0-1)/R2)
                           tmp13=tmp1*(-AKD**2*R2-tmp0+1)
                           G1=(a1*tmp12+tmp13)
                           G2=(a2*tmp12)
                           G3=(a3*tmp12)
                           G5=(a5*tmp12+tmp13)
                           G6=(a6*tmp12)
                           G9=(a9*tmp12+tmp13)
                         ELSE ! external point too near to dipole
                           G1=CXZERO
                           G2=CXZERO
                           G3=CXZERO
                           G5=CXZERO
                           G6=CXZERO
                           G9=CXZERO
                         ENDIF ! R<0.5
                         G4=G2
                         G7=G3
                         G8=G6
            end subroutine

