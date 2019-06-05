         PROGRAM GENDIP
!-------------------------------------------------------------------------
!     Purpose:  Creates the DIPY.SAV and/or DIPZ.SAV binary file 
!               it reads the files 'target' and 
!                'DipPol-Y' and/or 'DipPol-X
!               generated by adda -save_geom  -store_dip_pol
!               rename <type>.geom to target
!               if 
!     Author: Fabio Della Sala, 2009
!     e-mail: fabio.dellasala@unisalento.it
!    
!------------------------------
         implicit none 
         integer iouty,ioutx
         write(*,*) 
         write(*,'(A)') 'gendip :'
         write(*,*)      
         call createdip('DipPol-Y','DIPY.SAV',
     &                  iouty,'IncBeam-Y','IntField-Y')
         write(*,*) 
         call createdip('DipPol-X','DIPX.SAV',
     &                  ioutx,'IncBeam-X','IntField-X')
         write(*,*)   
         if (iouty.ne.0.and.ioutx.ne.0) then
          write(*,*) 'ERROR ! : Cannot found DipPol-Y or DipPol-X'
         else
           write(*,*) 'Bye Bye from gendip' 
         endif
         end program
  
!------------------------------------------------
         subroutine createdip(fdipin,fdipout,iout,fbeam,ffie)
         implicit none
!..................input................
         character*(*) fdipin,fdipout,ffie,fbeam
         integer iout 
!..................local................
         real*8  xin,yjn,zkn,qua,
     &      CXPL1R,CXPL1I,CXPL2R,
     &      CXPL2I,CXPL3R,CXPL3I
         real*8  fxin,fyjn,fzkn,fqua,
     &      fCXPL1R,fCXPL1I,fCXPL2R,
     &      fCXPL2I,fCXPL3R,fCXPL3I
         real*8  bxin,byjn,bzkn,bqua,
     &      bCXPL1R,bCXPL1I,bCXPL2R,
     &      bCXPL2I,bCXPL3R,bCXPL3I

         integer xi,yj,zk
         integer OS
         complex*16, allocatable :: CXPOL(:,:)
         complex*16, allocatable :: INCBEAM(:,:)
         complex*16, allocatable :: FIELD(:,:)
         real*8,    allocatable :: POS(:,:)
         real*8,    allocatable :: RPOS(:,:)
         complex*16 CXUNIT
         integer imax,imay,imaz
         integer imix,imiy,imiz
         real*8 rmax,rmay,rmaz
         real*8 rmix,rmiy,rmiz
         integer MAXOS,j,i
         real*8 mindist(3),dist(3)
         integer*8 sumx,sumy,sumz
         real*8 rcx,rcy,rcz ! center of the computational domain (in 'target' reference frame)
         integer pcx,pcy,pcz ! 2*symmetry center of the particle (in 'target' reference frame)
         real*8 diffx,diffy,diffz,DX(3),bxx,byy,bzz,mad
         integer iosb,iosf,iext,ierr
         integer*8 mem,amem


         write(*,*) 'Reading ',fdipin
         open(UNIT=3113,FILE=fdipin,STATUS='OLD',ERR=901)
         iout=0
         read(3113,*) ! x y z |P|^2 Px.r Px.i Py.r Py.i Pz.r Pz.i

         write(*,*) 'Reading target ...'
!---------------------get number of dipoles-------------------------
         OS=0
         open(UNIT=3114,FILE='target',STATUS='OLD',ERR=900)
         read(3114,*) 
         read(3114,*)      
         read(3114,*)
 399     read(3114,*,ERR=400,END=200) xi,yj,zk
         goto 401
 400          write(*,*) 'NMAT>1'
         goto 399
 401          OS=OS+1
         goto 399
 200     write(*,*) 'Readed Dipoles',OS
         close(3114)

!-------------------------------------------------------------------
         write(*,*) 
         open(UNIT=3114,FILE='target',STATUS='OLD',ERR=900)
         read(3114,*) 
         read(3114,*)
         read(3114,*)
         read(3114,*)
         
         iext=0
         open(UNIT=3115,FILE=fbeam,STATUS='OLD',iostat=iosb)
         if (iosb.eq.0) then
          iext=iext+1 
          write(*,*) 'Reading ',fbeam 
          read(3115,*) ! x y z |P|^2 Px.r Px.i Py.r Py.i Pz.r Pz.i 
         endif    
        

         open(UNIT=3116,FILE=ffie,STATUS='OLD',iostat=iosf)
         if (iosf.eq.0) then
         if (iosb.ne.0) then
          iosf=-1
!              for field you need also beam
         else 
          iext=iext+1
          write(*,*) 'Reading ',ffie
           read(3116,*) ! x y z |P|^2 Px.r Px.i Py.r Py.i Pz.r Pz.i
         endif
         endif

         MAXOS=OS
         write(*,*) 'Maximum number of Dipole',MAXOS
         mem=(MAXOS*8)
         amem=(1+1+2)*3
         if (iosf.eq.0) amem=amem+2*3
         if (iosb.eq.0) amem=amem+2*3
         mem=mem*amem

         write(*,*) 'Memory to be allocated',
     &               mem*1.d0/1024.d0/1024.d0,'MB'

         allocate ( POS  (1:MAXOS,1:3),stat=ierr)
         if (ierr.ne.0) stop 'cannot allocate POS' 
         allocate ( CXPOL(1:MAXOS,1:3),stat=ierr)
         if (ierr.ne.0) stop 'cannot allocate CXPOL'  
         allocate ( RPOS (1:MAXOS,1:3),stat=ierr)
         if (ierr.ne.0) stop 'cannot allocate RPOS'
         if (iosf.eq.0) then
           allocate ( FIELD  (1:MAXOS,1:3),stat=ierr)
           if (ierr.ne.0) stop 'cannot allocate FIELD'
         endif
         if (iosb.eq.0) then
           allocate ( INCBEAM (1:MAXOS,1:3),stat=ierr)
           if (ierr.ne.0) stop 'cannot allocate INCBEAM'
         endif

         POS=0.d0
         CXPOL=(0.d0,0.d0)
         RPOS=0.d0
         CXUNIT=(0.d0,1.d0)
         if (iosf.eq.0) FIELD=(0.d0,0.d0)
         if (iosb.eq.0) INCBEAM=(0.d0,0.d0)


         rmax=-1d9
         rmay=-1d9
         rmaz=-1d9
         rmix= 1d9
         rmiy= 1d9
         rmiz= 1d9

         sumx=0
         sumy=0
         sumz=0
         imax=-100000000
         imay=-100000000
         imaz=-100000000
         imix= 100000000
         imiy= 100000000
         imiz =100000000
         mindist(1:3)=1d9

         do i=1,OS !-------------OS loop------------------
          
          
            read(3113,*,END=920) xin,yjn,zkn,qua,  ! read DipPol-Y
     &      CXPL1R,CXPL1I,
     &      CXPL2R,CXPL2I,
     &      CXPL3R,CXPL3I

            read(3114,*) xi,yj,zk

            if (iosb.eq.0) then 
            read(3115,*,END=920) bxin,byjn,bzkn,bqua,  
     &      bCXPL1R,bCXPL1I,
     &      bCXPL2R,bCXPL2I,
     &      bCXPL3R,bCXPL3I     
            endif

            if (iosf.eq.0) then 
            read(3116,*,END=920) fxin,fyjn,fzkn,fqua,  
     &      fCXPL1R,fCXPL1I,
     &      fCXPL2R,fCXPL2I,
     &      fCXPL3R,fCXPL3I    
            endif
      
         sumx=sumx+xi
         sumy=sumy+yj
         sumz=sumz+zk
        
         if (mod(i,100000).eq.0) write(*,*) 'Readed',i

         CXPOL(i,1)=  CXPL1R+CXUNIT*CXPL1I
         CXPOL(i,2)=  CXPL2R+CXUNIT*CXPL2I
         CXPOL(i,3)=  CXPL3R+CXUNIT*CXPL3I

         if (iosf.eq.0) then
          FIELD(i,1)=  fCXPL1R+CXUNIT*fCXPL1I
          FIELD(i,2)=  fCXPL2R+CXUNIT*fCXPL2I
          FIELD(i,3)=  fCXPL3R+CXUNIT*fCXPL3I
         endif

         if (iosb.eq.0) then
          INCBEAM(i,1)=  bCXPL1R+CXUNIT*bCXPL1I
          INCBEAM(i,2)=  bCXPL2R+CXUNIT*bCXPL2I
          INCBEAM(i,3)=  bCXPL3R+CXUNIT*bCXPL3I
         endif

         POS(i,1)= xi
         POS(i,2)= yj
         POS(i,3)= zk

         RPOS(i,1)= xin
         RPOS(i,2)= yjn
         RPOS(i,3)= zkn

        
         if (i.gt.1) then
          do j=1,3
          if (j.eq.1) dist(j)= abs(xin-   RPOS(i-1,1))
          if (j.eq.2) dist(j)= abs(yjn-   RPOS(i-1,2))
          if (j.eq.3) dist(j)= abs(zkn-   RPOS(i-1,3))
          if (dist(j).gt.1e-9) mindist(j)=min(mindist(j),dist(j))
          enddo
         endif

         rmax=max(xin,rmax)
         rmay=max(yjn,rmay)
         rmaz=max(zkn,rmaz)
         rmix=min(xin,rmix)
         rmiy=min(yjn,rmiy)
         rmiz=min(zkn,rmiz)
      
     
         imax=max(xi,imax)
         imay=max(yj,imay)
         imaz=max(zk,imaz)       
         imix=min(xi,imix)
         imiy=min(yj,imiy)
         imiz=min(zk,imiz)
!          if (OS.ge.MAXOS) then
!           write(*,*) 'ERROR: too many dipoles !!'
!           goto 200
!          endif
        
         enddo !-----------------  OS loop -----------------------
         
         close(3114)
         close(3113)
         if (iosb.eq.0) close(3115)
         if (iosf.eq.0) close(3116)

        
         write(*,'(A,2F15.8,2I6)') ' Dimension X',rmix,rmax,imix,imax
         write(*,'(A,2F15.8,2I6)') ' Dimension Y',rmiy,rmay,imiy,imay
         write(*,'(A,2F15.8,2I6)') ' Dimension Z',rmiz,rmaz,imiz,imaz
!          baricenter
!         rcx=1.d0*sumx/OS                                                                                  
!         rcy=1.d0*sumy/OS                                                                                                
!         rcz=1.d0*sumz/OS  
!         icx=sumx/OS 
!         icy=sumy/OS
!         icz=sumz/OS
!        These are exact coordinates of the origing of the particle reference frame (used in 
!        DipPol) relative to the 'target' reference frame for any possible situation. I.e. it
!        is more general than the formula used before.
         rcx=(rmax*imix-rmix*imax)/(rmax-rmix)
         rcy=(rmay*imiy-rmiy*imay)/(rmay-rmiy)
         rcz=(rmaz*imiz-rmiz*imaz)/(rmaz-rmiz)
!        These are two times exact coordinates of the particle symmetry center (assuming it exists), 
!        which is used in symmetry-testing functions. We use double the value to keep it integer and 
!        avoid potential problems during rounding.
         pcx= imax+imix
         pcy= imay+imiy
         pcz= imaz+imiz

         write(*,'(A,3F10.5)') 'domain center xyz',rcx,rcy,rcz
         write(*,*) 'particle center xyz',pcx/2.d0,pcy/2.d0,pcz/2.d0

!         write(*,*) ((imax-imix)+2)/2+0.5d0
!         write(*,*) ((imay-imiy)+2)/2+0.5d0
!         write(*,*) ((imaz-imiz)+2)/2+0.5d0

         write(*,*) 'DX:',(rmax-rmix)/(imax-imix),mindist(1)
         write(*,*) 'DY:',(rmay-rmiy)/(imay-imiy),mindist(2)
         write(*,*) 'DZ:',(rmaz-rmiz)/(imaz-imiz),mindist(3)

         write(*,*) 
         write(*,*) 'Creating ', fdipout 
         open(UNIT=738,FILE=fdipout,FORM='UNFORMATTED')
         write(738)  OS
         write(738)  rmix,rmax,imix,imax
         write(738)  rmiy,rmay,imiy,imay
         write(738)  rmiz,rmaz,imiz,imaz
         write(738)  rcx,rcy,rcz
         write(738)  pcx,pcy,pcz
         write(738)  mindist(1),mindist(2),mindist(3)
         write(738) (CXPOL(i,1),i=1,OS)
         write(738) (CXPOL(i,2),i=1,OS)
         write(738) (CXPOL(i,3),i=1,OS)

         write(738) (POS(i,1),i=1,OS)
         write(738) (POS(i,2),i=1,OS)
         write(738) (POS(i,3),i=1,OS)
         write(738) iext
         write(*,*) 'Extensions',iext 
         if (iext.ge.1) then 
          write(*,*) ' Saving IncBeam' 
          write(738)  (INCBEAM(i,1),i=1,OS)
          write(738)  (INCBEAM(i,2),i=1,OS)
          write(738)  (INCBEAM(i,3),i=1,OS)
         endif
!          write(*,*) iext 
         if (iext.ge.2) then 
          write(*,*) ' Saving IntField' 
          write(738)  (FIELD(i,1),i=1,OS)
          write(738)  (FIELD(i,2),i=1,OS)
          write(738)  (FIELD(i,3),i=1,OS)
         endif        
         write(*,*)

!         do i=1,10
!         write(*,*)  POS(i,1),POS(i,2), POS(i,3)    
!         enddo
!         do i=1,10
!         write(*,*)  RPOS(i,1),RPOS(i,2), RPOS(i,3)    
!         enddo

         close(UNIT=738)

        
         write(*,*) 'Checking ORDER'
          DX(1)=MINDIST(1)
          DX(2)=MINDIST(2)
          DX(3)=MINDIST(3)
!          write(*,*)  'Delta',DX(1:3)
         mad=0.d0 
         do i=1,OS
!          external target : start from 1
!          bxx=(POS(i,1)-((imax-imix)+2)/2-0.5d0)*DX(1)
!          byy=(POS(i,2)-((imay-imiy)+2)/2-0.5d0)*DX(2)
!          bzz=(POS(i,3)-((imaz-imiz)+2)/2-0.5d0)*DX(3)
!          read.geom : start from 0

!          bxx=(POS(i,1)-((imax-imix)+2)/2+0.5d0)*DX(1)
!          byy=(POS(i,2)-((imay-imiy)+2)/2+0.5d0)*DX(2)
!          bzz=(POS(i,3)-((imaz-imiz)+2)/2+0.5d0)*DX(3)
          bxx=(POS(i,1)-rcx)*DX(1)                                                                                                
          byy=(POS(i,2)-rcy)*DX(2)      
          bzz=(POS(i,3)-rcz)*DX(3) 
          diffx=dabs(bxx-RPOS(i,1))
          diffy=dabs(byy-RPOS(i,2))
          diffz=dabs(bzz-RPOS(i,3))
          mad=max(mad,diffx,diffy,diffz)
!          write(*,'(I4,2X,9F10.5)') i,
!     &   POS(i,1),POS(i,2),POS(i,3),
!     &   bxx,byy,bzz,
!     &   RPOS(i,1),RPOS(i,2),RPOS(i,3)
         enddo
         write(*,*) 'MAX ERROR',mad
         if (mad.gt.1.d-6) stop 'error order'
         write(*,*) 'Ok!'
!         write(*,*) 
!         write(*,*) 
!         do i=1,2
!         write(*,'(I8,6D18.10)')  i, 
!     & CXPOL(i,1),CXPOL(i,2), CXPOL(i,3)    
!         enddo
!         write(*,*) '....'
!         do i=OS-1,OS
!         write(*,'(I8,6D18.10)')  i, 
!     & CXPOL(i,1),CXPOL(i,2), CXPOL(i,3)    
!         enddo

         

         deallocate(POS)
         deallocate(RPOS)
         deallocate(CXPOL)
         if (iosb.eq.0) deallocate(INCBEAM)
         if (iosf.eq.0) deallocate(FIELD)

         return
 901     write(*,*) ' Cannot Open file ',fdipin
         iout=-1
         return

 900     write(*,*) ' ERROR!: Cannot Open file target'
         stop

 920     write(*,*) 'ERROR!: Unexepcted end of file reading dipoles'
         stop
         end
