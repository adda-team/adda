!*********************************************************************
      subroutine checksymm(imix,imax,imiy,imay,imiz,imaz,
     &                     pcx,pcy,pcz,OS,CXPOL,POS)
      implicit none
!********************************************************************
!      prepare array and check symmetry planes for dipoles
!      here dipole polarization is checked for symmetry
!      thus the input parameters refer to the particle (set of dipoles) 
!............param............................
      integer    imix,imax,imiy,imay,imiz,imaz
      integer    pcx,pcy,pcz
      integer    OS
      complex*16 CXPOL(OS,3)
      real*8     POS(OS,3)    
!.............local.....................
!      real*8     mazdev,tmpval
      integer   i,ix,iy,iz
      real*8 , allocatable :: DIPMAT(:,:,:)
!----------------         
      write(*,*)
      write(*,*) 'Checking Dipole Symmetry { '
      allocate(DIPMAT(imix:imax,imiy:imay,imiz:imaz))     
      write(*,*) imix,imax,imiy,imay,imiz,imaz
      DIPMAT=0.d0
      do i=1,OS
       ix=POS(i,1)
       iy=POS(i,2)  
       iz=POS(i,3)
       DIPMAT(ix,iy,iz)= (cdabs(CXPOL(i,1))**2+
     &                          cdabs(CXPOL(i,2))**2+
     &                          cdabs(CXPOL(i,3))**2)    *1.d18    
      enddo
      call  checkplanes(DIPMAT,imix,imax,imiy,imay,imiz,imaz,
     &                         pcx,pcy,pcz)
      deallocate(DIPMAT)        
      write(*,*) '} Checking Dipole Symmetry '     
      end subroutine 
!*****************************************************
!****************************************************
!*****************************************************
       subroutine checkplanes(MAT,imix,imax,imiy,imay,imiz,imaz,
     &                        pcx,pcy,pcz)
       implicit none
       real*8 compare
       external compare
!****************************************************8
!     CHECK for symmetry planes:
!      for dipoles imix-imaz refer to the particle (set of dipoles)
!      for fields imix-imaz refer to the domain in which they are calculated 
!      pcx-pcz always refer to the particle center (doubled) 
!.............param...............
       integer imix,imax,imiy,imay,imiz,imaz 
       real*8 MAT(imix:imax,imiy:imay,imiz:imaz)
       integer    pcx,pcy,pcz
!............local...............
       real*8 mazdev,tmpval
       integer ix,iy,iz,d1,d2,d3,off 
       integer icx,icy,icz
!----------- info --------------------
! pari:  [0....14][15   29]
! pcx/2:         14.5
! icx:           14
!        d1=15   d2=16 =>  d3=15  off=1
!
! dispari: [0....14][15][17....30]
! pcx/2:             15
! icx:               15
!         d1=16      d2=16> d3=16 off=0
!--------------------------------------
      write(*,*) ' CHECKING DIP SYMMETRY X -> -X' 
      icx=pcx/2;
      off=pcx-2*icx;
      mazdev=0.d0
      d1=icx-imix
      d2=imax-(icx+off)
      d3=min(d1,d2)
      write(*,*) 'offset',d1,d2,d3, off  
      if (d3.gt.0) then
       write(*,*) icx-(1-off)     ,'<->', icx-(d3)    
       write(*,*) icx+(1-off)+off ,'<->', icx+(d3)+off 
       do iz=imiz,imaz
       do iy=imiy,imay
! off=1  0,d3
! off=0  1,d3-1
       do ix=1-off,d3 
        tmpval= compare(MAT(icx-ix    ,iy,iz),
     &                  MAT(icx+ix+off,iy,iz)) 
        mazdev=max(mazdev,tmpval)
       enddo
       enddo
       enddo
       write(*,* ) ' ERROR_SX : ',mazdev
      endif
!--------------------check symm---------------------------
      write(*,*) ' CHECKING DIP SYMMETRY Y -> -Y' 
      icy=pcy/2;
      off=pcy-2*icy;
      mazdev=0.d0
      d1=icy-imiy
      d2=imay-(icy+off)
      d3=min(d1,d2)
      write(*,*) 'offset', d1,d2,d3,off
      if (d3.gt.0) then 
       write(*,*) icy-(1-off)     ,'<->', icy-(d3)     
       write(*,*) icy+(1-off)+off ,'<->', icy+(d3)+off 
       do iz=imiz,imaz
       do ix=imix,imax
       do iy=1-off,d3
       tmpval=compare(MAT(ix,icy-iy    ,iz),
     &                MAT(ix,icy+iy+off,iz))
       mazdev=max(mazdev,tmpval)
       enddo
       enddo
       enddo
       write(*,* )  ' ERROR_SY :',mazdev
      endif  
!--------------------check symm---------------------------
      write(*,*) ' CHECKING DIP SYMMETRY Z -> -Z' 
      icz=pcz/2;
      off=pcz-2*icz;
      mazdev=0.d0
      d1=icz-imiz
      d2=imaz-(icz+off)
      d3=min(d1,d2)
      write(*,*) 'offset', d1,d2,d3,off
      if (d3.gt.0) then
       write(*,*) icz-(1-off)    ,'<->', icz-(d3)                                                                         
       write(*,*) icz+(1-off)+off,'<->', icz+(d3)+off 
       do iy=imiy,imay
       do ix=imix,imax
       do iz=1-off,d3
        tmpval= compare(MAT(ix,iy,icz-iz    ),
     &                  MAT(ix,iy,icz+iz+off))
        mazdev=max(mazdev,tmpval)
       enddo
       enddo
       enddo
       write(*,* ) ' ERROR_SZ :',mazdev
      endif
!------------------------------------
      end subroutine 


!************************************************************
!************************************************************
!************************************************************

       subroutine testbeam(rcx,rcy,rcz,debug,
     &                     POS,OS,INCBEAM,dir,pol,AKD,DX) 
       implicit none
       external comparec
       real*8   comparec
!...............param...............   
!       integer    imax,imix,imay,imiy,imaz,imiz
       integer    debug
       real*8     rcx,rcy,rcz
       integer    OS
       real*8     POS(OS,3)  
       complex*16 INCBEAM(OS,3)
       real*8     dir(3),pol(3)
       real*8     AKD,DX(3)
!...............local..................
       real*8 mazdev
       integer T
       complex*16 CXUNIT
       real*8 bxx,byy,bzz
       complex*16 EE(3)
!............................... 
       CXUNIT=(0.d0,1.d0)
       write(*,*)  
       write(*,*) 'Test Beam {' 
       write(*,*) ' dir',dir
       write(*,*) ' pol',pol
       mazdev=0.d0
       if (debug.ge.3) 
     & open(UNIT=2311,FILE='BEAM_CHECK',STATUS='UNKNOWN')
       DO T=1,OS
    
         call calcbeam(POS(T,1),POS(T,2),POS(T,3),DX,
     &                 rcx,rcy,rcz,dir,pol,AKD,EE)

       bxx=(POS(T,1)-rcx)*DX(1)
       byy=(POS(T,2)-rcy)*DX(2)
       bzz=(POS(T,3)-rcz)*DX(3)  
         if (debug.ge.3)  write(2311,'(15D12.3)')  
     &                    bxx,
     &                    byy,
     &                    bzz,
     &                    EE(1),INCBEAM(T,1),
     &                    EE(2),INCBEAM(T,2),
     &                    EE(3),INCBEAM(T,3)

          mazdev=max(mazdev, comparec(EE(1),INCBEAM(T,1))+ 
     &                       comparec(EE(2),INCBEAM(T,2))+ 
     &                       comparec(EE(3),INCBEAM(T,3))  )


       enddo
       write(*,*) ' ERROR_Beam : ', mazdev
       if (debug.ge.2) close(2311)
       write(*,*) '} Test Beam'
       end subroutine

!***********************************************************
!***********************************************************
       subroutine calcbeam(x,y,z,dx,rcx,rcy,rcz,dir,pol,AKD,
     &  EE) 
       implicit none
!...............param................        
       real*8     rcx,rcy,rcz
       real*8     x,y,z 
       real*8     dir(3),pol(3)
       real*8     AKD,DX(3)
       complex*16 EE(3) ! output
!...............local..................
       real*8 bxx,byy,bzz,mydot  
       complex*16 CXUNIT,tmp  
       CXUNIT=(0.d0,1.d0)
!         bxx=(POS(T,1)-((imax-imix)+2)/2-0.5d0)*DX(1)
!         byy=(POS(T,2)-((imay-imiy)+2)/2-0.5d0)*DX(2)
!         bzz=(POS(T,3)-((imaz-imiz)+2)/2-0.5d0)*DX(3)
       bxx=(x-rcx)*DX(1)
       byy=(y-rcy)*DX(2)
       bzz=(z-rcz)*DX(3)     
       mydot=bxx*dir(1)+byy*dir(2)+bzz*dir(3)
       tmp=cdexp(CXUNIT*AKD*mydot)
       EE(1)=pol(1)*tmp
       EE(2)=pol(2)*tmp
       EE(3)=pol(3)*tmp
       end subroutine


!**********************************************************
!**********************************************************
!**********************************************************

        subroutine testinternalfield(FIELD,CXPOL,INCBEAM,OS,
     &                               AKD,alpha,succe,vol)
        implicit none
        real*8 FOURPI
        parameter(fourpi=12.566370614359172953850573533118d0)
!...............param..............
        integer OS      
        complex*16 CXPOL(OS,3)
        complex*16 INCBEAM(OS,3)
        complex*16 FIELD(OS,3)
        complex*16 alpha,succe
        real*8 vol,AKD
!..................local...........
        real*8   normsquared
        external normsquared
        complex*16 qqq,ttt
        integer T
        complex*16 MQX,MQY,MQZ
        complex*16 VVX,VVY,VVZ
!        complex*16 DDX,DDY,DDZ
        complex*16 QQZ,QQY,QQX
        real*8 CABS,CABSt,CEXT,CEXTt
        real*8 CABS2,CABS2t,CABS3t,p2
        real*8 errx1,
     &         errx2,
     &         errx3,
     &         errx4
        write(*,*) 
        write(*,*) 'Debug InternalField {'
        qqq=  - FOURPI/3.d0 * succe
        ttt= 1+ FOURPI/3.d0 * succe


        write(*,*) ' AKD      :',AKD
        write(*,*) ' vol      :',vol
        write(*,*) ' alpha^-1 :',1.d0/alpha
        write(*,*) ' chi      :',succe
        write(*,*) ' epsilon  :',ttt


        errx1=-1.d9
        errx2=-1.d9
        errx3=-1.d9
        errx4=-1.d9

!        open(UNIT=2307,FILE='EXC',STATUS='UNKNOWN')          ! exc
!        open(UNIT=2308,FILE='ESELF',STATUS='UNKNOWN')         ! eself
!!        open(UNIT=2309,FILE='ETOT2_CHECK',STATUS='UNKNOWN')   ! (ESELF+EXC)/ETOT  must be one
!        open(UNIT=2301,FILE='DIP',STATUS='UNKNOWN')           ! CONFRONTO
!        open(UNIT=2303,FILE='DIPSELF',STATUS='UNKNOWN')       ! dipsel
!!        open(UNIT=2304,FILE='ESELF_CHECK',STATUS='UNKNOWN')   ! file ESELFCHECK  ! must be one
!!        open(UNIT=2305,FILE='INTFIE_CHECK',STATUS='UNKNOWN')  ! must be one
!!        open(UNIT=2306,FILE='ETOT_CHECK',STATUS='UNKNOWN')    ! must be one

        CABS2=0.d0
        CABS=0.d0
        CEXT=0.d0
        DO T=1,OS
         p2= normsquared( CXPOL(T,1), CXPOL(T,2),CXPOL(T,3)) 

! ABS_FINITE: IM[P E*] = IM [P invchiV* P* ] = IM[invchiV*] P^2 = -IM[invchiV] P^2 
         CABSt=dimag( CXPOL(T,1)*dconjg(FIELD(T,1))+
     &                CXPOL(T,2)*dconjg(FIELD(T,2))+
     &                CXPOL(T,3)*dconjg(FIELD(T,3)) )
         CABSt=-dimag(1.d0/(succe*vol))*p2
 
         CABS=CABS+ CABSt

         CEXTt=dimag( CXPOL(T,1)*dconjg(INCBEAM(T,1))+
     &                CXPOL(T,2)*dconjg(INCBEAM(T,2))+
     &                CXPOL(T,3)*dconjg(INCBEAM(T,3)) )        
         CEXT=CEXT + CEXTt 

! ABS_DRAINE: IM[P Exc*] = IM [P invalpha* P* ] = IM[invalpha*] P^2 = -IM[invalpha] P^2    
         VVX=1.d0*(CXPOL(T,1))/alpha ! EXC=  a^-1 P 
         VVY=1.d0*(CXPOL(T,2))/alpha
         VVZ=1.d0*(CXPOL(T,3))/alpha
!         write(2307,'(6F15.8)') VVX,VVY,VVZ  ! EXC 

         CABS2t=dimag(dconjg(VVX)*CXPOL(T,1)+
     &                dconjg(VVY)*CXPOL(T,2)+
     &                dconjg(VVZ)*CXPOL(T,3) )

         CABS2t= -dimag(1.d0/(alpha))  * p2  
         CABS3t= -((2.d0/3.d0)*(AKD**3)) * p2 
                                   
         CABS2=CABS2 +CABS2t+CABS3t
!         write(989,*) T,CABS2*FOURPI*AKD,CABS2t+CABS3t           
        
!         ETOT= CXPOL/(V*succe)
!         EINC-a^-1 P = -EDIP= -sum_i<>j Gij Pj
!         ETOT = EINC + EDIP + ESELF
        

!                      ETOT / [CXPOL/(V*succe)]
         MQX=1.d0*FIELD(T,1) /   ( (CXPOL(T,1)/(succe*vol)) ) 
         MQY=1.d0*FIELD(T,2) /   ( (CXPOL(T,2)/(succe*vol)) )
         MQZ=1.d0*FIELD(T,3) /   ( (CXPOL(T,3)/(succe*vol)) )


         errx1=max(errx1,dabs(cdabs(MQX)+cdabs(MQY)+cdabs(MQZ)-3.d0))
!!        write(2305,'(6F15.8)') MQX,MQY,MQZ  ! MUST BE ONE => INTFIE_CHECK



      

 !        QQX=qqq* 1.d0*FIELD(T,1)   ! ESELF= -L x Etot
 !        QQY=qqq* 1.d0*FIELD(T,2)
 !        QQZ=qqq* 1.d0*FIELD(T,3)
 !        write(2308,'(6F15.8)') QQX,QQY,QQZ  ! ESELF

         QQX=(qqq* 1.d0*FIELD(T,1)+ 1.d0*(CXPOL(T,1))/alpha)/
     &     (1.d0*FIELD(T,1))
         QQY=(qqq* 1.d0*FIELD(T,2)+ 1.d0*(CXPOL(T,2))/alpha)/
     &     (1.d0*FIELD(T,2))
         QQZ=(qqq* 1.d0*FIELD(T,3)+ 1.d0*(CXPOL(T,3))/alpha)/
     &     (1.d0*FIELD(T,3)) 
!!        write(2309,'(6F15.8)') QQX,QQY,QQZ  ! (ESELF+EXC)/ETOT => ETOT2_CHECK

         errx2=max(errx2,dabs(cdabs(QQX)+cdabs(QQY)+cdabs(QQZ)-3.d0))

!         VVX=1.d0*INCBEAM(T,1) -   CXPOL(T,1)/alpha ! EINC - a^-1 P  = - EDIP
!         VVY=1.d0*INCBEAM(T,2) -   CXPOL(T,2)/alpha
!         VVZ=1.d0*INCBEAM(T,3)  -  CXPOL(T,3)/alpha
!         write(2301,'(6F15.8)') VVX,VVY,VVZ  ! CONFRONTO dip ! dip is compute by near


!         DDX=FIELD(T,1)-INCBEAM(T,1) !ETOT-EINC
!         DDY=FIELD(T,2)-INCBEAM(T,2) !ETOT-EINC
!         DDZ=FIELD(T,3)-INCBEAM(T,3) !ETOT-EINC
!         write(2303,'(6F15.8)') DDX,DDY,DDZ !! file='dipself' ETOT-EINC

!                           a^-1 P =  (1+L*x) * Etot
         QQX =  (CXPOL(T,1)/alpha) / ( ttt * FIELD(T,1)) 
         QQY =  (CXPOL(T,2)/alpha) / ( ttt * FIELD(T,2))  
         QQZ =  (CXPOL(T,3)/alpha) / ( ttt * FIELD(T,3))  
!!         write(2306,'(6F15.8)') QQX,QQY,QQZ  ! ETOTCHECK

         errx3=max(errx3,dabs(cdabs(QQX)+cdabs(QQY)+cdabs(QQZ)-3.d0))

         QQX=(1.d0*FIELD(T,1) -   CXPOL(T,1)/alpha)
     &       /( qqq* 1.d0*FIELD(T,1) )  ! ETOT-a^-1 P = ESELF 
         QQY=(1.d0*FIELD(T,2) -   CXPOL(T,2)/alpha) 
     &       /( qqq* 1.d0*FIELD(T,2) )  ! ETOT-a^-1 P = ESELF 
         QQZ=(1.d0*FIELD(T,3) -   CXPOL(T,3)/alpha)
     &       /( qqq* 1.d0*FIELD(T,3) )  ! ETOT-a^-1 P = ESELF 
!!          write(2304,'(6F15.8)') QQX,QQY,QQZ  ! file ESELFCHECK  => ESELF_CHECK'

         errx4=max(errx4,dabs(cdabs(QQX)+cdabs(QQY)+cdabs(QQZ)-3.d0))

! a^-1 P = Exc = Etot- Eself = Etot - ( -L)*x*Etot = (1+L*x) * Etot
!  a^-1 P = (1+L*x) Etot
!
! -  a^-1 P = -  (1+L*x) Etot
! Etot -  a^-1 P = -L*x Etot
! etot -a^-1 P  = Eself = - L x Etot


        
        ENDDO
!        close(2301)
!        close(2303)
!        close(2304)
!        close(2305)
!        close(2306)
!        close(2307)
!        close(2308)
!        close(2309)
        write(*,*) ' ERROR_F1 :',errx1
        write(*,*) ' ERROR_F2 :',errx2
        write(*,*) ' ERROR_F3 :',errx3
        write(*,*) ' ERROR_F4 :',errx4

        CABS =CABS*FOURPI*AKD
        CEXT =CEXT*FOURPI*AKD 
        CABS2=CABS2*FOURPI*AKD
        write(*,*) ' CABS_finite : ',CABS
        write(*,*) ' CABS_draine : ',CABS2
        write(*,*) '        CEXT : ',CEXT    
        write(*,*) '} Debug InternalField '
        write(*,*) 
        end subroutine
!**********************************************
!**********************************************

       real*8 function normsquared(ex,ey,ez)
       implicit none
       complex*16 ex,ey,ez 
       normsquared=cdabs(ex)**2+cdabs(ey)**2+cdabs(ez)**2
       end function

!**************************************************
       real*8 function compare(x,y)
       implicit none
       real*8 x,y,comparet
       if (dabs(x).lt.1.d-99) then
        if (dabs(y).lt.1.d-99) then
!       x zero y zero
        comparet=0.d0
        else
!       x zero y not zero
        comparet=dabs(x)
        endif
       else if (dabs(y).lt.1.d-99) then
!      x not zero y zero
        comparet=dabs(y)
       else 
!      x not zero y not zero
       comparet=(x-y)/((dabs(x)+dabs(y))/2.d0)
       endif
       compare=comparet
       end function

!**************************************************
       real*8 function comparec(x,y)
       implicit none
       complex*16 x,y
       real*8 comparet
       if (cdabs(x).lt.1.d-99) then
        if (cdabs(y).lt.1.d-99) then
!       x zero y zero
        comparet=0.d0
        else
!       x zero y not zero
        comparet=cdabs(x)
        endif
       else if (cdabs(y).lt.1.d-99) then
!      x not zero y zero
        comparet=cdabs(y)
       else 
!      x not zero y not zero
       comparet=cdabs(x-y)/((cdabs(x)+cdabs(y))/2.d0)
       endif
       comparec=comparet
       end function
