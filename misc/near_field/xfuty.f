!-------------------------------------------------
! Support routines to plot 2D data using XFarbe
!
! Author: Fabio Della Sala, 2009
! e-mail: fabio.dellasala@unisalento.it
!...................................................
        subroutine prepxf(iunit,namefull,minv,maxv,shift)
        implicit none
        double precision maxv,minv
        integer            iunit,shift,pal
        character*(*)      namefull
        
        maxv=-1.d50
        minv=1.d50
        shift=0
        open(UNIT=iunit,FILE=namefull,status='UNKNOWN') 
        write(iunit,'(A70)') namefull    
        end subroutine

!...................................................
       subroutine postxf(iunit,namefull,minv,maxv,shift,pal,levels,
     & d1,d2,d3,d4)
        implicit none 
        double precision maxv,minv
        integer            iunit,shift,pal,levels
        character*(*)      namefull
!        -------------------------------
         double precision count(2*levels+1),reap 
         double precision d1,d2,d3,d4
         integer j
         write(*,'(A,A30,2D18.10,I2)')
     &  'saved file ',namefull ,minv,maxv,pal

        call enlarge(maxv,minv)

        if (pal.eq.0) then
!        if (maxv.gt.0.d0.and.minv.lt.0.d0) then
         maxv=max(maxv,dabs(minv))
         minv=-maxv
!        if (maxv.le.0.d0) maxv=-minv  ! only negative values
!        if (minv.ge.0.d0) minv=-maxv  ! only positive values 
           
!        write(*,*) minv,maxv

!        if (maxr.gt.1.and.maxz.gt.1) then  
         
           do j=1,levels
            reap=dble(j-1)/dble(levels)
            count(j)= minv+dabs(minv)*reap ! da -vmax a o   
           enddo

           count(levels+1)=0.0

           do j=(levels+1)+1,2*levels+1
            reap=dble(j-(levels+1))/dble(levels) 
            count(j)= reap*dabs(maxv)  ! da o a vmax
           enddo
         
       else
!            write(*,*) minv,maxv
           do j=1,2*levels+1
            reap=dble(j-1)/dble(2*levels)
            count(j)= minv+dabs(maxv-minv)*reap    
           enddo
       endif

 
           write(iunit,*) 2*levels+1
           do j=1,2*levels+1
           write(iunit,'(I5,E20.10)') j+shift,count(j)
           enddo
 

          write(iunit,*) 2*levels+1+shift+1
          write(iunit,*) d1,d2!-dr*bohr/10.0*maxrn,+dr*bohr/10.0*maxrn
          write(iunit,*) d3,d4!dz/2*bohr/10.d0,(dz/2+(maxz-1)*dz)*bohr/10.d0 
!        endif

         close(iunit)

        end subroutine

!---------------------------------------------------------------
         subroutine enlarge(zmaxv,zminv)
         implicit none
         double precision zmaxv,zminv
         double precision maxv,minv

!.........increase a lillte.................
         if (zmaxv.gt.0.d0)    then
           maxv=zmaxv*1.15d0
         else if (zmaxv.eq.0.d0) then
           maxv=-zminv*0.1d0
         else
           maxv=zmaxv*0.85d0
         endif

         if (zminv.lt.0.d0) then
          minv=zminv*1.15d0
         else if (zminv.eq.0.d0) then
          minv=-zmaxv*0.1d0
         else
          minv=zminv*0.85d0
         endif

         zmaxv=maxv
         zminv=minv

         end subroutine

