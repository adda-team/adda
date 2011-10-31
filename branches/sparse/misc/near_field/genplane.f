      program genplane
!--------------------------------------------
!     Program to generate points on a plane with
!     equation x=val     (whatis=1)
!              y=val     (whatis=2)
!              z=val     (whatis=3)
!
!     input : 
!       val, whatis
!       xmin,xmax,NGPX
!       ymin,ymax,NGPY
!
!   xmin,xmax,ymin,ymax are the initial and final xy coordinate
!   NGPX NGPY are the number of points.
!
!   Author: Fabio Della Sala, 2009
!   e-mail: fabio.dellasala@unisalento.it
!
!-------------------------------------------
      implicit none
      real xpos,ypos,z
      real xmin,xmax,stepx,ymin,ymax,stepy
      integer NGPX,NGPY,i,j,whatis
      read(*,*) z, whatis
      read(*,*) xmin,xmax,NGPX
      read(*,*) ymin,ymax,NGPY
      write(*,*) 2
      write(*,*) NGPX+1,xmin,xmax
      write(*,*) NGPY+1,ymin,ymax
      write(*,*)  (NGPX+1)*(NGPY+1)
      stepx=(xmax-xmin)/NGPX
      stepy=(ymax-ymin)/NGPY
      do i=0,NGPX
       do j=0,NGPY
       xpos=xmin+i*stepx
       ypos=ymin+j*stepy
       if (whatis.eq.3) write(*,'(3F10.3)') xpos,ypos,z
       if (whatis.eq.2) write(*,'(3F10.3)') xpos,z,ypos
       if (whatis.eq.1) write(*,'(3F10.3)') z,xpos,ypos
      enddo
      enddo
      end
