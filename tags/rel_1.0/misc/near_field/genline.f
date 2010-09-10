      program genline
!--------------------------------------------
!     Program to generate points on a line
!     input : 
!         x0 y0 z0
!         x1 y1 z1
!         min max NGP 
!
!    where x0,y0,z0 is the first point
!          x1,y1,z1 is the last point
!          min,max  define min,max values to label the points over a line
!          NGP is the number of point

!
!   Author: Fabio Della Sala, 2009
!   e-mail: fabio.dellasala@unisalento.it
!
!-------------------------------------------
      implicit none
      real*8 x,y,z,xmin,xmax
      real*8 x0,y0,z0,x1,y1,z1
      integer NGP,i
!      write(*,*) 'input x0,y0,z0'
      read(*,*) x0,y0,z0
!       write(*,*) 'input x1,y1,z1'
      read(*,*) x1,y1,z1
!      write(*,*)
      read(*,*) xmin,xmax,NGP
      write(*,*) 1
      write(*,*) xmin,xmax,NGP+1
      do i=0,NGP
       x=x0+(x1-x0)*i*1.d0/NGP
       y=y0+(y1-y0)*i*1.d0/NGP
       z=z0+(z1-z0)*i*1.d0/NGP
       write(*,'(3F15.5)') x,y,z
      enddo
      end
