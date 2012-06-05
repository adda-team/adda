subroutine cor3_limits( cor3_max, cor3_num, cor3, xmin, xave, xmax, ymin, yave, ymax, zmin, zave, zmax)
!
!*******************************************************************************
!
!! COR3_RANGE computes and prints the coordinate minima and maxima.
!
!
!  Modified:
!
!    22 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer COR3_MAX, the maximum number of points.
!
!    Input, integer COR3_NUM, the number of points.
!
!    Input, real COR3(3,COR3_MAX), the coordinates of points.
!
  integer cor3_max
!
  real cor3(3,cor3_max)
  integer cor3_num
  integer i
  real xave
  real xmax
  real xmin
  real yave
  real ymax
  real ymin
  real zave
  real zmax
  real zmin
!
  if ( cor3_num <= 0 ) then

    xave = 0.0
    xmax = 0.0
    xmin = 0.0

    yave = 0.0
    ymax = 0.0
    ymin = 0.0

    zave = 0.0
    zmax = 0.0
    zmin = 0.0

  else

    xave = cor3(1,1)
    xmax = cor3(1,1)
    xmin = cor3(1,1)

    yave = cor3(2,1)
    ymax = cor3(2,1)
    ymin = cor3(2,1)

    zave = cor3(3,1)
    zmax = cor3(3,1)
    zmin = cor3(3,1)

    do i = 2, cor3_num

      xave = xave + cor3(1,i)
      xmin = min ( xmin, cor3(1,i) )
      xmax = max ( xmax, cor3(1,i) )
  
      yave = yave + cor3(2,i)
      ymin = min ( ymin, cor3(2,i) )
      ymax = max ( ymax, cor3(2,i) )

      zave = zave + cor3(3,i)
      zmin = min ( zmin, cor3(3,i) )
      zmax = max ( zmax, cor3(3,i) )

    end do

    xave = xave / real ( cor3_num )
    yave = yave / real ( cor3_num )
    zave = zave / real ( cor3_num )

  end if

  write ( *, * ) ' '
  write ( *, * ) 'COR3_RANGE - Nodal coordinate range:'
  write ( *, * ) ' '
  write ( *, * ) '     Minimum     Average    Maximum       Range'
  write ( *, * ) ' '
  write ( *, '(a1,4g12.4)' ) 'X', xmin, xave, xmax, xmax-xmin
  write ( *, '(a1,4g12.4)' ) 'Y', ymin, yave, ymax, ymax-ymin
  write ( *, '(a1,4g12.4)' ) 'Z', zmin, zave, zmax, zmax-zmin
 
  return
end
