function arc_cosine ( c )

!*******************************************************************************
!
!! ARC_COSINE computes the arc cosine function, with argument truncation.
!
!  Discussion:
!
!    If you call geo_your system ACOS routine with an input argument that is
!    even slightly outside the range [-1.0, 1.0 ], you may get an unpleasant 
!    surprise (I did).
!
!    This routine simply truncates arguments outside the range.
!
!  Modified:
!
!    02 December 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) C, the argument.
!
!    Output, real ( kind = 8 ) ARC_COSINE, an angle whose cosine is C.
!
  implicit none

  real ( kind = 8 ) arc_cosine
  real ( kind = 8 ) c
  real ( kind = 8 ) c2

  c2 = c
  c2 = max ( c2, -1.0D+00 )
  c2 = min ( c2, +1.0D+00 )

  arc_cosine = acos ( c2 )

  return
end
subroutine geo_dvec_cross_3d ( v1, v2, v3 )

!*******************************************************************************
!
!! DVEC_CROSS_3D computes the cross product of two vectors in 3D.
!
!  Discussion:
!
!    The cross product in 3D can be regarded as the determinant of the
!    symbolic matrix:
!
!          |  i  j  k |
!      det | x1 y1 z1 |
!          | x2 y2 z2 |
!
!      = ( y1 * z2 - z1 * y2 ) * i
!      + ( z1 * x2 - x1 * z2 ) * j
!      + ( x1 * y2 - y1 * x2 ) * k
!
!  Modified:
!
!    07 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) V1(3), V2(3), the two vectors.
!
!    Output, real ( kind = 8 ) V3(3), the cross product vector.
!
  implicit none

  integer, parameter :: dim_num = 3

  real ( kind = 8 ) v1(dim_num)
  real ( kind = 8 ) v2(dim_num)
  real ( kind = 8 ) v3(dim_num)

  v3(1) = v1(2) * v2(3) - v1(3) * v2(2)
  v3(2) = v1(3) * v2(1) - v1(1) * v2(3)
  v3(3) = v1(1) * v2(2) - v1(2) * v2(1)

  return
end
function dvec_length ( dim_num, x )

!*******************************************************************************
!
!! DVEC_LENGTH returns the Euclidean length of a vector.
!
!  Modified:
!
!    08 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer DIM_NUM, the spatial dimension.
!
!    Input, real ( kind = 8 ) X(DIM_NUM), the vector.
!
!    Output, real ( kind = 8 ) DVEC_LENGTH, the Euclidean length of the vector.
!
  implicit none

  integer dim_num

  real ( kind = 8 ) dvec_length
  real ( kind = 8 ) x(dim_num)

  dvec_length = sqrt ( sum ( ( x(1:dim_num) )**2 ) )

  return
end
function dvec_triple_product ( v1, v2, v3 )

!*******************************************************************************
!
!! DVEC_TRIPLE_PRODUCT finds the triple product in 3D.
!
!  Discussion:
!
!    [A,B,C] = A dot ( B cross C ) 
!            = B dot ( C cross A )
!            = C dot ( A cross B )
!
!    The volume of a parallelepiped, whose sides are given by
!    vectors A, B, and C, is abs ( A dot ( B cross C ) ).
!
!    Three vectors are coplanar if and only if their triple product vanishes.
!
!  Modified:
!
!    07 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Eric Weisstein,
!    "Scalar Triple Product",
!    CRC Concise Encyclopedia of Mathematics, 1999
!
!  Parameters:
!
!    Input, real ( kind = 8 ) V1(3), V2(3), V3(3), the vectors.
!
!    Output, real ( kind = 8 ) DVEC_TRIPLE_PRODUCT, the triple product.
!
  implicit none

  integer, parameter :: dim_num = 3

  real ( kind = 8 ) dvec_triple_product
  real ( kind = 8 ) v1(dim_num)
  real ( kind = 8 ) v2(dim_num)
  real ( kind = 8 ) v3(dim_num)
  real ( kind = 8 ) v4(dim_num)

  call geo_dvec_cross_3d ( v2, v3, v4 )

  dvec_triple_product = dot_product ( v1(1:dim_num), v4(1:dim_num) )

  return
end
function geo_i_modp ( i, j )

!*******************************************************************************
!
!! geo_i_modp returns the nonnegative remainder of integer division.
!
!  Discussion:
!
!    If
!      NREM = geo_i_modp ( I, J )
!      NMULT = ( I - NREM ) / J
!    then
!      I = J * NMULT + NREM
!    where NREM is always nonnegative.
!
!    The MOD function computes a result with the same sign as the
!    quantity being divided.  Thus, suppose you had an angle A,
!    and you wanted to ensure that it was between 0 and 360.
!    Then mod(A,360) would do, if A was positive, but if A
!    was negative, your result would be between -360 and 0.
!
!    On the other hand, geo_i_modp(A,360) is between 0 and 360, always.
!
!  Examples:
!
!        I     J     MOD  geo_i_modp    Factorization
!
!      107    50       7       7    107 =  2 *  50 + 7
!      107   -50       7       7    107 = -2 * -50 + 7
!     -107    50      -7      43   -107 = -3 *  50 + 43
!     -107   -50      -7      43   -107 =  3 * -50 + 43
!
!  Modified:
!
!    02 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer I, the number to be divided.
!
!    Input, integer J, the number that divides I.
!
!    Output, integer geo_i_modp, the nonnegative remainder when I is
!    divided by J.
!
  implicit none

  integer i
  integer geo_i_modp
  integer j

  if ( j == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'geo_i_modp - Fatal error!'
    write ( *, '(a,i8)' ) '  geo_i_modp ( I, J ) called with J = ', j
    stop
  end if

  geo_i_modp = mod ( i, j )

  if ( geo_i_modp < 0 ) then
    geo_i_modp = geo_i_modp + abs ( j )
  end if

  return
end
function i_wrap ( ival, ilo, ihi )

!*******************************************************************************
!
!! I_WRAP forces an integer to lie between given limits by wrapping.
!
!  Example:
!
!    ILO = 4, IHI = 8
!
!    I  I_WRAP
!
!    -2     8
!    -1     4
!     0     5
!     1     6
!     2     7
!     3     8
!     4     4
!     5     5
!     6     6
!     7     7
!     8     8
!     9     4
!    10     5
!    11     6
!    12     7
!    13     8
!    14     4
!
!  Modified:
!
!    19 August 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer IVAL, an integer value.
!
!    Input, integer ILO, IHI, the desired bounds for the integer value.
!
!    Output, integer I_WRAP, a "wrapped" version of IVAL.
!
  implicit none

  integer geo_i_modp
  integer i_wrap
  integer ihi
  integer ilo
  integer ival
  integer jhi
  integer jlo
  integer wide

  jlo = min ( ilo, ihi )
  jhi = max ( ilo, ihi )

  wide = jhi - jlo + 1

  if ( wide == 1 ) then
    i_wrap = jlo
  else
    i_wrap = jlo + geo_i_modp ( ival - jlo, wide )
  end if

  return
end
subroutine geo_polygon_normal_3d ( n, v, normal ) 

!*********************************************************************
!
!! POLYGON_NORMAL_3D computes the normal vector to a polygon in 3D.
!
!  Discussion:
!
!    If the polygon is planar, then this calculation is correct.
!
!    Otherwise, the normal vector calculated is the simple average
!    of the normals defined by the planes of successive triples
!    of vertices.
!
!    If the polygon is "almost" planar, this is still acceptable.
!    But as the polygon is less and less planar, so this averaged normal
!    vector becomes more and more meaningless.
!
!  Modified:
!
!    12 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paulo Cezar Pinto Carvalho, Paulo Roma Cavalcanti,
!    Point in Polyhedron Testing Using Spherical Polygons,
!    in Graphics Gems V,
!    edited by Alan Paeth,
!    Academic Press, 1995, T385.G6975.
!
!  Parameters:
!
!    Input, integer N, the number of vertices.
!
!    Input, real ( kind = 8 ) V(3,N), the coordinates of the vertices.
!
!    Output, real ( kind = 8 ) NORMAL(3), the averaged normal vector
!    to the polygon. 
!
  implicit none

  integer, parameter :: dim_num = 3
  integer n

  real ( kind = 8 ) dvec_length
  integer i
  integer j
  real ( kind = 8 ) normal(dim_num)
  real ( kind = 8 ) normal_norm
  real ( kind = 8 ) p(dim_num)
  real ( kind = 8 ) v(dim_num,n)
  real ( kind = 8 ) v1(dim_num)
  real ( kind = 8 ) v2(dim_num)

  normal(1:dim_num) = 0.0D+00

  v1(1:dim_num) = v(1:dim_num,2) - v(1:dim_num,1)

  do j = 3, n

    v2(1:dim_num) = v(1:dim_num,j) - v(1:dim_num,1)

    call geo_dvec_cross_3d ( v1, v2, p )

    normal(1:dim_num) = normal(1:dim_num) + p(1:dim_num)

    v1(1:dim_num) = v2(1:dim_num)

  end do
!
!  Normalize.
!
  normal_norm = dvec_length ( dim_num, normal )

  if ( normal_norm == 0.0D+00 ) then
    return
  end if

  normal(1:dim_num) = normal(1:dim_num) / normal_norm

  return
end
subroutine geo_polygon_solid_angle_3d ( n, v, p, solid_angle )

!*********************************************************************
!
!! POLYGON_SOLID_ANGLE_3D: projected solid angle of a 3D plane polygon.
!
!  Discussion:
!
!    A point P is at the center of the unit sphere.  A planar polygon
!    is to be projected onto the surface of the unit sphere, by drawing
!    the ray from P to each polygonal vertex, and noting where this ray
!    intersects the unit sphere.  The area of the projected polygon is
!    equal to the solid angle, since we are considering the unit sphere.
!
!  Modified:
!
!    30 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paulo Cezar Pinto Carvalho, Paulo Roma Cavalcanti,
!    Point in Polyhedron Testing Using Spherical Polygons,
!    in Graphics Gems V,
!    edited by Alan Paeth,
!    Academic Press, 1995, T385.G6975.
!
!  Parameters:
!
!    Input, integer N, the number of vertices.
!
!    Input, real ( kind = 8 ) V(3,N), the coordinates of the vertices.
!
!    Input, real ( kind = 8 ) P(3), the point at the center of the unit sphere.
!
!    Output, double SOLID_ANGLE, the solid angle subtended
!    by the polygon, as projected onto the unit sphere around the point P.
!
  implicit none

  integer, parameter :: dim_num = 3
  integer n

  real ( kind = 8 ) a(dim_num)
  real ( kind = 8 ) angle
  real ( kind = 8 ) arc_cosine
  real ( kind = 8 ) area
  real ( kind = 8 ) b(dim_num)
  real ( kind = 8 ) dvec_length
  real ( kind = 8 ) dvec_triple_product
  integer i
  integer i_wrap
  integer j
  integer jp1
  real ( kind = 8 ) normal1(dim_num)
  real ( kind = 8 ) normal1_norm
  real ( kind = 8 ) normal2(dim_num)
  real ( kind = 8 ) normal2_norm
  real ( kind = 8 ) p(dim_num)
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) plane(dim_num)
  real ( kind = 8 ) r1(dim_num)
  real ( kind = 8 ) s
  real ( kind = 8 ) solid_angle
  real ( kind = 8 ) v(dim_num,n)

  if ( n < 3 ) then
    solid_angle = 0.0D+00
    return
  end if

  call geo_polygon_normal_3d ( n, v, plane )
 
  a(1:dim_num) = v(1:dim_num,n) - v(1:dim_num,1)

  area = 0.0D+00

  do j = 1, n

    r1(1:dim_num) = v(1:dim_num,j) - p(1:dim_num)

    jp1 = i_wrap ( j + 1, 1, n )

    b(1:dim_num) = v(1:dim_num,jp1) - v(1:dim_num,j)

    call geo_dvec_cross_3d ( a, r1, normal1 )

    normal1_norm = dvec_length ( dim_num, normal1 )

    call geo_dvec_cross_3d ( r1, b, normal2 )

    normal2_norm = dvec_length ( dim_num, normal2 )
    
    s = dot_product ( normal1(1:dim_num), normal2(1:dim_num) ) &
      / ( normal1_norm * normal2_norm )

    angle = arc_cosine ( s )

    s = dvec_triple_product ( b, a, plane )

    if ( 0.0D+00 < s ) then
      area = area + pi - angle
    else
      area = area + pi + angle
    end if

    a(1:dim_num) = -b(1:dim_num)

  end do

  area = area - pi * real ( n - 2, kind = 8 )

  if ( 0.0D+00 < dot_product ( plane(1:dim_num), r1(1:dim_num) ) ) then
    solid_angle = -area
  else
    solid_angle = area
  end if

  return
end
