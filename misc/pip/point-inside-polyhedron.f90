subroutine polyhedron_contains_point_3d ( node_num, face_num, face_order_max, &
  v, face_order, face_point, p, face_material, material_num, inside_mat )
!*********************************************************************
!
!! POLYHEDRON_CONTAINS_POINT_3D determines if a point is inside a polyhedron.
!  R. Schuh - October 2006
!
!  PIP (point inside polyhedron)
!
!  Discussion:
!
!    The reference states that the polyhedron should be simple (that
!    is, the faces should form a single connected surface), and that 
!    the individual faces should be consistently oriented.
!
!    However, the polyhedron does not, apparently, need to be convex.
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
!    Input, integer NODE_NUM, the number of vertices.
!
!    Input, integer FACE_NUM, the number of faces.
!
!    Input, integer FACE_ORDER_MAX, the maximum order of any face.
!
!    Input, real ( kind = 8 ) V(3,NODE_NUM), the coordinates of the vertices.
!
!    Input, integer FACE_ORDER(FACE_NUM), the order of each face.
!
!    Input, integer FACE_POINT(FACE_ORDER_MAX,FACE_NUM), the indices of the
!    nodes that make up each face.
!
!    Input, real ( kind = 8 ) P(3), the point to be tested.
!
!    Input, integer FACE_MATERIAL(FACE_NUM), the material of each face.
!
!    Input, integer MATERIAL_NUM, the number of materials.
!
!    Output, integer INSIDE_MAT is the component (domain) in which the point falls; zero if 
!    outside all domains
!
  implicit none

  integer, parameter :: dim_num = 3
  integer face_num
  integer face_order_max
  integer node_num

  real ( kind = 8 ) area(material_num)
  real ( kind = 8 ) wind_num
  integer face
  integer face_order(face_num)
  integer face_point(face_order_max,face_num)
  integer face_material(face_num)
  integer material_num
  integer mat
  integer inside_mat
  integer k
  integer node
  integer node_num_face
  real ( kind = 8 ) p(dim_num)
  real ( kind = 8 ), parameter :: four_pi = 12.566370614359173D+00
  real ( kind = 8 ) solid_angle(material_num)
  real ( kind = 8 ) v(dim_num,node_num)
  real ( kind = 8 ) v_face(dim_num,face_order_max)

  area = 0.0D+00

  do face = 1, face_num

    node_num_face = face_order(face)
    mat = face_material(face)

    do k = 1, node_num_face

      node = face_point(k,face)

      v_face(1:dim_num,k) = v(1:dim_num,node)

    end do

    call geo_polygon_solid_angle_3d ( node_num_face, v_face, p, solid_angle(mat) )

    area(mat) = area(mat) + solid_angle(mat)

  end do

  inside_mat = 0
  do mat = 1, material_num
!   winding number should be close to 0 (outside) or to 1 (inside)
    wind_num = dabs(area(mat)/four_pi)
    if ( wind_num > 0.01 ) then
      if ( wind_num < 0.99 ) then
        write(*,*)'WARNING: the surface with material #',mat,' seems to be not closed ', &
          'or the alignment of normals is inconsistent (when testing point ',p,')'
      else
        if (inside_mat /= 0) then
          write(*,*)'WARNING: interiors of the surfaces with materials #',inside_mat,' and #',mat, &
          ' intersect. Using the latter material for point',p
        endif
        inside_mat = mat
        if ( wind_num > 1.01 ) then
        write(*,*)'WARNING: the surface with material #',mat,' seems to consist of ', &
          'several inclosed parts with inconsistent alignment of normals (when testing point ',p,')'
        end if
      end if
    end if
  end do
  
  return
end