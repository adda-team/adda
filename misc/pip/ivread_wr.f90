subroutine ivread_wr(filein_name,face_normal,face_area,face_num,cor3_num, &
  cor3,face,order_max,face_order,face_material,material_num,vertex_range) 
! 
!  ivread.f90  Wriedt 04 July 2000
!
!*******************************************************************************
!
!! IVREAD is the main program for converting a graphics file.
!
!
!  Introduction:
!
!    This program was written for the Studio for Creative Inquiry.
!
!    This program was originally intended to read an Inventor 3D graphics
!    file, and write a corresponding VLA graphics file.  
!
!    However, it can also try to read files of several different types,
!    convert the data internally to line data, and write out files
!    of various types.
!
!    DXF and Inventor file formats are sophisticated and extensive.
!    This program was only designed to work on a few simple cases, and
!    may easily be confused by unexpected (though legal) data.
!
!    A sketch of the file formats and data structures is included
!    in most of the READ and WRITE routines.
!
!    The sizes of COR3_MAX, FACE_MAX and LINE_MAX control how much 
!    information this program can handle. 
!
!    The "head.25.iv" file has 240146 faces and 119,736 points.
!    The "fish.iv" file has about 76000 faces.
!    The "brain.iv" file has between 40000 and 75000 faces.
!
!  Development:
!
!    26 May 1999: Added LINE_PRUNE switch, which will try to cut down
!    (by about half) the number of superfluous lines created when
!    faces are turned into lines by FACE_TO_LINE for VLA_WRITE output.
!
!    22 May 1999: For VLA output files, the program now will automatically 
!    try to temporarily convert all face information to line information.  No
!    sophisticated attempt is made to delete superfluous lines (the
!    way FACE_TO_EDGE tries.)
!
!    The "<<" merge command:
!
!      On 20 April 1999, the "<<" command was added, to allow data
!      from two or more files to be merged.  It works, on a simple example,
!      for OBJ files.  However, some tuning of OBJ_READ was necessary.
!      Similar testing and tuning must be done to the other READ routines
!      before they will work with this option.
!      On 21 April 1999, the "<<" command worked on a simple example
!      using two IV files as input.  SMF_READ was also updated for the "<<"
!      command, but not tested.  ASE_READ, DXF_READ, HRC_READ,
!      STLA_READ and VLA_READ may already be OK.
!      On 22 April 1999, "<<" command works with ASE_READ, HRC_READ,
!      SMF_READ and STLA_READ.
!
!    The "MatrixTransform" field in IV_READ/IV_WRITE.
!      I'm having problems because I am reading in an IV file that has
!      a matrix transform that is not the identity.  I just ignore it,
!      and so my data is not rotated and scaled, when it should be.
!      As a start toward addressing this issue, I have IV_WRITE 
!      writing out the current transform matrix.  One problem with
!      Inventor is that the transform matrix can be specified on
!      every level, and the actual transform matrix that applies
!      has to be deduced from where you are in the tree.
!      Right now, all I've done is have IV_READ read the matrix,
!      and IV_WRITE write it out.  No concatenation is possible right now,
!      but my kludgy code will at least apply ONE transformation matrix
!      to the data, in IV_READ, anyway...
!
!    Adding material/normal/texture binding stubs, because 
!    A) new SMF format allows it;
!    B) Inventor uses it;
!    C) SCI wants to do textures eventually.
!
!    SMF_READ and SMF_WRITE can now read and write face and vertex colors 
!    of SMF2.0 files.
!
!  Modified:
!
!    26 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Parameter, integer COR3_MAX, the maximum number of points.
!
!    Parameter, integer EDGE_MAX, the maximum number of edges.
!
!    Parameter, integer FACE_MAX, the maximum number of faces.
!
!    Parameter, integer LINE_MAX, the maximum number of line definition items.
!
!    Parameter, integer MATERIAL_MAX, the maximum number of materials.
!
!    Parameter, integer ORDER_MAX, the maximum number of vertices per face.
!
!    Parameter, integer TEXTURE_MAX, the maximum number of textures.
!
  implicit none
  integer, parameter :: cor3_max = 100000
  integer, parameter :: edge_max = 100
  integer, parameter :: face_max = 100000
  integer, parameter :: line_max = 100000
  integer, parameter :: material_max = 200
! yurkin - order_max is now an argument
  integer order_max
  integer, parameter :: texture_max = 10
!
  real cor3(3,cor3_max)
  integer cor3_material(cor3_max)
  real cor3_normal(3,cor3_max)
  integer cor3_num

  real cor3_tex_uv(2,cor3_max)
  logical debug
  integer face(order_max,face_max)
  real face_area(face_max)
  integer face_material(face_max)
  real face_normal(3,face_max)
  integer face_order(face_max)
  real face_tex_uv(2,face_max)
  character ( len = 100 ) filein_name
  character ( len = 100 ) fileout_name
!*thw 29.11.00! 
  integer ierror
  integer line_dex(line_max)
  integer line_material(line_max)
  integer line_prune
  character ( len = 100 ) material_name(material_max)
  real material_rgba(4,material_max)
  real normal_temp(3,order_max*face_max)
  character ( len = 100 ) object_name
  character ( len = 100 ) texture_name(texture_max)
  real texture_temp(2,order_max*face_max)
  real transform_matrix(4,4)
  integer vertex_material(order_max,face_max)
  real vertex_normal(3,order_max,face_max)
  real vertex_tex_uv(2,order_max,face_max)
  real vertex_range(6)
!*thw 29.11.00! 
  integer face_num
  character ( len = 10 ) filein_type
  integer material_num
!
!  Initialize a few program variables.
!
  debug = .false.
! filein_name = ' '
  fileout_name = ' '
  ierror = 0
  line_prune = 1
!
  call command_line ( cor3, cor3_material, cor3_normal, cor3_tex_uv, &
    debug, face, face_area, face_material, face_normal, face_order, face_tex_uv, &
    filein_name, ierror, line_dex, line_material, material_name, material_rgba, &
    cor3_max, face_max, line_max, material_max, order_max, texture_max, &
    normal_temp, object_name, texture_name, texture_temp, transform_matrix, &
    vertex_material, vertex_normal, vertex_tex_uv, face_num, filein_type, &
    cor3_num, material_num, vertex_range)

  return
end
subroutine c_cap ( c )
!
!*******************************************************************************
!
!! C_CAP capitalizes a single character.
!
!
!  Modified:
!
!    19 July 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character C, the character to capitalize.
!
  implicit none
  character c
  integer itemp
!
  itemp = ichar ( c )

  if ( 97 <= itemp .and. itemp <= 122 ) then
    c = char ( itemp - 32 )
  end if

  return
end
function c_eqi ( c1, c2 )
!
!*******************************************************************************
!
!! C_EQI is a case insensitive comparison of two characters for equality.
!
!
!  Examples:
!
!    C_EQI ( 'A', 'a' ) is .TRUE.
!
!  Modified:
!
!    14 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character C1, C2, the characters to compare.
!
!    Output, logical C_EQI, the result of the comparison.
!
  implicit none
  logical c_eqi
  character c1
  character c2
  character cc1
  character cc2
!
  cc1 = c1
  cc2 = c2

  call c_cap ( cc1 )
  call c_cap ( cc2 )

  if ( cc1 == cc2 ) then
    c_eqi = .true.
  else
    c_eqi = .false.
  end if

  return
end
function c_is_control ( c )
!
!*******************************************************************************
!
!! C_IS_CONTROL reports whether a character is a control character or not.
!
!
!  Definition:
!
!    A "control character" has ASCII code <= 31 or ASCII code => 127.
!
!  Modified:
!
!    14 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character C, the character to be tested.
!
!    Output, logical C_IS_CONTROL, TRUE if C is a control character, and
!    FALSE otherwise.
!
  implicit none
  character c
  logical c_is_control
!
  if ( ichar ( c ) <= 31 .or. ichar ( c ) >= 127 ) then
    c_is_control = .true.
  else
    c_is_control = .false.
  end if

  return
end
subroutine c_to_digit ( c, digit )
!
!*******************************************************************************
!
!! C_TO_DIGIT returns the integer value of a base 10 digit.
!
!
!  Example:
!
!     C   DIGIT
!    ---  -----
!    '0'    0
!    '1'    1
!    ...  ...
!    '9'    9
!    ' '    0
!    'X'   -1
!
!  Modified:
!
!    04 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character C, the decimal digit, '0' through '9' or blank
!    are legal.
!
!    Output, integer DIGIT, the corresponding integer value.  If C was
!    'illegal', then DIGIT is -1.
!
  implicit none
  character c
  integer digit
!
  if ( lge ( c, '0' ) .and. lle ( c, '9' ) ) then

    digit = ichar ( c ) - 48

  else if ( c == ' ' ) then

    digit = 0

  else

    digit = -1

  end if

  return
end
subroutine command_line ( cor3, cor3_material, cor3_normal, cor3_tex_uv, &
  debug, face, face_area, face_material, face_normal, face_order, face_tex_uv, &
  filein_name, ierror, line_dex, line_material, material_name, material_rgba, &
  cor3_max, face_max, line_max, material_max, order_max, texture_max, &
  normal_temp, object_name, texture_name, texture_temp, transform_matrix, &
  vertex_material, vertex_normal, vertex_tex_uv, face_num, filein_type, &
  cor3_num, material_num, vertex_range)
!
!*******************************************************************************
!
!! COMMAND_LINE works with command line parameters.
!
!
!  Discussion:
!
!    This routine is invoked when the user command is something like
!
!      ivread filein_name fileout_name
!
!    or
!
!      ivread -rn filein_name fileout_name
!
!    where "-rn" signals the "reverse normals" option, or
!
!      ivread -rf filein_name fileout_name
!
!    where "-rf" signals the "reverse faces" option.
!
!  Modified:
!
!    26 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Workspace, real COR3(3,COR3_MAX), the coordinates of points.
!
!    Workspace, integer COR3_MATERIAL(COR3_MAX), the material of each node.
!
!    Workspace, real COR3_NORMAL(3,COR3_MAX), normals at nodes.
!
!    Workspace, real COR3_TEMP(COR3_MAX).
!
!    Workspace, real COR3_TEX_UV(2,COR3_MAX), UV texture coordinates for nodes.
!
!    Input, logical DEBUG, debugging switch.
!
!    Workspace, integer FACE(ORDER_MAX,FACE_MAX), the nodes making faces.
!
!    Workspace, real FACE_AREA(FACE_MAX), area of each face.
!
!    Workspace, integer FACE_MATERIAL(FACE_MAX), the material of each face.
!
!    Workspace, real FACE_NORMAL(3,FACE_MAX), the normal vector at each face.
!
!    Workspace, integer FACE_ORDER(FACE_MAX), the number of vertices per face.
!
!    Workspace, real FACE_TEX_UV(2,FACE_MAX), face texture coordinates.
!
!    Workspace, character ( len = 100 ) FILEIN_NAME, the name of the input file.
!
!    Output, integer IERROR, error flag.
!
!    Workspace, integer LINE_DEX(LINE_MAX), nodes forming a line, terminated by -1.
!
!    Workspace, integer LINE_MAT(LINE_MAX), material index for line.
!
!    Workspace, character ( len = 100 ) MATERIAL_NAME(MATERIAL_MAX), material names.
!
!    Workspace, real MATERIAL_RGBA(4,MATERIAL_MAX), material R, G, B and A values.
!
!    Input, integer COR3_MAX, the maximum number of points.
!
!    Input, integer FACE_MAX, the maximum number of faces.
!
!    Input, integer LINE_MAX, the maximum number of line definition items.
!
!    Input, integer MATERIAL_MAX, the maximum number of materials.
!
!    Input, integer ORDER_MAX, the maximum number of vertices per face.
!
!    Input, integer TEXTURE_MAX, the maximum number of textures.
!
!    Workspace, real NORMAL_TEMP(3,ORDER_MAX*FACE_MAX).
!!
!    Workspace, character ( len = 100 ) OBJECT_NAME, object name.
!
!    Workspace, character ( len = 100 ) TEXTURE_NAME(TEXTURE_MAX), texture names.
!
!    Workspace, real TRANSFORM_MATRIX(4,4), the transformation matrix.
!
!    Workspace, integer VERTEX_MAT(ORDER_MAX,FACE_MAX), vertex materials.
!
!    Workspace, real VERTEX_NORMAL(3,ORDER_MAX,FACE_MAX), normals at vertices.  
!
!    Workspace, real VERTEX_TEX_UV(2,ORDER_MAX,FACE_MAX), vertex texture coordinates.
!
!    ....
!
!    Output, real VERTEX_RANGE(6), min & max values for vertex coordinate
!
!*thw 29.11.00! 
! use DFLIB
! use DFPORT 

  implicit none
  integer cor3_max
  integer face_max
  integer line_max
  integer material_max
  integer order_max
  integer texture_max
!
  integer color_num
  real cor3(3,cor3_max)
  integer cor3_material(cor3_max)
  real cor3_normal(3,cor3_max)
  integer cor3_num
  real cor3_tex_uv(2,cor3_max)
  logical debug
  integer face(order_max,face_max)
  real face_area(face_max)
  integer face_material(face_max)
  real face_normal(3,face_max)
  integer face_num
  integer face_order(face_max)
  real face_tex_uv(2,face_max)
  character ( len = 100 ) filein_name
  character ( len = 10 ) filein_type
  integer group_num
  integer ierror
  integer iseed
  integer line_dex(line_max)
  integer line_material(line_max)
  integer line_num
  character ( len = 100 ) material_name(material_max)
  integer material_num
  real material_rgba(4,material_max)
  real normal_temp(3,order_max*face_max)
  character ( len = 100 ) object_name
  integer object_num
  logical reverse_faces
  logical reverse_normals
  character ( len = 100 ) texture_name(texture_max)
  integer texture_num
  real texture_temp(2,order_max*face_max)
  real transform_matrix(4,4)
  integer vertex_material(order_max,face_max)
  real vertex_normal(3,order_max,face_max)
  real vertex_tex_uv(2,order_max,face_max)
  real vertex_range(6)
!
  reverse_faces = .false.
  reverse_normals = .false.
!
!  Initialize the graphics data.
!
  call data_init ( cor3, cor3_material, cor3_normal, cor3_tex_uv, face, face_area, &
    face_material, face_normal, face_order, face_tex_uv, iseed, line_dex, line_material, &
    material_name, material_rgba, cor3_max, face_max, line_max, material_max, &
    order_max, texture_max, normal_temp, color_num, cor3_num, face_num, group_num, &
    line_num, material_num, object_num, texture_num, object_name, texture_name, &
    texture_temp, transform_matrix, vertex_material, vertex_normal, vertex_tex_uv )
!
!
! yurkin - here the code for reading the command line was removed, since it was skipped
! anyway. It is OK since this routine is not used as an entry point.

!  Check the input file name.
!
  call infile ( filein_name, ierror, filein_type )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'COMMAND_LINE - Fatal error!'
    write ( *, * ) '  Improper input file name.'
    return
  end if
!
!  Read the input file.
!
  call data_read ( cor3, cor3_material, cor3_normal, cor3_tex_uv, &
    debug, face, face_area, face_material, face_normal, face_order, face_tex_uv, &
    filein_name, filein_type, ierror, iseed, line_dex, line_material, &
    material_name, material_rgba, cor3_max, face_max, line_max, material_max, &
    order_max, texture_max, normal_temp, color_num, cor3_num, face_num, &
    group_num, line_num, material_num, object_num, texture_num, texture_name, &
    texture_temp, vertex_material, vertex_normal, vertex_tex_uv, vertex_range )
!
!  -RN: Reverse the normal vectors at points, vertices, and faces.
!
  if ( reverse_normals ) then

    cor3_normal(1:3,1:cor3_num) = - cor3_normal(1:3,1:cor3_num)
    vertex_normal(1:3,1:order_max,1:face_num) = &
      - vertex_normal(1:3,1:order_max,1:face_num)

    face_normal(1:3,1:face_num) = - face_normal(1:3,1:face_num)
!
!  -RF: Reverse the order of the nodes defining each face.
!
  else if ( reverse_faces ) then
 
    call face_reverse_order ( cor3_normal, face, face_normal, face_order, &
      cor3_max, face_max, order_max, cor3_num, face_num, vertex_material, &
      vertex_normal, vertex_tex_uv )

  end if

  return
end
subroutine cor3_normal_set ( cor3_normal, face, face_area, &
  face_order, cor3_max, face_max, order_max, face_num, vertex_normal )
!
!*******************************************************************************
!
!! COR3_NORMAL_SET recomputes zero node normal vectors.
!
!
!  Modified:
!
!    04 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real COR3_NORMAL(3,COR3_MAX), normals at nodes.  
!
!    Input, integer FACE(ORDER_MAX,FACE_MAX), the nodes making faces.
!
!    Input, integer FACE_ORDER(FACE_MAX), the number of vertices per face.
!
!    Input, integer COR3_MAX, the maximum number of points.
!
!    Input, integer FACE_MAX, the maximum number of faces.
!
!    Input, integer ORDER_MAX, the maximum number of vertices per face.
!
!    Input, integer FACE_NUM, the number of faces.
!
!    Input, real VERTEX_NORMAL(3,ORDER_MAX,FACE_MAX), normals at vertices.  
!
  implicit none
  integer cor3_max
  integer face_max
  integer order_max
!
  real cor3_normal(3,cor3_max)
  real cor3_temp(cor3_max)
  integer face(order_max,face_max)
  real face_area(face_max)
  integer face_num
  integer face_order(face_max)
  integer icor3
  integer iface
  integer ivert
  integer j
  real norm
  real vertex_normal(3,order_max,face_max)
!
!  Determine which nodes need to have their normals recomputed.
!
  do icor3 = 1, cor3_max

    norm = 0.0
    do j = 1, 3
      norm = norm + cor3_normal(j,icor3)**2
    end do

    cor3_temp(icor3) = norm

  end do
!
!  Add up the vertex normals from all the faces to which the node belongs,
!  weighted by area.
!
  do iface = 1, face_num

    do ivert = 1, face_order(iface)

      icor3 = face(ivert,iface)

      if ( cor3_temp(icor3) == 0.0 ) then

        do j = 1, 3
          cor3_normal(j,icor3) = cor3_normal(j,icor3) &
               + face_area(iface) * vertex_normal(j,ivert,iface)
        end do

      end if

    end do

  end do
!
!  Renormalize.
!
  do icor3 = 1, cor3_max

    norm = 0.0
    do j = 1, 3
      norm = norm + cor3_normal(j,icor3)**2
    end do

    if ( norm == 0.0 ) then

      norm = 3.0
      cor3_normal(1:3,icor3) = 1.0

    end if

    norm = sqrt ( norm )

    cor3_normal(1:3,icor3) = cor3_normal(1:3,icor3) / norm

  end do

  return
end
subroutine cor3_range ( cor3_max, cor3_num, cor3, vertex_range )
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
!    Output, real VERTEX_RANGE(6), min & max values for vertex coordinate
!
  implicit none
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
  real vertex_range(6)
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

  vertex_range(1) = xmin
  vertex_range(2) = xmax
  vertex_range(3) = ymin
  vertex_range(4) = ymax
  vertex_range(5) = zmin
  vertex_range(6) = zmax

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
subroutine cross0_3d ( x0, y0, z0, x1, y1, z1, x2, y2, z2, x3, y3, z3 )
!
!*******************************************************************************
!
!! CROSS0_3D computes the cross product of (P1-P0) and (P2-P0) in 3D.
!
!
!  Discussion:
!
!    The vectors are specified with respect to a basis point P0.
!    We are computing the normal to the triangle (P0,P1,P2).
!
!  Modified:
!
!    19 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X0, Y0, Z0, X1, Y1, Z1, X2, Y2, Z2, the coordinates of
!    three points.  The basis point is (X0,Y0,Z0).
!
!    Output, real X3, Y3, Z3, the cross product (P1-P0) x (P2-P0), or
!    (X1-X0,Y1-Y0,Z1-Z0) x (X2-X0,Y2-Y0,Z2-Z0).
!
  implicit none
  real x0
  real x1
  real x2
  real x3
  real y0
  real y1
  real y2
  real y3
  real z0
  real z1
  real z2
  real z3
!
  x3 = ( y1 - y0 ) * ( z2 - z0 ) - ( z1 - z0 ) * ( y2 - y0 )
  y3 = ( z1 - z0 ) * ( x2 - x0 ) - ( x1 - x0 ) * ( z2 - z0 )
  z3 = ( x1 - x0 ) * ( y2 - y0 ) - ( y1 - y0 ) * ( x2 - x0 )

  return
end
subroutine data_check ( face_order, material_name, cor3_max, face_max, &
  line_max, material_max, order_max, texture_max, cor3_num, face_num, line_num, &
  material_num, texture_num, texture_name )
!
!*******************************************************************************
!
!! DATA_CHECK checks the input data, and enforces limits.
!
!
!  Modified:
!
!    05 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer FACE_ORDER(FACE_MAX), the number of vertices per face.
!
!    Input/output, character ( len = 100 ) MATERIAL_NAME(MATERIAL_MAX), material names.
!
!    Input, integer COR3_MAX, the maximum number of points.
!
!    Input, integer FACE_MAX, the maximum number of faces.
!
!    Input, integer LINE_MAX, the maximum number of line definition items.
!
!    Input, integer MATERIAL_MAX, the maximum number of materials.
!
!    Input, integer ORDER_MAX, the maximum number of vertices per face.
!
!    Input, integer TEXTURE_MAX, the maximum number of textures.
!
!    Input/output, integer COR3_NUM, the number of points.
!
!    Input/output, integer FACE_NUM, the number of faces.
!
!    Input/output, integer LINE_NUM, the number of line definition items.
!
!    Input/output, integer MATERIAL_NUM, the number of materials.
!
!    Input/output, integer TEXTURE_NUM, the number of textures.
!
!    Input/output, character ( len = 100 ) TEXTURE_NAME(TEXTURE_MAX), texture names.
!
  implicit none
  integer cor3_max
  integer face_max
  integer line_max
  integer material_max
  integer order_max
  integer texture_max
!
  character ( len = 4 ) char4
  integer cor3_num
  integer face_num
  integer face_order(face_max)
  integer i
  integer iface
  integer line_num
  integer material_num
  character ( len = 100 ) material_name(material_max)
  integer nfix
  character ( len = 100 ) texture_name(texture_max)
  integer texture_num
!
  if ( cor3_num > cor3_max ) then
    write ( *, * ) ' '
    write ( *, * ) 'DATA_CHECK - Warning!'
    write ( *, * ) '  The input data requires ', cor3_num, ' points.'
    write ( *, * ) '  There was only room for ', cor3_max
    cor3_num = cor3_max
  else if ( cor3_num < 0 ) then
    cor3_num = 0
  end if
 
  if ( face_num > face_max ) then
    write ( *, * ) ' '
    write ( *, * ) 'DATA_CHECK - Warning!'
    write ( *, * ) '  The input data requires ', face_num, ' faces.'
    write ( *, * ) '  There was only room for ', face_max
    face_num = face_max
  else if ( face_num <0 ) then
    face_num = 0
  end if

  if ( line_num > line_max ) then
    write ( *, * ) ' '
    write ( *, * ) 'DATA_CHECK - Warning!'
    write ( *, * ) '  The input data requires ', line_num, ' line items.'
    write ( *, * ) '  There was only room for ', line_max
    line_num = line_max
  else if ( line_num < 0 ) then 
    line_num = 0
  end if

  if ( material_num > material_max ) then
    write ( *, * ) ' '
    write ( *, * ) 'DATA_CHECK - Warning!'
    write ( *, * ) '  The input data requires ', material_num, ' materials.'
    write ( *, * ) '  There was only room for ', material_max
    material_num = material_max
  else if ( material_num < 0 ) then 
    material_num = 0
  end if

  nfix = 0

  do iface = 1, face_num

    if ( face_order(iface) < 3 ) then
      nfix = nfix + 1
      face_order(iface) = 0
    end if

  end do

  if ( nfix > 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'DATA_CHECK - Warning!'
    write ( *, * ) '  Corrected ', nfix, &
      ' faces using less than 3 vertices per face.'
  end if
  nfix = 0

  do iface = 1, face_num

    if ( face_order(iface) > order_max ) then
      nfix = nfix + 1
      face_order(iface) = order_max
    end if

  end do

  if ( nfix > 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'DATA_CHECK - Warning!'
    write ( *, * ) '  Corrected ', nfix, ' faces using more than ', &
      order_max, ' vertices per face.'
  end if

  do i = 1, material_num
    if ( material_name(i) == ' ' ) then
      call i_to_s_zero ( i, char4 )
      material_name(i) = 'Material_' // char4
    end if
  end do

  do i = 1, texture_num
    if ( texture_name(i) == ' ' ) then
      call i_to_s_zero ( i, char4 )
      texture_name(i) = 'Texture_' // char4
    end if
  end do

  return
end
subroutine data_init ( cor3, cor3_material, cor3_normal, cor3_tex_uv, face, &
  face_area, face_material, face_normal, face_order, face_tex_uv, iseed, &
  line_dex, line_material, material_name, material_rgba, cor3_max, face_max, &
  line_max, material_max, order_max, texture_max, normal_temp, color_num, cor3_num, &
  face_num, group_num, line_num, material_num, object_num, texture_num, object_name, &
  texture_name, texture_temp, transform_matrix, vertex_material, vertex_normal, &
  vertex_tex_uv )
!
!*******************************************************************************
!
!! DATA_INIT initializes internal graphics data.
!
!
!  Modified:
!
!    02 June 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real COR3(3,COR3_MAX), the coordinates of points.
!
!    Output, integer COR3_MATERIAL(COR3_MAX), the material of each node.
!
!    Output, real COR3_TEX_UV(2,COR3_MAX), UV texture coordinates for nodes.
!
!    Output, integer FACE(ORDER_MAX,FACE_MAX), the nodes making faces.
!
!    Output, integer FACE_MATERIAL(FACE_MAX), the material of each face.
!
!    Output, real FACE_NORMAL(3,FACE_MAX), the normal vector at each face.
!
!    Output, integer FACE_ORDER(FACE_MAX), the number of vertices per face.
!
!    Output, integer LINE_DEX(LINE_MAX), nodes forming a line, terminated by -1.
!
!    Output, integer LINE_MAT(LINE_MAX), material index for line.
!
!    Output, character ( len = 100 ) MATERIAL_NAME(MATERIAL_MAX), material names.
!
!    Input, integer COR3_MAX, the maximum number of points.
!
!    Input, integer FACE_MAX, the maximum number of faces.
!
!    Input, integer LINE_MAX, the maximum number of line definition items.
!
!    Input, integer MATERIAL_MAX, the maximum number of materials.
!
!    Input, integer ORDER_MAX, the maximum number of vertices per face.
!
!    Workspace, real NORMAL_TEMP(3,ORDER_MAX*FACE_MAX).
!
!    Output, integer COR3_NUM, the number of points.
!
!    Output, integer FACE_NUM, the number of faces.
!
!    Output, integer LINE_NUM, the number of line data items.
!
!    Output, integer MATERIAL_NUM, the number of materials.
!
!    Output, character ( len = 100 ) TEXTURE_NAME(TEXTURE_MAX), texture names.
!
!    Output, integer VERTEX_MAT(ORDER_MAX,FACE_MAX), vertex materials.
!
!    Output, real VERTEX_NORMAL(3,ORDER_MAX,FACE_MAX), normals at vertices.  
!
!    Output, real VERTEX_TEX_UV(2,ORDER_MAX,FACE_MAX), vertex texture coordinates.
!
  implicit none
  integer, parameter :: OFFSET = 1
!
  integer cor3_max
  integer face_max
  integer line_max
  integer material_max
  integer order_max
  integer texture_max
!
  character ( len = 4 ) char4
  integer color_num
  real cor3(3,cor3_max)
  integer cor3_material(cor3_max)
  real cor3_normal(3,cor3_max)
  integer cor3_num
  real cor3_tex_uv(2,cor3_max)
  integer face(order_max,face_max)
  real face_area(face_max)
  integer face_material(face_max)
  real face_normal(3,face_max)
  integer face_num
  integer face_order(face_max)
  real face_tex_uv(2,face_max)
  integer group_num
  integer i
  integer iseed
  integer j
  integer line_dex(line_max)
  integer line_material(line_max)
  integer line_num
  character ( len = 100 ) material_name(material_max)
  integer material_num
  real material_rgba(4,material_max)
  real normal_temp(3,order_max*face_max)
  character ( len = 100 ) object_name
  integer object_num
  character ( len = 100 ) texture_name(texture_max)
  integer texture_num
  real texture_temp(2,order_max*face_max)
  real transform_matrix(4,4)
  integer vertex_material(order_max,face_max)
  real vertex_normal(3,order_max,face_max)
  real vertex_tex_uv(2,order_max,face_max)
!
  cor3(1:3,1:cor3_max) = 0.0
  cor3_material(1:cor3_max) = -1
  cor3_normal(1:3,1:cor3_max) = 0.0
  cor3_tex_uv(1:2,1:cor3_max) = 0.0
  face(1:order_max,1:face_max) = -1 + OFFSET
  face_area(1:face_max) = 0.0
  face_material(1:face_max) = 0 + OFFSET
  face_normal(1:3,1:face_max) = 0.0
  face_order(1:face_max) = 0
  face_tex_uv(1:2,1:face_max) = 0.0

  call get_seed ( iseed )

  line_dex(1:line_max) = -1 + OFFSET
  line_material(1:line_max) = 0 + OFFSET

  do j = 1, material_max
    material_name(j) = ' '
  end do

  material_rgba(1:4,1:material_max) = 0.0

  do j = 1, order_max*face_max
    normal_temp(1:3,j) = 0.0
  end do

  color_num = 0
  cor3_num = 0
  face_num = 0
  group_num = 0
  line_num = 0
  material_num = 0
  object_num = 0
  texture_num = 0

  object_name = 'IVCON'

  do i = 1, texture_max
    call i_to_s_zero ( i, char4 )
    texture_name(i) = 'Texture_' // char4
  end do

  do j = 1, order_max*face_max
    texture_temp(1:2,j) = 0.0
  end do

  call tmat_init ( transform_matrix )

  vertex_material(1:order_max,1:face_max) = 0 + OFFSET
  vertex_normal(1:3,1:order_max,1:face_max) = 0.0
  vertex_tex_uv(1:2,1:order_max,1:face_max) = 0.0

  write ( *, * ) ' '
  write ( *, * ) 'DATA_INIT: Initialized graphics data.'

  return
end
subroutine data_read ( cor3, cor3_material, cor3_normal, cor3_tex_uv, &
  debug, face, face_area, face_material, face_normal, face_order, face_tex_uv, &
  filein_name, filein_type, ierror, iseed, line_dex, line_material, material_name, &
  material_rgba, cor3_max, face_max, line_max, material_max, order_max, texture_max, &
  normal_temp, color_num, cor3_num, face_num, group_num, line_num, material_num, &
  object_num, texture_num, texture_name, texture_temp, vertex_material, vertex_normal, &
  vertex_tex_uv, vertex_range )
!
!*******************************************************************************
!
!! DATA_READ reads a file into internal graphics data.
!
!
!  Modified:
!
!    26 September 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real COR3(3,COR3_MAX), the coordinates of points.
!
!    Output, integer COR3_MATERIAL(COR3_MAX), the material of each node.
!
!    Output, real COR3_TEX_UV(2,COR3_MAX), UV texture coordinates for nodes.
!
!    Output, integer FACE(ORDER_MAX,FACE_MAX), the nodes making faces.
!
!    Output, integer FACE_MATERIAL(FACE_MAX), the material of each face.
!
!    Output, real FACE_NORMAL(3,FACE_MAX), the normal vector at each face.
!
!    Output, integer FACE_ORDER(FACE_MAX), the number of vertices per face.
!
!    Input/output, character ( len = 100 ) FILEIN_NAME, the name of the input file.
!
!    Output, integer IERROR, an error flag.
!
!    Output, integer LINE_DEX(LINE_MAX), nodes forming a line, terminated by -1.
!
!    Output, integer LINE_MAT(LINE_MAX), material index for line.
!
!    Output, character ( len = 100 ) MATERIAL_NAME(MATERIAL_MAX), material names.
!
!    Input, integer COR3_MAX, the maximum number of points.
!
!    Input, integer FACE_MAX, the maximum number of faces.
!
!    Input, integer LINE_MAX, the maximum number of line definition items.
!
!    Input, integer MATERIAL_MAX, the maximum number of materials.
!
!    Input, integer ORDER_MAX, the maximum number of vertices per face.
!
!    Workspace, real NORMAL_TEMP(3,ORDER_MAX*FACE_MAX).
!
!    Output, integer COR3_NUM, the number of points.
!
!    Output, integer FACE_NUM, the number of faces.
!
!    Output, integer LINE_NUM, the number of line definition items.
!
!    Output, integer MATERIAL_NUM, the number of materials.
!
!    Output, character ( len = 100 ) TEXTURE_NAME(TEXTURE_MAX), texture names.
!
!    Output, integer VERTEX_MAT(ORDER_MAX,FACE_MAX), vertex materials.
!
!    Output, real VERTEX_NORMAL(3,ORDER_MAX,FACE_MAX), normals at vertices.  
!
!    Output, real VERTEX_TEX_UV(2,ORDER_MAX,FACE_MAX), vertex texture coordinates.
!
!    Output, real VERTEX_RANGE(6), min & max values for vertex coordinate
!
  implicit none
  integer, parameter :: iunit = 1
  integer, parameter :: OFFSET = 1
!
  integer cor3_max
  integer face_max
  integer line_max
  integer material_max
  integer order_max
  integer texture_max
!
  integer bad_num
  integer color_num
  real cor3(3,cor3_max)
  integer cor3_material(cor3_max)
  real cor3_normal(3,cor3_max)
  integer cor3_num
  real cor3_tex_uv(2,cor3_max)
  logical debug
  integer dup_num
  integer face(order_max,face_max)
  real face_area(face_max)
  integer face_material(face_max)
  real face_normal(3,face_max)
  integer face_num
  integer face_order(face_max)
  real face_tex_uv(2,face_max)
  character ( len = 100 ) filein_name
  character ( len = 10 ) filein_type
  integer group_num
  integer i
  integer ierror
  integer iline
  integer imax
  integer ios
  integer iseed
  integer j
  integer line_dex(line_max)
  integer line_material(line_max)
  integer line_num
  character ( len = 100 ) material_name(material_max)
  integer material_num
  real material_rgba(4,material_max)
  integer ncol_oogl
  integer nrow_oogl
  real normal_temp(3,order_max*face_max)
  integer ntemp
  integer object_num
  logical s_eqi
  integer text_num
  character ( len = 100 ) texture_name(texture_max)
  integer texture_num
  real texture_temp(2,order_max*face_max)
  integer vertex_material(order_max,face_max)
  real vertex_normal(3,order_max,face_max)
  real vertex_tex_uv(2,order_max,face_max)
  real vertex_range(6)
!
  ierror = 0
  bad_num = 0
  dup_num = 0
  text_num = 0
!
!  Open the file.
!
  open ( unit = iunit, file = filein_name, status = 'old', iostat = ios )
 
  if ( ios /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'DATA_READ - Error!'
    write ( *, * ) '  The input file:'
    write ( *, '(a)' ) trim ( filein_name )
    write ( *, * ) '  could not be opened!'
    ierror = ios
    return
  end if
!
!  Read the information in the file.
!
  if ( s_eqi ( filein_type, 'DXF' ) ) then
 
    call dxf_read ( cor3, face, face_material, face_order, ierror, iunit, &
      line_dex, line_material, cor3_max, face_max, line_max, order_max, bad_num, &
      cor3_num, dup_num, face_num, line_num, material_num, text_num )

  else if ( s_eqi ( filein_type, 'OBJ' ) ) then

    call obj_read ( cor3, face, face_material, face_order, ierror, iseed, iunit, &
      line_dex, line_material, material_name, material_rgba, cor3_max, face_max, &
      line_max, material_max, order_max, normal_temp, bad_num, cor3_num, face_num, &
      group_num, line_num, material_num, object_num, text_num, vertex_material, &
      vertex_normal )

  else if ( s_eqi ( filein_type, 'STL' ) .or. &
            s_eqi ( filein_type, 'STLA' ) ) then

    call stla_read ( cor3, face, face_material, face_normal, face_order, ierror, &
      iunit, cor3_max, face_max, order_max, bad_num, cor3_num, dup_num, &
      face_num, material_num, object_num, text_num, vertex_material, vertex_normal )

  else if ( s_eqi ( filein_type, 'WRL' ) ) then

    call vrml_read ( cor3, cor3_material, face, face_material, face_order, &
      ierror, iunit, material_name, material_rgba, cor3_max, face_max, &
      material_max, order_max, bad_num, cor3_num, face_num, line_num, &
      material_num, text_num, vertex_material )

  else 

    write ( *, * ) ' '
    write ( *, * ) 'DATA_READ - Fatal error!'
    write ( *, * ) '  Unrecognized input file type:'
    write ( *, '(a)' ) trim ( filein_type )
    ierror = 1
    return

  end if
 
  close ( unit = iunit )
!
!  Make a report on what we saw in the file.
!
  ntemp = min ( face_num, face_max )
  call ivec_max ( ntemp, face_order, imax )

  call data_report ( imax, bad_num, color_num, cor3_num, dup_num, face_num, &
    group_num, line_num, material_num, object_num, text_num )
!
!  Warn about any errors that occurred during reading.
!
  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'DATA_READ - Warning!'
    write ( *, * ) '  An error occurred reading the input file.'
    return
  end if
!
!  Check the data.
!
!  You MUST wait until after this check before doing other computations,
!  since COR3_NUM, and the other variables could be much larger than
!  the legal maximums, until corrected by this routine.
!
  call data_check ( face_order, material_name, cor3_max, face_max, line_max, &
    material_max, order_max, texture_max, cor3_num, face_num, line_num, material_num, &
    texture_num, texture_name )
!
!  MATERIALS FIXUPS:
!
!  If there are no materials, define one.
!
  if ( material_num <= 0 ) then

    material_num = 1

    material_name(material_num) = 'Default_Material'

    material_rgba(1,material_num) = 0.7
    material_rgba(2,material_num) = 0.7
    material_rgba(3,material_num) = 0.7
    material_rgba(4,material_num) = 1.0

  end if
!
!  If a point hasn't been assigned a material, set it to material 1.
!
  do i = 1, cor3_num
    if ( cor3_material(i) < 1 .or. cor3_material(i) > material_num ) then
      cor3_material(i) = 1
    end if
  end do
!
!  If a vertex hasn't been assigned a material, set it to material 1.
!
  do i = 1, face_num
    do j = 1, face_order(i)
      if ( vertex_material(j,i) < 1 .or. vertex_material(j,i) > material_num ) then
        vertex_material(j,i) = 1
      end if
    end do
  end do
!
!  If a face hasn't been assigned a material, set it to material 1.
!
  do i = 1, face_num
    if ( face_material(i) < 1 .or. face_material(i) > material_num ) then
      face_material(i) = 1
    end if
  end do
!
!  If a line item hasn't been assigned a material, set it to material 1.
!
  do iline = 1, line_num
    if ( line_dex(iline) == -1 + OFFSET ) then
      line_material(iline) = -1 + OFFSET
    else if ( line_material(iline) < 1 .or. line_material(iline) > material_num ) then
      line_material(iline) = material_num
    end if
  end do
!
!  NULL EDGE DELETION.
!
  call edge_null_delete ( cor3, face, face_order, cor3_max, face_max, &
    order_max, face_num, vertex_normal )
!
!  COMPUTE FACE AREAS.
!
  call face_area_set ( cor3, face, face_area, face_order, cor3_max, face_max, &
    order_max, face_num )
!
!  NULL FACE DELETION.
!
  call face_null_delete ( face, face_area, face_material, face_order, face_max, &
    order_max, face_num, vertex_material, vertex_normal )
!
!  NORMAL VECTOR FIXUPS:
!
!  Recompute zero vertex normals from vertex positions.
!
  call vertex_normal_set ( cor3, face, face_order, cor3_max, face_max, &
    order_max, face_num, vertex_normal )
!
!  Recompute zero node normals by averaging vertex normals.
!
  call cor3_normal_set ( cor3_normal, face, face_area, face_order, &
    cor3_max, face_max, order_max, face_num, vertex_normal )
!
!  Recompute zero face normals by averaging vertex normals.
!
  call face_normal_ave ( face_normal, face_order, face_max, order_max, &
    face_num, vertex_normal )
!
!  Report the range of the nodal coordinates.
!
  call cor3_range ( cor3_max, cor3_num, cor3, vertex_range )

  return
end
subroutine data_report ( imax, bad_num, color_num, cor3_num, dup_num, &
  face_num, group_num, line_num, material_num, object_num, text_num )
!
!*******************************************************************************
!
!! DATA_REPORT gives a summary of the contents of the data file.
!
!
!  Modified:
!
!    09 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer COR3_NUM, the number of points.
!
!    Input, integer FACE_NUM, the number of faces.
!
!    Input, integer LINE_NUM, the number of line definition items.
!
!    Input, integer MATERIAL_NUM, the number of materials.
!
  implicit none
  integer imax
  integer bad_num
  integer color_num
  integer cor3_num
  integer dup_num
  integer face_num
  integer group_num
  integer line_num
  integer material_num
  integer object_num
  integer text_num
!
  write ( *, * ) ' '
  write ( *, * ) 'DATA_REPORT - The input file contains:'
  write ( *, * ) ' '
  write ( *, * ) '  Bad data items             ', bad_num
  write ( *, * ) '  Text lines                 ', text_num
  write ( *, * ) '  Colors                     ', color_num
  write ( *, * ) '  Duplicate points           ', dup_num
  write ( *, * ) '  Faces                      ', face_num
  write ( *, * ) '  Groups                     ', group_num
  write ( *, * ) '  Vertices per face, maximum ', imax
  write ( *, * ) '  Line items                 ', line_num
  write ( *, * ) '  Materials                  ', material_num
  write ( *, * ) '  Points                     ', cor3_num
  write ( *, * ) '  Objects                    ', object_num

  return
end
function degrees_to_radians ( angle )
!
!*******************************************************************************
!
!! DEGREES_TO_RADIANS converts an angle from degrees to radians.
!
!
!  Modified:
!
!    10 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ANGLE, an angle in degrees.
!
!    Output, real DEGREES_TO_RADIANS, the equivalent angle
!    in radians.
!
  implicit none
  real, parameter :: pi = 3.14159265358979323846264338327950288419716939937510

  real angle
  real degrees_to_radians

  degrees_to_radians = ( angle / 180.0 ) * pi

  return
end
subroutine digit_to_c ( digit, c )
!
!*******************************************************************************
!
!! DIGIT_TO_C returns the character representation of a decimal digit.
!
!
!  Example:
!
!    DIGIT   C
!    -----  ---
!      0    '0'
!      1    '1'
!    ...    ...
!      9    '9'
!     17    '*'
!
!  Modified:
!
!    04 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer DIGIT, the digit value between 0 and 9.
!
!    Output, character C, the corresponding character, or '*' if DIGIT
!    was illegal.
!
  implicit none
  character c
  integer digit
!
  if ( 0 <= digit .and. digit <= 9 ) then

    c = char ( digit + 48 )

  else

    c = '*'

  end if

  return
end
subroutine dxf_read ( cor3, face, face_material, face_order, ierror, iunit, &
  line_dex, line_material, cor3_max, face_max, line_max, order_max, bad_num, &
  cor3_num, dup_num, face_num, line_num, material_num, text_num )
!
!*******************************************************************************
!
!! DXF_READ reads graphics information from an AutoCAD DXF file.
!
!
!  Discussion:
!
!    It is intended that the information read from the file can
!    either start a whole new graphics object, or simply be added
!    to a current graphics object via the '<<' command.  
!
!    This is controlled by whether the input values have been zeroed 
!    out or not.  This routine simply tacks on the information it 
!    finds to the current graphics object.
!
!  Example:
!
!      0
!    SECTION
!      2
!    HEADER
!    999
!    diamond.dxf created by IVREAD.
!    999
!    Original data in diamond.obj.
!      0
!    ENDSEC
!      0
!    SECTION
!      2
!    TABLES
!      0
!    ENDSEC
!      0
!    SECTION
!      2
!    BLOCKS
!      0
!    ENDSEC
!      0
!    SECTION
!      2
!    ENTITIES
!      0
!    LINE
!      8
!    0
!     10
!      0.00  (X coordinate of beginning of line.)
!     20
!      0.00  (Y coordinate of beginning of line.)
!     30
!      0.00  (Z coordinate of beginning of line.)
!     11
!      1.32  (X coordinate of end of line.)
!     21
!      1.73  (Y coordinate of end of line.)
!     31
!      2.25  (Z coordinate of end of line.)
!      0
!    3DFACE
!      8
!     Cube
!    10
!    -0.50  (X coordinate of vertex 1)
!    20
!     0.50  (Y coordinate of vertex 1)   
!    30
!      1.0  (Z coordinate of vertex 1)  
!    11
!     0.50  (X coordinate of vertex 2)  
!    21
!     0.50  (Y coordinate of vertex 2)
!    31
!      1.0  (Z coordinate of vertex 2)
!    12
!     0.50  (X coordinate of vertex 3) 
!    22
!     0.50  (Y coordinate of vertex 3)
!    32
!     0.00  (Z coordinate of vertex 3)
!     0
!    ENDSEC
!      0
!    EOF
!
!  Modified:
!
!    02 June 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real COR3(3,COR3_MAX), the coordinates of points.
!
!    Input/output, integer FACE(ORDER_MAX,FACE_MAX), the nodes making faces.
!
!    Input/output, integer FACE_MATERIAL(FACE_MAX), the material of each face.
!
!    Input/output, integer FACE_ORDER(FACE_MAX), the number of vertices per face.
!
!    Output, integer IERROR, an error flag.
!
!    Input, integer IUNIT, the FORTRAN unit from which data is read.
!
!    Input/output, integer LINE_DEX(LINE_MAX), nodes forming a line, terminated by -1.
!
!    Input/output, integer LINE_MAT(LINE_MAX), material index for line.
!
!    Input, integer COR3_MAX, the maximum number of points.
!
!    Input, integer FACE_MAX, the maximum number of faces.
!
!    Input, integer LINE_MAX, the maximum number of line items.
!
!    Input, integer ORDER_MAX, the maximum number of vertices per face.
!
!    Output, integer BAD_NUM, the number of "bad" lines of input text.
!
!    Input/output, integer COR3_NUM, the number of points.
!
!    Output, integer DUP_NUM, the number of duplicate nodes discovered.
!
!    Input/output, integer FACE_NUM, the number of faces.
!
!    Input/output, integer LINE_NUM, the number of line definition items.
!
!    Output, integer TEXT_NUM, the number of lines of input text.
!
  implicit none
  integer, parameter :: OFFSET = 1

  integer cor3_max
  integer face_max
  integer line_max
  integer order_max
!
  integer bad_num
  real cor3(3,cor3_max)
  integer cor3_num
  real cvec(3)
  integer dup_num
  integer face(order_max,face_max)
  integer face_material(face_max)
  integer face_num
  integer face_order(face_max)
  integer icode
  integer icor3
  integer ierror
  integer ios
  integer iunit
  integer ival
  integer ivert
  integer ixyz
  integer lchar
  integer line_dex(line_max)
  integer line_material(line_max)
  integer line_num
  integer material_num
  character ( len = 256 ) mode
  real rval
  logical s_eqi
  integer text_num
  character ( len = 256 ) text1
  character ( len = 256 ) text2
!
  ierror = 0
  mode = ' '
!
!  Read the next pair of lines.  
!  TEXT1 is a numeric tag, TEXT2 contains data.
!
  do
 
    read ( iunit, '(a)', iostat = ios ) text1

    if ( ios /= 0 ) then
      exit
    end if

    text_num = text_num + 1
 
    call s_to_i ( text1, icode, ierror, lchar )
 
    if ( ierror /= 0 ) then
      ierror = 1
      bad_num = bad_num + 1
      write ( *, * ) ' '
      write ( *, * ) 'DXF_READ - Fatal error!'
      write ( *, * ) '  Could not interpret line:'
      write ( *, '(a)' ) trim ( text1 )
      return
    end if
!
!  Read the second item, which might be a label or numeric value.
!
    read ( iunit, '(a)', iostat = ios ) text2

    if ( ios /= 0 ) then

      if ( text1 /= ' ' ) then
        ierror = 3
        write ( *, * ) ' '
        write ( *, * ) 'DXF_READ - Warning!'
        write ( *, * ) '  The last code was not followed by data:'
        write ( *, '(a)' ) trim ( text1 )
      end if

      exit

    end if

    text_num = text_num + 1
!
!  Codes 0 through 9 are followed by a string.
!
!  All we want to do is know when we can expect to be reading
!  LINE data and when we can expect to be reading FACE data.
!
    if ( icode >= 0 .and. icode <= 9 ) then

      if ( s_eqi ( text2(1:6), '3DFACE' ) ) then
        mode = '3DFACE'
        ivert = 0
      else if ( s_eqi ( text2(1:4), 'LINE' ) ) then
        mode = 'LINE'
      end if
!
!  Codes 10 through 59 are followed by a real value.
!
!    10, 11, 12, ... are followed by a line of X data;
!    20, 21, 22, ... are followed by a line of Y data;
!    30, 31, 33, ... are followed by a line of Z data.
!
    else if ( icode >= 10 .and. icode <= 59 ) then
 
      call s_to_r ( text2, rval, ierror, lchar )
 
      if ( ierror /= 0 ) then

        ierror = 2
        rval = 0.0
        bad_num = bad_num + 1

        if ( bad_num == 1 ) then
          write ( *, * ) ' '
          write ( *, * ) 'DXF_READ - Fatal error!'
          write ( *, * ) '  Could not interpret line:'
          write ( *, '(a)' ) trim ( text2 )
        end if

      end if
 
      if ( mode == 'LINE' ) then

        if ( icode == 10 .and. line_num > 0 ) then

          line_num = line_num + 1
          if ( line_num <= line_max ) then
            line_dex(line_num) = -1 + OFFSET
            line_material(line_num) = -1 + OFFSET
          end if

        end if

      else if ( mode == '3DFACE' ) then

        if ( icode == 10 ) then

          face_num = face_num + 1
          face_order(face_num) = 0
          face_material(face_num) = material_num

        end if

      end if

      if ( icode >= 10 .and. icode <= 19 ) then
        ixyz = 1
      else if ( icode >= 20 .and. icode <= 29 ) then
        ixyz = 2
      else if ( icode >= 30 .and. icode <= 39 ) then
        ixyz = 3
      end if

      cvec(ixyz) = rval
!
!  Once the entire (X,Y,Z) triple has been read, check to see if the
!  values in CVEC already exist in COR3.  If so, save space by using 
!  the index of a previous copy.
!
!  Otherwise, add CVEC to COR3, and increment COR3_NUM.
!
      if ( ixyz == 3 ) then

        if ( cor3_num <= 1000 ) then
          call rcol_find ( 3, cor3_num, cor3, cvec, icor3 )
        else
          icor3 = 0
        end if
 
        if ( icor3 == 0 ) then

          cor3_num = cor3_num + 1
          icor3 = cor3_num

          if ( cor3_num <= cor3_max ) then
            cor3(1:3,cor3_num) = cvec(1:3)
          end if

        else

          dup_num = dup_num + 1

        end if

        if ( mode == 'LINE' ) then

          line_num = line_num + 1

          if ( line_num <= line_max ) then
            line_dex(line_num) = icor3 - 1 + OFFSET
            line_material(line_num) = material_num
          end if

        else if ( mode == '3DFACE' ) then

          ivert = ivert + 1
          face(ivert,face_num) = icor3
          face_order(face_num) = face_order(face_num) + 1

        end if

      end if
!
!  Codes 60 through 79 are followed by an integer.
!
    else if ( icode >= 60 .and. icode <= 79 ) then

      call s_to_i ( text2, ival, ierror, lchar )
!
!  Codes 140 through 147 are followed by a real.
!
    else if ( icode >= 140 .and. icode <= 147 ) then

      call s_to_r ( text2, rval, ierror, lchar )
!
!  Codes 170 through 175 are followed by an integer.
!
    else if ( icode >= 170 .and. icode <= 175 ) then

      call s_to_i ( text2, ival, ierror, lchar )
!
!  Codes 210 through 239 are followed by a real.
!
    else if ( icode >= 210 .and. icode <= 239 ) then

      call s_to_r ( text2, rval, ierror, lchar )
!
!  Code 999 is followed by a (comment) string.
!
    else if ( icode == 999 ) then
!
!  Codes 1000 through 1009 are followed by a string.
!
    else if ( icode >= 1000 .and. icode <= 1009 ) then
!
!  Codes 1010 through 1059 are followed by a real.
!
    else if ( icode >= 1010 .and. icode <= 1059 ) then

      call s_to_r ( text2, rval, ierror, lchar )
!
!  Codes 1060 through 1079 are followed by an integer.
!
    else if ( icode >= 1060 .and. icode <= 1079 ) then

      call s_to_i ( text2, ival, ierror, lchar )
!
!  Unrecognized code.
!
    else

      bad_num = bad_num + 1

    end if
 
  end do
!
!  END OF INPUT.
!
!  Slap a trailing "-1" on the end of the line indices.
!
  if ( line_num > 0 ) then
    line_num = line_num + 1
    if ( line_num <= line_max ) then
      line_dex(line_num) = -1 + OFFSET
      line_material(line_num) = -1 + OFFSET
    end if
  end if
 
  return
end
subroutine edge_null_delete ( cor3, face, face_order, cor3_max, face_max, &
  order_max, face_num, vertex_normal )
!
!*******************************************************************************
!
!! EDGE_NULL_DELETE deletes face edges with zero length.
!
!
!  Modified:
!
!    12 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real COR3(3,COR3_MAX), the coordinates of points.
!
!    Input, integer FACE(ORDER_MAX,FACE_MAX), the nodes making faces.
!
!    Input, integer FACE_ORDER(FACE_MAX), the number of vertices per face.
!
!    Input, integer COR3_MAX, the maximum number of points.
!
!    Input, integer FACE_MAX, the maximum number of faces.
!
!    Input, integer ORDER_MAX, the maximum number of vertices per face.
!
!    Input, integer FACE_NUM, the number of faces.
!
  implicit none
  integer cor3_max
  integer face_max
  integer order_max
!
  real cor3(3,cor3_max)
  real distsq
  integer edge_num
  integer edge_num_del
  integer face(order_max,face_max)
  integer face_num
  integer face_order(face_max)
  integer face_order2
  integer face2(100)
  integer iface
  integer inode
  integer ivert
  integer jnode
  integer jvert
  real vertex_normal(3,order_max,face_max)
  real vertex_normal2(3,100)
!
  edge_num = 0
  edge_num_del = 0
!
!  Consider each face.
!
  do iface = 1, face_num
!
!  Consider each pair of consecutive vertices.
!
    face_order2 = 0

    do ivert = 1, face_order(iface)

      edge_num = edge_num + 1

      jvert = ivert + 1
      if ( jvert > face_order(iface) ) then
        jvert = 1
      end if

      inode = face(ivert,iface)
      jnode = face(jvert,iface)

      distsq =   ( cor3(1,inode) - cor3(1,jnode) )**2 &
               + ( cor3(2,inode) - cor3(2,jnode) )**2 &
               + ( cor3(3,inode) - cor3(3,jnode) )**2

      if ( distsq /= 0.0 ) then
        face_order2 = face_order2 + 1
        face2(face_order2) = face(ivert,iface)
        vertex_normal2(1,face_order2) = vertex_normal(1,ivert,iface)
        vertex_normal2(2,face_order2) = vertex_normal(2,ivert,iface)
        vertex_normal2(3,face_order2) = vertex_normal(3,ivert,iface)
      else
        edge_num_del = edge_num_del + 1
      end if

    end do

    face_order(iface) = face_order2
    do ivert = 1, face_order(iface)
      face(ivert,iface) = face2(ivert)
      vertex_normal(1:3,ivert,iface) = vertex_normal2(1:3,ivert)
    end do

  end do

  write ( *, * ) ' '
  write ( *, * ) 'EDGE_NULL_DELETE:'
  write ( *, * ) '  There are a total of ', edge_num, ' edges.'
  write ( *, * ) '  Of these, ', edge_num_del, &
    ' were of zero length, and deleted.'

  return
end
subroutine face_area_set ( cor3, face, face_area, face_order, cor3_max, &
  face_max, order_max, face_num )
!
!*******************************************************************************
!
!! FACE_AREA_SET computes the area of the faces.
!
!
!  Formula:
!
!    The area is the sum of the areas of the triangles formed by
!    node N with consecutive pairs of nodes.
!
!  Reference:
!
!    Adrian Bowyer and John Woodwark,
!    A Programmer's Geometry,
!    Butterworths, 1983.
!
!  Modified:
!
!    13 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real COR3(3,COR3_MAX), the coordinates of points.
!
!    Input, integer FACE(ORDER_MAX,FACE_MAX), the nodes making faces.
!
!    Output, real FACE_AREA(FACE_MAX), the area of each face.
!
!    Input, integer FACE_ORDER(FACE_MAX), the number of vertices per face.
!
!    Input, integer COR3_MAX, the maximum number of points.
!
!    Input, integer FACE_MAX, the maximum number of faces.
!
!    Input, integer ORDER_MAX, the maximum number of vertices allowed per face.
!
!    Input, integer FACE_NUM, the number of faces.
!
  implicit none
  integer cor3_max
  integer face_max
  integer order_max
!
  real alpha
  real area_max
  real area_min
  real area_tri
  real base
  real cor3(3,cor3_max)
  real dot
  integer face(order_max,face_max)
  real face_area(face_max)
  integer face_num
  integer face_num_del
  integer face_order(face_max)
  real height
  integer i
  integer i1
  integer i2
  integer i3
  integer iface
  real tol
  real x1
  real x2
  real x3
  real y1
  real y2
  real y3
  real z1
  real z2
  real z3
!
  do iface = 1, face_num

    face_area(iface) = 0.0

    do i = 1, face_order(iface) - 2

      i1 = face(i,iface)
      i2 = face(i+1,iface)
      i3 = face(i+2,iface)

      x1 = cor3(1,i1)
      y1 = cor3(2,i1)
      z1 = cor3(3,i1)

      x2 = cor3(1,i2)
      y2 = cor3(2,i2)
      z2 = cor3(3,i2)

      x3 = cor3(1,i3)
      y3 = cor3(2,i3)
      z3 = cor3(3,i3)
!
!  Find the projection of (P3-P1) onto (P2-P1).
!
      dot = ( x2 - x1 ) * ( x3 - x1 ) + &
            ( y2 - y1 ) * ( y3 - y1 ) + &
            ( z2 - z1 ) * ( z3 - z1 )

      base = sqrt ( ( x2 - x1 )**2 + ( y2 - y1 )**2 + ( z2 - z1 )**2 )
!
!  The height of the triangle is the length of (P3-P1) after its
!  projection onto (P2-P1) has been subtracted.
!
      if ( base == 0.0 ) then

        height = 0.0

      else

        alpha = dot / base**2
  
        height = sqrt ( ( x3 - x1 - alpha * ( x2 - x1 ) )**2 + &
                        ( y3 - y1 - alpha * ( y2 - y1 ) )**2 + &
                        ( z3 - z1 - alpha * ( z2 - z1 ) )**2 )

      end if

      area_tri = 0.5 * base * height

      face_area(iface) = face_area(iface) + area_tri

    end do

  end do

  area_min = face_area(1)
  area_max = face_area(1)
  do iface = 2, face_num
    area_min = min ( area_min, face_area(iface) )
    area_max = max ( area_max, face_area(iface) )
  end do

  write ( *, * ) ' '
  write ( *, * ) 'FACE_AREA_SET:'
  write ( *, * ) '  Minimum face area is ', area_min
  write ( *, * ) '  Maximum face area is ', area_max

  tol = area_max / 10000.0
  if ( area_min < tol ) then

    face_num_del = 0

    do iface = 1, face_num
      if ( face_area(iface) < tol ) then
        face_order(iface) = 0
        face_num_del = face_num_del + 1
      end if
    end do

    write ( *, * ) '  Marked ', face_num_del, ' tiny faces for deletion.'

  end if

  return
end
subroutine face_normal_ave ( face_normal, face_order, face_max, order_max, &
  face_num, vertex_normal )
!
!*******************************************************************************
!
!! FACE_NORMAL_AVE sets face normals as average of face vertex normals.
!
!
!  Modified:
!
!    09 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real FACE_NORMAL(3,FACE_MAX), the normal vector at each face.
!
!    Input, integer FACE_ORDER(FACE_MAX), the number of vertices per face.
!
!    Input, integer FACE_MAX, the maximum number of faces.
!
!    Input, integer ORDER_MAX, the maximum number of vertices per face.
!
!    Input, integer FACE_NUM, the number of faces.
!
!    Input, real VERTEX_NORMAL(3,ORDER_MAX,FACE_MAX), normals at vertices.  
!
  implicit none
  integer face_max
  integer order_max
!
  real face_normal(3,face_max)
  integer face_num
  integer face_order(face_max)
  integer i
  integer iface
  integer ivert
  integer nfix
  real norm
  real vertex_normal(3,order_max,face_max)
!
  if ( face_num <= 0 ) then
    return
  end if

  nfix = 0

  do iface = 1, face_num

    norm = 0.0
    do i = 1, 3
      norm = norm + face_normal(i,iface)**2
    end do
    norm = sqrt ( norm )

    if ( norm == 0.0 ) then

      nfix = nfix + 1

      face_normal(1:3,iface) = 0.0

      do ivert = 1, face_order(iface)
        face_normal(1:3,iface) = face_normal(1:3,iface) + &
          vertex_normal(1:3,ivert,iface)
      end do

      norm = 0.0
      do i = 1, 3
        norm = norm + face_normal(i,iface)**2
      end do
      norm = sqrt ( norm )

      if ( norm == 0.0 ) then
        face_normal(1:3,iface) = 1.0 / sqrt ( 3.0 )
      else
        face_normal(1:3,iface) = face_normal(1:3,iface) / norm
      end if

    end if

  end do

  if ( nfix > 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'FACE_NORMAL_AVE: Recomputed ', nfix, &
      ' face normals by averaging face vertex normals.'
  end if

  return
end
subroutine face_null_delete ( face, face_area, face_material, face_order, &
  face_max, order_max, face_num, vertex_material, vertex_normal )
!
!*******************************************************************************
!
!! FACE_NULL_DELETE deletes faces of order less than 3.
!
!
!  Modified:
!
!    23 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer FACE(ORDER_MAX,FACE_MAX), the nodes making faces.
!
!    Input, integer FACE_MATERIAL(FACE_MAX), the material of each face.
!
!    Input, integer FACE_ORDER(FACE_MAX), the number of vertices per face.
!
!    Input, integer COR3_MAX, the maximum number of points.
!
!    Input, integer FACE_MAX, the maximum number of faces.
!
!    Input, integer ORDER_MAX, the maximum number of vertices per face.
!
!    Input, integer FACE_NUM, the number of faces.
!
!    Input/output, real VERTEX_NORMAL(3,ORDER_MAX,FACE_MAX), normal vectors
!    at vertices.  
!
  implicit none
  integer face_max
  integer order_max
!
  integer face(order_max,face_max)
  real face_area(face_max)
  integer face_material(face_max)
  integer face_num
  integer face_num2
  integer face_order(face_max)
  integer iface
  integer vertex_material(order_max,face_max)
  real vertex_normal(3,order_max,face_max)
!
!  Drop faces of order 0, 1 or 2.
!
  face_num2 = 0

  do iface = 1, face_num

    if ( face_order(iface) >= 3 ) then

      face_num2 = face_num2 + 1

      if ( face_num2 /= iface ) then

        face_area(face_num2) = face_area(iface)
        face_material(face_num2) = face_material(iface)
        face_order(face_num2) = face_order(iface)
        face(1:order_max,face_num2) = face(1:order_max,iface)
        vertex_material(1:order_max,face_num2) = vertex_material(1:order_max,iface)
        vertex_normal(1:3,1:order_max,face_num2) = &
          vertex_normal(1:3,1:order_max,iface)

      end if

    end if

  end do

  write ( *, * ) ' '
  write ( *, * ) 'FACE_NULL_DELETE'
  write ( *, * ) '  There are a total of ', face_num, ' faces.'
  write ( *, * ) '  Of these, ', face_num2, ' passed the order test.'

  face_num = face_num2 

  return
end
subroutine face_reverse_order ( cor3_normal, face, face_normal, face_order, &
  cor3_max, face_max, order_max, cor3_num, face_num, vertex_material, &
  vertex_normal, vertex_tex_uv )
!
!*******************************************************************************
!
!! FACE_REVERSE_ORDER reverses the order of the nodes in each face.
!
!
!  Discussion:
!
!    Reversing the order of the nodes requires that the normal vectors
!    be reversed as well, so this routine will automatically reverse
!    the normals associated with nodes, vertices and faces.
!
!  Modified:
!
!    25 June 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real COR3_NORMAL(3,COR3_MAX), normals at nodes.
!
!    Input/output, integer FACE(ORDER_MAX,FACE_MAX), the nodes making faces.
!
!    Input/output, real FACE_NORMAL(3,FACE_MAX), normals at faces.
!
!    Input, integer FACE_ORDER(FACE_MAX), the number of vertices per face.
!
!    Input, integer COR3_MAX, the maximum number of nodes.
!
!    Input, integer FACE_MAX, the maximum number of faces.
!
!    Input, integer ORDER_MAX, the maximum number of vertices per face.
!
!    Input, integer COR3_NUM, the number of nodes.
!
!    Input, integer FACE_NUM, the number of faces.
!
!    Input/output, integer VERTEX_MAT(ORDER_MAX,FACE_MAX), vertex materials.
!
!    Input/output, real VERTEX_NORMAL(3,ORDER_MAX,FACE_MAX), normals at vertices.  
!
!    Input/output, real VERTEX_TEX_UV(2,ORDER_MAX,FACE_MAX), vertex textures.
!
  implicit none
  integer cor3_max
  integer face_max
  integer order_max
!
  real cor3_normal(3,cor3_max)
  integer cor3_num
  integer face(order_max,face_max)
  real face_normal(3,face_max)
  integer face_num
  integer face_order(face_max)
  integer iface
  integer ivert
  integer j
  integer m
  integer vertex_material(order_max,face_max)
  real vertex_normal(3,order_max,face_max)
  real vertex_tex_uv(2,order_max,face_max)
!
  do iface = 1, face_num

    m = face_order(iface)

    do ivert = 1, m/2

      call i_swap ( face(ivert,iface), face(m+1-ivert,iface) )
      call i_swap ( vertex_material(ivert,iface), vertex_material(m+1-ivert,iface) )

      do j = 1, 3
        call r_swap ( vertex_normal(j,ivert,iface), &
                      vertex_normal(j,m+1-ivert,iface) )
      end do

      do j = 1, 2
        call r_swap ( vertex_tex_uv(j,ivert,iface), &
                      vertex_tex_uv(j,m+1-ivert,iface) )
      end do

    end do

  end do

  cor3_normal(1:3,1:cor3_num) = - cor3_normal(1:3,1:cor3_num)
  face_normal(1:3,1:face_num) = - face_normal(1:3,1:face_num)

  write ( *, * ) ' '
  write ( *, * ) 'FACE_REVERSE_ORDER'
  write ( *, * ) '  Each list of nodes defining a face'
  write ( *, * ) '  has been reversed; related information,'
  write ( *, * ) '  including normal vectors, was also updated.'

  return
end
subroutine file_get_next_word ( iunit, word, text, text_num, ierror )
!
!*******************************************************************************
!
!! FILE_GET_NEXT_WORD returns the next word and trailing context from a file.
!
!
!  Discussion:
!
!    The file should have been opened before calling this routine.
!    The file should contain ASCII text, which can be thought of as
!    words separated by one or more blanks or commas.
!
!  Modified:
!
!    12 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters
!
!    Input, integer IUNIT, the unit number associated with the file.
!
!    Output, character ( len = * ) WORD, the next word in the file.  If the
!    current line of the file is blank, or if the file has been exhausted, 
!    WORD will be set to ' '.
!
!    Input/output, character ( len = * ) TEXT, the remaining text of the line 
!    that contains the information in WORD.  On each call, the next word
!    in TEXT is extracted until TEXT is empty, when it is refilled by
!    reading another line from the file.  Because TEXT contains information
!    needed by this routine, it should not be altered by the user
!    between calls.
!
!    Input/output, integer TEXT_NUM, the number of lines read from the file.
!    Before the first call to this routine, the user should set TEXT_NUM
!    to 0.
!
!    Output, integer IERROR, error flag. 
!    0, no error, another word was read, and returned in WORD.
!    1, end of file.  WORD and TEXT were set to ' '.
!
  implicit none
  integer ierror
  integer ihi
  integer ilo
  integer ios
  integer iunit
  integer lenc
  integer text_num
  character ( len = * ) text
  character ( len = * ) word
!
  ierror = 0

  if ( text_num <= 0 ) then
    text_num = 0
    text = ' '
  end if
!
!  If TEXT is blank, try to read a new line from the file.
!
  if ( text == ' ' ) then

    read ( iunit, '(a)', iostat = ios ) text

    if ( ios /= 0 ) then
      ierror = ios
      word = ' '
      text = ' '
      return
    end if

    text_num = text_num + 1

    if ( text == ' ' ) then
      word = ' '
      return
    end if

  end if
!
!  Extract the next word from TEXT into WORD and return.
!
  lenc = len_trim ( text )
!
!  Find ILO, the index of the first nonblank in TEXT.
!
  ilo = 1

10    continue

  if ( text(ilo:ilo) == ' ' ) then
    ilo = ilo + 1
    go to 10
  end if
!
!  Find IHI, the index of the last consecutive nonblank after the
!  one at ILO.
!
  ihi = ilo

20    continue

  if ( text(ihi:ihi) /= ',' ) then

    if ( ihi+1 <= lenc ) then
      if ( text(ihi+1:ihi+1) /= ' ' ) then
        ihi = ihi + 1
        go to 20
      end if
    end if

  end if
!
!  Set WORD.
!
  word = text(ilo:ihi)
!
!  Slide TEXT to the left.
!
  if ( ihi+1 <= lenc ) then
    text = text(ihi+1:)
  else
    text = ' '
  end if

  return
end
subroutine file_name_ext_get ( filnam, i, j )
!
!*******************************************************************************
!
!! FILE_NAME_EXT_GET determines the "extension" of a file name.
!
!
!  Definition:
!
!    The "extension" of a filename is the string of characters
!    that appears after the LAST period in the name.  A file
!    with no period, or with a period as the last character
!    in the name, has a "null" extension.
!
!  Note:
!
!    Blanks are unusual in filenames.  This routine ignores all
!    trailing blanks, but will treat initial or internal blanks
!    as regular characters acceptable in a file name.
!
!  Example:
!
!    FILNAM   I  J
!
!    bob.for    5  7
!    N.B.C.D    7  7
!    Naomi.     0  0
!    Arthur     0  0
!
!  Modified:
!
!    17 July 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILNAM, a file name to be examined.
!
!    Output, integer I, J, the indices of the first and last characters
!    in the file extension.
!
!    If at least one period occurs in the filename, and at least one
!    nonblank character follows that period, then I will be the index
!    of the first character after the period, and J the index of the
!    last nonblank character after the period.  The extension is
!    therefore equal to FILNAM(I:J).
!
!    Otherwise, I and J will be returned as 0, indicating that the file
!    has no extension.
!
  implicit none
  character ( len = * ) filnam
  integer i
  integer j
  integer s_index_last
!
  i = s_index_last ( filnam, '.' )

  if ( i /= 0 ) then

    j = len_trim ( filnam )

    if ( i == j ) then
      i = 0
      j = 0
    else
      i = i + 1
    end if

  else

    j = 0

  end if

  return
end
subroutine get_seed ( iseed )
!
!*******************************************************************************
!
!! GET_SEED returns a seed for the random number generator.
!
!
!  Discussion:
!
!    The seed depends on the current time, and ought to be (slightly)
!    different every millisecond.  Once the seed is obtained, a random
!    number generator should be called a few times to further process
!    the seed.
!
!  Modified:
!
!    23 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ISEED, a pseudorandom seed value.
!
  implicit none
  integer, parameter :: I_MAX = 2147483647
!
  integer iseed
  double precision temp
  integer values(8)
!
  character ( len = 10 ) time
  character ( len = 8 ) today
  character ( len = 5 ) zone
!
  call date_and_time ( today, time, zone, values )

  temp = 0.0

  temp = temp + dble ( values(2) - 1 ) / 11.0
  temp = temp + dble ( values(3) - 1 ) / 30.0
  temp = temp + dble ( values(5) ) / 23.0
  temp = temp + dble ( values(6) ) / 59.0
  temp = temp + dble ( values(7) ) / 59.0
  temp = temp + dble ( values(8) ) / 999.0
  temp = temp / 6.0

  if ( temp <= 0.0 ) then
    temp = 1.0 / 3.0
  else if ( temp >= 1.0 ) then
    temp = 2.0 / 3.0
  end if

  iseed = int ( dble ( I_MAX ) * temp )
!
!  Never use a seed of 0 or I_MAX.
!
  if ( iseed == 0 ) then
    iseed = 1
  end if

  if ( iseed == I_MAX ) then
    iseed = I_MAX - 1
  end if

  return
end
subroutine i_swap ( i, j )
!
!*******************************************************************************
!
!! I_SWAP switches two integer values.
!
!
!  Modified:
!
!    30 November 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer I, J.  On output, the values of I and
!    J have been interchanged.
!
  implicit none
  integer i
  integer j
  integer k
!
  k = i
  i = j
  j = k

  return
end
subroutine i_to_s_zero ( intval, s )
!
!*******************************************************************************
!
!! I_TO_S_ZERO converts an integer to a string, with zero padding.
!
!
!  Examples:
!
!    Assume that S is 6 characters long:
!
!    INTVAL  S
!
!         1  000001
!        -1  -00001
!         0  000000
!      1952  001952
!    123456  123456
!   1234567  ******  <-- Not enough room!
!
!  Modified:
!
!    04 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer INTVAL, an integer to be converted.
!
!    Output, character ( len = * ) S, the representation of the integer.
!    The integer will be right justified, and zero padded.
!    If there is not enough space, the string will be filled with stars.
!
  implicit none
  character c
  integer i
  integer idig
  integer ihi
  integer ilo
  integer intval
  integer ipos
  integer ival
  character ( len = * ) s
!
  s = ' '

  ilo = 1
  ihi = len ( s )

  if ( ihi <= 0 ) then
    return
  end if
!
!  Make a copy of the integer.
!
  ival = intval
!
!  Handle the negative sign.
!
  if ( ival < 0 ) then

    if ( ihi <= 1 ) then
      s(1:1) = '*'
      return
    end if

    ival = - ival
    s(1:1) = '-'
    ilo = 2

  end if
!
!  The absolute value of the integer goes into S(ILO:IHI).
!
  ipos = ihi

10    continue
!
!  Find the last digit of IVAL, strip it off, and stick it into
!  the string.
!
  idig = mod ( ival, 10 )
  ival = ival / 10

  if ( ipos < ilo ) then
    do i = 1, ihi
      s(i:i) = '*'
    end do
    return
  end if

  call digit_to_c ( idig, c )

  s(ipos:ipos) = c
  ipos = ipos - 1

  if ( ival /= 0 ) then
    go to 10
  end if
!
!  Fill the empties with zeroes.
!
  do i = ilo, ipos
    s(i:i) = '0'
  end do

  return
end
subroutine infile ( filein_name, ierror, filein_type )
!
!*******************************************************************************
!
!! INFILE determines the input filename and type.
!
!
!  Modified:
!
!    26 June 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = 100 ) FILEIN_NAME, the name of the input file.
!
!    Output, integer IERROR, an error flag.
!
!    Output, character ( len = 10 ) FILEIN_TYPE, the type of the file, which is
!    set to the filename extension. 
!
  implicit none
  character ( len = 100 ) filein_name
  character ( len = 10 ) filein_type
  integer i1
  integer i2
  integer ierror
  integer ios
  logical s_eqi
!
  ierror = 0

  if ( filein_name == ' ' ) then
 
    write ( *, * ) ' '
    write ( *, * ) 'INFILE:'
    write ( *, * ) '  Enter the name of a graphics file to be read:'
 
    read ( *, '(a)', iostat = ios ) filein_name
 
    if ( ios /= 0 ) then
      write ( *, * ) ' '
      write ( *, * ) 'INFILE - Error!'
      write ( *, * ) '  The input file was not specified correctly.'
      ierror = ios
      return
    end if
 
  end if
!
!  Determine the input file type.
!
  call file_name_ext_get ( filein_name, i1, i2 )

  if ( i1 /= 0 ) then
    filein_type = filein_name(i1:i2)
  else
    filein_type = ' '
  end if
 
  if ( filein_type == ' ' ) then

    write ( *, * ) ' '
    write ( *, * ) 'INFILE - Warning!'
    write ( *, * ) '  Could not the file type of the input file.'
    write ( *, * ) '  The input file name is:'
    write ( *, '(a)' ) trim ( filein_name )
    write ( *, * ) ' '
    write ( *, * ) '  The file type should occur after the period.'
    write ( *, * ) '  Please specify the file type you are using:'
    write ( *, * ) '    ase, byu, dxf, hrc, iv, ' // &
      'obj, oogl, smf, stl, stla, tri, tria, or vla:'
    read ( *, '(a)', iostat = ios ) filein_type

    if ( ios /= 0 ) then
      ierror = ios
      write ( *, * ) ' '
      write ( *, * ) 'INFILE - Fatal error!'
      write ( *, * ) '  Input failure while reading input file type.'
      stop
    end if

  end if
 
  if ( .not. ( &
    s_eqi ( filein_type, 'ASE' ) .or. &
    s_eqi ( filein_type, 'BYU' ) .or. &
    s_eqi ( filein_type, 'DXF' ) .or. &
    s_eqi ( filein_type, 'HRC' ) .or. &
    s_eqi ( filein_type, 'IV' )  .or. &
    s_eqi ( filein_type, 'OBJ' ) .or. &
    s_eqi ( filein_type, 'OOGL' ) .or. &
    s_eqi ( filein_type, 'SMF' ) .or. &
    s_eqi ( filein_type, 'STL' ) .or. &
    s_eqi ( filein_type, 'STLA' )  .or. &
    s_eqi ( filein_type, 'TRI' )  .or. &
    s_eqi ( filein_type, 'TRIA' )  .or. &
    s_eqi ( filein_type, 'VLA' ) .or. &
    s_eqi ( filein_type, 'WRL' ) ) ) then
    ierror = 1
    write ( *, * ) ' '
    write ( *, * ) 'INFILE - Error!'
    write ( *, * ) '  The input file type was not acceptable!'
    return
  end if

  return
end
subroutine ivec_max ( nval, iarray, imax )
!
!*******************************************************************************
!
!! IVEC_MAX computes the maximum element of an integer array.
!
!
!  Modified:
!
!    09 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NVAL, the number of entries in the array.
!
!    Input, integer IARRAY(NVAL), the array.
!
!    Output, integer IMAX, the value of the largest entry in the array.
!
  implicit none
  integer nval
!
  integer i
  integer iarray(nval)
  integer imax
!
  if ( nval <= 0 ) then

    imax = 0

  else

    imax = iarray(1)

    do i = 2, nval
      imax = max ( imax, iarray(i) )
    end do

  end if
 
  return
end
subroutine obj_read ( cor3, face, face_material, face_order, ierror, iseed, iunit, &
  line_dex, line_material, material_name, material_rgba, cor3_max, face_max, &
  line_max, material_max, order_max, normal_temp, bad_num, cor3_num, face_num, &
  group_num, line_num, material_num, object_num, text_num, vertex_material, &
  vertex_normal )
!
!*******************************************************************************
!
!! OBJ_READ reads graphics information from a Wavefront OBJ file.
!
!
!  Discussion:
!
!    It is intended that the information read from the file can
!    either start a whole new graphics object, or simply be added
!    to a current graphics object via the '<<' command.  
!
!    This is controlled by whether the input values have been zeroed 
!    out or not.  This routine simply tacks on the information it 
!    finds to the current graphics object.
!
!  Example:
!
!    #  magnolia.obj
!
!    mtllib ./vp.mtl
!
!    g
!    v -3.269770 -39.572201 0.876128
!    v -3.263720 -39.507999 2.160890
!    ...
!    v 0.000000 -9.988540 0.000000
!    vn 1.0 0.0 0.0
!    ...
!    vn 0.0 1.0 0.0
!    g stem
!    s 1
!    usemtl brownskn
!    f 8 9 11 10
!    f 12 13 15 14
!    ...
!    f 788 806 774
!
!  Modified:
!
!    02 June 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real COR3(3,COR3_MAX), the coordinates of points.
!
!    Input/output, integer FACE(ORDER_MAX,FACE_MAX), the nodes making faces.
!
!    Input/output, integer FACE_MATERIAL(FACE_MAX), the material of each face.
!
!    Input/output, integer FACE_ORDER(FACE_MAX), the number of vertices per face.
!
!    Output, integer IERROR, an error flag.
!
!    Input, integer IUNIT, the FORTRAN unit from which data is read.
!
!    Input/output, integer LINE_DEX(LINE_MAX), nodes forming a line, terminated 
!    by -1.
!
!    Input/output, integer LINE_MAT(LINE_MAX), material index for each line.
!
!    Input/output, character ( len = 100 ) MATERIAL_NAME(MATERIAL_MAX), material names.
!
!    Input, integer COR3_MAX, the maximum number of points.
!
!    Input, integer FACE_MAX, the maximum number of faces.
!
!    Input, integer LINE_MAX, the maximum number of line definition items.
!
!    Input, integer ORDER_MAX, the maximum number of vertices per face.
!
!    Workspace, real NORMAL_TEMP(3,ORDER_MAX*FACE_MAX).
!
!    Output, integer BAD_NUM, the number of bad text lines encountered.
!
!    Input/output, integer COR3_NUM, the number of points.
!
!    Input/output, integer FACE_NUM, the number of faces.
!
!    Input/output, integer GROUP_NUM, the number of groups.
!
!    Input/output, integer LINE_NUM, the number of line definition items.
!
!    Input/output, integer OBJECT_NUM, the number of objects.
!
!    Output, integer TEXT_NUM, the number of lines of text read from
!    the file.
!
!    Input/output, integer VERTEX_MAT(ORDER_MAX,FACE_MAX), vertex materials.
!
!    Output, real VERTEX_NORMAL(3,ORDER_MAX,FACE_MAX), normals at vertices.  
!
  implicit none
  integer, parameter :: OFFSET = 1
!
  integer cor3_max
  integer face_max
  integer line_max
  integer material_max
  integer order_max
!
  integer bad_num
  real cor3(3,cor3_max)
  integer cor3_num
  integer cor3_num_base
  logical done
  integer face(order_max,face_max)
  integer face_material(face_max)
  integer face_num
  integer face_order(face_max)
  integer group_num
  integer i
  integer i1
  integer i2
  integer ierror
  integer iseed
  integer itemp
  integer iunit
  integer ivert
  integer iword
  integer lchar
  character ( len = 256 ) line
  integer line_dex(line_max)
  integer line_material(line_max)
  integer line_num
  character ( len = 100 ) material_name(material_max)
  integer material_num
  integer imat
  real material_rgba(4,material_max)
  real normal_temp(3,order_max*face_max)
  integer object_num
  logical s_eqi
  integer text_num
  real temp
  real uniform_01_sample
  integer vertex_material(order_max,face_max)
  real vertex_normal(3,order_max,face_max)
  integer vertex_normal_num
  character ( len = 256 ) word
  character ( len = 256 ) word1
!
  ierror = 0
!
!  Save a copy of the input value of COR3_NUM to use as a base.
!
  cor3_num_base = cor3_num

  bad_num = 0
  text_num = 0

  vertex_normal_num = 0
  word = ' '
  imat = material_num
!
!  Read a line of text from the file.
!
10    continue
 
  read ( iunit, '(a)', end = 30 ) line
  text_num = text_num + 1
!
!  Replace any control characters (in particular, TAB's) by blanks.
!
  call s_control_blank ( line )
 
  done = .true.
  iword = 0
!
!  Read a word from the line.
!
  call word_nexrd ( line, word, done )
!
!  If no more words in this line, read a new line.
!
  if ( done ) then
    go to 10
  end if
!
!  If this word begins with '#' or '$', then it's a comment.  Read a new line.
!
  if ( word(1:1) == '#' .or. word(1:1) == '$' ) then
    go to 10
  end if

  iword = iword + 1
  if ( iword == 1 ) then
    word1 = word
  end if
!
!  BEVEL
!  Bevel interpolation.
!
  if ( s_eqi ( word1, 'BEVEL' ) ) then
    go to 10
!
!  BMAT
!  Basis matrix.
!
  else if ( s_eqi ( word1, 'BMAT' ) ) then
    go to 10
!
!  C_INTERP
!  Color interpolation.
!
  else if ( s_eqi ( word1, 'C_INTERP' ) ) then
    go to 10
!
!  CON
!  Connectivity between free form surfaces.
!
  else if ( s_eqi ( word1, 'CON' ) ) then
    go to 10
!
!  CSTYPE
!  Curve or surface type.
!
  else if ( s_eqi ( word1, 'CSTYPE' ) ) then
    go to 10
!
!  CTECH
!  Curve approximation technique.
!
  else if ( s_eqi ( word1, 'CTECH' ) ) then
    go to 10
!
!  CURV
!  Curve.
!
  else if ( s_eqi ( word1, 'CURV' ) ) then
    go to 10
!
!  CURV2
!  2D curve.
!
  else if ( s_eqi ( word1, 'CURV2' ) ) then
    go to 10
!
!  D_INTERP
!  Dissolve interpolation.
!
  else if ( s_eqi ( word1, 'D_INTERP' ) ) then
    go to 10
!
!  DEG
!  Degree.
!
  else if ( s_eqi ( word1, 'DEG' ) ) then
    go to 10
!
!  END
!  End statement.
!
  else if ( s_eqi ( word1, 'END' ) ) then
    go to 10
!
!  F V1 V2 V3 ...
!    or
!  F V1/VT1/VN1 V2/VT2/VN2 ...
!    or
!  F V1//VN1 V2//VN2 ...
!
!  Face.
!  A face is defined by the vertices.
!  Optionally, slashes may be used to include the texture vertex
!  and vertex normal indices.
!
  else if ( s_eqi ( word1, 'F' ) ) then

    face_num = face_num + 1

    ivert = 0

15      continue

    ivert = ivert + 1
 
    call word_nexrd ( line, word, done )

    if ( done ) then
      go to 10
    end if
!
!  Locate the slash characters in the word, if any.
!
    i1 = index ( word, '/' )
    if ( i1 > 0 ) then
      i2 = index ( word(i1+1:), '/' ) + i1
    else
      i2 = 0
    end if
!
!  Read the vertex index.
!
    call s_to_i ( word, itemp, ierror, lchar )

    if ( ierror /= 0 ) then
      itemp = -1
      ierror = 0
      write ( *, * ) 'OBJ_READ - Error!'
      write ( *, * ) '  Bad FACE field.'
      write ( *, '(a)' ) trim ( word )
    end if

    if ( ivert <= order_max .and. face_num <= face_max ) then
      face(ivert,face_num) = itemp + cor3_num_base
      vertex_material(ivert,face_num) = imat
    end if

    if ( face_num <= face_max ) then

      face_material(face_num) = imat
      face_order(face_num) = ivert

    end if
!
!  If there are two slashes, then read the data following the second one.
!
    if ( i2 > 0 ) then

      call s_to_i ( word(i2+1:), itemp, ierror, lchar )

      if ( 1 <= itemp .and. itemp <= vertex_normal_num ) then
        vertex_normal(1:3,ivert,face_num) = normal_temp(1:3,itemp)
      end if

    end if

    go to 15
!
!  G
!  Group name.
!
  else if ( s_eqi ( word1, 'G' ) ) then
    group_num = group_num + 1
    go to 10
!
!  HOLE
!  Inner trimming loop.
!
  else if ( s_eqi ( word1, 'HOLE' ) ) then
    go to 10
!
!  L
!  A line, described by a sequence of vertex indices.
!  Are the vertex indices 0 based or 1 based?
!
  else if ( s_eqi ( word1, 'L' ) ) then

25      continue
 
    call word_nexrd ( line, word, done )
!
!  If no more indices, tack a "-1" on the end.
!
    if ( done ) then

      line_num = line_num + 1

      if ( line_num <= line_max ) then
        line_dex(line_num) = -1 + OFFSET
        line_material(line_num) = -1 + OFFSET
      end if

      go to 10

    end if
!
!  Otherwise, extract the node index and add it to the line list.
!
    call s_to_i ( word, itemp, ierror, lchar )

    line_num = line_num + 1

    if ( line_num <= line_max ) then
      line_dex(line_num) = itemp
      line_material(line_num) = imat
    end if

    go to 25
!
!  LOD
!  Level of detail.
!
  else if ( s_eqi ( word1, 'LOD' ) ) then
    go to 10
!
!  MG
!  Merging group.
!
  else if ( s_eqi ( word1, 'MG' ) ) then
    go to 10
!
!  MTLLIB
!  Material library.
!
  else if ( s_eqi ( word1, 'MTLLIB' ) ) then
    go to 10
!
!  O
!  Object name.
!
  else if ( s_eqi ( word1, 'O' ) ) then
    object_num = object_num + 1
    go to 10
!
!  P
!  Point.
!
  else if ( s_eqi ( word1, 'P' ) ) then
    go to 10
!
!  PARM
!  Parameter values.
!
  else if ( s_eqi ( word1, 'PARM' ) ) then
    go to 10
!
!  S
!  Smoothing group.
!
  else if ( s_eqi ( word1, 'S' ) ) then
    go to 10
!
!  SCRV
!  Special curve.
!
  else if ( s_eqi ( word1, 'SCRV' ) ) then
    go to 10
!
!  SHADOW_OBJ
!  Shadow casting.
!
  else if ( s_eqi ( word1, 'SHADOW_OBJ' ) ) then
    go to 10
!
!  SP
!  Special point.
!
  else if ( s_eqi ( word1, 'SP' ) ) then
    go to 10
!
!  STECH
!  Surface approximation technique.
!
  else if ( s_eqi ( word1, 'STECH' ) ) then
    go to 10
!
!  STEP
!  Stepsize.
!
  else if ( s_eqi ( word1, 'STEP' ) ) then
    go to 10
!
!  SURF
!  Surface.
!
  else if ( s_eqi ( word1, 'SURF' ) ) then
    go to 10
!
!  TRACE_OBJ
!  Ray tracing.
!
  else if ( s_eqi ( word1, 'TRACE_OBJ' ) ) then
    go to 10
!
!  TRIM
!  Outer trimming loop.
!
  else if ( s_eqi ( word1, 'TRIM' ) ) then
    go to 10
!
!  USEMTL
!  Material name.
!
  else if ( s_eqi ( word1, 'USEMTL' ) ) then

    call word_nexrd ( line, word, done )

!   search for previous occurence of the same material
    imat = 0
    do i=1,min(material_num,material_max)
      if ( word == material_name(i) ) then
        imat = i
        exit
      end if
    end do

    if ( imat == 0 ) then
      material_num = material_num + 1
      imat = material_num
      if ( material_num <= material_max ) then
        material_name(material_num) = word
        material_rgba(1,material_num) = uniform_01_sample ( iseed )
        material_rgba(2,material_num) = uniform_01_sample ( iseed )
        material_rgba(3,material_num) = uniform_01_sample ( iseed )
        material_rgba(4,material_num) = 1.0
      end if
    end if

    go to 10
!
!  V X Y Z W
!  Geometric vertex.
!
!  (X,Y,Z) is the coordinate of the vertex.
!  W is optional, a weight used for rational curves and surfaces.
!  The default for W is 1.
!
  else if ( s_eqi ( word1, 'V' ) ) then

    cor3_num = cor3_num + 1

    do i = 1, 3
      call word_nexrd ( line, word, done )
      call s_to_r ( word, temp, ierror, lchar )
      if ( cor3_num <= cor3_max ) then
        cor3(i,cor3_num) = temp
      end if
    end do

    go to 10
!
!  VN
!  Vertex normals.
!
  else if ( s_eqi ( word1, 'VN' ) ) then

    vertex_normal_num = vertex_normal_num + 1

    if ( vertex_normal_num <= order_max*face_max ) then

      do i = 1, 3
        call word_nexrd ( line, word, done )
        call s_to_r ( word, temp, ierror, lchar )
        normal_temp(i,vertex_normal_num) = temp
      end do

    end if

    go to 10
!
!  VT
!  Vertex texture.
!
  else if ( s_eqi ( word1, 'VT' ) ) then
    go to 10
!
!  VP
!  Parameter space vertices.
!
  else if ( s_eqi ( word1, 'VP' ) ) then
    go to 10
!
!  Unrecognized keyword.
!
  else

    bad_num = bad_num + 1

    if ( bad_num <= 10 ) then
      write ( *, * ) ' '
      write ( *, * ) 'OBJ_READ: Bad data on line ', text_num
      write ( *, '(a)' ) '  Bad word: ' // trim ( word )
    end if

    go to 10

  end if
!
!  End of information in file.
!
30    continue
 
  return
end
function pi ( )
!
!*******************************************************************************
!
!! PI returns the value of pi.
!
!
!  Modified:
!
!    04 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real PI, the value of pi.
!
  implicit none
  real pi
!
  pi = 3.14159265358979323846264338327950288419716939937510

  return
end
subroutine r_swap ( x, y )
!
!*******************************************************************************
!
!! R_SWAP switches two real values.
!
!
!  Modified:
!
!    01 May 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real X, Y.  On output, the values of X and
!    Y have been interchanged.
!
  implicit none
  real x
  real y
  real z
!
  z = x
  x = y
  y = z

  return
end
subroutine rcol_find ( nrow, ncol, a, rvec, icol )
!
!*******************************************************************************
!
!! RCOL_FIND seeks a table column equal to a real vector.
!
!
!  Example:
!
!    Input:
!
!      NROW = 3, 
!      NCOL = 4,
!
!      A = (
!        1.  2.  3.  4.
!        5.  6.  7.  8.
!        9. 10. 11. 12. )
!
!      RVEC = ( 3., 
!               7., 
!              11. )
!
!    Output:
!
!      ICOL = 3
!
!  Modified:
!
!    09 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NROW, NCOL, the number of rows and columns in
!    the table.  NROW is also the length of RVEC.
!
!    Input, real A(NROW,NCOL), a table of numbers, regarded as
!    NCOL columns of vectors of length NROW.
!
!    Input, real RVEC(NROW), a vector to be matched with the data
!    in A.
!
!    Output, integer ICOL, the index of the first column of A
!    which exactly matches every entry of RVEC, or 0 if no match
!    could be found.
!
  implicit none
  integer ncol
  integer nrow
!
  real a(nrow,ncol)
  integer i
  integer icol
  integer j
  real rvec(nrow)
!
  icol = 0

  do j = 1, ncol

    do i = 1, nrow
      if ( rvec(i) /= a(i,j) ) then
        go to 10
      end if
    end do

    icol = j
    return

10      continue

  end do

  return
end
subroutine s_cat ( s1, s2, s3 )
!
!*******************************************************************************
!
!! S_CAT concatenates two strings to make a third string.
!
!
!  Modified:
!
!    11 May 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S1, the "prefix" string.
!
!    Input, character ( len = * ) S2, the "postfix" string.
!
!    Output, character ( len = * ) S3, the string made by
!    concatenating S1 and S2, ignoring any trailing blanks.
!
  implicit none
  character ( len = * ) s1
  character ( len = * ) s2
  character ( len = * ) s3
!
  s3 = trim ( s1 ) // trim ( s2 )

  return
end
subroutine s_control_blank ( s )
!
!*******************************************************************************
!
!! S_CONTROL_BLANK replaces control characters with blanks.
!
!
!  Definition:
!
!    A "control character" has ASCII code <= 31 or ASCII code => 127.
!
!  Modified:
!
!    14 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character ( len = * ) S, the string to be transformed.
!
  implicit none
  logical c_is_control
  integer i
  integer nchar
  character ( len = * ) s
!
  nchar = len_trim ( s )

  do i = 1, nchar
    if ( c_is_control ( s(i:i) ) ) then
      s(i:i) = ' '
    end if
  end do

  return
end
function s_eqi ( strng1, strng2 )
!
!*******************************************************************************
!
!! S_EQI is a case insensitive comparison of two strings for equality.
!
!
!  Examples:
!
!    S_EQI ( 'Anjana', 'ANJANA' ) is .TRUE.
!
!  Modified:
!
!    14 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) STRNG1, STRNG2, the strings to compare.
!
!    Output, logical S_EQI, the result of the comparison.
!
  implicit none
  integer i
  integer len1
  integer len2
  integer lenc
  logical s_eqi
  character s1
  character s2
  character ( len = * ) strng1
  character ( len = * ) strng2
!
  len1 = len ( strng1 )
  len2 = len ( strng2 )
  lenc = min ( len1, len2 )

  s_eqi = .false.

  do i = 1, lenc

    s1 = strng1(i:i)
    s2 = strng2(i:i)
    call c_cap ( s1 )
    call c_cap ( s2 )

    if ( s1 /= s2 ) then
      return
    end if

  end do

  do i = lenc + 1, len1
    if ( strng1(i:i) /= ' ' ) then
      return
    end if
  end do

  do i = lenc + 1, len2
    if ( strng2(i:i) /= ' ' ) then
      return
    end if
  end do

  s_eqi = .true.

  return
end
function s_index_last ( s, sub )
!
!*******************************************************************************
!
!! S_INDEX_LAST finds the LAST occurrence of a given substring.
!
!
!  Discussion:
!
!    It returns the location in the string at which the substring SUB is
!    first found, or 0 if the substring does not occur at all.
!
!    The routine is also trailing blank insensitive.  This is very
!    important for those cases where you have stored information in
!    larger variables.  If S is of length 80, and SUB is of
!    length 80, then if S = 'FRED' and SUB = 'RED', a match would
!    not be reported by the standard FORTRAN INDEX, because it treats
!    both variables as being 80 characters long!  This routine assumes that
!    trailing blanks represent garbage!
!
!    This means that this routine cannot be used to find, say, the last
!    occurrence of a substring 'A ', since it assumes the blank space
!    was not specified by the user, but is, rather, padding by the
!    system.  However, as a special case, this routine can properly handle
!    the case where either S or SUB is all blanks.
!
!  Modified:
!
!    14 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S, the string to be searched.
!
!    Input, character ( len = * ) SUB, the substring to search for.
!
!    Output, integer S_INDEX_LAST.  0 if SUB does not occur in
!    the string.  Otherwise S_INDEX_LAST = I, where S(I:I+LENS-1) = SUB,
!    where LENS is the length of SUB, and is the last place
!    this happens.
!
  implicit none
  integer i
  integer j
  integer llen1
  integer llen2
  character ( len = * ) s
  integer s_index_last
  character ( len = * ) sub
!
  s_index_last = 0

  llen1 = len_trim ( s )
  llen2 = len_trim ( sub )
!
!  In case S or SUB is blanks, use LEN
!
  if ( llen1 == 0 ) then
    llen1 = len ( s )
  end if

  if ( llen2 == 0 ) then
    llen2 = len ( sub )
  end if

  if ( llen2 > llen1 ) then
    return
  end if

  do j = 1, llen1+1-llen2

    i = llen1 + 2 - llen2 - j

    if ( s(i:i+llen2-1) == sub ) then
      s_index_last = i
      return
    end if

  end do

  return
end
subroutine s_to_i ( s, ival, ierror, last )
!
!*******************************************************************************
!
!! S_TO_I reads an integer value from a string.
!
!
!  Modified:
!
!    28 June 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S, a string to be examined.
!
!    Output, integer IVAL, the integer value read from the string.
!    If STRING is blank, then IVAL will be returned 0.
!
!    Output, integer IERROR, an error flag.
!    0, no error.
!    1, an error occurred.
!
!    Output, integer LAST, the last character of S used to make IVAL.
!
  implicit none
  character c
  integer i
  integer ierror
  integer isgn
  integer istate
  integer ival
  integer last
  character ( len = * ) s
!
  ierror = 0
  istate = 0
  isgn = 1
  ival = 0

  do i = 1, len_trim ( s )

    c = s(i:i)
!
!  Haven't read anything.
!
    if ( istate == 0 ) then

      if ( c == ' ' ) then

      else if ( c == '-' ) then
        istate = 1
        isgn = -1
      else if ( c == '+' ) then
        istate = 1
        isgn = + 1
      else if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
        istate = 2
        ival = ichar ( c ) - ichar ( '0' )
      else
        ierror = 1
        return
      end if
!
!  Have read the sign, expecting digits.
!
    else if ( istate == 1 ) then

      if ( c == ' ' ) then

      else if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
        istate = 2
        ival = ichar ( c ) - ichar ( '0' )
      else
        ierror = 1
        return
      end if
!
!  Have read at least one digit, expecting more.
!
    else if ( istate == 2 ) then

      if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
        ival = 10 * ival + ichar ( c ) - ichar ( '0' )
      else
        ival = isgn * ival
        last = i - 1
        return
      end if

    end if

  end do
!
!  If we read all the characters in the string, see if we're OK.
!
  if ( istate == 2 ) then
    ival = isgn * ival
    last = len_trim ( s )
  else
    ierror = 1
    last = 0
  end if

  return
end
subroutine s_to_r ( s, r, ierror, lchar )
!
!*******************************************************************************
!
!! S_TO_R reads a real number from a string.
!
!
!  Discussion:
!
!    This routine will read as many characters as possible until it reaches
!    the end of the string, or encounters a character which cannot be
!    part of the real number.
!
!    Legal input is:
!
!       1 blanks,
!       2 '+' or '-' sign,
!       2.5 spaces
!       3 integer part,
!       4 decimal point,
!       5 fraction part,
!       6 'E' or 'e' or 'D' or 'd', exponent marker,
!       7 exponent sign,
!       8 exponent integer part,
!       9 exponent decimal point,
!      10 exponent fraction part,
!      11 blanks,
!      12 final comma or semicolon.
!
!    with most quantities optional.
!
!  Examples:
!
!    S                 R
!
!    '1'               1.0
!    '     1   '       1.0
!    '1A'              1.0
!    '12,34,56'        12.0
!    '  34 7'          34.0
!    '-1E2ABCD'        -100.0
!    '-1X2ABCD'        -1.0
!    ' 2E-1'           0.2
!    '23.45'           23.45
!    '-4.2E+2'         -420.0
!    '17d2'            1700.0
!    '-14e-2'         -0.14
!    'e2'              100.0
!    '-12.73e-9.23'   -12.73 * 10.0**(-9.23)
!
!  Modified:
!
!    14 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S, the string containing the
!    data to be read.  Reading will begin at position 1 and
!    terminate at the end of the string, or when no more
!    characters can be read to form a legal real.  Blanks,
!    commas, or other nonnumeric data will, in particular,
!    cause the conversion to halt.
!
!    Output, real R, the real value that was read from the string.
!
!    Output, integer IERROR, error flag.
!
!    0, no errors occurred.
!
!    1, 2, 6 or 7, the input number was garbled.  The
!    value of IERROR is the last type of input successfully
!    read.  For instance, 1 means initial blanks, 2 means
!    a plus or minus sign, and so on.
!
!    Output, integer LCHAR, the number of characters read from
!    the string to form the number, including any terminating
!    characters such as a trailing comma or blanks.
!
  implicit none
  logical c_eqi
  character c
  integer ierror
  integer ihave
  integer isgn
  integer iterm
  integer jbot
  integer jsgn
  integer jtop
  integer lchar
  integer nchar
  integer ndig
  real r
  real rbot
  real rexp
  real rtop
  character ( len = * ) s
  character, parameter :: TAB = char ( 9 )
!
  nchar = len ( s )
  ierror = 0
  r = 0.0
  lchar = - 1
  isgn = 1
  rtop = 0.0
  rbot = 1.0
  jsgn = 1
  jtop = 0
  jbot = 1
  ihave = 1
  iterm = 0

10    continue

  lchar = lchar + 1
  c = s(lchar+1:lchar+1)
!
!  Blank or TAB character.
!
  if ( c == ' ' .or. c == TAB ) then

    if ( ihave == 2 ) then

    else if ( ihave == 6 .or. ihave == 7 ) then
      iterm = 1
    else if ( ihave > 1 ) then
      ihave = 11
    end if
!
!  Comma.
!
  else if ( c == ',' .or. c == ';' ) then

    if ( ihave /= 1 ) then
      iterm = 1
      ihave = 12
      lchar = lchar + 1
    end if
!
!  Minus sign.
!
  else if ( c == '-' ) then

    if ( ihave == 1 ) then
      ihave = 2
      isgn = - 1
    else if ( ihave == 6 ) then
      ihave = 7
      jsgn = - 1
    else
      iterm = 1
    end if
!
!  Plus sign.
!
  else if ( c == '+' ) then

    if ( ihave == 1 ) then
      ihave = 2
    else if ( ihave == 6 ) then
      ihave = 7
    else
      iterm = 1
    end if
!
!  Decimal point.
!
  else if ( c == '.' ) then

    if ( ihave < 4 ) then
      ihave = 4
    else if ( ihave >= 6 .and. ihave <= 8 ) then
      ihave = 9
    else
      iterm = 1
    end if
!
!  Exponent marker.
!
  else if ( c_eqi ( c, 'E' ) .or. c_eqi ( c, 'D' ) ) then

    if ( ihave < 6 ) then
      ihave = 6
    else
      iterm = 1
    end if
!
!  Digit.
!
  else if ( ihave < 11 .and. lge ( c, '0' ) .and. lle ( c, '9' ) ) then

    if ( ihave <= 2 ) then
      ihave = 3
    else if ( ihave == 4 ) then
      ihave = 5
    else if ( ihave == 6 .or. ihave == 7 ) then
      ihave = 8
    else if ( ihave == 9 ) then
      ihave = 10
    end if

    call c_to_digit ( c, ndig )

    if ( ihave == 3 ) then
      rtop = 10.0 * rtop + real ( ndig )
    else if ( ihave == 5 ) then
      rtop = 10.0 * rtop + real ( ndig )
      rbot = 10.0 * rbot
    else if ( ihave == 8 ) then
      jtop = 10 * jtop + ndig
    else if ( ihave == 10 ) then
      jtop = 10 * jtop + ndig
      jbot = 10 * jbot
    end if
!
!  Anything else is regarded as a terminator.
!
  else
    iterm = 1
  end if
!
!  If we haven't seen a terminator, and we haven't examined the
!  entire string, go get the next character.
!
  if ( iterm /= 1 .and. lchar+1 < nchar ) then
    go to 10
  end if
!
!  If we haven't seen a terminator, and we have examined the
!  entire string, then we're done, and LCHAR is equal to NCHAR.
!
  if ( iterm /= 1 .and. lchar+1 == nchar ) then
    lchar = nchar
  end if
!
!  Number seems to have terminated.  Have we got a legal number?
!  Not if we terminated in states 1, 2, 6 or 7!
!
  if ( ihave == 1 .or. ihave == 2 .or. ihave == 6 .or. ihave == 7 ) then

    ierror = ihave

    return
  end if
!
!  Number seems OK.  Form it.
!
  if ( jtop == 0 ) then
    rexp = 1.0
  else

    if ( jbot == 1 ) then
      rexp = 10.0**( jsgn * jtop )
    else
      rexp = jsgn * jtop
      rexp = rexp / jbot
      rexp = 10.0**rexp
    end if

  end if

  r = isgn * rexp * rtop / rbot

  return
end
subroutine stla_read ( cor3, face, face_material, face_normal, face_order, &
  ierror, iunit, cor3_max, face_max, order_max, bad_num, cor3_num, &
  dup_num, face_num, material_num, object_num, text_num, vertex_material, vertex_normal )
!
!*******************************************************************************
!
!! STLA_READ reads graphics information from an ASCII StereoLithography file.
!
!
!  Discussion:
!
!    It is intended that the information read from the file can
!    either start a whole new graphics object, or simply be added
!    to a current graphics object via the '<<' command.  
!
!    This is controlled by whether the input values have been zeroed 
!    out or not.  This routine simply tacks on the information it 
!    finds to the current graphics object.
!
!  Example:
!
!    solid MYSOLID
!      facet normal 0.4 0.4 0.2
!        outerloop
!          vertex  1.0 2.1 3.2
!          vertex  2.1 3.7 4.5
!          vertex  3.1 4.5 6.7
!      endloop
!    endfacet
!      ...
!      facet normal 0.2 0.2 0.4
!        outerloop
!          vertex  2.0 2.3 3.4
!          vertex  3.1 3.2 6.5
!          vertex  4.1 5.5 9.0
!      endloop
!    endfacet
!  endsolid MYSOLID
!
!  Modified:
!
!    02 June 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real COR3(3,COR3_MAX), the coordinates of points.
!
!    Input/output, integer FACE(ORDER_MAX,FACE_MAX), the nodes making faces.
!
!    Input/output, integer FACE_MATERIAL(FACE_MAX), the material of each face.
!
!    Input/output, real FACE_NORMAL(3,FACE_MAX), the normal vector at each face.
!
!    Input/output, integer FACE_ORDER(FACE_MAX), the number of vertices per face.
!
!    Output, integer IERROR, an error flag.
!
!    Input, integer IUNIT, the FORTRAN unit from which data is read.
!
!    Input, integer COR3_MAX, the maximum number of points.
!
!    Input, integer FACE_MAX, the maximum number of faces.
!
!    Input, integer ORDER_MAX, the maximum number of vertices per face.
!
!    Input/output, integer COR3_NUM, the number of points.
!
!    Input/output, integer FACE_NUM, the number of faces.
!
!    Input/output, integer LINE_NUM, the number of line items.
!
!    Input/output, integer VERTEX_MAT(ORDER_MAX,FACE_MAX), vertex materials.
!
!    Input/output, real VERTEX_NORMAL(3,ORDER_MAX,FACE_MAX), normals at vertices.  
!
  implicit none
  integer cor3_max
  integer face_max
  integer order_max
!
  integer bad_num
  real cor3(3,cor3_max)
  integer cor3_num
  logical done
  integer dup_num
  integer face(order_max,face_max)
  integer face_material(face_max)
  real face_normal(3,face_max)
  integer face_num
  integer face_order(face_max)
  integer i
  integer icor3
  integer ierror
  integer istate
  integer iunit
  integer ivert
  integer lchar
  integer material_num
  integer object_num
  real rval
  logical s_eqi
  real temp(3)
  character ( len = 256 ) text
  integer text_num
  integer vertex_material(order_max,face_max)
  real vertex_normal(3,order_max,face_max)
  character ( len = 256 ) word1
  character ( len = 256 ) word2
!
  ierror = 0
  istate = 0
!
!  Read the next line of text.
!
10    continue
 
  read ( iunit, '(a)', end = 30 ) text
  text_num = text_num + 1
  done = .true.
!
!  Read the first word in the line.
!
  call word_nexrd ( text, word1, done )
!
!  "Doctor" the text, changing a beginning occurrence of:
!
!      END FACET to ENDFACET
!      END LOOP to ENDLOOP
!      END SOLID to ENDSOLID
!      FACET NORMAL to FACETNORMAL
!      OUTER LOOP to OUTERLOOP
!
  if ( s_eqi ( word1, 'END' ) .or. s_eqi ( word1, 'FACET' ) .or. &
       s_eqi ( word1, 'OUTER' ) ) then

    call word_nexrd ( text, word2, done )
    call s_cat ( word1, word2, word1 )

  end if
!
!  This first word tells us what to do.
!
!  SOLID - begin a new solid.
!    Valid in state 0, moves to state 1.
!  ENDSOLID - end current solid.
!    Valid in state 1, moves to state 0.
!
!  FACETNORMAL - begin a new facet.
!    Valid in state 0 or 1, moves to state 2.
!  ENDFACET - end current facet.
!    Valid in state 2, moves to state 1.
!
!  OUTERLOOP - begin a list of vertices.
!    Valid in state 2, moves to state 3.
!  ENDLOOP - end vertex list.
!    Valid in state 3, moves to state 2.
!
!  VERTEX - give coordinates of next vertex.
!    Valid in state 3.
!
!  End of file - 
!    Valid in state 0 or 1.
!
  if ( s_eqi ( word1, 'SOLID' ) ) then

    if ( istate == 0 ) then
      istate = 1
      object_num = object_num + 1
    else
      write ( *, * ) ' '
      write ( *, * ) 'STLA_READ - Warning!'
      write ( *, * ) '  Model not in right state for SOLID.'
      bad_num = bad_num + 1
      ierror = 1
      go to 40
    end if

  else if ( s_eqi ( word1, 'ENDSOLID' ) ) then

    if ( istate == 1 ) then
      istate = 0
    else
      write ( *, * ) ' '
      write ( *, * ) 'STLA_READ - Warning!'
      write ( *, * ) '  Model not in right state for ENDSOLID.'
      bad_num = bad_num + 1
      ierror = 1
      go to 40
    end if

  else if ( s_eqi ( word1, 'FACETNORMAL' ) ) then

    if ( istate == 0 .or. istate == 1 ) then

      istate = 2
      face_num = face_num + 1

      if ( face_num <= face_max ) then

        face_material(face_num) = material_num
        face_order(face_num) = 0

        do i = 1, 3
          face_normal(i,face_num) = 0.0
          call word_nexrd ( text, word2, done )
          if ( .not. done ) then
            call s_to_r ( word2, rval, ierror, lchar )
            if ( ierror == 0 ) then
              face_normal(i,face_num) = rval
            end if
          end if
        end do

      end if

    else

      write ( *, * ) ' '
      write ( *, * ) 'STLA_READ - Warning!'
      write ( *, * ) '  Model not in right state for FACET.'
      bad_num = bad_num + 1
      ierror = 1
      go to 40

    end if

  else if ( s_eqi ( word1, 'ENDFACET' ) ) then

    if ( istate == 2 ) then
      istate = 1
    else
      write ( *, * ) ' '
      write ( *, * ) 'STLA_READ - Warning!'
      write ( *, * ) '  Model not in right state for ENDFACET.'
      bad_num = bad_num + 1
      ierror = 1
      go to 40
    end if

  else if ( s_eqi ( word1, 'OUTERLOOP' ) ) then

    if ( istate == 2 ) then
      istate = 3
      ivert = 0
    else
      write ( *, * ) ' '
      write ( *, * ) 'STLA_READ - Warning!'
      write ( *, * ) '  Model not in right state for OUTERLOOP.'
      bad_num = bad_num + 1
      ierror = 1
      go to 40
    end if

  else if ( s_eqi ( word1, 'ENDLOOP' ) ) then

    if ( istate == 3 ) then
      istate = 2
    else
      write ( *, * ) ' '
      write ( *, * ) 'STLA_READ - Warning!'
      write ( *, * ) '  Model not in right state for ENDLOOP.'
      bad_num = bad_num + 1
      ierror = 1
      go to 40
    end if

  else if ( s_eqi ( word1, 'VERTEX' ) ) then

    if ( istate == 3 ) then

      do i = 1, 3
        call word_nexrd ( text, word2, done )
        call s_to_r ( word2, rval, ierror, lchar )
        temp(i) = rval
      end do
!
!  If the coordinate values already exist in COR3, then  
!  save space by using the index of a previous copy.
!
      if ( cor3_num <= 1000 ) then
        call rcol_find ( 3, cor3_num, cor3, temp, icor3 )
      else
        icor3 = 0
      end if

      if ( icor3 == 0 ) then
        cor3_num = cor3_num + 1
        icor3 = cor3_num
        if ( cor3_num <= cor3_max ) then
          cor3(1:3,cor3_num) = temp(1:3)
        end if
      else
        dup_num = dup_num + 1
      end if

      ivert = ivert + 1

      if ( ivert <= order_max .and. face_num <= face_max ) then
        face(ivert,face_num) = icor3
        vertex_material(ivert,face_num) = material_num
        vertex_normal(1:3,ivert,face_num) = face_normal(1:3,face_num)
      end if

      if ( face_num <= face_max .and. face_order(face_num) < order_max ) then

        face_order(face_num) = face_order(face_num) + 1

      end if

    else

      write ( *, * ) ' '
      write ( *, * ) 'STLA_READ - Warning!'
      write ( *, * ) '  Model not in right state for VERTEX.'
      bad_num = bad_num + 1
      ierror = 1
      go to 40

    end if

  else

    write ( *, * ) ' '
    write ( *, * ) 'STLA_READ - Warning!'
    write ( *, * ) '  Unrecognized line in file.'
    bad_num = bad_num + 1
    ierror = 1
    go to 40

  end if
 
  go to 10
!
!  Come here on end of file.
!
30    continue

  if ( istate /= 0 .and. istate .ne. 1 ) then
    ierror = 1
    write ( *, * ) ' '
    write ( *, * ) 'STLA_READ - Warning.'
    write ( *, * ) '  End-of-file, but model not finished.'
  end if

40    continue

  return
end
subroutine tmat_init ( a )
!
!*******************************************************************************
!
!! TMAT_INIT initializes the geometric transformation matrix.
!
!
!  Definition:
!
!    The geometric transformation matrix can be thought of as a 4 by 4
!    matrix "A" having components:
!
!      r11 r12 r13 t1
!      r21 r22 r23 t2
!      r31 r32 r33 t3
!        0   0   0  1
!
!    This matrix encodes the rotations, scalings and translations that
!    are applied to graphical objects.
!
!    A point P = (x,y,z) is rewritten in "homogeneous coordinates" as 
!    PH = (x,y,z,1).  Then to apply the transformations encoded in A to 
!    the point P, we simply compute A * PH.
!
!    Individual transformations, such as a scaling, can be represented
!    by simple versions of the transformation matrix.  If the matrix
!    A represents the current set of transformations, and we wish to 
!    apply a new transformation B, then the original points are
!    transformed twice:  B * ( A * PH ).  The new transformation B can
!    be combined with the original one A, to give a single matrix C that
!    encodes both transformations: C = B * A.
!
!  Modified:
!
!    19 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Foley, van Dam, Feiner, Hughes,
!    Computer Graphics, Principles and Practice,
!    Addison Wesley, Second Edition, 1990.
!
!  Parameters:
!
!    Input, real A(4,4), the geometric transformation matrix.
!
  implicit none
  real a(4,4)
  integer i
  integer j
!
  do i = 1, 4
    do j = 1, 4
      if ( i == j ) then
        a(i,j) = 1.0
      else
        a(i,j) = 0.0
      end if
    end do
  end do

  return
end
subroutine tmat_mxm ( a, b, c )
!
!*******************************************************************************
!
!! TMAT_MXM multiplies two geometric transformation matrices.
!
!
!  Note:
!
!    The product is accumulated in a temporary array, and then assigned
!    to the result.  Therefore, it is legal for any two, or all three,
!    of the arguments to share memory.
!
!  Modified:
!
!    19 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Foley, van Dam, Feiner, Hughes,
!    Computer Graphics, Principles and Practice,
!    Addison Wesley, Second Edition, 1990.
!
!  Parameters:
!
!    Input, real A(4,4), the first geometric transformation matrix.
!
!    Input, real B(4,4), the second geometric transformation matrix.
!
!    Output, real C(4,4), the product A * B.
!
  implicit none
  real a(4,4)
  real b(4,4)
  real c(4,4)
  real d(4,4)
  integer i
  integer j
  integer k
!
  do i = 1, 4
    do k = 1, 4
      d(i,k) = 0.0
      do j = 1, 4
        d(i,k) = d(i,k) + a(i,j) * b(j,k)
      end do
    end do
  end do

  c(1:4,1:4) = d(1:4,1:4)

  return
end
subroutine tmat_mxv ( a, x, y )
!
!*******************************************************************************
!
!! TMAT_MXV multiplies a geometric transformation matrix times a vector.
!
!
!  Modified:
!
!    12 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Foley, van Dam, Feiner, Hughes,
!    Computer Graphics, Principles and Practice,
!    Addison Wesley, Second Edition, 1990.
!
!  Parameters:
!
!    Input, real A(4,4), the geometric transformation matrix.
!
!    Input, real X(3), the vector to be multiplied.  The fourth component
!    of X is implicitly assigned the value of 1.
!
!    Output, real Y(3), the result of A*X.  The product is accumulated in 
!    a temporary vector, and then assigned to the result.  Therefore, it 
!    is legal for X and Y to share memory.
!
  implicit none
  real a(4,4)
  integer i
  integer j
  real x(3)
  real y(3)
  real z(3)
!
  do i = 1, 3
    z(i) = 0.0
    do j = 1, 3
      z(i) = z(i) + a(i,j) * x(j)
    end do
    z(i) = z(i) + a(i,4)
  end do

  y(1:3) = z(1:3)

  return
end
subroutine tmat_rot_vector ( a, b, angle, axis )
!
!*******************************************************************************
!
!! TMAT_ROT_VECTOR applies an arbitrary axis rotation to the geometric transformation matrix.
!
!
!  Modified:
!
!    11 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Foley, van Dam, Feiner, Hughes,
!    Computer Graphics, Principles and Practice,
!    Addison Wesley, Second Edition, 1990.
!
!  Parameters:
!
!    Input, real A(4,4), the current geometric transformation matrix.
!
!    Output, real B(4,4), the modified geometric transformation matrix.
!    A and B may share the same memory.
!
!    Input, real ANGLE, the angle, in degrees, of the rotation.
!
!    Input, real AXIS(3), the axis vector about which rotation occurs.
!    AXIS may not be the zero vector.
!
  implicit none
  real a(4,4)
  real angle
  real angle_rad
  real axis(3)
  real b(4,4)
  real c(4,4)
  real ca
  real d(4,4)
  real degrees_to_radians
  real norm
  real sa
  real v1
  real v2
  real v3
!
  v1 = axis(1)
  v2 = axis(2)
  v3 = axis(3)

  norm = sqrt ( v1 * v1 + v2 * v2 + v3 * v3 )

  if ( norm == 0.0 ) then
    return
  end if

  v1 = v1 / norm
  v2 = v2 / norm
  v3 = v3 / norm

  angle_rad = degrees_to_radians ( angle )
  ca = cos ( angle_rad )
  sa = sin ( angle_rad )

  call tmat_init ( c )

  c(1,1) =                v1 * v1 + ca * ( 1.0 - v1 * v1 )
  c(1,2) = ( 1.0 - ca ) * v1 * v2 - sa * v3
  c(1,3) = ( 1.0 - ca ) * v1 * v3 + sa * v2

  c(2,1) = ( 1.0 - ca ) * v2 * v1 + sa * v3
  c(2,2) =                v2 * v2 + ca * ( 1.0 - v2 * v2 )
  c(2,3) = ( 1.0 - ca ) * v2 * v3 - sa * v1

  c(3,1) = ( 1.0 - ca ) * v3 * v1 - sa * v2
  c(3,2) = ( 1.0 - ca ) * v3 * v2 + sa * v1
  c(3,3) =                v3 * v3 + ca * ( 1.0 - v3 * v3 )

  call tmat_mxm ( c, a, d )

  b(1:4,1:4) = d(1:4,1:4)

  return
end
subroutine tmat_scale ( a, b, sx, sy, sz )
!
!*******************************************************************************
!
!! TMAT_SCALE applies a scaling to the geometric transformation matrix.
!
!
!  Modified:
!
!    19 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Foley, van Dam, Feiner, Hughes,
!    Computer Graphics, Principles and Practice,
!    Addison Wesley, Second Edition, 1990.
!
!  Parameters:
!
!    Input, real A(4,4), the current geometric transformation matrix.
!
!    Output, real B(4,4), the modified geometric transformation matrix.
!    A and B may share the same memory.
!
!    Input, real SX, SY, SZ, the scalings to be applied to the X, Y and
!    Z coordinates.
!
  implicit none
  real a(4,4)
  real b(4,4)
  real c(4,4)
  real d(4,4)
  real sx
  real sy
  real sz
!
  call tmat_init ( c )

  c(1,1) = sx
  c(2,2) = sy
  c(3,3) = sz

  call tmat_mxm ( c, a, d )

  b(1:4,1:4) = d(1:4,1:4)

  return
end
subroutine tmat_trans ( a, b, x, y, z )
!
!*******************************************************************************
!
!! TMAT_TRANS applies a translation to the geometric transformation matrix.
!
!
!  Modified:
!
!    19 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Foley, van Dam, Feiner, Hughes,
!    Computer Graphics, Principles and Practice,
!    Addison Wesley, Second Edition, 1990.
!
!  Parameters:
!
!    Input, real A(4,4), the current geometric transformation matrix.
!
!    Output, real B(4,4), the modified transformation matrix.
!    A and B may share the same memory.
!
!    Input, real X, Y, Z, the translation.  This may be thought of as the
!    point that the origin moves to under the translation.
!
  implicit none
  real a(4,4)
  real b(4,4)
  integer i
  integer j
  real x
  real y
  real z
!
  do i = 1, 4
    do j = 1, 4

      if ( i == 1 .and. j == 4 ) then
        b(1,4) = a(1,4) + x
      else if ( i == 2 .and. j == 4 ) then
        b(2,4) = a(2,4) + y
      else if ( i == 3 .and. j == 4 ) then
        b(3,4) = a(3,4) + z
      else
        b(i,j) = a(i,j)
      end if

    end do
  end do

  return
end
function uniform_01_sample ( iseed )
!
!*******************************************************************************
!
!! UNIFORM_01_SAMPLE is a portable random number generator.
!
!
!  Formula:
!
!    ISEED = ISEED * (7**5) mod (2**31 - 1)
!    RANDOM = ISEED * / ( 2**31 - 1 )
!
!  Modified:
!
!    01 March 1999
!
!  Parameters:
!
!    Input/output, integer ISEED, the integer "seed" used to generate
!    the output random number, and updated in preparation for the
!    next one.  ISEED should not be zero.
!
!    Output, real UNIFORM_01_SAMPLE, a random value between 0 and 1.
!

!
!  IA = 7**5
!  IB = 2**15
!  IB16 = 2**16
!  IP = 2**31-1
!
  implicit none
  integer, parameter :: ia = 16807
  integer, parameter :: ib15 = 32768
  integer, parameter :: ib16 = 65536
  integer, parameter :: ip = 2147483647
!
  integer iprhi
  integer iseed
  integer ixhi
  integer k
  integer leftlo
  integer loxa
  real uniform_01_sample
!
!  Don't let ISEED be 0.
!
  if ( iseed == 0 ) then
    iseed = ip
  end if
!
!  Get the 15 high order bits of ISEED.
!
  ixhi = iseed / ib16
!
!  Get the 16 low bits of ISEED and form the low product.
!
  loxa = ( iseed - ixhi * ib16 ) * ia
!
!  Get the 15 high order bits of the low product.
!
  leftlo = loxa / ib16
!
!  Form the 31 highest bits of the full product.
!
  iprhi = ixhi * ia + leftlo
!
!  Get overflow past the 31st bit of full product.
!
  k = iprhi / ib15
!
!  Assemble all the parts and presubtract IP.  The parentheses are
!  essential.
!
  iseed = ( ( ( loxa - leftlo * ib16 ) - ip ) + ( iprhi - k * ib15 ) * ib16 ) &
    + k
!
!  Add IP back in if necessary.
!
  if ( iseed < 0 ) then
    iseed = iseed + ip
  end if
!
!  Multiply by 1 / (2**31-1).
!
  uniform_01_sample = real ( iseed ) * 4.656612875e-10

  return
end
subroutine vertex_normal_set ( cor3, face, face_order, cor3_max, face_max, &
  order_max, face_num, vertex_normal )
!
!*******************************************************************************
!
!! VERTEX_NORMAL_SET recomputes the face vertex normal vectors.
!
!
!  Modified:
!
!    12 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real COR3(3,COR3_MAX), the coordinates of points.
!
!    Input, integer FACE(ORDER_MAX,FACE_MAX), the nodes making faces.
!
!    Input, integer FACE_ORDER(FACE_MAX), the number of vertices per face.
!
!    Input, integer COR3_MAX, the maximum number of points.
!
!    Input, integer FACE_MAX, the maximum number of faces.
!
!    Input, integer ORDER_MAX, the maximum number of vertices per face.
!
!    Input, integer FACE_NUM, the number of faces.
!
!    Input/output, real VERTEX_NORMAL(3,ORDER_MAX,FACE_MAX), normal vectors
!    at vertices.  
!
  implicit none
  integer cor3_max
  integer face_max
  integer order_max
!
  real cor3(3,cor3_max)
  integer face(order_max,face_max)
  integer face_order(face_max)
  integer i
  integer i0
  integer i1
  integer i2
  integer iface
  integer ivert
  integer jp1
  integer jp2
  real norm
  integer face_num
  integer fix_num
  integer vec_num
  integer zero_num
  real vertex_normal(3,order_max,face_max)
  real x0
  real x1
  real x2
  real xc
  real y0
  real y1
  real y2
  real yc
  real z0
  real z1
  real z2
  real zc
!
  if ( face_num <= 0 ) then
    return
  end if

  fix_num = 0
  vec_num = 0
  zero_num = 0
!
!  Consider each face.
!
  do iface = 1, face_num
!
!  Consider each vertex.
!
    do ivert = 1, face_order(iface)

      vec_num = vec_num + 1

      norm = 0.0
      do i = 1, 3
        norm = norm + vertex_normal(i,ivert,iface)**2
      end do
      norm = sqrt ( norm )

      if ( norm == 0.0 ) then

        fix_num = fix_num + 1

        i0 = face(ivert,iface)
        x0 = cor3(1,i0)
        y0 = cor3(2,i0)
        z0 = cor3(3,i0)

        jp1 = ivert + 1
        if ( jp1 > face_order(iface) ) then
          jp1 = jp1 - face_order(iface)
        end if
        i1 = face(jp1,iface)
        x1 = cor3(1,i1)
        y1 = cor3(2,i1)
        z1 = cor3(3,i1)

        jp2 = ivert + 2
        if ( jp2 > face_order(iface) ) then
          jp2 = jp2 - face_order(iface)
        end if
        i2 = face(jp2,iface)
        x2 = cor3(1,i2)
        y2 = cor3(2,i2)
        z2 = cor3(3,i2)

        call cross0_3d ( x0, y0, z0, x1, y1, z1, x2, y2, z2, xc, yc, zc )

        norm = sqrt ( xc * xc + yc * yc + zc * zc )

        if ( norm == 0.0 ) then
          zero_num = zero_num + 1
          xc = 1.0 / sqrt ( 3.0 )
          yc = 1.0 / sqrt ( 3.0 )
          zc = 1.0 / sqrt ( 3.0 )
        else
          xc = xc / norm
          yc = yc / norm
          zc = zc / norm
        end if

        vertex_normal(1,ivert,iface) = xc
        vertex_normal(2,ivert,iface) = yc
        vertex_normal(3,ivert,iface) = zc

      end if

    end do

  end do

  write ( *, * ) ' '
  write ( *, * ) 'VERTEX_NORMAL_SET - Note:'
  write ( *, * ) '  There are a total of ', vec_num, ' vertex normal vectors.'
  write ( *, * ) '  Initially, ', fix_num, ' were zero.'
  write ( *, * ) '  Of these, ', zero_num, ' represent zero areas.'

  return
end
subroutine vrml_read ( cor3, cor3_material, face, face_material, face_order, &
  ierror, iunit, material_name, material_rgba, cor3_max, face_max, &
  material_max, order_max, bad_num, cor3_num, face_num, line_num, &
  material_num, text_num, vertex_material )
!
!*******************************************************************************
!
!! VRML_READ reads graphics information from a VRML file.
!
!
!  Modified:
!
!    12 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
  implicit none
  logical, parameter :: debug = .false.
  integer, parameter :: ivec_max = 20
  integer, parameter :: level_max = 10
  integer, parameter :: rvec_max = 20
  integer, parameter :: offset = 1
!
  integer cor3_max
  integer face_max
  integer order_max
  integer material_max

  real angle
  real axis(3)
  integer bad_num
  character ( len = 4 ) char4
  real cor3(3,cor3_max)
  integer cor3_material(cor3_max)
  integer color_mat(cor3_max)
  real rgba(4)
  integer cor3_num
  integer cor3_num_old
  integer face(order_max,face_max)
  integer face_material(face_max)
  integer face_num
  integer face_num_old
  integer face_order(face_max)
  integer i
  integer icor3
  integer ierror
  integer iface
  integer imat
  integer iunit
  integer ival
  integer ivec(ivec_max)
  integer ivec_num
  integer ivert
  integer j
  integer k
  integer lchar
  integer level
  character ( len = 256 ) level_name(0:level_max)
  integer line_num
  integer line_num_old
  character ( len = 100 ) material_binding
  character ( len = 100 ) material_name(material_max)
  integer material_num
  integer mat_stored
  real material_rgba(4,material_max)
  integer new_color
  integer new_color_index
  integer cor3_new
  integer new_face
  integer node
  integer nlbrack
  integer nrbrack
  integer overall_mat
  real r01
  real r02
  real r03
  real r04
  real r05
  real r06
  real r07
  real r08
  real r09
  real r10
  real r11
  real r12
  real rval
  real rvec(rvec_max)
  integer rvec_num
  real rx
  real ry
  real rz
  logical s_eqi
  real sx
  real sy
  real sz
  character ( len = 256 ) text
  integer text_num
  real transform_matrix(4,4)
  real tx
  real ty
  real tz
  integer vertex_material(order_max,face_max)
  character ( len = 256 ) word
  character ( len = 256 ) wordm1
!
  ierror = 0
  level = 0
  level_name(0) = 'Top'
  nlbrack = 0
  nrbrack = 0
  ivec_num = 0
  rvec_num = 0
!
!  Save old counts in order to be able to add a new object to a set
!  of existing ones.
!
  cor3_num_old = cor3_num
  face_num_old = face_num
  line_num_old = line_num
!
!  Initialize the transformation matrix.
!
  call tmat_init ( transform_matrix )

  word = ' '
  text = ' '
!
!  Read the next word and its trailing context, from the input file.
!
10    continue

  wordm1 = word

11    continue

  call file_get_next_word ( iunit, word, text, text_num, ierror )

  if ( ierror /= 0 ) then
    go to 50
  end if

  if ( debug ) then
    write ( *, '(a)' ) trim ( word )
  end if
!
!  The first line of the file must be the header.
!
  if ( text_num == 1 ) then

    if ( s_eqi ( word, '#VRML' ) ) then
      call file_get_next_word ( iunit, word, text, text_num, ierror )
      if ( s_eqi ( word, 'V2.0' ) ) then
        call file_get_next_word ( iunit, word, text, text_num, ierror )
        if ( s_eqi ( word, 'utf8' ) ) then
          go to 10
        end if
      end if
    end if

    ierror = 1
    write ( *, * ) ' '
    write ( *, * ) 'VRML_READ - Fatal error!'
    write ( *, * ) '  The input file has a bad header.'
    write ( *, '(a)' ) trim ( text )
    return

  end if
!
!  Skip a comment and all text following it on a line.
!
  if ( word(1:1) == '#' ) then
    text = ' '
    go to 10
  end if
!
!  Ignore a blank line.
!
  if ( word == ' ' ) then
    go to 10
  end if
!
!  Ignore an isolated comma.
!
  if ( word == ',' ) then
    go to 10
  end if
!
!  If the word is a curly or square bracket, count it.
!
  if ( word == '{' .or. word == '[' ) then

    nlbrack = nlbrack + 1

  else if ( word .eq. '}' .or. word == ']' ) then

    nrbrack = nrbrack + 1

    if ( nlbrack < nrbrack ) then
      write ( *, * ) ' '
      write ( *, * ) 'VRML_READ - Fatal error!'
      write ( *, * ) '  Extraneous right bracket, line ', text_num
      write ( *, '(a)' ) trim ( text )
      write ( *, * ) 'Currently processing field:'
      write ( *, '(a)' ) trim ( level_name(level) )
      ierror = 1
      return
    end if

  end if
!
!  If the word is DEF, then read the next word right now,
!  and DON'T copy WORD into WORDM1.  
!
  if ( s_eqi ( word, 'DEF' ) ) then
    call file_get_next_word ( iunit, word, text, text_num, ierror )
    write ( *, * ) 'Skipping DEF ' // trim ( word )
    go to 11
  end if
!
!  If the word is a left bracket, then the previous word
!  is the name of a node.
!
  if ( word == '{' .or. word == '[' ) then

    level = nlbrack - nrbrack
    if ( level < 0 ) then
      write ( *, * ) 'Too many right brackets!'
      level = 0
    else if ( level > level_max ) then
      write ( *, * ) 'Too many left brackets!'
      level = level_max
    end if

    level_name(level) = wordm1

    if ( debug ) then
      write ( *, * ) ' '
      do i = 0, level
        write ( *, '(i3,2x,a)' ) i, trim ( level_name(i) )
      end do
    end if

  end if
!
!  ANCHOR
!
  if ( s_eqi ( level_name(level), 'Anchor' ) ) then

    if ( word == '{' ) then

    else if ( word == '}' ) then

      level = nlbrack - nrbrack

    end if
!
!  APPEARANCE
!
  else if ( s_eqi ( level_name(level), 'Appearance' ) ) then

    if ( word == '{' ) then

    else if ( word == '}' ) then

      level = nlbrack - nrbrack

    else if ( s_eqi ( word, 'material' ) ) then
    else if ( s_eqi ( word, 'texture' ) ) then
    else if ( s_eqi ( word, 'textureTransform' ) ) then

    end if
!
!  AUDIOCLIP
!
  else if ( s_eqi ( level_name(level), 'AudioClip' ) ) then

    if ( word == '{' ) then

    else if ( word == '}' ) then

      level = nlbrack - nrbrack

    end if
!
!  BACKGROUND
!
  else if ( s_eqi ( level_name(level), 'Background' ) ) then

    if ( word == '{' ) then

    else if ( word == '}' ) then

      level = nlbrack - nrbrack

    else if ( s_eqi ( word, 'skyColor' ) ) then
    else if ( s_eqi ( word, 'skyAngle' ) ) then
    else if ( s_eqi ( word, 'groundColor' ) ) then
    else if ( s_eqi ( word, 'groundAngle' ) ) then

    end if
!
!  BILLBOARD
!
  else if ( s_eqi ( level_name(level), 'Billboard' ) ) then

    if ( word == '{' ) then

    else if ( word == '}' ) then

      level = nlbrack - nrbrack

    end if
!
!  BOX
!
  else if ( s_eqi ( level_name(level), 'Box' ) ) then

    if ( word == '{' ) then

    else if ( word == '}' ) then

      level = nlbrack - nrbrack

    end if
!
!  COLLISION
!
  else if ( s_eqi ( level_name(level), 'Collision' ) ) then

    if ( word == '{' ) then

    else if ( word == '}' ) then

      level = nlbrack - nrbrack

    end if
!
!  COLOR
!
  else if ( level_name(level) == 'Color' ) then

    if ( word == '{' ) then

    else if ( word == '}' ) then

      level = nlbrack - nrbrack

    else if ( word == 'color' ) then

      write ( *, * ) 'DEBUG: COLOR saw color.'

    end if
!
!  COLOR { COLOR [] }
!
  else if ( level_name(level) == 'color' ) then

    if ( level_name(level-1) == 'Color' ) then

      if ( word == '[' ) then

        rvec_num = 0
        new_color = 0

      else if ( word == ']' ) then

        if ( rvec_num == 3 ) then

          write (*,*) 'ERROR during reading of WRL format'
          stop
      
        end if

        rvec_num = 0
        level = nlbrack - nrbrack

      else

        call s_to_r ( word, rval, ierror, lchar )
        rvec_num = rvec_num + 1
        rvec(rvec_num) = rval

        if ( rvec_num == 3 ) then
!         if the following is not satisfied, it will lead to silent erroneous behavior 
          if ( new_color <= cor3_max ) then
            new_color = new_color + 1
          
            rgba(1:3) = rvec(1:3)
            rgba(4) = 1.0
!           Search current color among previous ones
            mat_stored = min(material_num,material_max)
            call rcol_find ( 4, mat_stored, material_rgba, rgba, imat )
      
            if ( imat == 0 ) then
              material_num = material_num + 1
              imat = material_num
              
              if ( material_num <= material_max ) then
                material_rgba(1:4,material_num) = rgba(1:4)
                call i_to_s_zero ( material_num, char4 )
                material_name(material_num) = 'Material_' // char4
              end if
            end if
            color_mat(new_color) = imat
            
          end if
          rvec_num = 0

        end if

      end if

    end if
!
!  COLORINDEX
!
  else if ( s_eqi ( level_name(level), 'colorIndex' ) ) then

    if ( word == '[' ) then

      ivec_num = 0

    else if ( word == ']' ) then

      write ( *, * ) 'Hey, IVEC_NUM is ', ivec_num
      new_color_index = ivec_num

!         if ( ivec_num /= 0 ) then

!           face_num = face_num + 1

!           if ( face_num <= face_max ) then

!             if ( ivec_num > order_max ) then
!               ivec_num = order_max
!             end if

!             do i = 1, ivec_num
!               face(i,face_num) = ivec(i)
!             end do
!             face_order(face_num) = ivec_num
!             ivec_num = 0

!           end if

!         end if

      level = nlbrack - nrbrack

    else

      call s_to_i ( word, ival, ierror, lchar )

      if ( ivec_num < ivec_max ) then
        ivec_num = ivec_num + 1
        ivec(ivec_num) = ival + OFFSET
      end if

    end if
!
!  COLORINTERPOLATOR
!
  else if ( s_eqi ( level_name(level), 'ColorInterpolator' ) ) then

    if ( word == '{' ) then

    else if ( word == '}' ) then

      level = nlbrack - nrbrack

    end if
!
!  CONE
!
  else if ( s_eqi ( level_name(level), 'CONE' ) ) then

    if ( word == '{' ) then

    else if ( word == '}' ) then

      level = nlbrack - nrbrack

    end if
!
!  COORDINATE
!
  else if ( s_eqi ( level_name(level), 'Coordinate' ) ) then

    if ( word == '{' ) then

    else if ( word == '}' ) then

      level = nlbrack - nrbrack

    else if ( s_eqi ( word, 'point' ) ) then

    end if
!
!  COORDINATE POINT
!
  else if ( s_eqi ( level_name(level), 'point' ) ) then

    if ( s_eqi ( level_name(level-1), 'Coordinate' ) ) then

      if ( word == '[' ) then

        rvec_num = 0

      else if ( word == ']' ) then

        if ( rvec_num == 3 ) then

          cor3_num = cor3_num + 1

          if ( cor3_num <= cor3_max ) then

            call tmat_mxv ( transform_matrix, rvec, rvec )

            cor3(1:3,cor3_num) = rvec(1:3)

          end if

        end if

        rvec_num = 0
        level = nlbrack - nrbrack

      else

        call s_to_r ( word, rval, ierror, lchar )
        rvec_num = rvec_num + 1
        rvec(rvec_num) = rval

        if ( rvec_num == 3 ) then

          cor3_num = cor3_num + 1

          if ( cor3_num <= cor3_max ) then
            call tmat_mxv ( transform_matrix, rvec, rvec )
            cor3(1:3,cor3_num) = rvec(1:3)
          end if

          rvec_num = 0

        end if

      end if

    end if
!
!  COORDINATEINTERPOLATOR
!
  else if ( s_eqi ( level_name(level), 'CoordinateInterpolator' ) ) then

    if ( word == '{' ) then

    else if ( word == '}' ) then

      level = nlbrack - nrbrack

    end if
!
!  COORDINDEX
!
  else if ( s_eqi ( level_name(level), 'coordIndex' ) ) then

    if ( word == '[' ) then

      ivec_num = 0

    else if ( word == ']' ) then

      if ( ivec_num /= 0 ) then

        face_num = face_num + 1

        if ( face_num <= face_max ) then

          if ( ivec_num > order_max ) then
            ivec_num = order_max
          end if

          do i = 1, ivec_num
            face(i,face_num) = ivec(i)
          end do
          face_order(face_num) = ivec_num
          ivec_num = 0

        end if

      end if

      level = nlbrack - nrbrack

    else

      call s_to_i ( word, ival, ierror, lchar )

      if ( ival /= -1 ) then

        if ( ivec_num < ivec_max ) then
          ivec_num = ivec_num + 1
          ivec(ivec_num) = ival + cor3_num_old + OFFSET
        end if

      else

        face_num = face_num + 1

        if ( face_num <= face_max ) then

          if ( ivec_num > order_max ) then
            ivec_num = order_max
          end if

          do i = 1, ivec_num
            face(i,face_num) = ivec(i)
          end do
          face_order(face_num) = ivec_num

          ivec_num = 0

        end if
      end if

    end if
!
!  CYLINDER
!
  else if ( s_eqi ( level_name(level), 'Cylinder' ) ) then

    if ( word == '{' ) then

    else if ( word == '}' ) then

      level = nlbrack - nrbrack

    end if
!
!  CYLINDERSENSOR
!
  else if ( s_eqi ( level_name(level), 'CylinderSensor' ) ) then

    if ( word == '{' ) then

    else if ( word == '}' ) then

      level = nlbrack - nrbrack

    end if
!
!  DIRECTIONALLIGHT
!
  else if ( s_eqi ( level_name(level), 'DirectionalLight' ) ) then

    if ( word == '{' ) then

    else if ( word == '}' ) then

      level = nlbrack - nrbrack

    end if
!
!  ELEVATIONGRID
!
  else if ( s_eqi ( level_name(level), 'ElevationGrid' ) ) then

    if ( word == '{' ) then

    else if ( word == '}' ) then

      level = nlbrack - nrbrack

    end if
!
!  EXTRUSION
!
  else if ( s_eqi ( level_name(level), 'Extrusion' ) ) then

    if ( word == '{' ) then

    else if ( word == '}' ) then

      level = nlbrack - nrbrack

    end if
!
!  FOG
!
  else if ( s_eqi ( level_name(level), 'Fog' ) ) then

    if ( word == '{' ) then

    else if ( word == '}' ) then

      level = nlbrack - nrbrack

    end if
!
!  FONTSTYLE
!
  else if ( s_eqi ( level_name(level), 'FontStyle' ) ) then

    if ( word == '{' ) then

    else if ( word == '}' ) then

      level = nlbrack - nrbrack

    end if
!
!  GROUP
!
  else if ( s_eqi ( level_name(level), 'Group' ) ) then

    if ( word == '{' ) then

    else if ( word == '}' ) then

      level = nlbrack - nrbrack

    end if
!
!  IMAGETEXTURE
!
  else if ( s_eqi ( level_name(level), 'ImageTexture' ) ) then

    if ( word == '{' ) then

    else if ( word == '}' ) then

      level = nlbrack - nrbrack

    end if
!
!  INDEXEDFACESET
!
  else if ( s_eqi ( level_name(level), 'IndexedFaceSet' ) ) then

    if ( word == '{' ) then

      material_binding = 'PerVertex'
      cor3_num_old = cor3_num
      face_num_old = face_num
      new_color = 0
      new_color_index = 0

!     default material is not used here, rather index 0 is assigned
!     if needed, default material is assigned after postprocessing

      overall_mat = material_num

    else if ( word == '}' ) then

      level = nlbrack - nrbrack

      imat = overall_mat

      write ( *, * ) 'New_Color = ', new_color
      write ( *, * ) 'New_Color_Index = ', new_color_index
      write ( *, * ) 'Material binding is ' // material_binding(1:9)

      if ( material_binding == 'PerVertex' ) then

        cor3_new = min ( cor3_num, cor3_max ) - min ( cor3_num_old, cor3_max )

        do i = 1, cor3_new
          if ( new_color /= 0 ) then
            if ( new_color_index == 0 ) then
              j = mod ( i-1, new_color )
            else
              k = mod ( i-1, new_color_index ) + 1
              j = ivec(k)
            end if
            imat = overall_mat + j + 1
          end if
          icor3 = cor3_num_old + i
          cor3_material(icor3) = color_mat(imat)
        end do

        new_face = min ( face_num, face_max ) - min ( face_num_old, face_max )

        do i = 1, new_face
          iface = face_num_old + i
          do ivert = 1, face_order(iface)
            node = face(ivert,iface)
            vertex_material(ivert,iface) = cor3_material(node)
          end do
        end do

        ivert = 1
        do i = 1, new_face
          iface = face_num_old + i
          face_material(iface) = vertex_material(ivert,iface)
        end do

      else if ( material_binding == 'PerFace' ) then

        new_face = min ( face_num, face_max ) - min ( face_num_old, face_max )

        do i = 1, new_face

          if ( new_color /= 0 ) then
            if ( new_color_index == 0 ) then
              j = mod ( i-1, new_color )
            else
              k = mod ( i-1, new_color_index ) + 1
              j = ivec(k)
            end if
            imat = overall_mat + j + 1
          end if

          iface = face_num_old + i
          face_material(iface) = color_mat(imat)
        end do 

        do i = 1, new_face
          iface = face_num_old + i
          do ivert = 1, face_order(iface)
            vertex_material(ivert,iface) = face_material(iface)
          end do
        end do

        do i = 1, new_face
          iface = face_num_old + i
          do ivert = 1, face_order(iface)
            node = face(ivert,iface)
            cor3_material(node) = vertex_material(ivert,iface)
          end do
        end do

      else

        write ( *, * ) 'Cannot decide what material binding is...'

      end if

      cor3_num_old = cor3_num
      face_num_old = face_num

    else if ( s_eqi ( word, 'ccw' ) ) then
    else if ( s_eqi ( word, 'color' ) ) then
    else if ( s_eqi ( word, 'colorIndex' ) ) then
    else if ( s_eqi ( word, 'colorPerVertex' ) ) then

      call file_get_next_word ( iunit, word, text, text_num, ierror )

      if ( s_eqi ( word, 'TRUE' ) ) then
        material_binding = 'PerVertex'
      else if ( s_eqi ( word, 'FALSE' ) ) then
        material_binding = 'PerFace'
      end if

    else if ( s_eqi ( word, 'convex' ) ) then
    else if ( s_eqi ( word, 'coord' ) ) then
    else if ( s_eqi ( word, 'coordIndex' ) ) then
    else if ( s_eqi ( word, 'creaseAngle' ) ) then
    else if ( s_eqi ( word, 'normal' ) ) then
    else if ( s_eqi ( word, 'normalIndex' ) ) then
    else if ( s_eqi ( word, 'normalPerVertex' ) ) then
    else if ( s_eqi ( word, 'solid' ) ) then
    else if ( s_eqi ( word, 'texCoord' ) ) then
    else if ( s_eqi ( word, 'texCoordIndex' ) ) then

    end if
!
!  INDEXEDLINESET
!
  else if ( s_eqi ( level_name(level), 'IndexedLineSet' ) ) then

    if ( word == '{' ) then

    else if ( word == '}' ) then

      level = nlbrack - nrbrack

    end if
!
!  INLINE
!
  else if ( s_eqi ( level_name(level), 'Inline' ) ) then

    if ( word == '{' ) then

    else if ( word == '}' ) then

      level = nlbrack - nrbrack

    end if
!
!  LOD
!
  else if ( s_eqi ( level_name(level), 'LOD' ) ) then

    if ( word == '{' ) then

    else if ( word == '}' ) then

      level = nlbrack - nrbrack

    end if
!
!  MATERIAL
!
  else if ( s_eqi ( level_name(level), 'Material' ) ) then

    if ( word == '{' ) then

      r01 = 0.2

      r02 = 0.8
      r03 = 0.8
      r04 = 0.8

      r05 = 0.0
      r06 = 0.0
      r07 = 0.0

      r08 = 0.2

      r09 = 0.0
      r10 = 0.0
      r11 = 0.0

      r12 = 0.0

    else if ( word == '}' ) then

!     The following is commented out, since parsed materials in this part of the file
!     are not used in other parts anyway

!     material_num = material_num + 1
!
!     if ( material_num <= material_max ) then
!       material_rgba(1,material_num) = r02
!       material_rgba(2,material_num) = r03
!       material_rgba(3,material_num) = r04
!       material_rgba(4,material_num) = 1.0 - r12
!       call i_to_s_zero ( material_num, char4 )
!       material_name(material_num) = 'Material_' // char4
!     end if

      level = nlbrack - nrbrack

    else if ( s_eqi ( word, 'ambientIntensity' ) ) then

      call file_get_next_word ( iunit, word, text, text_num, ierror )
      call s_to_r ( word, r01, ierror, lchar )

    else if ( s_eqi ( word, 'diffuseColor' ) ) then

      call file_get_next_word ( iunit, word, text, text_num, ierror )
      call s_to_r ( word, r02, ierror, lchar )
      call file_get_next_word ( iunit, word, text, text_num, ierror )
      call s_to_r ( word, r03, ierror, lchar )
      call file_get_next_word ( iunit, word, text, text_num, ierror )
      call s_to_r ( word, r04, ierror, lchar )

    else if ( s_eqi ( word, 'emissiveColor' ) ) then

      call file_get_next_word ( iunit, word, text, text_num, ierror )
      call s_to_r ( word, r05, ierror, lchar )
      call file_get_next_word ( iunit, word, text, text_num, ierror )
      call s_to_r ( word, r06, ierror, lchar )
      call file_get_next_word ( iunit, word, text, text_num, ierror )
      call s_to_r ( word, r07, ierror, lchar )

    else if ( s_eqi ( word, 'shininess' ) ) then

      call file_get_next_word ( iunit, word, text, text_num, ierror )
      call s_to_r ( word, r08, ierror, lchar )

    else if ( s_eqi ( word, 'specularColor' ) ) then

      call file_get_next_word ( iunit, word, text, text_num, ierror )
      call s_to_r ( word, r09, ierror, lchar )
      call file_get_next_word ( iunit, word, text, text_num, ierror )
      call s_to_r ( word, r10, ierror, lchar )
      call file_get_next_word ( iunit, word, text, text_num, ierror )
      call s_to_r ( word, r11, ierror, lchar )

    else if ( s_eqi ( word, 'transparency' ) ) then

      call file_get_next_word ( iunit, word, text, text_num, ierror )
      call s_to_r ( word, r12, ierror, lchar )

    end if
!
!  MOVIETEXTURE
!
  else if ( s_eqi ( level_name(level), 'MovieTexture' ) ) then

    if ( word == '{' ) then

    else if ( word == '}' ) then

      level = nlbrack - nrbrack

    end if
!
!  NAVIGATIONINFO
!
  else if ( s_eqi ( level_name(level), 'NavigationInfo' ) ) then

    if ( word == '{' ) then

    else if ( word == '}' ) then

      level = nlbrack - nrbrack

    else

    end if
!
!  NORMAL
!
  else if ( s_eqi ( level_name(level), 'Normal' ) ) then

    if ( word == '{' ) then

    else if ( word == '}' ) then

      level = nlbrack - nrbrack

    end if
!
!  NORMALINTERPOLATOR
!
  else if ( s_eqi ( level_name(level), 'NormalInterpolator' ) ) then

    if ( word == '{' ) then

    else if ( word == '}' ) then

      level = nlbrack - nrbrack

    end if
!
!  ORIENTATIONINTERPOLATOR
!
  else if ( s_eqi ( level_name(level), 'OrientationInterpolator' ) ) then

    if ( word == '{' ) then

    else if ( word == '}' ) then

      level = nlbrack - nrbrack

    end if
!
!  PIXELTEXTURE
!
  else if ( s_eqi ( level_name(level), 'PixelTexture' ) ) then

    if ( word == '{' ) then

    else if ( word == '}' ) then

      level = nlbrack - nrbrack

    end if
!
!  PLANESENSOR
!
  else if ( s_eqi ( level_name(level), 'PlaneSensor' ) ) then

    if ( word == '{' ) then

    else if ( word == '}' ) then

      level = nlbrack - nrbrack

    end if
!
!  POINTLIGHT
!
  else if ( s_eqi ( level_name(level), 'PointLight' ) ) then

    if ( word == '{' ) then

    else if ( word == '}' ) then

      level = nlbrack - nrbrack

    end if
!
!  POINTSET
!
  else if ( s_eqi ( level_name(level), 'PointSet' ) ) then

    if ( word == '{' ) then

    else if ( word == '}' ) then

      level = nlbrack - nrbrack

    end if
!
!  POSITIONINTERPOLATOR
!
  else if ( s_eqi ( level_name(level), 'PositionInterpolator' ) ) then

    if ( word == '{' ) then

    else if ( word == '}' ) then

      level = nlbrack - nrbrack

    end if
!
!  PROXIMITYSENSOR
!
  else if ( s_eqi ( level_name(level), 'ProximitySensor' ) ) then

    if ( word == '{' ) then

    else if ( word == '}' ) then

      level = nlbrack - nrbrack

    end if
!
!  SCALARINTERPOLATOR
!
  else if ( s_eqi ( level_name(level), 'ScalarInterpolator' ) ) then

    if ( word == '{' ) then

    else if ( word == '}' ) then

      level = nlbrack - nrbrack

    end if
!
!  SCRIPT
!
  else if ( s_eqi ( level_name(level), 'Script' ) ) then

    if ( word == '{' ) then

    else if ( word == '}' ) then

      level = nlbrack - nrbrack

    end if
!
!  SHAPE
!
  else if ( s_eqi ( level_name(level), 'Shape' ) ) then

    if ( word == '{' ) then

    else if ( word == '}' ) then

      level = nlbrack - nrbrack

    else if ( s_eqi ( word, 'appearance' ) ) then
    else if ( s_eqi ( word, 'geometry' ) ) then

    end if
!
!  SOUND
!
  else if ( s_eqi ( level_name(level), 'Sound' ) ) then

    if ( word == '{' ) then

    else if ( word == '}' ) then

      level = nlbrack - nrbrack

    end if
!
!  SPHERE
!
  else if ( s_eqi ( level_name(level), 'Sphere' ) ) then

    if ( word == '{' ) then

    else if ( word == '}' ) then

      level = nlbrack - nrbrack

    end if
!
!  SPHERESENSOR
!
  else if ( s_eqi ( level_name(level), 'SphereSensor' ) ) then

    if ( word == '{' ) then

    else if ( word == '}' ) then

      level = nlbrack - nrbrack

    end if
!
!  SPOTLIGHT
!
  else if ( s_eqi ( level_name(level), 'SpotLight' ) ) then

    if ( word == '{' ) then

    else if ( word == '}' ) then

      level = nlbrack - nrbrack

    end if
!
!  SWITCH
!
  else if ( s_eqi ( level_name(level), 'Switch' ) ) then

    if ( word == '{' ) then

    else if ( word == '}' ) then

      level = nlbrack - nrbrack

    end if
!
!  TEXT
!
  else if ( s_eqi ( level_name(level), 'Text' ) ) then

    if ( word == '{' ) then

    else if ( word == '}' ) then

      level = nlbrack - nrbrack

    end if
!
!  TEXTURECOORDINATE
!
  else if ( s_eqi ( level_name(level), 'TextureCoordinate' ) ) then

    if ( word == '{' ) then

    else if ( word == '}' ) then

      level = nlbrack - nrbrack

    end if
!
!  TEXTURETRANSFORM
!
  else if ( s_eqi ( level_name(level), 'TextureTransform' ) ) then

    if ( word == '{' ) then

    else if ( word == '}' ) then

      level = nlbrack - nrbrack

    end if
!
!  TIMESENSOR
!
  else if ( s_eqi ( level_name(level), 'TimeSensor' ) ) then

    if ( word == '{' ) then

    else if ( word == '}' ) then

      level = nlbrack - nrbrack

    end if
!
!  TOP
!
  else if ( s_eqi ( level_name(level), 'Top' ) ) then

    if ( word == 'Anchor' ) then
    else if ( word == 'Appearance' ) then
    else if ( word == 'AudioClip' ) then
    else if ( word == 'Background' ) then
    else if ( word == 'Billboard' ) then
    else if ( word == 'Box' ) then
    else if ( word == 'Collision' ) then
    else if ( word == 'Color' ) then
    else if ( word == 'ColorInterpolator' ) then
    else if ( word == 'Cone' ) then
    else if ( word == 'CoordinateInterpolator' ) then
    else if ( word == 'Cylinder' ) then
    else if ( word == 'CylinderSensor' ) then

    else if ( word == 'DirectionalLight' ) then
    else if ( word == 'ElevationGrid' ) then
    else if ( word == 'Extrusion' ) then
    else if ( word == 'Fog' ) then
    else if ( word == 'Fontstyle' ) then
    else if ( word == 'Group' ) then
    else if ( word == 'ImageTexture' ) then
    else if ( word == 'IndexedFaceSet' ) then
    else if ( word == 'IndexedLineSet' ) then
    else if ( word == 'Inline' ) then
    else if ( word == 'LOD' ) then
    else if ( word == 'Material' ) then
    else if ( word == 'MovieTexture' ) then
    else if ( word == 'NavigationInfo' ) then
    else if ( word == 'Normal' ) then
    else if ( word == 'NormalInterpolator' ) then
    else if ( word == 'OrientationInterpolator' ) then
    else if ( word == 'PixelTexture' ) then
    else if ( word == 'PlaneSensor' ) then
    else if ( word == 'PointLight' ) then
    else if ( word == 'PointSet' ) then
    else if ( word == 'PositionInterpolator' ) then
    else if ( word == 'ProximitySensor' ) then
    else if ( word == 'ScalarInterpolator' ) then
    else if ( word == 'Script' ) then
    else if ( word == 'Sound' ) then
    else if ( word == 'Sphere' ) then
    else if ( word == 'SphereSensor' ) then
    else if ( word == 'SpotLight' ) then
    else if ( word == 'Switch' ) then
    else if ( word == 'Text' ) then
    else if ( word == 'TextureCoordinate' ) then
    else if ( word == 'TextureTransform' ) then
    else if ( word == 'TimeSensor' ) then
    else if ( word == 'TouchSensor' ) then
    else if ( word == 'Transform' ) then
    else if ( word == 'Viewpoint' ) then
    else if ( word == 'VisibilitySensor' ) then
    else if ( word == 'WorldInfo' ) then
    end if
!
!  TOUCHSENSOR
!
  else if ( s_eqi ( level_name(level), 'TouchSensor' ) ) then

    if ( word == '{' ) then

    else if ( word == '}' ) then

      level = nlbrack - nrbrack

    end if
!
!  TRANSFORM
!
  else if ( s_eqi ( level_name(level), 'Transform' ) ) then

    if ( word == '{' ) then

      angle = 0.0

      rx = 1.0
      ry = 1.0
      rz = 1.0

      sx = 1.0
      sy = 1.0
      sz = 1.0

      tx = 1.0
      ty = 1.0
      tz = 1.0

    else if ( word == '}' ) then

      level = nlbrack - nrbrack

    else if ( s_eqi ( word, 'rotation' ) ) then

      call file_get_next_word ( iunit, word, text, text_num, ierror )
      call s_to_r ( word, rx, ierror, lchar )
      call file_get_next_word ( iunit, word, text, text_num, ierror )
      call s_to_r ( word, ry, ierror, lchar )
      call file_get_next_word ( iunit, word, text, text_num, ierror )
      call s_to_r ( word, rz, ierror, lchar )
      call file_get_next_word ( iunit, word, text, text_num, ierror )
      call s_to_r ( word, angle, ierror, lchar )

    else if ( s_eqi ( word, 'scale' ) ) then

      call file_get_next_word ( iunit, word, text, text_num, ierror )
      call s_to_r ( word, sx, ierror, lchar )
      call file_get_next_word ( iunit, word, text, text_num, ierror )
      call s_to_r ( word, sy, ierror, lchar )
      call file_get_next_word ( iunit, word, text, text_num, ierror )
      call s_to_r ( word, sz, ierror, lchar )

    else if ( s_eqi ( word, 'translation' ) ) then

      call file_get_next_word ( iunit, word, text, text_num, ierror )
      call s_to_r ( word, tx, ierror, lchar )
      call file_get_next_word ( iunit, word, text, text_num, ierror )
      call s_to_r ( word, ty, ierror, lchar )
      call file_get_next_word ( iunit, word, text, text_num, ierror )

      call s_to_r ( word, tz, ierror, lchar )

    else if ( s_eqi ( word, 'children' ) ) then

      call tmat_init ( transform_matrix )

      call tmat_scale ( transform_matrix, transform_matrix, sx, sy, sz )

      axis(1) = rx
      axis(2) = ry
      axis(3) = rz

      call tmat_rot_vector ( transform_matrix, transform_matrix, angle, axis )

      call tmat_trans ( transform_matrix, transform_matrix, tx, ty, tz )

    end if
!
!  VIEWPOINT
!
  else if ( s_eqi ( level_name(level), 'Viewpoint' ) ) then

    if ( word == '{' ) then

    else if ( word == '}' ) then

      level = nlbrack - nrbrack

    end if
!
!  VISIBILITYSENSOR
!
  else if ( s_eqi ( level_name(level), 'VisibilitySensor' ) ) then

    if ( word == '{' ) then

    else if ( word == '}' ) then

      level = nlbrack - nrbrack

    end if
!
!  WORLDINFO
!
  else if ( s_eqi ( level_name(level), 'WorldInfo' ) ) then

    if ( word == '{' ) then

    else if ( word == '}' ) then

      level = nlbrack - nrbrack

    end if
!
!  Any other word:
!
  else

  end if

  go to 10
!
!  Bad data
!

  bad_num = bad_num + 1

  if ( bad_num <= 10 ) then
    write ( *, * ) ' '
    write ( *, * ) 'VRML_READ - Warning!'
    write ( *, * ) '  Bad data on level ' // trim ( level_name(level) )
    write ( *, '(a)' ) '  Bad word: ' // trim ( word )
    write ( *, * ) '  Line number: ', text_num
    write ( *, '(a)' ) trim ( text )
  else
    write ( *, * ) ' '
    write ( *, * ) 'VRML_READ - Fatal error!'
    write ( *, * ) '  Too many warnings!'
    return
  end if

  go to 10
!
!  End of information in file.
!
50    continue

  ierror = 0

  return
end
subroutine word_nexrd ( line, word, done )
!
!*******************************************************************************
!
!! WORD_NEXRD "reads" words from a string, one at a time.
!
!
!  Special cases:
!
!    The following characters are considered to be a single word,
!    whether surrounded by spaces or not:
!
!      " ( ) { } [ ]
!
!  Modified: 
!
!    16 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) LINE, a string, presumably containing words
!    separated by spaces. 
!
!    Output, character ( len = * ) WORD.  
!
!    If DONE is FALSE, then WORD contains the "next" word read from LINE.  
!    If DONE is TRUE, then WORD is blank, because there was no more to read.
!
!    Input/output, logical DONE.
!
!    On input with a fresh value of LINE, set DONE to TRUE.
!
!    On output, the routine sets DONE:
!      FALSE if another word was read from LINE,
!      TRUE if no more words could be read (LINE is exhausted).
!
  implicit none
  character, parameter :: TAB = char(9)
!
  logical done
  integer ilo
  integer, save :: lenc = 0
  character ( len = * ) line
  integer, save :: next = 1
  character ( len = * ) word
!
!  An input value of DONE = TRUE signals a new line of text to examine.
!
  if ( done ) then

    next = 1
    done = .false.
    lenc = len_trim ( line )

    if ( lenc <= 0 ) then
      done = .true.
      word = ' '
      return
    end if

  end if
!
!  Beginning at index NEXT, search LINE for the next nonblank,
!  which signals the beginning of a word.
!
  ilo = next
 
10    continue
!
!  ...LINE(NEXT:) is blank.  Return with WORD = ' ' and DONE = TRUE.
!
  if ( ilo > lenc ) then
    word = ' '
    done = .true.
    next = lenc + 1
    return
  end if
!
!  If the current character is blank, skip to the next one.
!
  if ( line(ilo:ilo) == ' ' .or. line(ilo:ilo) == TAB ) then
    ilo = ilo + 1
    go to 10
  end if
!
!  ILO is the index of the next nonblank character in the string.
!
!  If this initial nonblank is a special character, 
!  then that's the whole word as far as we're concerned,
!  so return immediately.
!
  if ( line(ilo:ilo) == '"' .or. line(ilo:ilo) == '(' .or. &
       line(ilo:ilo) == ')' .or. line(ilo:ilo) == '{' .or. &
       line(ilo:ilo) == '}' .or. line(ilo:ilo) == '[' .or. &
       line(ilo:ilo) == ']' ) then

    word = line(ilo:ilo)
    next = ilo + 1
    return

  end if     
!
!  Now search for the last contiguous character that is not a
!  blank, TAB, or special character.
!
  next = ilo + 1
 
20    continue
 
  if ( next > lenc ) then
    word = line(ilo:next-1)
    return
  end if
 
  if ( line(next:next) /= ' ' .and. line(next:next) /= TAB .and. &
       line(next:next) /= '"' .and. line(next:next) /= '(' .and. &
       line(next:next) /= ')' .and. line(next:next) /= '{' .and. &
       line(next:next) /= '}' .and. line(next:next) /= '[' .and. &
       line(next:next) /= ']' ) then

    next = next + 1
    go to 20

  end if
 
  word = line(ilo:next-1)
 
  return
end
