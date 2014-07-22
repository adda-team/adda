subroutine ivread_wr(filein_name,face_point,face_normal,face_area,face_num,cor3_num,cor3,face) 
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

  integer, parameter :: cor3_max = 100000
  integer, parameter :: edge_max = 100
  integer, parameter :: face_max = 100000
  integer, parameter :: line_max = 100000
  integer, parameter :: material_max = 200
  integer, parameter :: order_max = 3
  integer, parameter :: texture_max = 10
!
  integer arg_num
  logical byte_swap
  real cor3(3,cor3_max)
  integer cor3_material(cor3_max)
  real cor3_new(3,cor3_max)
  real cor3_normal(3,cor3_max)
  integer cor3_num

  real cor3_tex_uv(2,cor3_max)
  logical debug
  integer edge(4,edge_max)
  integer face(order_max,face_max)
  real face_area(face_max)
  integer face_material(face_max)
  real face_normal(3,face_max)
  integer face_object(face_max)
  integer face_order(face_max)
  integer face_rank(face_max)
  real face_tex_uv(2,face_max)
  integer face_tier(face_max)
  character ( len = 100 ) filein_name
  character ( len = 100 ) fileout_name
!*thw 29.11.00! 
  integer iargc
  integer ierror
  integer line_dex(line_max)
  integer line_material(line_max)
  integer line_prune
  integer list(cor3_max)
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
!*thw 29.11.00! 
  integer na
  integer face_num
  real face_point(3,face_max)
  character ( len = 10 ) filein_type
!
!  Initialize a few program variables.
!
  byte_swap = .true.
  debug = .false.
! filein_name = ' '
  fileout_name = ' '
  ierror = 0
  line_prune = 1
!
!  Get the number of command line arguments.
!
!*thw 29.11.00! 
! arg_num = IARGC ( )
! if ( arg_num >= 2 ) then

    call command_line ( cor3, cor3_material, cor3_normal, &
      cor3_tex_uv, debug, face, face_area, face_material, face_normal, &
      face_order, face_tex_uv, filein_name, fileout_name, ierror, &
      line_dex, line_material, line_prune, material_name, material_rgba, &
      cor3_max, face_max, line_max, material_max, order_max, texture_max, &
      normal_temp, arg_num, object_name, texture_name, texture_temp, &
      transform_matrix, vertex_material, vertex_normal, vertex_tex_uv, &
	  face_num, filein_type, cor3_num)

!  else
!
!    call interact ( byte_swap, cor3, cor3_material, cor3_new, &
!     cor3_normal, cor3_num, cor3_tex_uv, debug, edge, &
!      face, face_area, face_material, face_normal, face_object, &
!      face_order, face_rank, face_tex_uv, face_tier, filein_name, &
!      fileout_name, ierror, line_dex, line_material, line_prune, list, &
!      material_name, material_rgba, cor3_max, edge_max, face_max, &
!      line_max, material_max, order_max, texture_max, normal_temp, &
!      object_name, texture_name, texture_temp, transform_matrix, &
!      vertex_material, vertex_normal, vertex_tex_uv)
!
! end if
!
! if ( ierror /= 0 ) then
!   write ( *, * ) ' '
!    write ( *, * ) 'IVREAD - Error exit.'
!  else
!    write ( *, * ) ' '
!   write ( *, * ) 'IVREAD - Normal exit.'
!  end if
!
	do 100 na=1, face_num
!
!   compute centroid of faces 
!
	face_point(1,na)= (cor3(1,face(1,na))+cor3(1,face(2,na))+cor3(1,face(3,na)))/3.0
	face_point(2,na)= (cor3(2,face(1,na))+cor3(2,face(2,na))+cor3(2,face(3,na)))/3.0
	face_point(3,na)= (cor3(3,face(1,na))+cor3(3,face(2,na))+cor3(3,face(3,na)))/3.0
!
100	continue
!
110 format (I8)
111 format (I8)
112 format (I8,1x,I8,1x,I8)
120	format (F9.6,1x,F9.6,1x,F9.6,1x)
130	format (F9.6)
!*!thw 30.11.00
  return
end
function angle_rad_3d ( x1, y1, z1, x2, y2, z2, x3, y3, z3 )
!
!*******************************************************************************
!
!! ANGLE_RAD_3D returns the angle in radians between two vectors in 3D.
!
!
!  Discussion:
!
!    The routine always computes the SMALLER of the two angles between
!    two vectors.  Thus, if the vectors make an (exterior) angle of 
!    1.5 radians, the (interior) angle of 0.5 radians will be reported.
!
!  Formula:
!
!    X dot Y = Norm(X) * Norm(Y) * Cos ( Angle(X,Y) )
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
!    Input, real X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3, are three points
!    which define the vectors.  The vectors are:
!    ( X1-X2, Y1-Y2, Z1-Z2 ) and ( X3-X2, Y3-Y2, Z3-Z2 ).
!
!    Output, real ANGLE_RAD_3D, the angle between the two vectors, in radians.
!    This value will always be between 0 and PI.  If either vector has 
!    zero length, then the angle is returned as zero.
!
  real angle_rad_3d
  real dot
  real dot0_3d
  real enorm0_3d
  real v1norm
  real v2norm
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
  dot = dot0_3d ( x2, y2, y2, x1, y1, z1, x3, y3, z3 )
  v1norm = enorm0_3d ( x1, y1, z1, x2, y2, z2 )
  v2norm = enorm0_3d ( x3, y3, z3, x2, y2, z2 )
 
  if ( v1norm == 0.0 .or. v2norm == 0.0 ) then
    angle_rad_3d = 0.0
  else
    angle_rad_3d = acos ( dot / ( v1norm * v2norm ) )
  end if
 
  return
end
subroutine ase_read ( cor3, cor3_material, face, face_material, face_normal, face_order, &
  ierror, iunit, material_name, material_rgba, cor3_max, face_max, material_max, &
  order_max, bad_num, cor3_num, face_num, material_num, text_num, vertex_material, &
  vertex_normal )
!
!*******************************************************************************
!
!! ASE_READ reads graphics information from an ASE file.
!
!
!  Problems:
!
!    Processing of the MESH_FACELIST assumes faces are always triangles
!    or quadrilaterals.
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
!    Input/output, integer COR3_MATERIAL(COR3_MAX), the material index of each node.
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
!    Input/output, character ( len = 100 ) MATERIAL_NAME(MATERIAL_MAX), material names.
!
!    Input, integer COR3_MAX, the maximum number of points.
!
!    Input, integer FACE_MAX, the maximum number of faces.
!
!    Input, integer MATERIAL_MAX, the maximum number of materials.
!
!    Input, integer ORDER_MAX, the maximum number of vertices per face.
!
!    Input/output, integer COR3_NUM, the number of points.
!
!    Input/output, integer FACE_NUM, the number of faces.
!
!    Input/output, integer MATERIAL_NUM, the number of materials.
!
!    Output, real TRANSFORM_MATRIX(4,4), the transformation matrix.
!
!    Input/output, integer VERTEX_MAT(ORDER_MAX,FACE_MAX), vertex materials.
!
!    Input/output, real VERTEX_NORMAL(3,ORDER_MAX,FACE_MAX), normals at vertices.  
!
  logical, parameter :: debug = .false.
  integer, parameter :: level_max = 10
  integer, parameter :: OFFSET = 1
!
  integer cor3_max
  integer face_max
  integer material_max
  integer order_max
!
  integer bad_num
  real bval
  character ( len = 4 ) char4
  real cor3(3,cor3_max)
  integer cor3_material(cor3_max)
  integer cor3_num
  integer cor3_num_old
  logical done
  integer face(order_max,face_max)
  integer face_material(face_max)
  real face_normal(3,face_max)
  integer face_num
  integer face_num_old
  integer face_order(face_max)
  real gval
  integer i
  integer ierror
  integer iface
  integer imat
  integer iunit
  integer ivert
  integer iword
  integer k
  integer lchar
  integer level
  character ( len = 256 ) level_name(0:level_max)
  character ( len = 256 ) line
  character ( len = 100 ) material_name(material_max)
  integer material_num
  real material_rgba(4,material_max)
  integer nlbrack
  integer node
  integer nrbrack
  real rgba(4)
  real rval
  logical s_eqi
  real temp
  integer text_num
  real transform_matrix(4,4)
  integer vertex_material(order_max,face_max)
  real vertex_normal(3,order_max,face_max)
  character ( len = 256 ) word
  character ( len = 256 ) word1
  character ( len = 256 ) wordm1
  real x
  real y
  real z
!
  ierror = 0
  level = 0
  level_name(0) = 'Top'
  cor3_num_old = cor3_num
  face_num_old = face_num
  nlbrack = 0
  nrbrack = 0
  call tmat_init ( transform_matrix )
 
  word = ' '
  wordm1 = ' '
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
20    continue
 
  if ( word /= ' ' ) then
    wordm1 = word
  end if
 
  call word_nexrd ( line, word, done )
!
!  If no more words in this line, read a new line.
!
  if ( done ) then
    go to 10
  end if
 
  iword = iword + 1
  if ( iword == 1 ) then
    word1 = word
  end if
!
!  In cases where the word is a left bracket, record the level name,
!  and for right brackets, do a parity check.
!
  if ( word == '{' ) then
 
    nlbrack = nlbrack + 1

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
      do i = 0, level
        write ( *, * ) i, trim ( level_name(i) )
      end do
    end if

  else if ( word == '}' ) then
 
    nrbrack = nrbrack + 1
 
    if ( nlbrack < nrbrack ) then
      write ( *, * ) ' '
      write ( *, * ) 'ASE_READ - Fatal error!'
      write ( *, * ) '  Extraneous right bracket on line ', text_num
      write ( *, '(a)' ) trim ( line )
      write ( *, * ) '  Currently processing field:'
      write ( *, '(a)' ) trim ( level_name(level) )
      ierror = 1
      return
    end if
 
  end if
!
!  *3DSMAX_ASCIIEXPORT  200
!
  if ( word1 == '*3DSMAX_ASCIIEXPORT' ) then
 
    go to 10
!
!  *COMMENT
!
  else if ( word1 == '*COMMENT' ) then
 
    go to 10
!
!  *GEOMOBJECT
!
  else if ( level_name(level) == '*GEOMOBJECT' ) then
 
    if ( word == '{' ) then
      go to 20
    else if ( word == '}' ) then

      level = nlbrack - nrbrack

    else if ( word == '*NODE_NAME' ) then
      go to 10
    else if ( word == '*NODE_TM' ) then
      go to 20
    else if ( word == '*MESH' ) then
      go to 20
    else if ( word == '*PROP_CASTSHADOW' ) then
      go to 10
    else if ( word == '*PROP_MOTIONBLUR' ) then
      go to 10
    else if ( word == '*PROP_RECVSHADOW' ) then
      go to 10
    else
      go to 99
    end if
!
!  *MESH
!
  else if ( level_name(level) == '*MESH' ) then
 
    if ( word == '{' ) then
      go to 20
    else if ( word == '}' ) then

      level = nlbrack - nrbrack

    else if ( word == '*MESH_CFACELIST' ) then
      go to 20
    else if ( word == '*MESH_CVERTLIST' ) then
      go to 20
    else if ( word == '*MESH_FACE_LIST' ) then
      go to 20
    else if ( word == '*MESH_NORMALS' ) then
      go to 20
    else if ( word == '*MESH_NUMCVERTEX' ) then
      go to 10
    else if ( word == '*MESH_NUMCVFACES' ) then
      go to 10
    else if ( word == '*MESH_NUMFACES' ) then
      go to 10
    else if ( word == '*MESH_NUMTVERTEX' ) then
      go to 10
    else if ( word == '*MESH_NUMTVFACES' ) then
      go to 10
    else if ( word == '*MESH_NUMVERTEX' ) then
      go to 10
    else if ( word == '*MESH_TFACELIST' ) then
      go to 20
    else if ( word == '*MESH_TVERTLIST' ) then
      go to 20
    else if ( word == '*MESH_VERTEX_LIST' ) then
      go to 20
    else if ( word == '*TIMEVALUE' ) then
      go to 10
    else
      bad_num = bad_num + 1
      if ( bad_num < 10 ) then
        write ( *, * ) 'Bad data while reading *MESH.'
        write ( *, * ) trim ( line )
      end if
      go to 10
    end if
!
!  *MESH_CFACELIST
!
  else if ( level_name(level) == '*MESH_CFACELIST' ) then
 
    if ( word == '{' ) then
      go to 20
    else if ( word == '}' ) then

      level = nlbrack - nrbrack

    else if ( word == '*MESH_CFACE' ) then
      go to 10
    else
      go to 99
    end if
!
!  *MESH_CVERTLIST
!
!  Mesh vertex indices must be incremented by COR3_NUM_OLD before 
!  being stored in the internal array.
!
  else if ( level_name(level) == '*MESH_CVERTLIST' ) then
 
    if ( word == '{' ) then
 
      go to 20
 
    else if ( word == '}' ) then
 
      level = nlbrack - nrbrack

    else if ( word == '*MESH_VERTCOL' ) then
 
      call word_nexrd ( line, word, done )
      call s_to_i ( word, i, ierror, lchar )
      i = i + cor3_num_old + OFFSET
 
      call word_nexrd ( line, word, done )
      call s_to_r ( word, rval, ierror, lchar )
 
      call word_nexrd ( line, word, done )
      call s_to_r ( word, gval, ierror, lchar )
 
      call word_nexrd ( line, word, done )
      call s_to_r ( word, bval, ierror, lchar )

      rgba(1) = rval
      rgba(2) = gval
      rgba(3) = bval
      rgba(4) = 1.0

      if ( material_num <= 1000 ) then
        call rcol_find ( 4, material_num, material_rgba, rgba, imat )
      else
        imat = 0
      end if

      if ( imat == 0 ) then

        material_num = material_num + 1

        if ( material_num <= material_max ) then

          call i_to_s_zero ( material_num, char4 )
          material_name(material_num) = 'Material_' // char4
          material_rgba(1:4,material_num) = rgba(1:4)
          imat = material_num

        else

          imat = 0

        end if

      end if

      cor3_material(i) = imat
 
    else
      go to 99
    end if
!
!  *MESH_FACE_LIST
!
!  WARNING:
!  The following coding assumes that the faces are always triangles
!  or quadrilaterals, but not higher order.
!
  else if ( level_name(level) == '*MESH_FACE_LIST' ) then
 
    if ( word == '{' ) then
 
      go to 20
 
    else if ( word == '}' ) then
 
      level = nlbrack - nrbrack

    else if ( word == '*MESH_FACE' ) then
 
      face_num = face_num + 1

      if ( face_num <= face_max ) then

        call word_nexrd ( line, word, done )
        call s_to_i ( word, i, ierror, lchar )
        face_material(face_num) = material_num
        face_order(face_num) = 0

        call word_nexrd ( line, word, done )
        call word_nexrd ( line, word, done )
        call s_to_i ( word, i, ierror, lchar )
        face(1,face_num) = i + cor3_num_old + OFFSET
        vertex_material(1,face_num) = material_num
        face_order(face_num) = face_order(face_num) + 1

        call word_nexrd ( line, word, done )
        call word_nexrd ( line, word, done )
        call s_to_i ( word, i, ierror, lchar )
        face(2,face_num) = i + cor3_num_old + OFFSET
        vertex_material(2,face_num) = material_num
        face_order(face_num) = face_order(face_num) + 1

        call word_nexrd ( line, word, done )
        call word_nexrd ( line, word, done )
        call s_to_i ( word, i, ierror, lchar )
        face(3,face_num) = i + cor3_num_old + OFFSET
        vertex_material(3,face_num) = material_num
        face_order(face_num) = face_order(face_num) + 1

        call word_nexrd ( line, word, done )
        if ( s_eqi ( word, 'D:' ) ) then
          call word_nexrd ( line, word, done )
          call s_to_i ( word, i, ierror, lchar )
          face(4,face_num) = i + cor3_num_old + OFFSET
          vertex_material(4,face_num) = material_num
          face_order(face_num) = face_order(face_num) + 1
        end if

      end if
 
      go to 10
 
    else
      go to 99
    end if
!
!  *MESH_NORMALS
!
  else if ( level_name(level) == '*MESH_NORMALS' ) then
 
    if ( word == '{' ) then
      go to 20
    else if ( word == '}' ) then
      level = nlbrack - nrbrack
    else if ( word == '*MESH_FACENORMAL' ) then

      call word_nexrd ( line, word, done )
      call s_to_i ( word, iface, ierror, lchar )
      iface = iface + face_num_old + OFFSET
      ivert = 0

      call word_nexrd ( line, word, done )
      call s_to_r ( word, x, ierror, lchar )

      call word_nexrd ( line, word, done )
      call s_to_r ( word, y, ierror, lchar )

      call word_nexrd ( line, word, done )
      call s_to_r ( word, z, ierror, lchar )

      face_normal(1,iface) = x
      face_normal(2,iface) = y
      face_normal(3,iface) = z

      go to 10

    else if ( word == '*MESH_VERTEXNORMAL' ) then

      call word_nexrd ( line, word, done )
      call s_to_i ( word, node, ierror, lchar )
      ivert = ivert + 1

      call word_nexrd ( line, word, done )
      call s_to_r ( word, x, ierror, lchar )

      call word_nexrd ( line, word, done )
      call s_to_r ( word, y, ierror, lchar )

      call word_nexrd ( line, word, done )
      call s_to_r ( word, z, ierror, lchar )

      vertex_normal(1,ivert,iface) = x
      vertex_normal(2,ivert,iface) = y
      vertex_normal(3,ivert,iface) = z

      go to 10

    else
      go to 99
    end if
!
!  *MESH_TFACELIST
!
  else if ( level_name(level) == '*MESH_TFACELIST' ) then
 
    if ( word == '{' ) then
      go to 20
    else if ( word == '}' ) then

      level = nlbrack - nrbrack

    else if ( word1 == '*MESH_TFACE' ) then
      go to 10
    else
      go to 99
    end if
!
!  *MESH_TVERTLIST
!
  else if ( level_name(level) == '*MESH_TVERTLIST' ) then
 
    if ( word == '{' ) then
      go to 20
    else if ( word == '}' ) then
      level = nlbrack - nrbrack
    else if ( word1 == '*MESH_TVERT' ) then
      go to 10
    else
      go to 99
    end if
!
!  *MESH_VERTEX_LIST
!
  else if ( level_name(level) == '*MESH_VERTEX_LIST' ) then
 
    if ( word == '{' ) then 
      cor3_num_old = cor3_num
      go to 20
     else if ( word == '}' ) then
      level = nlbrack - nrbrack
    else if ( word1 == '*MESH_VERTEX' ) then
 
      call word_nexrd ( line, word, done ) 
      call s_to_i ( word, i, ierror, lchar )
      call word_nexrd ( line, word, done )
      call s_to_r ( word, x, ierror, lchar )
      call word_nexrd ( line, word, done )
      call s_to_r ( word, y, ierror, lchar )
      call word_nexrd ( line, word, done )
      call s_to_r ( word, z, ierror, lchar )
 
      i = i + cor3_num_old + OFFSET
      cor3_num = max ( cor3_num, i )
      if ( i <= cor3_max ) then
        cor3(1,i) =  transform_matrix(1,1) * x + transform_matrix(1,2) * y &
                   + transform_matrix(3,1) * z + transform_matrix(4,1)

        cor3(2,i) =  transform_matrix(2,1) * x + transform_matrix(2,2) * y &
                   + transform_matrix(2,3) * z + transform_matrix(2,4)

        cor3(3,i) =  transform_matrix(3,1) * x + transform_matrix(3,2) * y &
                   + transform_matrix(3,3) * z + transform_matrix(3,4)
      end if

      go to 10

    else
      go to 99
    end if
!
!  *NODE_TM
!
!  Each node should start out with a default transformation matrix.
!
  else if ( level_name(level) == '*NODE_TM' ) then
 
    if ( word == '{' ) then

      call tmat_init ( transform_matrix )
      go to 20

    else if ( word == '}' ) then
      level = nlbrack - nrbrack
    else if ( word == '*INHERIT_POS' ) then
      go to 10
    else if ( word == '*INHERIT_ROT' ) then
      go to 10
    else if ( word == '*INHERIT_SCL' ) then
      go to 10
    else if ( word == '*NODE_NAME' ) then
      go to 10
    else if ( word == '*TM_POS' ) then
      go to 10
    else if ( word == '*TM_ROTANGLE' ) then
      go to 10
    else if ( word == '*TM_ROTAXIS' ) then
      go to 10
    else if ( word == '*TM_ROW0' ) then
      call word_nexrd ( line, word, done )
      call s_to_r ( word, temp, ierror, lchar )
      transform_matrix(1,1) = temp
      call word_nexrd ( line, word, done )
      call s_to_r ( word, temp, ierror, lchar )
      transform_matrix(2,1) = temp
      call word_nexrd ( line, word, done )
      call s_to_r ( word, temp, ierror, lchar )
      transform_matrix(3,1) = temp
      go to 10
    else if ( word == '*TM_ROW1' ) then
      call word_nexrd ( line, word, done )
      call s_to_r ( word, temp, ierror, lchar )
      transform_matrix(1,2) = temp
      call word_nexrd ( line, word, done )
      call s_to_r ( word, temp, ierror, lchar )
      transform_matrix(2,2) = temp
      call word_nexrd ( line, word, done )
      call s_to_r ( word, temp, ierror, lchar )
      transform_matrix(3,2) = temp
      go to 10
    else if ( word == '*TM_ROW2' ) then
      call word_nexrd ( line, word, done )
      call s_to_r ( word, temp, ierror, lchar )
      transform_matrix(1,3) = temp
      call word_nexrd ( line, word, done )
      call s_to_r ( word, temp, ierror, lchar )
      transform_matrix(2,3) = temp
      call word_nexrd ( line, word, done )
      call s_to_r ( word, temp, ierror, lchar )
      transform_matrix(3,3) = temp
      go to 10
    else if ( word == '*TM_ROW3' ) then
      call word_nexrd ( line, word, done )
      call s_to_r ( word, temp, ierror, lchar )
      transform_matrix(1,4) = temp
      call word_nexrd ( line, word, done )
      call s_to_r ( word, temp, ierror, lchar )
      transform_matrix(2,4) = temp
      call word_nexrd ( line, word, done )
      call s_to_r ( word, temp, ierror, lchar )
      transform_matrix(3,4) = temp
      go to 10
    else if ( word == '*TM_SCALE' ) then
      go to 10
    else if ( word == '*TM_SCALEAXIS' ) then
      go to 10
    else if ( word == '*TM_SCALEAXISANG' ) then
      go to 10
    else
      go to 99
    end if
!
!  *SCENE
!
  else if ( level_name(level) == '*SCENE' ) then

    if ( word == '{' ) then
      go to 20
    else if ( word == '}' ) then
      level = nlbrack - nrbrack
    else if ( word == '*SCENE_AMBIENT_STATIC' ) then
      go to 10
    else if ( word == '*SCENE_BACKGROUND_STATIC' ) then
      go to 10
    else if ( word == '*SCENE_FILENAME' ) then
      go to 10
    else if ( word == '*SCENE_FIRSTFRAME' ) then
      go to 10
    else if ( word == '*SCENE_FRAMESPEED' ) then
      go to 10
    else if ( word == '*SCENE_LASTFRAME' ) then
      go to 10
    else if ( word == '*SCENE_TICKSPERFRAME' ) then
      go to 10
    else
      go to 99
    end if
 
  end if
 
  go to 20
!
!  Bad data
!
99    continue

  bad_num = bad_num + 1

  if ( bad_num <= 10 ) then
    write ( *, * ) ' '
    write ( *, * ) 'ASE_READ - Warning!'
    write ( *, * ) '  Bad data on level ' // trim ( level_name(level) )
    write ( *, '(a)' ) '  Bad word: ' // trim ( word )
    write ( *, * ) '  Line number: ', text_num
    write ( *, '(a)' ) trim ( line )
  end if

  go to 10
!
!  End of information in file.
!
30    continue

  return
end
subroutine ase_write ( cor3, face, face_normal, face_order, filein_name, &
  fileout_name, iunit, cor3_max, face_max, order_max, cor3_num, face_num, &
  vertex_normal )
!
!*******************************************************************************
!
!! ASE_WRITE writes graphics information to an ASE file.
!
!
!  Modified:
!
!    22 April 1999
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
!    Input, real FACE_NORMAL(3,FACE_MAX), the normal vector at each face.
!
!    Input, integer FACE_ORDER(FACE_MAX), the number of vertices per face.
!
!    Input, character ( len = 100 ) FILEIN_NAME, the name of the input file.
!
!    Input, character ( len = 100 ) FILEOUT_NAME, the name of the output file.
!
!    Input, integer IUNIT, the FORTRAN unit to which output is written.
!
!    Input, integer COR3_MAX, the maximum number of points.
!
!    Input, integer FACE_MAX, the maximum number of faces.
!
!    Input, integer ORDER_MAX, the maximum number of vertices per face.
!
!    Input, integer COR3_NUM, the number of points.
!
!    Input, integer FACE_NUM, the number of faces.
!
!    Input, real VERTEX_NORMAL(3,ORDER_MAX,FACE_MAX), normals at vertices.  
!
  integer, parameter :: OFFSET = 1
!
  integer cor3_max
  integer face_max
  integer order_max
!
  character ( len = 10 ) chrtmp
  real cor3(3,cor3_max)
  integer cor3_num
  integer face(order_max,face_max)
  real face_normal(3,face_max)
  integer face_num
  integer face_order(face_max)
  character ( len = 100 ) filein_name
  character ( len = 100 ) fileout_name
  integer i
  integer i1
  integer i2
  integer i3
  integer i4
  integer iface
  integer iunit
  integer ivert
  integer j
  integer text_num
  character ( len = 200 ) text
  real vertex_normal(3,order_max,face_max)
  real x
  real y
  real z
!
  text_num = 0
!
!  Write the header.
!
  write ( iunit, '(a)' ) '*3DSMAX_ASCIIEXPORT 200'
  write ( iunit, '(a)' ) '*COMMENT "' // trim ( fileout_name ) // &
    ', created by IVREAD."'
  write ( iunit, '(a)' ) '*COMMENT "Original data in ' // trim ( filein_name ) &
    // '."'

  text_num = text_num + 3
!
!  Write the scene block.
!
  write ( iunit, '(a)' ) '*SCENE {'
  write ( iunit, '(a)' ) '  *SCENE_FILENAME ""'
  write ( iunit, '(a)' ) '  *SCENE_FIRSTFRAME 0'
  write ( iunit, '(a)' ) '  *SCENE_LASTFRAME 100'
  write ( iunit, '(a)' ) '  *SCENE_FRAMESPEED 30'
  write ( iunit, '(a)' ) '  *SCENE_TICKSPERFRAME 160'
  write ( iunit, '(a)' ) '  *SCENE_BACKGROUND_STATIC 0.0000 0.0000 0.0000'
  write ( iunit, '(a)' ) '  *SCENE_AMBIENT_STATIC 0.0431 0.0431 0.0431'
  write ( iunit, '(a)' ) '}'

  text_num = text_num + 9
!
!  Begin the big geometry block.
!
  write ( iunit, '(a)' ) '*GEOMOBJECT {'
  write ( iunit, '(a)' ) '  *NODE_NAME "Object01"'

  text_num = text_num + 2
!
!  Sub block NODE_TM:
!
  write ( iunit, '(a)' ) '  *NODE_TM {'
  write ( iunit, '(a)' ) '    *NODE_NAME "Object01"'
  write ( iunit, '(a)' ) '    *INHERIT_POS 0 0 0'
  write ( iunit, '(a)' ) '    *INHERIT_ROT 0 0 0'
  write ( iunit, '(a)' ) '    *INHERIT_SCL 0 0 0'
  write ( iunit, '(a)' ) '    *TM_ROW0 1.0000 0.0000 0.0000'
  write ( iunit, '(a)' ) '    *TM_ROW1 0.0000 1.0000 0.0000'
  write ( iunit, '(a)' ) '    *TM_ROW2 0.0000 0.0000 1.0000'
  write ( iunit, '(a)' ) '    *TM_ROW3 0.0000 0.0000 0.0000'
  write ( iunit, '(a)' ) '    *TM_POS 0.0000 0.0000 0.0000'
  write ( iunit, '(a)' ) '    *TM_ROTAXIS 0.0000 0.0000 0.0000'
  write ( iunit, '(a)' ) '    *TM_ROTANGLE 0.0000'
  write ( iunit, '(a)' ) '    *TM_SCALE 1.0000 1.0000 1.0000'
  write ( iunit, '(a)' ) '    *TM_SCALEAXIS 0.0000 0.0000 0.0000'
  write ( iunit, '(a)' ) '    *TM_SCALEAXISANG 0.0000'
  write ( iunit, '(a)' ) '  }'

  text_num = text_num + 16
!
!  Sub block MESH:
!    Items
!
  write ( iunit, '(a)' ) '  *MESH {'
  write ( iunit, '(a)' ) '    *TIMEVALUE 0'
  write ( chrtmp, '(i8)' ) cor3_num
  write ( iunit, '(a)' ) '    *MESH_NUMVERTEX ' // trim ( chrtmp )
  write ( chrtmp, '(i8)' ) face_num
  write ( iunit, '(a)' ) '    *MESH_NUMFACES ' // trim ( chrtmp )

  text_num = text_num + 4
!
!  Sub sub block MESH_VERTEX_LIST
!
  write ( iunit, '(a)' ) '    *MESH_VERTEX_LIST {'

  do j = 1, cor3_num
    write ( text, '(a,i8,3g12.4)' ) '*MESH_VERTEX ', j - OFFSET, cor3(1:3,j)
    call s_blanks_delete ( text )
    write ( iunit, '(6x,a)' ) trim ( text )
  end do

  write ( iunit, '(a)' ) '    }'

  text_num = text_num + cor3_num + 2
!
!  Sub sub block MESH_FACE_LIST
!    Items MESH_FACE
!
!  ???  What do you do when the face has more than 4 vertices?
!
  write ( iunit, '(a)' ) '    *MESH_FACE_LIST {'

  do iface = 1, face_num

    i1 = face(1,iface) - OFFSET
    i2 = face(2,iface) - OFFSET
    i3 = face(3,iface) - OFFSET

    if ( face_order(iface) == 3 ) then

      write ( text, '(a,i8,a,i8,a,i8,a,i8,a)' )  '*MESH_FACE ', &
        iface - OFFSET, ': A: ', i1, ' B: ', i2, ' C: ', i3, &
        ' AB: 1 BC: 1 CA: 1 *MESH_SMOOTHING *MESH_MTLID 1'

      call s_blanks_delete ( text )
      write ( iunit, '(6x,a)' ) trim ( text )

    else if ( face_order(iface) == 4 ) then

      i4 = face(4,iface) - OFFSET
      write ( text, '(a,i8,a,i8,a,i8,a,i8,a,i8,a)' ) '*MESH_FACE ', &
        iface - OFFSET, ': A: ', i1, ' B: ', i2, ' C: ', i3, ' D: ', i4, &
        ' AB: 1 BC: 1 CD: 1 DA: 1 *MESH_SMOOTHING *MESH_MTLID 1'
      call s_blanks_delete ( text )
      write ( iunit, '(6x,a)' ) trim ( text )

    end if

  end do

  write ( iunit, '(a)' ) '    }'

  text_num = text_num + face_num + 2
!
!  Item MESH_NUMTVERTEX
!
  write ( iunit, '(a)' ) '    *MESH_NUMTVERTEX 0'
  text_num = text_num + 1
!
!  Item NUMCVERTEX
!
  write ( iunit, '(a)' ) '    *MESH_NUMCVERTEX 0'
  text_num = text_num + 1
!
!  Sub block MESH_NORMALS
!    Items MESH_FACENORMAL, MESH_VERTEXNORMAL (repeated)
!
  write ( iunit, '(a)' ) '    *MESH_NORMALS {'
  text_num = text_num + 1

  do iface = 1, face_num

    x = face_normal(1,iface)
    y = face_normal(2,iface)
    z = face_normal(3,iface)

    write ( text, '(a,i8,3g12.4)' ) '*MESH_FACENORMAL ', iface-OFFSET, x, y, z
    call s_blanks_delete ( text )

    write ( iunit, '(6x,a)' ) trim ( text )
    text_num = text_num + 1

    do ivert = 1, face_order(iface)

      x = vertex_normal(1,ivert,iface)
      y = vertex_normal(2,ivert,iface)
      z = vertex_normal(3,ivert,iface)

      write ( text, '(a,i8,3g12.4)' ) '*MESH_VERTEXNORMAL ', ivert-OFFSET, &
        x, y, z 
      call s_blanks_delete ( text )
      write ( iunit, '(6x,a)' ) trim ( text )
      text_num = text_num + 1

    end do

  end do

  write ( iunit, '(a)' ) '    }'
  text_num = text_num + 1
!
!  Close the MESH object.
!
  write ( iunit, '(a)' ) '  }'
!
!  A few closing parameters.
!
  write ( iunit, '(a)' ) '  *PROP_MOTIONBLUR 0'
  write ( iunit, '(a)' ) '  *PROP_CASTSHADOW 1'
  write ( iunit, '(a)' ) '  *PROP_RECVSHADOW 1'
!
!  Close the GEOM object.
!
  write ( iunit, '(a)' ) '}'

  text_num = text_num + 5
!
!  Report.
!
  write ( *, * ) ' '
  write ( *, * ) 'ASE_WRITE - Wrote ', text_num, ' text lines.'

  return
end
subroutine byu_read ( cor3, face, face_order, ierror, iunit, cor3_max, &
  face_max, order_max, cor3_num, face_num )
!
!*******************************************************************************
!
!! BYU_READ reads graphics data from a Movie.BYU surface geometry file.
!
!
!  Discussion:
!
!    This code will certainly read a BYU file created by BYU_WRITE, but
!    it will not handle more general files.  In particular, an object
!    can have several parts, the coordinate data can be grouped so
!    that there are 2 sets of (x,y,z) data per line, and so on.
!
!  Example:
!
!          0       8       6       0
!          0       6
!    0.00000E+00 0.00000E+00 0.00000E+00
!    1.00000E+00 0.00000E+00 0.00000E+00
!    1.00000E+00 2.00000E+00 0.00000E+00
!    0.00000E+00 2.00000E+00 0.00000E+00
!    0.00000E+00 0.00000E+00 1.00000E+00
!    1.00000E+00 0.00000E+00 1.00000E+00
!    1.00000E+00 2.00000E+00 1.00000E+00
!    0.00000E+00 2.00000E+00 1.00000E+00
!          4       3       2      -1
!          5       6       7      -8
!          1       5       8      -4
!          4       8       7      -3
!          3       7       6      -2
!          2       6       5      -1
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
!    Input/output, real COR3(3,COR3_MAX), the coordinates of points.
!
!    Input/output, integer FACE(ORDER_MAX,FACE_MAX), the nodes making faces.
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
  integer cor3_max
  integer face_max
  integer order_max
!
  real cor3(3,cor3_max)
  integer cor3_num
  integer cor3_num_new
  logical done
  integer face(order_max,face_max)
  integer face_num
  integer face_num_new
  integer face_order(face_max)
  integer i
  integer ierror
  integer iface
  integer ios
  integer irow
  integer iunit
  integer ival
  integer ivert
  integer j
  integer np
  integer np1
  integer np2
  character ( len = 256 ) text
  integer text_num
!
  ierror = 0
  text_num = 0

  read ( iunit, *, iostat = ios ) np, cor3_num_new, face_num_new, irow

  if ( ios /= 0 ) then
    ierror = ios
    write ( *, * ) ' '
    write ( *, * ) 'BYU_READ - Fatal error!'
    write ( *, * ) '  Unexpected end of file.'
    return
  end if

  text_num = text_num + 1

  read ( iunit, *, iostat = ios ) np1, np2

  if ( ios /= 0 ) then
    ierror = ios
    write ( *, * ) ' '
    write ( *, * ) 'BYU_READ - Fatal error!'
    write ( *, * ) '  Unexpected end of file.'
    return
  end if

  text_num = text_num + 1

  do j = cor3_num + 1, cor3_num + cor3_num_new

    read ( iunit, *, iostat = ios ) cor3(1:3,j)

    if ( ios /= 0 ) then
      ierror = ios
      write ( *, * ) ' '
      write ( *, * ) 'BYU_READ - Fatal error!'
      write ( *, * ) '  Unexpected end of file.'
      return
    end if

    text_num = text_num + 1

  end do

  do iface = face_num + 1, face_num + face_num_new

    read ( iunit, '(a)', iostat = ios ) text
    text_num = text_num + 1

    ivert = 0
    done = .true.

10  continue

    call intnex ( text, ival, done )

    if ( .not. done ) then

      ivert = ivert + 1
      face(ivert,iface) = abs ( ival ) + cor3_num

      if ( ival > 0 ) then
        go to 10
      end if

    end if

    face_order(iface) = ivert

  end do

  cor3_num = cor3_num + cor3_num_new
  face_num = face_num + face_num_new
!
!  Report.
!
  write ( *, * ) ' '
  write ( *, * ) 'BYU_READ - Read ', text_num, ' text lines.'

  return
end
subroutine byu_write ( cor3, face, face_order, iunit, cor3_max, face_max, &
  order_max, cor3_num, face_num )
!
!*******************************************************************************
!
!! BYU_WRITE writes out the graphics data as a Movie.BYU surface geometry file.
!
!
!  Example:
!
!          0       8       6       0
!          0       6
!    0.00000E+00 0.00000E+00 0.00000E+00
!    1.00000E+00 0.00000E+00 0.00000E+00
!    1.00000E+00 2.00000E+00 0.00000E+00
!    0.00000E+00 2.00000E+00 0.00000E+00
!    0.00000E+00 0.00000E+00 1.00000E+00
!    1.00000E+00 0.00000E+00 1.00000E+00
!    1.00000E+00 2.00000E+00 1.00000E+00
!    0.00000E+00 2.00000E+00 1.00000E+00
!          4       3       2      -1
!          5       6       7      -8
!          1       5       8      -4
!          4       8       7      -3
!          3       7       6      -2
!          2       6       5      -1
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
!    Input, real COR3(3,COR3_MAX), the coordinates of points.
!
!    Input, integer FACE(ORDER_MAX,FACE_MAX), the nodes making faces.
!
!    Input, integer FACE_ORDER(FACE_MAX), the number of vertices per face.
!
!    Input, integer IUNIT, the FORTRAN unit from which data is read.
!
!    Input, integer COR3_MAX, the maximum number of points.
!
!    Input, integer FACE_MAX, the maximum number of faces.
!
!    Input, integer ORDER_MAX, the maximum number of vertices per face.
!
!    Input, integer COR3_NUM, the number of points.
!
!    Input, integer FACE_NUM, the number of faces.
!
  integer cor3_max
  integer face_max
  integer order_max
!
  real cor3(3,cor3_max)
  integer cor3_num
  integer face(order_max,face_max)
  integer face_num
  integer face_order(face_max)
  integer i
  integer iface
  integer ihi
  integer irow
  integer iunit
  integer ivert
  integer j
  integer jp(8)
  integer nedge
  integer np
  integer text_num
!
  text_num = 0
!
!  NEDGE is the total number of edges.
!
  nedge = 0
  do iface = 1, face_num
    nedge = nedge + face_order(iface)
  end do

  irow = 0
  np = 0
  write ( iunit, '(10i8)' ) np, cor3_num, face_num, irow
  text_num = text_num + 1

  write ( iunit, '(10i8)' ) 0, face_num
  text_num = text_num + 1

  do j = 1, cor3_num
    write ( iunit, '(1p6e12.5)' ) cor3(1:3,j)
    text_num = text_num + 1
  end do
!
!  It takes a little mangling in order to print out all the edges in a
!  single list, one face at a time, with the last node in each face
!  negative, and written in groups of 8.
!
  ihi = 0

  do iface = 1, face_num

    do ivert = 1, face_order(iface)

      ihi = ihi + 1
      jp(ihi) = face(ivert,iface)
      if ( ivert == face_order(iface) ) then
        jp(ihi) = - jp(ihi)
      end if

      if ( ihi == 8 .or. ivert == face_order(iface) ) then
        write ( iunit, '(10i8)' ) jp(1:ihi)
        text_num = text_num + 1
        ihi = 0
      end if

    end do
  end do
!
!  Report.
!
  write ( *, * ) ' '
  write ( *, * ) 'BYU_WRITE - Wrote ', text_num, ' text lines.'

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
subroutine command_line ( cor3, cor3_material, cor3_normal, &
  cor3_tex_uv, debug, face, face_area, face_material, face_normal, &
  face_order, face_tex_uv, filein_name, fileout_name, ierror, &
  line_dex, line_material, line_prune, material_name, material_rgba, &
  cor3_max, face_max, line_max, material_max, order_max, texture_max, &
  normal_temp, arg_num, object_name, texture_name, &
  texture_temp, transform_matrix, vertex_material, vertex_normal, &
  vertex_tex_uv, face_num, filein_type, cor3_num)
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
!    Workspace, character ( len = 100 ) FILEOUT_NAME, the name of the output file.
!
!    Output, integer IERROR, error flag.
!
!    Workspace, integer LINE_DEX(LINE_MAX), nodes forming a line, terminated by -1.
!
!    Workspace, integer LINE_MAT(LINE_MAX), material index for line.
!
!    Input, integer LINE_PRUNE, pruning option.
!    0, no pruning, draw every line.
!    nonzero, prune.  Only draw the line from node I to node J if I < J.
!    This should cut down on repeated drawing of lines in the common
!    case of a face mesh where each line is drawn twice, once with positive
!    and once with negative orientation.  In other cases, pruning
!    may omit some lines that only occur with negative orientation.
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
!
!    Input, integer ARG_NUM, the number of command-line arguments.
!
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
!*thw 29.11.00! 
! use DFLIB
! use DFPORT 

  integer cor3_max
  integer face_max
  integer line_max
  integer material_max
  integer order_max
  integer texture_max
!
  integer arg_num
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
  character ( len = 100 ) fileout_name
  character ( len = 10 ) fileout_type
  integer group_num
  integer i
  integer iarg
  integer icor3
  integer ierror
  integer iface
  integer ilen
  integer iseed
  integer ivert
  integer line_dex(line_max)
  integer line_material(line_max)
  integer line_num
  integer line_prune
  character ( len = 100 ) material_name(material_max)
  integer material_num
  real material_rgba(4,material_max)
  real normal_temp(3,order_max*face_max)
  character ( len = 100 ) object_name
  integer object_num
  logical reverse_faces
  logical reverse_normals
  logical s_eqi
  character ( len = 100 ) texture_name(texture_max)
  integer texture_num
  real texture_temp(2,order_max*face_max)
  real transform_matrix(4,4)
  integer vertex_material(order_max,face_max)
  real vertex_normal(3,order_max,face_max)
  real vertex_tex_uv(2,order_max,face_max)
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
!  Sort out the command line arguments.
!
!*thw 29.11.00!
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
    texture_temp, vertex_material, vertex_normal, vertex_tex_uv )
!
!  Check the output file name.
!
!*thw 29.11.00!
! call outfile ( filein_name, fileout_name, ierror, fileout_type )
!
! if ( ierror /= 0 ) then
!    write ( *, * ) ' '
!    write ( *, * ) 'COMMAND_LINE - Fatal error!'
!    write ( *, * ) '  Improper output file name.'
!    return
!  end if
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
!
!  Write the output file.
!
!*! thw 1.12.00
!
!  call data_write ( cor3, cor3_material, cor3_normal, cor3_tex_uv, debug, face, &
!   face_material, face_normal, face_order, face_tex_uv, filein_name, &
!    fileout_name, fileout_type, ierror, line_dex, line_material, line_prune, &
!    material_name, material_rgba, cor3_max, face_max, line_max, material_max, &
!    order_max, texture_max, cor3_num, face_num, line_num, material_num, texture_num, &
!    object_name, texture_name, vertex_material, vertex_normal, vertex_tex_uv )

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
subroutine cor3_range ( cor3_max, cor3_num, cor3 )
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
      texture_name(i) = 'Material_' // char4
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
  integer iface
  integer iseed
  integer ivert
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
    call i_to_s_zero ( i, char4 )
    material_name(j) = 'Material' // char4
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
  vertex_tex_uv )
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
  if ( s_eqi ( filein_type, 'ASE' ) ) then
 
    call ase_read ( cor3, cor3_material, face, face_material, face_normal, face_order, &
      ierror, iunit, material_name, material_rgba, cor3_max, face_max, &
      material_max, order_max, bad_num, cor3_num, face_num, material_num, text_num, &
      vertex_material, vertex_normal )

    call node_to_vertex_material ( cor3_material, face, face_order, cor3_max, face_max, &
      order_max, face_num, vertex_material )

    face_material(1:face_num) = vertex_material(1,1:face_num)

  else if ( s_eqi ( filein_type, 'BYU' ) ) then
 
    call byu_read ( cor3, face, face_order, ierror, iunit, cor3_max, &
      face_max, order_max, cor3_num, face_num )

  else if ( s_eqi ( filein_type, 'DXF' ) ) then
 
    call dxf_read ( cor3, face, face_material, face_order, ierror, iunit, &
      line_dex, line_material, cor3_max, face_max, line_max, order_max, bad_num, &
      cor3_num, dup_num, face_num, line_num, material_num, text_num )

  else if ( s_eqi ( filein_type, 'HRC' ) ) then
 
    call hrc_read ( cor3, face, face_material, face_order, ierror, iunit, line_dex, &
      line_material, material_name, material_rgba, cor3_max, face_max, line_max, &
      material_max, order_max, texture_max, bad_num, cor3_num, dup_num, face_num, &
      line_num, material_num, texture_num, text_num, texture_name, vertex_material, &
      vertex_normal, vertex_tex_uv )

  else if ( s_eqi ( filein_type, 'IV' ) ) then
 
    call iv_read ( cor3, debug, face, face_order, ierror, iunit, line_dex, &
      line_material, cor3_max, face_max, line_max, order_max, texture_max, &
      normal_temp, bad_num, color_num, cor3_num, face_num, line_num, material_num, &
      texture_num, text_num, texture_name, texture_temp, vertex_material, &
      vertex_normal, vertex_tex_uv )

  else if ( s_eqi ( filein_type, 'OBJ' ) ) then

    call obj_read ( cor3, face, face_material, face_order, ierror, iseed, iunit, &
      line_dex, line_material, material_name, material_rgba, cor3_max, face_max, &
      line_max, material_max, order_max, normal_temp, bad_num, cor3_num, face_num, &
      group_num, line_num, material_num, object_num, text_num, vertex_material, &
      vertex_normal )

  else if ( s_eqi ( filein_type, 'OOGL' ) ) then
 
    call oogl_read ( cor3, cor3_material, cor3_normal, &
      face, face_area, face_material, face_normal, face_order, ierror, &
      iunit, material_name, material_rgba, cor3_max, face_max, &
      material_max, order_max, cor3_num, face_num, material_num, text_num, &
      ncol_oogl, nrow_oogl, vertex_material, vertex_normal )

  else if ( s_eqi ( filein_type, 'SMF' ) ) then

    call smf_read ( cor3, cor3_material, cor3_normal, cor3_tex_uv, &
      debug, face, face_material, face_normal, face_order, face_tex_uv, &
      ierror, iunit, material_name, material_rgba, cor3_max, &
      face_max, material_max, order_max, texture_max, bad_num, &
      cor3_num, face_num, group_num, material_num, texture_num, text_num, &
      texture_name, vertex_material )

  else if ( s_eqi ( filein_type, 'STL' ) .or. &
            s_eqi ( filein_type, 'STLA' ) ) then

    call stla_read ( cor3, face, face_material, face_normal, face_order, ierror, &
      iunit, cor3_max, face_max, order_max, bad_num, cor3_num, dup_num, &
      face_num, material_num, object_num, text_num, vertex_material, vertex_normal )

  else if ( s_eqi ( filein_type, 'TRI' ) .or. &
            s_eqi ( filein_type, 'TRIA' ) ) then

    call tria_read ( cor3, face, face_material, face_order, ierror, iunit, &
      cor3_max, face_max, order_max, cor3_num, dup_num, face_num, text_num, &
      vertex_material, vertex_normal )

  else if ( s_eqi ( filein_type, 'VLA' ) ) then
 
    call vla_read ( cor3, ierror, iunit, line_dex, line_material, cor3_max, &
      line_max, bad_num, cor3_num, dup_num, line_num, material_num, text_num )

  else if ( s_eqi ( filein_type, 'WRL' ) ) then

    call vrml_read ( cor3, cor3_material, face, face_material, face_order, ierror, &
      iunit, line_dex, line_material, material_name, material_rgba, cor3_max, &
      face_max, line_max, material_max, order_max, texture_max, bad_num, cor3_num, &
      face_num, line_num, material_num, texture_num, text_num, texture_name, &
      vertex_material, vertex_normal, vertex_tex_uv )

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
  call cor3_range ( cor3_max, cor3_num, cor3 )

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
subroutine data_write ( cor3, cor3_material, cor3_normal, cor3_tex_uv, &
  debug, face, face_material, face_normal, face_order, face_tex_uv, &
  filein_name, fileout_name, fileout_type, ierror, line_dex, &
  line_material, line_prune, material_name, material_rgba, cor3_max, &
  face_max, line_max, material_max, order_max, texture_max, cor3_num, &
  face_num, line_num, material_num, texture_num, object_name, &
  texture_name, vertex_material, vertex_normal, vertex_tex_uv )
!
!*******************************************************************************
!
!! DATA_WRITE writes the internal graphics data to a file.
!
!
!  Modified:
!
!    14 October 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real COR3(3,COR3_MAX), the coordinates of points.
!
!    Input, integer COR3_MATERIAL(COR3_MAX), the material of each node.
!
!    Input, real COR3_TEX_UV(2,COR3_MAX), UV texture coordinates for nodes.
!
!    Input, integer FACE(ORDER_MAX,FACE_MAX), the nodes making faces.
!
!    Input, integer FACE_MATERIAL(FACE_MAX), the material of each face.
!
!    Input, real FACE_NORMAL(3,FACE_MAX), the normal vector at each face.
!
!    Input, integer FACE_ORDER(FACE_MAX), the number of vertices per face.
!
!    Input, character ( len = 100 ) FILEIN_NAME, the name of the input file.
!
!    Input/output, character ( len = 100 ) FILEOUT_NAME, the name of the output file.
!
!    Output, integer IERROR, an error flag.
!
!    Input, integer LINE_DEX(LINE_MAX), nodes forming a line, terminated by -1.
!
!    Input, integer LINE_MAT(LINE_MAX), material index for line.
!
!    Input, character ( len = 100 ) MATERIAL_NAME(MATERIAL_MAX), material names.
!
!    Input, integer COR3_MAX, the maximum number of 3D points.
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
!    Input, integer COR3_NUM, the number of points.
!
!    Input, integer FACE_NUM, the number of faces.
!
!    Input, integer LINE_NUM, the number of line definition items.
!
!    Input, integer MATERIAL_NUM, the number of materials.
!
!    Input, character ( len = 100 ) TEXTURE_NAME(TEXTURE_MAX), texture names.
!
!    Input, real TRANSFORM_MATRIX(4,4), the transformation matrix.
!
!    Input, integer VERTEX_MAT(ORDER_MAX,FACE_MAX), vertex materials.
!
!    Input, real VERTEX_NORMAL(3,ORDER_MAX,FACE_MAX), normals at vertices.  
!
!    Input, real VERTEX_TEX_UV(2,ORDER_MAX,FACE_MAX), vertex texture coordinates.
!
  integer, parameter :: iunit = 1
!
  integer cor3_max
  integer face_max
  integer line_max
  integer material_max
  integer order_max
  integer texture_max
!
  real cor2(2,cor3_max)
  integer cor2_max
  integer cor2_num
  real cor3(3,cor3_max)
  integer cor3_material(cor3_max)
  real cor3_normal(3,cor3_max)
  integer cor3_num
  real cor3_tex_uv(2,cor3_max)
  logical debug
  integer face(order_max,face_max)
  integer face_material(face_max)
  real face_normal(3,face_max)
  integer face_num
  integer face_order(face_max)
  real face_tex_uv(2,face_max)
  character ( len = 100 ) filein_name
  character ( len = 100 ) fileout_name
  character ( len = 10 ) fileout_type
  integer ierror
  integer ios
  integer line_dex(line_max)
  integer line_material(line_max)
  integer line_num
  integer line_num_save
  integer line_prune
  character ( len = 100 ) material_name(material_max)
  integer material_num
  real material_rgba(4,material_max)
  character ( len = 100 ) object_name
  logical s_eqi
  character ( len = 100 ) texture_name(texture_max)
  integer texture_num
  integer vertex_material(order_max,face_max)
  real vertex_normal(3,order_max,face_max)
  real vertex_tex_uv(2,order_max,face_max)
!
  ierror = 0
!
!  Open the file.
!
  open ( unit = iunit, file = fileout_name, status = 'replace', iostat = ios )
  
  if ( ios /= 0 ) then
    ierror = ios
    write ( *, * ) ' '
    write ( *, * ) 'DATA_WRITE - Warning!'
    write ( *, * ) '  Could not open the output file.'
    return
  end if
!
!  Write an Autodesk file...
!
  if ( s_eqi ( fileout_type, 'ASE' ) ) then
 
    call ase_write ( cor3, face, face_normal, face_order, &
      filein_name, fileout_name, iunit, cor3_max, face_max, &
      order_max, cor3_num, face_num, vertex_normal )
!
!  ...or a BYU file...
!
  else if ( s_eqi ( fileout_type, 'BYU' ) ) then

    call byu_write ( cor3, face, face_order, iunit, cor3_max, &
      face_max, order_max, cor3_num, face_num )
!
!  ...or a DXF file...
!
  else if ( s_eqi ( fileout_type, 'DXF' ) ) then
 
    call dxf_write ( cor3, face, face_order, filein_name, &
      fileout_name, iunit, line_dex, cor3_max, face_max, line_max, &
      order_max, face_num, line_num )
!
!  ...or an HRC SoftImage file...
!
  else if ( s_eqi ( fileout_type, 'HRC' ) ) then
 
    call hrc_write ( cor3, face, face_material, face_order, &
      fileout_name, iunit, line_dex, material_name, material_rgba, &
      cor3_max, face_max, line_max, material_max, order_max, texture_max, &
      cor3_num, face_num, line_num, material_num, texture_num, texture_name, &
      vertex_normal, vertex_tex_uv )
!
!  ...or an IV Inventor file...
!
  else if ( s_eqi ( fileout_type, 'IV' ) ) then
 
    call iv_write ( cor3, cor3_normal, face, face_order, &
      filein_name, fileout_name, iunit, line_dex, line_material, &
      material_rgba, cor3_max, face_max, line_max, material_max, &
      order_max, texture_max, cor3_num, face_num, line_num, material_num, &
      texture_num, texture_name, vertex_material, vertex_tex_uv )
!
!  ...or a WaveFront OBJ file...
!
  else if ( s_eqi ( fileout_type, 'OBJ' ) ) then
 
    call obj_write ( cor3, face, face_order, filein_name, &
      fileout_name, iunit, line_dex, cor3_max, face_max, line_max, &
      order_max, cor3_num, face_num, line_num, vertex_normal )
!
!  ...or a POV file...
!
  else if ( s_eqi ( fileout_type, 'POV' ) ) then
 
    call pov_write ( cor3, face, face_material, face_order, filein_name, &
      fileout_name, iunit, material_rgba, cor3_max, face_max, &
      material_max, order_max, face_num, material_num, vertex_normal )
!
!  ...or a PS file...
!
  else if ( s_eqi ( fileout_type, 'PS' ) ) then
 
    write ( *, * ) ' '
    write ( *, * ) 'WATCH OUT!'
    write ( *, * ) 'PS_WRITE not ready yet!'

    cor2_num = cor3_num
    cor2_max = cor3_max

    call project_2d ( cor2, cor3, ierror, cor2_max, cor3_max, cor2_num, &
      cor3_num )

    if ( ierror == 0 ) then

      call ps_write ( cor2, face, face_material, face_order, &
        fileout_name, iunit, line_dex, line_material, material_rgba, &
        cor2_max, face_max, line_max, material_max, order_max, cor2_num, &
        face_num, line_num )

    else

      write ( *, * ) ' '
      write ( *, * ) 'DATA_WRITE - Error!'
      write ( *, * ) '  2D projection canceled.'
      return

    end if
!
!  ...or an SMF file...
!
  else if ( s_eqi ( fileout_type, 'SMF' ) ) then
 
    call smf_write ( cor3, cor3_material, cor3_normal, cor3_tex_uv, &
      face, face_order, filein_name, fileout_name, iunit, &
      material_rgba, cor3_max, face_max, material_max, order_max, &
      texture_max, cor3_num, face_num, texture_num, texture_name )
!
!  ...or an ASCII STL file...
!
  else if ( s_eqi ( fileout_type, 'STL' ) .or. &
            s_eqi ( fileout_type, 'STLA' ) ) then
 
    call stla_write ( cor3, face, face_normal, face_order, &
      filein_name, iunit, cor3_max, face_max, order_max, face_num )
!
!  ...or a TEC file...
!
  else if ( s_eqi ( fileout_type, 'TEC' ) ) then
 
    call tec_write ( cor3, cor3_material, face, face_order, &
      fileout_name, iunit, material_rgba, cor3_max, face_max, &
      material_max, order_max, cor3_num, face_num )
!
!  ...or a TRI/TRIA file...
!
  else if ( s_eqi ( fileout_type, 'TRI' ) .or. &
            s_eqi ( fileout_type, 'TRIA' ) ) then

    call tria_write ( cor3, cor3_normal, face, face_order, &
      iunit, cor3_max, face_max, order_max, face_num )
!
!  ...or a TXT text file...
!
  else if ( s_eqi ( fileout_type, 'TXT' ) ) then

    call txt_write ( cor3, cor3_material, cor3_normal, cor3_tex_uv, &
      face, face_material, face_normal, face_order, face_tex_uv, &
      filein_name, fileout_name, iunit, line_dex, line_material, &
      material_name, material_rgba, cor3_max, face_max, line_max, &
      material_max, order_max, texture_max, cor3_num, face_num, &
      line_num, material_num, texture_num, object_name, texture_name, &
      vertex_material, vertex_normal, vertex_tex_uv )
!
!  ...or a UCD file...
!
  else if ( s_eqi ( fileout_type, 'UCD' ) ) then
 
    call ucd_write ( cor3, cor3_material, face, face_material, face_order, &
      fileout_name, iunit, material_rgba, cor3_max, &
      face_max, material_max, order_max, cor3_num, face_num, material_num )
!
!  ...or a VLA file...
!
  else if ( s_eqi ( fileout_type, 'VLA' ) ) then
 
    line_num_save = line_num

    if ( face_num > 0 ) then

      write ( *, * ) ' '
      write ( *, * ) 'DATA_WRITE - Note:'
      write ( *, * ) '  Face information will be temporarily'
      write ( *, * ) '  converted to line information for '
      write ( *, * ) '  the VLA output.'

      call face_to_line ( debug, face, face_order, line_dex, &
        line_material, line_prune, face_max, line_max, order_max, &
        face_num, line_num, vertex_material )

      if ( line_num > line_max ) then

        write ( *, * ) ' '
        write ( *, * ) 'DATA_WRITE - Note:'
        write ( *, * ) '  Some face information was lost.'
        write ( *, * ) '  The maximum number of lines is ', line_max
        write ( *, * ) '  but we would need at least ', line_num
        line_num = line_max

      end if

    end if

    call vla_write ( cor3, filein_name, fileout_name, iunit, line_dex, &
      cor3_max, line_max, line_num )

    line_num = line_num_save
!
!  ...or a VRML file...
!
  else if ( s_eqi ( fileout_type, 'WRL' ) ) then
 
    call vrml_write ( cor3, face, face_order, filein_name, fileout_name, &
      iunit, line_dex, line_material, material_rgba, cor3_max, face_max, &
      line_max, material_max, order_max, cor3_num, face_num, line_num, material_num, &
      vertex_material )
!
!  ...or an XGL file...
!
  else if ( s_eqi ( fileout_type, 'XGL' ) ) then
 
    call xgl_write ( cor3, cor3_max, cor3_normal, cor3_num, face, &
      face_material, face_max, face_num, face_order, iunit, material_max, &
      material_num, material_rgba, order_max )
!
!  ...or what the hell happened?
!
  else

    write ( *, * ) ' '
    write ( *, * ) 'DATA_WRITE - Warning!'
    write ( *, * ) '  Unrecognized output file type.'
    ierror = 1
 
  end if
!
!  Close the file.
!
  close ( unit = iunit )
 
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

function dot0_3d ( x0, y0, z0, x1, y1, z1, x2, y2, z2 )
!
!*******************************************************************************
!
!! DOT0_3D computes the dot product of (P1-P0) and (P2-P0) in 3D.
!
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
!    Input, real X0, Y0, Z0, the coordinates of the point P0.
!
!    Input, real X1, Y1, Z1, the coordinates of the point P1.
!
!    Input, real X2, Y2, Z2, the coordinates of the point P2.
!
!    Output, real DOT0_3D, the dot product of (P1-P0) and (P2-P0).
!
  real dot0_3d
  real x0
  real x1
  real x2
  real y0
  real y1
  real y2
  real z0
  real z1
  real z2
!
  dot0_3d = ( x1 - x0 ) * ( x2 - x0 ) + ( y1 - y0 ) * ( y2 - y0 ) + &
    ( z1 - z0 ) * ( z2 - z0 )
 
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
  integer i
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
subroutine dxf_write ( cor3, face, face_order, filein_name, fileout_name, &
  iunit, line_dex, cor3_max, face_max, line_max, order_max, face_num, &
  line_num )
!
!*******************************************************************************
!
!! DXF_WRITE writes graphics data to an AutoCAD DXF file.
!
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
!
!  Modified:
!
!    23 May 1999
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
!    Input, character ( len = 100 ) FILEIN_NAME, the name of the input file.
!
!    Input, character ( len = 100 ) FILEOUT_NAME, the name of the output file.
!
!    Input, integer IUNIT, the FORTRAN unit to which output is written.
!
!    Input, integer LINE_DEX(LINE_MAX), nodes forming a line, terminated by -1.
!
!    Input, integer COR3_MAX, the maximum number of points.
!
!    Input, integer FACE_MAX, the maximum number of faces.
!
!    Input, integer LINE_MAX, the maximum number of line definition items.
!
!    Input, integer ORDER_MAX, the maximum number of vertices per face.
!
!    Input, integer FACE_NUM, the number of faces.
!
!    Input, integer LINE_NUM, the number of line definition items.
!
  integer, parameter :: OFFSET = 1
!
  integer cor3_max
  integer face_max
  integer line_max
  integer order_max
!
  real cor3(3,cor3_max)
  integer face(order_max,face_max)
  integer face_num
  integer face_order(face_max)
  character ( len = 100 ) filein_name
  character ( len = 100 ) fileout_name
  integer i
  integer icor3
  integer iface
  integer iunit
  integer ivert
  integer jcor3
  integer line_dex(line_max)
  integer line_num
  logical newline
  integer text_num
!
  text_num = 0

  write ( iunit, '(i3)' ) 0
  write ( iunit, '(a)' ) 'SECTION'
  write ( iunit, '(i3)' ) 2
  write ( iunit, '(a)' ) 'HEADER'
  write ( iunit, '(i3)' ) 999
  write ( iunit, '(a,a)' ) trim ( fileout_name ), ' created by IVREAD.'
  write ( iunit, '(i3)' ) 999
  write ( iunit, '(a,a)' ) 'Original data in ', trim ( filein_name ) // '.'
  write ( iunit, '(i3)' ) 0
  write ( iunit, '(a)' ) 'ENDSEC'
  text_num = text_num + 10

  write ( iunit, '(i3)' ) 0
  write ( iunit, '(a)' ) 'SECTION'
  write ( iunit, '(i3)' ) 2
  write ( iunit, '(a)' ) 'TABLES'
  write ( iunit, '(i3)' ) 0
  write ( iunit, '(a)' ) 'ENDSEC'
  text_num = text_num + 6

  write ( iunit, '(i3)' ) 0
  write ( iunit, '(a)' ) 'SECTION'
  write ( iunit, '(i3)' ) 2
  write ( iunit, '(a)' ) 'BLOCKS'
  write ( iunit, '(i3)' ) 0
  write ( iunit, '(a)' ) 'ENDSEC'
  text_num = text_num + 6

  write ( iunit, '(i3)' ) 0
  write ( iunit, '(a)' ) 'SECTION'
  write ( iunit, '(i3)' ) 2
  write ( iunit, '(a)' ) 'ENTITIES'
  text_num = text_num + 4

  jcor3 = 0
  newline = .true.

  do i = 1, line_num

    icor3 = line_dex(i)

    if ( line_dex(i) - OFFSET == -1 ) then

      newline = .true.
!
!  LINE_DEX(I) is the index of a new point that begins or continues a line.
!
    else
!
!  LINE_DEX(I) is the index of a new point that continues a line. 
!    Output the pair of points that define this segment of the line.
!
      if ( .not. newline ) then

        write ( iunit, '(i3)' ) 0
        write ( iunit, '(a)' ) 'LINE'
        write ( iunit, '(i3)' ) 8
        write ( iunit, '(i3)' ) 0

        write ( iunit, '(i3)' ) 10
        write ( iunit, '(g12.4)' ) cor3(1,jcor3)
        write ( iunit, '(i3)' ) 20
        write ( iunit, '(g12.4)' ) cor3(2,jcor3)
        write ( iunit, '(i3)' ) 30
        write ( iunit, '(g12.4)' ) cor3(3,jcor3)

        write ( iunit, '(i3)' ) 11
        write ( iunit, '(g12.4)' ) cor3(1,icor3)
        write ( iunit, '(i3)' ) 21
        write ( iunit, '(g12.4)' ) cor3(2,icor3)
        write ( iunit, '(i3)' ) 31
        write ( iunit, '(g12.4)' ) cor3(3,icor3)

        text_num = text_num + 16
 
      end if
!
!  Save the index of this new point, and note that a line is in progress.
!
      jcor3 = icor3
      newline = .false.

    end if
 
  end do
!
!  Handle faces.
!  This is going to fail bigtime if FACE_ORDER is larger than 9
!
  do iface = 1, face_num
    
    write ( iunit, '(a)' ) '  0'
    write ( iunit, '(a)' ) '3DFACE'
    write ( iunit, '(a)' ) '  8'
    write ( iunit, '(a)' ) '  Cube'
    text_num = text_num + 4
    
    do ivert = 1, face_order(iface)

      icor3 = face(ivert,iface)
  
      write ( iunit, '(i3)' ) 10 + ivert - 1
      write ( iunit, '(g12.4)' ) cor3(1,icor3)
      write ( iunit, '(i3)' ) 20 + ivert - 1
      write ( iunit, '(g12.4)' ) cor3(2,icor3)
      write ( iunit, '(i3)' ) 30 + ivert - 1
      write ( iunit, '(g12.4)' ) cor3(3,icor3)

      text_num = text_num + 6

    end do
  end do

  write ( iunit, '(i3)' ) 0
  write ( iunit, '(a)' ) 'ENDSEC'
  write ( iunit, '(i3)' ) 0
  write ( iunit, '(a)' ) 'EOF'
 
  text_num = text_num + 4
!
!  Report.
!
  write ( *, * ) ' '
  write ( *, * ) 'DXF_WRITE - Wrote ', text_num, ' text lines.'
 
  return
end
subroutine edge_add_nodes ( edge, edge_max, edge_num, iface, n1, n2, ierror )
!
!*******************************************************************************
!
!! EDGE_ADD_NODES adds the edge defined by two nodes to the edge list.
!
!
!  Modified:
!
!    12 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer EDGE(4,EDGE_MAX), contains edge information.
!    EDGE(1,I) is the starting node of edge I;
!    EDGE(2,I) is the ending node of edge I;
!    EDGE(3,I) is the positive face;
!    EDGE(4,I) is the negative face, if any.
!
!    Input, integer EDGE_MAX, the maximum number of edges.
!
!    Input/output, integer EDGE_NUM, the number of edges.
!
!    Input, integer IFACE, the face to which the nodes belong.
!
!    Input, integer N1, N2, two nodes which form an edge.
!
!    Output, integer IERROR, error flag, 0 = no error, nonzero = error.
!
  integer edge_max
!
  integer edge(4,edge_max)
  integer edge_num
  integer ierror
  integer iface
  integer n1
  integer n2
!
  if ( edge_num < edge_max ) then
    edge_num = edge_num + 1
    edge(1,edge_num) = n1
    edge(2,edge_num) = n2
    edge(3,edge_num) = iface
    edge(4,edge_num) = 0
    ierror = 0
  else
    write ( *, * ) ' '
    write ( *, * ) 'EDGE_ADD_NODES - Fatal error!'
    write ( *, * ) '  Exceeding EDGE_MAX = ', edge_max
    ierror = 1
  end if

  return
end
subroutine edge_bound ( edge, edge_max, edge_num )
!
!*******************************************************************************
!
!! EDGE_BOUND reports the edges which are part of the boundary.
!
!
!  Modified:
!
!    12 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer EDGE(4,EDGE_MAX), contains edge information.
!    EDGE(1,I) is the starting node of edge I;
!    EDGE(2,I) is the ending node of edge I;
!    EDGE(3,I) is the positive face;
!    EDGE(4,I) is the negative face, if any.
!
!    Input, integer EDGE_MAX, the maximum number of edges.
!
!    Input, integer EDGE_NUM, the number of edges.
!
  logical, parameter :: DEBUG = .false.
!
  integer edge_max
!
  integer bound_num
  integer edge(4,edge_max)
  integer edge_num
  integer iedge
!
  if ( DEBUG ) then
    write ( *, * ) ' '
    write ( *, * ) 'Boundary edges:'
    write ( *, * ) ' '
  end if

  bound_num = 0

  do iedge = 1, edge_num
    if ( edge(4,iedge) == 0 ) then
      bound_num = bound_num + 1
      if ( DEBUG ) then
        write ( *, * ) edge(2,iedge), edge(1,iedge)
      end if
    end if
  end do

  write ( *, * ) 'EDGE_BOUND found ', bound_num, ' boundary edges.'

  return
end
subroutine edge_match_face ( edge, edge_max, edge_num, facelist, n, index )
!
!*******************************************************************************
!
!! EDGE_MATCH_FACE seeks an edge common to a face and the edge list.
!
!
!  Note:
!
!    If a common edge is found, then the information in the face node
!    list is adjusted so that the first two entries correspond to the
!    matching edge in EDGE, but in reverse order.
!
!  Modified:
!
!    12 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer EDGE(4,EDGE_MAX), contains edge information.
!    EDGE(1,I) is the starting node of edge I;
!    EDGE(2,I) is the ending node of edge I;
!    EDGE(3,I) is the positive face;
!    EDGE(4,I) is the negative face, if any.
!
!    Input, integer EDGE_MAX, the maximum number of edges.
!
!    Input, integer EDGE_NUM, the number of edges.
!
!    Input/output, integer FACELIST(N), the list of nodes making a face.
!
!    Input, integer N, the number of nodes in the face.
!
!    Output, integer INDEX, the results of the search.
!    0, there is no edge common to the face and the EDGE array.
!    nonzero, edge INDEX is common to the face and the EDGE array.
!
  integer n
  integer edge_max
!
  integer edge(4,edge_max)
  integer edge_num
  integer facelist(n)
  integer iedge
  integer index
  integer j
  integer jp1
  integer n1
  integer n2
!
  index = 0

  if ( n <= 0 ) then
    return
  end if

  if ( edge_num <= 0 ) then
    return
  end if

  do j = 1, n

    if ( j == n ) then
      jp1 = 1
    else
      jp1 = j + 1
    end if

    n1 = facelist(j)
    n2 = facelist(jp1)

    do iedge = 1, edge_num

      if ( edge(1,iedge) == n2 .and. edge(2,iedge) == n1 ) then

        call ivec_rotate ( n, 1 - j, facelist )

        index = iedge
        return

      else if ( edge(1,iedge) == n1 .and. edge(2,iedge) == n2 ) then

        call ivec_rotate ( n, n - jp1, facelist )

        call ivec_reverse ( n, facelist )

        index = iedge
        return

      end if

    end do
   
  end do

  return
end
subroutine edge_match_nodes ( edge, edge_max, edge_num, n1, n2, iedge )
!
!*******************************************************************************
!
!! EDGE_MATCH_NODES seeks an edge of the form (N1,N2) or (N2,N1) in EDGE.
!
!
!  Modified:
!
!    12 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer EDGE(4,EDGE_MAX), contains edge information.
!    EDGE(1,I) is the starting node of edge I;
!    EDGE(2,I) is the ending node of edge I;
!    EDGE(3,I) is the positive face;
!    EDGE(4,I) is the negative face, if any.
!
!    Input, integer EDGE_MAX, the maximum number of edges.
!
!    Input, integer EDGE_NUM, the number of edges.
!
!    Input, integer N1, N2, two nodes that form an edge.
!
!    Output, integer IEDGE, the results of the search.
!    0, no matching edge was found.
!    nonzero, edge IEDGE of the EDGE array matches (N1,N2) or (N2,N1).
!
  integer edge_max
!
  integer edge(4,edge_max)
  integer edge_num
  integer i
  integer iedge
  integer n1
  integer n2
!
  iedge = 0
  do i = 1, edge_num

    if ( ( n1 == edge(1,i) .and. n2 == edge(2,i) ) .or. &
         ( n2 == edge(1,i) .and. n1 == edge(2,i) ) ) then
      iedge = i
      return
    end if

  end do

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
  integer j
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
function enorm_nd ( n, x )
!
!*******************************************************************************
!
!! ENORM_ND computes the Euclidean norm of a vector in ND.
!
!
!  Modified:
!
!    27 July 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the dimension of the space.
!
!    Input, real X(N), the coordinates of the vector.
!
!    Output, real ENORM_ND, the Euclidean norm of the vector.
!
  integer n
!
  real enorm_nd
  integer i
  real x(n)
!
  enorm_nd = 0.0

  do i = 1, n
    enorm_nd = enorm_nd + x(i) * x(i)
  end do

  enorm_nd = sqrt ( enorm_nd )
 
  return
end
function enorm0_3d ( x0, y0, z0, x1, y1, z1 )
!
!*******************************************************************************
!
!! ENORM0_3D computes the Euclidean norm of (P1-P0) in 3D.
!
!
!  Modified:
!
!    16 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X0, Y0, Z0, X1, Y1, Z1, the coordinates of the points 
!    P0 and P1.
!
!    Output, real ENORM0_3D, the Euclidean norm of (P1-P0).
!
  real enorm0_3d
  real x0
  real x1
  real y0
  real y1
  real z0
  real z1
!
  enorm0_3d = sqrt ( ( x1 - x0 )**2 + ( y1 - y0 )**2 + ( z1 - z0 )**2 )
 
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
subroutine face_check ( edge, face, face_material, face_normal, face_object, &
  face_order, face_rank, face_tier, edge_max, face_max, order_max, edge_num, &
  face_num, object_num, vertex_material, vertex_normal )
!
!*******************************************************************************
!
!! FACE_CHECK checks and analyzes a set of faces.
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
!    Output, integer EDGE(4,EDGE_MAX), contains edge information.
!    EDGE(1,I) is the starting node of edge I;
!    EDGE(2,I) is the ending node of edge I;
!    EDGE(3,I) is the positive face;
!    EDGE(4,I) is the negative face, or 0 if the edge is used once.
!
!    Input, integer FACE(ORDER_MAX,FACE_NUM), the nodes making faces
!
!    Input, integer FACE_MATERIAL(FACE_MAX), the material of each face.
!
!    Output, integer FACE_OBJECT(FACE_NUM), describes the objects.
!    FACE_OBJECT(I) is the index of the edge-connected "object" to 
!    which face I belongs.
!
!    Input, integer FACE_ORDER(FACE_NUM), the number of nodes per face.
!
!    Output, integer FACE_RANK(FACE_NUM), is an ordered list of faces.
!    FACE_RANK(1) is the index of the face in the first tier of the 
!    first object, followed by second tier faces, and so on until
!    object one is complete.  Object two follows, and so on.
!
!    Output, integer FACE_TIER(FACE_NUM).  FACE_TIER(I) is the "tier"
!    of face I in its object.  The seed of the object has tier 1,
!    the neighbors of the seed have tier 2, and so on.
!
!    Input, integer EDGE_MAX, the maximum number of edges.
!
!    Input, integer ORDER_MAX, is the maximum number of nodes that can
!    make up a face, required to dimension FACE.
!
!    Output, integer EDGE_NUM, the number of edges.
!
!    Input, integer FACE_NUM, the number of faces.
!
!    Output, integer OBJECT_NUM, the number of objects.
!
!    Input, integer VERTEX_MAT(ORDER_MAX,FACE_MAX), vertex materials.
!
  integer edge_max
  integer order_max
  integer face_max
!
  integer edge(4,edge_max)
  integer edge_num
  integer face(order_max,face_max)
  integer face_material(face_max)
  real face_normal(3,face_max)
  integer face_num
  integer face_object(face_max)
  integer face_order(face_max)
  integer face_rank(face_max)
  integer face_tier(face_max)
  integer i
  integer ierror
  integer j
  integer nfix
  integer object_num
  integer vertex_material(order_max,face_max)
  real vertex_normal(3,order_max,face_max)
!
!  Organize the faces into layered objects.
!
  write ( *, * ) ' '
  write ( *, * ) 'Determine edge-connected objects.'

  call object_build ( face, face_object, face_order, face_rank, face_tier, &
    order_max, face_num, object_num )

  write ( *, * ) ' '
  write ( *, * ) 'Number of objects = ', object_num

  if ( face_num <= 20 ) then

    write ( *, * ) ' '
    write ( *, * ) 'Face, Object, Tier'
    write ( *, * ) ' '

    do i = 1, face_num
      write ( *, * ) i, face_object(i), face_tier(i)
    end do

  end if

  if ( face_num <= 20 ) then
    write ( *, * ) ' '
    write ( *, * ) 'Preferred order:'
    write ( *, * ) '  Order, Face'
    write ( *, * ) ' '
    do i = 1, face_num
      write ( *, * ) i, face_rank(i)
    end do
  end if
!
!  Reorder the faces by object and tier.
!
  write ( *, * ) ' '
  write ( *, * ) 'Reorder the faces.'

  call face_sort ( face, face_material, face_normal, face_object, &
    face_order, face_tier, face_max, order_max, face_num, &
    vertex_material, vertex_normal )

  if ( face_num <= 20 ) then
    write ( *, * ) ' '
    write ( *, * ) 'Face, Label, Object, Tier'
    write ( *, * ) ' '
    do i = 1, face_num
      write ( *, * ) i, face_rank(i), face_object(i), face_tier(i)
    end do
  end if
!
!  Construct the edge list.
!
  write ( *, * ) ' '
  write ( *, * ) 'Construct the edge list.'
  write ( *, * ) '(While doing so, check for edges used more'
  write ( *, * ) 'than twice.)'

  call face_to_edge ( edge, face, face_order, ierror, edge_max, &
    order_max, edge_num, face_num )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'FACE_CHECK - Fatal error!'
    write ( *, * ) '  FACE_TO_EDGE failed.'
    return
  end if

  if ( face_num <= 20 ) then
    write ( *, * ) ' '
    write ( *, * ) 'Edge, Node1, Node2, Face1, Face2, Tier, Object'
    write ( *, * ) ' '
    write ( *, * ) ' I, node1(i), node2(i), face1(i), face2(i)'
    write ( *, * ) ' '

    do i = 1, edge_num
      write ( *, '(10i3)' ) i, edge(1:4,i)
    end do

    write ( *, * ) ' '
    write ( *, * ) 'Face, Order, Nodes'
    write ( *, * ) ' '
    do i = 1, face_num
      write ( *, '(10i3)' ) i, face_order(i), ( face(j,i), j = 1, face_order(i) )
    end do
  end if
!
!  Now force faces to have a consistent orientation.
!
  write ( *, * ) ' '
  write ( *, * ) 'Force faces to consistent orientation.'
  
  call face_flip ( edge, face, face_order, edge_max, &
    order_max, nfix, edge_num, face_num )

  if ( face_num <= 20 ) then
    write ( *, * ) ' '
    write ( *, * ) 'Face, Order, Nodes'
    write ( *, * ) ' '
    do i = 1, face_num
      write ( *, '(10i3)' ) i, face_order(i), ( face(j,i), j = 1, face_order(i) )
    end do
  end if

  write ( *, * ) ' '
  write ( *, * ) 'List boundary edges.'

  call edge_bound ( edge, edge_max, edge_num )

  return
end
subroutine face_flip ( edge, face, face_order, edge_max, order_max, nfix, &
  edge_num, face_num )
!
!*******************************************************************************
!
!! FACE_FLIP flips faces to achieve a consistent orientation.
!
!
!  Modified:
!
!    12 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer EDGE(4,EDGE_MAX), contains edge information.
!    EDGE(1,I) is the starting node of edge I;
!    EDGE(2,I) is the ending node of edge I;
!    EDGE(3,I) is the positive face;
!    EDGE(4,I) is the negative face, if any.
!
!    Input, integer FACE(ORDER_MAX,FACE_NUM), the nodes making faces.
!
!    Input, integer FACE_ORDER(FACE_NUM), the number of nodes per face.
!
!    Input, integer EDGE_MAX, the maximum number of edges.
!
!    Input, integer ORDER_MAX, the maximum number of nodes that can
!    make up a face, required to dimension FACE.
!
!    Output, integer NFIX, the number of bad faces that were found.
!
!    Input, integer EDGE_NUM, the number of edges.
!
!    Input, integer FACE_NUM, the number of faces.
!
  integer edge_max
  integer order_max
  integer face_num
!
  integer edge(4,edge_max)
  integer edge_num
  integer f1
  integer f2
  integer face(order_max,face_num)
  integer face_order(face_num)
  integer iedge
  integer j
  integer jp1
  integer m1
  integer m2
  integer n1
  integer n2
  integer nfix
!
  nfix = 0

  do iedge = 1, edge_num

    n1 = edge(1,iedge)
    n2 = edge(2,iedge)
    f1 = edge(3,iedge)
    f2 = edge(4,iedge)
!
!  For now, just whine unless (N1,N2) is positive in F1 and negative in F2.
!
    if ( f1 /= 0 ) then

      do j = 1, face_order(f1)

        if ( j < face_order(f1) ) then
          jp1 = j + 1
        else
          jp1 = j
        end if

        m1 = face(j,f1)
        m2 = face(jp1,f1)

        if ( m1 == n1 .and. m2 == n2 ) then
          go to 10
        else if ( m1 == n2 .and. m2 == n1 ) then
          nfix = nfix + 1
          write ( *, * ) 'Bad orientation, face ', f1, ' side ', j
          go to 10
        end if

      end do

    end if

10      continue

    if ( f2 /= 0 ) then

    do j = 1, face_order(f2)

      if ( j < face_order(f2) ) then
        jp1 = j + 1
      else
        jp1 = j
      end if

      m1 = face(j,f2)
      m2 = face(jp1,f2)

      if ( m1 == n2 .and. m2 == n1 ) then
        go to 20
      else if ( m1 == n1 .and. m2 == n2 ) then
        nfix = nfix + 1
        write ( *, * ) 'Bad orientation, face ', f2, ' side ', j
        go to 20
      end if

    end do

    end if

20      continue

  end do

  write ( *, * ) 'FACE_FLIP found ', nfix, ' badly oriented faces.'

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
  integer ivert
  integer j
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
subroutine face_print ( cor3, face, face_material, face_normal, face_order, &
  iface, cor3_max, face_max, order_max, face_num, vertex_material, vertex_normal )
!
!*******************************************************************************
!
!! FACE_PRINT prints out information about a face.
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
!    Input, real COR3(3,COR3_MAX), the coordinates of points.
!
!    Input, integer FACE(ORDER_MAX,FACE_MAX), the nodes making faces.
!
!    Input, integer FACE_MATERIAL(FACE_MAX), the material of each face.
!
!    Input, real FACE_NORMAL(3,FACE_MAX), the normal vector at each face.
!
!    Input, integer FACE_ORDER(FACE_MAX), the number of vertices per face.
!
!    Input, integer IFACE, the face about which information is desired.
!
!    Input, integer COR3_MAX, the maximum number of points.
!
!    Input, integer FACE_MAX, the maximum number of faces.
!
!    Input, integer ORDER_MAX, the maximum number of vertices per face.
!
!    Input, integer FACE_NUM, the number of faces.
!
!    Input, integer VERTEX_MAT(ORDER_MAX,FACE_MAX), vertex materials.
!
!    Input, real VERTEX_NORMAL(3,ORDER_MAX,FACE_MAX), normals at vertices.  
!
  integer cor3_max
  integer face_max
  integer order_max
!
  real cor3(3,cor3_max)
  integer face(order_max,face_max)
  integer face_material(face_max)
  real face_normal(3,face_max)
  integer face_num
  integer face_order(order_max)
  integer i
  integer iface
  integer ivert
  integer j
  integer k
  integer vertex_material(order_max,face_max)
  real vertex_normal(3,order_max,face_max)
!
  if ( iface < 1 .or. iface > face_num ) then
    write ( *, * ) ' '
    write ( *, * ) 'FACE_PRINT - Error!'
    write ( *, * ) '  Face indices must be between 1 and ', face_num
    write ( *, * ) '  But your requested value was ', iface
    return
  end if

  write ( *, * ) ' '
  write ( *, * ) 'FACE_PRINT'
  write ( *, * ) '  Information about face ', iface
  write ( *, * ) ' '
  write ( *, * ) '  The number of vertices is ', face_order(iface)
  write ( *, * ) '  Face material is ', face_material(iface)
  write ( *, * ) ' '
  write ( *, * ) '  Vertex list:'
  write ( *, * ) '    Vertex #, Node #, Material #, X, Y, Z'
  write ( *, * ) ' '
  do ivert = 1, face_order(iface)
    j = face(ivert,iface)
    k = vertex_material(ivert,iface)
    write ( *, '(3i8,3f10.4)' ) ivert, j, k, cor3(1,j), cor3(2,j), cor3(3,j)
  end do

  write ( *, * ) ' '
  write ( *, * ) '  Face normal vector:'
  write ( *, * ) ' '
  write ( *, '(3f10.4)' ) face_normal(1:3,iface)

  write ( *, * ) ' '
  write ( *, * ) '  Vertex face normals:'
  write ( *, * ) ' '
  do ivert = 1, face_order(iface)
    write ( *, '(i8,3f10.4)' ) ivert, vertex_normal(1:3,ivert,iface)
  end do

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
  integer i
  integer iface
  integer itemp
  integer ivert
  integer j
  integer m
  real temp
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
subroutine face_sort ( face, face_material, face_normal, face_object, face_order, &
  face_tier, face_max, order_max, face_num, vertex_material, vertex_normal )
!
!*******************************************************************************
!
!! FACE_SORT renumbers the faces in order of object and tier.
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
!    Input/output, integer FACE(ORDER_MAX,FACE_NUM), the nodes making faces.
!
!    Input/output, integer FACE_MATERIAL(FACE_MAX), the material of each face.
!
!    Input/output, integer FACE_OBJECT(FACE_NUM), describes the objects.
!    FACE_OBJECT(I) is the index of the edge-connected "object" to 
!    which face I belongs.
!
!    Input/output, integer FACE_ORDER(FACE_NUM), the number of nodes per face.
!
!    Input/output, integer FACE_TIER(FACE_NUM).  FACE_TIER(I) is the "tier"
!    of face I in its object.  The seed of the object has tier 1,
!    the neighbors of the seed have tier 2, and so on.
!
!    Input, integer ORDER_MAX, is the maximum number of nodes that can
!    make up a face, required to dimension FACE.
!
!    Input, integer FACE_NUM, the number of faces.
!
!    Input/output, integer VERTEX_MAT(ORDER_MAX,FACE_MAX), vertex materials.
!
  integer order_max
  integer face_max
!
  integer face(order_max,face_max)
  integer face_material(face_max)
  real face_normal(3,face_max)
  integer face_num
  integer face_object(face_max)
  integer face_order(face_max)
  integer face_tier(face_max)
  integer i
  integer iface
  integer indx
  integer isgn
  integer itemp
  integer ivert
  integer jface
  real temp
  integer vertex_material(order_max,face_max)
  real vertex_normal(3,order_max,face_max)
!
  iface = 0
  jface = 0
  indx = 0
  isgn = 0

  do

    call sort_heap_external ( face_num, indx, iface, jface, isgn )
!
!  Interchange faces IFACE and JFACE.
!
    if ( indx > 0 ) then

      do i = 1, order_max
        call i_swap ( face(i,iface), face(i,jface) )
      end do

      call i_swap ( face_material(iface),    face_material(jface) )
      call i_swap ( face_object(iface), face_object(jface) )
      call i_swap ( face_order(iface),  face_order(jface) )
      call i_swap ( face_tier(iface),   face_tier(jface) )

      do i = 1, 3
        call r_swap ( face_normal(i,iface), face_normal(i,jface) )
      end do

      do ivert = 1, order_max
        call i_swap ( vertex_material(ivert,iface), vertex_material(ivert,jface) )
      end do

      do i = 1, 3
        do ivert = 1, order_max
          call r_swap ( vertex_normal(i,ivert,iface), vertex_normal(i,ivert,jface) )
        end do
      end do
!
!  Compare faces IFACE and JFACE.
!
    else if ( indx < 0 ) then

      if ( ( face_object(iface) < face_object(jface) ) .or. &
           ( face_object(iface) == face_object(jface) .and. &
             face_tier(iface) < face_tier(jface) ) ) then
        isgn = -1
      else
        isgn = +1
      end if

    else if ( indx == 0 ) then

      exit

    end if

  end do

  return
end
subroutine face_subset ( cor3, face, face_material, face_normal, &
  face_order, ierror, list, cor3_max, face_max, order_max, &
  cor3_num, face_num, line_num, vertex_material, vertex_normal )
!
!*******************************************************************************
!
!! FACE_SUBSET selects a subset of the current faces as the new object.
!
!
!  Warning:
!
!    The original graphic object is overwritten by the new one.
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
!    Input, real COR3(3,COR3_MAX), the coordinates of points.
!
!    Input/output, integer FACE(ORDER_MAX,FACE_MAX), the nodes making faces.
!
!    Input/output, integer FACE_MATERIAL(FACE_MAX), the material of each face.
!
!    Input/output, real FACE_NORMAL(3,FACE_MAX), the normal vector at each face.
!
!    Input, integer FACE_ORDER(FACE_MAX), the number of vertices per face.
!
!    Workspace, integer LIST(COR3_MAX), contains the indices of the points
!    to be copied from the old graphics object to the new one.
!
!    Input, integer COR3_MAX, the maximum number of points.
!
!    Input, integer FACE_MAX, the maximum number of faces.
!
!    Input, integer ORDER_MAX, the maximum number of vertices per face.
!
!    Input/output, integer LINE_NUM, the number of lines.  
!    This routine resets LINE_NUM to zero, since we will be dropping
!    as many points as possible.
!
!    Input/output, integer VERTEX_MAT(ORDER_MAX,FACE_MAX), vertex materials.
!
!    Input/output, real VERTEX_NORMAL(3,ORDER_MAX,FACE_MAX), normals at vertices.  
!
  integer cor3_max
  integer face_max
  integer order_max
!
  real cor3(3,cor3_max)
  integer cor3_num
  integer cor3_num2
  integer face(order_max,face_max)
  integer face_material(face_max)
  real face_normal(3,face_max)
  integer face_num
  integer face_order(face_max)
  integer i
  integer ierror
  integer iface
  integer iface1
  integer iface2
  integer inc
  integer ivert
  integer j
  integer k
  integer line_num
  integer list(cor3_max)
  integer vertex_material(order_max,face_max)
  real vertex_normal(3,order_max,face_max)
!
  ierror = 0

  line_num = 0
!
!  Get the first and last faces to save, IFACE1 and IFACE2.
!
  write ( *, * ) ' '
  write ( *, * ) 'Enter lowest face number to save,'
  write ( *, * ) 'between 1 and ', face_num
  read ( *, * ) iface1
  if ( iface1 < 1 .or. iface1 > face_num ) then
    write ( *, * ) 'Illegal choice!'
    ierror = 1
    return
  end if

  write ( *, * ) ' '
  write ( *, * ) 'Enter highest face number to save'
  write ( *, * ) 'between ', iface1, ' and ', face_num
  read ( *, * ) iface2
  if ( iface2 < iface1 .or. iface2 > face_num ) then
    write ( *, * ) 'Illegal choice!'
    ierror = 1
    return
  end if

  inc = iface1 - 1
!
!  "Slide" the data for the saved faces down the face arrays.
!
  do iface = 1, iface2 + 1 - iface1

    face_material(iface) = face_material(iface+inc)
    face_order(iface) = face_order(iface+inc)

    do ivert = 1, order_max
      face(ivert,iface) = face(ivert,iface+inc)
      vertex_material(ivert,iface) = vertex_material(ivert,iface+inc)
      vertex_normal(1:3,ivert,iface) = vertex_normal(1:3,ivert,iface+inc)
    end do

    face_normal(1:3,iface) = face_normal(1:3,iface+inc)

  end do
!
!  Now reset the number of faces.
!
  face_num = iface2 + 1 - iface1
!
!  Now, for each point I, set LIST(I) = J if point I is the J-th
!  point we are going to save, and 0 otherwise.  Then J will be
!  the new label of point I.
!
  list(1:cor3_num) = 0
  cor3_num2 = 0

  do iface = 1, face_num

    do ivert = 1, face_order(iface)

      j = face(ivert,iface)

      if ( list(j) == 0 ) then
        cor3_num2 = cor3_num2 + 1
        list(j) = cor3_num2
      end if

    end do

  end do
!
!  Now make the nonzero list entries rise in order, so that
!  we can compress the COR3 data in a minute.
!
  cor3_num2 = 0
  do i = 1, cor3_num
    if ( list(i) /= 0 ) then
      cor3_num2 = cor3_num2 + 1
      list(i) = cor3_num2
    end if
  end do
!
!  Relabel the FACE array with the new node indices.
!
  do iface = 1, face_num
    do ivert = 1, face_order(iface)
      j = face(ivert,iface)
      face(ivert,iface) = list(j)
    end do
  end do
!
!  Rebuild the COR3 array by sliding data down.
!
  do i = 1, cor3_num
    k = list(i)
    if ( k /= 0 ) then
      cor3(1:3,k) = cor3(1:3,i)
    end if
  end do

  cor3_num = cor3_num2

  return
end
subroutine face_to_edge ( edge, face, face_order, ierror, edge_max, order_max, &
  edge_num, face_num )
!
!*******************************************************************************
!
!! FACE_TO_EDGE converts face data to edge data.
!
!
!  Discussion:
!
!    The computation will fail if:
!
!      More than two faces claim to share an edge (Node1,Node2).
!      Not enough storage is set aside by EDGE_MAX.
!
!    If is expected that the edge (Node1,Node2) in Face1 is traversed in
!    the opposite sense, as (Node2,Node1), in Face2.  If this is not the
!    case, then some faces may need to be reoriented, but that will not
!    affect the computation.
!    
!  Modified:
!
!    12 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer EDGE(4,EDGE_MAX), contains edge information.
!    EDGE(1,I) is the starting node of edge I;
!    EDGE(2,I) is the ending node of edge I;
!    EDGE(3,I) is the positive face;
!    EDGE(4,I) is the negative face, or 0 if the edge is used once.
!
!    Input, integer FACE(ORDER_MAX,FACE_NUM), the nodes making faces.
!
!    Input, integer FACE_ORDER(FACE_NUM), the number of nodes per face.
!
!    Output, integer IERROR, error flag: 0 = no error, nonzero = error.
!
!    Input, integer EDGE_MAX, the maximum number of edges.
!
!    Input, integer ORDER_MAX, the maximum number of nodes that can
!    make up a face, required to dimension FACE.
!
!    Output, integer EDGE_NUM, the number of edges.
!
!    Input, integer FACE_NUM, the number of faces.
!
  integer edge_max
  integer order_max
  integer face_num
!
  integer edge(4,edge_max)
  integer edge_num
  integer face(order_max,face_num)
  integer face_order(face_num)
  integer i
  integer iedge
  integer ierror
  integer iface
  integer index
  integer j
  integer jp1
  integer n1
  integer n2
!
!  Initialize.
!
  ierror = 0

  edge(1:4,1:edge_max) = 0

  edge_num = 0
!
!  Consider face #I.
!
  do iface = 1, face_num
!
!  Seek an edge of face IFACE that already occurs in the edge list.
!  If there is one, then slide and reverse the entries in FACE(*,IFACE)
!  so that that edge occurs first, and in the opposite sense to its
!  occurrence in the edge list.
!
    call edge_match_face ( edge, edge_max, edge_num, face(1,iface), &
      face_order(iface), index )
!
!  Now, in any case, we know that the first two nodes in FACE(*,IFACE)
!  are the negative of an existing edge, or no nodes in FACE(*,IFACE)
!  occur in any existing edge.
!
    do j = 1, face_order(iface)

      n1 = face(j,iface)
      if ( j == face_order(iface) ) then
        jp1 = 1
      else
        jp1 = j + 1
      end if

      n2 = face(jp1,iface)

      call edge_match_nodes ( edge, edge_max, edge_num, n1, n2, iedge )

      if ( iedge == 0 ) then

        call edge_add_nodes ( edge, edge_max, edge_num, iface, n1, n2, ierror )

        if ( ierror /= 0 ) then
          write ( *, * ) ' '
          write ( *, * ) 'FACE_TO_EDGE - Fatal error!'
          write ( *, * ) '  EDGE_ADD_NODES failed.'
          ierror = 1
          return
        end if

      else if ( edge(4,iedge) == 0 ) then

        edge(4,iedge) = iface

      else 

        write ( *, * ) ' '
        write ( *, * ) 'FACE_TO_EDGE - Fatal error!'
        write ( *, * ) '  Edge between nodes ', edge(1,iedge), edge(2,iedge)
        write ( *, * ) '  is used at least 3 times, by faces:'
        write ( *, * ) edge(3,iedge), edge(4,iedge), iface
        ierror = 1
        return

      end if

    end do
  end do

  return
end
subroutine face_to_line ( debug, face, face_order, line_dex, line_material, &
  line_prune, face_max, line_max, order_max, face_num, line_num, vertex_material )
!
!*******************************************************************************
!
!! FACE_TO_LINE converts face information to line information.
!
!
!  Discussion:
!
!    In some cases, the graphic information represented by polygonal faces
!    must be converted to a representation based solely on line segments.
!    This is particularly true if a VLA file is being written.
!
!  Modified:
!
!    29 September 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer FACE(ORDER_MAX,FACE_MAX), the nodes making faces.
!
!    Input, integer FACE_ORDER(FACE_MAX), the number of vertices per face.
!
!    Input/output, integer LINE_DEX(LINE_MAX), nodes forming a line, 
!    terminated by -1.
!
!    Output, integer LINE_MAT(LINE_MAX), material index for line.
!
!    Input, integer LINE_PRUNE, pruning option.
!    0, no pruning, draw every line.
!    nonzero, prune.  Only draw the line from node I to node J if I < J.
!    This should cut down on repeated drawing of lines in the common
!    case of a face mesh where each line is drawn twice, once with positive
!    and once with negative orientation.  In other cases, pruning
!    may omit some lines that only occur with negative orientation.
!
!    Input, integer COR3_MAX, the maximum number of points.
!
!    Input, integer FACE_MAX, the maximum number of faces.
!
!    Input, integer LINE_MAX, the maximum number of line definition items.
!
!    Input, integer ORDER_MAX, the maximum number of vertices per face.
!
!    Input/output, integer LINE_NUM, the number of line data items.
!
!    Input, integer VERTEX_MAT(ORDER_MAX,FACE_MAX), vertex materials.
!
  integer, parameter :: OFFSET = 1

  integer face_max
  integer line_max
  integer order_max
!
  logical debug
  integer face(order_max,face_max)
  integer face_num
  integer face_order(face_max)
  integer i
  integer icor3
  integer iface
  integer ivert
  integer jcor3
  integer jvert
  integer line_dex(line_max)
  integer line_material(line_max)
  integer line_num
  integer line_prune
  integer vertex_material(order_max,face_max)
!
!  Case 1: 
!  No line pruning.
!
  if ( line_prune == 0 ) then

    do iface = 1, face_num

      do ivert = 1, face_order(iface)

        icor3 = face(ivert,iface)
 
        line_num = line_num + 1
        if ( line_num <= line_max ) then
          line_dex(line_num) = icor3
          line_material(line_num) = vertex_material(ivert,iface)
        end if

      end do

      ivert = 1
      icor3 = face(ivert,iface)

      line_num = line_num + 1
      if ( line_num <= line_max ) then
        line_dex(line_num) = icor3
        line_material(line_num) = vertex_material(ivert,iface)
      end if

      line_num = line_num + 1
      if ( line_num <= line_max ) then
        line_dex(line_num) = -1 + OFFSET
        line_material(line_num) = -1 + OFFSET
      end if

    end do
!
!  Case 2: 
!    Simple-minded line pruning.
!    Only draw line (I,J) if I < J.
!
  else

    do iface = 1, face_num

      do ivert = 1, face_order(iface)

        icor3 = face(ivert,iface)

        if ( ivert < face_order(iface) ) then
          jvert = ivert + 1
        else
          jvert = 1
        end if

        jcor3 = face(jvert,iface)

        if ( icor3 < jcor3 ) then

          if ( line_num + 3 <= line_max ) then

            line_num = line_num + 1
            line_dex(line_num) = icor3
            line_material(line_num) = vertex_material(ivert,iface)
 
            line_num = line_num + 1
            line_dex(line_num) = jcor3
            line_material(line_num) = vertex_material(jvert,iface)

            line_num = line_num + 1
            line_dex(line_num) = -1 + OFFSET
            line_material(line_num) = -1 + OFFSET

          end if

        end if

      end do

    end do

  end if

  if ( debug ) then

    write ( *, * ) ' '
    write ( *, * ) 'FACE_TO_LINE:'
    write ( *, * ) ' '
    write ( *, * ) 'I, LINE_DEX(I), LINE_MAT(I)'
    write ( *, * ) ' '

    do i = 1, line_num
      write ( *, '(i6,2x,i6,2x,i6)' ) i, line_dex(i), line_material(i)
    end do

  end if

  return
end
subroutine face_touch ( face, face_order, order_max, face_num, iface, jface, &
  touch )
!
!*******************************************************************************
!
!! FACE_TOUCH reports whether two polygonal faces touch.
!
!
!  Modified:
!
!    06 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer FACE(ORDER_MAX,FACE_NUM), the nodes making faces.
!
!    Input, integer FACE_ORDER(FACE_NUM), the number of nodes per face.
!
!    Input, integer ORDER_MAX, is the maximum number of nodes that can
!    make up a face, required to dimension FACE.
!
!    Input, integer FACE_NUM, the number of faces.
!
!    Input, integer IFACE, JFACE, the faces to be checked.
!
!    Output, integer TOUCH:
!     0, the faces do not touch;
!    +1, the faces touch, both using an arc in the same direction;
!    -1, the faces touch, using an arc in opposite directions.
!
  integer order_max
  integer face_num
!
  integer face(order_max,face_num)
  integer face_order(face_num)
  integer i
  integer iface
  integer j
  integer jface
  integer m
  integer mp1
  integer mm1
  integer n
  integer np1
  integer touch
!
  touch = 0
!
!  Arc N1-N2 on IFACE must be matched by arc N1-N2 or N2-N1 on JFACE.
!
  do i = 1, face_order(iface)

    n = face(i,iface)
    if ( i < face_order(iface) ) then
      np1 = face(i+1,iface)
    else
      np1 = face(1,iface)
    end if

    do j = 1, face_order(jface)

      m = face(j,jface)
      if ( j < face_order(jface) ) then
        mp1 = face(j+1,jface)
      else
        mp1 = face(1,jface)
      end if
      if ( j > 1 ) then
        mm1 = face(j-1,jface)
      else
        mm1 = face(face_order(jface),jface)
      end if

      if ( n == m ) then
        if ( np1 == mp1 ) then
          touch = + 1
          return
        else if ( np1 == mm1 ) then
          touch = - 1
          return
        end if
      end if

    end do
  end do

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
subroutine hello ( cor3_max, face_max, line_max, material_max, order_max, texture_max )
!
!*******************************************************************************
!
!! HELLO prints out a message about the program.
!
!
!  Modified:
!
!    12 April 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer COR3_MAX, the maximum number of points.
!
!    Input, integer FACE_MAX, the maximum number of faces.
!
!    Input, integer LINE_MAX, the maximum number of lines.
!
!    Input, integer MATERIAL_MAX, the maximum number of materials.
!
!    Input, integer ORDER_MAX, the maximum vertices per face.
!
!    Input, integer TEXTURE_MAX, the maximum number of textures.
!
  integer cor3_max
  integer face_max
  integer line_max
  integer material_max
  integer order_max
  integer texture_max
!
  write ( *, * ) ' '
  write ( *, * ) 'Hello:  This is IVRead,'
  write ( *, * ) '  a program which can convert some files from'
  write ( *, * ) '  some 3D graphics format to some others:'
  write ( *, * ) ' '
  write ( *, * ) '    ".ase"  3D Studio Max ASCII export;'
  write ( *, * ) '    ".byu"  Movie.BYU surface geometry;'
  write ( *, * ) '    ".dxf"  AutoCAD DXF;'
  write ( *, * ) '    ".hrc"  SoftImage hierarchy;'
  write ( *, * ) '    ".iv"   SGI Open Inventor;'
  write ( *, * ) '    ".obj"  WaveFront Advanced Visualizer;'
  write ( *, * ) '    ".oogl" OOGL file (input only);'  
  write ( *, * ) '    ".pov"  Persistence of Vision (output only);'
  write ( *, * ) '    ".ps"   PostScript (output only)(NOT READY);'
  write ( *, * ) '    ".smf"  Michael Garland''s format;'
  write ( *, * ) '    ".stl"  ASCII StereoLithography;'
  write ( *, * ) '    ".stla" ASCII StereoLithography;'
  write ( *, * ) '    ".tec"  TECPLOT (output only);'
  write ( *, * ) '    ".tri"  [Greg Hood triangles];'
  write ( *, * ) '    ".tria" [Greg Hood triangles];'
  write ( *, * ) '    ".txt"  Text (output only);'
  write ( *, * ) '    ".ucd"  AVS unstructured cell data (output only);'
  write ( *, * ) '    ".vla"  VLA;' 
  write ( *, * ) '    ".wrl"  VRML;'
  write ( *, * ) '    ".xgl"  XGL (output only) (DEVELOPMENT)'

  write ( *, * ) ' '
  write ( *, * ) '  Current limits:'
  write ( *, * ) ' '
  write ( *, * ) cor3_max,  ' points;'
  write ( *, * ) line_max,  ' line items;'
  write ( *, * ) face_max,  ' faces.'
  write ( *, * ) ' '
  write ( *, * ) order_max, ' vertices per face;'
  write ( *, * ) material_max,   ' materials;'
  write ( *, * ) texture_max,   ' textures.'
  write ( *, * ) ' '
  write ( *, * ) '  Last modifed: 14 June 2000.'
  write ( *, * ) ' '
  write ( *, * ) '  Send problem reports to burkardt@psc.edu.'
 
  return
end
subroutine help
!
!*******************************************************************************
!
!! HELP prints out a help message about the interactive commands.
!
!
!  Modified:
!
!    29 June 1999
!
!  Author:
!
!    John Burkardt
!
  write ( *, * ) ' '
  write ( *, * ) 'HELP:'
  write ( *, * ) '  Batch commands to convert IN_FILE to OUT_FILE:'
  write ( *, * ) ' '
  write ( *, * ) '  ivread  in_file  out_file'
  write ( *, * ) ' '
  write ( *, * ) '  ivread  -rn  in_file  out_file'
  write ( *, * ) '    Reverse normals before output.'
  write ( *, * ) ' '
  write ( *, * ) '  ivread  -rf  in_file  out_file'
  write ( *, * ) '    Reverse faces and normals before output.'
  write ( *, * ) ' '
  write ( *, * ) ' '
  write ( *, * ) 'HELP:'
  write ( *, * ) '  These are legal interactive commands:'
  write ( *, * ) ' '
  write ( *, * ) '  < in_file    Read data from a file.'
  write ( *, * ) '  << out_file  Append data from a file'
  write ( *, * ) '                 to current information.'
  write ( *, * ) '  > out_file   Write data out to a file.'
  write ( *, * ) '  B            Switch byte swapping option.'
  write ( *, * ) '  D            Switch debug option.'
  write ( *, * ) '  F            Print info about one face.'
  write ( *, * ) '  H            Print this help message.'
  write ( *, * ) '  I            Info, print out recent changes.'
  write ( *, * ) '  LINE_PRUNE   Set FACE_TO_LINE pruning option.'
  write ( *, * ) '  LINES        Convert faces to lines.'
  write ( *, * ) '  N            Recompute normal vectors.'
  write ( *, * ) '  O            Use an average node normal.'
  write ( *, * ) '  Q            Quit.'
  write ( *, * ) '  REVERSE      Reverse the normal vectors.'
  write ( *, * ) '  RELAX        Smooth surface via relaxation.'
  write ( *, * ) '  S            Select face subset.'
  write ( *, * ) '  T            Transform data.'
  write ( *, * ) '  U            Renumber faces and analyze.'
  write ( *, * ) '  V            Convert polygons to triangles.'
  write ( *, * ) '  W            Reverse faces and normals.'
  return
end
subroutine hrc_read ( cor3, face, face_material, face_order, ierror, iunit, &
  line_dex, line_material, material_name, material_rgba, cor3_max, face_max, &
  line_max, material_max, order_max, texture_max, bad_num, cor3_num, dup_num, &
  face_num, line_num, material_num, texture_num, text_num, texture_name, vertex_material, &
  vertex_normal, vertex_tex_uv )
!
!*******************************************************************************
!
!! HRC_READ reads graphics information from a SoftImage HRC file.
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
!    HRCH: Softimage 4D Creative Environment v3.00
!
!
!    model
!    {
!      name         "cube_10x10"
!      scaling      1.000 1.000 1.000
!      rotation     0.000 0.000 0.000
!      translation  0.000 0.000 0.000
!
!      mesh
!      {
!        flag    ( PROCESS )
!        discontinuity  60.000
!
!        vertices   8
!        {
!          [0] position  -5.000  -5.000  -5.000
!          [1] position  -5.000  -5.000  5.000
!          [2] position  -5.000  5.000  -5.000
!          [3] position  -5.000  5.000  5.000
!          [4] position  5.000  -5.000  -5.000
!          [5] position  5.000  -5.000  5.000
!          [6] position  5.000  5.000  -5.000
!          [7] position  5.000  5.000  5.000
!        }
!
!        polygons   6
!        {
!          [0] nodes  4
!              {
!                [0] vertex  0
!                    normal  -1.000  0.000  0.000
!                    uvTexture  0.000  0.000
!                    vertexColor 255 178 178 178
!                [1] vertex  1
!                    normal  -1.000  0.000  0.000
!                    uvTexture  0.000  0.000
!                    vertexColor 255 178 178 178
!                [2] vertex  3
!                    normal  -1.000  0.000  0.000
!                    uvTexture  0.000  0.000
!                    vertexColor 255 178 178 178
!                [3] vertex  2
!                    normal  -1.000  0.000  0.000
!                    uvTexture  0.000  0.000
!                    vertexColor 255 178 178 178
!              }
!              material  0
!          [1] nodes  4
!             {
!                [0] vertex  1
!                    normal  0.000  0.000  1.000
!                    uvTexture  0.000  0.000
!                    vertexColor 255 178 178 178
!                [1] vertex  5
!
!    ...etc.....
!
!          [5] nodes  4
!              {
!                [0] vertex  2
!                    normal  0.000  1.000  0.000
!                    uvTexture  0.000  0.000
!                    vertexColor 255 178 178 178
!                [1] vertex  3
!                    normal  0.000  1.000  0.000
!                    uvTexture  0.000  0.000
!                    vertexColor 255 178 178 178
!                [2] vertex  7
!                    normal  0.000  1.000  0.000
!                    uvTexture  0.000  0.000
!                    vertexColor 255 178 178 178
!                [3] vertex  6
!                    normal  0.000  1.000  0.000
!                    uvTexture  0.000  0.000
!                    vertexColor 255 178 178 178
!              }
!              material  0
!        }
!
!        edges   12
!        {
!          [1] vertices  3  2
!          [2] vertices  2  0
!          [3] vertices  0  1
!          [4] vertices  1  3
!          [5] vertices  7  3
!          [6] vertices  1  5
!          [7] vertices  5  7
!          [8] vertices  6  7
!          [9] vertices  5  4
!          [10] vertices  4  6
!          [11] vertices  2  6
!          [12] vertices  4  0
!        }
!      }
!
!      material [0]
!      {
!      name           "kazoo"
!      type           PHONG
!      ambient        0.0  1.0  0.0
!      diffuse        1.0  0.0  0.0
!      specular       0.0  0.0  1.0
!      exponent      50.0
!      reflectivity   0.0
!      transparency   0.0
!      refracIndex    1.0
!      glow           0
!      coc            0.0
!      }
!
!      texture [0]
!      {
!      name          "/usr/users/foss/HOUSE/PICTURES/mellon"
!      glbname       "t2d1"
!      anim          STATIC
!      method        XY
!      repeat        1  1
!      scaling       1.000  1.000
!      offset        0.000  0.000
!      pixelInterp
!      effect        INTENSITY
!      blending      1.000
!      ambient       0.977
!      diffuse       1.000
!      specular      0.966
!      reflect       0.000
!      transp        0.000
!      roughness     0.000
!      reflMap       1.000
!      rotation      0.000
!      txtsup_rot    0.000  0.000  0.000
!      txtsup_trans  0.000  0.000  0.000
!      txtsup_scal   1.000  1.000  1.000
!      }
!
!    }
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
!    Input/output, integer LINE_MAT(LINE_MAX), material index for line.
!
!    Input/output, character ( len = 100 ) MATERIAL_NAME(MATERIAL_MAX), material names.
!
!    Input/output, real MATERIAL_RGBA(4,MATERIAL_MAX), material R, G, B and A values.
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
!    Output, integer BAD_NUM, the number of bad text lines.
!
!    Input/output, integer COR3_NUM, the number of points.
!
!    Output, integer DUP_NUM, the number of duplicate points that were dropped.
!
!    Input/output, integer FACE_NUM, the number of faces.
!
!    Input/output, integer LINE_NUM, the number of line definition items.
!
!    Input/output, integer MATERIAL_NUM, the number of materials.
!
!    Input/output, integer TEXTURE_NUM, the number of textures.
!
!    Output, integer TEXT_NUM, the number of text lines.
!
!    Input/output, character ( len = 100 ) TEXTURE_NAME(TEXTURE_MAX), texture names.
!
!    Input/output, integer VERTEX_MAT(ORDER_MAX,FACE_MAX), vertex materials.
!
!    Input/output, real VERTEX_NORMAL(3,ORDER_MAX,FACE_MAX), normals at vertices.
!
!    Input/output, real VERTEX_TEX_UV(2,ORDER_MAX,FACE_MAX), vertex texture coordinates.
!
  logical, parameter :: debug = .FALSE.
  integer, parameter :: level_max = 10
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
  character ( len = 4 ) char4
  real cor3(3,cor3_max)
  integer cor3_num
  integer cor3_num_old
  logical done
  integer dup_num
  integer face(order_max,face_max)
  integer face_material(face_max)
  integer face_num
  integer face_order(face_max)
  integer i
  integer icor3
  integer ierror
  integer iunit
  integer ival1
  integer ival2
  integer ival3
  integer ival4
  integer ivert
  integer iword
  integer jval
  integer lchar
  integer level
  character ( len = 256 ) level_name(0:level_max)
  integer line_dex(line_max)
  integer line_material(line_max)
  integer line_num
  logical lval
  character ( len = 100 ) material_name(material_max)
  integer material_num
  integer material_num_old
  real material_rgba(4,material_max)
  integer nlbrack
  integer nrbrack
  real rval
  logical s_eqi
  logical s_is_i
  real temp(3)
  character ( len = 256 ) text
  integer text_num
  character ( len = 100 ) texture_name(texture_max)
  integer texture_num
  integer vertex_material(order_max,face_max)
  real vertex_normal(3,order_max,face_max)
  real vertex_tex_uv(2,order_max,face_max)
  character ( len = 256 ) word
  character ( len = 256 ) word2
  character ( len = 256 ) wordm1
  real x
  real y
  real z
!
  ierror = 0
  ival1 = 0
  ival2 = 0
  ival3 = 0
  ival4 = 0
  jval = 0
  level = 0
  level_name(0) = 'Top'
  nlbrack = 0
  nrbrack = 0
  cor3_num_old = cor3_num
  material_num_old = material_num
  word = ' '
  wordm1 = ' '
!
!  Read a line of text from the file.
!
10    continue

  read ( iunit, '(a)', end = 50 ) text
  text_num = text_num + 1

  if ( text == ' ' ) then
    go to 10
  end if
!
!  The first line of the file must be the header.
!
  if ( text_num == 1 ) then

    if ( .not. s_eqi ( text(1:5), 'HRCH:' ) ) then
      ierror = 1
      write ( *, * ) ' '
      write ( *, * ) 'HRC_READ - Fatal error!'
      write ( *, * ) '  The input file has a bad header.'
      write ( *, '(a)' ) trim ( text )
      return
    else
      go to 10
    end if

  end if

  done = .true.
  iword = 0
!
!  Save the previous word read.  
!  It helps when a word depends on its cotext_num.
!
20    continue

  if ( word /= ' ' .and. word .ne. ',' ) then
    wordm1 = word
  end if
!
!  Read a word from the line.
!
  call word_nexrd ( text, word, done )
!
!  If no more words in this line, read in a whole new line.
!
  if ( done ) then
    go to 10
  end if
!
!  Ignore blanks and commas.
!
  if ( word == ' ' .or. word == ',' ) then
    go to 20
  end if
!
!  Count the words in the current line, and the total.
!
  iword = iword + 1
!
!  If the word is a curly bracket, count it.
!
  if ( word == '{' ) then

    nlbrack = nlbrack + 1

  else if ( word .eq. '}' ) then

    nrbrack = nrbrack + 1

    if ( nlbrack < nrbrack ) then
      write ( *, * ) ' '
      write ( *, * ) 'HRC_READ - Fatal error!'
      write ( *, * ) '  Extraneous right bracket, line ', text_num
      write ( *, '(a)' ) trim ( text )
      write ( *, * ) 'Currently processing field:'
      write ( *, '(a)' ) trim ( level_name(level) )
      ierror = 1
      return
    end if

  end if
!
!  If the word is a left bracket, the previous word is the name of a node.
!
  if ( word == '{' ) then

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
!  CONTROLPOINTS
!
  if ( s_eqi ( level_name(level), 'CONTROLPOINTS' ) ) then

    if ( word == '{' ) then

    else if ( word == '}' ) then

      line_num = line_num + 1
      if ( line_num <= line_max ) then
        line_dex(line_num) = -1 + OFFSET
        line_material(line_num) = material_num
      end if

      level = nlbrack - nrbrack

    else if ( word(1:1) == '[' ) then

      call word_nexrd ( text, word2, done )
      call word_nexrd ( text, word2, done )

    else if ( s_eqi ( word, 'POSITION' ) ) then

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
        if ( cor3_num <= cor3_max ) then
          cor3(1,cor3_num) = temp(1)
          cor3(2,cor3_num) = temp(2)
          cor3(3,cor3_num) = temp(3)
        end if

        icor3 = cor3_num 

      else

        dup_num = dup_num + 1

      end if

      line_num = line_num + 1
      if ( line_num <= line_max ) then
        line_dex(line_num) = icor3
        line_material(line_num) = material_num
      end if

    else

      go to 99

    end if
!
!  EDGES
!
  else if ( s_eqi ( level_name(level), 'EDGES' ) ) then

    if ( word == '{' ) then

    else if ( word == '}' ) then

      level = nlbrack - nrbrack

    else if ( word(1:1) == '[' ) then

      call word_nexrd ( text, word2, done )
      call word_nexrd ( text, word2, done )

    else if ( s_eqi ( word, 'VERTICES' ) ) then

      call word_nexrd ( text, word2, done )
      lval = s_is_i ( word2, jval )
      line_num = line_num + 1

      if ( line_num <= line_max ) then
        line_dex(line_num) = jval + cor3_num_old + OFFSET
        line_material(line_num) = material_num
      end if

      call word_nexrd ( text, word2, done )
      lval = s_is_i ( word2, jval )
      line_num = line_num + 1

      if ( line_num <= line_max ) then
        line_dex(line_num) = jval + cor3_num_old + OFFSET
        line_material(line_num) = material_num
      end if

      line_num = line_num + 1

      if ( line_num <= line_max ) then
        line_dex(line_num) = -1 + OFFSET
        line_material(line_num) = -1 + OFFSET
      end if

    else

      go to 99

    end if
!
!  MATERIAL
!
  else if ( s_eqi ( level_name(level), 'MATERIAL' ) ) then

    if ( word == '{' ) then

      material_num = material_num + 1
      if ( material_num <= material_max ) then
        call i_to_s_zero ( material_num, char4 )
        material_name(material_num) = 'Material_' // char4
      end if

    else if ( word == '}' ) then

      level = nlbrack - nrbrack

    else if ( word(1:1) == '[' ) then

    else if ( s_eqi ( word, 'AMBIENT' ) ) then

      call word_nexrd ( text, word2, done )
      call word_nexrd ( text, word2, done )
      call word_nexrd ( text, word2, done )

    else if ( s_eqi ( word, 'COC' ) ) then

      call word_nexrd ( text, word2, done )

    else if ( s_eqi ( word, 'DIFFUSE' ) ) then

      call word_nexrd ( text, word2, done )
      call s_to_r ( word2, rval, ierror, lchar )
      material_rgba(1,material_num) = rval

      call word_nexrd ( text, word2, done )
      call s_to_r ( word2, rval, ierror, lchar )
      material_rgba(2,material_num) = rval

      call word_nexrd ( text, word2, done )
      call s_to_r ( word2, rval, ierror, lchar )
      material_rgba(3,material_num) = rval

    else if ( s_eqi ( word, 'EXPONENT' ) ) then

      call word_nexrd ( text, word2, done )

    else if ( s_eqi ( word, 'GLOW' ) ) then

      call word_nexrd ( text, word2, done )

    else if ( s_eqi ( word, 'NAME' ) ) then

      call word_nexrd ( text, word2, done )
      call word_nexrd ( text, word2, done )
      material_name(material_num) = word2
      call word_nexrd ( text, word2, done )

    else if ( s_eqi ( word, 'REFLECTIVITY' ) ) then

      call word_nexrd ( text, word2, done )

    else if ( s_eqi ( word, 'REFRACINDEX' ) ) then

      call word_nexrd ( text, word2, done )

    else if ( s_eqi ( word, 'SPECULAR' ) ) then

      call word_nexrd ( text, word2, done )
      call word_nexrd ( text, word2, done )
      call word_nexrd ( text, word2, done )

    else if ( s_eqi ( word, 'TRANSPARENCY' ) ) then

      call word_nexrd ( text, word2, done )
      call s_to_r ( word2, rval, ierror, lchar )
      material_rgba(4,material_num) = 1.0 - rval

    else if ( s_eqi ( word, 'TYPE' ) ) then

      call word_nexrd ( text, word2, done )

    else

      go to 99

    end if
!
!  MESH
!
  else if ( s_eqi ( level_name(level), 'MESH' ) ) then

    if ( word == '{' ) then

    else if ( word == '}' ) then
      level = nlbrack - nrbrack
    else if ( s_eqi ( word, 'DISCONTINUITY' ) ) then
      go to 10
    else if ( s_eqi ( word, 'EDGES' ) ) then
      go to 10
    else if ( s_eqi ( word, 'FLAG' ) ) then
      go to 10
    else if ( s_eqi ( word, 'POLYGONS' ) ) then
      go to 10
    else if ( s_eqi ( word, 'VERTICES' ) ) then
      go to 10
    else
      go to 99
    end if
!
!  MODEL
!
  else if ( s_eqi ( level_name(level), 'MODEL' ) ) then

    if ( word == '{' ) then

    else if ( word == '}' ) then

      level = nlbrack - nrbrack

    else if ( s_eqi ( word, 'MATERIAL' ) ) then
      go to 10
    else if ( s_eqi ( word, 'MESH' ) ) then
      go to 10
    else if ( s_eqi ( word, 'NAME' ) ) then
      go to 10
    else if ( s_eqi ( word, 'PATCH' ) ) then
      go to 10
    else if ( s_eqi ( word, 'ROTATION' ) ) then
      go to 10
    else if ( s_eqi ( word, 'SCALING' ) ) then
      go to 10
    else if ( s_eqi ( word, 'SPLINE' ) ) then
      go to 10
    else if ( s_eqi ( word, 'TRANSLATION' ) ) then
      go to 10
    else
      go to 99
    end if
!
!  NODES
!
  else if ( s_eqi ( level_name(level), 'NODES' ) ) then

    if ( word == '{' ) then

      face_num = face_num + 1
      ivert = 0
      face_order(face_num) = 0

    else if ( word == '}' ) then

      level = nlbrack - nrbrack

    else if ( word(1:1) == '[' ) then

      call word_nexrd ( text, word2, done )
      call word_nexrd ( text, word2, done )

    else if ( s_eqi ( word, 'NORMAL' ) ) then

      call word_nexrd ( text, word2, done )
      call s_is_r ( word2, x, lval )

      call word_nexrd ( text, word2, done )
      call s_is_r ( word2, y, lval )

      call word_nexrd ( text, word2, done )
      call s_is_r ( word2, z, lval )

      if ( face_num <= face_max .and. ivert <= order_max ) then
        vertex_normal(1,ivert,face_num) = x
        vertex_normal(2,ivert,face_num) = y
        vertex_normal(3,ivert,face_num) = z
      end if

    else if ( s_eqi ( word, 'UVTEXTURE' ) ) then

      call word_nexrd ( text, word2, done )
      call s_is_r ( word2, x, lval )

      call word_nexrd ( text, word2, done )
      call s_is_r ( word2, y, lval )

      if ( face_num <= face_max .and. ivert <= order_max ) then
        vertex_tex_uv(1,ivert,face_num) = x
        vertex_tex_uv(2,ivert,face_num) = y
      end if

    else if ( s_eqi ( word, 'VERTEX' ) ) then

      call word_nexrd ( text, word2, done )
      lval = s_is_i ( word2, jval )
      ivert = ivert + 1

      if ( ivert <= order_max .and. face_num <= face_max ) then
        face_order(face_num) = face_order(face_num) + 1
        face(ivert,face_num) = jval + cor3_num_old + OFFSET
      end if
!
!  What do I do with this?  Define a vertex material?
!
    else if ( s_eqi ( word, 'VERTEXCOLOR' ) ) then

      call word_nexrd ( text, word2, done )
      lval = s_is_i ( word2, ival1 )

      call word_nexrd ( text, word2, done )
      lval = s_is_i ( word2, ival2 )

      call word_nexrd ( text, word2, done )
      lval = s_is_i ( word2, ival3 )

      call word_nexrd ( text, word2, done )
      lval = s_is_i ( word2, ival4 )

    else
      go to 99
    end if
!
!  PATCH
!
!  JVB: I don't know what to do with this yet.
!
  else if ( s_eqi ( level_name(level), 'PATCH' ) ) then

    if ( word == '{' ) then

    else if ( word == '}' ) then
      level = nlbrack - nrbrack
    else if ( s_eqi ( word, 'APPROX_TYPE' ) ) then
      call word_nexrd ( text, word2, done )
    else if ( s_eqi ( word, 'CONTROLPOINTS' ) ) then
      go to 10
    else if ( s_eqi ( word, 'CURV_U' ) ) then
      call word_nexrd ( text, word2, done )
    else if ( s_eqi ( word, 'CURV_V' ) ) then
      call word_nexrd ( text, word2, done )
    else if ( s_eqi ( word, 'RECMIN' ) ) then
      call word_nexrd ( text, word2, done )
    else if ( s_eqi ( word, 'RECMAX' ) ) then
      call word_nexrd ( text, word2, done )
    else if ( s_eqi ( word, 'RECURSION' ) ) then
      call word_nexrd ( text, word2, done )
    else if ( s_eqi ( word, 'SPACIAL' ) ) then
      call word_nexrd ( text, word2, done )
    else if ( s_eqi ( word, 'TAGGEDPOINTS' ) ) then
      go to 10
    else if ( s_eqi ( word, 'UCURVE' ) ) then
      call word_nexrd ( text, word2, done )
    else if ( s_eqi ( word, 'UPOINT' ) ) then
      call word_nexrd ( text, word2, done )
    else if ( s_eqi ( word, 'USTEP' ) ) then
      call word_nexrd ( text, word2, done )
    else if ( s_eqi ( word, 'UTENSION' ) ) then
      call word_nexrd ( text, word2, done )
    else if ( s_eqi ( word, 'UTYPE' ) ) then
      call word_nexrd ( text, word2, done )
    else if ( s_eqi ( word, 'VCLOSE' ) ) then

    else if ( s_eqi ( word, 'VCURVE' ) ) then
      call word_nexrd ( text, word2, done )
    else if ( s_eqi ( word, 'VIEWDEP' ) ) then
      call word_nexrd ( text, word2, done )
    else if ( s_eqi ( word, 'VPOINT' ) ) then
      call word_nexrd ( text, word2, done )
    else if ( s_eqi ( word, 'VSTEP' ) ) then
      call word_nexrd ( text, word2, done )
    else if ( s_eqi ( word, 'VTENSION' ) ) then
      call word_nexrd ( text, word2, done )
    else if ( s_eqi ( word, 'VTYPE' ) ) then
      call word_nexrd ( text, word2, done )
    else
      go to 99
    end if
!
!  POLYGONS
!
  else if ( s_eqi ( level_name(level), 'POLYGONS' ) ) then

    if ( word == '{' ) then

    else if ( word == '}' ) then

      level = nlbrack - nrbrack

    else if ( word == '[' ) then

      call word_nexrd ( text, word2, done )
      call word_nexrd ( text, word2, done )

    else if ( s_eqi ( word, 'MATERIAL' ) ) then

      call word_nexrd ( text, word2, done )
      lval = s_is_i ( word2, jval )

      face_material(face_num) = jval + material_num_old + OFFSET

      do i = 1, order_max
        vertex_material(i,face_num) = jval + material_num_old + OFFSET
      end do

    else if ( s_eqi ( word, 'NODES' ) ) then
      call word_nexrd ( text, word2, done )
    else
      go to 99
    end if
!
!  SPLINE
!
  else if ( s_eqi ( level_name(level), 'SPLINE' ) ) then

    if ( word == '{' ) then

    else if ( word == '}' ) then

      level = nlbrack - nrbrack

    else if ( s_eqi ( word, 'CONTROLPOINTS' ) ) then
      go to 10
    else if ( s_eqi ( word, 'NAME' ) ) then
      go to 10
    else if ( s_eqi ( word, 'NBKEYS' ) ) then
      go to 10
    else if ( s_eqi ( word, 'STEP' ) ) then
      go to 10
    else if ( s_eqi ( word, 'TENSION' ) ) then
      go to 10
    else if ( s_eqi ( word, 'TYPE' ) ) then
      go to 10
    else
      go to 99
    end if
!
!  TAGGEDPOINTS
!
  else if ( s_eqi ( level_name(level), 'TAGGEDPOINTS' ) ) then

    if ( word == '{' ) then

    else if ( word == '}' ) then

      level = nlbrack - nrbrack

    else if ( word(1:1) == '[' ) then

      call word_nexrd ( text, word2, done )
      call word_nexrd ( text, word2, done )

    else if ( s_eqi ( word, 'TAGGED' ) ) then
      call word_nexrd ( text, word2, done )
    else
      go to 99
    end if
!
!  TEXTURE
!
  else if ( s_eqi ( level_name(level), 'TEXTURE' ) ) then

    if ( word == '{' ) then

      texture_num = texture_num + 1

      if ( texture_num <= texture_max ) then
        call i_to_s_zero ( texture_num, char4 )
        texture_name(texture_num) = 'Texture_' // char4
      end if

    else if ( word == '}' ) then

      level = nlbrack - nrbrack

    else if ( word(1:1) == '[' ) then

    else if ( s_eqi ( word, 'AMBIENT' ) ) then

      call word_nexrd ( text, word2, done )

    else if ( s_eqi ( word, 'ANIM' ) ) then

      call word_nexrd ( text, word2, done )

    else if ( s_eqi ( word, 'BLENDING' ) ) then

      call word_nexrd ( text, word2, done )

    else if ( s_eqi ( word, 'DIFFUSE' ) ) then

      call word_nexrd ( text, word2, done )

    else if ( s_eqi ( word, 'EFFECT' ) ) then

      call word_nexrd ( text, word2, done )

    else if ( s_eqi ( word, 'GLBNAME' ) ) then

      call word_nexrd ( text, word2, done )

    else if ( s_eqi ( word, 'METHOD' ) ) then

      call word_nexrd ( text, word2, done )
!
!  (I assume there are initial and trailing quotes in the NAME field,
!  which are treated as separate words.)
!
    else if ( s_eqi ( word, 'NAME' ) ) then

      call word_nexrd ( text, word2, done )
      call word_nexrd ( text, word2, done )
      texture_name(texture_num) = word2
      call word_nexrd ( text, word2, done )

    else if ( s_eqi ( word, 'OFFSET' ) ) then

      call word_nexrd ( text, word2, done )
      call word_nexrd ( text, word2, done )

    else if ( s_eqi ( word, 'PIXELINTERP' ) ) then

    else if ( s_eqi ( word, 'REFLECT' ) ) then

      call word_nexrd ( text, word2, done )

    else if ( s_eqi ( word, 'REFLMAP' ) ) then

      call word_nexrd ( text, word2, done )

    else if ( s_eqi ( word, 'REPEAT' ) ) then

      call word_nexrd ( text, word2, done )
      call word_nexrd ( text, word2, done )

    else if ( s_eqi ( word, 'ROTATION' ) ) then

      call word_nexrd ( text, word2, done )

    else if ( s_eqi ( word, 'ROUGHNESS' ) ) then

      call word_nexrd ( text, word2, done )

    else if ( s_eqi ( word, 'SCALING' ) ) then

      call word_nexrd ( text, word2, done )
      call word_nexrd ( text, word2, done )

    else if ( s_eqi ( word, 'SPECULAR' ) ) then

      call word_nexrd ( text, word2, done )

    else if ( s_eqi ( word, 'TRANSP' ) ) then

      call word_nexrd ( text, word2, done )

    else if ( s_eqi ( word, 'TXTSUP_ROT' ) ) then

      call word_nexrd ( text, word2, done )
      call word_nexrd ( text, word2, done )
      call word_nexrd ( text, word2, done )

    else if ( s_eqi ( word, 'TXTSUP_SCAL' ) ) then

      call word_nexrd ( text, word2, done )
      call word_nexrd ( text, word2, done )
      call word_nexrd ( text, word2, done )

    else if ( s_eqi ( word, 'TXTSUP_TRANS' ) ) then

      call word_nexrd ( text, word2, done )
      call word_nexrd ( text, word2, done )
      call word_nexrd ( text, word2, done )

    else

      go to 99

    end if
!
!  VERTICES
!
  else if ( s_eqi ( level_name(level), 'VERTICES' ) ) then

    if ( word == '{' ) then

    else if ( word == '}' ) then

      level = nlbrack - nrbrack

    else if ( word(1:1) == '[' ) then

      call word_nexrd ( text, word2, done )
      call word_nexrd ( text, word2, done )

    else if ( s_eqi ( word, 'POSITION' ) ) then

      cor3_num = cor3_num + 1

      call word_nexrd ( text, word2, done )
      call s_is_r ( word2, x, lval )

      call word_nexrd ( text, word2, done )
      call s_is_r ( word2, y, lval )

      call word_nexrd ( text, word2, done )
      call s_is_r ( word2, z, lval )

      if ( cor3_num <= cor3_max ) then
        cor3(1,cor3_num) = x
        cor3(2,cor3_num) = y
        cor3(3,cor3_num) = z
      end if

    else
      go to 99
    end if
!
!  Any other word:
!
  else

  end if

  go to 20
!
!  Bad data
!
99    continue

  bad_num = bad_num + 1

  if ( bad_num <= 10 ) then
    write ( *, * ) ' '
    write ( *, * ) 'HRC_READ - Warning!'
    write ( *, * ) '  Bad data on level ' // trim ( level_name(level) )
    write ( *, '(a)' ) '  Bad word: ' // trim ( word )
    write ( *, * ) '  Line number: ', text_num
    write ( *, '(a)' ) trim ( text )
  else
    write ( *, * ) ' '
    write ( *, * ) 'HRC_READ - Fatal error!'
    write ( *, * ) '  Too many warnings!'
    return
  end if

  go to 10
!
!  Normal end of information in file.
!
50    continue
!
!  Check the "materials" defining a line.
!
!  If COORDINDEX is -1, so should be the MATERIALINDEX.
!  If COORDINDEX is not -1, then the MATERIALINDEX shouldn't be either.
!
  do i = 1, line_num

    if ( line_dex(i) == -1 + OFFSET ) then
      line_material(i) = -1 + OFFSET
    else if ( line_material(i) == -1 + OFFSET ) then
      line_material(i) = material_num
    end if

  end do

  return
end
subroutine hrc_write ( cor3, face, face_material, face_order, fileout_name, &
  iunit, line_dex, material_name, material_rgba, cor3_max, face_max, &
  line_max, material_max, order_max, texture_max, cor3_num, face_num, line_num, &
  material_num, texture_num, texture_name, vertex_normal, vertex_tex_uv )
!
!*******************************************************************************
!
!! HRC_WRITE writes graphics data to an HRC SoftImage file.
!
!
!  Example:
!
!    HRCH: Softimage 4D Creative Environment v3.00
!
!
!    model
!    {
!      name         "cube_10x10"
!      scaling      1.000 1.000 1.000
!      rotation     0.000 0.000 0.000
!      translation  0.000 0.000 0.000
!
!      mesh
!      {
!        flag    ( PROCESS )
!        discontinuity  60.000
!
!        vertices   8
!        {
!          [0] position  -5.000  -5.000  -5.000
!          [1] position  -5.000  -5.000  5.000
!          [2] position  -5.000  5.000  -5.000
!          [3] position  -5.000  5.000  5.000
!          [4] position  5.000  -5.000  -5.000
!          [5] position  5.000  -5.000  5.000
!          [6] position  5.000  5.000  -5.000
!          [7] position  5.000  5.000  5.000
!        }
!
!        polygons   6
!        {
!          [0] nodes  4
!              {
!                [0] vertex  0
!                    normal  -1.000  0.000  0.000
!                    uvTexture  0.000  0.000
!                    vertexColor 255 178 178 178
!                [1] vertex  1
!                    normal  -1.000  0.000  0.000
!                    uvTexture  0.000  0.000
!                    vertexColor 255 178 178 178
!                [2] vertex  3
!                    normal  -1.000  0.000  0.000
!                    uvTexture  0.000  0.000
!                    vertexColor 255 178 178 178
!                [3] vertex  2
!                    normal  -1.000  0.000  0.000
!                    uvTexture  0.000  0.000
!                    vertexColor 255 178 178 178
!              }
!              material  0
!          [1] nodes  4
!             {
!                [0] vertex  1
!                    normal  0.000  0.000  1.000
!                    uvTexture  0.000  0.000
!                    vertexColor 255 178 178 178
!                [1] vertex  5
!
!    ...etc.....
!
!          [5] nodes  4
!              {
!                [0] vertex  2
!                    normal  0.000  1.000  0.000
!                    uvTexture  0.000  0.000
!                    vertexColor 255 178 178 178
!                [1] vertex  3
!                    normal  0.000  1.000  0.000
!                    uvTexture  0.000  0.000
!                    vertexColor 255 178 178 178
!                [2] vertex  7
!                    normal  0.000  1.000  0.000
!                    uvTexture  0.000  0.000
!                    vertexColor 255 178 178 178
!                [3] vertex  6
!                    normal  0.000  1.000  0.000
!                    uvTexture  0.000  0.000
!                    vertexColor 255 178 178 178
!              }
!              material  0
!        }
!
!        edges   12
!        {
!          [1] vertices  3  2
!          [2] vertices  2  0
!          [3] vertices  0  1
!          [4] vertices  1  3
!          [5] vertices  7  3
!          [6] vertices  1  5
!          [7] vertices  5  7
!          [8] vertices  6  7
!          [9] vertices  5  4
!          [10] vertices  4  6
!          [11] vertices  2  6
!          [12] vertices  4  0
!        }
!      }
!
!      material [0]
!      {
!      name           "kazoo"
!      type           PHONG
!      ambient        0.0  1.0  0.0
!      diffuse        1.0  0.0  0.0
!      specular       0.0  0.0  1.0
!      exponent      50.0
!      reflectivity   0.0
!      transparency   0.0
!      refracIndex    1.0
!      glow           0
!      coc            0.0
!      }
!
!      texture [0]
!      {
!      name          "/usr/users/foss/HOUSE/PICTURES/mellon"
!      glbname       "t2d1"
!      anim          STATIC
!      method        XY
!      repeat        1  1
!      scaling       1.000  1.000
!      offset        0.000  0.000
!      pixelInterp
!      effect        INTENSITY
!      blending      1.000
!      ambient       0.977
!      diffuse       1.000
!      specular      0.966
!      reflect       0.000
!      transp        0.000
!      roughness     0.000
!      reflMap       1.000
!      rotation      0.000
!      txtsup_rot    0.000  0.000  0.000
!      txtsup_trans  0.000  0.000  0.000
!      txtsup_scal   1.000  1.000  1.000
!      }
!
!    }
!
!  Modified:
!
!    24 June 1999
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
!    Input, integer FACE_MATERIAL(FACE_MAX), face materials.
!
!    Input, integer FACE_ORDER(FACE_MAX), the number of vertices per face.
!
!    Input, character ( len = 100 ) FILEOUT_NAME, the name of the output file.
!
!    Input, integer IUNIT, the FORTRAN unit to which output is written.
!
!    Input, integer LINE_DEX(LINE_MAX), nodes forming a line, terminated by -1.
!
!    Input, character ( len = 100 ) MATERIAL_NAME(MATERIAL_MAX), material names.
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
!    Input, integer COR3_NUM, the number of points.
!
!    Input, integer FACE_NUM, the number of faces.
!
!    Input, integer LINE_NUM, the number of line definition items.
!
!    Input, integer MATERIAL_NUM, the number of materials.
!
!    Input, integer TEXTURE_NUM, the number of textures.
!
!    Input, character ( len = 100 ) TEXTURE_NAME(TEXTURE_MAX), texture names.
!
!    Input, real VERTEX_NORMAL(3,ORDER_MAX,FACE_MAX), normals at vertices.  
!
!    Input, real VERTEX_TEX_UV(2,ORDER_MAX,FACE_MAX), vertex texture coordinates.
!
  integer, parameter :: OFFSET = 1
!
  integer cor3_max
  integer face_max
  integer line_max
  integer material_max
  integer order_max
  integer texture_max
!
  real cor3(3,cor3_max)
  integer cor3_num
  integer face(order_max,face_max)
  integer face_material(face_max)
  integer face_num
  integer face_order(face_max)
  character ( len = 100 ) fileout_name
  integer i
  integer iface
  integer iunit
  integer ivert
  integer j
  integer jhi
  integer jlo
  integer jrel
  integer k
  integer line_dex(line_max)
  integer line_num
  character ( len = 100 ) material_name(material_max)
  integer material_num
  real material_rgba(4,material_max)
  integer npts
  integer nseg
  character ( len = 100 ) text
  integer text_num
  character ( len = 100 ) texture_name(texture_max)
  integer texture_num
  real vertex_normal(3,order_max,face_max)
  real vertex_tex_uv(2,order_max,face_max)
  character ( len = 10 ) word
  real x
  real y
  real z
!
  nseg = 0
  text_num = 0

  write ( iunit, '(a)' ) 'HRCH: Softimage 4D Creative Environment v3.00'
  write ( iunit, '(a)' ) ' '
  write ( iunit, '(a)' ) ' '
  text_num = text_num + 3

  write ( iunit, '(a)' ) 'model'
  write ( iunit, '(a)' ) '{'
  write ( iunit, '(a)' ) '  name         "' // trim ( fileout_name ) // '"'
  write ( iunit, '(a)' ) '  scaling      1.000 1.000 1.000'
  write ( iunit, '(a)' ) '  rotation     0.000 0.000 0.000'
  write ( iunit, '(a)' ) '  translation  0.000 0.000 0.000'
  text_num = text_num + 6

  if ( face_num > 0 ) then

    write ( iunit, '(a)' ) ' '
    write ( iunit, '(a)' ) '  mesh'
    write ( iunit, '(a)' ) '  {'
    write ( iunit, '(a)' ) '    flag    ( PROCESS )'
    write ( iunit, '(a)' ) '    discontinuity  60.000'
    text_num = text_num + 5
!
!  Point coordinates.
!
    if ( cor3_num > 0 ) then

      write ( iunit, '(a)' ) ' '
      write ( text, '(a, i8)' ) 'vertices ', cor3_num
      call s_blanks_delete ( text )
      write ( iunit, '(4x,a)' ) trim ( text )
      write ( iunit, '(a)' )     '    {'
      text_num = text_num + 3
 
      do j = 1, cor3_num

        write ( word, '( ''['', i8, '']'' )' ) j-OFFSET
        call s_blank_delete ( word )

        write ( text, '(a,'' position '',3f12.3)' ) trim ( word ), cor3(1:3,j)
        call s_blanks_delete ( text )
        write ( iunit, '(6x,a)' ) trim ( text )
        text_num = text_num + 1
      end do

      write ( iunit, '(a)' )     '    }'
      text_num = text_num + 1

    end if
!
!  Faces.
!
    write ( iunit, '(a)' ) ' '
    write ( text, '(a,i8)' ) 'polygons ', face_num
    call s_blanks_delete ( text )
    write ( iunit, '(4x,a)' ) trim ( text )
    write ( iunit, '(a)' ) '    {'
    text_num = text_num + 3

    do iface = 1, face_num

      write ( word, '( ''['', i8, '']'' )' ) iface-OFFSET
      call s_blank_delete ( word )
      write ( text, '(a,'' nodes '',i8 )' ) trim ( word ), face_order(iface)
      call s_blanks_delete ( text )
      write ( iunit, '(6x,a)' ) trim ( text )
      write ( iunit, '(a)' ) '      {'
      text_num = text_num + 2

      do ivert = 1, face_order(iface)

        write ( word, '( ''['', i8, '']'' )' ) ivert-OFFSET
        call s_blank_delete ( word )
        write ( text, '( a,'' vertex '',i8 )' ) trim ( word ), &
          face(ivert,iface) - OFFSET
        call s_blanks_delete ( text )
        write ( iunit, '(8x,a)' ) trim ( text )

        x = vertex_normal(1,ivert,iface)
        y = vertex_normal(2,ivert,iface)
        z = vertex_normal(3,ivert,iface)
        write ( text, '(a,3f12.3)' ) 'normal ', x, y, z
        call s_blanks_delete ( text )
        write ( iunit, '(12x,a)' ) trim ( text )

        x = vertex_tex_uv(1,ivert,iface)
        y = vertex_tex_uv(2,ivert,iface)
        write ( text, '(a,2f12.3)' ) 'uvTexture ', x, y
        call s_blanks_delete ( text )
        write ( iunit, '(12x,a)' ) trim ( text )

        text = 'vertexColor 255 178 178 178'
        call s_blanks_delete ( text )
        write ( iunit, '(12x,a)' ) trim ( text )

        text_num = text_num + 4

      end do

      write ( iunit, '(a)' ) '      }'
      write ( text, '(''material '',i8)' ) face_material(iface) - OFFSET
      call s_blanks_delete ( text )
      write ( iunit, '(6x,a)' ) trim ( text )
      text_num = text_num + 2

    end do
 
    write ( iunit, '(a)' ) '    }'
    write ( iunit, '(a)' ) '  }'
 
    text_num = text_num + 2

  end if
!
!  IndexedLineSet.
!
  if ( line_num > 0 ) then

    nseg = 0

    jhi = 0

10      continue

    jlo = jhi
!
!  Look for the next index JLO that is not -1.
!
11      continue

    jlo = jlo + 1

    if ( jlo > line_num ) then
      go to 20
    end if

    if ( line_dex(jlo) == -1+OFFSET ) then
      go to 11
    end if
!
!  Look for the highest following index JHI that is not -1.
!
    jhi = jlo + 1

    if ( jhi > line_num ) then
      go to 20
    end if

    if ( line_dex(jhi) == -1+OFFSET ) then
      go to 10
    end if

15      continue
    
    if ( jhi < line_num ) then
      if ( line_dex(jhi+1) /= -1+OFFSET ) then
        jhi = jhi + 1
        go to 15
      end if
    end if
!
!  Our next line segment involves LINE_DEX indices JLO through JHI.
!
    nseg = nseg + 1
    write ( text, '(''spl'', i8 )' ) nseg
    call s_blank_delete ( text )

    npts = jhi + 1 - jlo

    write ( iunit, '(a)' ) ' '
    write ( iunit, '(a)' ) '  spline'
    write ( iunit, '(a)' ) '  {'
    write ( iunit, '(a)' ) '    name     "' // trim ( text ) // '"'
    write ( iunit, '(a)' ) '    type     LINEAR'
    write ( iunit, '(a,i8)' ) '    nbKeys   ', npts
    write ( iunit, '(a)' ) '    tension  0.000'
    write ( iunit, '(a)' ) '    step     1'
    write ( iunit, '(a)' ) ' '
    text_num = text_num + 9

    write ( iunit, '(a)' ) '    controlpoints'
    write ( iunit, '(a)' ) '    {'
    text_num = text_num + 2

    do j = jlo, jhi
 
      jrel = j - jlo
      k = line_dex(j)
      write ( word, '( ''['', i8, '']'')' ) jrel
      call s_blank_delete ( word )
      write ( text, '( a, '' position '', 3f12.4)' ) &
        trim ( word ), cor3(1,k), cor3(2,k), cor3(3,k)
      call s_blanks_delete ( text )
      write ( iunit, '(6x,a)' ) trim ( text )
      text_num = text_num + 1

    end do

    write ( iunit, '(a)' ) '    }'
    write ( iunit, '(a)' ) '  }'
    text_num = text_num + 2

    go to 10

20      continue

  end if
!
!  MATERIALS
!
  do i = 1, material_num

    write ( text, '(''['', i8, '']'' )' ) i-OFFSET
    call s_blank_delete ( text )
    write ( iunit, '(a)' ) '  material ' // trim ( text )

    write ( iunit, '(a)' ) '  {'
    write ( iunit, '(a)' ) 'name          "' // trim ( material_name(i) ) // '"'
    write ( iunit, '(a)' ) '    type           PHONG'
    write ( iunit, '(a,3f10.4)' ) '    ambient        ', &
      material_rgba(1,i), material_rgba(2,i), material_rgba(3,i)
    write ( iunit, '(a,3f10.4)' ) '    diffuse        ', &
      material_rgba(1,i), material_rgba(2,i), material_rgba(3,i)
    write ( iunit, '(a,3f10.4)' ) '    specular       ', &
      material_rgba(1,i), material_rgba(2,i), material_rgba(3,i)
    write ( iunit, '(a)' ) '    exponent      50.0'
    write ( iunit, '(a)' ) '    reflectivity   0.0'
    write ( iunit, '(a,f10.4)' ) '    transparency   ', 1.0 - material_rgba(4,i)
    write ( iunit, '(a)' ) '    refracIndex    1.0'
    write ( iunit, '(a)' ) '    glow           0'
    write ( iunit, '(a)' ) '    coc            0.0'
    write ( iunit, '(a)' ) '  }'

    text_num = text_num + 14

  end do
!
!  TEXTURES
!
  do i = 1, texture_num

    write ( text, '(''['', i8, '']'' )' ) i-OFFSET
    call s_blank_delete ( text )
    write ( iunit, '(a)' ) '  texture [' // trim ( text ) //']'

    write ( iunit, '(a)' ) '{'
    write ( iunit, '(a)' ) 'name          "' // trim ( texture_name(i) ) // '"'
    write ( iunit, '(a)' ) 'glbname       "t2d1"'
    write ( iunit, '(a)' ) 'anim          STATIC'
    write ( iunit, '(a)' ) 'method        XY'
    write ( iunit, '(a)' ) 'repeat        1  1'
    write ( iunit, '(a)' ) 'scaling       1.000  1.000'
    write ( iunit, '(a)' ) 'offset        0.000  0.000'
    write ( iunit, '(a)' ) 'pixelInterp'
    write ( iunit, '(a)' ) 'effect        INTENSITY'
    write ( iunit, '(a)' ) 'blending      1.000'
    write ( iunit, '(a)' ) 'ambient       0.977'
    write ( iunit, '(a)' ) 'diffuse       1.000'
    write ( iunit, '(a)' ) 'specular      0.966'
    write ( iunit, '(a)' ) 'reflect       0.000'
    write ( iunit, '(a)' ) 'transp        0.000'
    write ( iunit, '(a)' ) 'roughness     0.000'
    write ( iunit, '(a)' ) 'reflMap       1.000'
    write ( iunit, '(a)' ) 'rotation      0.000'
    write ( iunit, '(a)' ) 'txtsup_rot    0.000  0.000  0.000'
    write ( iunit, '(a)' ) 'txtsup_trans  0.000  0.000  0.000'
    write ( iunit, '(a)' ) 'txtsup_scal   1.000  1.000  1.000'
    write ( iunit, '(a)' ) '}'
    write ( iunit, '(1x)' ) 

    text_num = text_num + 25

  end do

  write ( iunit, '(a)' ) '}'
  text_num = text_num + 1
!
!  Report.
!
  write ( *, * ) ' '
  write ( *, * ) 'HRC_WRITE - Wrote ', text_num, ' text lines.'

  return
end
function i_modp ( i, j )
!
!*******************************************************************************
!
!! I_MODP returns the nonnegative remainder of integer division.
!
!
!  Formula:
!
!    If
!      NREM = I_MODP ( I, J )
!      NMULT = ( I - NREM ) / J
!    then
!      I = J * NMULT + NREM
!    where NREM is always nonnegative.
!
!  Comments:
!
!    The MOD function computes a result with the same sign as the
!    quantity being divided.  Thus, suppose you had an angle A,
!    and you wanted to ensure that it was between 0 and 360.
!    Then mod(A,360) would do, if A was positive, but if A
!    was negative, your result would be between -360 and 0.
!
!    On the other hand, I_MODP(A,360) is between 0 and 360, always.
!
!  Example:
!
!        I         J     MOD  I_MODP   I_MODP Factorization
!
!      107        50       7       7    107 =  2 *  50 + 7
!      107       -50       7       7    107 = -2 * -50 + 7
!     -107        50      -7      43   -107 = -3 *  50 + 43
!     -107       -50      -7      43   -107 =  3 * -50 + 43
!
!  Modified:
!
!    02 March 199
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
!    Output, integer I_MODP, the nonnegative remainder when I is
!    divided by J.
!
  integer i
  integer i_modp
  integer j
!
  if ( j == 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'I_MODP - Fatal error!'
    write ( *, * ) '  I_MODP ( I, J ) called with J = ', j
    stop
  end if

  i_modp = mod ( i, j )

  if ( i_modp < 0 ) then
    i_modp = i_modp + abs ( j )
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
subroutine interact ( byte_swap, cor3, cor3_material, cor3_new, cor3_normal, &
  cor3_num, cor3_tex_uv, debug, edge, face, face_area, face_material, &
  face_normal, face_object, face_order, face_rank, face_tex_uv, face_tier, &
  filein_name, fileout_name, ierror, line_dex, line_material, line_prune, list, &
  material_name, material_rgba, cor3_max, edge_max, face_max, line_max, &
  material_max, order_max, texture_max, normal_temp, object_name, texture_name, &
  texture_temp, transform_matrix, vertex_material, vertex_normal, vertex_tex_uv )
!
!*******************************************************************************
!
!! INTERACT interacts with the user to specify input and output files.
!
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
!    Workspace, real COR3(3,COR3_MAX), the coordinates of points.
!
!    Workspace, integer COR3_MATERIAL(COR3_MAX), the material of each node.
!
!    Workspace, real COR3_TEX_UV(2,COR3_MAX), UV texture coordinates for nodes.
!
!    Workspace, integer FACE(ORDER_MAX,FACE_MAX), the nodes making faces.
!
!    Workspace, integer FACE_MATERIAL(FACE_MAX), the material of each face.
!
!    Workspace, real FACE_NORMAL(3,FACE_MAX), the normal vector at each face.
!
!    Workspace, integer FACE_ORDER(FACE_MAX), the number of vertices per face.
!
!    Output, integer IERROR, error flag.  0 = no error, nonzero = error.
!
!    Workspace, integer LINE_DEX(LINE_MAX), nodes forming a line, terminated by -1.
!
!    Workspace, integer LINE_MAT(LINE_MAX), material index for line.
!
!    Workspace, character ( len = 100 ) MATERIAL_NAME(MATERIAL_MAX), material names.
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
!    Workspace, character ( len = 100 ) TEXTURE_NAME(TEXTURE_MAX), texture names.
!
!    Workspace, integer VERTEX_MAT(ORDER_MAX,FACE_MAX), vertex materials.
!
!    Workspace, real VERTEX_NORMAL(3,ORDER_MAX,FACE_MAX), normals at vertices.  
!
!    Workspace, real VERTEX_TEX_UV(2,ORDER_MAX,FACE_MAX), vertex texture coordinates.
!
  integer cor3_max
  integer edge_max
  integer face_max
  integer line_max
  integer material_max
  integer order_max
  integer texture_max
!
  logical byte_swap
  character ( len = 100 ) command
  real cor3(3,cor3_max)
  integer cor3_material(cor3_max)
  real cor3_new(3,cor3_max)
  real cor3_normal(3,cor3_max)
  integer cor3_num
  integer cor3_number(cor3_max)
  real cor3_tex_uv(2,cor3_max)
  logical debug
  integer edge(4,edge_max)
  integer face(order_max,face_max)
  real face_area(face_max)
  integer face_material(face_max)
  real face_normal(3,face_max)
  integer face_object(face_max)
  integer face_order(face_max)
  integer face_rank(face_max)
  real face_tex_uv(2,face_max)
  integer face_tier(face_max)
  character ( len = 100 ) filein_name
  character ( len = 10 ) filein_type
  character ( len = 100 ) fileout_name
  character ( len = 10 ) fileout_type
  integer i
  integer icor3
  integer ierror
  integer iface
  integer ios
  integer iseed
  integer ivert
  integer j
  integer line_dex(line_max)
  integer line_material(line_max)
  integer line_prune
  integer list(cor3_max)
  character ( len = 100 ) material_name(material_max)
  real material_rgba(4,material_max)
  real normal_temp(3,order_max*face_max)
  integer color_num
  integer edge_num
  integer face_num
  integer group_num
  integer line_num
  integer material_num
  integer object_num
  integer texture_num
  character ( len = 100 ) object_name
  logical s_eqi
  character ( len = 100 ) texture_name(texture_max)
  real texture_temp(2,order_max*face_max)
  real transform_matrix(4,4)
  integer vertex_material(order_max,face_max)
  real vertex_normal(3,order_max,face_max)
  real vertex_tex_uv(2,order_max,face_max)
  real x
  real y
  real z
!
!  Print an introductory message.
!
  call hello ( cor3_max, face_max, line_max, material_max, order_max, texture_max )
!
!  Get the next user command.
!
10    continue
 
  ierror = 0

  write ( *, * ) ' '
  write ( *, * )'Enter command ( or "Help" )'

  read ( *, '(a)', iostat = ios ) command

  if ( ios /= 0 ) then
    ierror = ios
    write ( *, * ) ' '
    write ( *, * ) 'INTERACT - Fatal error!'
    write ( *, * ) '  Unexpected end of input.'
    return
  end if
!
!  << means read a new file, and add it to the current information.
!
  if ( command(1:2) == '<<' ) then
 
    filein_name = command(3:)
    call s_blank_delete ( filein_name )
 
    call infile ( filein_name, ierror, filein_type )

    if ( ierror /= 0 ) then
      write ( *, * ) ' '
      write ( *, * ) 'INTERACT - Warning!'
      write ( *, * ) '  The input file name was unacceptable!'
      go to 10
    end if

    call data_read ( cor3, cor3_material, cor3_normal, cor3_tex_uv, &
      debug, face, face_area, face_material, face_normal, face_order, face_tex_uv, &
      filein_name, filein_type, ierror, iseed, line_dex, line_material, &
      material_name, material_rgba, cor3_max, face_max, line_max, material_max, &
      order_max, texture_max, normal_temp, color_num, cor3_num, face_num, &
      group_num, line_num, material_num, object_num, texture_num, texture_name, &
      texture_temp, vertex_material, vertex_normal, vertex_tex_uv )

    if ( ierror /= 0 ) then
      write ( *, * ) ' ' 
      write ( *, * ) 'INTERACT - Warning!'
      write ( *, * ) '  The input data was not read properly.'
      ierror = 0
    end if
!
!  Get the next input file to be examined.
!
  else if ( command(1:1) == '<' ) then
 
    filein_name = command(2:)
    call s_blank_delete ( filein_name )
 
    call infile ( filein_name, ierror, filein_type )

    if ( ierror /= 0 ) then
      write ( *, * ) ' '
      write ( *, * ) 'INTERACT - Warning!'
      write ( *, * ) '  The input file name was unacceptable!'
      go to 10
    end if

    call data_init ( cor3, cor3_material, cor3_normal, cor3_tex_uv, face, &
      face_area, face_material, face_normal, face_order, face_tex_uv, iseed, &
      line_dex, line_material, material_name, material_rgba, cor3_max, face_max, &
      line_max, material_max, order_max, texture_max, normal_temp, color_num, &
      cor3_num, face_num, group_num, line_num, material_num, object_num, texture_num, &
      object_name, texture_name, texture_temp, transform_matrix, vertex_material, &
      vertex_normal, vertex_tex_uv )

    call data_read ( cor3, cor3_material, cor3_normal, cor3_tex_uv, &
      debug, face, face_area, face_material, face_normal, face_order, face_tex_uv, &
      filein_name, filein_type, ierror, iseed, line_dex, line_material, &
      material_name, material_rgba, cor3_max, face_max, line_max, material_max, &
      order_max, texture_max, normal_temp, color_num, cor3_num, face_num, &
      group_num, line_num, material_num, object_num, texture_num, texture_name, &
      texture_temp, vertex_material, vertex_normal, vertex_tex_uv )

    if ( ierror /= 0 ) then
      write ( *, * ) ' ' 
      write ( *, * ) 'INTERACT - Warning!'
      write ( *, * ) '  The input data was not read properly.'
      ierror = 0
    end if
!
!  The ">" command specifies output.
!
  else if ( command(1:1) == '>' ) then

    fileout_name = command(2:)

    call s_blank_delete ( fileout_name )
!
!  Check the output filename.
!
    call outfile ( filein_name, fileout_name, ierror, fileout_type )

    if ( ierror /= 0 ) then
      write ( *, * ) ' '
      write ( *, * ) 'INTERACT - Warning!'
      write ( *, * ) '  Improper output file name.'
      ierror = 0
      go to 10
    end if
!
!  Write the output file.
!
    call data_write ( cor3, cor3_material, cor3_normal, cor3_tex_uv, &
      debug, face, face_material, face_normal, face_order, face_tex_uv, &
      filein_name, fileout_name, fileout_type, ierror, line_dex, &
      line_material, line_prune, material_name, material_rgba, cor3_max, &
      face_max, line_max, material_max, order_max, texture_max, cor3_num, &
      face_num, line_num, material_num, texture_num, object_name, &
      texture_name,vertex_material, vertex_normal, vertex_tex_uv )

    if ( ierror /= 0 ) then
      write ( *, * ) ' '
      write ( *, * ) 'INTERACT - Warning!'
      write ( *, * ) '  An error occurred during output.'
      ierror = 0
      go to 10
    end if
!
!  B: Switch byte swapping option:
!
  else if ( s_eqi ( command(1:1), 'B' ) ) then

    byte_swap = .not. byte_swap

    if ( byte_swap ) then
      write ( *, * ) 'Byte swapping set to TRUE.'
    else
      write ( *, * ) 'Byte swapping set to FALSE.'
    end if
!
!  D: Switch debug option:
!
  else if ( s_eqi ( command(1:1), 'D' ) ) then

    debug = .not. debug

    if ( debug ) then
      write ( *, * ) 'Debug option set to TRUE.'
    else
      write ( *, * ) 'Debug option set to FALSE.'
    end if
!
!  F: Check a face:
!
  else if ( s_eqi ( command(1:1), 'F' ) ) then

    if ( face_num <= 0 ) then
      write ( *, * ) ' '
      write ( *, * ) 'Input a graphical object with faces first!'
      go to 10
    end if

    write ( *, * ) ' '
    write ( *, * ) 'Enter a face between 1 and ', face_num
    read ( *, * ) iface

    call face_print ( cor3, face, face_material, face_normal, face_order, iface, &
      cor3_max, face_max, order_max, face_num, vertex_material, vertex_normal )
!
!  HELP:
!
  else if ( s_eqi ( command(1:1), 'H' ) ) then
 
    call help
!
!  INFO:
!
  else if ( s_eqi ( command, 'INFO' ) ) then

    call news
!
!  INVERT:
!  Make an inverted copy of the object to give it thickness.
!
  else if ( s_eqi ( command(1:3), 'INV' ) ) then

    call object_invert ( cor3, cor3_material, cor3_normal, face, &
      face_material, face_normal, face_order, material_name, &
      material_rgba, cor3_max, face_max, material_max, order_max, &
      cor3_num, face_num, material_num, vertex_material, vertex_normal )
!
!  LINE_PRUNE:
!  Set line pruning option.
!
  else if ( s_eqi ( command, 'LINE_PRUNE' ) ) then

    write ( *, * ) ' '
    write ( *, * ) 'SET THE LINE PRUNING OPTION:'
    write ( *, * ) ' '
    write ( *, * ) '  0 means no pruning;'
    write ( *, * ) '  nonzero means only generate line (I,J)'
    write ( *, * ) '    if I < J.'
    write ( *, * ) ' '
    write ( *, * ) '  The line pruning option is now ', line_prune
    write ( *, * ) ' '
    write ( *, * ) '  Enter your line_pruning option:'

    read ( *, * ) line_prune
!
!  LINES:
!  Convert all faces to lines.
!
  else if ( s_eqi ( command, 'LINES' ) ) then
 
    if ( face_num > 0 ) then

      write ( *, * ) ' '
      write ( *, * ) 'INTERACT - Note:'
      write ( *, * ) '  Face information will be converted'
      write ( *, * ) '  to line information.'

      call face_to_line ( debug, face, face_order, line_dex, line_material, &
        line_prune, face_max, line_max, order_max, face_num, line_num, &
        vertex_material )

      if ( line_num > line_max ) then

        write ( *, * ) ' '
        write ( *, * ) 'INTERACT - Note:'
        write ( *, * ) '  Some face information was lost.'
        write ( *, * ) '  The maximum number of lines is ', line_max
        write ( *, * ) '  but we would need at least ', line_num

        line_num = line_max

      end if

      face_num = 0

    else

      write ( *, * ) ' '
      write ( *, * ) 'INTERACT - Note:'
      write ( *, * ) '  There were no faces to convert.'

    end if
!
!  NORMALS: recompute the normal vectors for vertices on faces,
!  and average these to get face normal vectors.
!
  else if ( s_eqi ( command(1:1), 'N' ) ) then

    do iface = 1, face_num
      do ivert = 1, face_order(iface)
        vertex_normal(1:3,ivert,iface) = 0.0
      end do
    end do

    call vertex_normal_set ( cor3, face, face_order, cor3_max, &
      face_max, order_max, face_num, vertex_normal )

    cor3_normal(1:3,1:cor3_num) = 0.0

    call cor3_normal_set ( cor3_normal, face, face_area, &
      face_order, cor3_max, face_max, order_max, face_num, vertex_normal )

    face_normal(1:3,1:face_num) = 0.0

    call face_normal_ave ( face_normal, face_order, face_max, &
      order_max, face_num, vertex_normal )
!
!  OHELL: 
!    Use the node normals.
!    Set the vertex normals equal to the node normals.
!    Set the face normals by averaging vertex normals.
!
  else if ( s_eqi ( command(1:1), 'O' ) ) then
!
!  Recompute any zero vertex normals from vertex positions.
!
    write ( *, * ) ' '
    write ( *, * ) 'Making sure all vertex normals are defined.'

    call vertex_normal_set ( cor3, face, face_order, cor3_max, &
      face_max, order_max, face_num, vertex_normal )
!
!  Compute node normals by averaging vertex normals.
!
    write ( *, * ) ' '
    write ( *, * ) 'Averaging vertex normals at each node.'

    cor3_normal(1:3,1:cor3_num) = 0.0

    call cor3_normal_set ( cor3_normal, face, face_area, &
      face_order, cor3_max, face_max, order_max, face_num, vertex_normal )
!
!  Copy node normals into vertex normals.
!
    write ( *, * ) ' '
    write ( *, * ) 'Replacing vertex normals by average.'

    do iface = 1, face_num
      do ivert = 1, face_order(iface)
        icor3 = face(ivert,iface)
        vertex_normal(1:3,ivert,iface) = cor3_normal(1:3,icor3)
      end do
    end do
!
!  Recompute zero face normals by averaging vertex normals.
!
    write ( *, * ) ' '
    write ( *, * ) 'Averaging vertex normals to get face normals.'

    call face_normal_ave ( face_normal, face_order, face_max, &
      order_max, face_num, vertex_normal )
!
!  QUIT:
!
  else if ( s_eqi ( command(1:1), 'Q' ) ) then
 
    write ( *, * ) ' '
    write ( *, * ) 'INTERACT - End of interaction.'
    ierror = 0
    return
!
!  REVERSE: Reverse normal vectors.
!
  else if ( s_eqi ( command(1:3), 'REV' ) ) then
 
    cor3_normal(1:3,1:cor3_num) = - cor3_normal(1:3,1:cor3_num)

    face_normal(1:3,1:face_num) = - face_normal(1:3,1:face_num)

    do iface = 1, face_num
      do ivert = 1, face_order(iface)
        vertex_normal(1:3,ivert,iface) = - vertex_normal(1:3,ivert,iface)
      end do
    end do

    write ( *, * ) ' '
    write ( *, * ) 'INTERACT - Note:'
    write ( *, * ) '  Reversed node, face and vertex normals.'
!
!  RELAX: Smooth the surface via relaxation.
!
  else if ( s_eqi ( command(1:3), 'REL' ) ) then

    call node_relax ( cor3, cor3_new, cor3_number, face, face_order, cor3_max, &
      face_max, order_max, cor3_num, face_num )

    cor3(1:3,1:cor3_num) = cor3_new(1:3,1:cor3_num)
!
!  S: Select a few faces, discard rest:
!
  else if ( s_eqi ( command(1:1), 'S' ) ) then
 
    call face_subset ( cor3, face, face_material, face_normal, &
      face_order, ierror, list, cor3_max, face_max, order_max, &
      cor3_num, face_num, line_num, vertex_material, vertex_normal )
!
!  T: Transform data.
!
  else if ( s_eqi ( command(1:1), 'T' ) ) then

    write ( *, * ) ' '
    write ( *, * ) 'For now, we only offer point scaling.'
    write ( *, * ) 'Enter X, Y, Z scale factors:'
    read ( *, * ) x, y, z
    do j = 1, cor3_num
      cor3(1,j) = cor3(1,j) * x
      cor3(2,j) = cor3(2,j) * y
      cor3(3,j) = cor3(3,j) * z
    end do

    call cor3_range ( cor3_max, cor3_num, cor3 )

    do iface = 1, face_num
      do ivert = 1, face_order(iface)
        vertex_normal(1:3,ivert,iface) = 0.0
      end do
    end do

    call vertex_normal_set ( cor3, face, face_order, cor3_max, face_max, &
      order_max, face_num, vertex_normal )

    cor3_normal(1:3,1:cor3_num) = 0.0

    call cor3_normal_set ( cor3_normal, face, face_area, &
      face_order, cor3_max, face_max, order_max, face_num, vertex_normal )

    face_normal(1:3,1:face_num) = 0.0

    call face_normal_ave ( face_normal, face_order, face_max, order_max, &
      face_num, vertex_normal )
!
!  U: Renumber faces, count objects:
!
  else if ( s_eqi ( command(1:1), 'U' ) ) then
 
    call face_check ( edge, face, face_material, face_normal, face_object, &
      face_order, face_rank, face_tier, edge_max, face_max, order_max, &
      edge_num, face_num, object_num, vertex_material, vertex_normal )
!
!  V: Convert polygons to triangles.
!
  else if ( s_eqi ( command(1:1), 'V' ) ) then

    write ( *, * ) ' '
    write ( *, * ) 'Convert polygonal faces to triangles.'

    call poly_2_tri ( face, face_material, face_order, ierror, face_max, &
      order_max, face_num, vertex_material )

    if ( ierror /= 0 ) then

      write ( *, * ) ' '
      write ( *, * ) 'Conversion attempt abandoned.'

    else

      write ( *, * ) ' '
      write ( *, * ) 'Number of faces is now ', face_num

      do iface = 1, face_num
        do ivert = 1, face_order(iface)
          vertex_normal(1:3,ivert,iface) = 0.0
        end do
      end do

      call vertex_normal_set ( cor3, face, face_order, cor3_max, &
        face_max, order_max, face_num, vertex_normal )

      cor3_normal(1:3,1:cor3_num) = 0.0

      call cor3_normal_set ( cor3_normal, face, face_area, &
        face_order, cor3_max, face_max, order_max, face_num, vertex_normal )

      face_normal(1:3,1:face_num) = 0.0

      call face_normal_ave ( face_normal, face_order, face_max, &
        order_max, face_num, vertex_normal )

    end if
!
!  W: Reverse the order of the nodes that define each face.
!
  else if ( s_eqi ( command(1:1), 'W' ) ) then
 
    call face_reverse_order ( cor3_normal, face, face_normal, &
      face_order, cor3_max, face_max, order_max, cor3_num, face_num, &
      vertex_material, vertex_normal, vertex_tex_uv )
!
!  Unintelligible!
!
  else
 
    write ( *, * ) ' '
    write ( *, * ) 'INTERACT - Warning!'
    write ( *, * ) '  Your command was not recognized:'
    write ( *, '(a)' ) trim ( command )
 
  end if
 
  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'INTERACT - Warning!'
    write ( *, * ) '  An error occurred during this action.'
    ierror = 0
  end if

  go to 10
end
subroutine intnex ( line, ival, done )
!
!*******************************************************************************
!
!! INTNEX "reads" integers from a string, one at a time.
!
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
!    Input, character ( len = * ) LINE, a string, presumably containing
!    integers.  These may be separated by spaces or commas.
!
!    Output, integer IVAL.  If DONE is FALSE, then IVAL contains the
!    "next" integer read from LINE.  If DONE is TRUE, then
!    IVAL is zero.
!
!    Input/output, logical DONE.
!    On input with a fresh value of LINE, the user should set
!    DONE to TRUE.
!    On output, the routine sets DONE to FALSE if another integer
!    was read, or TRUE if no more integers could be read.
!
  logical done
  integer ierror
  integer ival
  integer lchar
  character ( len = * ) line
  integer, save :: next = 1
!
  ival = 0

  if ( done ) then
    next = 1
    done = .false.
  end if

  if ( next > len(line) ) then
    done = .true.
    return
  end if

  call s_to_i ( line(next:), ival, ierror, lchar )

  if ( ierror /= 0 .or. lchar == 0 ) then
    done = .true.
    next = 1
  else
    done = .false.
    next = next + lchar
  end if

  return
end
subroutine iv_read ( cor3, debug, face, face_order, ierror, iunit, line_dex, &
  line_material, cor3_max, face_max, line_max, order_max, texture_max, normal_temp, &
  bad_num, color_num, cor3_num, face_num, line_num, material_num, texture_num, &
  text_num, texture_name, texture_temp, vertex_material, vertex_normal, &
  vertex_tex_uv )
!
!*******************************************************************************
!
!! IV_READ reads graphics information from an Inventor file.
!
!
!  Diagnostics:
!
!    For now, we are going to apply the following kludge, which is an
!    improvement over the previous situation.  The transform matrix
!    will be initialized to the identity; every time a new transform
!    matrix is specified, it will OVERWRITE the old one.  Every point
!    and vector that is read in will be multiplied by the current
!    transform matrix.  That's it for now.  We need to start using
!    the transform matrix, and eventually, we need to start using
!    it more accurately than this.
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
!     #Inventor V2.0 ascii
!
!     Separator {
!       Info {
!         string "Inventor file generated by IVREAD.
!         Original data in file cube.iv."
!       }
!       Separator {
!         LightModel {
!           model PHONG
!         }
!         MatrixTransform { matrix
!           0.9  0.0  0.0  0.0
!           0.0 -0.9  0.0  0.0
!           0.0  0.0 -1.5  0.0
!           0.0  0.0  0.0  1.0
!         }
!         Material {
!           ambientColor  0.2 0.2 0.2
!           diffuseColor  [
!             0.8 0.8 0.8,
!             0.7 0.1 0.1,
!             0.1 0.8 0.2,
!           ]
!           emissiveColor 0.0 0.0 0.0
!           specularColor 0.0 0.0 0.0
!           shininess     0.2
!           transparency  [
!             0.0, 0.5, 1.0,
!           ]
!         }
!         Texture2 {
!           filename      "fred.rgb"
!           wrapS         REPEAT
!           wrapT         REPEAT
!           model         MODULATE
!           blendColor    0.0 0.0 0.0
!         }
!         TextureCoordinateBinding {
!           value PER_VERTEX_INDEXED
!         }
!         MaterialBinding {
!           value PER_VERTEX_INDEXED
!         }
!         NormalBinding {
!           value PER_VERTEX_INDEXED
!         }
!         ShapeHints {
!           vertexOrdering COUNTERCLOCKWISE
!           shapeType UNKNOWN_SHAPE_TYPE
!           faceType CONVEX
!           creaseAngle 6.28319
!         }
!         Coordinate3 {
!           point [
!                8.59816       5.55317      -3.05561,
!                8.59816       2.49756      0.000000E+00,
!                ...etc...
!                2.48695       2.49756      -3.05561,
!           ]
!         }
!         TextureCoordinate2 {
!           point [
!                0.0  1.0,
!                0.1, 0.8,
!                ...etc...
!                0.4  0.7,
!           ]
!         }
!         Normal {
!           vector [
!             0.71 0.71 0.0,
!             ...etc...
!             0.32 0.32 0.41,
!           ]
!         }
!
!         IndexedLineSet {
!           coordIndex [
!              0,    1,    2,   -1,    3,    4,    5,   -1,    7,    8,    9,
!            ...etc...
!            191,   -1,
!           ]
!           materialIndex [
!              0,    0,    0,   -1,    1,    1,    1,   -1,    2,    2,    2,
!            ...etc...
!             64,   -1,
!           ]
!         }
!
!         IndexedFaceSet {
!           coordIndex [
!              0,    1,    2,   -1,    3,    4,    5,   -1,    7,    8,    9,
!            ...etc...
!            191,   -1,
!           ]
!           materialIndex [
!              0,    0,    0,   -1,    1,    1,    1,   -1,    2,    2,    2,
!            ...etc...
!             64,   -1,
!           ]
!           normalIndex [
!              0,    0,    0,   -1,    1,    1,    1,   -1,    2,    2,    2,
!            ...etc...
!             64,   -1,
!           ]
!           textureCoordIndex [
!              0,    0,    0,   -1,    1,    1,    1,   -1,    2,    2,    2,
!            ...etc...
!             64,   -1,
!           ]
!         }
!
!         IndexedTriangleStripSet {
!           vertexProperty VertexProperty {
!             vertex [ x y z,
!                      ...
!                      x y z ]
!             normal [ x y z,
!                      ...
!                      x y z ]
!             materialBinding OVERALL
!             normalBinding PER_VERTEX_INDEXED
!           }
!           coordIndex [ i, j, k, l, m, -1, n, o, p, q, r, s, t, u, -1,
!             v, w, x, -1 ..., -1 ]
!           normalIndex -1
!         }
!
!       }
!     }
!
!  Modified:
!
!    28 June 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real COR3(3,COR3_MAX), the coordinates of points.
!
!    Input, logical DEBUG, debugging switch.
!
!    Input/output, integer FACE(ORDER_MAX,FACE_MAX), the nodes making faces.
!
!    Input/output, integer FACE_ORDER(FACE_MAX), the number of vertices 
!    per face.
!
!    Output, integer IERROR, an error flag.
!
!    Input, integer IUNIT, the FORTRAN unit from which data is read.
!
!    Input/output, integer LINE_DEX(LINE_MAX), nodes forming a line, 
!    terminated by -1.
!
!    Input/output, integer LINE_MAT(LINE_MAX), material index for line.
!
!    Input, integer COR3_MAX, the maximum number of points.
!
!    Input, integer FACE_MAX, the maximum number of faces.
!
!    Input, integer LINE_MAX, the maximum number of line definition items.
!
!    Input, integer ORDER_MAX, the maximum number of vertices per face.
!
!    Input, integer TEXTURE_MAX, the maximum number of textures.
!
!    Workspace, NORMAL_TEMP(3,ORDER_MAX*MAXFACE).
!
!    Output, integer BAD_NUM, number of bad lines of text read.
!
!    Output, integer COLOR_NUM, number of colors.
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
!    Output, integer TEXT_NUM, number of lines in the file.
!
!    Input/output, character ( len = 100 ) TEXTURE_NAME(TEXTURE_MAX), texture names.
!
!    Workspace, real TEXTURE_TEMP(2,ORDER_MAX*FACE_MAX), texture coordinates.
!
!    Output, real TRANSFORM_MATRIX(4,4), the transformation matrix.
!
!    Input/output, integer VERTEX_MAT(ORDER_MAX,FACE_MAX), vertex materials.
!
!    Output, real VERTEX_NORMAL(3,ORDER_MAX,FACE_MAX), normals at vertices.  
!
!    Input/output, real VERTEX_TEX_UV(2,ORDER_MAX,FACE_MAX), vertex texture coordinates.
!
  integer, parameter :: level_max = 10
  integer, parameter :: OFFSET = 1
!
  integer cor3_max
  integer face_max
  integer line_max
  integer order_max
  integer texture_max
!
  real b
  character ( len = 4 ) char4
  real cor3(3,cor3_max)
  logical debug
  logical done
  integer face(order_max,face_max)
  integer face_order(face_max)
  real g
  integer i
  integer icol
  integer icolor
  integer ierror
  integer iface
  integer iface_num
  integer iprint
  integer irow
  integer iunit
  integer iuv
  integer ivert
  integer iword
  integer ix
  integer ixyz
  integer iy
  integer iz
  integer jval
  integer k
  integer level
  character ( len = 256 ) level_name(0:level_max)
  integer line_dex(line_max)
  integer line_material(line_max)
  logical lval
!     character ( len = 30 ) material_binding
  integer nlbrack
!     character ( len = 30 ) normal_binding
  real normal_temp(3,order_max*face_max)
  integer nrbrack
  integer nu
  integer bad_num
  integer color_num
  integer cor3_num
  integer cor3_num_old
  integer face_num
  integer face_num2
  integer line_num
  integer line_num2
  integer material_num
  integer normal_bad_num
  integer normal_temp_num
  integer texture_num
  integer text_num
  integer text_numure_temp
  integer nv
  real r
  real r3vec(3)
  real rval
  logical s_eqi
  logical s_is_i
  character ( len = 256 ) text
! character ( len = 30 ) texture_coordinate_binding
  character ( len = 100 ) texture_name(texture_max)
  real texture_temp(2,order_max*face_max)
  real transform_matrix(4,4)
  integer vertex_material(order_max,face_max)
  real vertex_normal(3,order_max,face_max)
  real vertex_tex_uv(2,order_max,face_max)
  character ( len = 256 ) word
  character ( len = 256 ) wordm1
!
  icol = 0
  ierror = 0
  iface_num = face_num
  iprint = 1
  irow = 0
  iuv = 0
  ix = 0
  ixyz = 0
  iy = 0
  iz = 0
  jval = 0
  level = 0
  level_name(0) = 'Top'
!     material_binding = 'PER_VERTEX_INDEXED'
  cor3_num_old = cor3_num
  nlbrack = 0
!     normal_binding = 'PER_VERTEX_INDEXED'
  nrbrack = 0
  nu = 0
  face_num2 = face_num
  line_num2 = line_num
  normal_bad_num = 0
  normal_temp_num = 0
  text_numure_temp = 0
  nv = 0
  rval = 0.0
  call tmat_init ( transform_matrix )
  word = ' '
  wordm1 = ' '
!
!  Read a line of text from the file.
!
10    continue
 
  read ( iunit, '(a)', end = 50 ) text
  text_num = text_num + 1

  if ( debug ) then

    if ( text_num == iprint ) then
      write ( *, * ) 'Line ', iprint
      iprint = 2 * iprint
    end if

  end if
 
  if ( text == ' ' ) then
    go to 10
  end if
!
!  The first line of the file must be the header.
!
  if ( text_num == 1 ) then

    if ( .not. s_eqi ( text(1:9), '#Inventor' ) ) then
      ierror = 1
      write ( *, * ) ' '
      write ( *, * ) 'IV_READ - Fatal error!'
      write ( *, * ) '  The input file has a bad header.'
      write ( *, '(a)' ) trim ( text )
      return
    else
      go to 10
    end if

  end if
 
  done = .true.
  iword = 0
!
!  Save the previous word read.  It helps when a word depends
!  on its context.
!
20    continue

  if ( word /= ' ' .and. word .ne. ',' ) then
    wordm1 = word
  end if
!
!  Read a word from the line.
!
  call word_nexrd ( text, word, done )
!
!  If no more words in this line, read in a whole new line.
!
  if ( done ) then
    go to 10
  end if
!
!  Skip over comments.
!
  if ( word(1:1) == '#' ) then
    go to 10
  end if
!
!  Ignore blanks and commas.
!
  if ( word == ' ' .or. word == ',' ) then
    go to 20
  end if
!
!  Count the words in the current line, and the total.
!
  iword = iword + 1
!
!  If the word is a curly or square bracket, count it.
!
  if ( word == '{' .or. word == '[' ) then

    nlbrack = nlbrack + 1

    if ( debug ) then
      write ( *, '(a)' ) word(1:1)
    end if

  else if ( word .eq. '}' .or. word == ']' ) then

    if ( debug ) then
      write ( *, '(a)' ) word(1:1)
    end if

    nrbrack = nrbrack + 1

    if ( nlbrack < nrbrack ) then
      write ( *, * ) ' '
      write ( *, * ) 'IV_READ - Fatal error!'
      write ( *, * ) '  Extraneous right bracket, line ', text_num
      write ( *, '(a)' ) trim ( text )
      write ( *, * ) 'Currently processing field:'
      write ( *, '(a)' ) trim ( level_name(level) )
      ierror = 1
      return
    end if

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
!  BASECOLOR
!
  if ( s_eqi ( level_name(level), 'BASECOLOR' ) ) then
 
    if ( word == '{' ) then

    else if ( word == '}' ) then

      level = nlbrack - nrbrack

    else if ( s_eqi ( word, 'RGB' ) ) then

    else
      go to 99
    end if
!
!  COORDINATE3
!
  else if ( s_eqi ( level_name(level), 'COORDINATE3' ) ) then
 
    if ( word == '{' ) then

    else if ( word == '}' ) then

      level = nlbrack - nrbrack

    else if ( s_eqi ( word, 'POINT' ) ) then

    else 
      go to 99
    end if
!
!  COORDINATE4
!
  else if ( s_eqi ( level_name(level), 'COORDINATE4' ) ) then
 
    if ( word == '{' ) then
 
    else if ( word == '}' ) then

      level = nlbrack - nrbrack

    else if ( s_eqi ( word, 'POINT' ) ) then

    else
      go to 99
    end if
!
!  COORDINDEX
!
  else if ( s_eqi ( level_name(level), 'COORDINDEX' ) ) then 
 
    if ( word == '[' ) then

      ivert = 0

    else if ( word == ']' ) then

      level = nlbrack - nrbrack
!
!  (indexedlineset) COORDINDEX
!
    else if ( s_eqi ( level_name(level-1), 'INDEXEDLINESET' ) ) then

      lval = s_is_i ( word, jval )
 
      if ( lval ) then

        if ( jval < -1 ) then

          bad_num = bad_num + 1

        else

          line_num = line_num + 1

          if ( line_num <= line_max ) then

            if ( jval == -1 ) then
              line_material(line_num) = -1 + OFFSET
            else
              jval = jval + cor3_num_old
            end if

            line_dex(line_num) = jval + OFFSET

          end if

        end if

      else
        go to 99
      end if
!
!  (indexedfaceset) COORDINDEX
!
    else if ( s_eqi ( level_name(level-1), 'INDEXEDFACESET' ) ) then
 
      lval = s_is_i ( word, jval )
 
      if ( lval ) then

        if ( jval == -1 ) then

          ivert = 0

        else

          ivert = ivert + 1

          if ( ivert == 1 ) then
            face_num = face_num + 1
            if ( face_num <= face_max ) then
              face_order(face_num) = 0
            end if
          end if

          if ( face_num <= face_max ) then
            face_order(face_num) = face_order(face_num) + 1
            if ( ivert <= order_max ) then
              face(ivert,face_num) = jval + cor3_num_old + OFFSET
            end if
          end if

        end if

      end if
!
!  (indexednurbssurface) COORDINDEX
!
    else if ( s_eqi ( level_name(level-1), 'INDEXEDNURBSSURFACE' ) ) then

      lval = s_is_i ( word, jval )

      if ( lval ) then

      else if ( word == '[' ) then

      else if ( word == ']' ) then
 
        do i = 1, nu-1
          line_num = line_num + 1
          if ( line_num <= line_max ) then
            line_dex(line_num) = i - 1 + cor3_num_old + OFFSET
          end if
        end do
 
        do i = nu, nu*(nv-1), nu
          line_num = line_num + 1
          if ( line_num <= line_max ) then
            line_dex(line_num) = i - 1 + cor3_num_old + OFFSET
          end if
        end do
 
        do i = nu*nv, nu*nv-nu+2, -1
          line_num = line_num + 1
          if ( line_num <= line_max ) then
            line_dex(line_num) = i - 1 + cor3_num_old + OFFSET
          end if
        end do
 
        do i = nu*(nv-1)+1, 1, -nu
          line_num = line_num + 1
          if ( line_num <= line_max ) then
            line_dex(line_num) = i - 1 + cor3_num_old + OFFSET
          end if
        end do
 
        line_num = line_num + 1
        if ( line_num <= line_max ) then
          line_dex(line_num) = -1 + OFFSET
        end if

      end if
!
!  (indexedtrianglestripset) COORDINDEX
!
!  First three coordinate indices I1, I2, I3 define a triangle.
!  Next triangle is defined by I2, I3, I4 (actually, I4, I3, I2
!  to stay with same counterclockwise sense).
!  Next triangle is defined by I3, I4, I5 ( don't need to reverse
!  odd numbered triangles) and so on.
!  List is terminated with -1.
!
    else if ( s_eqi ( level_name(level-1), 'INDEXEDTRIANGLESTRIPSET' ) ) then

      lval = s_is_i ( word, jval )
 
      if ( lval ) then

        if ( jval == -1 ) then

          ivert = 0

        else

          ivert = ivert + 1

          ix = iy
          iy = iz 
          iz = jval + cor3_num_old

          if ( ivert == 1 ) then

            face_num = face_num + 1
            if ( face_num <= face_max ) then
              face(ivert,face_num) = jval + cor3_num_old + OFFSET
              face_order(face_num) = 3
            end if

          else if ( ivert == 2 .or. ivert == 3 ) then

            if ( face_num <= face_max ) then
              face(ivert,face_num) = jval + cor3_num_old + OFFSET
            end if

          else

            face_num = face_num + 1
            if ( face_num <= face_max ) then
              face_order(face_num) = 3
              if ( mod ( ivert, 2 ) == 1 ) then
                face(1,face_num) = ix + OFFSET
                face(2,face_num) = iy + OFFSET
                face(3,face_num) = iz + OFFSET
              else
                face(1,face_num) = iz + OFFSET
                face(2,face_num) = iy + OFFSET
                face(3,face_num) = ix + OFFSET
              end if
            end if

          end if
!
!  ??? This can't be right.???
!
!  Very very tentative guess as to how indices into the normal
!  vector array are set up...
!
!            if ( face_num <= face_max .and. ivert >= 3 ) then
!               do j = 1, order_max
!                 vertex_normal(1:3,j,face_num) = normal(1:3,ix + OFFSET)
!               end do
!             end if

        end if

      end if

    end if
!
!  INDEXEDFACESET
!
  else if ( s_eqi ( level_name(level), 'INDEXEDFACESET' ) ) then
 
    if ( word == '{' ) then
 
    else if ( word == '}' ) then

      level = nlbrack - nrbrack

    else if ( s_eqi ( word, 'COORDINDEX' ) ) then

      ivert = 0
 
    else if ( s_eqi ( word, 'MATERIALINDEX' ) ) then

    else if ( s_eqi ( word, 'NORMALINDEX' ) ) then

    else if ( s_eqi ( word, 'TEXTURECOORDINDEX' ) ) then

      if ( texture_num <= 0 ) then
        texture_num = 1
        call i_to_s_zero ( texture_num, char4 )
        texture_name(texture_num) = 'Texture_' // char4
      end if

    else
      go to 99
    end if
!
!  INDEXEDLINESET
!
  else if ( s_eqi ( level_name(level), 'INDEXEDLINESET' ) ) then 
 
    if ( word == '{' ) then

    else if ( word == '}' ) then

      level = nlbrack - nrbrack
 
    else if ( s_eqi ( word, 'COORDINDEX' ) ) then

    else if ( s_eqi ( word, 'MATERIALINDEX' ) ) then

    else
      go to 99
    end if
!
!  INDEXEDNURBSSURFACE
!
  else if ( s_eqi ( level_name(level), 'INDEXEDNURBSSURFACE' ) ) then
 
    if ( word == '{' ) then

    else if ( word == '}' ) then

      level = nlbrack - nrbrack
 
    else if ( s_eqi ( word, 'NUMUCONTROLPOINTS') ) then

      call word_nexrd ( text, word, done )
      lval = s_is_i ( word, jval )

      if ( lval ) then
        nu = jval
      else
        nu = 0
        go to 99
      end if

    else if ( s_eqi ( word, 'NUMVCONTROLPOINTS' ) ) then

      call word_nexrd ( text, word, done)
      lval = s_is_i ( word, jval ) 

      if ( lval ) then
        nv = jval
      else
        nv = 0
        go to 99
      end if

    else if ( s_eqi ( word, 'COORDINDEX' ) ) then

    else if ( s_eqi ( word, 'UKNOTVECTOR' ) ) then

    else if ( s_eqi ( word, 'VKNOTVECTOR' ) ) then

    else
      go to 99
    end if
!
!  INDEXEDTRIANGLESTRIPSET
!
  else if ( s_eqi ( level_name(level), 'INDEXEDTRIANGLESTRIPSET' ) ) then 
 
    if ( word == '{' ) then

    else if ( word == '}' ) then

      level = nlbrack - nrbrack
 
    else if ( s_eqi ( word, 'VERTEXPROPERTY' ) ) then

      call word_nexrd ( text, word, done )

    else if ( s_eqi ( word, 'COORDINDEX' ) ) then

      ivert = 0

    else if ( s_eqi ( word, 'NORMALINDEX' ) ) then

      call word_nexrd ( text, word, done )

      if ( word == '[' ) then

        nlbrack = nlbrack + 1
        level = level + 1
        level_name(level) = 'NORMALINDEX'

      else if ( word == '-1' ) then

        do iface = 1, face_num
          do ivert = 1, face_order(iface)
            k = face(ivert,iface)
            vertex_normal(1:3,ivert,iface) = normal_temp(1:3,k)
          end do
        end do

      end if

    else
      go to 99
    end if
!
!  INFO
!
  else if ( s_eqi ( level_name(level), 'INFO' ) ) then

    if ( word == '{' ) then
 
    else if ( word == '}' ) then

      level = nlbrack - nrbrack

    else if ( s_eqi ( word, 'STRING' ) ) then

    else if ( word == '"' ) then

    else

    end if
!
!  LIGHTMODEL
!
!  Read, but ignore.
!
  else if ( s_eqi ( level_name(level), 'LIGHTMODEL' ) ) then

    if ( word == '{' ) then

    else if ( word == '}' ) then

      level = nlbrack - nrbrack

    else if ( s_eqi ( word, 'MODEL' ) ) then

    else

    end if
!
!  MATERIAL
!
!  Read, but ignore for now.
!
  else if ( s_eqi ( level_name(level), 'MATERIAL' ) ) then

    if ( word == '{' ) then

    else if ( word == '}' ) then

      level = nlbrack - nrbrack

    else if ( s_eqi ( word, 'AMBIENTCOLOR' ) ) then

      call word_nexrd ( text, word, done )
      call s_is_r ( word, r, lval )
      call word_nexrd ( text, word, done )
      call s_is_r ( word, g, lval )
      call word_nexrd ( text, word, done )
      call s_is_r ( word, b, lval )

    else if ( s_eqi ( word, 'EMISSIVECOLOR' ) ) then

      call word_nexrd ( text, word, done )
      call s_is_r ( word, r, lval )
      call word_nexrd ( text, word, done )
      call s_is_r ( word, g, lval )
      call word_nexrd ( text, word, done )
      call s_is_r ( word, b, lval )

    else if ( s_eqi ( word, 'DIFFUSECOLOR' ) ) then

      call word_nexrd ( text, word, done )
      call s_is_r ( word, r, lval )
      call word_nexrd ( text, word, done )
      call s_is_r ( word, g, lval )
      call word_nexrd ( text, word, done )
      call s_is_r ( word, b, lval )

    else if ( s_eqi ( word, 'SHININESS' ) ) then

      call word_nexrd ( text, word, done )
      call s_is_r ( word, r, lval )

    else if ( s_eqi ( word, 'SPECULARCOLOR' ) ) then

      call word_nexrd ( text, word, done )
      call s_is_r ( word, r, lval )
      call word_nexrd ( text, word, done )
      call s_is_r ( word, g, lval )
      call word_nexrd ( text, word, done )
      call s_is_r ( word, b, lval )

    else if ( s_eqi ( word, 'TRANSPARENCY' ) ) then

      call word_nexrd ( text, word, done )
      call s_is_r ( word, r, lval )

    else

    end if
!
!  MATERIALBINDING
!
!    OVERALL: Whole object has same material;
!    PER_FACE: One material for each face;
!    PER_FACE_INDEXED: One material for each face, indexed;
!    PER_PART: One material for each part;
!    PER_PART_INDEXED: One material for each part, indexed;
!    PER_VERTEX: One material for each vertex;
!    PER_VERTEX_INDEXED: One material for each vertex, indexed.
!
  else if ( s_eqi ( level_name(level), 'MATERIALBINDING' ) ) then
  
    if ( word == '{' ) then

    else if ( word == '}' ) then
      level = nlbrack - nrbrack
    else if ( s_eqi ( word, 'VALUE' ) ) then
      call word_nexrd ( text, word, done )
!         material_binding = word
    else
      bad_num = bad_num + 1
      write ( *, * ) '  Bad data on level ' // trim ( level_name(level) )
      write ( *, '(a)' ) '  Bad word: ' // trim ( word )
      write ( *, * ) '  Line number: ', text_num
      write ( *, '(a)' ) trim ( text )
    end if
!
!  MATERIALINDEX
!
  else if ( s_eqi ( level_name(level), 'MATERIALINDEX' ) ) then
  
    if ( word == '[' ) then

      ivert = 0

    else if ( word == ']' ) then

      level = nlbrack - nrbrack
!
!  (indexedfaceset) MATERIALINDEX
!
    else if ( s_eqi ( level_name(level-1), 'INDEXEDFACESET' ) ) then

      lval = s_is_i ( word, jval )
 
      if ( lval ) then

        if ( jval == -1 ) then

          ivert = 0

        else

          ivert = ivert + 1

          if ( ivert == 1 ) then
            face_num2 = face_num2 + 1
          end if

          if ( face_num2 <= face_max ) then
            if ( jval /= -1 ) then
              jval = jval + cor3_num_old
            end if
            vertex_material(ivert,face_num2) = jval + OFFSET
          end if

        end if

      else
        go to 99
      end if
!
!  (indexedlineset) MATERIALINDEX
!
    else if ( s_eqi ( level_name(level-1), 'INDEXEDLINESET' ) ) then

      lval = s_is_i ( word, jval )
 
      if ( lval ) then

        line_num2 = line_num2 + 1
        if ( line_num2 <= line_max ) then
          if ( jval /= -1 ) then
            jval = jval + cor3_num_old
          end if
          line_material(line_num2) = jval + OFFSET
        end if

      else
        go to 99
      end if

    else
 
      lval = s_is_i ( word, jval )
 
      if ( lval ) then

      else
        go to 99
      end if

    end if
!
!  MATRIXTRANSFORM
!
  else if ( s_eqi ( level_name(level), 'MATRIXTRANSFORM' ) ) then

      if ( word == '{' ) then

      else if ( word == '}' ) then

        if ( irow /= 4 .or. icol .ne. 4 ) then
          write ( *, * ) 'Incomplete MatrixTransform!'
        end if

        level = nlbrack - nrbrack

      else if ( s_eqi ( word, 'MATRIX' ) ) then

        irow = 1
        icol = 0

      else

        call s_is_r ( word, rval, lval )

        if ( lval ) then

          icol = icol + 1
          if ( icol > 4 ) then
            icol = 1
            irow = irow + 1
            if ( irow > 4 ) then
              go to 99
            end if
          end if

          transform_matrix(irow,icol) = rval

        else
          go to 99
        end if

      end if
!
!  NORMAL
!
!  The field "VECTOR" may be followed by three numbers,
!  (handled here),  or by a square bracket, and sets of three numbers.
!
  else if ( s_eqi ( level_name(level), 'NORMAL' ) ) then
!
!  (vertexproperty) NORMAL
!  Just stick the normal vectors into NORMAL_TEMP for now,
!  retrieve them at COORDINDEX time.
!
    if ( s_eqi ( level_name(level-1), 'VERTEXPROPERTY' ) ) then

      if ( word == '[' ) then

        ixyz = 0

      else if ( word == ']' ) then

        level = nlbrack - nrbrack

      else

        call s_is_r ( word, rval, lval )
 
        if ( lval ) then
 
          ixyz = ixyz + 1
 
          if ( ixyz > 3 ) then
            ixyz = 1
          end if

          if ( ixyz == 1 ) then
            normal_temp_num = normal_temp_num + 1
          end if
 
          if ( normal_temp_num <= order_max*face_max ) then

            r3vec(ixyz) = rval

            if ( ixyz == 3 ) then
              call tmat_mxv ( transform_matrix, r3vec, r3vec )
              normal_temp(1:3,normal_temp_num) = r3vec(1:3)
            end if

          end if
 
        else
          go to 99
        end if

      end if
!
!  (anythingelse) NORMAL.
!
    else

      if ( word == '{' ) then

        ixyz = 0
        normal_temp_num = 0

      else if ( word == '}' ) then

        level = nlbrack - nrbrack

      else if ( s_eqi ( word, 'VECTOR' ) ) then

      else

        go to 99

      end if

    end if
!
!  NORMALBINDING 
!    OVERALL: Whole object has same normal.
!    PER_FACE: One normal per face;
!    PER_FACE_INDEXED: one normal per face, indexed;
!    PER_FACE: one normal per part;
!    PER_FACE_INDEXED: one normal per part, indexed.
!    PER_VERTEX: one normal per vertex;
!    PER_VERTEX_INDEXED: one normal per vertex, indexed.
!
  else if ( s_eqi ( level_name(level), 'NORMALBINDING' ) ) then
  
    if ( word == '{' ) then

    else if ( word == '}' ) then
      level = nlbrack - nrbrack
    else if ( s_eqi ( word, 'VALUE' ) ) then
      call word_nexrd ( text, word, done )
!         normal_binding = word
    else
      bad_num = bad_num + 1
      write ( *, * ) '  Bad data on level ' // trim ( level_name(level) )
      write ( *, '(a)' ) '  Bad word: ' // trim ( word )
      write ( *, * ) '  Line number: ', text_num
      write ( *, '(a)' ) trim ( text )
    end if
!
!  NORMALINDEX
!
  else if ( s_eqi ( level_name(level), 'NORMALINDEX' ) ) then
!
!  (indexedtrianglestripset) NORMALINDEX
!
    if ( s_eqi ( level_name(level-1), 'INDEXEDTRIANGLESTRIPSET' ) ) then

      lval = s_is_i ( word, jval )

      if ( lval ) then

      else if ( word == '[' ) then

      else if ( word == ']' ) then
        level = nlbrack - nrbrack
      end if
!
!  (indexedFaceSet) NORMALINDEX
!
    else 

      if ( word == '[' ) then
        ivert = 0
      else if ( word == ']' ) then
        level = nlbrack - nrbrack
      else
 
        lval = s_is_i ( word, jval )
 
        if ( lval ) then

          if ( jval == -1 ) then

            ivert = 0

          else

            ivert = ivert + 1

            if ( ivert == 1 ) then
              iface_num = iface_num + 1
            end if

            if ( iface_num <= face_max .and. jval <= normal_temp_num ) then
              vertex_normal(1:3,ivert,iface_num) = normal_temp(1:3,jval + OFFSET)
            end if

          end if

        else
          go to 99
        end if

      end if

    end if
!
!  (coordinate3) POINT
!
  else if ( s_eqi ( level_name(level), 'POINT' ) ) then
 
    if ( s_eqi ( level_name(level-1), 'COORDINATE3' ) ) then

      if ( word == '[' ) then

        ixyz = 0
        cor3_num_old = cor3_num
 
      else if ( word == ']' ) then

        level = nlbrack - nrbrack

      else
 
        call s_is_r ( word, rval, lval )

        if ( lval ) then
 
          ixyz = ixyz + 1
 
          if ( ixyz > 3 ) then
            ixyz = 1
          end if
 
          r3vec(ixyz) = rval

          if ( ixyz == 3 ) then

            cor3_num = cor3_num + 1
            call tmat_mxp ( transform_matrix, r3vec, r3vec )

            if ( cor3_num <= cor3_max ) then
              cor3(1:3,cor3_num) = r3vec(1:3)
            end if

          end if
 
        else
          go to 99
        end if
 
      end if
!
!  (texturecoordinate2) POINT
!
    else if ( s_eqi ( level_name(level-1), 'TEXTURECOORDINATE2' ) ) then

      if ( word == '[' ) then

        iuv = 0
        text_numure_temp = 0
 
      else if ( word == ']' ) then

        level = nlbrack - nrbrack

      else
 
        call s_is_r ( word, rval, lval )

        if ( lval ) then
 
          iuv = iuv + 1
 
          if ( iuv > 2 ) then
            iuv = 1
          end if
 
          texture_temp(iuv,text_numure_temp+1) = rval

          if ( iuv == 2 ) then
            text_numure_temp = text_numure_temp + 1
          end if
 
        else
          go to 99
        end if
 
      end if

    end if
!
!  RGB
!
  else if ( s_eqi ( level_name(level), 'RGB' ) ) then

    if ( s_eqi ( level_name(level-1), 'BASECOLOR' ) ) then

      if ( word == '[' ) then
 
        icolor = 0

      else if ( word == ']' ) then

        level = nlbrack - nrbrack

      else 

        call s_is_r ( word, rval, lval )
 
        if ( lval ) then
 
          icolor = icolor + 1
 
          if ( icolor > 3 ) then
            icolor = 1
          end if

          if ( icolor == 1 ) then
            color_num = color_num + 1
          end if

        else
          go to 99
        end if

      end if

    end if
!
!  SEPARATOR
!
  else if ( s_eqi ( level_name(level), 'SEPARATOR' ) ) then
 
    if ( word == '{' ) then
 
    else if ( word == '}' ) then

      level = nlbrack - nrbrack

    else
 
    end if
!
!  SHAPEHINTS
!
!  Read, but ignore.
!
  else if ( s_eqi ( level_name(level), 'SHAPEHINTS' ) ) then
  
    if ( word == '{' ) then

    else if ( word == '}' ) then

      level = nlbrack - nrbrack

    else if ( s_eqi ( word, 'CREASEANGLE' ) ) then

      call word_nexrd ( text, word, done )
      call s_is_r ( word, rval, lval )
 
    else if ( s_eqi ( word, 'FACETYPE' ) ) then

      call word_nexrd ( text, word, done )

    else if ( s_eqi ( word, 'SHAPETYPE' ) ) then

      call word_nexrd ( text, word, done )

    else if ( s_eqi ( word, 'VERTEXORDERING' ) ) then

      call word_nexrd ( text, word, done )

    else
      go to 99
    end if
!
!  TEXTURE2
!
  else if ( s_eqi ( level_name(level), 'TEXTURE2' ) ) then

    if ( word == '{' ) then

      texture_num = texture_num + 1

    else if ( word == '}' ) then

      level = nlbrack - nrbrack

    else if ( s_eqi ( word, 'BLENDCOLOR' ) ) then

      call word_nexrd ( text, word, done )
      call s_is_r ( word, r, lval )
      call word_nexrd ( text, word, done )
      call s_is_r ( word, g, lval )
      call word_nexrd ( text, word, done )
      call s_is_r ( word, b, lval )

    else if ( s_eqi ( word, 'FILENAME' ) ) then

      call word_nexrd ( text, word, done )
      call word_nexrd ( text, word, done )
      texture_name(texture_num) = word
      call word_nexrd ( text, word, done )

    else if ( s_eqi ( word, 'IMAGE' ) ) then

      go to 99

    else if ( s_eqi ( word, 'MODEL' ) ) then

      call word_nexrd ( text, word, done )

    else if ( s_eqi ( word, 'WRAPS' ) ) then

      call word_nexrd ( text, word, done )

    else if ( s_eqi ( word, 'WRAPT' ) ) then

      call word_nexrd ( text, word, done )

    else

    end if
!
!  TEXTURECOORDINATE2
!
  else if ( s_eqi ( level_name(level), 'TEXTURECOORDINATE2' ) ) then
 
    if ( word == '{' ) then

    else if ( word == '}' ) then

      level = nlbrack - nrbrack

    else if ( s_eqi ( word, 'POINT' ) ) then

    else 
      go to 99
    end if
!
!  TEXTURECOORDINATEBINDING
!
!    PER_VERTEX: One texture coordinate for each vertex;
!    PER_VERTEX_INDEXED: One texture coordinate for each vertex, indexed.
!
  else if ( s_eqi ( level_name(level), 'TEXTURECOORDINATEBINDING' ) ) then
  
    if ( word == '{' ) then

    else if ( word == '}' ) then
      level = nlbrack - nrbrack
    else if ( s_eqi ( word, 'VALUE' ) ) then
      call word_nexrd ( text, word, done )
!         texture_coordinate_binding = word
    else
      bad_num = bad_num + 1
      write ( *, * ) '  Bad data on level ' // trim ( level_name(level) )
      write ( *, '(a)' ) '  Bad word: ' // trim ( word )
      write ( *, * ) '  Line number: ', text_num
      write ( *, '(a)' ) trim ( text )
    end if
!
!  TEXTURECOORDINDEX
!
  else if ( s_eqi ( level_name(level), 'TEXTURECOORDINDEX' ) ) then

    if ( word == '[' ) then
      ivert = 0
      iface = 0
    else if ( word == ']' ) then
      level = nlbrack - nrbrack
    else
 
      lval = s_is_i ( word, jval )
 
      if ( lval ) then

        if ( jval == -1 ) then

          ivert = 0

        else

          ivert = ivert + 1

          if ( ivert == 1 ) then
            iface = iface + 1
          end if

          if ( iface <= face_max ) then
            vertex_tex_uv(1:2,ivert,iface) = texture_temp(1:2,jval + OFFSET)
          end if

        end if

      else
        go to 99
      end if

    end if
!
!  UKNOTVECTOR
!
  else if ( s_eqi ( level_name(level), 'UKNOTVECTOR' ) ) then

    if ( word == '[' ) then
      go to 20
    else if ( word == ']' ) then
      level = nlbrack - nrbrack
      go to 20
    else
      lval = s_is_i ( word, jval )
    end if
!
!  VECTOR
!
  else if ( s_eqi ( level_name(level), 'VECTOR' ) ) then

    if ( word == '[' ) then

    else if ( word == ']' ) then

      level = nlbrack - nrbrack
!
!  (normal) VECTOR.
!  For some reason, Joel's program is spewing out occasional
!  NAN normal vectors.  This should not halt the program.
!
    else if ( s_eqi ( level_name(level-1), 'NORMAL' ) ) then

      if ( word(1:13) == 'nan0x7ffffe00' ) then

        lval = .true.
        rval = 1.0 / sqrt ( 3.0 )
        normal_bad_num = normal_bad_num + 1

      else

        call s_is_r ( word, rval, lval )
 
      end if

      if ( lval ) then
 
        ixyz = ixyz + 1
 
        if ( ixyz > 3 ) then
          ixyz = 1
        end if

        if ( ixyz == 1 ) then
          normal_temp_num = normal_temp_num + 1
        end if
 
        r3vec(ixyz) = rval

        if ( ixyz == 3 ) then

          if ( normal_temp_num <= order_max * face_max ) then
            normal_temp(1:3,normal_temp_num) = r3vec(1:3)
          end if

        end if
 
      else
        go to 99
      end if

    end if
!
!  (vertexproperty) VERTEX
!
  else if ( s_eqi ( level_name(level), 'VERTEX' ) ) then

    if ( s_eqi ( level_name(level-1), 'VERTEXPROPERTY' ) ) then

      if ( word == '[' ) then

        ixyz = 0
        cor3_num_old = cor3_num

      else if ( word == ']' ) then

        level = nlbrack - nrbrack

      else

        call s_is_r ( word, rval, lval )
 
        if ( lval ) then
 
          ixyz = ixyz + 1
 
          if ( ixyz > 3 ) then
            ixyz = 1
          end if
 
          if ( cor3_num+1 <= cor3_max ) then
            cor3(ixyz,cor3_num+1) = rval
          end if
 
          if ( ixyz == 3 ) then
            cor3_num = cor3_num + 1
          end if
 
        else
          go to 99
        end if

      end if

    end if
!
!  (indexedtrianglestripset) VERTEXPROPERTY
!
  else if ( s_eqi ( level_name(level), 'VERTEXPROPERTY' ) ) then

    if ( word == '{' ) then

    else if ( word == '}' ) then

      level = nlbrack - nrbrack

    else if ( s_eqi ( word, 'VERTEX' ) ) then

    else if ( s_eqi ( word, 'MATERIALBINDING' ) ) then

      call word_nexrd ( text, word, done )
!         material_binding = word

    else if ( s_eqi ( word, 'NORMAL' ) ) then

      ixyz = 0

    else if ( s_eqi ( word, 'NORMALBINDING' ) ) then

      call word_nexrd ( text, word, done )
!         normal_binding = word

    else
      go to 99
    end if
!
!  VKNOTVECTOR
!
  else if ( s_eqi ( level_name(level), 'VKNOTVECTOR' ) ) then
    
    if ( word == '[' ) then
      go to 20
    else if ( word == ']' ) then
      level = nlbrack - nrbrack
      go to 20
    else
      lval = s_is_i ( word, jval )
    end if
!
!  Any other word:
!
  else
 
  end if
 
  go to 20
!
!  Bad data
!
99    continue

  bad_num = bad_num + 1

  if ( bad_num <= 10 ) then
    write ( *, * ) ' '
    write ( *, * ) 'IV_READ - Warning!'
    write ( *, * ) '  Bad data on level ' // trim ( level_name(level) )
    write ( *, '(a)' ) '  Bad word: ' // trim ( word )
    write ( *, * ) '  Line number: ', text_num
    write ( *, '(a)' ) trim ( text )
  else
    write ( *, * ) ' '
    write ( *, * ) 'IV_READ - Fatal error!'
    write ( *, * ) '  Too many warnings!'
    return
  end if

  go to 10
!
!  End of information in file.
!
50    continue
!
!  Reset the transformation matrix to the identity, cause we
!  went ahead and applied it.
!
  call tmat_init ( transform_matrix )
!
!  Check the "materials" defining a line.  
!
!  If COORDINDEX is -1, so should be the MATERIALINDEX.  
!  If COORDINDEX is not -1, then the MATERIALINDEX shouldn't be either.
!
  do i = 1, line_num

    if ( line_dex(i) == -1 + OFFSET ) then
      line_material(i) = -1 + OFFSET
    else if ( line_material(i) == -1 + OFFSET ) then
      line_material(i) = material_num
    end if

  end do
!
!  Complain once about bad entries in normal vectors.
!
  if ( normal_bad_num > 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'IV_READ - Warning!'
    write ( *, * ) normal_bad_num, ' bad normal vector entries.'
  end if

  return
end
subroutine iv_write ( cor3, cor3_normal, face, face_order, filein_name, &
  fileout_name, iunit, line_dex, line_material, material_rgba, cor3_max, face_max, &
  line_max, material_max, order_max, texture_max, cor3_num, face_num, line_num, &
  material_num, texture_num, texture_name, vertex_material, vertex_tex_uv )
!
!*******************************************************************************
!
!! IV_WRITE writes graphics data to an Inventor file.
!
!
!  Example:
!
!     #Inventor V2.0 ascii
!
!     Separator {
!       Info {
!         string "Inventor file generated by IVREAD.
!         Original data in file cube.iv."
!       }
!       Separator {
!         LightModel {
!           model PHONG
!         }
!         MatrixTransform { matrix
!           0.9  0.0  0.0  0.0
!           0.0 -0.9  0.0  0.0
!           0.0  0.0 -1.5  0.0
!           0.0  0.0  0.0  1.0
!         }
!         Material {
!           ambientColor  0.2 0.2 0.2
!           diffuseColor  [
!             0.8 0.8 0.8,
!             0.7 0.1 0.1,
!             0.1 0.8 0.2,
!           ]
!           emissiveColor 0.0 0.0 0.0
!           specularColor 0.0 0.0 0.0
!           shininess     0.2
!           transparency  [
!             0.0, 0.5, 1.0,
!           ]
!         }
!         Texture2 {
!           filename      "fred.rgb"
!           wrapS         REPEAT
!           wrapT         REPEAT
!           model         MODULATE
!           blendColor    0.0 0.0 0.0
!         }
!         TextureCoordinateBinding {
!           value PER_VERTEX_INDEXED
!         }
!         MaterialBinding {
!           value PER_VERTEX_INDEXED
!         }
!         NormalBinding {
!           value PER_VERTEX_INDEXED
!         }
!         ShapeHints {
!           vertexOrdering COUNTERCLOCKWISE
!           shapeType UNKNOWN_SHAPE_TYPE
!           faceType CONVEX
!           creaseAngle 6.28319
!         }
!         Coordinate3 {
!           point [
!                8.59816       5.55317      -3.05561,
!                8.59816       2.49756      0.000000E+00,
!                ...etc...
!                2.48695       2.49756      -3.05561,
!           ]
!         }
!         TextureCoordinate2 {
!           point [
!                0.0  1.0,
!                0.1, 0.8,
!                ...etc...
!                0.4  0.7,
!           ]
!         }
!         Normal {
!           vector [
!             0.71 0.71 0.0,
!             ...etc...
!             0.32 0.32 0.41,
!           ]
!         }
!
!         IndexedLineSet {
!           coordIndex [
!              0,    1,    2,   -1,    3,    4,    5,   -1,    7,    8,    9,
!            ...etc...
!            191,   -1,
!           ]
!           materialIndex [
!              0,    0,    0,   -1,    1,    1,    1,   -1,    2,    2,    2,
!            ...etc...
!             64,   -1,
!           ]
!         }
!
!         IndexedFaceSet {
!           coordIndex [
!              0,    1,    2,   -1,    3,    4,    5,   -1,    7,    8,    9,
!            ...etc...
!            191,   -1,
!           ]
!           materialIndex [
!              0,    0,    0,   -1,    1,    1,    1,   -1,    2,    2,    2,
!            ...etc...
!             64,   -1,
!           ]
!           normalIndex [
!              0,    0,    0,   -1,    1,    1,    1,   -1,    2,    2,    2,
!            ...etc...
!             64,   -1,
!           ]
!           textureCoordIndex [
!              0,    0,    0,   -1,    1,    1,    1,   -1,    2,    2,    2,
!            ...etc...
!             64,   -1,
!           ]
!         }
!
!         IndexedTriangleStripSet {
!           vertexProperty VertexProperty {
!             vertex [ x y z,
!                      ...
!                      x y z ]
!             normal [ x y z,
!                      ...
!                      x y z ]
!             materialBinding OVERALL
!             normalBinding PER_VERTEX_INDEXED
!           }
!           coordIndex [ i, j, k, l, m, -1, n, o, p, q, r, s, t, u, -1,
!             v, w, x, -1 ..., -1 ]
!           normalIndex -1
!         }
!
!       }
!     }
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
!    Input, real COR3(3,COR3_MAX), the coordinates of points.
!
!    Input, real COR3_NORMAL(3,COR3_MAX), normals at nodes.  
!
!    Input, integer FACE(ORDER_MAX,FACE_MAX), the nodes making faces.
!
!    Input, integer FACE_ORDER(FACE_MAX), the number of vertices per face.
!
!    Input, character ( len = 100 ) FILEIN_NAME, the name of the input file.
!
!    Input, character ( len = 100 ) FILEOUT_NAME, the name of the output file.
!
!    Input, integer IUNIT, the FORTRAN unit to which output is written.
!
!    Input, integer LINE_DEX(LINE_MAX), nodes forming a line, terminated by -1.
!
!    Input, integer LINE_MAT(LINE_MAX), the material of each line.
!
!    Input, real MATERIAL_RGBA(4,MATERIAL_MAX), material R, G, B and A values.
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
!    Input, integer COR3_NUM, the number of points.
!
!    Input, integer FACE_NUM, the number of faces.
!
!    Input, integer LINE_NUM, the number of line definition items.
!
!    Input, integer MATERIAL_NUM, the number of materials.
!
!    Input, integer TEXTURE_NUM, the number of textures.
!
!    Input, real TRANSFORM_MATRIX(4,4), the transformation matrix.
!
!    Input, integer VERTEX_MAT(ORDER_MAX,FACE_MAX), vertex materials.
!
!    Input, real VERTEX_TEX_UV(2,ORDER_MAX,FACE_MAX), vertex texture coordinates.
!
  integer, parameter :: OFFSET = 1

  integer cor3_max
  integer face_max
  integer line_max
  integer material_max
  integer order_max
  integer texture_max
!
  real cor3(3,cor3_max)
  real cor3_normal(3,cor3_max)
  integer cor3_num
  integer face(order_max,face_max)
  integer face_num
  integer face_order(face_max)
  character ( len = 100 ) filein_name
  character ( len = 100 ) fileout_name
  integer i
  integer icor3
  integer iface
  integer ihi
  integer ilo
  integer itemp
  integer iunit
  integer ivert
  integer j
  integer jtemp
  integer length
  integer line_dex(line_max)
  integer line_material(line_max)
  integer line_num
  integer material_num
  real material_rgba(4,material_max)
  integer ndx
  character ( len = 200 ) text
  integer text_num
  character ( len = 100 ) texture_name(texture_max)
  integer texture_num
  real transform_matrix(4,4)
  integer vertex_material(order_max,face_max)
  real vertex_tex_uv(2,order_max,face_max)
  character ( len = 20 ) word
!
  text_num = 0

  write ( iunit, '(a)' ) '#Inventor V2.0 ascii'
  write ( iunit, '(a)' ) ' '
  write ( iunit, '(a)' ) 'Separator {'
  write ( iunit, '(a)' ) '  Info {'
  write ( iunit, '(a)' ) '    string "' // trim ( fileout_name ) &
    // ' generated by IVREAD."'
  write ( iunit, '(a)' ) '    string "Original data in file ' // &
    trim ( filein_name ) // '."'
  write ( iunit, '(a)' ) '  }'
  write ( iunit, '(a)' ) '  Separator {'
  text_num = text_num + 8
!
!  LightModel:
!
!    BASE_COLOR ignores light sources, and uses only diffuse color
!      and transparency.  Even without normal vector information,
!      the object will show up.  However, you won't get shadow
!      and lighting effects.
!
!    PHONG uses the Phong lighting model, accounting for light sources
!      and surface orientation.  This is the default.  I believe
!      you need accurate normal vector information in order for this
!      option to produce nice pictures.
!
!    DEPTH ignores light sources, and calculates lighting based on
!      the location of the object within the near and far planes
!      of the current camera's view volume.
!
  write ( iunit, '(a)' ) '    LightModel {'
  write ( iunit, '(a)' ) '      model PHONG'
  write ( iunit, '(a)' ) '    }'
  text_num = text_num + 3
!
!  Transformation matrix.
!
  call tmat_init ( transform_matrix )

  write ( iunit, '(a)' ) '    MatrixTransform { matrix'
  do i = 1, 4
    write ( iunit, '(8x,4g14.6)' ) transform_matrix(i,1:4)
  end do
  write ( iunit, '(a)' ) '    }'
!
!  Material
!
  write ( iunit, '(a)' ) '    Material {'
  write ( iunit, '(a)' ) '      ambientColor  0.2 0.2 0.2'
  write ( iunit, '(a)' ) '      diffuseColor  ['
  text_num = text_num + 3

  do i = 1, material_num
    write ( iunit, '(8x,3f8.4,'','')' ) material_rgba(1:3,i)
    text_num = text_num + 1
  end do

  write ( iunit, '(a)' ) '      ]'
  write ( iunit, '(a)' ) '      emissiveColor 0.0 0.0 0.0'
  write ( iunit, '(a)' ) '      specularColor 0.0 0.0 0.0'
  write ( iunit, '(a)' ) '      shininess     0.2'
  write ( iunit, '(a)' ) '      transparency  ['
  text_num = text_num + 4

  do ilo = 1, material_num, 10
    ihi = min ( ilo + 9, material_num )
    write ( iunit, '(8x,10(f7.3,'',''))' ) 1.0 - material_rgba(4,ilo:ihi)
    text_num = text_num + 1
  end do

  write ( iunit, '(a)' ) '      ]'
  write ( iunit, '(a)' ) '    }'
  text_num = text_num + 2
!
!  TEXTURE2
!
!  FLAW: We can only handle one texture right now.
!
  if ( texture_num > 0 ) then
    write ( iunit, '(a)' ) '    Texture2 {'
    write ( iunit, '(a)' ) '      filename    "' // trim ( texture_name(1) ) &
      // '"'
    write ( iunit, '(a)' ) '      wrapS       REPEAT'
    write ( iunit, '(a)' ) '      wrapT       REPEAT'
    write ( iunit, '(a)' ) '      model       MODULATE'
    write ( iunit, '(a)' ) '      blendColor  0.0 0.0 0.0'
    write ( iunit, '(a)' ) '    }'
    text_num = text_num + 7
  end if
!
!  TextureCoordinateBinding
!
  write ( iunit, '(a)' ) '    TextureCoordinateBinding {'
  write ( iunit, '(a)' ) '      value PER_VERTEX_INDEXED'
  write ( iunit, '(a)' ) '    }'
  text_num = text_num + 3
!
!  MaterialBinding
!
  write ( iunit, '(a)' ) '    MaterialBinding {'
  write ( iunit, '(a)' ) '      value PER_VERTEX_INDEXED'
  write ( iunit, '(a)' ) '    }'
  text_num = text_num + 3
!
!  NormalBinding
!
!    PER_VERTEX promises that we will write a list of normal vectors
!    in a particular order, namely, the normal vectors for the vertices
!    of the first face, then the second face, and so on.
!
!    PER_VERTEX_INDEXED promises that we will write a list of normal vectors,
!    and then, as part of the IndexedFaceSet, we will give a list of
!    indices referencing this normal vector list.
!
  write ( iunit, '(a)' ) '    NormalBinding {'
  write ( iunit, '(a)' ) '      value PER_VERTEX_INDEXED'
  write ( iunit, '(a)' ) '    }'
  text_num = text_num + 3
!
!  ShapeHints
!
  write ( iunit, '(a)' ) '    ShapeHints {'
  write ( iunit, '(a)' ) '      vertexOrdering COUNTERCLOCKWISE'
  write ( iunit, '(a)' ) '      shapeType UNKNOWN_SHAPE_TYPE'
  write ( iunit, '(a)' ) '      faceType CONVEX'
  write ( iunit, '(a)' ) '      creaseAngle 6.28319'
  write ( iunit, '(a)' ) '    }'
  text_num = text_num + 6
!
!  Point coordinates.
!
  write ( iunit, '(a)' ) '    Coordinate3 {'
  write ( iunit, '(a)' ) '      point ['
  text_num = text_num + 2
 
  do icor3 = 1, cor3_num
    write ( text, '(3f12.4,'','')' ) cor3(1:3,icor3)
    call s_blanks_delete ( text )
    write ( iunit, '(8x,a)' ) trim ( text )
    text_num = text_num + 1
  end do

  write ( iunit, '(a)' ) '      ]'
  write ( iunit, '(a)' ) '    }'
  text_num = text_num + 2
!
!  Texture coordinates.
!
  if ( texture_num > 0 ) then

    write ( iunit, '(a)' ) '    TextureCoordinate2 {'
    write ( iunit, '(a)' ) '      point ['
    text_num = text_num + 2
 
    do iface = 1, face_num
      do ivert = 1, face_order(iface)
        write ( text, '(2f12.4,'','')' ) vertex_tex_uv(1:2,ivert,iface)
        call s_blanks_delete ( text )
        write ( iunit, '(8x,a)' ) trim ( text )
        text_num = text_num + 1
      end do
    end do

    write ( iunit, '(a)' ) '      ]'
    write ( iunit, '(a)' ) '    }'
    text_num = text_num + 2

  end if
!
!  Normal vectors.
!    Use the normal vectors associated with nodes.
!
  if ( face_num > 0 ) then

    write ( iunit, '(a)' ) '    Normal { '
    write ( iunit, '(a)' ) '      vector ['
    text_num = text_num + 2

    do icor3 = 1, cor3_num
      write ( text, '(3f12.4,'','')' ) cor3_normal(1:3,icor3)
      call s_blanks_delete ( text )
      write ( iunit, '(8x,a)' ) trim ( text )
      text_num = text_num + 1
    end do

    write ( iunit, '(a)' ) '      ]'
    write ( iunit, '(a)' ) '    }'
    text_num = text_num + 2

  end if
!
!  IndexedLineSet
!
  if ( line_num > 0 ) then

    write ( iunit, '(a)' ) '    IndexedLineSet {'
!
!  IndexedLineSet coordIndex
!
    write ( iunit, '(a)' ) '      coordIndex ['
    text_num = text_num + 2

    text = ' '
    length = 0

    do j = 1, line_num
 
      write ( word, '(i8,'','')' ) line_dex(j) - OFFSET
      call s_cat ( text, word, text )
      length = length + 1

      if ( line_dex(j)-OFFSET == -1 .or. length >= 10 .or. j >= line_num ) then

        call s_blanks_delete ( text )
        write ( iunit, '(8x,a)' ) trim ( text )
        text_num = text_num + 1
        text = ' '
        length = 0

      end if

    end do
 
    write ( iunit, '(a)' ) '      ]'
    text_num = text_num + 1
!
!  IndexedLineSet materialIndex.
!
    write ( iunit, '(a)' ) '      materialIndex ['
    text_num = text_num + 1

    text = ' '
    length = 0

    do j = 1, line_num
 
      write ( word, '(i8,'','')' ) line_material(j) - OFFSET

      call s_cat ( text, word, text )
      length = length + 1

      if ( line_material(j)-OFFSET == -1 .or. length >= 10 .or. j >= line_num ) then

        call s_blanks_delete ( text )
        write ( iunit, '(8x,a)' ) trim ( text )
        text_num = text_num + 1
        text = ' '
        length = 0

      end if
 
    end do
 
    write ( iunit, '(a)' ) '      ]'
    write ( iunit, '(a)' ) '    }'
    text_num = text_num + 2

  end if
!
!  IndexedFaceSet.
!
  if ( face_num > 0 ) then

    write ( iunit, '(a)' ) '    IndexedFaceSet {'
!
!  IndexedFaceSet coordIndex
!
    write ( iunit, '(a)' ) '      coordIndex ['
    text_num = text_num + 2

    text = ' '
    length = 0
 
    do iface = 1, face_num

      do ivert = 1, face_order(iface) + 1

        if ( ivert <= face_order(iface) ) then
          itemp = face(ivert,iface) - OFFSET
        else
          itemp = 0 - OFFSET
        end if

        write ( word, '(i8,'','')' ) itemp
  
        call s_cat ( text, word, text )
        length = length + 1

        if ( itemp == -1 .or. length >= 10 .or. &
           ( iface == face_num .and. ivert == face_order(iface) + 1 )  ) then

          call s_blanks_delete ( text )
          write ( iunit, '(8x,a)' ) trim ( text )
          text_num = text_num + 1
          text = ' '
          length = 0

        end if

      end do

    end do

    write ( iunit, '(a)' ) '      ]'
    text_num = text_num + 1
!
!  IndexedFaceSet normalIndex
!
    if ( texture_num > 0 ) then

      write ( iunit, '(a)' ) '      normalIndex ['
      text_num = text_num + 1

      text = ' '
      length = 0

      do iface = 1, face_num
  
        do ivert = 1, face_order(iface) + 1

          if ( ivert <= face_order(iface) ) then
            itemp = face(ivert,iface) - OFFSET
          else
            itemp = 0 - OFFSET
          end if
  
          write ( word, '(i8,'','')' ) itemp

          call s_cat ( text, word, text )
          length = length + 1

          if ( itemp == -1 .or. length >= 10 .or. &
             ( iface == face_num .and. ivert == face_order(iface) + 1 )  ) then

            call s_blanks_delete ( text )
            write ( iunit, '(8x,a)' ) trim ( text )
            text_num = text_num + 1
            text = ' '
            length = 0

          end if

        end do

      end do

      write ( iunit, '(a)' ) '      ]'
      text_num = text_num + 1

    end if
!
!  IndexedFaceSet textureCoordIndex
!
    write ( iunit, '(a)' ) '      textureCoordIndex ['
    text_num = text_num + 1

    text = ' '
    length = 0
    itemp = 0

    do iface = 1, face_num

      do ivert = 1, face_order(iface) + 1

        if ( ivert <= face_order(iface) ) then
          jtemp = itemp
          itemp = itemp + 1
        else
          jtemp = - 1
        end if
  
        write ( word, '(i8,'','')' ) jtemp

        call s_cat ( text, word, text )
        length = length + 1

        if ( jtemp == -1 .or. length >= 10 .or. &
           ( iface == face_num .and. ivert == face_order(iface) + 1 )  ) then

          call s_blanks_delete ( text )
          write ( iunit, '(8x,a)' ) trim ( text )
          text_num = text_num + 1
          text = ' '
          length = 0

        end if

      end do

    end do

    write ( iunit, '(a)' ) '      ]'
    text_num = text_num + 1
!
!  IndexedFaceSet materialIndex
!
    write ( iunit, '(a)' ) '      materialIndex ['
    text_num = text_num + 1

    text = ' '
    length = 0
    ndx = 0
 
    do iface = 1, face_num
 
      do ivert = 1, face_order(iface) + 1

        if ( ivert <= face_order(iface) ) then
          itemp = vertex_material(ivert,iface) - OFFSET
          ndx = ndx + 1
        else
          itemp = 0 - OFFSET
        end if

        write ( word, '(i8,'','')' ) itemp

        call s_cat ( text, word, text )
        length = length + 1

        if ( itemp == -1 .or. length >= 10 .or. &
           ( iface == face_num .and. ivert == face_order(iface) + 1 )  ) then

          call s_blanks_delete ( text )
          write ( iunit, '(8x,a)' ) trim ( text )
          text_num = text_num + 1
          text = ' '
          length = 0
  
        end if
 
      end do

    end do

    write ( iunit, '(a)' ) '      ]'
    text_num = text_num + 1

    write ( iunit, '(a)' ) '    }'
    text_num = text_num + 1

  end if
!
!  Close up the Separator node.
!
  write ( iunit, '(a)' ) '  }'
!
!  Close up the Separator node.
!
  write ( iunit, '(a)' ) '}'
 
  text_num = text_num + 2
!
!  Report.
!
  write ( *, * ) ' '
  write ( *, * ) 'IV_WRITE - Wrote ', text_num, ' text lines.'

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
subroutine ivec_reverse ( n, x )
!
!*******************************************************************************
!
!! IVEC_REVERSE reverses the elements of an integer vector.
!
!
!  Example:
!
!    Input:
!
!      N = 5, X = ( 11, 12, 13, 14, 15 ).
!
!    Output:
!
!      X = ( 15, 14, 13, 12, 11 ).
!
!  Modified:
!
!    06 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries in the array.
!
!    Input/output, integer X(N), the array to be reversed.
!
  integer n
!
  integer i
  integer temp
  integer x(n)
!
  do i = 1, n/2
    temp = x(i)
    x(i) = x(n+1-i)
    x(n+1-i) = temp
  end do

  return
end
subroutine ivec_rotate ( n, m, x )
!
!*******************************************************************************
!
!! IVEC_ROTATE rotates an object in place.
!
!
!  Example:
!
!    Input:
!
!      N = 5, M = 2
!      X    = ( 1, 2, 3, 4, 5 )
!
!    Output:
!
!      X    = ( 4, 5, 1, 2, 3 ).
!
!  Modified:
!
!    07 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of objects.
!
!    Input, integer M, the number of positions to the right that
!    each element should be moved.  Elements that shift pass position
!    N "wrap around" to the beginning of the array.
!
!    Input/output, integer X(N), the array to be rotated.
!
  integer n
!
  integer iget
  integer i_modp
  integer iput
  integer istart
  integer m
  integer mcopy
  integer nset
  integer temp
  integer x(n)
!
!  Force M to be positive, between 0 and N-1.
!
  mcopy = i_modp ( m, n )

  if ( mcopy == 0 ) then
    return
  end if

  istart = 0
  nset = 0

10    continue

  istart = istart + 1

  if ( istart > n ) then
    return
  end if

  temp = x(istart)
  iget = istart
!
!  Copy the new value into the vacated entry.
!
20    continue

  iput = iget

  iget = iget - mcopy
  if ( iget < 1 ) then
    iget = iget + n
  end if

  if ( iget /= istart ) then
    x(iput) = x(iget)
    nset = nset + 1
    go to 20
  end if

  x(iput) = temp
  nset = nset + 1

  if ( nset < n ) then
    go to 10
  end if

  return
end
function lcon ( chr )
!
!*******************************************************************************
!
!! LCON reports whether a character is a control character or not.
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
!    Input, character CHR, is the character to be tested.
!
!    Output, logical LCON, TRUE if CHR is a control character, and
!    FALSE otherwise.
!
  character chr
  logical lcon
!
  if ( ichar ( chr ) <= 31 .or. ichar ( chr ) >= 127 ) then
    lcon = .true.
  else
    lcon = .false.
  end if

  return
end
subroutine news
!
!*******************************************************************************
!
!! NEWS prints out news (old and new) about the program.
!
!
!  Modified:
!
!    06 October 1999
!
!  Author:
!
!    John Burkardt
!
  write ( *, * ) ' '
  write ( *, * ) 'NEWS:'
  write ( *, * ) '  This is a list of changes to the program:'
  write ( *, * ) ' '
  write ( *, * ) '  14 June 2000:'
  write ( *, * ) '    FORTRAN90 conversion.'
  write ( *, * ) '  06 October 1999:'
  write ( *, * ) '    Trying to retrieve PS_WRITE from old code.'
  write ( *, * ) '  03 August 1999:'
  write ( *, * ) '    Corrected TMAT_ROT_VECTOR.'
  write ( *, * ) '  29 June 1999:'
  write ( *, * ) '    IV_READ and IV_WRITE support UV textures.'
  write ( *, * ) '  24 June 1999:'
  write ( *, * ) '    HRC_WRITE and TXT_WRITE write texture info.'
  write ( *, * ) '    BYU_WRITE writes Movie.BYU format.'
  write ( *, * ) '  23 June 1999:'
  write ( *, * ) '    HRC_READ and HRC_WRITE use vertex textures.'
  write ( *, * ) '  08 June 1999:'
  write ( *, * ) '    Added simple TECPLOT output.'
  write ( *, * ) '    Added Greg Hood "TRI" triangle output.'
  write ( *, * ) '  02 June 1999:'
  write ( *, * ) '    Adding material names.'
  write ( *, * ) '  26 May 1999:'
  write ( *, * ) '    Internal LINE_PRUNE option added.'
  write ( *, * ) '  24 May 1999:'
  write ( *, * ) '    STLA_WRITE automatically decomposes any'
  write ( *, * ) '    non-triangular faces before writing them.'
  write ( *, * ) '  23 May 1999:'
  write ( *, * ) '    DXF_WRITE can output faces.'
  write ( *, * ) '  22 May 1999:'
  write ( *, * ) '    VLA output includes line versions of faces.'
  write ( *, * ) '  04 May 1999:'
  write ( *, * ) '    SMF_READ/WRITE support face/node material,'
  write ( *, * ) '    normal, texture coordinates.'
  write ( *, * ) '  01 May 1999:'
  write ( *, * ) '    Relaxation smoothing option added.'
  write ( *, * ) '  27 April 1999:'
  write ( *, * ) '    SMF_READ and SMF_WRITE handle SMF2.0 colors.'
  write ( *, * ) '  23 April 1999:'
  write ( *, * ) '    Better POV_WRITE material treatment.'
  write ( *, * ) '    FACE_MATERIAL vs VERTEX_MAT fixup.'
  write ( *, * ) '  21 April 1999:'
  write ( *, * ) '    Trying to get IV_READ to merge two files.'
  write ( *, * ) '  20 April 1999:'
  write ( *, * ) '    Added << option, trying to get OBJ_READ'
  write ( *, * ) '      to merge two files.'
  write ( *, * ) '  26 March 1999'
  write ( *, * ) '    Added RGB_TO_HUE routine;'
  write ( *, * ) '    Adding UCD_WRITE.'
  write ( *, * ) '  23 February 1999'
  write ( *, * ) '    In HRC_WRITE, made specular and ambient'
  write ( *, * ) '    material colors the same as diffuse.'
  write ( *, * ) '  15 February 1999'
  write ( *, * ) '    Trying to get grid lines in OOGL.'
  write ( *, * ) '  13 February 1999'
  write ( *, * ) '    Added face area calculation.'
  write ( *, * ) '  12 February 1999'
  write ( *, * ) '    Added new color scheme to IV_WRITE.'
  write ( *, * ) '  10 February 1999'
  write ( *, * ) '    HRC_READ should now be able to read '
  write ( *, * ) '    material data.'
  write ( *, * ) '  09 February 1999'
  write ( *, * ) '    HRC_WRITE now writes out material data.'
  write ( *, * ) '    Moving all RGB material information'
  write ( *, * ) '      into MATERIAL_RGBA, with other items'
  write ( *, * ) '      using pointers to index it.'
  write ( *, * ) '  08 February 1999'
  write ( *, * ) '    Adding "OOGL" file format for Greg Foss.'
  write ( *, * ) '  02 December 1998'
  write ( *, * ) '    Restored VRML write.'
  write ( *, * ) '    Set up simple hooks for texture map names.'
  write ( *, * ) '  19 November 1998'
  write ( *, * ) '    IV_WRITE uses PER_VERTEX normal binding.'
  write ( *, * ) '  18 November 1998'
  write ( *, * ) '    Added node-based normal vectors.'
  write ( *, * ) '  17 November 1998'
  write ( *, * ) '    Added face node ordering reversal option.'
  write ( *, * ) '  23 October 1998'
  write ( *, * ) '    Added polygon to triangle option.'
  write ( *, * ) '  20 October 1998'
  write ( *, * ) '    Added interactive scaling patch.'
  write ( *, * ) '    Inserted TMAT routines.'
  write ( *, * ) '  19 October 1998'
  write ( *, * ) '    SMF_READ and SMF_WRITE added.'
  write ( *, * ) '  12 October 1998'
  write ( *, * ) '    Added FACE_CHECK code.'
  write ( *, * ) '  08 October 1998'
  write ( *, * ) '    Added POV_WRITE;'
  write ( *, * ) '    Added SET_VERTEX_NORMAL;'
  write ( *, * ) '    Modified normal vector computation.'
  write ( *, * ) '  30 August 1998'
  write ( *, * ) '    Still trying to fix up normals, because'
  write ( *, * ) '    of OBJ_READ and OBJ_WRITE complications.'
  write ( *, * ) '  29 August 1998'
  write ( *, * ) '    OBJ_READ and OBJ_WRITE now handle normals,'
  write ( *, * ) '    and read and write normals to file.'
  write ( *, * ) '    OBJ_READ can handle // face format.'
  write ( *, * ) '  28 August 1998'
  write ( *, * ) '    STLA_READ and STLA_WRITE seem OK after'
  write ( *, * ) '    the normal changes.'
  write ( *, * ) '  27 August 1998'
  write ( *, * ) '    Trying better NORMAL storage approach.'
  write ( *, * ) '  21 August 1998'
  write ( *, * ) '    Trying to add HRC_READ.'
  write ( *, * ) '    TXT_WRITE improved.'
  write ( *, * ) '  20 August 1998'
  write ( *, * ) '    Adding linear splines to HRC_WRITE.'
  write ( *, * ) '  19 August 1998'
  write ( *, * ) '    Automatic normal computation for OBJ files.'
  write ( *, * ) '    SoftImage HRC output added.'
  write ( *, * ) '    Normal vector computation improved.'
  write ( *, * ) '  18 August 1998'
  write ( *, * ) '    Improved treatment of face materials and normals.'
  write ( *, * ) '  17 August 1998'
  write ( *, * ) '    The maximum number of vertices per face'
  write ( *, * ) '    was increased to 35.'
  write ( *, * ) '    The maximum input line length was increased'
  write ( *, * ) '    to 256 characters.'
  write ( *, * ) '  10 August 1998:'
  write ( *, * ) '    Output DXF files have a comment now.'
  write ( *, * ) '    OBJ_READ corrected line indexing problem.'
  write ( *, * ) '  24 July 1998:'
  write ( *, * ) '    INCHECK checks the input data.'
  write ( *, * ) '    DXF_READ suppresses duplicate points.'
  write ( *, * ) '    Removed grid routines.'
  write ( *, * ) '    LINES(2,*) -> LINE_DEX(), LINE_MAT().'
  write ( *, * ) '    PS and VRML output dropped.'
  write ( *, * ) '    OBJ_WRITE line output tightened up.'
  write ( *, * ) '  22 July 1998:'
  write ( *, * ) '    STLA_READ suppresses duplicate points.'
  write ( *, * ) '  21 July 1998:'
  write ( *, * ) '    VLA_READ suppresses duplicate points.'
  write ( *, * ) '    OBJ_WRITE outputs line data now.'
  write ( *, * ) '  15 July 1998:'
  write ( *, * ) '    Added STLA_READ and STLA_WRITE.'
  write ( *, * ) '  11 July 1998:'
  write ( *, * ) '    DXF_READ and DXF_WRITE use IV data.'
  write ( *, * ) '    Dropped XYZ data structures.'
  write ( *, * ) '  10 July 1998:'
  write ( *, * ) '    Dropped XYZ input/output option.'
  write ( *, * ) '    VLA_READ and VLA_WRITE use IV data.'
  write ( *, * ) '  08 July 1998:'
  write ( *, * ) '    Added OBJ_READ and OBJ_WRITE.'
  write ( *, * ) '    Set ORDER_MAX=4, to allow for quad faces.'
  write ( *, * ) '  05 July 1998:'
  write ( *, * ) '    Added RF command to reverse faces.'
  write ( *, * ) '    Fixed 0/1 index based problem for FACE.'
  write ( *, * ) '  04 July 1998:'
  write ( *, * ) '    Added CHECK command to examine a face.'
  write ( *, * ) '  03 July 1998:'
  write ( *, * ) '    Only converting data when necessary.'
  write ( *, * ) '    Alternate IV triangles have opposite sense.'
  write ( *, * ) '    NORMALS command will recompute normals.'
  write ( *, * ) '  02 July 1998:'
  write ( *, * ) '    Trying to write simple ASE files.'
  write ( *, * ) '  01 July 1998:'
  write ( *, * ) '    Tentative attempts to read new IV data,'
  write ( *, * ) '    INDEXEDTRIANGLESTRIP.'
  write ( *, * ) '  03 June 1998:'
  write ( *, * ) '    MATERIALINDEX works for IV faces AND lines.'
  write ( *, * ) '  02 June 1998:'
  write ( *, * ) '    Trying to sort out -1/0 LINES confusion.'
  write ( *, * ) '    VRML_WRITE uses Inventor data.'
  write ( *, * ) '    I got VRML_WRITE to do color lines.'
  write ( *, * ) '    I need to reconcile RGBCOLOR/FACE/NODE.'
  write ( *, * ) '  15 May 1998:'
  write ( *, * ) '    Set up RGBCOLOR for new color handling.'
  write ( *, * ) '  08 May 1998:'
  write ( *, * ) '    Preparing for IV PER_VERTEX_INDEXED.'
  write ( *, * ) '  06 May 1998:'
  write ( *, * ) '    Added "reverse normal" option.'
  write ( *, * ) '    Added 3D projected plane grid.'
  write ( *, * ) '  05 May 1998:'
  write ( *, * ) '    Added projection into 3D plane.'
  write ( *, * ) '    Added spherical grid lines.'
  write ( *, * ) '  04 May 1998:'
  write ( *, * ) '    Sphere projection set up.'
  write ( *, * ) '    The PostScript output seems OK.'
  write ( *, * ) '  30 April 1998:'
  write ( *, * ) '    Adding 2D PostScript output option.'
  write ( *, * ) '  23 April 1998:'
  write ( *, * ) '    ASE->IV surface information sorta works.'
  write ( *, * ) '  22 April 1998:'
  write ( *, * ) '    IVREAD accepts command line arguments.'
  write ( *, * ) '  21 April 1998:'
  write ( *, * ) '    IV_WRITE now writes a default material.'
  write ( *, * ) '  20 April 1998:'
  write ( *, * ) '    ASE_READ tries to read vertex color.'
  write ( *, * ) '  17 April 1998:'
  write ( *, * ) '    ASE_READ reads transform matrix.'
  write ( *, * ) '    Overhauled ReadIV routine.'
  write ( *, * ) '  16 April 1998:'
  write ( *, * ) '    Increased big array limits.'
  write ( *, * ) '    Adding ASE_READ.'
  write ( *, * ) '  15 April 1998:'
  write ( *, * ) '    VLA intensities parameterized in INTENSE.'
  write ( *, * ) '    Got VRML option to work.'
  write ( *, * ) '  14 April 1998:'
  write ( *, * ) '    Added VRML_WRITE.'
  write ( *, * ) '  13 April 1998:'
  write ( *, * ) '    Made program command driven.'
  write ( *, * ) '    Started projection option.'
  write ( *, * ) '  10 April 1998:'
  write ( *, * ) '    Added Min/Max coordinate print.'
  write ( *, * ) '    Compressed IV output.'

  return
end
subroutine node_relax ( cor3, cor3_new, cor3_number, face, face_order, &
  cor3_max, face_max, order_max, cor3_num, face_num )
!
!*******************************************************************************
!
!! NODE_RELAX smooths a shape by an averaging operation on the node positions.
!
!  
!  Modified:
!
!    30 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real COR3(3,MAXCOR3), the coordinates of the nodes.
!
!    Output, real COR3_NEW(3,MAXCOR3), the new, averaged coordinates of 
!    the nodes.
!
!    Workspace, integer COR3_NUMBER(MAXCOR3).  On output, COR3_NUMBER(I)
!    will contain the number of node neighbors of node I.
!
!    Input, integer FACE(ORDER_MAX,FACE_MAX), the nodes making faces. 
!
!    Input, integer FACE_ORDER(FACE_MAX), is the number of nodes
!    making up each face.
!
!    Input, integer FACE_MAX, the maximum number of faces.
!
!    Input, integer ORDER_MAX, is the maximum number of nodes that can
!    make up a face, required to dimension FACE.
!
!    Input, integer FACE_NUM, the number of faces.
!
  integer cor3_max
  integer face_max
  integer order_max
!
  real cor3(3,cor3_max)
  real cor3_new(3,cor3_max)
  integer cor3_num
  integer cor3_number(cor3_max)
  integer face(order_max,face_max)
  integer face_num
  integer face_order(face_max)
  integer i
  integer icor3
  integer iface
  integer inode
  integer ivert
  integer jnode
!
!  COR3_NEW will contain the new averaged coordinates.
!
  cor3_number(1:cor3_num) = 0
  cor3_new(1:3,1:cor3_num) = 0.0
!
!  Consider each edge.  Essentially, the edge (I,J) is a signal to
!  add the old coordinates of I to the new J coordinates, and vice versa.
!
!  Because we are using a face representation, many, perhaps all the
!  edges, will show up repeatedly, probably twice.  To keep the algorithm
!  simple, for now we will simply use an edge every time it shows up
!  in a face, which means that edges that occur in multiple faces
!  will be weighted more.
!
  do iface = 1, face_num

    inode = face(face_order(iface),iface)

    do ivert = 1, face_order(iface)

      jnode = inode
      inode = face(ivert,iface)

      cor3_number(inode) = cor3_number(inode) + 1
      cor3_number(jnode) = cor3_number(jnode) + 1

      cor3_new(1:3,jnode) = cor3_new(1:3,jnode) + cor3(1:3,inode)
      cor3_new(1:3,inode) = cor3_new(1:3,inode) + cor3(1:3,jnode)

    end do

  end do
!
!  Copy the new into the old.
!
  do icor3 = 1, cor3_num

    if ( cor3_number(icor3) /= 0 ) then
      cor3_new(1:3,icor3) = cor3_new(1:3,icor3) / real ( cor3_number(icor3) )
    end if

  end do

  return
end
subroutine node_to_vertex_material ( cor3_material, face, face_order, cor3_max, &
  face_max, order_max, face_num, vertex_material )
!
!*******************************************************************************
!
!! NODE_TO_VERTEX_MAT extends node material definitions to vertices.
!
!
!  Discussion:
!
!    A NODE is a point in space.
!    A VERTEX is a node as used in a particular face.
!    One node may be used as a vertex in several faces, or none.
!
!  Modified:
!
!    27 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer COR3_MATERIAL(COR3_MAX), the material index of each node.
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
!    Input/output, integer FACE_NUM, the number of faces.
!
!    Output, integer VERTEX_MAT(ORDER_MAX,FACE_MAX), vertex materials.
!
  integer cor3_max
  integer face_max
  integer order_max
!
  integer cor3_material(cor3_max)
  integer face(order_max,face_max)
  integer face_order(face_max)
  integer iface
  integer ivert
  integer node
  integer face_num
  integer vertex_material(order_max,face_max)
!
  do iface = 1, face_num
    do ivert = 1, face_order(iface)
      node = face(ivert,iface)
      vertex_material(ivert,iface) = cor3_material(node)
    end do
  end do

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
      vertex_material(ivert,face_num) = material_num
    end if

    if ( face_num <= face_max ) then

      face_material(face_num) = material_num
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
      line_material(line_num) = material_num
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

    material_num = material_num + 1

    if ( material_num <= material_max ) then
      material_name(material_num) = word
      material_rgba(1,material_num) = uniform_01_sample ( iseed )
      material_rgba(2,material_num) = uniform_01_sample ( iseed )
      material_rgba(3,material_num) = uniform_01_sample ( iseed )
      material_rgba(4,material_num) = 1.0
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
subroutine obj_write ( cor3, face, face_order, filein_name, fileout_name, &
  iunit, line_dex, cor3_max, face_max, line_max, order_max, cor3_num, &
  face_num, line_num, vertex_normal )
!
!*******************************************************************************
!
!! OBJ_WRITE writes graphics information to a WaveFront OBJ file.
!
!
!  Example:
!
!    #  magnolia.obj
!
!    mtllib ./vp.mtl
!
!    g Group001
!    v -3.269770 -39.572201 0.876128
!    v -3.263720 -39.507999 2.160890
!    ...
!    v 0.000000 -9.988540 0.000000
!    vn 0.0 1.0 0.0
!    vn 1.0 0.0 0.0
!    ...
!    vn 0.0 0.0 1.0
!    s 1
!    usemtl brownskn
!    f 8//1 9//2 11//3 10//4
!    f 12//5 13//6 15//7 14//8
!    ...
!    f 788//800 806//803 774//807
!
!  Modified:
!
!    30 August 1998
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
!    Input, character ( len = 100 ) FILEIN_NAME, the name of the input file.
!
!    Input, character ( len = 100 ) FILEOUT_NAME, the name of the output file.
!
!    Input, integer IUNIT, the FORTRAN unit to which output is written.
!
!    Input, integer LINE_DEX(LINE_MAX), nodes forming a line, terminated by -1.
!
!    Input, integer COR3_MAX, the maximum number of points.
!
!    Input, integer FACE_MAX, the maximum number of faces.
!
!    Input, integer LINE_MAX, the maximum number of line definition items.
!
!    Input, integer ORDER_MAX, the maximum number of vertices per face.
!
!    Input, integer COR3_NUM, the number of points.
!
!    Input, integer FACE_NUM, the number of faces.
!
!    Input, integer LINE_NUM, the number of line definition items.
!
!    Input, real VERTEX_NORMAL(3,ORDER_MAX,FACE_MAX), normals at vertices.  
!
  integer, parameter :: OFFSET = 1
!
  integer cor3_max
  integer face_max
  integer line_max
  integer order_max
!
  real cor3(3,cor3_max)
  integer cor3_num
  integer face(order_max,face_max)
  integer face_num
  integer face_order(face_max)
  character ( len = 100 ) filein_name
  character ( len = 100 ) fileout_name
  integer i
  integer iface
  integer ihi
  integer ilo
  integer indexvn
  integer iunit
  integer ivert
  integer j
  integer line_dex(line_max)
  integer line_num
  integer nl
  character ( len = 256 ) text
  integer text_num
  character ( len = 256 ) text2
  real vertex_normal(3,order_max,face_max)
  real w
!
  text_num = 0
  write ( iunit, '(a)' ) '# ' // trim ( fileout_name ) // ' created by IVREAD.'
  write ( iunit, '(a)' ) '# Original data in ' // trim ( filein_name ) // '.'
  write ( iunit, '(a)' ) ' '
  write ( iunit, '(a)' ) 'g Group001'
  write ( iunit, '(a)' ) ' '

  text_num = text_num + 5
!
!  V: vertex coordinates.
!
  w = 1.0
  do j = 1, cor3_num
    write ( text, '(a1,2x,4g14.6)' ) 'v', cor3(1:3,j), w
    call s_blanks_delete ( text )
    write ( iunit, '(a)' ) trim ( text )
    text_num = text_num + 1
  end do
!
!  VN: vertex face normal vectors.
!
  if ( face_num > 0 ) then
    write ( iunit, '(a)' ) ' '
    text_num = text_num + 1
  end if

  do iface = 1, face_num

    do ivert = 1, face_order(iface)

      write ( text, '(a2,2x,3f7.3)' ) 'vn', vertex_normal(1:3,ivert,iface)
      call s_blanks_delete ( text )
      write ( iunit, '(a)' ) trim ( text )
      text_num = text_num + 1

    end do

  end do
!
!  F: Faces, specified as 
!    vertex index/texture vertex index/normal index
!
  if ( face_num > 0 ) then
    write ( iunit, '(a)' ) ' '
    text_num = text_num + 1
  end if

  indexvn = 0

  do iface = 1, face_num

    text = 'f'

    do ivert = 1, face_order(iface)
      indexvn = indexvn + 1
      text2 = ' '     
      write ( text2(2:), '(i8, ''//'', i8 )' ) face(ivert,iface), indexvn
      call s_blank_delete ( text2(2:) )
      call s_cat ( text, text2, text )
    end do

    write ( iunit, '(a)' ) trim ( text )
    text_num = text_num + 1

  end do
!
!  L: lines, specified as a sequence of vertex indices.
!
  nl = 0

  if ( line_num > 0 ) then

    write ( iunit, '(a)' ) ' '
    text_num = text_num + 1
  
    ihi = 0

    ilo = ihi

20      continue

    ihi = ihi + 1

    if ( ihi < line_num ) then
      if ( line_dex(ihi) /= -1 + OFFSET ) then
        go to 20
      end if
    end if

    write ( iunit, '(a,20i8)' ) 'l', ( line_dex(i), i = ilo+1, ihi-1 )

    text_num = text_num + 1
    nl = nl + 1

    if ( ihi < line_num ) then
      ilo = ihi
      go to 20
    end if

  end if
!
!  Report.
!
  write ( *, * ) ' '
  write ( *, * ) 'OBJ_WRITE - Wrote ', text_num, ' text lines.'

  return
end
subroutine object_build ( face, face_object, face_order, face_rank, &
  face_tier, order_max, face_num, object_num )
!
!*******************************************************************************
!
!! OBJECT_BUILD builds edge-connected "objects" out of polygonal faces.
!
!
!  Modified:
!
!    14 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer FACE(ORDER_MAX,FACE_NUM), the nodes making faces.
!
!    Output, integer FACE_OBJECT(FACE_NUM), describes the objects.
!    FACE_OBJECT(I) is the index of the edge-connected "object" to 
!    which face I belongs.
!
!    Input, integer FACE_ORDER(FACE_NUM), the number of nodes per face.
!
!    Output, integer FACE_RANK(FACE_NUM), is an ordered list of faces.
!    FACE_RANK(1) is the index of the face in the first tier of the 
!    first object, followed by second tier faces, and so on until
!    object one is complete.  Object two follows, and so on.
!
!    Output, integer FACE_TIER(FACE_NUM).  FACE_TIER(I) is the "tier"
!    of face I in its object.  The seed of the object has tier 1,
!    the neighbors of the seed have tier 2, and so on.
!
!    Input, integer ORDER_MAX, is the maximum number of nodes that can
!    make up a face, required to dimension FACE.
!
!    Input, integer FACE_NUM, the number of faces.
!
!    Output, integer OBJECT_NUM, the number of objects.
!
  logical, parameter :: DEBUG = .false.
!
  integer order_max
  integer face_num
!
  integer face(order_max,face_num)
  integer face_object(face_num)
  integer face_order(face_num)
  integer face_rank(face_num)
  integer face_tier(face_num)
  integer i
  integer iface
  integer ihi
  integer ihi_next
  integer ilo
  integer ilo_next
  integer iprint
  integer irank
  integer jface
  integer jrank
  integer object_num
  integer seed
  integer tier
  integer touch
!
!  Initialization.
!
  iprint = 1
  object_num = 0

  if ( face_num <= 0 ) then
    return
  end if

  do iface = 1, face_num
    face_object(iface) = 0
    face_rank(iface) = 0
    face_tier(iface) = 0
  end do

  irank = 0

  seed = 1
!
!  Begin the next object, seeded with face SEED.
!
10    continue

  tier = 1

  object_num = object_num + 1
  irank = irank + 1
  jrank = irank
  if ( irank == iprint .or. irank == face_num ) then
    write ( *, * ) irank, seed
    iprint = 2 * iprint
  end if

  face_rank(irank) = seed
  face_tier(seed) = tier
  face_object(seed) = object_num

  ilo = irank
  ihi = irank
!
!  Begin the next tier of faces, which are neighbors of faces we
!  found in the previous tier.
!
20    continue

  tier = tier + 1

  ilo_next = ihi + 1
  ihi_next = ihi

  do jface = 1, face_num

    if ( face_tier(jface) == 0 ) then

      do i = ilo, ihi

        iface = face_rank(i)

        call face_touch ( face, face_order, order_max, face_num, iface, &
          jface, touch )

        if ( touch /= 0 ) then
          if ( DEBUG ) then
            write ( *, * ) 'Touching faces: ', iface, jface
          end if
          ihi_next = ihi_next + 1
          irank = irank + 1
          if ( irank == iprint .or. irank == face_num ) then
            write ( *, * ) irank, jface
            iprint = 2 * iprint
          end if
          face_rank(irank) = jface
          face_tier(jface) = tier
          face_object(jface) = object_num
          go to 30
        end if

      end do

    end if

30      continue

  end do

  if ( ihi_next >= ilo_next ) then
    ilo = ilo_next
    ihi = ihi_next
    go to 20
  end if

  write ( *, * ) 'Object ', object_num, ' uses ', irank + 1 - jrank, ' faces.'
  jrank = irank
!
!  No neighbors were found, so this object is complete.  
!  Search for an unused face, which will be the seed of the next object.
!
  do iface = 1, face_num

    if ( face_tier(iface) == 0 ) then
      seed = iface
      go to 10
    end if

  end do

  return
end
subroutine object_invert ( cor3, cor3_material, cor3_normal, face, face_material, &
  face_normal, face_order, material_name, material_rgba, cor3_max, face_max, &
  material_max, order_max, cor3_num, face_num, material_num, vertex_material, vertex_normal )
!
!*******************************************************************************
!
!! OBJECT_INVERT makes an inverted duplicate of the object.
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
!    Input/output, integer FACE(ORDER_MAX,FACE_MAX), the nodes making faces.
!
!    Input/output, integer FACE_MATERIAL(FACE_MAX), the material of each face.
!
!    Input/output, integer FACE_ORDER(FACE_MAX), the number of vertices per face.
!
!    Input/output, character ( len = 100 ) MATERIAL_NAME(MATERIAL_MAX), material names.
!
!    Input, integer FACE_MAX, the maximum number of faces.
!
!    Input, integer ORDER_MAX, the maximum number of vertices per face.
!
!    Output, integer FACE_NUM, the number of faces.
!
!    Output, integer MATERIAL_NUM, the number of materials.
!
!    Input/output, integer VERTEX_MAT(ORDER_MAX,FACE_MAX), vertex materials.
!
  real, parameter :: EPS = 0.01
!
  integer cor3_max
  integer face_max
  integer material_max
  integer order_max
!
  real cor3(3,cor3_max)
  integer cor3_material(cor3_max)
  real cor3_normal(3,cor3_max)
  integer cor3_num
  integer face(order_max,face_max)
  integer face_material(face_max)
  real face_normal(3,face_max)
  integer face_num
  integer face_order(face_max)
  integer icor3
  integer icor32
  integer iface
  integer iface2
  integer ivert
  integer ivert2
  integer j
  character ( len = 100 ) material_name(material_max)
  integer material_num
  real material_rgba(4,material_max)
  integer vertex_material(order_max,face_max)
  real vertex_normal(3,order_max,face_max)
!
!  Check.
!
  if ( 2 * face_num > face_max ) then
    write ( *, * ) ' '
    write ( *, * ) 'OBJECT_INVERT - Fatal error!'
    write ( *, * ) '  2 * FACE_NUM exceeds FACE_MAX.'
    return
  end if

  if ( 2 * cor3_num > cor3_max ) then
    write ( *, * ) ' '
    write ( *, * ) 'OBJECT_INVERT - Fatal error!'
    write ( *, * ) '  2 * COR3_NUM exceeds COR3_MAX.'
    return
  end if
!
!  If there aren't 5 materials, add some.
!
  if ( material_num < 5 ) then

    material_num = material_num + 1
    if ( material_num <= material_max ) then
      material_name(material_num) = 'Red_Material'
      material_rgba(1,material_num) = 1.0
      material_rgba(2,material_num) = 0.0
      material_rgba(3,material_num) = 0.0
      material_rgba(4,material_num) = 1.0
      write ( *, * ) 'OBJECT_INVERT - Adding dummy red material'
    end if
  end if

  if ( material_num < 5 ) then

    material_num = material_num + 1
    if ( material_num <= material_max ) then
      material_name(material_num) = 'Green_Material'
      material_rgba(1,material_num) = 0.0
      material_rgba(2,material_num) = 1.0
      material_rgba(3,material_num) = 0.0
      material_rgba(4,material_num) = 1.0
      write ( *, * ) 'OBJECT_INVERT - Adding dummy green material'
    end if
  end if
!
!  Generate new points, displaced by EPS in the negative direction.
!
  do icor3 = 1, cor3_num

    icor32 = icor3 + cor3_num

    if ( cor3_material(icor3) == 1 ) then
      cor3_material(icor32) = 2
    else if ( cor3_material(icor3) == 3 ) then
      cor3_material(icor32) = 4
    else if ( cor3_material(icor3) == 4 ) then
      cor3_material(icor3) = 3
      cor3_material(icor32) = 4
    else if ( cor3_material(icor3) == 5 ) then
      cor3_material(icor3) = 3
      cor3_material(icor32) = 5
    end if

    cor3(1:3,icor32) = cor3(1:3,icor3) - EPS * cor3_normal(1:3,icor3)
    cor3_normal(1:3,icor32) = - cor3_normal(1:3,icor3)

  end do
!
!  Generate new faces.
!
  do iface = 1, face_num

    iface2 = face_num + iface
    face_order(iface2) = face_order(iface)

    if ( face_material(iface) == 1 ) then
      face_material(iface2) = 2
    else if ( face_material(iface) == 3 ) then
      face_material(iface2) = 4
    else if ( face_material(iface) == 4 ) then
      face_material(iface) = 3
      face_material(iface2) = 4
    else if ( face_material(iface) == 5 ) then
      face_material(iface) = 3
      face_material(iface2) = 5
    end if

    do ivert = 1, face_order(iface)

      ivert2 = face_order(iface) + 1 - ivert
      face(ivert2,iface2) = face(ivert,iface) + cor3_num

      if ( vertex_material(ivert,iface) == 1 ) then
        vertex_material(ivert2,iface2) = 2
      else if ( vertex_material(ivert,iface) == 3 ) then
        vertex_material(ivert2,iface2) = 4
      else if ( vertex_material(ivert,iface) == 4 ) then
        vertex_material(ivert,iface) = 3
        vertex_material(ivert2,iface2) = 4
      else if ( vertex_material(ivert,iface) == 5 ) then
        vertex_material(ivert,iface) = 3
        vertex_material(ivert2,iface2) = 5
      end if

      vertex_normal(1:3,ivert2,iface2) = - vertex_normal(1:3,ivert,iface)

    end do

    face_normal(1:3,iface2) = - face_normal(1:3,iface)

  end do

  cor3_num = 2 * cor3_num
  face_num = 2 * face_num

  write ( *, * ) ' '
  write ( *, * ) 'OBJECT_INVERT: Information:'
  write ( *, * ) '  Number of points = ', cor3_num
  write ( *, * ) '  Number of faces =  ', face_num

  return
end
subroutine oogl_read ( cor3, cor3_material, cor3_normal, face, &
  face_area, face_material, face_normal, face_order, ierror, iunit, material_name, &
  material_rgba, cor3_max, face_max, material_max, order_max, cor3_num, face_num, &
  material_num, text_num, ncol_oogl, nrow_oogl, vertex_material, vertex_normal )
!
!*******************************************************************************
!
!! OOGL_READ reads graphics information from a OOGL file.
!
!
!  Diagnostics:
!
!    Note that raw READ statements are used.  As written, the
!    code can't handle a blank line, or deal with a case where
!    information runs over to a new line.
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
!    {CMESH
!      5 3
!      0.0  0.0  0.0  0.94  0.70  0.15  1.00
!      1.0  0.0  0.0  0.94  0.70  0.15  1.00
!      2.0  0.0  0.0  0.94  0.70  0.15  1.00
!      3.0  0.0  0.0  0.94  0.70  0.15  1.00
!      4.0  0.0  0.0  0.94  0.70  0.15  1.00
!      0.0  1.0  0.0  0.94  0.70  0.15  1.00
!      1.0  1.0  0.0  0.94  0.70  0.15  1.00
!      2.0  1.0  0.0  0.94  0.70  0.15  1.00
!      3.0  1.0  0.0  0.94  0.70  0.15  1.00
!      4.0  1.0  0.0  0.94  0.70  0.15  1.00
!      0.0  2.0  0.0  0.94  0.70  0.15  1.00
!      1.0  2.0  0.0  0.94  0.70  0.15  1.00
!      2.0  2.0  0.0  0.94  0.70  0.15  1.00
!      3.0  2.0  0.0  0.94  0.70  0.15  1.00
!      4.0  2.0  0.0  0.94  0.70  0.15  1.00
!    }
!
!  Modified:
!
!    27 April 1999
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
!    Workspace, real FACE_AREA(FACE_MAX), the area of each face.
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
!    Input/output, character ( len = 100 ) MATERIAL_NAME(MATERIAL_MAX), material names.
!
!    Input/output, real MATERIAL_RGBA(4,MATERIAL_MAX), material R, G, B and A values.
!
!    Input, integer COR3_MAX, the maximum number of points.
!
!    Input, integer FACE_MAX, the maximum number of faces.
!
!    Input, integer MATERIAL_MAX, the maximum number of materials.
!
!    Input, integer ORDER_MAX, the maximum number of vertices per face.
!
!    Input/output, integer COR3_NUM, the number of points.
!
!    Input/output, integer FACE_NUM, the number of faces.
!
!    Input/output, integer MATERIAL_NUM, the number of materials.
!
!    Output, integer TEXT_NUM, the number of lines of text read from
!    the file.
!
!    ?, integer NCOL_OOGL, ?
!
!    ?, integer NROW_OOGL, ?
!
!    Input/output, integer VERTEX_MAT(ORDER_MAX,FACE_MAX), vertex materials.
!
!    ?, real VERTEX_NORMAL(3,ORDER_MAX,FACE_MAX), ?
!
  integer cor3_max
  integer face_max
  integer material_max
  integer order_max
!
  real a
  real b
  integer b2l2
  integer b2r2
  integer black
  character ( len = 4 ) char4
  logical clipping
  real cor3(3,cor3_max)
  integer cor3_material(cor3_max)
  real cor3_normal(3,cor3_max)
  integer cor3_num
  real dist
  integer face(order_max,face_max)
  real face_area(face_max)
  integer face_material(face_max)
  real face_normal(3,face_max)
  integer face_num
  integer face_order(face_max)
  real g
  logical grid
  integer i
  integer icor3
  logical identify
  integer ierror
  integer iface
  integer imat
  logical invert
  integer itemp
  integer iunit
  integer j
  integer jvert
  integer k
  character ( len = 100 ) material_name(material_max)
  integer material_num
  real material_rgba(4,material_max)
  integer nclip
  integer ncol_oogl
  integer ngold
  integer nrow_oogl
  real r
  real rgba(4)
  integer t2l2
  integer t2r2
  integer text_num
  integer vertex_material(order_max,face_max)
  real vertex_normal(3,order_max,face_max)
  integer white
  real x
  real y
  real z
!
  ierror = 0
  text_num = 0
!
!  30 March 1999
!  KLUDGE WARNING:
!  For our current purposes, we need to have the following materials defined.
!
!  Set up material #1 = GOLD.
!
  material_num = material_num + 1
  if ( material_num <= material_max ) then
    material_name(material_num) = 'Gold_Material'
    material_rgba(1,material_num) = 0.94
    material_rgba(2,material_num) = 0.70
    material_rgba(3,material_num) = 0.15
    material_rgba(4,material_num) = 1.0
  end if
!
!  Set up material #2 = DARK BLUE.
!
  material_num = material_num + 1
  if ( material_num <= material_max ) then
    material_name(material_num) = 'Dark_Blue_Material'
    material_rgba(1,material_num) = 0.24
    material_rgba(2,material_num) = 0.00
    material_rgba(3,material_num) = 0.85
    material_rgba(4,material_num) = 1.0
  end if
!
!  Set up material #3 = LIGHT BLUE.
!
  material_num = material_num + 1
  if ( material_num <= material_max ) then
    material_name(material_num) = 'Light_Blue_Material'
    material_rgba(1,material_num) = 0.24
    material_rgba(2,material_num) = 0.70
    material_rgba(3,material_num) = 0.85
    material_rgba(4,material_num) = 1.0
  end if
!
!  A node with grid coordinates (I,J) will be mapped to a
!  node number ( I - 1 ) * NCOL + J
!
  read ( iunit, *, end = 50 )
  text_num = text_num + 1

  read ( iunit, *, end = 50 ) ncol_oogl, nrow_oogl
  text_num = text_num + 1

  do i = 1, nrow_oogl
    do j = 1, ncol_oogl

      read ( iunit, *, end = 50 ) x, y, z, r, g, b, a
      text_num = text_num + 1

      cor3_num = cor3_num + 1
!
!  25 February 1999:
!  In the data we've seen, the first and last columns have almost
!  identical X,Y,Z coordinates and it might help if they were identical.
!
      if ( cor3_num <= cor3_max ) then

        cor3(1,cor3_num) = x
        cor3(2,cor3_num) = y
        cor3(3,cor3_num) = z

        if ( j == ncol_oogl ) then

          itemp = cor3_num + 1 - ncol_oogl

          dist = sqrt ( ( cor3(1,cor3_num) - cor3(1,itemp) )**2 &
                      + ( cor3(2,cor3_num) - cor3(2,itemp) )**2 &
                      + ( cor3(3,cor3_num) - cor3(3,itemp) )**2 )

          if ( dist > 0.0 .and. dist < 0.000001 ) then
            cor3(1,cor3_num) = cor3(1,itemp)
            cor3(2,cor3_num) = cor3(2,itemp)
            cor3(3,cor3_num) = cor3(3,itemp)
          end if

        end if

      end if
!
!  Some of the input data has had RGBA values that are negative.
!  Do not allow this.
!
      rgba(1) = min ( max ( 0.0, r ), 1.0 )
      rgba(2) = min ( max ( 0.0, g ), 1.0 )
      rgba(3) = min ( max ( 0.0, b ), 1.0 )
      rgba(4) = min ( max ( 0.0, a ), 1.0 )
!
!  See if the RGBA values of this material match those of a material
!  that has already been defined.
!
      if ( material_num <= 1000 ) then
        call rcol_find ( 4, material_num, material_rgba, rgba, imat )
      else
        imat = 0
      end if

      if ( imat == 0 ) then

        material_num = material_num + 1

        if ( material_num <= material_max ) then

          call i_to_s_zero ( material_num, char4 )

          material_name(material_num) = 'Material_' // char4              
          material_rgba(1:4,material_num) = rgba(1:4)
          imat = material_num

        else

          imat = 0

        end if

      end if

      if ( cor3_num <= cor3_max ) then
        cor3_material(cor3_num) = imat
      end if
  
    end do
  end do

  read ( iunit, *, end = 50 )
  text_num = text_num + 1
!
!  Now set up the faces from the grid of points that were defined.
!
  clipping = .false.
  nclip = 0

  do i = 2, nrow_oogl
    do j = 2, ncol_oogl

      b2l2 = ( i - 2 ) * ncol_oogl + j - 1
      b2r2 = ( i - 2 ) * ncol_oogl + j
      t2l2 = ( i - 1 ) * ncol_oogl + j - 1
      t2r2 = ( i - 1 ) * ncol_oogl + j

      if ( .not. clipping  ) then

        face_num = face_num + 1
        face_material(face_num) = cor3_material(t2r2)
        face_order(face_num) = 3

        k = t2r2
        face(1,face_num) = k
        vertex_material(1,face_num) = cor3_material(t2r2)

        k = b2r2
        face(2,face_num) = k
        vertex_material(2,face_num) = cor3_material(b2r2)

        k = b2l2
        face(3,face_num) = k
        vertex_material(3,face_num) = cor3_material(b2l2)

        face_num = face_num + 1
        face_material(face_num) = cor3_material(t2l2)
        face_order(face_num) = 3

        k = t2l2
        face(1,face_num) = k
        vertex_material(1,face_num) = cor3_material(t2l2)

        k = t2r2
        face(2,face_num) = k
        vertex_material(2,face_num) = cor3_material(t2r2)

        k = b2l2
        face(3,face_num) = k
        vertex_material(3,face_num) = cor3_material(b2l2)

      else 

        ngold = 0
        if ( cor3_material(t2r2) == 3 ) then
          ngold = ngold + 1
        end if
        if ( cor3_material(b2r2) == 3 ) then
          ngold = ngold + 1
        end if
        if ( cor3_material(b2l2) == 3 ) then
          ngold = ngold + 1
        end if
        if ( cor3_material(t2l2) == 3 ) then
          ngold = ngold + 1
        end if

        if ( ngold >= 3 ) then

          face_num = face_num + 1
          face_material(face_num) = 3
          face_order(face_num) = 3

          k = t2r2
          face(1,face_num) = k
          vertex_material(1,face_num) = 3

          k = b2r2
          face(2,face_num) = k
          vertex_material(2,face_num) = 3

          k = b2l2
          face(3,face_num) = k
          vertex_material(3,face_num) = 3

          face_num = face_num + 1
          face_material(face_num) = 3
          face_order(face_num) = 3

          k = t2l2
          face(1,face_num) = k
          vertex_material(1,face_num) = 3

          k = t2r2
          face(2,face_num) = k
          vertex_material(2,face_num) = 3

          k = b2l2
          face(3,face_num) = k
          vertex_material(3,face_num) = 3

        else

          nclip = nclip + 2

        end if

      end if

    end do
  end do

  if ( nclip > 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'OOGL_READ:'
    write ( *, * ) '  "Clipped" ', nclip, ' faces.'
  end if
!
!  Set "identify" TRUE to identify the repeated points at the first
!  and last, and along the sides.
!
  identify = .true.

  if ( identify ) then

    do iface = 1, face_num
      do jvert = 1, face_order(iface)
 
        i = face(jvert,iface)
        if ( i <= ncol_oogl ) then
          face(jvert,iface) = 1
        else if ( i <= nrow_oogl * ncol_oogl .and. &
                  i > ( nrow_oogl - 1 ) * ncol_oogl ) then
          face(jvert,iface) = ( nrow_oogl - 1 ) * ncol_oogl + 1
        end if

      end do
    end do

    do iface = 1, face_num
      do jvert = 1, face_order(iface)
        i = face(jvert,iface)
        if ( i <= nrow_oogl * ncol_oogl .and. mod ( i, ncol_oogl ) == 0 ) then
          face(jvert,iface) = face(jvert,iface) - ncol_oogl + 1
        end if

      end do
    end do

  end if

  write ( *, * ) ' '
  write ( *, * ) 'OOGL_READ - Information:'
  write ( *, * ) '  Initial number of points = ', cor3_num
  write ( *, * ) '  Initial number of faces =  ', face_num

  invert = .false.
  grid = .false.

  if ( invert .or. grid) then
!
!  Set up the normal vector information.
!
    call vertex_normal_set ( cor3, face, face_order, cor3_max, &
      face_max, order_max, face_num, vertex_normal )

    call face_area_set ( cor3, face, face_area, face_order, &
      cor3_max, face_max, order_max, face_num )

    cor3_normal(1:3,1:cor3_num) = 0.0

    call cor3_normal_set ( cor3_normal, face, face_area, &
      face_order, cor3_max, face_max, order_max, face_num, vertex_normal )

  end if
!
!  Make the other side of the surface.
!
  if ( invert ) then
    call object_invert ( cor3, cor3_material, cor3_normal, face, face_material, &
      face_normal, face_order, material_name, material_rgba, cor3_max, &
      face_max, material_max, order_max, cor3_num, face_num, material_num, vertex_material, &
      vertex_normal )

  end if
!
!  Make the grid.
!
  if ( grid ) then

    black = 1
    white = 2

    call oogl_grid ( cor3, cor3_material, cor3_normal, face, face_material, face_order, &
      cor3_max, face_max, order_max, cor3_num, face_num, ncol_oogl, nrow_oogl, &
      black, white, invert, vertex_material )

  end if
!
!  06 April 1999
!  KLUDGE #2 WARNING:
!  For our current purposes, we need to guarantee that there 
!  is at least one node of each material.
!
  cor3_num = cor3_num + 1
  cor3(1,cor3_num) = cor3(1,1)
  cor3(2,cor3_num) = cor3(1,1)
  cor3(3,cor3_num) = cor3(1,1)
  cor3_material(cor3_num) = 1
  cor3_normal(1,cor3_num) = 1.0
  cor3_normal(2,cor3_num) = 0.0
  cor3_normal(3,cor3_num) = 0.0

  cor3_num = cor3_num + 1
  cor3(1,cor3_num) = cor3(1,1)
  cor3(2,cor3_num) = cor3(1,1)
  cor3(3,cor3_num) = cor3(1,1)
  cor3_material(cor3_num) = 2
  cor3_normal(1,cor3_num) = 1.0
  cor3_normal(2,cor3_num) = 0.0
  cor3_normal(3,cor3_num) = 0.0

  cor3_num = cor3_num + 1
  cor3(1,cor3_num) = cor3(1,1)
  cor3(2,cor3_num) = cor3(1,1)
  cor3(3,cor3_num) = cor3(1,1)
  cor3_material(cor3_num) = 3
  cor3_normal(1,cor3_num) = 1.0
  cor3_normal(2,cor3_num) = 0.0
  cor3_normal(3,cor3_num) = 0.0

  return
!
!  Unexpected end of information.
!
50    continue

  write ( *, * ) ' '
  write ( *, * ) 'OOGL_READ - Fatal error!'
  write ( *, * ) '  Unexpected end of information!'

  return
end
subroutine oogl_grid ( cor3, cor3_material, cor3_normal, face, face_material, &
  face_order, cor3_max, face_max, order_max, cor3_num, face_num, ncol_oogl, &
  nrow_oogl, black, white, invert, vertex_material )
!
!*******************************************************************************
!
!! OOGL_GRID adds a grid to an OOGL data file.
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
!    Output, real COR3(3,COR3_MAX), the coordinates of points.
!
!    Output, integer FACE(ORDER_MAX,FACE_MAX), the nodes making faces.
!
!    Output, integer FACE_MATERIAL(FACE_MAX), the material of each face.
!
!    Output, integer FACE_ORDER(FACE_MAX), the number of vertices per face.
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
!    Output, integer COR3_NUM, the number of points.
!
!    Output, integer FACE_NUM, the number of faces.
!
!    Output, integer MATERIAL_NUM, the number of materials.
!
!    Output, integer VERTEX_MAT(ORDER_MAX,FACE_MAX), vertex materials.
!
  real, parameter :: EPS = 0.0025
  integer, parameter :: GRID_NX_NUM = 20
  integer, parameter :: GRID_NY_NUM = 20
  real, parameter :: GRID_WIDTH = 0.008
!
  integer cor3_max
  integer face_max
  integer order_max
!
  integer base
  integer b2l2
  integer b2r1u
  integer b2r2
  integer b2r2u
  integer black
  real cor3(3,cor3_max)
  integer cor3_material(cor3_max)
  real cor3_normal(3,cor3_max)
  integer cor3_num
  integer face(order_max,face_max)
  integer face_material(face_max)
  integer face_num
  integer face_order(face_max)
  integer grid_nx
  integer grid_ny
  integer i
  logical invert
  integer j
  integer k
  integer ncol_oogl
  real norm
  integer nrow_oogl
  integer nvert
  integer t1l2u
  integer t1r2u
  integer t2l2
  integer t2l2u
  integer t2r1u
  integer t2r2
  integer t2r2u
  integer vertex_material(order_max,face_max)
  integer white
!
!  Determine the number of grid lines.
!
  grid_nx = ncol_oogl / GRID_NX_NUM
  grid_nx = min ( grid_nx, ncol_oogl - 1 )
  grid_nx = max ( grid_nx, 1 )

  grid_ny = nrow_oogl / GRID_NY_NUM
  grid_ny = min ( grid_ny, nrow_oogl - 1 )
  grid_ny = max ( grid_ny, 1 )

  write ( *, * ) ' '
  write ( *, * ) 'OOGL_GRID:'
  write ( *, * ) '  Grid spacing is ', grid_nx, ' by ', grid_ny
!
!  Do the grid lines that I think of as running along the "top" of
!  the affected faces.
!
!    T2L2U---- ---- T2R2U
!    T1L2U---- ---- T1R2U
!    |    .... ....    |
!    +--- ---- ---- ---+
!
!    T2L2 ---- ---- T2R2
!    |    .... ....    |
!    |    .... ....    |
!    B2L2 ---- ---- B2R2
!
  do i = grid_ny+1, nrow_oogl, grid_ny
    do j = 2, ncol_oogl

      b2l2 = ( i - 2 ) * ncol_oogl + j - 1
      b2r2 = ( i - 2 ) * ncol_oogl + j
      t2l2 = ( i - 1 ) * ncol_oogl + j - 1
      t2r2 = ( i - 1 ) * ncol_oogl + j

      cor3_num = cor3_num + 1
      t2l2u = cor3_num
      cor3(1:3,t2l2u) = cor3(1:3,t2l2) + EPS * cor3_normal(1:3,t2l2)
      cor3_material(cor3_num) = black

      cor3_num = cor3_num + 1
      t2r2u = cor3_num
      do k = 1, 3
        cor3(k,t2r2u) = cor3(k,t2r2) + EPS * cor3_normal(k,t2r2)
      end do
      cor3_material(cor3_num) = black

      cor3_num = cor3_num + 1
      t1l2u = cor3_num

      norm = 0.0
      do k = 1, 3
        norm = norm + ( cor3(k,b2l2) + EPS * cor3_normal(k,b2l2) &
          - cor3(k,t2l2) )**2
      end do
      norm = sqrt ( norm )
      if ( norm /= 0.0 ) then
        norm = 1.0 / norm
      end if

      do k = 1, 3
        cor3(k,t1l2u) = cor3(k,t2l2u) + GRID_WIDTH * ( cor3(k,b2l2) &
          + EPS * cor3_normal(k,b2l2) - cor3(k,t2l2) ) * norm
      end do

      cor3_material(cor3_num) = black


      cor3_num = cor3_num + 1
      t1r2u = cor3_num

      norm = 0.0
      do k = 1, 3
        norm = norm + ( cor3(k,b2r2) + EPS * cor3_normal(k,b2r2) &
          - cor3(k,t2r2) )**2
      end do
      norm = sqrt ( norm )
      if ( norm /= 0.0 ) then
        norm = 1.0 / norm
      end if

      do k = 1, 3
        cor3(k,t1r2u) = cor3(k,t2r2u) &
          + GRID_WIDTH * ( cor3(k,b2r2) + EPS * cor3_normal(k,b2r2) &
          - cor3(k,t2r2) ) * norm
      end do

      cor3_material(cor3_num) = black

      face_num = face_num + 1
      face_material(face_num) = black
      face_order(face_num) = 3
      nvert = 0

      nvert = nvert + 1
      k = t2l2u
      face(nvert,face_num) = k
      vertex_material(nvert,face_num) = black

      nvert = nvert + 1
      k = t2r2u
      face(nvert,face_num) = k
      vertex_material(nvert,face_num) = black

      nvert = nvert + 1
      k = t1l2u
      face(nvert,face_num) = k
      vertex_material(nvert,face_num) = black

      face_num = face_num + 1
      face_material(face_num) = black
      face_order(face_num) = 3
      nvert = 0

      nvert = nvert + 1
      k = t1l2u
      face(nvert,face_num) = k
      vertex_material(nvert,face_num) = black

      nvert = nvert + 1
      k = t2r2u
      face(nvert,face_num) = k
      vertex_material(nvert,face_num) = black

      nvert = nvert + 1
      k = t1r2u
      face(nvert,face_num) = k
      vertex_material(nvert,face_num) = black

    end do
  end do

  write ( *, * ) 'OOGL_GRID2: Done "top" grid.'
!
!  Do the grid lines that I think of as running along the "right" of
!  the affected faces.
!
!    +--- ---- T2R1U-T2R2U
!    |    ....    |    |
!    |    ....    |    |
!    +--- ---- B2R1U-B2R2U
!
!    T2L2 ---- ---- T2R2
!    |    .... ....    |
!    |    .... ....    |
!    B2L2 ---- ---- B2R2
!
  do j = grid_nx+1, ncol_oogl, grid_nx
    do i = 2, nrow_oogl

      b2l2 = ( i - 2 ) * ncol_oogl + j - 1
      b2r2 = ( i - 2 ) * ncol_oogl + j
      t2l2 = ( i - 1 ) * ncol_oogl + j - 1
      t2r2 = ( i - 1 ) * ncol_oogl + j

      cor3_num = cor3_num + 1
      b2r2u = cor3_num
      do k = 1, 3
        cor3(k,b2r2u) = cor3(k,b2r2) + EPS * cor3_normal(k,b2r2)
      end do
      cor3_material(cor3_num) = black

      cor3_num = cor3_num + 1
      t2r2u = cor3_num
      do k = 1, 3
        cor3(k,t2r2u) = cor3(k,t2r2) + EPS * cor3_normal(k,t2r2)
      end do
      cor3_material(cor3_num) = black


      cor3_num = cor3_num + 1
      t2r1u = cor3_num

      norm = 0.0
      do k = 1, 3
        norm = norm + ( cor3(k,t2l2) + EPS * cor3_normal(k,t2l2) &
          - cor3(k,t2r2) )**2
      end do
      norm = sqrt ( norm )
      if ( norm /= 0.0 ) then
        norm = 1.0 / norm
      end if

      do k = 1, 3
        cor3(k,t2r1u) = cor3(k,t2r2u) + GRID_WIDTH * ( cor3(k,t2l2) + &
          EPS * cor3_normal(k,t2l2) - cor3(k,t2r2) ) * norm
      end do

      cor3_material(cor3_num) = black

      cor3_num = cor3_num + 1
      b2r1u = cor3_num

      norm = 0.0
      do k = 1, 3
        norm = norm + ( cor3(k,b2l2) + EPS * cor3_normal(k,b2l2) &
          - cor3(k,b2r2) )**2
      end do
      norm = sqrt ( norm )
      if ( norm /= 0.0 ) then
        norm = 1.0 / norm
      end if

      do k = 1, 3
        cor3(k,b2r1u) = cor3(k,b2r2u) + GRID_WIDTH * ( cor3(k,b2l2) &
          + EPS * cor3_normal(k,b2l2) - cor3(k,b2r2) ) * norm
      end do

      cor3_material(cor3_num) = black



      face_num = face_num + 1
      face_material(face_num) = black
      face_order(face_num) = 3
      nvert = 0

      nvert = nvert + 1
      k = b2r1u
      face(nvert,face_num) = k
      vertex_material(nvert,face_num) = black

      nvert = nvert + 1
      k = t2r1u
      face(nvert,face_num) = k
      vertex_material(nvert,face_num) = black

      nvert = nvert + 1
      k = t2r2u
      face(nvert,face_num) = k
      vertex_material(nvert,face_num) = black

      face_num = face_num + 1
      face_material(face_num) = black
      face_order(face_num) = 3
      nvert = 0

      nvert = nvert + 1
      k = b2r1u
      face(nvert,face_num) = k
      vertex_material(nvert,face_num) = black

      nvert = nvert + 1
      k = t2r2u
      face(nvert,face_num) = k
      vertex_material(nvert,face_num) = black

      nvert = nvert + 1
      k = b2r2u
      face(nvert,face_num) = k
      vertex_material(nvert,face_num) = black

    end do
  end do

  write ( *, * ) ' '
  write ( *, * ) 'OOGL_GRID - Added black "grid lines".'
  write ( *, * ) '  Number of points = ', cor3_num
  write ( *, * ) '  Number of faces =  ', face_num
!
!  Do the grid lines that I think of as running along the "top" of
!  the affected faces.
!
!    T2L2U---- ---- T2R2U
!    T1L2U---- ---- T1R2U
!    |    .... ....    |
!    +--- ---- ---- ---+
!
!    T2L2 ---- ---- T2R2
!    |    .... ....    |
!    |    .... ....    |
!    B2L2 ---- ---- B2R2
!
  if ( .not. invert ) then
    return
  end if

  base = ncol_oogl * nrow_oogl
  write ( *, * ) 'BASE = ', base

  do i = grid_ny+1, nrow_oogl, grid_ny
    do j = 2, ncol_oogl

      b2l2 = base + ( i - 2 ) * ncol_oogl + j - 1
      b2r2 = base + ( i - 2 ) * ncol_oogl + j
      t2l2 = base + ( i - 1 ) * ncol_oogl + j - 1
      t2r2 = base + ( i - 1 ) * ncol_oogl + j

      cor3_num = cor3_num + 1
      t2l2u = cor3_num
      do k = 1, 3
        cor3(k,t2l2u) = cor3(k,t2l2) + EPS * cor3_normal(k,t2l2)
      end do
      cor3_material(cor3_num) = white

      cor3_num = cor3_num + 1
      t2r2u = cor3_num
      do k = 1, 3
        cor3(k,t2r2u) = cor3(k,t2r2) + EPS * cor3_normal(k,t2r2)
      end do
      cor3_material(cor3_num) = white

      cor3_num = cor3_num + 1
      t1l2u = cor3_num

      norm = 0.0
      do k = 1, 3
        norm = norm + ( cor3(k,b2l2) + EPS * cor3_normal(k,b2l2) &
          - cor3(k,t2l2) )**2
      end do
      norm = sqrt ( norm )
      if ( norm /= 0.0 ) then
        norm = 1.0 / norm
      end if

      do k = 1, 3
        cor3(k,t1l2u) = cor3(k,t2l2u) + GRID_WIDTH * ( cor3(k,b2l2) &
          + EPS * cor3_normal(k,b2l2) - cor3(k,t2l2) ) * norm
      end do

      cor3_material(cor3_num) = white


      cor3_num = cor3_num + 1
      t1r2u = cor3_num

      norm = 0.0
      do k = 1, 3
        norm = norm + ( cor3(k,b2r2) + EPS * cor3_normal(k,b2r2) &
          - cor3(k,t2r2) )**2
      end do
      norm = sqrt ( norm )
      if ( norm /= 0.0 ) then
        norm = 1.0 / norm
      end if

      do k = 1, 3
        cor3(k,t1r2u) = cor3(k,t2r2u) + GRID_WIDTH * ( cor3(k,b2r2) &
          + EPS * cor3_normal(k,b2r2) - cor3(k,t2r2) ) * norm
      end do

      cor3_material(cor3_num) = white

      face_num = face_num + 1
      face_material(face_num) = white
      face_order(face_num) = 3
      nvert = 0

      nvert = nvert + 1
      k = t2l2u
      face(nvert,face_num) = k
      vertex_material(nvert,face_num) = white

      nvert = nvert + 1
      k = t2r2u
      face(nvert,face_num) = k
      vertex_material(nvert,face_num) = white

      nvert = nvert + 1
      k = t1l2u
      face(nvert,face_num) = k
      vertex_material(nvert,face_num) = white

      face_num = face_num + 1
      face_material(face_num) = white
      face_order(face_num) = 3
      nvert = 0

      nvert = nvert + 1
      k = t1l2u
      face(nvert,face_num) = k
      vertex_material(nvert,face_num) = white

      nvert = nvert + 1
      k = t2r2u
      face(nvert,face_num) = k
      vertex_material(nvert,face_num) = white

      nvert = nvert + 1
      k = t1r2u
      face(nvert,face_num) = k
      vertex_material(nvert,face_num) = white

    end do
  end do

  write ( *, * ) ' '
  write ( *, * ) 'OOGL_GRID - Added white top "grid lines".'
  write ( *, * ) '  Number of points = ', cor3_num
  write ( *, * ) '  Number of faces =  ', face_num
!
!  Do the grid lines that I think of as running along the "right" of
!  the affected faces.
!
!    +--- ---- T2R1U-T2R2U
!    |    ....    |    |
!    |    ....    |    |
!    +--- ---- B2R1U-B2R2U
!
!    T2L2 ---- ---- T2R2
!    |    .... ....    |
!    |    .... ....    |
!    B2L2 ---- ---- B2R2
!
  do j = grid_nx+1, ncol_oogl, grid_nx
    do i = 2, nrow_oogl

      b2l2 = base + ( i - 2 ) * ncol_oogl + j - 1
      b2r2 = base + ( i - 2 ) * ncol_oogl + j
      t2l2 = base + ( i - 1 ) * ncol_oogl + j - 1
      t2r2 = base + ( i - 1 ) * ncol_oogl + j

      cor3_num = cor3_num + 1
      b2r2u = cor3_num

      do k = 1, 3
        cor3(k,b2r2u) = cor3(k,b2r2) + EPS * cor3_normal(k,b2r2)
      end do
      cor3_material(cor3_num) = white

      cor3_num = cor3_num + 1
      t2r2u = cor3_num
      do k = 1, 3
        cor3(k,t2r2u) = cor3(k,t2r2) + EPS * cor3_normal(k,t2r2)
      end do
      cor3_material(cor3_num) = white

      cor3_num = cor3_num + 1
      t2r1u = cor3_num

      norm = 0.0
      do k = 1, 3
        norm = norm + ( cor3(k,t2l2) + EPS * cor3_normal(k,t2l2) &
          - cor3(k,t2r2) )**2
      end do
      norm = sqrt ( norm )
      if ( norm /= 0.0 ) then
        norm = 1.0 / norm
      end if

      do k = 1, 3
        cor3(k,t2r1u) = cor3(k,t2r2u) + GRID_WIDTH * ( cor3(k,t2l2) &
          + EPS * cor3_normal(k,t2l2) - cor3(k,t2r2) ) * norm
      end do

      cor3_material(cor3_num) = white


      cor3_num = cor3_num + 1
      b2r1u = cor3_num

      norm = 0.0
      do k = 1, 3
        norm = norm + ( cor3(k,b2l2) + EPS * cor3_normal(k,b2l2) &
          - cor3(k,b2r2) )**2
      end do
      norm = sqrt ( norm )
      if ( norm /= 0.0 ) then
        norm = 1.0 / norm
      end if

      do k = 1, 3
        cor3(k,b2r1u) = cor3(k,b2r2u) + GRID_WIDTH * ( cor3(k,b2l2) &
          + EPS * cor3_normal(k,b2l2) - cor3(k,b2r2) ) * norm
      end do

      cor3_material(cor3_num) = white



      face_num = face_num + 1
      face_material(face_num) = white
      face_order(face_num) = 3
      nvert = 0

      nvert = nvert + 1
      k = b2r1u
      face(nvert,face_num) = k
      vertex_material(nvert,face_num) = white

      nvert = nvert + 1
      k = t2r1u
      face(nvert,face_num) = k
      vertex_material(nvert,face_num) = white

      nvert = nvert + 1
      k = t2r2u
      face(nvert,face_num) = k
      vertex_material(nvert,face_num) = white

      face_num = face_num + 1
      face_material(face_num) = white
      face_order(face_num) = 3
      nvert = 0

      nvert = nvert + 1
      k = b2r1u
      face(nvert,face_num) = k
      vertex_material(nvert,face_num) = white

      nvert = nvert + 1
      k = t2r2u
      face(nvert,face_num) = k
      vertex_material(nvert,face_num) = white

      nvert = nvert + 1
      k = b2r2u
      face(nvert,face_num) = k
      vertex_material(nvert,face_num) = white

    end do
  end do

  write ( *, * ) ' '
  write ( *, * ) 'OOGL_GRID - Added white "grid lines".'
  write ( *, * ) '  Number of points = ', cor3_num
  write ( *, * ) '  Number of faces =  ', face_num

  return
end
subroutine outfile ( filein_name, fileout_name, ierror, fileout_type )
!
!*******************************************************************************
!
!! OUTFILE determines the output filename and type.
!
!
!  Modified:
!
!    02 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = 100 ) FILEIN_NAME, the name of the input file.
!
!    Output, character ( len = 100 ) FILEOUT_NAME, the name of the output file.
!
!    Output, integer IERROR, an error flag.
!
!    Output, character ( len = 10 ) FILEOUT_TYPE, the type of the file, which is
!    set to the filename extension.  Typical values include
!    'ase', 'dxf', 'iv', 'obj', 'pov', 'ps', 'smf', 'stl', 'stla', 
!    'tec', 'tri', 'txt', 'vla', or 'wrl'.
!
  character ( len = 100 ) filein_name
  character ( len = 100 ) fileout_name
  character ( len = 10 ) fileout_type
  integer i1
  integer i2
  integer ierror
  integer ios
  logical s_eqi
!
  ierror = 0

  if ( filein_name == ' ' ) then
    ierror = 1
    write ( *, * ) 'OUTFILE - Error!'
    write ( *, * ) '  You must read a file IN before you can'
    write ( *, * ) '  write a file OUT.'
    return
  end if
 
  if ( fileout_name == ' ' ) then
 
    write ( *, * ) ' '
    write ( *, * ) 'OUTFILE:'
    write ( *, * ) '  Enter the output filename to be created,'
    write ( *, * ) '  or hit return if done.'
 
    read ( *, '(a)', iostat = ios ) fileout_name
 
    if ( ios /= 0 ) then
      write ( *, * ) ' '
      write ( *, * ) 'OUTFILE - Error!'
      write ( *, * ) '  The output file was not specified correctly.'
      ierror = ios
      fileout_name = ' '
      return
    end if
 
  end if
!
!  Determine the output file type.
!
  call file_name_ext_get ( fileout_name, i1, i2 )

  if ( i1 /= 0 ) then
    fileout_type = fileout_name(i1:i2)
  else
    fileout_type = ' '
  end if
 
  if ( .not. ( &
    s_eqi ( fileout_type, 'ASE'  ) .or. s_eqi ( fileout_type, 'BYU'  ) .or. &
    s_eqi ( fileout_type, 'DXF'  ) .or. s_eqi ( fileout_type, 'HRC'  ) .or. &
    s_eqi ( fileout_type, 'IV'   ) .or. s_eqi ( fileout_type, 'OBJ'  ) .or. &
    s_eqi ( fileout_type, 'POV'  ) .or. s_eqi ( fileout_type, 'PS'   ) .or. &
    s_eqi ( fileout_type, 'SMF'  ) .or. s_eqi ( fileout_type, 'STL'  ) .or. &
    s_eqi ( fileout_type, 'STLA' ) .or. s_eqi ( fileout_type, 'TEC'  ) .or. &
    s_eqi ( fileout_type, 'TRI'  ) .or. s_eqi ( fileout_type, 'TRIA' ) .or. &
    s_eqi ( fileout_type, 'TXT'  ) .or. s_eqi ( fileout_type, 'UCD'  ) .or. &
    s_eqi ( fileout_type, 'VLA'  ) .or. s_eqi ( fileout_type, 'WRL'  ) .or. &
    s_eqi ( fileout_type, 'XGL'  ) ) ) then

    write ( *, * ) ' '
    write ( *, * ) 'OutFile could not determine the file type.'
    write ( *, * ) '  The output file name is:'
    write ( *, '(a)' ) trim ( fileout_name )
    write ( *, * ) ' '
    write ( *, * ) '  The file type should occur after the period.'
    write ( *, * ) '  Please specify the file type you are using:'
    write ( *, * ) '  Enter "ase", "byu", "dxf", "hrc", "iv",'
    write ( *, * ) '  "obj", "pov", "ps", "smf", "stl", "tec",'
    write ( *, * ) '  "tri", "txt", "ucd", "vla", "wrl", or "xgl":'

    read ( *, '(a)' ) fileout_type
    call s_cap ( fileout_type )

    if ( .not. ( &
      s_eqi ( fileout_type, 'ASE'  ) .or. s_eqi ( fileout_type, 'BYU'  ) .or. &
      s_eqi ( fileout_type, 'DXF'  ) .or. s_eqi ( fileout_type, 'HRC'  ) .or. &
      s_eqi ( fileout_type, 'IV'   ) .or. s_eqi ( fileout_type, 'OBJ'  ) .or. &
      s_eqi ( fileout_type, 'POV'  ) .or. s_eqi ( fileout_type, 'PS'   ) .or. &
      s_eqi ( fileout_type, 'SMF'  ) .or. s_eqi ( fileout_type, 'STL'  ) .or. &
      s_eqi ( fileout_type, 'STLA' ) .or. s_eqi ( fileout_type, 'TEC'  ) .or. &
      s_eqi ( fileout_type, 'TRI'  ) .or. s_eqi ( fileout_type, 'TRIA' ) .or. &
      s_eqi ( fileout_type, 'TXT'  ) .or. s_eqi ( fileout_type, 'UCD'  ) .or. &
      s_eqi ( fileout_type, 'VLA'  ) .or. s_eqi ( fileout_type, 'WRL'  ) .or. &
      s_eqi ( fileout_type, 'XGL'  ) ) ) then
      ierror = 1
      write ( *, * ) ' '
      write ( *, * ) 'OUTFILE - Error!'
      write ( *, * ) '  The file type was not acceptable!'
      return
    end if

  end if

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
  real pi
!
  pi = 3.14159265358979323846264338327950288419716939937510

  return
end
subroutine plane_exp2imp_3d ( x1, y1, z1, x2, y2, z2, x3, y3, z3, a, b, c, d )
!
!*******************************************************************************
!
!! PLANE_EXP2IMP_3D converts an explicit plane to implicit form in 3D.
!
!
!  Definition:
!
!    The explicit form of a plane in 3D is
!
!      (X1,Y1,Z1), (X2,Y2,Z2), (X3,Y3,Z3).
!
!    The implicit form of a plane in 3D is
!
!      A * X + B * Y + C * Z + D = 0
!
!  Reference:
!
!    Adrian Bowyer and John Woodwark,
!    A Programmer's Geometry,
!    Butterworths, 1983.
!
!  Modified:
!
!    04 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X1, Y1, Z1, X2, Y2, X2, X3, Y3, Z3, are three points 
!    on the plane, which must be distinct, and not collinear.
!
!    Output, real A, B, C, D, coefficients which describe the plane.
!
  real a
  real b
  real c
  real d
  real x1
  real x2
  real x3
  real y1
  real y2
  real y3
  real z1
  real z2
  real z3

  a = ( y2 - y1 ) * ( z3 - z1 ) - ( z2 - z1 ) * ( y3 - y1 )
  b = ( z2 - z1 ) * ( x3 - x1 ) - ( x2 - x1 ) * ( z3 - z1 )
  c = ( x2 - x1 ) * ( y3 - y1 ) - ( y2 - y1 ) * ( x3 - x1 ) 
  d = - x2 * a - y2 * b - z2 * c 

  return
end
subroutine plane_imp_point_nearest_3d ( a, b, c, d, x, y, z, xn, yn, zn )
!
!*******************************************************************************
!
!! PLANE_IMP_POINT_NEAREST_3D: nearest point on a implicit plane to a point in 3D.
!
!
!  Definition:
!
!    The implicit form of a plane in 3D is:
!
!      A * X + B * Y + C * Z + D = 0
!
!  Modified:
!
!    16 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A, B, C, D, coefficients that define the plane as
!    the set of points for which A*X+B*Y+C*Z+D = 0.
!
!    Input, real X, Y, Z, the coordinates of the point.
!
!    Output, real XN, YN, ZN, the coordinates of the nearest point on
!    the plane.
! 
  real a
  real b
  real c
  real d
  real t
  real x
  real xn
  real y
  real yn
  real z
  real zn
!
  if ( a == 0.0 .and. b == 0.0 .and. c == 0.0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'PLANE_IMP_POINT_NEAREST_3D - Fatal error!'
    write ( *, * ) '  A = B = C = 0.'
    stop
  end if
!
!  The normal N to the plane is (A,B,C).
!
!  The line defined by (XN-X)/A = (YN-Y)/B = (ZN-Z)/C = T
!  goes through (X,Y,Z) and is parallel to N.
!
!  Solving for the point (XN,YN,ZN) we get
!
!    XN = A*T+X
!    YN = B*T+Y
!    ZN = C*T+Z
!
!  Now place these values in the equation for the plane:
!
!    A*(A*T+X) + B*(B*T+Y) + C*(C*T+Z) + D = 0
!
!  and solve for T:
!
!    T = (-A*X-B*Y-C*Z-D) / (A * A + B * B + C * C )
!
  t = - ( a * x + b * y + c * z + d ) / ( a * a + b * b + c * c )
 
  xn = x + a * t
  yn = y + b * t
  zn = z + c * t
 
  return
end
subroutine points_distance_3d ( dis, x1, y1, z1, x2, y2, z2 )
!
!*******************************************************************************
!
!! POINTS_DISTANCE_3D finds the distance between two points in 3D.
!
!
!  Modified:
!
!    27 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real DIS, the distance between the points.
!
!    Input, real X1, Y1, Z1, X2, Y2, Z2, determines the pair of points
!    (X1,Y1,Z1) and (X2,Y2,Z2) whose distance apart is be determined.
!
  real dis
  real x1
  real x2
  real y1
  real y2
  real z1
  real z2
!
  dis = sqrt ( ( x1 - x2 )**2 + ( y1 - y2 )**2 + ( z1 - z2 )**2 )
 
  return
end
subroutine poly_2_tri ( face, face_material, face_order, ierror, face_max, &
  order_max, face_num, vertex_material )
!
!*******************************************************************************
!
!! POLY_2_TRI converts a collection of polygons into a collection of triangles.
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
!    Input/output, integer FACE(ORDER_MAX,FACE_MAX), the nodes making faces.
!
!    Input/output, integer FACE_MATERIAL(FACE_MAX), the material of each face.
!
!    Output, integer IERROR, error flag.
!    0, no error.
!    1, the algorithm failed because FACE_MAX was too small.
!    2, the algorithm failed because there were faces of order < 3.
!    3, the algorithm failed because there were faces of order > ORDER_MAX.
!
!    Input, integer FACE_MAX, the maximum number of faces allowed.
!
!    Input, integer ORDER_MAX, the maximum number of vertices allowed per face.
!
!    Input/output, integer FACE_NUM, the number of faces.  This value is updated
!    on return.
!
!    Input/output, integer VERTEX_MAT(ORDER_MAX,FACE_MAX), vertex materials.
!
  integer face_max
  integer order_max
!
  integer face(order_max,face_max)
  integer face_material(face_max)
  integer face_num
  integer face_num2
  integer face_order(face_max)
  integer ierror
  integer iface
  integer iface_old
  integer ivert
  integer k
  integer vertex_material(order_max,face_max)
!
  ierror = 0
  face_num2 = 0

  do iface = 1, face_num

    if ( face_order(iface) < 3 ) then
      write ( *, * ) ' '
      write ( *, * ) 'POLY_2_TRI - Fatal error!'
      write ( *, * ) '  Face ', iface, ' is illegal.'
      write ( *, * ) '  Number of vertices is ', face_order(iface)
      ierror = 2
      return
    else if ( face_order(iface) > order_max ) then
      write ( *, * ) ' '
      write ( *, * ) 'POLY_2_TRI - Fatal error!'
      write ( *, * ) '  Face ', iface, ' is illegal.'
      write ( *, * ) '  Number of vertices is ', face_order(iface)
      write ( *, * ) '  ORDER_MAX is ', order_max
      ierror = 3
      return
    end if

    do ivert = 3, face_order(iface)
      face_num2 = face_num2 + 1
    end do

  end do

  if ( face_num2 > face_max ) then
    write ( *, * ) ' '
    write ( *, * ) 'POLY_2_TRI - Fatal error!'
    write ( *, * ) '  FACE_MAX is too small to replace all faces'
    write ( *, * ) '  by triangles.'
    write ( *, * ) '  FACE_MAX = ', face_max
    write ( *, * ) '  FACE_NUM2 = ', face_num2
    ierror = 1
    return
  end if

  iface_old = face_num
  k = face_order(iface_old)

  do iface = face_num2, 1, -1

    if ( k < 3 ) then
      iface_old = iface_old - 1
      k = face_order(iface_old)
    end if

    face_material(iface) = face_material(iface_old)
    face_order(iface) = 3
    face(1,iface) = face(1,iface_old)
    vertex_material(1,iface) = vertex_material(1,iface_old)
    do ivert = 2, 3
      face(ivert,iface) = face(k+ivert-3,iface_old)
      vertex_material(ivert,iface) = vertex_material(k+ivert-3,iface_old)
    end do

    k = k - 1

  end do

  face_num = face_num2

  return
end
subroutine pov_write ( cor3, face, face_material, face_order, filein_name, &
  fileout_name, iunit, material_rgba, cor3_max, face_max, material_max, &
  order_max, face_num, material_num, vertex_normal )
!
!*******************************************************************************
!
!! POV_WRITE writes graphics information to a POV file.
!
!
!  Example:
!
!    // cone.pov created by IVREAD.
!    // Original data in cone.iv.
!
!    #version 3.0
!    #include "colors.inc"
!    #include "shapes.inc"
!    global_settings { assumed_gamma 2.2 }
!
!    camera { 
!     right < 4/3, 0, 0>
!     up < 0, 1, 0 >
!     sky < 0, 1, 0 >
!     angle 20
!     location < 0, 0, -300 >
!     look_at < 0, 0, 0>
!    }
!
!    light_source { < 20, 50, -100 > color White }
!
!    background { color SkyBlue }
!
!    #declare Material001 = texture { 
!      pigment { color rgb < 0.8, 0.2, 0.2> } 
!      finish { ambient 0.2 diffuse 0.5 }
!    }
! 
!    #declare Material002 = texture { 
!      pigment { color rgb < 0.2, 0.2, 0.8> } 
!      finish { ambient 0.2 diffuse 0.5 }
!    }
!
!    mesh {
!      smooth_triangle { 
!        < 0.29, -0.29, 0.0>, < 0.0, 0.0, -1.0 >,
!        < 38.85, 10.03, 0.0>, < 0.0, 0.0, -1.0 >,
!        < 40.21, -0.29, 0.0>, <  0.0, 0.0, -1.0 > 
!        texture { Material002 } }
!        ...
!      smooth_triangle { 
!        <  0.29, -0.29, 70.4142 >, < 0.0,  0.0, 1.0 >,
!        <  8.56,  -2.51, 70.4142 >, < 0.0,  0.0, 1.0 >,
!        <  8.85, -0.29, 70.4142 >, < 0.0,  0.0, 1.0 > 
!        texture { Material001 } }
!    }
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
!    Input, real COR3(3,COR3_MAX), the coordinates of points.
!
!    Input, integer FACE(ORDER_MAX,FACE_MAX), the nodes making faces.
!
!    Input, integer FACE_MATERIAL(FACE_MAX), face materials.
!
!    Input, integer FACE_ORDER(FACE_MAX), the number of vertices per face.
!
!    Input, character ( len = 100 ) FILEIN_NAME, the name of the input file.
!
!    Input, character ( len = 100 ) FILEOUT_NAME, the name of the output file.
!
!    Input, integer IUNIT, the FORTRAN unit to which output is written.
!
!    Input, real MATERIAL_RGBA(4,MATERIAL_MAX), material R, G, B and A values.
!
!    Input, integer COR3_MAX, the maximum number of points.
!
!    Input, integer FACE_MAX, the maximum number of faces.
!
!    Input, integer MATERIAL_MAX, the maximum number of materials.
!
!    Input, integer ORDER_MAX, the maximum number of vertices per face.
!
!    Input, integer FACE_NUM, the number of faces.
!
!    Input, integer MATERIAL_NUM, the number of materials.
!
!    Input, real VERTEX_NORMAL(3,ORDER_MAX,FACE_MAX), normals at vertices.  
!
  integer cor3_max
  integer face_max
  integer material_max
  integer order_max
!
  real b
  character ( len = 4 ) char4
  character comma
  real cor3(3,cor3_max)
  integer face(order_max,face_max)
  integer face_material(face_max)
  integer face_num
  integer face_order(face_max)
  character ( len = 100 ) filein_name
  character ( len = 100 ) fileout_name
  real g
  integer i
  integer i_mat
  integer iunit
  integer j
  integer jj
  integer jlo
  integer k
  real material_rgba(4,material_max)
  integer material_num
  real r
  character ( len = 100 ) text
  integer text_num
  real vertex_normal(3,order_max,face_max)
!
  text_num = 0

  write ( iunit, '(a)' ) '// ' // trim ( fileout_name ) // ' created by IVREAD.'
  write ( iunit, '(a)' ) '// Original data in ' // trim ( filein_name ) // '.'

  text_num = text_num + 2
!
!  Initial declarations.
!
  write ( iunit, '(a)' ) ' '
  write ( iunit, '(a)' ) '#version 3.0'
  write ( iunit, '(a)' ) '#include "colors.inc"'
  write ( iunit, '(a)' ) '#include "shapes.inc"'
  write ( iunit, '(a)' ) 'global_settings { assumed_gamma 2.2 }'
  write ( iunit, '(a)' ) ' '
  write ( iunit, '(a)' ) 'camera { '
  write ( iunit, '(a)' ) ' right < 4/3, 0, 0>'
  write ( iunit, '(a)' ) ' up < 0, 1, 0 >'
  write ( iunit, '(a)' ) ' sky < 0, 1, 0 >'
  write ( iunit, '(a)' ) ' angle 20'
  write ( iunit, '(a)' ) ' location < 0, 0, -300 >'
  write ( iunit, '(a)' ) ' look_at < 0, 0, 0>'
  write ( iunit, '(a)' ) '}'
  write ( iunit, '(a)' ) ' '
  write ( iunit, '(a)' ) 'light_source { < 20, 50, -100 > color White }'
  write ( iunit, '(a)' ) ' '
  write ( iunit, '(a)' ) 'background { color SkyBlue }'

  text_num = text_num + 15
!
!  Declare RGB textures.
!
  do i = 1, material_num

    write ( iunit, '(a)' ) ' '

    call i_to_s_zero ( i, char4 )
    write ( iunit, '(a)' ) '#declare Material_' // char4 // ' = texture { '

    r = material_rgba(1,i)
    g = material_rgba(2,i)
    b = material_rgba(3,i)
    write ( iunit, '(a,f4.2,a,f4.2,a,f4.2,a)' ) '  pigment { color rgb < ', &
      r, ',', g, ',', b ,' > } '

    write ( iunit, '(a)' ) '  finish { ambient 0.2 diffuse 0.5 }'

    write ( iunit, '(a)' ) '}'

  end do
!
!  Write one big object.
!
  write ( iunit, '(a)' ) 'mesh {'
  text_num = text_num + 1
!
!  Do the next face.
!
  do i = 1, face_num
!
!  Break the face up into triangles, anchored at node 1.
!
    do jlo = 1, face_order(i) - 2

      write ( iunit, '(a)' )  '  smooth_triangle { '
      text_num = text_num + 1

      do j = jlo, jlo + 2

        if ( j == jlo ) then
          jj = 1
        else
          jj = j
        end if

        k = face(jj,i)

        if ( j < jlo + 2 ) then
          comma = ','
        else
          comma = ' '
        end if
 
        write ( text, '(a,3(f10.3,a),3(f6.2,a),a )' ) &
          '<', cor3(1,k), ',', cor3(2,k), ',', cor3(3,k), '>, <', &
          vertex_normal(1,jj,i), ',', vertex_normal(2,jj,i), ',', &
          vertex_normal(3,jj,i), '>', comma

        call s_blanks_delete ( text )
        write ( iunit, '(a)' ) trim ( text )
        text_num = text_num + 1

      end do

      i_mat = face_material(i)
      call i_to_s_zero ( i_mat, char4 )
      write ( iunit, '(a)' ) 'texture { Material_' // char4 // ' } }'
      text_num = text_num + 1

    end do

  end do

  write ( iunit, '(a)' ) '}'
  text_num = text_num + 1
!
!  Report.
!
  write ( *, * ) ' '
  write ( *, * ) 'POV_WRITE - Wrote ', text_num, ' text lines.'

  return
end
subroutine project_2d ( cor2, cor3, ierror, cor2_max, cor3_max, cor2_num, &
  cor3_num )
!
!*******************************************************************************
!
!! PROJECT_2D projects 3D data to 2D based on user choices.
!
!
!  Discussion:
!
!    Projections include:
!
!      drop X coordinate, display YZ;
!      drop Y coordinate, display XZ;
!      drop Z coordinate, display XY;
!      orthographic projection into a plane;
!      perspective projection into a plane through a focus point;
!      project X into YZ using an angle THETA;
!      project Y into XZ using an angle THETA;
!      project Z into XY using an angle THETA.
!
!  Modified:
!
!    14 October 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real COR2(2,COR2_MAX), the projected 2D data.
!
!    Input, real COR3(3,COR3_MAX), the data to project.
!
!    Input, integer COR2_MAX, COR3_MAX, the maximum number of 2D and
!    3D points that can be handled.
!
!    Output, integer COR2_NUM, the number of 2D points.
!
!    Input, integer COR3_NUM, the number of 3D points.
!
  integer cor2_max
  integer cor3_max
!
  real cor2(2,cor2_max)
  integer cor2_num
  real cor3(3,cor3_max)
  integer cor3_num
  integer i
  integer ierror
  integer ios
  character ( len = 20 ) isay
  logical s_eqi
  real theta
  real x1
  real x2
  real x3
  real xf
  real y1
  real y2
  real y3
  real yf
  real z1
  real z2
  real z3
  real zf
!
  ierror = 0

  if ( cor3_num <= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'PROJECT_2D - Fatal error!'
    write ( *, * ) '  Input COR3_NUM <= 0.'
    ierror = 1
    return
  end if

  cor2_num = cor3_num

  write ( *, * ) ' '
  write ( *, * ) 'Choose a projection from 3D -> 2D:'
  write ( *, * ) ' '
  write ( *, * ) '-X     drop X coordinate, display YZ;'
  write ( *, * ) '-Y     drop Y coordinate, display XZ;'
  write ( *, * ) '-Z     drop Z coordinate, display XY;'
  write ( *, * ) 'OPLANE orthographic projection into plane.'
  write ( *, * ) 'PPLANE perspective projection into plane.'
  write ( *, * ) 'PX     project X into YZ using THETA;'
  write ( *, * ) 'PY     project Y into XZ using THETA;'
  write ( *, * ) 'PZ     project Z into XY using THETA;'
 
  read ( *, '(a)', iostat = ios ) isay

  if ( ios /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'PROJECT_2D - Unexpected end of input.'
    ierror = ios
    return
  end if
 
  if ( s_eqi ( isay, '-X' ) ) then
 
    do i = 1, cor3_num
      cor2(1,i) = cor3(2,i)
      cor2(2,i) = cor3(3,i)
    end do
 
  else if ( s_eqi ( isay, '-Y' ) ) then
 
    do i = 1, cor3_num
      cor2(1,i) = cor3(1,i)
      cor2(2,i) = cor3(3,i)
    end do
 
  else if ( s_eqi ( isay, '-Z' ) ) then
 
    do i = 1, cor3_num
      cor2(1,i) = cor3(1,i)
      cor2(2,i) = cor3(2,i)
    end do
 
  else if ( s_eqi ( isay, 'OPLANE' ) ) then
 
    write ( *, * ) ' '
    write ( *, * ) 'Enter 3 (X,Y,Z) points on the plane:'
    read ( *, *, iostat = ios ) x1, y1, z1, x2, y2, z2, x3, y3, z3

    if ( ios /= 0 ) then
      write ( *, * ) ' '
      write ( *, * ) 'PROJECT_2D - Unexpected end of input.'
      ierror = ios
      return
    end if
 
    call project_oplane ( x1, y1, z1, x2, y2, z2, x3, y3, z3, &
      cor2, cor3, cor2_max, cor3_max, cor3_num )

  else if ( s_eqi ( isay, 'PPLANE' ) ) then
 
    write ( *, * ) ' '
    write ( *, * ) 'Enter 3 (X,Y,Z) points on the plane:'
    read ( *, *, iostat = ios ) x1, y1, z1, x2, y2, z2, x3, y3, z3

    if ( ios /= 0 ) then
      write ( *, * ) ' '
      write ( *, * ) 'PROJECT_2D - Unexpected end of input.'
      ierror = ios
      return
    end if

    write ( *, *) 'Enter focus point (X,Y,Z):'
    read ( *, *, iostat = ios ) xf, yf, zf

    if ( ios /= 0 ) then
      write ( *, * ) ' '
      write ( *, * ) 'PROJECT_2D - Unexpected end of input.'
      ierror = ios
      return
    end if
 
    call project_pplane ( x1, y1, z1, x2, y2, z2, x3, y3, z3, xf, yf, zf, &
      cor2, cor3, cor2_max, cor3_max, cor2_num, cor3_num )
 
  else if ( s_eqi ( isay, 'PX' ) ) then
 
    write ( *, * ) ' '
    write ( *, * ) 'Enter projection angle THETA in degrees:'
    read ( *, *, iostat = ios ) theta

    if ( ios /= 0 ) then
      write ( *, * ) ' '
      write ( *, * ) 'PROJECT_2D - Unexpected end of input.'
      ierror = ios
      return
    end if

    call project_angle ( 'X', cor2, cor3, cor2_max, cor3_max, cor2_num, theta )

  else if ( s_eqi ( isay, 'PY' ) ) then
 
    write ( *, * ) ' '
    write ( *, * ) 'Enter projection angle THETA in degrees:'
    read ( *, *, iostat = ios ) theta

    if ( ios /= 0 ) then
      write ( *, * ) ' '
      write ( *, * ) 'PROJECT_2D - Unexpected end of input.'
      ierror = ios
      return
    end if

    call project_angle ( 'Y', cor2, cor3, cor2_max, cor3_max, cor2_num, theta )

  else if ( s_eqi ( isay, 'PZ' ) ) then
 
    write ( *, * ) ' '
    write ( *, * ) 'Enter projection angle THETA in degrees:'
    read ( *, *, iostat = ios ) theta

    if ( ios /= 0 ) then
      write ( *, * ) ' '
      write ( *, * ) 'PROJECT_2D - Unexpected end of input.'
      ierror = ios
      return
    end if

    call project_angle ( 'Z', cor2, cor3, cor2_max, cor3_max, cor2_num, theta )

  else
 
    write ( *, * ) ' '
    write ( *, * ) 'PROJECT_2D - Error!'
    write ( *, * ) '  Unrecognized projection option!'
    cor2_num = 0
    ierror = 1
    
  end if
 
  return
end
subroutine project_angle ( cor, cor2, cor3, cor2_max, cor3_max, cor2_num, &
  theta )
!
!*******************************************************************************
!
!! PROJECT_ANGLE converts 3D data to 2D using a presentation angle.
!
!
!  Discussion:
!
!    A "presentation angle" THETA is used to project the 3D point
!    (X3D, Y3D, Z3D) to the 2D projection (XVAL,YVAL).
!
!    The formula used if COR = 'X' is
!
!      X2D = Y3D - sin ( THETA ) * X3D
!      Y2D = Z3D - sin ( THETA ) * X3D
!
!  Modified:
!
!    14 October 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character COR, the coordinate to be projected.
!    COR should be 'X', 'Y', or 'Z'.
!
!    Output, real COR2(2,COR2_MAX), the 2D projections.
!
!    Input, real COR3(3,COR3_MAX), the 3D points to be projected.
!
!    Input, integer COR2_MAX, COR3_MAX, the maximum number of 2D and
!    3D points allowed.
!
!    Input, integer COR2_NUM, the number of 2D values to be computed.
!
!    Input, real THETA, the presentation angle in degrees.
!
  integer cor2_max
  integer cor3_max
!
  character cor
  real cor2(2,cor2_max)
  integer cor2_num
  real cor3(3,cor3_max)
  integer i
  real pi
  logical s_eqi
  real stheta
  real theta
!
  stheta = sin ( pi() * theta / 180.0 )

  if ( cor2_num > cor2_max ) then
    write ( *, * ) ' '
    write ( *, * ) 'PROJECT_ANGLE - Fatal error!'
    write ( *, * ) '  COR2_NUM is greater than COR2_MAX.'
    stop
  end if

  if ( cor2_num > cor3_max ) then
    write ( *, * ) ' '
    write ( *, * ) 'PROJECT_ANGLE - Fatal error!'
    write ( *, * ) '  COR2_NUM is greater than COR3_MAX.'
    stop
  end if

  if ( s_eqi ( cor, 'X' ) ) then

    do i = 1, cor2_num
      cor2(1,i) = cor3(2,i) - stheta * cor3(1,i)
      cor2(2,i) = cor3(3,i) - stheta * cor3(1,i)
    end do

  else if ( s_eqi ( cor, 'Y' ) ) then

    do i = 1, cor2_num
      cor2(1,i) = cor3(1,i) - stheta * cor3(2,i)
      cor2(2,i) = cor3(3,i) - stheta * cor3(2,i)
    end do

  else if ( s_eqi ( cor, 'Z' ) ) then

    do i = 1, cor2_num
      cor2(1,i) = cor3(1,i) - stheta * cor3(3,i)
      cor2(2,i) = cor3(2,i) - stheta * cor3(3,i)
    end do

  else

    write ( *, * ) ' '
    write ( *, * ) 'PROJECT_ANGLE - Fatal error!'
    write ( *, * ) '  Unrecognized axis.'
    stop

  end if

  return
end
subroutine project_oplane ( x1, y1, z1, x2, y2, z2, x3, y3, z3, cor2, cor3, &
  cor2_max, cor3_max, cor3_num )
!
!*******************************************************************************
!
!! PROJECT_OPLANE projects 3D points onto an orthographic plane.
!
!
!  Modified:
!
!    14 October 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3, are the 
!    coordinates of three points on the plane.
!
!    Output, real COR2(2,COR2_MAX), the "local" in-plane coordinates
!    of the projections of the object points.
!
!    Input, real COR3(3,COR3_MAX), the (X,Y,Z) coordinates of the object 
!     points.
!
!    Input, integer COR2_MAX, COR3_MAX, the maximum number of 2D and
!    3D points that can be handled.
!
!    Input, integer COR3_NUM, the number of points to project.
!
  integer cor2_max
  integer cor3_max
!
  real a
  real b
  real c
  real cor2(2,cor2_max)
  real cor3(3,cor3_max)
  integer cor3_num
  real d
  real dot
  integer i
  real v1(3)
  real v2(3)
  real x1
  real x2
  real x3
  real xn
  real y1
  real y2
  real y3
  real yn
  real z1
  real z2
  real z3
  real zn
!
!  Put the plane into ABCD form.
!
  call plane_exp2imp_3d ( x1, y1, z1, x2, y2, z2, x3, y3, z3, a, b, c, d )
!
!  Compute the two basis vectors for the affine plane.
!
  v1(1) = x2 - x1
  v1(2) = y2 - y1
  v1(3) = z2 - z1

  call vector_unit_nd ( 3, v1 )

  v2(1) = x3 - x1
  v2(2) = y3 - y1
  v2(3) = z3 - z1

  dot = v1(1) * v2(1) + v1(2) * v2(2) + v1(3) * v2(3)

  v2(1:3) = v2(1:3) - dot * v1(1:3)

  call vector_unit_nd ( 3, v2 )
!
!  For each point, its image in the plane is the nearest point 
!  in the plane.
!
  do i = 1, min ( cor3_num, cor3_max )
 
    call plane_imp_point_nearest_3d ( a, b, c, d, cor3(1,i), cor3(2,i), &
      cor3(3,i), xn, yn, zn )

    cor2(1,i) = ( xn - x1 ) * v1(1) + ( yn - y1 ) * v1(2) + ( zn - z1 ) * v1(3)
    cor2(2,i) = ( xn - x1 ) * v2(1) + ( yn - y1 ) * v2(2) + ( zn - z1 ) * v2(3)
 
  end do

  return
end
subroutine project_pplane ( x1, y1, z1, x2, y2, z2, x3, y3, z3, xf, yf, zf, &
  cor2, cor3, cor2_max, cor3_max, cor2_num, cor3_num )
!
!*******************************************************************************
!
!! PROJECT_PPLANE projects a point through a focus point onto a perspective plane.
!
!
!  Modified:
!
!    14 October 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3, are the 
!    coordinates of three points on the plane.
!
!    Input, real XF, YF, ZF, are the coordinates of the focus point.
!
!    Output, real COR2(2,COR2_MAX), the "local" in-plane coordinates
!    of the projections of the object points.
!
!    Input, real COR3(3,COR3_MAX), the (X,Y,Z) coordinates of points
!    to be projected.
!
!    Input, integer COR2_MAX, COR3_MAX, the maximum number of 2D and
!    3D points that can be handled.
!
!    Output, integer COR2_NUM, the number of projected points.
!
!    Input, integer COR3_NUM, the number of points to project.
!
  integer cor2_max
  integer cor3_max
!
  real a
  real alpha
  real angle_rad_3d
  real b
  real beta
  real c
  real cor2(2,cor2_max)
  integer cor2_num
  real cor3(3,cor3_max)
  integer cor3_num
  real d
  real disfo
  real disfn
  real dot
  integer i
  real v1(3)
  real v2(3)
  real x1
  real x2
  real x3
  real xf
  real xn
  real xp
  real y1
  real y2
  real y3
  real yf
  real yn
  real yp
  real z1
  real z2
  real z3
  real zf
  real zn
  real zp
!
!  Put the plane into ABCD form.
!
  call plane_exp2imp_3d ( x1, y1, z1, x2, y2, z2, x3, y3, z3, a, b, c, d )
!
!  Get the nearest point on the plane to the focus.
!
  call plane_imp_point_nearest_3d ( a, b, c, d, xf, yf, zf, xn, yn, zn )
!
!  Get the distance from the focus to the plane.
!
  call points_distance_3d ( disfn, xf, yf, zf, xn, yn, zn )
!
!  If the focus lies in the plane, this is bad.  We could still
!  project points that actually lie in the plane, but we'll
!  just bail out.
!
  if ( disfn == 0.0 ) then

    cor2_num = 0
    return

  end if
!
!  Compute the two basis vectors for the affine plane.
!
  v1(1) = x2 - x1
  v1(2) = y2 - y1
  v1(3) = z2 - z1

  call vector_unit_nd ( 3, v1 )

  v2(1) = x3 - x1
  v2(2) = y3 - y1
  v2(3) = z3 - z1

  dot = v1(1) * v2(1) + v1(2) * v2(2) + v1(3) * v2(3)

  v2(1:3) = v2(1:3) - dot * v1(1:3)

  call vector_unit_nd ( 3, v2 )
!
!  Process the points.
!
  do i = 1, min ( cor3_num, cor3_max )
!
!  Get the distance from the focus to the object.
!
    call points_distance_3d ( disfo, xf, yf, zf, cor3(1,i), cor3(2,i), &
      cor3(3,i) )
 
    if ( disfo == 0.0 ) then
 
      xp = xn
      yp = yn
      zp = zn 

    else
!
!  Compute ALPHA, the angle between (OBJECT-FOCUS) and (NEAREST-FOCUS).
!
      alpha = angle_rad_3d ( cor3(1,i), cor3(2,i), cor3(3,i), &
        xf, yf, zf, xn, yn, zn )
 
      if ( cos ( alpha ) == 0.0 ) then
 
        xp = xn
        yp = yn
        zp = zn
 
      else
!
!  Multiplier BETA is Dist(NEAREST-FOCUS) / ( Cos(ALPHA)*Dist(OBJECT-FOCUS) )
!
        beta = disfn / ( cos ( alpha ) * disfo )
!
!  Set the projected point.
!
        xp = xf + beta * ( cor3(1,i) - xf )
        yp = yf + beta * ( cor3(2,i) - yf )
        zp = zf + beta * ( cor3(3,i) - zf )
 
      end if
 
    end if
 
    cor2(1,i) = ( xp - x1 ) * v1(1) + ( yp - y1 ) * v1(2) + ( zp - z1 ) * v1(3)
    cor2(2,i) = ( xp - x1 ) * v2(1) + ( yp - y1 ) * v2(2) + ( zp - z1 ) * v2(3)

  end do
 
  return
end
subroutine ps_write ( cor2, face, face_material, face_order, fileout_name, iunit, &
  line_dex, line_material, material_rgba, cor2_max, face_max, line_max, material_max, &
  order_max, cor2_num, face_num, line_num )
!
!*******************************************************************************
!
!! PS_WRITE writes 2D face and line information to a PostScript file.
!
!
!  Discussion:
!
!    The intent is that a 3D model will be projected in some way
!    to a 2D model that can be printed out as a standard PostScript object.
!
!  Modified:
!
!    06 October 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real COR2(2,COR2_MAX), the X and Y components of 2D points.
!
!    Input, integer FACE(ORDER_MAX,FACE_MAX), the nodes making faces.
!
!    Input, integer FACE_MATERIAL(FACE_MAX), the material of each face.
!
!    Input, integer FACE_ORDER(FACE_MAX), the number of vertices per face.
!
!    Input, character ( len = 80 ) FILEOUT_NAME, the name of the output file.
!
!    Input, integer IUNIT, the FORTRAN unit to which data is written.
!
!    Input, integer LINE_DEX(LINE_MAX), nodes forming a line, terminated by -1.
!
!    Input, integer LINE_MAT(LINE_MAX), material index for line.
!
!    Input, real MATERIAL_RGBA(4,MATERIAL_MAX), material R, G, B and A values.
!
!    Input, integer COR2_MAX, the maximum number of 2D points.
!
!    Input, integer FACE_MAX, the maximum number of faces.
!
!    Input, integer LINE_MAX, the maximum number of line definition items.
!
!    Input, integer MATERIAL_MAX, the maximum number of materials.
!
!    Input, integer ORDER_MAX, the maximum number of vertices per face.
!
!    Input, integer COR2_NUM, the number of 2D points.
!
!    Input, integer FACE_NUM, the number of faces defined.
!
!    Input, integer LINE_NUM, the number of line definition items.
!
  integer cor2_max
  integer face_max
  integer line_max
  integer material_max
  integer order_max
!
  real alpha
  real blue
  real cor2(2,cor2_max)
  integer cor2_num
  integer face(order_max,face_max)
  integer face_material(face_max)
  integer face_num
  integer face_order(face_max)
  character ( len = 80 ) fileout_name
  real green
  integer i
  integer iface
  integer imat
  integer imat_old
  integer iunit
  integer j
  integer jhi
  integer k
  logical lineopen
  integer line_dex(line_max)
  integer line_material(line_max)
  integer line_num
  integer margin
  real material_rgba(4,material_max)
  integer page_num
  integer pagexmax
  integer pagexmin
  integer pageymax
  integer pageymin
  integer plotxmax
  integer plotxmin
  integer plotxmin2
  integer plotymax
  integer plotymin
  integer plotymin2
  integer px
  integer py
  real red
  integer text_num
  real xmax
  real xmin
  real ymax
  real ymin
!
!  Initialization
!
  imat_old = -999
  page_num = 1
  text_num = 0
!
!  Comput the bounding box.
!
  xmin = cor2(1,1)
  xmax = cor2(1,1)
  ymin = cor2(2,1)
  ymax = cor2(2,1)
  do i = 2, cor2_num
    xmin = min ( xmin, cor2(1,i) )
    xmax = max ( xmax, cor2(1,i) )
    ymin = min ( ymin, cor2(2,i) )
    ymax = max ( ymax, cor2(2,i) )
  end do

  if ( xmin == xmax ) then
    xmin = cor2(1,1) - 0.5
    xmax = cor2(1,1) + 0.5
  end if

  if ( ymin == ymax ) then
    ymin = cor2(2,1) - 0.5
    ymax = cor2(2,1) + 0.5
  end if
!
!  Compute the scale factor.
!
  pagexmax = 612
  pagexmin = 0
  pageymax = 792
  pageymin = 0

  margin = 36

  plotxmax = pagexmax - margin
  plotxmin = pagexmin + margin
  plotymax = pageymax - margin
  plotymin = pageymin + margin

  alpha = min ( ( plotxmax - plotxmin ) / ( xmax - xmin ), &
                ( plotymax - plotymin ) / ( ymax - ymin ) )
!
!  Adjust PLOTXMIN and PLOTYMIN to center the image.
!
  plotxmin2 = 0.5 * ( plotxmin + plotxmax - alpha * ( xmax - xmin ) )
  plotymin2 = 0.5 * ( plotymin + plotymax - alpha * ( ymax - ymin ) )
!
!  Prolog
!
  write ( iunit, '(a)' ) '%!PS-Adobe-3.0'
  write ( iunit, '(a)' ) '%%Document-Fonts: Times-Roman'
  write ( iunit, '(a)' ) '%%Title: ' // trim ( fileout_name )
  write ( iunit, '(a)' ) '%%Creator: IVREAD/PS_WRITE'
  write ( iunit, '(a)' ) '%%CreationDate: Tue Apr 4 15:57:00 1997'
  write ( iunit, '(a,4i5)' ) '%%BoundingBox', plotxmin, plotymin, plotxmax, &
    plotymax
  write ( iunit, '(a)' ) '%%LanguageLevel: 2'
  write ( iunit, '(a)' ) '%%EndComments'
  write ( iunit, '(a)' ) '%%BeginProlog'
  write ( iunit, '(a)' ) '%%EndProlog'

  text_num = text_num + 10
!
!  Fill the faces.
!
  red = 0.7
  green = 0.7
  blue = 0.0
  write ( iunit, '(3f7.4,a)' ) red, green, blue, ' setrgbcolor'

  text_num = text_num + 1

  do iface = 1, face_num

    imat = face_material(iface)

    if ( imat /= imat_old ) then
      imat_old = imat
      red = material_rgba(1,imat)
      green = material_rgba(2,imat)
      blue = material_rgba(3,imat)
      write ( iunit, '(3f7.4,a)' ) red, green, blue, ' setrgbcolor'
    end if

    jhi = face_order(iface)

    do j = 1, jhi + 1

      if ( j <= face_order(iface) ) then
        k = face(j,iface)
      else
        k = face(1,iface)
      end if

      px = plotxmin2 + nint ( alpha * ( cor2(1,k) - xmin ) )
      py = plotymin2 + nint ( alpha * ( cor2(2,k) - ymin ) )

      if ( j == 1 ) then
        write ( iunit, '(a)' ) ' newpath'
        write ( iunit, '(2i4,a)' ) px, py, ' moveto'
        text_num = text_num + 2
      else
        write ( iunit, '(2i4,a)' ) px, py, ' lineto'
        text_num = text_num + 1
      end if

    end do

    write ( iunit, '(a)' ) ' fill'
    text_num = text_num + 1

  end do
!
!  Draw the boundaries of the faces as black lines.
!
  red = 0.0
  green = 0.0
  blue = 0.0
  write ( iunit, '(3f7.4,a)' ) red, green, blue, ' setrgbcolor'
  text_num = text_num + 1

  do iface = 1, face_num

    jhi = face_order(iface)

    do j = 1, jhi + 1

      if ( j <= face_order(iface) ) then
        k = face(j,iface)
      else
        k = face(1,iface)
      end if

      px = plotxmin2 + nint ( alpha * ( cor2(1,k) - xmin ) )
      py = plotymin2 + nint ( alpha * ( cor2(2,k) - ymin ) )

      if ( j == 1 ) then
        write ( iunit, '(a)' ) ' newpath'
        write ( iunit, '(2i4,a)' ) px, py, ' moveto'
        text_num = text_num + 2
      else
        write ( iunit, '(2i4,a)' ) px, py, ' lineto'
        text_num = text_num + 1
      end if

    end do

    write ( iunit, '(a)' ) ' stroke'
    text_num = text_num + 1

  end do
!
!  Draw other lines.
!
!  We need to set the color of the lines as specified by the user using LINE_MAT.
!
  lineopen = .false.

  do i = 1, line_num

    j = line_dex(i)

    if ( j <= 0 ) then

      write ( iunit, '(a)' ) ' stroke'
      text_num = text_num + 1
      lineopen = .false.

    else

      imat = line_material(i)

      if ( imat /= imat_old ) then
        imat_old = imat
        red = material_rgba(1,imat)
        green = material_rgba(2,imat)
        blue = material_rgba(3,imat)
        write ( iunit, '(3f7.4,a)' ) red, green, blue, ' setrgbcolor'
      end if

      px = plotxmin2 + nint ( alpha * ( cor2(1,j) - xmin ) )
      py = plotymin2 + nint ( alpha * ( cor2(2,j) - ymin ) )

      if ( lineopen ) then
        write ( iunit, '(2i4,a)' ) px, py, ' lineto'
        text_num = text_num + 1
      else
        write ( iunit, '(a)' ) ' newpath'
        write ( iunit, '(2i4,a)' ) px, py, ' moveto'
        text_num = text_num + 2
        lineopen = .true.
      end if

    end if

  end do
!
!  Write the epilog.
!
  write ( iunit, '(a)' ) 'showpage'
  write ( iunit, '(a)' ) 'grestore'
  write ( iunit, '(a)' ) '%%Trailer'
  write ( iunit, '(a,i2)' ) '%%Pages: ', page_num
  text_num = text_num + 4
!
!  Report.
!
  write ( *, * ) ' '
  write ( *, * ) 'PS_WRITE - Wrote ', text_num, ' text lines.'

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
subroutine rgb_to_hue ( r, g, b, h )
!
!*******************************************************************************
!
!! RGB_TO_HUE converts (R,G,B) colors to a hue value between 0 and 1.
!
!
!  Discussion:
!
!    The hue computed here should be the same as the H value computed
!    for HLS and HSV, except that it ranges from 0 to 1 instead of
!    0 to 360.
!
!    A monochromatic color ( white, black, or a shade of gray) does not
!    have a hue.  This routine will return a special value of H = -1
!    for such cases.
!
!  Example:
!
!    Color    R    G    B     H
!
!    red      1.0  0.0  0.0   0.00
!    yellow   1.0  1.0  0.0   0.16
!    green    0.0  1.0  0.0   0.33
!    cyan     0.0  1.0  1.0   0.50
!    blue     0.0  0.0  1.0   0.67
!    magenta  1.0  0.0  1.0   0.83
!
!    black    0.0  0.0  0.0  -1.00
!    gray     0.5  0.5  0.5  -1.00
!    white    1.0  1.0  1.0  -1.00
!
!  Modified:
!
!    25 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real R, G, B, the red, green and blue values of the color.
!    These values should be between 0 and 1.
!
!    Output, real H, the corresponding hue of the color, or -1.0 if
!    the color is monochromatic.
!
  real b
  real b2
  real g
  real g2
  real h
  real r
  real r2
  real rgbmax
  real rgbmin
!
!  Make sure the colors are between 0 and 1.
!
  r2 = min ( max ( r, 0.0 ), 1.0 )
  g2 = min ( max ( g, 0.0 ), 1.0 )
  b2 = min ( max ( b, 0.0 ), 1.0 )
!
!  Compute the minimum and maximum of R, G and B.
!
  rgbmax = r2
  rgbmax = max ( rgbmax, g2 )
  rgbmax = max ( rgbmax, b2 )

  rgbmin = r2
  rgbmin = min ( rgbmin, g2 )
  rgbmin = min ( rgbmin, b2 )
!
!  If RGBMAX = RGBMIN, then the color has no hue.
!
  if ( rgbmax == rgbmin ) then

    h = - 1.0
!
!  Otherwise, we need to determine the dominant color.
!
  else

    if ( r2 == rgbmax ) then
      h = ( g2 - b2 ) / ( rgbmax - rgbmin )
    else if ( g2 == rgbmax ) then
      h = 2.0 + ( b2 - r2 ) / ( rgbmax - rgbmin )
    else if ( b2 == rgbmax ) then
      h = 4.0 + ( r2 - g2 ) / ( rgbmax - rgbmin )
    end if

    h = h / 6.0
!
!  Make sure H lies between 0 and 1.0.
!
    if ( h < 0.0 ) then
      h = h + 1.0
    else if ( h > 1.0 ) then
      h = h - 1.0
    end if

  end if

  return
end
subroutine relnex ( line, rval, done )
!
!*******************************************************************************
!
!! RELNEX "reads" real numbers from a string, one at a time.
!
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
!    Input, character ( len = * ) LINE, a string, presumably containing real
!    numbers.  These may be separated by spaces or commas.
!
!    Output, real RVAL.  If DONE is FALSE, then RVAL contains the
!    "next" real value read from LINE.  If DONE is TRUE, then
!    RVAL is zero.
!
!    Input/output, logical DONE.
!    On input with a fresh value of LINE, the user should set
!    DONE to TRUE.
!    On output, the routine sets DONE to FALSE if another real
!    value was read, or TRUE if no more reals could be read.
!
  logical done
  integer ierror
  integer lchar
  character ( len = * ) line
  integer, save :: next = 1
  real rval
!
  rval = 0.0

  if ( done ) then
    next = 1
    done = .false.
  end if

  if ( next > len_trim ( line ) ) then
    done = .true.
    return
  end if

  call s_to_r ( line(next:), rval, ierror, lchar )

  if ( ierror /= 0 .or. lchar == 0 ) then
    done = .true.
    next = 1
  else
    done = .false.
    next = next + lchar
  end if

  return
end
subroutine rvec_to_s ( n, x, s )
!
!*******************************************************************************
!
!! RVEC_TO_S "writes" a real vector into a string.
!
!
!  Discussion:
!
!    The values will be separated by commas and a single space.
!    If the string is too short, then data will be lost.
!
!  Modified:
!
!    30 June 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the dimension of X.
!
!    Input, real X(N), a vector to be written to a string.
!
!    Output, character ( len = * ) S, a string to which the real vector
!    has been written.
!
  integer n
!
  integer i
  character ( len = * ) s
  character ( len = 14 ) s2
  real x(n)
!
  do i = 1, n

    if ( x(i) == 0.0 ) then
      s2 = '0'
    else if ( abs ( x(i) ) >= 1.0E+10 ) then
      write ( s2, '(g14.6)' ) x(i)
      call s_trim_zeros ( s2 )
    else if ( real ( int ( x(i) ) ) == x(i) ) then
      write ( s2, '(i12)' ) int ( x(i) )
    else
      write ( s2, '(g14.6)' ) x(i)
      call s_trim_zeros ( s2 )
    end if

    if ( i == 1 ) then
      s = adjustl ( s2 )
    else
      s = trim ( s ) // ', ' // adjustl ( s2 )
    end if

  end do

  return
end
subroutine s_blank_delete ( s )
!
!*******************************************************************************
!
!! S_BLANK_DELETE removes blanks from a string, left justifying the remainder.
!
!
!  Comment:
!
!    All TAB characters are also removed.
!
!  Modified:
!
!    26 July 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character ( len = * ) S, the string to be transformed.
!
  character c
  integer iget
  integer iput
  integer nchar
  character ( len = * ) s
  character, parameter :: TAB = char ( 9 )
!
  iput = 0
  nchar = len_trim ( s )

  do iget = 1, nchar

    c = s(iget:iget)

    if ( c /= ' ' .and. c /= TAB ) then
      iput = iput + 1
      s(iput:iput) = c
    end if

  end do

  s(iput+1:) = ' '

  return
end
subroutine s_blanks_delete ( s )
!
!*******************************************************************************
!
!! S_BLANKS_DELETE replaces consecutive blanks by one blank.
!
!
!  Discussion:
!
!    The remaining characters are left justified and right padded with blanks.
!    TAB characters are converted to spaces.
!
!  Modified:
!
!    26 July 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character ( len = * ) S, the string to be transformed.
!
  integer i
  integer j
  integer nchar
  character newchr
  character oldchr
  character ( len = * ) s
  character, parameter :: TAB = char ( 9 )
!
  j = 0
  newchr = ' '
  nchar = len_trim ( s )

  do i = 1, nchar

    oldchr = newchr
    newchr = s(i:i)

    if ( newchr == TAB ) then
      newchr = ' '
    end if

    s(i:i) = ' '

    if ( oldchr /= ' ' .or. newchr /= ' ' ) then
      j = j + 1
      s(j:j) = newchr
    end if

  end do

  return
end
subroutine s_cap ( string )
!
!*******************************************************************************
!
!! S_CAP replaces any lowercase letters by uppercase ones in a string.
!
!
!  Modified:
!
!    16 May 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character ( len = * ) STRING, the string to be transformed.
!
  character c
  integer i
  integer nchar
  character ( len = * ) string
!
  nchar = len_trim ( string )

  do i = 1, nchar

    c = string(i:i)
    call c_cap ( c )
    string(i:i) = c

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
function s_is_i ( string, ival )
!
!*******************************************************************************
!
!! S_IS_I returns .TRUE. if STRING represents an integer.
!
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
!    Input, character ( len = * ) STRING, the string to be checked.
!
!    Output, integer IVAL.  If S_IS_INT is TRUE, then IVAL is the
!    integer represented.  Otherwise IVAL is 0.
!
!    Output, logical S_IS_I, .TRUE. if STRING represents an integer.
!
  integer ierror
  integer ival
  integer lchar
  integer lenc
  logical s_is_i
  character ( len = * ) string
!
  lenc = len_trim ( string )

  call s_to_i ( string, ival, ierror, lchar )

  if ( ierror == 0 .and. lchar >= lenc ) then
    s_is_i = .true.
  else
    s_is_i = .false.
    ival = 0
  end if

  return
end
subroutine s_is_r ( string, rval, lval )
!
!*******************************************************************************
!
!! S_IS_R returns .TRUE. if STRING represents a real number.
!
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
!    Input, character ( len = * ) STRING, the string to be checked.
!
!    Output, real RVAL.  If ISREAL is TRUE, then RVAL is the real
!    number represented.  Otherwise RVAL is 0.
!
!    Output, logical LVAL, .TRUE. if STRING represents a real number.
!
  integer ierror
  integer lchar
  logical lval
  real rval
  character ( len = * ) string
!
  call s_to_r ( string, rval, ierror, lchar )

  if ( ierror == 0 .and. lchar >= len_trim ( string ) ) then
    lval = .true.
  else
    lval = .false.
    rval = 0.0
  end if

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
subroutine s_trim_zeros ( s )
!
!*******************************************************************************
!
!! S_TRIM_ZEROS removes trailing zeros from a string.
!
!
!  Example:
!
!    Input:
!
!      S = '1401.072500'
!
!    Output:
!
!      S = '1401.0725'
!
!  Modified:
!
!    30 June 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character ( len = * ) S, the string to be operated on.
!
  integer i
  character ( len = * ) s
!
  i = len_trim ( s )

  do while ( i > 0 .and. s(i:i) == '0' )
    s(i:i) = ' '
    i = i - 1
  end do

  return
end
subroutine smf_read ( cor3, cor3_material, cor3_normal, cor3_tex_uv, debug, face, &
  face_material, face_normal, face_order, face_tex_uv, ierror, iunit, &
  material_name, material_rgba, cor3_max, face_max, material_max, order_max, &
  texture_max, bad_num, cor3_num, face_num, group_num, material_num, texture_num, &
  text_num, texture_name, vertex_material )
!
!*******************************************************************************
!
!! SMF_READ reads graphics information from an SMF file.
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
!    #SMF2.0
!    #  cube_face.smf
!    #  This example demonstrates how an RGB color can be assigned to
!    #  each face of an object.
!    #    
!    # First, define the geometry of the cube.
!    #
!    v 0.0  0.0  0.0
!    v 1.0  0.0  0.0
!    v 0.0  1.0  0.0
!    v 1.0  1.0  0.0
!    v 0.0  0.0  1.0
!    v 1.0  0.0  1.0
!    v 0.0  1.0  1.0
!    v 1.0  1.0  1.0
!    f 1 4 2
!    f 1 3 4
!    f 5 6 8
!    f 5 8 7
!    f 1 2 6
!    f 1 6 5
!    f 2 4 8
!    f 2 8 6
!    f 4 3 7
!    f 4 7 8
!    f 3 1 5
!    f 3 5 7
!    #
!    #  Colors will be bound 1 per face.
!    #
!    bind c face
!    c 1.0  0.0  0.0
!    c 1.0  0.0  0.0
!    c 0.0  1.0  0.0
!    c 0.0  1.0  0.0
!    c 0.0  0.0  1.0
!    c 0.0  0.0  1.0
!    c 1.0  1.0  0.0
!    c 1.0  1.0  0.0
!    c 0.0  1.0  1.0
!    c 0.0  1.0  1.0
!    c 1.0  0.0  1.0
!    c 1.0  0.0  1.0
!    #
!    #  Normal vectors will be bound 1 per face.
!    #
!    bind n face
!    n  0.0   0.0  -1.0
!    n  0.0   0.0  -1.0
!    n  0.0   0.0   1.0
!    n  0.0   0.0   1.0
!    n  0.0  -1.0   0.0
!    n  0.0  -1.0   0.0
!    n  1.0   0.0   0.0
!    n  1.0   0.0   0.0
!    n  0.0   1.0   0.0
!    n  0.0   1.0   0.0
!    n -1.0   0.0   0.0
!    n -1.0   0.0   0.0
!    #
!    #  Texture coordinate pairs will be bound 1 per face.
!    #
!    bind r face
!    r  0.0   0.0
!    r  0.0   0.1
!    r  0.0   0.2
!    r  0.0   0.3
!    r  0.1   0.0
!    r  0.1   0.1
!    r  0.1   0.2
!    r  0.1   0.3
!    r  0.2   0.0
!    r  0.2   0.1
!    r  0.2   0.2
!    r  0.2   0.3
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
!    Input/output, real COR3(3,COR3_MAX), the coordinates of points.
!
!    Input/output, integer COR3_MATERIAL(COR3_MAX), the material index of each node.
!
!    Input/output, real COR3_TEX_UV(2,COR3_MAX), UV texture coordinates for nodes.
!
!    Input, logical DEBUG, debugging switch.
!
!    Input/output, integer FACE(ORDER_MAX,FACE_MAX), the nodes making faces.
!
!    Input/output, integer FACE_MATERIAL(FACE_MAX), the material of each face.
!
!    Input/output, real FACE_NORMAL(3,FACE_MAX), the normal vector at each face.
!
!    Input/output, integer FACE_ORDER(FACE_MAX), the number of 
!    vertices per face.
!
!    Output, integer IERROR, an error flag.
!
!    Input, integer IUNIT, the FORTRAN unit from which data is read.
!
!    Input/output, character ( len = 100 ) MATERIAL_NAME(MATERIAL_MAX), material names.
!
!    Input/output, real MATERIAL_RGBA(4,MATERIAL_MAX), material R, G, B and A values.
!
!    Input, integer COR3_MAX, the maximum number of points.
!
!    Input, integer FACE_MAX, the maximum number of faces.
!
!    Input, integer MATERIAL_MAX, the maximum number of materials.
!
!    Input, integer ORDER_MAX, the maximum number of vertices per face.
!
!    Output, integer BAD_NUM, the number of bad lines of text in the file.
!
!    Input/output, integer COR3_NUM, the number of points.
!
!    Input/output, integer FACE_NUM, the number of faces.
!
!    Input/output, integer GROUP_NUM, the number of groups.
!
!    Input/output, integer MATERIAL_NUM, the number of materials.
!
!    Output, integer TEXT_NUM, the number of lines of text in the file.
!
!    Input/output, character ( len = 100 ) TEXTURE_NAME(TEXTURE_MAX), texture names.
!
!    Output, real TRANSFORM_MATRIX(4,4), the transformation matrix.
!
!    Input/output, integer VERTEX_MAT(ORDER_MAX,FACE_MAX), vertex materials.
!
  integer cor3_max
  integer face_max
  integer material_max
  integer order_max
  integer texture_max
!
  real angle
  character axis
  real b
  integer bad_num
  character ( len = 4 ) char4
  character cnr
  real cor3(3,cor3_max)
  integer cor3_material(cor3_max)
  real cor3_normal(3,cor3_max)
  integer cor3_num
  real cor3_tex_uv(2,cor3_max)
  logical debug
  logical done
  real dx
  real dy
  integer face(order_max,face_max)
  integer face_count
  integer face_material(face_max)
  real face_normal(3,face_max)
  integer face_num
  integer face_order(face_max)
  real face_tex_uv(2,face_max)
  real g
  integer group_num
  integer i
  integer icor3_normal
  integer icor3_tex_uv
  integer ierror
  integer iface_normal
  integer iface_tex_uv
  integer imat
  integer itemp
  integer iunit
  integer ivert
  integer iword
  integer k
  integer lchar
  integer level
  character ( len = 256 ) line
  character ( len = 30 ) material_binding
  character ( len = 100 ) material_name(material_max)
  integer material_num
  real material_rgba(4,material_max)
  integer node_count
  character ( len = 30 ) normal_binding
  real r
  real rgba(4)
  logical s_eqi
  real sx
  real sy
  real sz
  real temp
  integer text_num
  character ( len = 30 ) texture_binding
  character ( len = 100 ) texture_name(texture_max)
  integer texture_num
  real transform_matrix(4,4)
  character ( len = 30 ) type
  real u
  real v
  integer vertex_base
  integer vertex_correction
  integer vertex_material(order_max,face_max)
  character ( len = 256 ) word
  character ( len = 256 ) word1
  real x
  real xvec(3)
  real y
  real z
!
  face_count = 0
  ierror = 0
  icor3_normal = 0
  icor3_tex_uv = 0
  iface_normal = 0
  iface_tex_uv = 0
  level = 0
  material_binding = 'UNDEFINED'
  normal_binding = 'UNDEFINED'
  node_count = 0
  texture_binding = 'UNDEFINED'
  vertex_base = cor3_num
  vertex_correction = 0
  word = ' '

  call tmat_init ( transform_matrix )
!
!  Read a line of text from the file.
!
10    continue
 
  read ( iunit, '(a)', end = 30 ) line

  if ( debug ) then
    write ( *, '(a)' ) trim ( line )
  end if

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

  if ( debug ) then
    write ( *, '(a)' ) trim ( word )
  end if
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
!  BEGIN
!  Reset the transformation matrix to identity.
!  Node numbering starts at zero again.  (Really, this is level based)
!  (Really should define a new transformation matrix, and concatenate.)
!  (Also, might need to keep track of level.)
!
  if ( s_eqi ( word1, 'BEGIN' ) ) then

    level = level + 1

    vertex_base = cor3_num
    group_num = group_num + 1
    call tmat_init ( transform_matrix )
!
!  BIND [c|n|r] [vertex|face]
!  Specify the binding for RGB color, Normal, or Texture.
!  Options are "vertex" or "face"
!
  else if ( s_eqi ( word1, 'BIND' ) ) then

    call word_nexrd ( line, cnr, done )

    call word_nexrd ( line, type, done )

    if ( s_eqi ( cnr, 'C' ) ) then

      if ( s_eqi ( type, 'VERTEX' ) ) then
        material_binding = 'PER_VERTEX'
      else if ( s_eqi ( type, 'FACE' ) ) then
        material_binding = 'PER_FACE'
      end if

    else if ( s_eqi ( cnr, 'N' ) ) then

      if ( s_eqi ( type, 'VERTEX' ) ) then
        normal_binding = 'PER_VERTEX'
      else if ( s_eqi ( type, 'FACE' ) ) then
        normal_binding = 'PER_FACE'
      end if

    else if ( s_eqi ( cnr, 'R' ) ) then

      if ( s_eqi ( type, 'VERTEX' ) ) then
        texture_binding = 'PER_VERTEX'
      else if ( s_eqi ( type, 'FACE' ) ) then
        texture_binding = 'PER_FACE'
      end if

    end if
!
!  C <r> <g> <b>
!  Specify an RGB color, with R, G, B between 0.0 and 1.0.
!
  else if ( s_eqi ( word1, 'C' ) ) then

    call word_nexrd ( line, word, done )
    call s_to_r ( word, r, ierror, lchar )
    call word_nexrd ( line, word, done )
    call s_to_r ( word, g, ierror, lchar )
    call word_nexrd ( line, word, done )
    call s_to_r ( word, b, ierror, lchar )
!
!    Set up a temporary material (R,G,B,1.0).
!    Add the material to the material database, or find the index of
!      a matching material already in.
!    Assign the material of the node or face to this index.
!
    rgba(1) = r
    rgba(2) = g
    rgba(3) = b
    rgba(4) = 1.0

    if ( material_num <= 1000 ) then
      call rcol_find ( 4, material_num, material_rgba, rgba, imat )
    else
      imat = 0
    end if

    if ( imat == 0 ) then

      material_num = material_num + 1

      if ( material_num <= material_max ) then

        call i_to_s_zero ( material_num, char4 )

        material_name(material_num) = 'Material_' // char4
        material_rgba(1:4,material_num) = rgba(1:4)
        imat = material_num

      else

        imat = 0

      end if

    end if

    if ( material_binding == 'PER_FACE' ) then

      face_count = face_count + 1
      face_material(face_count) = imat

    else if ( material_binding == 'PER_VERTEX' ) then

      node_count = node_count + 1
      cor3_material(node_count) = imat

    else

      write ( *, * ) ' '
      write ( *, * ) 'SMF_READ - Fatal error!'
      write ( *, * ) '  Material binding undefined!'
      stop

    end if
!
!  END
!  Drop down a level.
!
  else if ( s_eqi ( word1, 'END' ) ) then

    level = level - 1

    if ( level < 0 ) then
      write ( *, * ) ' '
      write ( *, * ) 'SMF_READ - Fatal error!'
      write ( *, * ) '  More END statements than BEGINs.'
      write ( *, * ) '  Stopping on line ', text_num
      stop
    end if
!
!  F V1 V2 V3 ...
!  Face.
!  A face is defined by the vertices.
!
  else if ( s_eqi ( word1, 'F' ) ) then

    face_num = face_num + 1
    face_material(face_num) = material_num

    ivert = 0

15      continue

    ivert = ivert + 1
 
    call word_nexrd ( line, word, done )

    if ( done ) then
      go to 10
    end if
!
!  Read the vertex index.
!  Note that vertex indices start back at 0 each time a BEGIN is entered.
!  The strategy here won't handle nested BEGIN's, just one at a time.
!
    call s_to_i ( word, itemp, ierror, lchar )

    if ( ierror /= 0 ) then
      itemp = -1
      ierror = 0
      write ( *, * ) 'SMF_READ - Error!'
      write ( *, * ) '  Bad FACE field.'
      write ( *, '(a)' ) trim ( word )
    end if

    if ( ivert <= order_max .and. face_num <= face_max ) then
      face(ivert,face_num) = itemp + vertex_base
      vertex_material(ivert,face_num) = material_num
    end if

    if ( face_num <= face_max ) then
      face_order(face_num) = ivert
    end if

    go to 15
!
!  N <x> <y> <z>
!  Specify a normal vector.
!
  else if ( s_eqi ( word1, 'N' ) ) then

    call word_nexrd ( line, word, done )
    call s_to_r ( word, x, ierror, lchar )
    call word_nexrd ( line, word, done )
    call s_to_r ( word, y, ierror, lchar )
    call word_nexrd ( line, word, done )
    call s_to_r ( word, z, ierror, lchar )

    if ( normal_binding == 'PER_FACE' ) then

      iface_normal = iface_normal + 1

      face_normal(1,iface_normal) = x
      face_normal(2,iface_normal) = y
      face_normal(3,iface_normal) = z

    else if ( normal_binding == 'PER_VERTEX' ) then

      icor3_normal = icor3_normal + 1

      cor3_normal(1,icor3_normal) = x
      cor3_normal(2,icor3_normal) = y
      cor3_normal(3,icor3_normal) = z

    else

      write ( *, * ) ' '
      write ( *, * ) 'SMF_READ - Fatal error!'
      write ( *, * ) '  Normal binding undefined!'
      stop
                
    end if
!
!  R <u> <v>
!  Specify a texture coordinate.
!  (Read, but ignore for now.)
!
  else if ( s_eqi ( word1, 'R' ) ) then

    call word_nexrd ( line, word, done )
    call s_to_r ( word, u, ierror, lchar )
    call word_nexrd ( line, word, done )
    call s_to_r ( word, v, ierror, lchar )

    if ( texture_binding == 'PER_FACE' ) then

      iface_tex_uv = iface_tex_uv + 1

      face_tex_uv(1,iface_tex_uv) = u
      face_tex_uv(2,iface_tex_uv) = v

    else if ( texture_binding == 'PER_VERTEX' ) then

      icor3_tex_uv = icor3_tex_uv + 1

      cor3_tex_uv(1,icor3_tex_uv) = u
      cor3_tex_uv(2,icor3_tex_uv) = v

    else
      write ( *, * ) ' '
      write ( *, * ) 'SMF_READ - Fatal error!'
      write ( *, * ) '  Texture binding undefined!'
      stop
    end if
!
!  ROT [x|y|z] <theta>
!  (Read, but ignore for now.)
!
  else if ( s_eqi ( word1, 'ROT' ) ) then

    call word_nexrd ( line, axis, done )

    call word_nexrd ( line, word, done )
    call s_to_r ( word, angle, ierror, lchar )

    call tmat_rot_axis ( transform_matrix, transform_matrix, angle, axis )
!
!  SCALE <sx> <sy> <sz>
!  (Read, but ignore for now.)
!
  else if ( s_eqi ( word1, 'SCALE' ) ) then

    call word_nexrd ( line, word, done )
    call s_to_r ( word, sx, ierror, lchar )
    call word_nexrd ( line, word, done )
    call s_to_r ( word, sy, ierror, lchar )
    call word_nexrd ( line, word, done )
    call s_to_r ( word, sz, ierror, lchar )

    call tmat_scale ( transform_matrix, transform_matrix, sx, sy, sz )
!
!  SET VERTEX_CORRECTION <i>
!  Specify increment to add to vertex indices in file.
!  (Read, but ignore for now.)
!
  else if ( s_eqi ( word1, 'SET' ) ) then

    call word_nexrd ( line, word, done )

    call word_nexrd ( line, word, done )
    call s_to_i ( word, vertex_correction, ierror, lchar )
!
!  T_SCALE <dx> <dy>
!  Specify a translation to texture coordinates.
!  (Read, but ignore for now.)
!
  else if ( s_eqi ( word1, 'T_SCALE' ) ) then

    call word_nexrd ( line, word, done )
    call s_to_r ( word, dx, ierror, lchar )
    call word_nexrd ( line, word, done )
    call s_to_r ( word, dy, ierror, lchar )
!
!  T_TRANS <dx> <dy>
!  Specify a translation to texture coordinates.
!  (Read, but ignore for now.)
!
  else if ( s_eqi ( word1, 'T_TRANS' ) ) then

    call word_nexrd ( line, word, done )
    call s_to_r ( word, dx, ierror, lchar )
    call word_nexrd ( line, word, done )
    call s_to_r ( word, dy, ierror, lchar )
!
!  TEX <filename>
!  Specify a filename containing the texture.
!
  else if ( s_eqi ( word1, 'TEX' ) ) then

    call word_nexrd ( line, word, done )

    texture_num = texture_num + 1
    texture_name(texture_num) = word
!
!  TRANS <dx> <dy> <dz>
!
  else if ( s_eqi ( word1, 'TRANS' ) ) then

    call word_nexrd ( line, word, done )
    call s_to_r ( word, x, ierror, lchar )
    call word_nexrd ( line, word, done )
    call s_to_r ( word, y, ierror, lchar )
    call word_nexrd ( line, word, done )
    call s_to_r ( word, z, ierror, lchar )

    call tmat_trans ( transform_matrix, transform_matrix, x, y, z )
!
!  V X Y Z
!  Geometric vertex.
!
  else if ( s_eqi ( word1, 'V' ) ) then

    cor3_num = cor3_num + 1

    do i = 1, 3
      call word_nexrd ( line, word, done )
      call s_to_r ( word, temp, ierror, lchar )
      xvec(i) = temp
    end do
!
!  Apply current transformation matrix.
!  Right now, we can only handle one matrix, not a stack of
!  matrices representing nested BEGIN/END's.
!
    call tmat_mxp ( transform_matrix, xvec, xvec )

    if ( cor3_num <= cor3_max ) then
      cor3(1:3,cor3_num) = xvec(1:3)
    end if
!
!  Unrecognized keyword.
!
  else

    bad_num = bad_num + 1

    if ( bad_num <= 10 ) then
      write ( *, * ) ' '
      write ( *, * ) 'SMF_READ: Bad data on line ', text_num
      write ( *, '(a)' ) '  Bad word: ' // trim ( word )
    end if

  end if

  go to 10
!
!  End of information in file.
!
30    continue
!
!  Extend the material definition 
!  * from the face to the vertices and nodes, or
!  * from the vertices to the faces and nodes.
!
  if ( material_binding == 'PER_FACE' ) then

    do ivert = 1, order_max
      vertex_material(ivert,1:face_num) = face_material(1:face_num)
    end do

    call vertex_to_node_material ( cor3_material, face, face_order, &
      cor3_max, face_max, order_max, face_num, vertex_material )

  else if ( material_binding == 'PER_VERTEX' ) then

    call node_to_vertex_material ( cor3_material, face, face_order, &
      cor3_max, face_max, order_max, face_num, vertex_material )

    face_material(1:face_num) = vertex_material(1,1:face_num)

  end if

  return
end
subroutine smf_write ( cor3, cor3_material, cor3_normal, cor3_tex_uv, face, &
  face_order, filein_name, fileout_name, iunit, material_rgba, cor3_max, &
  face_max, material_max, order_max, texture_max, cor3_num, face_num, texture_num, &
  texture_name )
!
!*******************************************************************************
!
!! SMF_WRITE writes graphics information to an SMF file.
!
!
!  Example:
!
!    #SMF2.0
!    #  cube_face.smf
!    #  This example demonstrates how an RGB color can be assigned to
!    #  each face of an object.
!    #    
!    # First, define the geometry of the cube.
!    #
!    v 0.0  0.0  0.0
!    v 1.0  0.0  0.0
!    v 0.0  1.0  0.0
!    v 1.0  1.0  0.0
!    v 0.0  0.0  1.0
!    v 1.0  0.0  1.0
!    v 0.0  1.0  1.0
!    v 1.0  1.0  1.0
!    f 1 4 2
!    f 1 3 4
!    f 5 6 8
!    f 5 8 7
!    f 1 2 6
!    f 1 6 5
!    f 2 4 8
!    f 2 8 6
!    f 4 3 7
!    f 4 7 8
!    f 3 1 5
!    f 3 5 7
!    #
!    #  Colors will be bound 1 per face.
!    #
!    bind c face
!    c 1.0  0.0  0.0
!    c 1.0  0.0  0.0
!    c 0.0  1.0  0.0
!    c 0.0  1.0  0.0
!    c 0.0  0.0  1.0
!    c 0.0  0.0  1.0
!    c 1.0  1.0  0.0
!    c 1.0  1.0  0.0
!    c 0.0  1.0  1.0
!    c 0.0  1.0  1.0
!    c 1.0  0.0  1.0
!    c 1.0  0.0  1.0
!    #
!    #  Normal vectors will be bound 1 per face.
!    #
!    bind n face
!    n  0.0   0.0  -1.0
!    n  0.0   0.0  -1.0
!    n  0.0   0.0   1.0
!    n  0.0   0.0   1.0
!    n  0.0  -1.0   0.0
!    n  0.0  -1.0   0.0
!    n  1.0   0.0   0.0
!    n  1.0   0.0   0.0
!    n  0.0   1.0   0.0
!    n  0.0   1.0   0.0
!    n -1.0   0.0   0.0
!    n -1.0   0.0   0.0
!    #
!    #  Texture coordinate pairs will be bound 1 per face.
!    #
!    bind r face
!    r  0.0   0.0
!    r  0.0   0.1
!    r  0.0   0.2
!    r  0.0   0.3
!    r  0.1   0.0
!    r  0.1   0.1
!    r  0.1   0.2
!    r  0.1   0.3
!    r  0.2   0.0
!    r  0.2   0.1
!    r  0.2   0.2
!    r  0.2   0.3
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
!    Input, real COR3(3,COR3_MAX), the coordinates of points.
!
!    Input, integer COR3_MATERIAL(COR3_MAX), the material index of each node.
!
!    Input, real COR3_NORMAL(3,COR3_MAX), normals at nodes.  
!
!    Input, real COR3_TEX_UV(2,COR3_MAX), UV texture coordinates for nodes.
!
!    Input, integer FACE(ORDER_MAX,FACE_MAX), the nodes making faces.
!
!    Input, integer FACE_ORDER(FACE_MAX), the number of vertices per face.
!
!    Input, character ( len = 100 ) FILEIN_NAME, the name of the input file.
!
!    Input, character ( len = 100 ) FILEOUT_NAME, the name of the output file.
!
!    Input, integer IUNIT, the FORTRAN unit to which output is written.
!
!    Input, real MATERIAL_RGBA(4,MATERIAL_MAX), material R, G, B and A values.
!
!    Input, integer COR3_MAX, the maximum number of points.
!
!    Input, integer FACE_MAX, the maximum number of faces.
!
!    Input, integer MATERIAL_MAX, the maximum number of materials.
!
!    Input, integer ORDER_MAX, the maximum number of vertices per face.
!
!    Input, integer TEXTURE_MAX, the maximum number of textures.
!
!    Input, integer COR3_NUM, the number of points.
!
!    Input, integer FACE_NUM, the number of faces.
!
!    Input, integer TEXTURE_NUM, the number of textures.
!
!    Input, character ( len = 100 ) TEXTURE_NAME(TEXTURE_MAX), texture names.
!
  integer cor3_max
  integer face_max
  integer material_max
  integer order_max
  integer texture_max
!
  real b
  real cor3(3,cor3_max)
  integer cor3_material(cor3_max)
  real cor3_normal(3,cor3_max)
  integer cor3_num
  real cor3_tex_uv(2,cor3_max)
  integer face(order_max,face_max)
  integer face_num
  integer face_order(face_max)
  character ( len = 100 ) filein_name
  character ( len = 100 ) fileout_name
  real g
  integer i
  integer icor3
  integer iface
  integer imat
  integer iunit
  integer ivert
  real material_rgba(4,material_max)
  real r
  character ( len = 256 ) text
  integer text_num
  character ( len = 100 ) texture_name(texture_max)
  integer texture_num
!
  text_num = 0

  write ( iunit, '(a)' ) '#$SMF 1.0'
  write ( iunit, '(a, i8)' ) '#$vertices ', cor3_num
  write ( iunit, '(a, i8)' ) '#$faces ', face_num
  write ( iunit, '(a)' ) '#'
  write ( iunit, '(a)' ) '# ' // trim ( fileout_name ) // ' created by IVREAD.'
  write ( iunit, '(a)' ) '# Original data in ' // trim ( filein_name ) // '.'
  write ( iunit, '(a)' ) '#'
  text_num = text_num + 7
!
!  V: vertex coordinates.
!
  do icor3 = 1, cor3_num
    write ( text, '(a1,2x,3g14.6)' ) 'v', cor3(1:3,icor3)
    call s_blanks_delete ( text )
    write ( iunit, '(a)' ) trim ( text )
    text_num = text_num + 1
  end do
!
!  F: Faces.
!
  if ( face_num > 0 ) then
    write ( iunit, '(a)' ) ' '
    text_num = text_num + 1
  end if

  do iface = 1, face_num
    write ( text, '(a1,2x,10i8)' ) 'f', &
      ( face(ivert,iface), ivert = 1, face_order(iface) )
    call s_blanks_delete ( text )
    write ( iunit, '(a)' ) trim ( text )
    text_num = text_num + 1
  end do
!
!  Material binding.
!
  write ( iunit, '(a)' ) 'bind c vertex'
  text_num = text_num + 1
!
!  Material RGB values at each node.
!
  do icor3 = 1, cor3_num
    imat = cor3_material(icor3)
    r = material_rgba(1,imat)
    g = material_rgba(2,imat)
    b = material_rgba(3,imat)
    write ( iunit, '(a,1x,3f6.2)' ) 'c', r, g, b
    text_num = text_num + 1
  end do
!
!  Normal binding.
!
  write ( iunit, '(a)' ) 'bind n vertex'
  text_num = text_num + 1
!
!  Normal vector at each node.
!
  do icor3 = 1, cor3_num
    write ( iunit, '(a,1x,3f6.2)' ) 'n', cor3_normal(1:3,icor3)
    text_num = text_num + 1
  end do

  if ( texture_num > 0 ) then
!
!  Texture filename
!
    write ( iunit, '(a)' ) 'tex ' // trim ( texture_name(1) )
    text_num = text_num + 1
!
!  Texture binding.
!
    write ( iunit, '(a)' ) 'bind r vertex'
    text_num = text_num + 1
!
!  Texture coordinates at each node.
!
    do icor3 = 1, cor3_num
      write ( iunit, '(a,1x,3f6.2)' ) 'r', cor3_tex_uv(1:2,icor3)
      text_num = text_num + 1
    end do

  end if
!
!  Report.
!
  write ( *, * ) ' '
  write ( *, * ) 'SMF_WRITE - Wrote ', text_num, ' text lines.'

  return
end
subroutine sort_heap_external ( n, indx, i, j, isgn )
!
!*******************************************************************************
!
!! SORT_HEAP_EXTERNAL externally sorts a list of items into linear order.
!
!
!  Discussion:
!
!    The actual list of data is not passed to the routine.  Hence this
!    routine may be used to sort integers, reals, numbers, names, 
!    dates, shoe sizes, and so on.  After each call, the routine asks
!    the user to compare or interchange two items, until a special
!    return value signals that the sorting is completed.
!
!  Modified:
!
!    19 May 1999
!
!  Reference:
!
!    A Nijenhuis and H Wilf,
!    Combinatorial Algorithms,
!    Academic Press, 1978, second edition,
!    ISBN 0-12-519260-6.
!
!  Parameters:
!
!    Input, integer N, the number of items to be sorted.
!
!    Input/output, integer INDX, the main communication signal.
!
!    The user must set INDX to 0 before the first call.
!    Thereafter, the user should not change the value of INDX until
!    the sorting is done.
!
!    On return, if INDX is
!
!      greater than 0, 
!      * interchange items I and J;
!      * call again.
!
!      less than 0, 
!      * compare items I and J;
!      * set ISGN = -1 if I < J, ISGN = +1 if I > J;
!      * call again.
!
!      equal to 0, the sorting is done.
!
!    Output, integer I, J, the indices of two items.
!    On return with INDX positive, elements I and J should be interchanged.  
!    On return with INDX negative, elements I and J should be compared, and
!    the result reported in ISGN on the next call.
!
!    Input, integer ISGN, results of comparison of elements I and J.
!    (Used only when the previous call returned INDX less than 0).
!    ISGN <= 0 means I is less than or equal to J;
!    ISGN => 0 means I is greater than or equal to J.
!
  integer i
  integer indx
  integer isgn
  integer j
  integer, save :: k = 0
  integer, save :: k1 = 0
  integer n
  integer, save :: n1 = 0
!
!  INDX = 0: This is the first call.
!
  if ( indx == 0 ) then

    n1 = n
    k = n / 2
    k1 = k
!
!  INDX < 0: The user is returning the results of a comparison.
!
  else if ( indx < 0 ) then

    if ( indx == -2 ) then

      if ( isgn < 0 ) then
        i = i + 1
      end if

      j = k1
      k1 = i
      indx = - 1
      return

    end if

    if ( isgn > 0 ) then
      indx = 2
      return
    end if

    if ( k <= 1 ) then

      if ( n1 == 1 ) then
        indx = 0
      else
        i = n1
        n1 = n1 - 1
        j = 1
        indx = 1
      end if

      return

    end if

    k = k - 1
    k1 = k
!
!  INDX > 0, the user was asked to make an interchange.
!
  else if ( indx == 1 ) then

    k1 = k

  end if

10    continue

  i = 2 * k1

  if ( i == n1 ) then
    j = k1
    k1 = i
    indx = - 1
    return
  else if ( i <= n1 ) then
    j = i + 1
    indx = - 2
    return
  end if

  if ( k > 1 ) then
    k = k - 1
    k1 = k
    go to 10
  end if

  if ( n1 == 1 ) then
    indx = 0
  else
    i = n1
    n1 = n1 - 1
    j = 1
    indx = 1
  end if

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
subroutine stla_write ( cor3, face, face_normal, face_order, filein_name, &
  iunit, cor3_max, face_max, order_max, face_num )
!
!*******************************************************************************
!
!! STLA_WRITE writes graphics information to an ASCII StereoLithography file.
!
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
!  Discussion:
!
!    The polygons in an STL file should only be triangular.  This routine 
!    will try to automatically decompose higher-order polygonal faces into 
!    suitable triangles, without actually modifying the internal graphics 
!    data.
!
!  Modified:
!
!    24 May 1999
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
!    Input, real FACE_NORMAL(3,FACE_MAX), the normal vector at each face.
!
!    Input, integer FACE_ORDER(FACE_MAX), the number of vertices per face.
!
!    Input, character ( len = 100 ) FILEIN_NAME, the name of the input file.
!
!    Input, integer IUNIT, the FORTRAN unit to which output is written.
!
!    Input, integer COR3_MAX, the maximum number of points.
!
!    Input, integer FACE_MAX, the maximum number of faces.
!
!    Input, integer ORDER_MAX, the maximum number of vertices per face.
!
!    Input, integer FACE_NUM, the number of faces.
!
  integer cor3_max
  integer face_max
  integer order_max
!
  real cor3(3,cor3_max)
  integer face(order_max,face_max)
  real face_normal(3,face_max)
  integer face_num
  integer face_num2
  integer face_order(face_max)
  character ( len = * ) filein_name
  integer i
  integer iface
  integer iunit
  integer jvert
  integer node
  integer text_num
!
  text_num = 0

  write ( iunit, '(a,a)' ) 'solid MYSOLID created by IVREAD, ' // &
    'original data in ', trim ( filein_name )

  text_num = text_num + 1
  face_num2 = 0

  do iface = 1, face_num

    do jvert = 3, face_order(iface)

      face_num2 = face_num2 + 1

      write ( iunit, '(''  facet normal '', 3g14.6)' ) face_normal(1:3,iface)
      text_num = text_num + 1

      write ( iunit, '(a)' ) '    outer loop'
      text_num = text_num + 1

      node = face(1,iface)
      write ( iunit, '(''      vertex '', 3g14.6)' ) cor3(1:3,node)
      text_num = text_num + 1

      node = face(jvert-1,iface)
      write ( iunit, '(''      vertex '', 3g14.6)' ) cor3(1:3,node)
      text_num = text_num + 1

      node = face(jvert,iface)
      write ( iunit, '(''      vertex '', 3g14.6)' ) cor3(1:3,node)
      text_num = text_num + 1

      write ( iunit, '(a)' ) '  endloop'
      write ( iunit, '(a)' ) 'endfacet'
      text_num = text_num + 2

    end do

  end do

  write ( iunit, '(a)' ) 'endsolid MYSOLID'
  text_num = text_num + 1
!
!  Report.
!
  write ( *, * ) ' '
  write ( *, * ) 'STLA_WRITE - Wrote ', text_num, ' text lines.'
  if ( face_num2 /= face_num ) then
    write ( *, * ) '  Number of faces in original data was ', face_num
    write ( *, * ) '  Number of triangular faces in decomposed data was ', &
      face_num2
  end if

  return
end
subroutine tec_write ( cor3, cor3_material, face, face_order, fileout_name, iunit, &
  material_rgba, cor3_max, face_max, material_max, order_max, cor3_num, face_num )
!
!*******************************************************************************
!
!! TEC_WRITE writes graphics information to a TECPLOT file.
!
!
!  Discussion:
!
!    The file format used is appropriate for 3D finite element
!    surface zone data.  Polygons are decomposed into triangles where
!    necessary.
!
!  Example:
!
!    TITLE = "cube.tec created by IVREAD."
!    VARIABLES = "X", "Y", "Z", "R", "G", "B"
!    ZONE T="TRIANGLES", N=8, E=12, F=FEPOINT, ET=TRIANGLE
!    0.0 0.0 0.0 0.0 0.0 0.0
!    1.0 0.0 0.0 1.0 0.0 0.0
!    1.0 1.0 0.0 1.0 1.0 0.0
!    0.0 1.0 0.0 0.0 1.0 0.0
!    0.0 0.0 1.0 0.0 0.0 1.0
!    1.0 0.0 1.0 1.0 0.0 1.0
!    1.0 1.0 1.0 1.0 1.0 1.0
!    0.0 1.0 1.0 0.0 1.0 1.0
!    1 4 2
!    2 4 3
!    1 5 8
!    1 2 5
!    2 6 5
!    2 3 6
!    3 7 6
!    3 4 7
!    4 8 7
!    4 1 8
!    5 6 8
!    6 7 8
!
!  Modified:
!
!    08 June 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real COR3(3,COR3_MAX), the coordinates of points.
!
!    Input, integer COR3_MATERIAL(COR3_MAX), the material of each node.
!
!    Input, integer FACE(ORDER_MAX,FACE_MAX), the nodes making faces.
!
!    Input, integer FACE_ORDER(FACE_MAX), the number of vertices per face.
!
!    Input, character ( len = 100 ) FILEOUT_NAME, the name of the output file.
!
!    Input, integer IUNIT, the FORTRAN unit to which output is written.
!
!    Input, real MATERIAL_RGBA(4,MATERIAL_MAX), material R, G, B and A values.
!
!    Input, integer COR3_MAX, the maximum number of points.
!
!    Input, integer FACE_MAX, the maximum number of faces.
!
!    Input, integer MATERIAL_MAX, the maximum number of materials.
!
!    Input, integer ORDER_MAX, the maximum number of vertices per face.
!
!    Input, integer COR3_NUM, the number of points.
!
!    Input, integer FACE_NUM, the number of faces.
!
  integer cor3_max
  integer face_max
  integer material_max
  integer order_max
!
  real b
  real cor3(3,cor3_max)
  integer cor3_material(cor3_max)
  integer cor3_num
  integer face(order_max,face_max)
  integer face2(3)
  integer face_num
  integer face_num2
  integer face_order(face_max)
  character ( len = 100 ) fileout_name
  real g
  integer i
  integer icor3
  integer iface
  integer imat
  integer iunit
  integer jlo
  real material_rgba(4,material_max)
  integer text_num
  real r
!
!  Determine the number of triangular faces.
!
  face_num2 = 0
  do i = 1, face_num
    do jlo = 1, face_order(i) - 2
      face_num2 = face_num2 + 1
    end do
  end do

  text_num = 0

  write ( iunit, '(a)' ) '"' // trim ( fileout_name ) // ' created by IVREAD."'
  write ( iunit, '(a)' ) 'VARIABLES = "X", "Y", "Z", "R", "G", "B"'
  write ( iunit, '(a,i6,a,i6,a)' ) 'ZONE T="TRIANGLES", N=', cor3_num, &
    ', E=', face_num2, ', F=FEPOINT, ET=TRIANGLE'

  text_num = text_num + 3
!
!  Write out X, Y, Z, R, G, B per node.
!
  do icor3 = 1, cor3_num
    imat = cor3_material(icor3)
    r = material_rgba(1,imat)
    g = material_rgba(2,imat)
    b = material_rgba(3,imat)
    write ( iunit, '(6g11.3)' ) cor3(1,icor3), cor3(2,icor3), cor3(3,icor3), &
      r, g, b
    text_num = text_num + 1
  end do
!
!  Do the next face.
!
  do iface = 1, face_num
!
!  Break the face up into triangles, anchored at node 1.
!
    do jlo = 1, face_order(iface) - 2

      face2(1) = face(    1,iface)
      face2(2) = face(jlo+1,iface)
      face2(3) = face(jlo+2,iface)

      write ( iunit, '(3i6)' ) face2(1), face2(2), face2(3)
      text_num = text_num + 1

    end do

  end do
!
!  Report.
!
  write ( *, * ) ' '
  write ( *, * ) 'TEC_WRITE - Wrote ', text_num, ' text lines.'

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
subroutine tmat_mxp ( a, x, y )
!
!*******************************************************************************
!
!! TMAT_MXP multiplies a geometric transformation matrix times a point.
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
!    Input, real A(4,4), the geometric transformation matrix.
!
!    Input, real X(3), the point to be multiplied.  The fourth component
!    of X is implicitly assigned the value of 1.
!
!    Output, real Y(3), the result of A*X.  The product is accumulated in 
!    a temporary vector, and then assigned to the result.  Therefore, it 
!    is legal for X and Y to share memory.
!
  real a(4,4)
  integer i
  integer j
  real x(3)
  real y(3)
  real z(3)
!
  do i = 1, 3
    z(i) = a(i,4)
    do j = 1, 3
      z(i) = z(i) + a(i,j) * x(j)
    end do
  end do

  y(1:3) = z(1:3)

  return
end
subroutine tmat_mxp2 ( a, x, y, n )
!
!*******************************************************************************
!
!! TMAT_MXP2 multiplies a geometric transformation matrix times N points.
!
!
!  Modified:
!
!    20 October 1998
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
!    Input, real X(3,N), the points to be multiplied.  
!
!    Output, real Y(3,N), the transformed points.  Each product is 
!    accumulated in a temporary vector, and then assigned to the
!    result.  Therefore, it is legal for X and Y to share memory.
!
  integer n
!
  real a(4,4)
  integer i
  integer j
  integer k
  real x(3,n)
  real y(3,n)
  real z(3)
!
  do k = 1, n

    do i = 1, 3
      z(i) = a(i,4)
      do j = 1, 3
        z(i) = z(i) + a(i,j) * x(j,k)
      end do
    end do

    y(1:3,k) = z(1:3)

  end do

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
subroutine tmat_rot_axis ( a, b, angle, axis )
!
!*******************************************************************************
!
!! TMAT_ROT_AXIS applies a coordinate axis rotation to the geometric transformation matrix.
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
!    Input, character AXIS, is 'X', 'Y' or 'Z', specifying the coordinate
!    axis about which the rotation occurs.
!
  real a(4,4)
  real angle
  real angle_rad
  character axis
  real b(4,4)
  real c(4,4)
  real d(4,4)
  real degrees_to_radians
  integer i
  integer j
!
  angle_rad = degrees_to_radians ( angle )

  call tmat_init ( c )

  if ( axis == 'X' .or. axis == 'x' ) then
    c(2,2) =   cos ( angle_rad )
    c(2,3) = - sin ( angle_rad )
    c(3,2) =   sin ( angle_rad )
    c(3,3) =   cos ( angle_rad )
  else if ( axis == 'Y' .or. axis == 'y' ) then
    c(1,1) =   cos ( angle_rad )
    c(1,3) =   sin ( angle_rad )
    c(3,1) = - sin ( angle_rad )
    c(3,3) =   cos ( angle_rad )
  else if ( axis == 'Z' .or. axis == 'z' ) then
    c(1,1) =   cos ( angle_rad )
    c(1,2) = - sin ( angle_rad )
    c(2,1) =   sin ( angle_rad )
    c(2,2) =   cos ( angle_rad )
  else
    write ( *, * ) ' '
    write ( *, * ) 'TMAT_ROT_AXIS - Fatal error!'
    write ( *, * ) '  Illegal rotation axis: ', axis
    write ( *, * ) '  Legal choices are ''X'', ''Y'', or ''Z''.'
    return
  end if

  call tmat_mxm ( c, a, d )

  b(1:4,1:4) = d(1:4,1:4)

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
  real a(4,4)
  real angle
  real angle_rad
  real axis(3)
  real b(4,4)
  real c(4,4)
  real ca
  real d(4,4)
  real degrees_to_radians
  integer i
  integer j
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
  real a(4,4)
  real b(4,4)
  real c(4,4)
  real d(4,4)
  integer i
  integer j
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
subroutine tmat_shear ( a, b, axis, s )
!
!*******************************************************************************
!
!! TMAT_SHEAR applies a shear to the geometric transformation matrix.
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
!    Input, character ( len = 2 ) AXIS, is 'XY', 'XZ', 'YX', 'YZ', 'ZX' or 'ZY',
!    specifying the shear equation:
!
!      XY:  x' = x + s * y;
!      XZ:  x' = x + s * z;
!      YX:  y' = y + s * x;
!      YZ:  y' = y + s * z;
!      ZX:  z' = z + s * x;
!      ZY:  z' = z + s * y.
!
!    Input, real S, the shear coefficient.
!
  real a(4,4)
  character ( len = 2 ) axis
  real b(4,4)
  real c(4,4)
  real d(4,4)
  integer i
  integer j
  real s
!
  call tmat_init ( c )

  if ( axis == 'XY' .or. axis == 'xy' ) then
    c(1,2) = s
  else if ( axis == 'XZ' .or. axis == 'xz' ) then
    c(1,3) = s
  else if ( axis == 'YX' .or. axis == 'yx' ) then
    c(2,1) = s
  else if ( axis == 'YZ' .or. axis == 'yz' ) then
    c(2,3) = s
  else if ( axis == 'ZX' .or. axis == 'zx' ) then
    c(3,1) = s
  else if ( axis == 'ZY' .or. axis == 'zy' ) then
    c(3,2) = s
  else
    write ( *, * ) ' '
    write ( *, * ) 'TMAT_SHEAR - Fatal error!'
    write ( *, * ) '  Illegal shear axis: ', axis
    write ( *, * ) '  Legal choices are XY, XZ, YX, YZ, ZX, or ZY.'
    return
  end if

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
subroutine tria_read ( cor3, face, face_material, face_order, ierror, iunit, &
  cor3_max, face_max, order_max, cor3_num, dup_num, face_num, text_num, &
  vertex_material, vertex_normal )
!
!*******************************************************************************
!
!! TRIA_READ reads graphics information from an ASCII triangle file.
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
!    12                    <-- Number of triangles
!
!                          (x,y,z) and (nx,ny,nz) of normal vector at:
!
!    0.0 0.0 0.0 0.3 0.3 0.3   node 1 of triangle 1.
!    1.0 0.0 0.0 0.3 0.1 0.3   node 2 of triangle 1,
!    0.0 1.0 0.0 0.3 0.1 0.3   node 3 of triangle 1,
!    1.0 0.5 0.0 0.3 0.1 0.3   node 1 of triangle 2,
!    ...
!    0.0 0.5 0.5 0.3 0.1 0.3   node 3 of triangle 12.
!
!  Modified:
!
!    22 June 1999
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
!    Input, integer COR3_MAX, the maximum number of points.
!
!    Input, integer FACE_MAX, the maximum number of faces.
!
!    Input, integer ORDER_MAX, the maximum number of vertices per face.
!
!    Input/output, integer COR3_NUM, the number of points.
!
!    Output, integer DUP_NUM, the number of duplicate nodes discovered.
!
!    Input/output, integer FACE_NUM, the number of faces.
!
!    Output, integer TEXT_NUM, the number of lines of input text.
!
!    Input/output, integer VERTEX_MAT(ORDER_MAX,FACE_MAX), vertex materials.
!
!    Input/output, real VERTEX_NORMAL(3,ORDER_MAX,FACE_MAX), normals at vertices.
!
  integer cor3_max
  integer face_max
  integer order_max
!
  real cor3(3,cor3_max)
  integer cor3_num
  real cvec(3)
  integer dup_num
  integer face(order_max,face_max)
  integer face_material(face_max)
  integer face_num
  integer face_num2
  integer face_order(face_max)
  integer i
  integer icor3
  integer ierror
  integer iface
  integer iface_hi
  integer iface_lo
  integer iunit
  integer ivert
  real r1
  real r2
  real r3
  real r4
  real r5
  real r6
  integer text_num
  integer vertex_material(order_max,face_max)
  real vertex_normal(3,order_max,face_max)
!
  ierror = 0
!
!  Read the number of (triangular) faces.
!  (This is added on to the current number, if any).
!
  read ( iunit, *, end = 10 ) face_num2
  text_num = text_num + 1
!
!  For each triangle.
!
  iface_lo = face_num + 1
  iface_hi = face_num + face_num2

  do iface = iface_lo, iface_hi

    if ( iface <= FACE_MAX ) then
      face_order(iface) = 3
      face_material(iface) = 1
    end if
!
!  For each face of a triangle:
!
    do ivert = 1, face_order(iface)

      read ( iunit, *, end = 10 ) r1, r2, r3, r4, r5, r6
      text_num = text_num + 1

      cvec(1) = r1
      cvec(2) = r2
      cvec(3) = r3

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

      if ( iface <= FACE_MAX ) then

        face(ivert,iface) = icor3

        vertex_material(ivert,iface) = 1
        vertex_normal(1,ivert,iface) = r4
        vertex_normal(2,ivert,iface) = r5
        vertex_normal(3,ivert,iface) = r6

      end if

    end do

    face_num = iface

  end do

  return
!
!  Come here on unexpected end of file.
!
10    continue

  ierror = 1
  write ( *, * ) ' '
  write ( *, * ) 'TRIA_READ - Warning.'
  write ( *, * ) '  End-of-file, but the model was not finished.'

  return
end
subroutine tria_write ( cor3, cor3_normal, face, face_order, iunit, cor3_max, &
  face_max, order_max, face_num )
!
!*******************************************************************************
!
!! TRIA_WRITE writes the graphics data to an ASCII "triangle" file.
!
!
!  Discussion:
!
!    This is just a private format that Greg Hood requested from me.
!
!  Example:
!
!    12                        <-- Number of triangles
!
!                              (x,y,z) and (nx,ny,nz) of normal vector at:
!
!    0.0 0.0 0.0 0.3 0.3 0.3   node 1 of triangle 1.
!    1.0 0.0 0.0 0.3 0.1 0.3   node 2 of triangle 1,
!    0.0 1.0 0.0 0.3 0.1 0.3   node 3 of triangle 1,
!    1.0 0.5 0.0 0.3 0.1 0.3   node 1 of triangle 2,
!    ...
!    0.0 0.5 0.5 0.3 0.1 0.3   node 3 of triangle 12.
!
!  Modified:
!
!    08 June 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real COR3(3,COR3_MAX), the coordinates of points.
!
!    Input, real COR3_NORMAL(3,COR3_MAX), the normal vector at each node.
!
!    Input, integer FACE(ORDER_MAX,FACE_MAX), the nodes making faces.
!
!    Input, integer FACE_ORDER(FACE_MAX), the number of vertices per face.
!
!    Input, integer IUNIT, the FORTRAN unit to which output is written.
!
!    Input, integer COR3_MAX, the maximum number of points.
!
!    Input, integer FACE_MAX, the maximum number of faces.
!
!    Input, integer ORDER_MAX, the maximum number of vertices per face.
!
!    Input, integer FACE_NUM, the number of faces.
!
  integer cor3_max
  integer face_max
  integer order_max
!
  real cor3(3,cor3_max)
  real cor3_normal(3,cor3_max)
  integer face(order_max,face_max)
  integer face_num
  integer face_num2
  integer face_order(face_max)
  integer face2(3)
  integer i
  integer icor3
  integer iface
  integer iunit
  integer jlo
  integer k
  real nx
  real ny
  real nz
  integer text_num
  real x
  real y
  real z
!
  text_num = 0
!
!  Determine the number of triangular faces.
!
  face_num2 = 0
  do i = 1, face_num
    do jlo = 1, face_order(i) - 2
      face_num2 = face_num2 + 1
    end do
  end do

  write ( iunit, '(i6)' ) face_num2
  text_num = text_num + 1
!
!  Do the next face.
!
  do iface = 1, face_num
!
!  Break the face up into triangles, anchored at node 1.
!
    do jlo = 1, face_order(iface) - 2

      face2(1) = face(    1,iface)
      face2(2) = face(jlo+1,iface)
      face2(3) = face(jlo+2,iface)

      do k = 1, 3

        icor3 = face2(k)

        x = cor3(1,icor3)
        y = cor3(2,icor3)
        z = cor3(3,icor3)

        nx = cor3_normal(1,icor3)
        ny = cor3_normal(2,icor3)
        nz = cor3_normal(3,icor3)

        write ( iunit, '(6f10.4)' ) x, y, z, nx, ny, nz

        text_num = text_num + 1

      end do

    end do

  end do
!
!  Report.
!
  write ( *, * ) ' '
  write ( *, * ) 'TRIA_WRITE - Wrote ', text_num, ' text lines.'
 
  return
end
subroutine txt_write ( cor3, cor3_material, cor3_normal, cor3_tex_uv, face, &
  face_material, face_normal, face_order, face_tex_uv, filein_name, fileout_name, &
  iunit, line_dex, line_material, material_name, material_rgba, cor3_max, face_max, &
  line_max, material_max, order_max, texture_max, cor3_num, face_num, line_num, &
  material_num, texture_num, object_name, texture_name, vertex_material, vertex_normal, &
  vertex_tex_uv )
!
!*******************************************************************************
!
!! TXT_WRITE writes the graphics data to a text file.
!
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
!    Input, real COR3(3,COR3_MAX), the coordinates of points.
!
!    Input, integer COR3_MATERIAL(COR3_MAX), the material of each node.
!
!    Input, real COR3_NORMAL(3,COR3_MAX), normals at nodes.  
!
!    Input, real COR3_TEX_UV(2,COR3_MAX), UV texture coordinates for nodes.
!
!    Input, integer FACE(ORDER_MAX,FACE_MAX), the nodes making faces.
!
!    Input, integer FACE_MATERIAL(FACE_MAX), the material of each face.
!
!    Input, real FACE_NORMAL(3,FACE_MAX), the normal vector at each face.
!
!    Input, integer FACE_ORDER(FACE_MAX), the number of vertices per face.
!
!    Input, character ( len = 100 ) FILEIN_NAME, the name of the input file.
!
!    Input, character ( len = 100 ) FILEOUT_NAME, the name of the output file.
!
!    Input, integer IUNIT, the FORTRAN unit to which output is written.
!
!    Input, integer LINE_DEX(LINE_MAX), nodes forming a line, terminated by -1.
!
!    Input, integer LINE_MAT(LINE_MAX), material index for each line.
!
!    Input, character ( len = 100 ) MATERIAL_NAME(MATERIAL_MAX), material names.
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
!    Input, integer COR3_NUM, the number of points.
!
!    Input, integer FACE_NUM, the number of faces.
!
!    Input, integer LINE_NUM, the number of line definition items.
!
!    Input, integer MATERIAL_NUM, the number of materials.
!
!    Input, integer TEXTURE_NUM, the number of textures.
!
!    Input, character ( len = 100 ) TEXTURE_NAME(TEXTURE_MAX), texture names.
!
!    Input, integer VERTEX_MAT(ORDER_MAX,FACE_MAX), vertex materials.
!
!    Input, real VERTEX_NORMAL(3,ORDER_MAX,FACE_MAX), normals at vertices.  
!
  integer, parameter :: OFFSET = 1
!
  integer cor3_max
  integer face_max
  integer line_max
  integer material_max
  integer order_max
  integer texture_max
!
  real cor3(3,cor3_max)
  integer cor3_material(cor3_max)
  real cor3_normal(3,cor3_max)
  integer cor3_num
  real cor3_tex_uv(2,cor3_max)
  integer face(order_max,face_max)
  integer face_material(face_max)
  real face_normal(3,face_max)
  integer face_num
  integer face_order(face_max)
  real face_tex_uv(2,face_max)
  character ( len = 100 ) filein_name
  character ( len = 100 ) fileout_name
  integer i
  integer iface
  integer imat
  integer iunit
  integer ivert
  integer j
  integer jhi
  integer jlo
  integer length
  integer line_dex(line_max)
  integer line_material(line_max)
  integer line_num
  character ( len = 100 ) material_name(material_max)
  integer material_num
  real material_rgba(4,material_max)
  character ( len = 100 ) object_name
  character ( len = 100 ) text
  integer text_num
  character ( len = 100 ) texture_name(texture_max)
  integer texture_num
  integer vertex_material(order_max,face_max)
  real vertex_normal(3,order_max,face_max)
  real vertex_tex_uv(2,order_max,face_max)
  character ( len = 10 ) word
!
  text_num = 0

  write ( iunit, * ) trim ( fileout_name ) // ' created by IVREAD.'
  write ( iunit, * ) 'Original data in ' // trim ( filein_name ) // '.'
  write ( iunit, * ) 'Object name is ' // trim ( object_name ) // '.'
  text_num = text_num + 3
!
!  NODES.
!
  write ( iunit, * ) ' '
  write ( iunit, * ) cor3_num, ' nodes:'
  text_num = text_num + 2

  if ( cor3_num > 0 ) then

    write ( iunit, * ) ' '
    write ( iunit, '(a)' ) '    Node  Coordinates:'
    write ( iunit, '(a)' ) '========  ============'
    write ( iunit, * ) ' '
    text_num = text_num + 4

    do j = 1, cor3_num
      write ( iunit, '(i8,2x,3g12.4)' ) j-OFFSET, cor3(1:3,j)
      text_num = text_num + 1
    end do

    write ( iunit, * ) ' '
    write ( iunit, '(a)' ) '    Node  Normal vector'
    write ( iunit, '(a)' ) '========  ============'
    write ( iunit, * ) ' '
    text_num = text_num + 4

    do j = 1, cor3_num
      write ( iunit, '(i8,2x,3g12.4)' ) j-OFFSET, cor3_normal(1:3,j)
      text_num = text_num + 1
    end do

    write ( iunit, * ) ' '
    write ( iunit, '(a)' ) '    Node  Material'
    write ( iunit, '(a)' ) '========  ============'
    write ( iunit, * ) ' '
    text_num = text_num + 4

    do j = 1, cor3_num
      write ( iunit, '(i8,2x,i8)' ) j-OFFSET, cor3_material(j)-OFFSET
      text_num = text_num + 1
    end do

    if ( texture_num > 0 ) then

    write ( iunit, * ) ' '
    write ( iunit, '(a)' ) '    Node  Texture coordinates:'
    write ( iunit, '(a)' ) '========  ===================='
    write ( iunit, * ) ' '
    text_num = text_num + 4

      do j = 1, cor3_num
        write ( iunit, '(i8,2x,3g12.4)' ) j-OFFSET, cor3_tex_uv(1:2,j)
        text_num = text_num + 1
      end do

    end if

  end if
!
!  LINES
!
  write ( iunit, * ) ' '
  write ( iunit, * ) line_num, ' line data items.'
  text_num = text_num + 2

  if ( line_num > 0 ) then

    write ( iunit, * ) ' '
    write ( iunit, * ) '  Line index data:'
    write ( iunit, * ) ' '
    text_num = text_num + 3

    text = ' '
    length = 0
      
    do j = 1, line_num
 
      write ( word, '(i8,'','')' ) line_dex(j) - OFFSET
      call s_cat ( text, word, text )
      length = length + 1

      if ( line_dex(j)-OFFSET == -1 .or. length >= 10 .or. j >= line_num ) then

        call s_blanks_delete ( text )
        write ( iunit, '(8x,a)' ) trim ( text )
        text_num = text_num + 1
        text = ' '
        length = 0

      end if

    end do

    write ( iunit, * ) ' '
    write ( iunit, * ) 'Line materials:'
    write ( iunit, * ) ' '
    text_num = text_num + 3

    text = ' '
    length = 0
     
    do j = 1, line_num
 
      write ( word, '(i8,'','')' ) line_material(j) - OFFSET
      call s_cat ( text, word, text )
      length = length + 1

      if ( line_material(j)-OFFSET == -1 .or. length >= 10 .or. j >= line_num ) then

        call s_blanks_delete ( text )
        write ( iunit, '(8x,a)' ) trim ( text )
        text_num = text_num + 1
        text = ' '
        length = 0

      end if

    end do

  end if
!
!  FACES
!
  write ( iunit, * ) ' '
  write ( iunit, * ) face_num, ' faces:'
  text_num = text_num + 2

  if ( face_num > 0 ) then

    write ( iunit, * ) ' '
    write ( iunit, '(a)' ) '    Face  Material     Order'
    write ( iunit, '(a)' ) '========  ========  ========'
    write ( iunit, * ) ' '
    text_num = text_num + 4

    do iface = 1, face_num
      write ( iunit, '(i8,2x,i8,2x,i8)' ) iface-OFFSET, &
        face_material(iface)-OFFSET, face_order(iface)
      text_num = text_num + 1
    end do

    write ( iunit, * ) ' '
    write ( iunit, '(a)' ) '    Face  Vertices:'
    write ( iunit, '(a)' ) '========  ========================'
    write ( iunit, * ) ' '
    text_num = text_num + 4

    do iface = 1, face_num

      do jlo = 1, face_order(iface), 10
        jhi = min ( jlo + 9, face_order(iface) )
        if ( jlo == 1 ) then
          write ( iunit, '(i8,2x,10i8)' ) iface-OFFSET, &
                ( face(ivert,iface)-OFFSET, ivert = jlo, jhi )
        else
          write ( iunit, '(10x,10i8)' ) &
                ( face(ivert,iface)-OFFSET, ivert = jlo, jhi )
        end if
        text_num = text_num + 1
      end do

    end do

    write ( iunit, * ) ' '
    write ( iunit, '(a)' ) '    Face  Vertex Material Indices'
    write ( iunit, '(a)' ) '========  ========================'
    write ( iunit, * ) ' '
    text_num = text_num + 4

    do iface = 1, face_num

      do jlo = 1, face_order(iface), 10
        jhi = min ( jlo + 9, face_order(iface) )
        if ( jlo == 1 ) then
          write ( iunit, '(i8,2x,10i8)' ) iface-OFFSET, &
                ( vertex_material(ivert,iface)-OFFSET, ivert = jlo, jhi )
        else
          write ( iunit, '(10x,10i8)' ) &
                ( vertex_material(ivert,iface)-OFFSET, ivert = jlo, jhi )
        end if
        text_num = text_num + 1
      end do

    end do

    write ( iunit, * ) ' '
    write ( iunit, '(a)' ) '    Face  Normal vector'
    write ( iunit, '(a)' ) '========  ========================'
    write ( iunit, * ) ' '
    text_num = text_num + 4

    do iface = 1, face_num
      write ( iunit, '(i8,3g12.4)' ) iface-OFFSET, face_normal(1:3,iface)
      text_num = text_num + 1
    end do

    if ( texture_num > 0 ) then

      write ( iunit, * ) ' '
      write ( iunit, '(a)' ) '    Face  Texture coordinates:'
      write ( iunit, '(a)' ) '========  ========================'
      write ( iunit, * ) ' '
      text_num = text_num + 4

      do iface = 1, face_num
        write ( iunit, '(i8,2x,3g12.4)' ) iface-OFFSET, face_tex_uv(1:2,iface)
        text_num = text_num + 1
      end do

    end if

  end if
!
!  VERTICES.
!
  if ( face_num > 0 ) then

    write ( iunit, * ) ' '
    write ( iunit, '(a)' ) '    Face    Vertex  Normal vector:'
    write ( iunit, '(a)' ) '========  ========  =============='
    write ( iunit, * ) ' '
    text_num = text_num + 4

    do iface = 1, face_num
      write ( iunit, '(a)' ) ' '
      text_num = text_num + 1
      do ivert = 1, face_order(iface)
        write ( iunit, '(i8,2x,i8,2x,3f12.4)' ) iface-OFFSET, ivert-OFFSET, &
          vertex_normal(1:3,ivert,iface)
        text_num = text_num + 1
      end do
    end do

    if ( texture_num > 0 ) then

    write ( iunit, * ) ' '
    write ( iunit, '(a)' ) '    Face    Vertex  Texture coordinates:'
    write ( iunit, '(a)' ) '========  ========  =================='
    write ( iunit, * ) ' '
    text_num = text_num + 4

      do iface = 1, face_num
        write ( iunit, '(a)' ) ' '
        text_num = text_num + 1
        do ivert = 1, face_order(iface)
          write ( iunit, '(i8,2x,i8,2x,3f12.4)' ) iface-OFFSET, ivert-OFFSET, &
            vertex_tex_uv(1:2,ivert,iface)
          text_num = text_num + 1
        end do
      end do

    end if

  end if
!
!  MATERIALS
!
  write ( iunit, * ) ' '
  write ( iunit, * ) material_num, ' materials:'
  write ( iunit, * ) ' '
  text_num = text_num + 3

  if ( material_num > 0 ) then
    write ( iunit, '(a)' ) 'Material        Name            ' // &
      '    R         G         B         A'
    write ( iunit, '(a)' ) '========  ====================  ' // &
      '========================================'
    write ( iunit, * ) ' '
    do imat = 1, material_num
      write ( iunit, '(i8,2x,a20,2x,4f10.6)' ) imat-OFFSET, &
        material_name(imat)(1:20), material_rgba(1:4,imat)
    end do
  end if
!
!  TEXTURES:
!
  write ( iunit, * ) ' '
  write ( iunit, * ) texture_num, ' textures:'
  text_num = text_num + 2

  if ( texture_num > 0 ) then

    write ( iunit, * ) ' '
    write ( iunit, '(a)' ) ' Texture    Name'
    write ( iunit, '(a)' ) '========  ========'
    write ( iunit, * ) ' '
    text_num = text_num + 4

    do i = 1, texture_num
      write ( iunit, '(i6,2x,a)' ) i-OFFSET, trim ( texture_name(i) )
      text_num = text_num + 1
    end do
  end if
!
!  Report.
!
  write ( *, * ) ' '
  write ( *, * ) 'TXT_WRITE - Wrote ', text_num, ' text lines.'
 
  return
end
subroutine ucd_write ( cor3, cor3_material, face, face_material, face_order, &
  fileout_name, iunit, material_rgba, cor3_max, face_max, material_max, &
  order_max, cor3_num, face_num, material_num )
!
!*******************************************************************************
!
!! UCD_WRITE writes graphics data to an AVS UCD file.
!
!
!  Example:
!
!    #  cube.ucd created by IVREAD.
!    #
!    #  Material RGB to hue map:
!    #
!    #  material    R    G      B   Alpha     Hue
!    #
!    #    0       0.94  0.70  0.15  1.000   0.116
!    #    1       0.24  0.70  0.85  1.000   0.541
!    #    2       0.24  0.00  0.85  1.000   0.666
!    #
!    #  The node data is
!    #    node # / material # / RGBA / Hue
!    #
!    8  6  6  0  0
!    0  0.0  0.0  0.0
!    1  1.0  0.0  0.0
!    2  1.0  1.0  0.0
!    3  0.0  1.0  0.0
!    4  0.0  0.0  1.0
!    5  1.0  0.0  1.0
!    6  1.0  1.0  1.0
!    7  0.0  1.0  1.0
!    0  0  quad  0  1  2  3
!    1  0  quad  0  4  5  1 
!    2  0  quad  1  5  6  2
!    3  0  quad  2  6  7  3
!    4  0  quad  3  7  4  0
!    5  0  quad  4  7  6  5
!    3  1 4 1
!    material, 0...2
!    RGBA, 0-1/0-1/0-1/0-1
!    Hue, 0-1
!    0  0  0.94  0.70  0.15  1.0  0.116
!    1  0  0.94  0.70  0.15  1.0  0.116
!    2  0  0.94  0.70  0.15  1.0  0.116
!    3  0  0.94  0.70  0.15  1.0  0.116
!    4  1  0.24  0.70  0.85  1.0  0.541
!    5  1  0.24  0.70  0.85  1.0  0.541
!    6  2  0.24  0.24  0.85  0.0  0.666
!    7  2  0.24  0.24  0.85  0.0  0.666
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
!    Input, real COR3(3,COR3_MAX), the coordinates of points.
!
!    Input, integer FACE(ORDER_MAX,FACE_MAX), the nodes making faces.
!
!    Input, integer FACE_MATERIAL(FACE_MAX), the material of each face.
!
!    Input, integer FACE_ORDER(FACE_MAX), the number of vertices per face.
!
!    Input, character ( len = 100 ) FILEOUT_NAME, the name of the output file.
!
!    Input, integer IUNIT, the FORTRAN unit to which output is written.
!
!    Input, integer LINE_DEX(LINE_MAX), nodes forming a line, terminated by -1.
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
!    Input, integer COR3_NUM, the number of points.
!
!    Input, integer FACE_NUM, the number of faces.
!
!    Input, integer LINE_NUM, the number of line definition items.
!
!    Input, integer MATERIAL_NUM, the number of materials.
!
!    Input, real VERTEX_NORMAL(3,ORDER_MAX,FACE_MAX), normals at vertices.
!
  integer, parameter :: OFFSET = 1
!
  integer cor3_max
  integer face_max
  integer material_max
  integer order_max
!
  real a
  real b
  character ( len = 4 ) cell_type
  real cor3(3,cor3_max)
  integer cor3_material(cor3_max)
  integer cor3_num
  integer face(order_max,face_max)
  integer face_material(face_max)
  integer face_num
  integer face_order(face_max)
  character ( len = 100 ) fileout_name
  real g
  real h
  integer i
  integer imat
  integer iunit
  integer j
  integer material_num
  real material_rgba(4,material_max)
  real r
  character ( len = 100 ) text
  integer text_num
!
  text_num = 0

  write ( iunit, '(a)' ) '#  ' // trim ( fileout_name ) // ' created by IVREAD.'
  write ( iunit, '(a)' ) '#'
  write ( iunit, '(a)' ) '#  Material RGB to Hue map:'
  write ( iunit, '(a)' ) '#'
  write ( iunit, '(a)' ) '#  material    R    G      B     Alpha  Hue'
  write ( iunit, '(a)' ) '#'
  text_num = text_num + 6

  do j = 1, material_num
    r = material_rgba(1,j)
    g = material_rgba(2,j)
    b = material_rgba(3,j)
    a = material_rgba(4,j)
    call rgb_to_hue ( r, g, b, h )
    write ( text, '(a,2x,i4,5f7.3)' ) '#  ', j-OFFSET, r, g, b, a, h
    call s_blanks_delete ( text )
    write ( iunit, '(a)' ) trim ( text )
    text_num = text_num + 1
  end do

  write ( iunit, '(a)' ) '#'
  write ( iunit, '(a)' ) '#  The node data is'
  write ( iunit, '(a)' ) '#    node # / material # / RGBA / Hue'
  write ( iunit, '(a)' ) '#'
  text_num = text_num + 4

  write ( text, '(i6,2x,i6,2x,a)' ) cor3_num, face_num, '6  0  0'
  call s_blanks_delete ( text )
  write ( iunit, '(a)' ) trim ( text )
  text_num = text_num + 1

  do j = 1, cor3_num
    write ( text, '(i6,2x,3f7.3)' ) j - OFFSET, cor3(1:3,j)
    call s_blanks_delete ( text )
    write ( iunit, '(a)' ) trim ( text )
    text_num = text_num + 1
  end do
!
!  NOTE:
!  UCD only accepts triangles and quadrilaterals, not higher order
!  polygons.  We would need to break polygons up to proceed.
!
!  Also, we just use the material of vertex 1 of the face as the
!  material of the face.
!
  do j = 1, face_num

    if ( face_order(j) == 3 ) then
      cell_type = 'tri'
    else if ( face_order(j) .eq. 4 ) then
      cell_type = 'quad'
    else
      cell_type = '???'
    end if

    imat = face_material(j) - OFFSET
    write ( text, '(i6,2x,i6,2x,a4,2x,10i6)' ) j-OFFSET, imat, &
      cell_type, ( face(i,j) - OFFSET, i = 1, face_order(j) )
    call s_blanks_delete ( text )
    write ( iunit, '(a)' ) trim ( text )
    text_num = text_num + 1

  end do

  write ( iunit, '(a)' )  '3  1  4  1'
  write ( iunit, '(a,i6)' )  'material, 0...', material_num - OFFSET
  write ( iunit, '(a)' )  'RGBA, 0-1/0-1/0-1/0-1'
  write ( iunit, '(a)' )  'Hue, 0-1'
  text_num = text_num + 4

  do j = 1, cor3_num
    imat = cor3_material(j)
    r = material_rgba(1,imat)
    g = material_rgba(2,imat)
    b = material_rgba(3,imat)
    a = material_rgba(4,imat)
    call rgb_to_hue ( r, g, b, h )
    write ( text, '(i6,2x,i6,2x,5f7.3)' ) j-OFFSET, imat-OFFSET, r, g, b, a, h
    call s_blanks_delete ( text )
    write ( iunit, '(a)' ) trim ( text )
    text_num = text_num + 1
  end do
!
!  Report.
!
  write ( *, * ) ' '
  write ( *, * ) 'UCD_WRITE - Wrote ', text_num, ' text lines.'

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
subroutine vector_unit_nd ( n, v )
!
!*******************************************************************************
!
!! VECTOR_UNIT_ND normalizes a vector in ND.
!
!
!  Modified:
!
!    07 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries in the vector.
!
!    Input/output, real V(N), the vector to be normalized.  On output,
!    V should have unit Euclidean norm.  However, if the input vector
!    has zero Euclidean norm, it is not altered.
!
  integer n
!
  real enorm_nd
  integer i
  real temp
  real v(n)
!
  temp = enorm_nd ( n, v )
 
  if ( temp /= 0.0 ) then
    v(1:n) = v(1:n) / temp
  end if
 
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
subroutine vertex_to_node_material ( cor3_material, face, face_order, cor3_max, &
  face_max, order_max, face_num, vertex_material )
!
!*******************************************************************************
!
!! VERTEX_TO_NODE_MATERIAL extends vertex material definitions to nodes.
!
!
!  Discussion:
!
!    A NODE is a point in space.
!    A VERTEX is a node as used in a particular face.
!    One node may be used as a vertex in several faces, or none.
!    This routine simply runs through all the vertices, and assigns
!    the material of the vertex to the corresponding node.  If a
!    node appears as a vertex several times, then the node will 
!  end up having the material of the vertex that occurs "last".
!
!  Modified:
!
!    27 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer COR3_MATERIAL(COR3_MAX), the material index of each node.
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
!    Input, integer VERTEX_MAT(ORDER_MAX,FACE_MAX), vertex materials.
!
  integer cor3_max
  integer face_max
  integer order_max
!
  integer cor3_material(cor3_max)
  integer face(order_max,face_max)
  integer face_order(face_max)
  integer iface
  integer ivert
  integer node
  integer face_num
  integer vertex_material(order_max,face_max)
!
  do iface = 1, face_num
    do ivert = 1, face_order(iface)
      node = face(ivert,iface)
      cor3_material(node) = vertex_material(ivert,iface)
    end do
  end do

  return
end
subroutine vla_read ( cor3, ierror, iunit, line_dex, line_material, cor3_max, &
  line_max, bad_num, cor3_num, dup_num, line_num, material_num, text_num ) 
!
!*******************************************************************************
!
!! VLA_READ reads graphics information from a VLA file.
!
!
!  Comments:
!
!    Internal comments begin with a semicolon in column 1.
!
!    The X, Y, Z coordinates of points begin with a "P" to
!    denote the beginning of a line, and "L" to denote the
!    continuation of a line.  The fourth entry is intensity, which 
!    should be between 0.0 and 1.0.
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
!     set comment cube.vla created by IVREAD
!     set comment from data in file cube.iv
!     set comment
!     set intensity EXPLICIT
!     set parametric NON_PARAMETRIC
!     set filecontent LINES
!     set filetype NEW
!     set depthcue 0
!     set defaultdraw stellar
!     set coordsys RIGHT
!     set author IVREAD
!     set site Buhl Planetarium
!     set library_id UNKNOWN
!     ; DXF LINE entity
!     P   8.59816       5.55317      -3.05561       1.00000
!     L   8.59816       2.49756      0.000000E+00   1.00000
!     L   8.59816       2.49756      -3.05561       1.00000
!     L   8.59816       5.55317      -3.05561       1.00000
!     ; DXF LINE entity
!     P   8.59816       5.55317      0.000000E+00   1.00000
!     ...etc...
!     L   2.48695       2.49756      -3.05561       1.00000
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
!    Output, integer IERROR, an error flag.
!
!    Input, integer IUNIT, the FORTRAN unit from which data is read.
!
!    Output, integer LINE_DEX(LINE_MAX), nodes forming a line, terminated by -1.
!
!    Output, integer LINE_MAT(LINE_MAX), material index for each line.
!
!    Input, integer COR3_MAX, the maximum number of points.
!
!    Input, integer LINE_MAX, the maximum number of line items.
!
!    Output, integer COR3_NUM, the number of points.
!
!    Output, integer LINE_NUM, the number of line definition items.
!
  integer, parameter :: OFFSET = 1

  integer cor3_max
  integer line_max
!
  integer bad_num
  real cor3(3,cor3_max)
  integer cor3_num
  real cvec(3)
  logical done
  integer dup_num
  integer i
  integer icor3
  integer ierror
  integer iunit
  integer iword
  integer lchar
  integer line_dex(line_max)
  integer line_material(line_max)
  integer line_num
  integer material_num
  character prevcode
  real rval
  logical s_eqi
  character ( len = 256 ) text
  integer text_num
  character ( len = 256 ) word
  character ( len = 256 ) word1
!
  ierror = 0
  prevcode = 'P'

10    continue
 
  read ( iunit, '(a)', end = 30 ) text
  text_num = text_num + 1
 
  done = .true.
  iword = 0
 
20    continue
 
  call word_nexrd ( text, word, done )
!
!  If no more words in this line, read a new line.
!
  if ( done ) then
    go to 10
  end if
 
  iword = iword + 1
!
!  The first word in the line tells us what's happening.
!
  if ( iword == 1 ) then
    word1 = word
  end if
!
!  If WORD1 is "SET", then we regard this line as comments.
!
  if ( s_eqi ( word1, 'set' ) ) then
!
!  If WORD1 is ";", then we regard this line as comments.
!
  else if ( word1 == ';' ) then
!
!  If WORD1 is "P", then this is the initial point on a line.
!  If WORD1 is "L", then this is a followup point on a line.
!
  else if ( s_eqi ( word1, 'P' ) .or. s_eqi ( word1, 'L' ) ) then
!
!  Terminate the current line if the new code is 'P'.
!
    if ( s_eqi ( prevcode, 'L' ) .and. s_eqi ( word1, 'P' ) ) then

      line_num = line_num + 1

      if ( line_num <= line_max ) then
        line_dex(line_num) = -1 + OFFSET
        line_material(line_num) = -1 + OFFSET
      end if

    end if

    prevcode = word1
!
!  Read in the coordinates of the point, stored temporarily in CVEC.
!
    do i = 1, 3

      call word_nexrd ( text, word, done )

      if ( done ) then
        bad_num = bad_num + 1
        go to 10
      end if

      call s_to_r ( word, rval, ierror, lchar )

      if ( ierror /= 0 ) then
        bad_num = bad_num + 1
        go to 10
      end if

      if ( cor3_num <= cor3_max ) then
        cvec(i) = rval
      end if

    end do
!
!  If the values in CVEC already exist in COR3, then  
!  save space by using the index of a previous copy.
!
!  Otherwise, add CVEC to COR3, and increment COR3_NUM.
!
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
!
!  Now define the line.
!
    line_num = line_num + 1

    if ( line_num <= line_max ) then
      line_dex(line_num) = icor3 - 1 + OFFSET
      line_material(line_num) = material_num
    end if

    go to 10
!
!  If the first word is unrecognized, then skip the whole line.
!
  else

    bad_num = bad_num + 1
    go to 10

  end if
 
  go to 20

30    continue
!
!  Terminate the very last line.
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
subroutine vla_write ( cor3, filein_name, fileout_name, iunit, line_dex, &
  cor3_max, line_max, line_num )
!
!*******************************************************************************
!
!! VLA_WRITE writes graphics data to a VLA file.
!
!
!  Discussion:
!
!    Comments begin with a semicolon in column 1.
!    The X, Y, Z coordinates of points begin with a "P" to
!    denote the beginning of a line, and "L" to denote the
!    continuation of a line.  The fourth entry is intensity, which 
!    should be between 0.0 and 1.0.
!
!  Example:
!
!     set comment cube.vla created by IVREAD
!     set comment Original data in cube.iv.
!     set comment
!     set intensity EXPLICIT
!     set parametric NON_PARAMETRIC
!     set filecontent LINES
!     set filetype NEW
!     set depthcue 0
!     set defaultdraw stellar
!     set coordsys RIGHT
!     set author IVREAD
!     set site Buhl Planetarium
!     set library_id UNKNOWN
!     P   8.59816       5.55317      -3.05561       1.00000
!     L   8.59816       2.49756      0.000000E+00   1.00000
!     L   8.59816       2.49756      -3.05561       1.00000
!     L   8.59816       5.55317      -3.05561       1.00000
!     P   8.59816       5.55317      0.000000E+00   1.00000
!     ...etc...
!     L   2.48695       2.49756      -3.05561       1.00000
!
!  Modified:
!
!    18 August 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real COR3(3,COR3_MAX), the coordinates of points.
!
!    Input, character ( len = 100 ) FILEIN_NAME, the name of the input file.
!
!    Input, character ( len = 100 ) FILEOUT_NAME, the name of the output file.
!
!    Input, integer IUNIT, the FORTRAN unit to which output is written.
!
!    Input, integer LINE_DEX(LINE_MAX), nodes forming a line, terminated by -1.
!
!    Input, integer COR3_MAX, the maximum number of points.
!
!    Input, integer LINE_MAX, the maximum number of line definition items.
!
!    Output, integer LINE_NUM, the number of line definition items.
!
  integer, parameter :: OFFSET = 1
!
  integer cor3_max
  integer line_max
!
  character code
  real cor3(3,cor3_max)
  character ( len = 100 ) filein_name
  character ( len = 100 ) fileout_name
  integer i
  real intense
  integer iunit
  integer j
  integer k
  integer line_dex(line_max)
  integer line_num
  logical newline
  integer nl
  integer np
  integer text_num
!
  intense = 1.00
  np = 0
  nl = 0
 
  write ( iunit, '(a)' ) 'set comment ' // trim ( fileout_name ) // &
    ' created by IVREAD'

  write ( iunit, '(a)' ) 'set comment Original data in ' &
    // trim ( filein_name ) // '.'

  write ( iunit, '(a)' ) 'set comment'
  write ( iunit, '(a)' ) 'set intensity EXPLICIT'
  write ( iunit, '(a)' ) 'set parametric NON_PARAMETRIC'
  write ( iunit, '(a)' ) 'set filecontent LINES'
  write ( iunit, '(a)' ) 'set filetype NEW'
  write ( iunit, '(a)' ) 'set depthcue 0'
  write ( iunit, '(a)' ) 'set defaultdraw stellar'
  write ( iunit, '(a)' ) 'set coordsys RIGHT'
  write ( iunit, '(a)' ) 'set author IVREAD'
  write ( iunit, '(a)' ) 'set site Buhl Planetarium'
  write ( iunit, '(a)' ) 'set library_id UNKNOWN'

  text_num = 13
!
!  Print the line data.
!
  newline = .TRUE.

  do j = 1, line_num
 
    k = line_dex(j) - OFFSET

    if ( k == - 1 ) then

      newline = .TRUE.

    else 

      if ( newline ) then
        code = 'P'
        np = np + 1
      else
        code = 'L'
        nl = nl + 1
      end if

      write ( iunit, '(a,4g12.4)' ) code, cor3(1:3,k+OFFSET), intense

      text_num = text_num + 1
      newline = .FALSE.

    end if

  end do
!
!  Report.
!
  write ( *, * ) ' '
  write ( *, * ) 'VLA_WRITE - Wrote ', text_num, ' text lines.'
 
  return
end
subroutine vrml_read ( cor3, cor3_material, face, face_material, face_order, ierror, &
  iunit, line_dex, line_material, material_name, material_rgba, cor3_max, &
  face_max, line_max, material_max, order_max, texture_max, bad_num, cor3_num, &
  face_num, line_num, material_num, texture_num, text_num, texture_name, vertex_material, &
  vertex_normal, vertex_tex_uv )
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
  logical, parameter :: debug = .false.
  integer, parameter :: ivec_max = 20
  integer, parameter :: level_max = 10
  integer, parameter :: rvec_max = 20
  integer, parameter :: offset = 1
!
  integer cor3_max
  integer face_max
  integer line_max
  integer material_max
  integer order_max
  integer texture_max
!
  real angle
  real axis(3)
  integer bad_num
  character ( len = 4 ) char4
  real cor3(3,cor3_max)
  integer cor3_material(cor3_max)
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
  integer line_dex(line_max)
  integer line_material(line_max)
  integer line_num
  integer line_num_old
  character ( len = 100 ) material_binding
  character ( len = 100 ) material_name(material_max)
  integer material_num
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
  character ( len = 100 ) texture_name(texture_max)
  integer texture_num
  real transform_matrix(4,4)
  real tx
  real ty
  real tz
  integer vertex_material(order_max,face_max)
  real vertex_normal(3,order_max,face_max)
  real vertex_tex_uv(2,order_max,face_max)
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

          if ( material_num <= material_max ) then
            new_color = new_color + 1
            material_rgba(1,material_num) = rvec(1)
            material_rgba(2,material_num) = rvec(2)
            material_rgba(3,material_num) = rvec(3)
            material_rgba(4,material_num) = 1.0
            call i_to_s_zero ( material_num, char4 )
            material_name(material_num) = 'Material_' // char4
          end if

        end if

        rvec_num = 0
        level = nlbrack - nrbrack

      else

        call s_to_r ( word, rval, ierror, lchar )
        rvec_num = rvec_num + 1
        rvec(rvec_num) = rval

        if ( rvec_num == 3 ) then

          material_num = material_num + 1

          if ( material_num <= material_max ) then
            new_color = new_color + 1
            material_rgba(1,material_num) = rvec(1)
            material_rgba(2,material_num) = rvec(2)
            material_rgba(3,material_num) = rvec(3)
            material_rgba(4,material_num) = 1.0
            call i_to_s_zero ( material_num, char4 )
            material_name(material_num) = 'Material_' // char4
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

      if ( material_num == 0 ) then

        material_num = material_num + 1

        material_rgba(1,material_num) = 1.0
        material_rgba(2,material_num) = 0.0
        material_rgba(3,material_num) = 0.0
        material_rgba(4,material_num) = 1.0
        call i_to_s_zero ( material_num, char4 )
        material_name(material_num) = 'Material_' // char4

      end if

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
          cor3_material(icor3) = imat
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
          face_material(iface) = imat
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

      material_num = material_num + 1

      if ( material_num <= material_max ) then
        material_rgba(1,material_num) = r02
        material_rgba(2,material_num) = r03
        material_rgba(3,material_num) = r04
        material_rgba(4,material_num) = 1.0 - r12
        call i_to_s_zero ( material_num, char4 )
        material_name(material_num) = 'Material_' // char4
      end if

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
99    continue

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
subroutine vrml_write ( cor3, face, face_order, filein_name, fileout_name, &
  iunit, line_dex, line_material,  material_rgba, cor3_max, face_max, line_max, &
  material_max, order_max, cor3_num, face_num, line_num, material_num, vertex_material )
!
!*******************************************************************************
!
!! VRML_WRITE writes graphics data to a VRML file.
!
!
!  The VRML files written by this routine have the form:
!
!
!     #VRML V2.0 utf8
!
!       WorldInfo {
!         title "cube.iv."
!         string "VRML file generated by IVREAD.
!       }
!
!       Group {
!         children [
!
!           Shape {
!
!             appearance Appearance {
!               material Material {
!                 diffuseColor   0.0 0.0 0.0
!                 emissiveColor  0.0 0.0 0.0
!                 shininess      1.0
!               }
!             } #end of appearance
!
!             geometry IndexedLineSet {
!
!               coord Coordinate {
!                 point [
!                   8.59816       5.55317      -3.05561
!                   8.59816       2.49756      0.000000E+00
!                   ...etc...
!                   2.48695       2.49756      -3.05561
!                 ]
!               }
!
!               coordIndex [
!                   0     1     2    -1     3     4     5     6     7     8    -1
!                   9    10    -1    11    12    -1    13    14    15    -1    16
!                 ...etc...
!                 191    -1
!               ]
!
!               colorPerVertex TRUE
!
!               colorIndex [
!                   0     0     0    -1     2     3     1     1     4     7    -1
!                  10     9    -1     7     7    -1     3     2     2    -1    12
!                 ...etc...
!                 180    -1
!               ]
!
!             }  #end of geometry
!
!           }  #end of Shape
!
!         ]  #end of children
!       }  #end of Group
!
!  Modified:
!
!    23 April 1998
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
!    Input, character ( len = 100 ) FILEIN_NAME, the name of the input file.
!
!    Input, character ( len = 100 ) FILEOUT_NAME, the name of the output file.
!
!    Input, integer IUNIT, the FORTRAN unit to which output is written.
!
!    Input, integer LINE_DEX(LINE_MAX), nodes forming a line, terminated by -1.
!
!    Input, integer LINE_MAT(LINE_MAX), material index for line.
!
!    Input, integer COR3_MAX, the maximum number of points.
!
!    Input, integer FACE_MAX, the maximum number of faces.
!
!    Input, integer LINE_MAX, the maximum number of line definition items.
!
!    Input, integer ORDER_MAX, the maximum number of vertices per face.
!
!    Input, integer COR3_NUM, the number of points.
!
!    Input, integer FACE_NUM, the number of faces.
!
!    Input, integer LINE_NUM, the number of line definition items.
!
  integer, parameter :: OFFSET = 1

  integer cor3_max
  integer face_max
  integer line_max
  integer material_max
  integer order_max
!
  real cor3(3,cor3_max)
  integer cor3_num
  integer face(order_max,face_max)
  integer face_num
  integer face_order(face_max)
  character ( len = 100 ) filein_name
  character ( len = 100 ) fileout_name
  integer i
  integer icor3
  integer iface
  integer itemp
  integer iunit
  integer ivert
  integer j
  integer length
  integer line_dex(line_max)
  integer line_material(line_max)
  integer line_num
  integer material_num
  real material_rgba(4,material_max)
  integer ndx
  character ( len = 100 ) text
  integer text_num
  integer vertex_material(order_max,face_max)
  character ( len = 10 ) word
!
  write ( iunit, '(a)' ) '#VRML V2.0 utf8'
  write ( iunit, '(a)' ) ' '
  write ( iunit, '(a)' ) '  WorldInfo {'
  write ( iunit, '(a)' ) '    title "' // trim ( fileout_name ) //'"'
  write ( iunit, '(a)' ) '    info "VRML file generated by IVREAD."'
  write ( iunit, '(a)' ) '    info "Original data in file ' // &
    trim ( filein_name ) // '."'
  write ( iunit, '(a)' ) '  }'
  write ( iunit, '(a)' ) ' '
  write ( iunit, '(a)' ) '  Group {'
  write ( iunit, '(a)' ) '    children ['
  write ( iunit, '(a)' ) '      Shape {'
  write ( iunit, '(a)' ) '        appearance Appearance {'
  write ( iunit, '(a)' ) '          material Material {'
  write ( iunit, '(a)' ) '            diffuseColor   0.0 0.0 0.0'
  write ( iunit, '(a)' ) '            emissiveColor  0.0 0.0 0.0'
  write ( iunit, '(a)' ) '            shininess      1.0'
  write ( iunit, '(a)' ) '          }'
  write ( iunit, '(a)' ) '        } '

  text_num = 18
!
!  IndexedLineSet
!
  if ( line_num > 0 ) then

    write ( iunit, '(a)' ) '        geometry IndexedLineSet {'
!
!  IndexedLineSet coord
!
    write ( iunit, '(a)' ) '          coord Coordinate {'
    write ( iunit, '(a)' ) '            point ['
 
    text_num = text_num + 3
 
    do icor3 = 1, cor3_num
      write ( text, '(3f12.4,'','')' ) cor3(1:3,icor3)
      call s_blanks_delete ( text )
      write ( iunit, '(8x,a)' ) trim ( text )
      text_num = text_num + 1
    end do

    write ( iunit, '(a)' ) '            ]'
    write ( iunit, '(a)' ) '          }' 
    text_num = text_num + 2
!
!  IndexedLineSet coordIndex.
!
    write ( iunit, '(a)' ) '          coordIndex ['

    text_num = text_num + 1

    text = ' '
    length = 0
    
    do j = 1, line_num
 
      write ( word, '(i8,'','')' ) line_dex(j) - OFFSET
      call s_cat ( text, word, text )
      length = length + 1

      if ( line_dex(j)-OFFSET == -1 .or. length >= 10 .or. j >= line_num ) then

        call s_blanks_delete ( text )
        write ( iunit, '(8x,a)' ) trim ( text )
        text_num = text_num + 1
        text = ' '
        length = 0

      end if

    end do
 
    write ( iunit, '(a)' ) '          ]'
    text_num = text_num + 1
!
!  Colors. (materials)
!
    write ( iunit, '(a)' ) '          color Color {'
    write ( iunit, '(a)' ) '            color ['
    text_num = text_num + 2

    do j = 1, material_num
      write ( iunit, '(3f12.4,'','')' ) material_rgba(1:3,j)
      text_num = text_num + 1
    end do

    write ( iunit, '(a)' ) '            ]'
    write ( iunit, '(a)' ) '          }'
    write ( iunit, '(a)' ) '          colorPerVertex TRUE'
!
!  IndexedLineset colorIndex
!
    write ( iunit, '(a)' ) '          colorIndex ['

    text_num = text_num + 4

    text = ' '
    length = 0
    
    do j = 1, line_num
 
      write ( word, '(i8,'','')' ) line_material(j) - OFFSET
      call s_cat ( text, word, text )
      length = length + 1

      if ( line_dex(j)-OFFSET == -1 .or. length >= 10 .or. j >= line_num ) then

        call s_blanks_delete ( text )
        write ( iunit, '(8x,a)' ) trim ( text )
        text_num = text_num + 1
        text = ' '
        length = 0

      end if

    end do

    if ( text /= ' ' ) then
      write ( iunit, '(a)' ) trim ( text )
      text_num = text_num + 1
      text = ' '
    end if

    write ( iunit, '(a)' ) '          ]'
    write ( iunit, '(a)' ) '        }'
    text_num = text_num + 2

  end if
!
!  End of IndexedLineSet
!
!
!  IndexedFaceSet
!
  if ( face_num > 0 ) then

    write ( iunit, '(a)' ) '        geometry IndexedFaceSet {'
!
!  IndexedFaceSet coord
!
    write ( iunit, '(a)' ) '          coord Coordinate {'
    write ( iunit, '(a)' ) '            point ['
 
    text_num = text_num + 3
 
    do icor3 = 1, cor3_num
      write ( text, '(3f12.4,'','')' ) cor3(1:3,icor3)
      call s_blanks_delete ( text )
      write ( iunit, '(8x,a)' ) trim ( text )
      text_num = text_num + 1
    end do

    write ( iunit, '(a)' ) '            ]'
    write ( iunit, '(a)' ) '          }' 
!
!  IndexedFaceSet coordIndex.
!
    write ( iunit, '(a)' ) '          coordIndex ['

    text_num = text_num + 3

    text = ' '
    length = 0
 
    do iface = 1, face_num

      do ivert = 1, face_order(iface) + 1

        if ( ivert <= face_order(iface) ) then
          itemp = face(ivert,iface) - OFFSET
        else
          itemp = 0 - OFFSET
        end if

        write ( word, '(i8,'','')' ) itemp

        call s_cat ( text, word, text )
        length = length + 1

        if ( itemp == -1 .or. length >= 10 .or. &
           ( iface == face_num .and. ivert == face_order(iface) + 1 ) ) then

          call s_blanks_delete ( text )
          write ( iunit, '(8x,a)' ) trim ( text )
          text_num = text_num + 1
          text = ' '
          length = 0

        end if
 
      end do

    end do
 
    write ( iunit, '(a)' ) '          ]'
    text_num = text_num + 1
!
!  IndexedFaceset colorIndex
!
    write ( iunit, '(a)' ) '          colorIndex ['

    text_num = text_num + 4

    text = ' '
    length = 0
    ndx = 0
 
    do iface = 1, face_num
 
      do ivert = 1, face_order(iface) + 1

        if ( ivert <= face_order(iface) ) then
          itemp = vertex_material(ivert,iface) - OFFSET
          ndx = ndx + 1
        else
          itemp = 0 - OFFSET
        end if

        write ( word, '(i8,'','')' ) itemp

        call s_cat ( text, word, text )
        length = length + 1

        if ( itemp == -1 .or.  length >= 10 .or. &
           ( iface == face_num .and.  ivert == face_order(iface) + 1 )  ) then

          call s_blanks_delete ( text )
          write ( iunit, '(8x,a)' ) trim ( text )
          text_num = text_num + 1
          text = ' '
          length = 0
  
        end if

      end do
 
    end do

    if ( text /= ' ' ) then
      write ( iunit, '(a)' ) trim ( text )
      text_num = text_num + 1
      text = ' '
    end if

    write ( iunit, '(a)' ) '          ]'
    write ( iunit, '(a)' ) '        }'
    text_num = text_num + 2

  end if
!
!  End of IndexedFaceSet
!
!  End of:
!        Shape
!      children
!    Group
!
  write ( iunit, '(a)' ) '      }'
  write ( iunit, '(a)' ) '    ]'
  write ( iunit, '(a)' ) '  }'
 
  text_num = text_num + 3
!
!  Report.
!
  write ( *, * ) ' '
  write ( *, * ) 'VRML_WRITE - Wrote ', text_num, ' text lines.'

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
subroutine xgl_write ( cor3, cor3_max, cor3_normal, cor3_num, face, &
  face_material, face_max, face_num, face_order, iunit, material_max, &
  material_num, material_rgba, order_max )
!
!*******************************************************************************
!
!! XGL_WRITE writes graphics data to an XGL file.
!
!
!  Discussion:
!
!    Only triangular faces are allowed.
!
!  Example:
!
!    <WORLD>
!
!      <BACKGROUND>
!        <BACKCOLOR>r,g,b</BACKCOLOR>
!      </BACKGROUND>
! 
!      <LIGHTING>
!        <AMBIENT>r,g,b</AMBIENT>
!        <DIRECTIONALLIGHT>
!          <DIFFUSE>r,g,b</DIFFUSE>
!          <DIRECTION>x,y,z</DIRECTION>
!          <SPECULAR>r,g,b</SPECULAR>
!        </DIRECTIONALLIGHT>
!      </LIGHTING>
!
!      <MESH ID="0">
!
!        <P ID = "0"> x, y, z </P>
!        ...
!        <P ID = "99"> x, y, z </P>
!
!        <N ID = "0"> nx, ny, nz </N>
!        ...
!        <N ID = "55"> nx, ny, nz </N>
!
!        <MAT ID = "0">
!          <ALPHA></ALPHA>
!          <AMB>r,g,b</AMB>
!          <DIFF>r,g,b</DIFF>
!          <EMISS>r,g,b</EMISS>
!          <SHINE></SHINE>
!          <SPEC>r,g,b</SPEC>
!        </MAT>
!
!        <F>
!          <MATREF>0</MATREF>
!          <FV1><PREF>0</PREF><NREF>0</NREF></FV1>
!          <FV2><PREF>1</PREF><NREF>1</NREF></FV2>
!          <FV3><PREF>2</PREF><NREF>2</NREF></FV3>
!        </F>
!        ...
!        <F>
!          <MATREF>0</MATREF>
!          <FV1><PREF>12</PREF><NREF>12</NREF></FV1>
!          <FV2><PREF>13</PREF><NREF>13</NREF></FV2>
!          <FV3><PREF>14</PREF><NREF>14</NREF></FV3>
!        </F>
!
!       </MESH>
!
!       <OBJECT>
!         <TRANSFORM>
!           <FORWARD>x,y,z</FORWARD>
!           <POSITION>x,y,z</POSITION>
!           <SCALE>x,y,z</SCALE>
!           <UP>x,y,z</UP>
!         </TRANSFORM>
!         <MESHREF>0</MESHREF>
!       </OBJECT>
!
!       <OBJECT>
!         <TRANSFORM>
!           <FORWARD>x,y,z</FORWARD>
!           <POSITION>x,y,z</POSITION>
!           <SCALE>x,y,z</SCALE>
!           <UP>x,y,z</UP>
!         </TRANSFORM>
!         <MESHREF>0</MESHREF>
!       </OBJECT>
!
!     </WORLD>
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
!    Input, real COR3(3,COR3_MAX), the coordinates of points.
!
!    Input, integer COR3_MAX, the maximum number of points.
!
!    Input, real COR3_NORMAL(3,COR3_MAX), normals at nodes.
!
!    Input, integer COR3_NUM, the number of points.
!
!    Input, integer FACE(ORDER_MAX,FACE_MAX), the nodes making faces.
!
!    Input, integer FACE_MAX, the maximum number of faces.
!
!    Input, integer FACE_NUM, the number of faces.
!
!    Input, integer FACE_ORDER(FACE_MAX), the number of vertices per face.
!
!    Input, integer IUNIT, the FORTRAN unit to which output is written.
!
!    Input, integer MATERIAL_MAX, the maximum number of materials.
!
!    Input, integer MATERIAL_NUM, the number of materials.
!
!    Input, real MATERIAL_RGBA(4,MATERIAL_MAX), material R, G, B and A values.
!
!    Input, integer ORDER_MAX, the maximum number of vertices per face.
!
  integer, parameter :: OFFSET = 1

  integer cor3_max
  integer face_max
  integer material_max
  integer order_max
!
  real background_rgb(3)
  real cor3(3,cor3_max)
  real cor3_normal(3,cor3_max)
  integer cor3_num
  integer face(order_max,face_max)
  integer face_material(face_max)
  integer face_num
  integer face_order(face_max)
  integer i
  integer iface
  integer iunit
  integer ivert
  integer j
  real light_ambient_rgb(3)
  real light_diffuse_rgb(3)
  real light_direction(3)
  real light_specular_rgb(3)
  real material_alpha(1)
  real material_amb_rgb(3)
  real material_diff_rgb(3)
  real material_emiss_rgb(3)
  real material_shine(1)
  real material_spec_rgb(3)
  integer material
  integer material_num
  real material_rgba(4,material_max)
  integer mesh
  integer, parameter :: mesh_num = 1
  integer, parameter :: object_num = 1
  character ( len = 100 ) s
  character s1
  character ( len = 6 ) s2
  real transform_forward(3)
  real transform_position(3)
  real transform_scale(3)
  real transform_up(3)
!
!  Define some dummy values.
!
  background_rgb(1) = 0.1
  background_rgb(2) = 0.1
  background_rgb(3) = 0.1

  light_ambient_rgb(1) = 0.2
  light_ambient_rgb(2) = 0.1
  light_ambient_rgb(3) = 0.1

  light_diffuse_rgb(1) = 0.1
  light_diffuse_rgb(2) = 0.2
  light_diffuse_rgb(3) = 0.1

  light_direction(1) =   0.0
  light_direction(2) =   0.0
  light_direction(3) = 100.0

  light_specular_rgb(1) = 0.1
  light_specular_rgb(2) = 0.1
  light_specular_rgb(3) = 0.2

  material_alpha = 0.9

  material_amb_rgb(1) = 0.1
  material_amb_rgb(2) = 0.1
  material_amb_rgb(3) = 0.1

  material_diff_rgb(1) = 0.2
  material_diff_rgb(2) = 0.1
  material_diff_rgb(3) = 0.1

  material_emiss_rgb(1) = 0.1
  material_emiss_rgb(2) = 0.2
  material_emiss_rgb(3) = 0.1

  material_shine = 0.8

  material_spec_rgb(1) = 0.1
  material_spec_rgb(2) = 0.1
  material_spec_rgb(3) = 0.2

  transform_forward(1) = 0.0
  transform_forward(2) = 0.0
  transform_forward(3) = 0.0

  transform_position(1) = 0.0
  transform_position(2) = 0.0
  transform_position(3) = 0.0

  transform_scale(1) = 1.0
  transform_scale(2) = 1.0
  transform_scale(3) = 1.0

  transform_up(1) = 1.0
  transform_up(2) = 1.0
  transform_up(3) = 1.0
!
  write ( iunit, '(a)' ) '<WORLD>'
  write ( iunit, '(a)' ) ' '

  write ( iunit, '(a)' ) '  <BACKGROUND>'
  call rvec_to_s ( 3, background_rgb, s )
  write ( iunit, '(a)' ) '    <BACKCOLOR> ' // trim ( s ) // ' </BACKCOLOR>'
  write ( iunit, '(a)' ) '  </BACKGROUND>'
  write ( iunit, '(a)' ) ' '
  write ( iunit, '(a)' ) '  <LIGHTING>'
  call rvec_to_s ( 3, light_ambient_rgb, s )
  write ( iunit, '(a)' ) '    <AMBIENT> ' // trim ( s ) // ' </AMBIENT>'
  write ( iunit, '(a)' ) '    <DIRECTIONALLIGHT>'
  call rvec_to_s ( 3, light_diffuse_rgb, s )
  write ( iunit, '(a)' ) '      <DIFFUSE> ' // trim ( s ) // ' </DIFFUSE>'
  call rvec_to_s ( 3, light_direction, s )
  write ( iunit, '(a)' ) '      <DIRECTION> ' // trim ( s ) // ' </DIRECTION>'
  call rvec_to_s ( 3, light_specular_rgb, s )
  write ( iunit, '(a)' ) '      <SPECULAR> ' // trim ( s ) // ' </SPECULAR>'
  write ( iunit, '(a)' ) '    </DIRECTIONALLIGHT>'
  write ( iunit, '(a)' ) '  </LIGHTING>'

  do mesh = 1, mesh_num

    write ( iunit, '(a)' ) ' '
    write ( s2, '(i6)' ) mesh - OFFSET
    s2 = adjustl ( s2 )
    write ( iunit, '(a)' ) '  <MESH ID = "' // trim ( s2 ) // '">'

    write ( iunit, '(a)' ) ' '
    do j = 1, cor3_num
      write ( s2, '(i6)' ) j - OFFSET
      s2 = adjustl ( s2 )
      call rvec_to_s ( 3, cor3(1,j), s )
      write ( iunit, '(a)' ) '    <P ID="' // trim ( s2 ) // '"> ' // &
        trim ( s ) // ' </P>'
    end do

    write ( iunit, '(a)' ) ' '
    do j = 1, cor3_num
      write ( s2, '(i6)' ) j - OFFSET
      s2 = adjustl ( s2 )
      call rvec_to_s ( 3, cor3_normal(1,j), s )
      write ( iunit, '(a)' ) '    <N ID="' // trim ( s2 ) // '"> ' // &
        trim ( s ) // ' </N>'
    end do

    do material = 1, material_num
      write ( iunit, '(a)' ) ' '
      write ( s2, '(i6)' ) material - OFFSET
      s2 = adjustl ( s2 )
      write ( iunit, '(a)' ) '    <MAT ID="' // trim ( s2 ) // '">'
      call rvec_to_s ( 1, material_alpha, s )
      write ( iunit, '(a)' ) '      <ALPHA> ' // trim ( s ) // ' </ALPHA>'
      call rvec_to_s ( 3, material_amb_rgb, s )
      write ( iunit, '(a)' ) '      <AMB> ' // trim ( s ) // ' </AMB>'
      call rvec_to_s ( 3, material_diff_rgb, s )
      write ( iunit, '(a)' ) '      <DIFF> ' // trim ( s ) // ' </DIFF>'
      call rvec_to_s ( 3, material_emiss_rgb, s )
      write ( iunit, '(a)' ) '      <EMISS> ' // trim ( s ) // ' </EMISS>'
      call rvec_to_s ( 1, material_shine, s )
      write ( iunit, '(a)' ) '      <SHINE> ' // trim ( s ) // ' </SHINE>'
      call rvec_to_s ( 3, material_spec_rgb, s )
      write ( iunit, '(a)' ) '      <SPEC> ' // trim ( s ) // ' </SPEC>'
      write ( iunit, '(a)' ) '    </MAT>'
    end do

    write ( iunit, '(a)' ) ' '

    do iface = 1, face_num
      write ( iunit, '(a)' ) '    <F>'
      write ( s2, '(i6)' ) face_material(iface) - OFFSET
      s2 = adjustl ( s2 )
      write ( iunit, '(a)' ) '      <MATREF> ' // trim ( s2 ) // ' </MATREF>'
      do ivert = 1, face_order(iface)
        write ( s1, '(i1)' ) ivert
        write ( s2, '(i6)' ) face(ivert,iface) - OFFSET
        s2 = adjustl ( s2 )
        write ( iunit, '(a)' ) '      <FV' // trim ( s1 ) // '><PREF> ' // &
          trim ( s2 ) // ' </PREF><NREF> ' // &
          trim ( s2 ) // ' </NREF></FV' // trim ( s1 ) // '>'
      end do
      write ( iunit, '(a)' ) '    </F>'
    end do

    write ( iunit, '(a)' ) '  </MESH>'

  end do

  write ( iunit, '(a)' ) ' '

  do i = 1, object_num

    write ( iunit, '(a)' ) '  <OBJECT>'
    write ( iunit, '(a)' ) '    <TRANSFORM>'
    call rvec_to_s ( 3, transform_forward, s )
    write ( iunit, '(a)' ) '      <FORWARD> ' // trim ( s ) // ' </FORWARD>'
    call rvec_to_s ( 3, transform_position, s )
    write ( iunit, '(a)' ) '      <POSITION> ' // trim ( s ) // ' </POSITION>'
    call rvec_to_s ( 3, transform_scale, s )
    write ( iunit, '(a)' ) '      <SCALE> ' // trim ( s ) // ' </SCALE>'
    call rvec_to_s ( 3, transform_up, s )
    write ( iunit, '(a)' ) '      <UP> ' // trim ( s ) // ' </UP>'
    write ( iunit, '(a)' ) '    </TRANSFORM>'
    mesh = 1
    write ( s2, '(i6)' ) mesh - OFFSET
    s2 = adjustl ( s2 )
    write ( iunit, '(a)' ) '    <MESHREF> ' // trim ( s2 ) // ' </MESHREF>'
    write ( iunit, '(a)' ) '  </OBJECT>'

  end do

  write ( iunit, '(a)' ) ' '
  write ( iunit, '(a)' ) '</WORLD>'

  return
end
