subroutine bjndd ( n, x, bj, dj, fj )

!*****************************************************************************80
!
!! BJNDD computes Bessel functions Jn(x) and first and second derivatives.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However,
!    they give permission to incorporate this routine into a user program
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    11 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order.
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Output, real ( kind = 8 ) BJ(N+1), DJ(N+1), FJ(N+1), the values of
!    Jn(x), Jn'(x) and Jn''(x) in the last entries.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) bj(n+1)
  real ( kind = 8 ) bs
  real ( kind = 8 ) dj(n+1)
  real ( kind = 8 ) f
  real ( kind = 8 ) f0
  real ( kind = 8 ) f1
  real ( kind = 8 ) fj(n+1)
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mt
  integer ( kind = 4 ) nt
  real ( kind = 8 ) x

  do nt = 1, 900
    mt = int ( 0.5D+00 * log10 ( 6.28D+00 * nt ) &
      - nt * log10 ( 1.36D+00 * abs ( x ) / nt ) )
    if ( 20 < mt ) then
      exit
    end if
  end do

  m = nt
  bs = 0.0D+00
  f0 = 0.0D+00
  f1 = 1.0D-35
  do k = m, 0, -1
    f = 2.0D+00 * ( k + 1.0D+00 ) * f1 / x - f0
    if ( k <= n ) then
      bj(k+1) = f
    end if
    if ( k == 2 * int ( k / 2 ) ) then
      bs = bs + 2.0D+00 * f
    end if
    f0 = f1
    f1 = f
  end do

  do k = 0, n
    bj(k+1) = bj(k+1) / ( bs - f )
  end do

  dj(1) = -bj(2)
  fj(1) = -1.0D+00 * bj(1) - dj(1) / x
  do k = 1, n
    dj(k+1) = bj(k) - k * bj(k+1) / x
    fj(k+1) = ( k * k / ( x * x ) - 1.0D+00 ) * bj(k+1) - dj(k+1) / x
  end do

  return
end
subroutine cik01 ( z, cbi0, cdi0, cbi1, cdi1, cbk0, cdk0, cbk1, cdk1 )

!*****************************************************************************80
!
!! CIK01: modified Bessel I0(z), I1(z), K0(z) and K1(z) for complex argument.
!
!  Discussion:
!
!    This procedure computes the modified Bessel functions I0(z), I1(z),
!    K0(z), K1(z), and their derivatives for a complex argument.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However,
!    they give permission to incorporate this routine into a user program
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    31 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, complex ( kind = 8 ) Z, the argument.
!
!    Output, complex ( kind = 8 ) CBI0, CDI0, CBI1, CDI1, CBK0, CDK0, CBK1,
!    CDK1, the values of I0(z), I0'(z), I1(z), I1'(z), K0(z), K0'(z), K1(z),
!    and K1'(z).
!
  implicit none

  real ( kind = 8 ), save, dimension ( 12 ) :: a = (/ &
    0.125D+00,           7.03125D-02,&
    7.32421875D-02,      1.1215209960938D-01,&
    2.2710800170898D-01, 5.7250142097473D-01,&
    1.7277275025845D+00, 6.0740420012735D+00,&
    2.4380529699556D+01, 1.1001714026925D+02,&
    5.5133589612202D+02, 3.0380905109224D+03 /)
  real ( kind = 8 ) a0
  real ( kind = 8 ), save, dimension ( 10 ) :: a1 = (/ &
    0.125D+00,            0.2109375D+00, &
    1.0986328125D+00,     1.1775970458984D+01, &
    2.1461706161499D+002, 5.9511522710323D+03, &
    2.3347645606175D+05,  1.2312234987631D+07, &
    8.401390346421D+08,   7.2031420482627D+10 /)
  real ( kind = 8 ), save, dimension ( 12 ) :: b = (/ &
   -0.375D+00,           -1.171875D-01, &
   -1.025390625D-01,     -1.4419555664063D-01, &
   -2.7757644653320D-01, -6.7659258842468D-01, &
   -1.9935317337513D+00, -6.8839142681099D+00, &
   -2.7248827311269D+01, -1.2159789187654D+02, &
   -6.0384407670507D+02, -3.3022722944809D+03 /)
  complex ( kind = 8 ) ca
  complex ( kind = 8 ) cb
  complex ( kind = 8 ) cbi0
  complex ( kind = 8 ) cbi1
  complex ( kind = 8 ) cbk0
  complex ( kind = 8 ) cbk1
  complex ( kind = 8 ) cdi0
  complex ( kind = 8 ) cdi1
  complex ( kind = 8 ) cdk0
  complex ( kind = 8 ) cdk1
  complex ( kind = 8 ) ci
  complex ( kind = 8 ) cr
  complex ( kind = 8 ) cs
  complex ( kind = 8 ) ct
  complex ( kind = 8 ) cw
  integer ( kind = 4 ) k
  integer ( kind = 4 ) k0
  real ( kind = 8 ) pi
  real ( kind = 8 ) w0
  complex ( kind = 8 ) z
  complex ( kind = 8 ) z1
  complex ( kind = 8 ) z2
  complex ( kind = 8 ) zr
  complex ( kind = 8 ) zr2

  pi = 3.141592653589793D+00
  ci = cmplx ( 0.0D+00, 1.0D+00, kind = 8 )
  a0 = abs ( z )
  z2 = z * z
  z1 = z

  if ( a0 == 0.0D+00 ) then
    cbi0 = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
    cbi1 = cmplx ( 0.0D+00, 0.0D+00, kind = 8 )
    cdi0 = cmplx ( 0.0D+00, 0.0D+00, kind = 8 )
    cdi1 = cmplx ( 0.5D+00, 0.0D+00, kind = 8 )
    cbk0 = cmplx ( 1.0D+30, 0.0D+00, kind = 8 )
    cbk1 = cmplx ( 1.0D+30, 0.0D+00, kind = 8 )
    cdk0 = - cmplx ( 1.0D+30, 0.0D+00, kind = 8 )
    cdk1 = - cmplx ( 1.0D+30, 0.0D+00, kind = 8 )
    return
  end if

  if ( real ( z, kind = 8 ) < 0.0D+00 ) then
    z1 = -z
  end if

  if ( a0 <= 18.0D+00 ) then

    cbi0 = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
    cr = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
    do k = 1, 50
      cr = 0.25D+00 * cr * z2 / ( k * k )
      cbi0 = cbi0 + cr
      if ( abs ( cr / cbi0 ) < 1.0D-15 ) then
        exit
      end if
    end do

    cbi1 = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
    cr = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
    do k = 1, 50
      cr = 0.25D+00 * cr * z2 / ( k * ( k + 1 ) )
      cbi1 = cbi1 + cr
      if ( abs ( cr / cbi1 ) < 1.0D-15 ) then
        exit
      end if
    end do

    cbi1 = 0.5D+00 * z1 * cbi1

  else

    if ( a0 < 35.0D+00 ) then
      k0 = 12
    else if ( a0 < 50.0D+00 ) then
      k0 = 9
    else
      k0 = 7
    end if

    ca = exp ( z1 ) / sqrt ( 2.0D+00 * pi * z1 )
    cbi0 = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
    zr = 1.0D+00 / z1
    do k = 1, k0
      cbi0 = cbi0 + a(k) * zr ** k
    end do
    cbi0 = ca * cbi0
    cbi1 = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
    do k = 1, k0
      cbi1 = cbi1 + b(k) * zr ** k
    end do
    cbi1 = ca * cbi1

  end if

  if ( a0 <= 9.0D+00 ) then

    cs = cmplx ( 0.0D+00, 0.0D+00, kind = 8 )
    ct = - log ( 0.5D+00 * z1 ) - 0.5772156649015329D+00
    w0 = 0.0D+00
    cr = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
    do k = 1, 50
      w0 = w0 + 1.0D+00 / k
      cr = 0.25D+00 * cr / ( k * k ) * z2
      cs = cs + cr * ( w0 + ct )
      if ( abs ( ( cs - cw ) / cs ) < 1.0D-15 ) then
        exit
      end if
      cw = cs
    end do

    cbk0 = ct + cs

  else

    cb = 0.5D+00 / z1
    zr2 = 1.0D+00 / z2
    cbk0 = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
    do k = 1, 10
      cbk0 = cbk0 + a1(k) * zr2 ** k
    end do
    cbk0 = cb * cbk0 / cbi0

  end if

  cbk1 = ( 1.0D+00 / z1 - cbi1 * cbk0 ) / cbi0

  if ( real ( z, kind = 8 ) < 0.0D+00 ) then

    if ( imag ( z ) < 0.0D+00 ) then
      cbk0 = cbk0 + ci * pi * cbi0
      cbk1 = - cbk1 + ci * pi * cbi1
    else
      cbk0 = cbk0 - ci * pi * cbi0
      cbk1 = - cbk1 - ci * pi * cbi1
    end if

    cbi1 = - cbi1

  end if

  cdi0 = cbi1
  cdi1 = cbi0 - 1.0D+00 / z * cbi1
  cdk0 = - cbk1
  cdk1 = - cbk0 - 1.0D+00 / z * cbk1

  return
end
subroutine jy01b ( x, bj0, dj0, bj1, dj1, by0, dy0, by1, dy1 )

!*****************************************************************************80
!
!! JY01B computes Bessel functions J0(x), J1(x), Y0(x), Y1(x) and derivatives.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    02 August 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Output, real ( kind = 8 ) BJ0, DJ0, BJ1, DJ1, BY0, DY0, BY1, DY1,
!    the values of J0(x), J0'(x), J1(x), J1'(x), Y0(x), Y0'(x), Y1(x), Y1'(x).
!
  implicit none

  real ( kind = 8 ) a0
  real ( kind = 8 ) bj0
  real ( kind = 8 ) bj1
  real ( kind = 8 ) by0
  real ( kind = 8 ) by1
  real ( kind = 8 ) dj0
  real ( kind = 8 ) dj1
  real ( kind = 8 ) dy0
  real ( kind = 8 ) dy1
  real ( kind = 8 ) p0
  real ( kind = 8 ) p1
  real ( kind = 8 ) pi
  real ( kind = 8 ) q0
  real ( kind = 8 ) q1
  real ( kind = 8 ) t
  real ( kind = 8 ) t2
  real ( kind = 8 ) ta0
  real ( kind = 8 ) ta1
  real ( kind = 8 ) x

  pi = 3.141592653589793D+00

  if ( x == 0.0D+00 ) then

    bj0 = 1.0D+00
    bj1 = 0.0D+00
    dj0 = 0.0D+00
    dj1 = 0.5D+00
    by0 = -1.0D+300
    by1 = -1.0D+300
    dy0 = 1.0D+300
    dy1 = 1.0D+300
    return

  else if ( x <= 4.0D+00 ) then

    t = x / 4.0D+00
    t2 = t * t

    bj0 = (((((( &
      - 0.5014415D-03 * t2 &
      + 0.76771853D-02 ) * t2 &
      - 0.0709253492D+00 ) * t2 &
      + 0.4443584263D+00 ) * t2 &
      - 1.7777560599D+00 ) * t2 &
      + 3.9999973021D+00 ) * t2 &
      - 3.9999998721D+00 ) * t2 &
      + 1.0D+00

    bj1 = t * ((((((( &
      - 0.1289769D-03 * t2 &
      + 0.22069155D-02 ) * t2 &
      - 0.0236616773D+00 ) * t2 &
      + 0.1777582922D+00 ) * t2 &
      - 0.8888839649D+00 ) * t2 &
      + 2.6666660544D+00 ) * t2 &
      - 3.9999999710D+00 ) * t2 &
      + 1.9999999998D+00 )

    by0 = ((((((( &
      - 0.567433D-04 * t2 &
      + 0.859977D-03 ) * t2 &
      - 0.94855882D-02 ) * t2 &
      + 0.0772975809D+00 ) * t2 &
      - 0.4261737419D+00 ) * t2 &
      + 1.4216421221D+00 ) * t2 &
      - 2.3498519931D+00 ) * t2 &
      + 1.0766115157D+00 ) * t2 &
      + 0.3674669052D+00

    by0 = 2.0D+00 / pi * log ( x / 2.0D+00 ) * bj0 + by0

    by1 = (((((((( &
        0.6535773D-03 * t2 &
      - 0.0108175626D+00 ) * t2 &
      + 0.107657606D+00 ) * t2 &
      - 0.7268945577D+00 ) * t2 &
      + 3.1261399273D+00 ) * t2 &
      - 7.3980241381D+00 ) * t2 &
      + 6.8529236342D+00 ) * t2 &
      + 0.3932562018D+00 ) * t2 &
      - 0.6366197726D+00 ) / x

    by1 = 2.0D+00 / pi * log ( x / 2.0D+00 ) * bj1 + by1

  else

    t = 4.0D+00 / x
    t2 = t * t
    a0 = sqrt ( 2.0D+00 / ( pi * x ) )

    p0 = (((( &
      - 0.9285D-05 * t2 &
      + 0.43506D-04 ) * t2 &
      - 0.122226D-03 ) * t2 &
      + 0.434725D-03 ) * t2 &
      - 0.4394275D-02 ) * t2 &
      + 0.999999997D+00

    q0 = t * ((((( &
        0.8099D-05 * t2 &
      - 0.35614D-04 ) * t2 &
      + 0.85844D-04 ) * t2 &
      - 0.218024D-03 ) * t2 &
      + 0.1144106D-02 ) * t2 &
      - 0.031249995D+00 )

    ta0 = x - 0.25D+00 * pi
    bj0 = a0 * ( p0 * cos ( ta0 ) - q0 * sin ( ta0 ) )
    by0 = a0 * ( p0 * sin ( ta0 ) + q0 * cos ( ta0 ) )

    p1 = (((( &
        0.10632D-04 * t2 &
      - 0.50363D-04 ) * t2 &
      + 0.145575D-03 ) * t2 &
      - 0.559487D-03 ) * t2 &
      + 0.7323931D-02 ) * t2 &
      + 1.000000004D+00

    q1 = t * ((((( &
      - 0.9173D-05      * t2 &
      + 0.40658D-04 )   * t2 &
      - 0.99941D-04 )   * t2 &
      + 0.266891D-03 )  * t2 &
      - 0.1601836D-02 ) * t2 &
      + 0.093749994D+00 )

    ta1 = x - 0.75D+00 * pi
    bj1 = a0 * ( p1 * cos ( ta1 ) - q1 * sin ( ta1 ) )
    by1 = a0 * ( p1 * sin ( ta1 ) + q1 * cos ( ta1 ) )

  end if

  dj0 = - bj1
  dj1 = bj0 - bj1 / x
  dy0 = - by1
  dy1 = by0 - by1 / x

  return
end
subroutine jyna ( n, x, nm, bj, dj, by, dy )

!*****************************************************************************80
!
!! JYNA computes Bessel functions Jn(x) and Yn(x) and derivatives.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    29 April 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order.
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Output, integer ( kind = 4 ) NM, the highest order computed.
!
!    Output, real ( kind = 8 ) BJ(0:N), DJ(0:N), BY(0:N), DY(0:N), the values
!    of Jn(x), Jn'(x), Yn(x), Yn'(x).
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) bj(0:n)
  real ( kind = 8 ) bj0
  real ( kind = 8 ) bj1
  real ( kind = 8 ) bjk
  real ( kind = 8 ) by(0:n)
  real ( kind = 8 ) by0
  real ( kind = 8 ) by1
  real ( kind = 8 ) cs
  real ( kind = 8 ) dj(0:n)
  real ( kind = 8 ) dj0
  real ( kind = 8 ) dj1
  real ( kind = 8 ) dy(0:n)
  real ( kind = 8 ) dy0
  real ( kind = 8 ) dy1
  real ( kind = 8 ) f
  real ( kind = 8 ) f0
  real ( kind = 8 ) f1
  real ( kind = 8 ) f2
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  integer ( kind = 4 ) msta1
  integer ( kind = 4 ) msta2
  integer ( kind = 4 ) nm
  real ( kind = 8 ) x

  nm = n

  if ( x < 1.0D-100 ) then

    do k = 0, n
      bj(k) = 0.0D+00
      dj(k) = 0.0D+00
      by(k) = -1.0D+300
      dy(k) = 1.0D+300
    end do
    bj(0) = 1.0D+00
    dj(1) = 0.5D+00
    return

  end if

  call jy01b ( x, bj0, dj0, bj1, dj1, by0, dy0, by1, dy1 )
  bj(0) = bj0
  bj(1) = bj1
  by(0) = by0
  by(1) = by1
  dj(0) = dj0
  dj(1) = dj1
  dy(0) = dy0
  dy(1) = dy1

  if ( n <= 1 ) then
    return
  end if

  if ( n < int ( 0.9D+00 * x) ) then

    do k = 2, n
      bjk = 2.0D+00 * ( k - 1.0D+00 ) / x * bj1 - bj0
      bj(k) = bjk
      bj0 = bj1
      bj1 = bjk
    end do

  else

    m = msta1 ( x, 200 )

    if ( m < n ) then
      nm = m
    else
      m = msta2 ( x, n, 15 )
    end if

    f2 = 0.0D+00
    f1 = 1.0D-100
    do k = m, 0, -1
      f = 2.0D+00 * ( k + 1.0D+00 ) / x * f1 - f2
      if ( k <= nm ) then
        bj(k) = f
      end if
      f2 = f1
      f1 = f
    end do

    if ( abs ( bj1 ) < abs ( bj0 ) ) then
      cs = bj0 / f
    else
      cs = bj1 / f2
    end if

    do k = 0, nm
      bj(k) = cs * bj(k)
    end do

  end if

  do k = 2, nm
    dj(k) = bj(k-1) - k / x * bj(k)
  end do

  f0 = by(0)
  f1 = by(1)
  do k = 2, nm
    f = 2.0D+00 * ( k - 1.0D+00 ) / x * f1 - f0
    by(k) = f
    f0 = f1
    f1 = f
  end do

  do k = 2, nm
    dy(k) = by(k-1) - k * by(k) / x
  end do

  return
end
function msta1 ( x, mp )

!*****************************************************************************80
!
!! MSTA1 determines a backward recurrence starting point for Jn(x).
!
!  Discussion:
!
!    This procedure determines the starting point for backward  
!    recurrence such that the magnitude of    
!    Jn(x) at that point is about 10^(-MP).
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    08 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Input, integer ( kind = 4 ) MP, the negative logarithm of the 
!    desired magnitude.
!
!    Output, integer ( kind = 4 ) MSTA1, the starting point.
!
  implicit none

  real ( kind = 8 ) a0
  real ( kind = 8 ) envj
  real ( kind = 8 ) f
  real ( kind = 8 ) f0
  real ( kind = 8 ) f1
  integer ( kind = 4 ) it
  integer ( kind = 4 ) mp
  integer ( kind = 4 ) msta1
  integer ( kind = 4 ) n0
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) nn
  real ( kind = 8 ) x

  a0 = abs ( x )
  n0 = int ( 1.1D+00 * a0 ) + 1
  f0 = envj ( n0, a0 ) - mp
  n1 = n0 + 5
  f1 = envj ( n1, a0 ) - mp
  do it = 1, 20       
    nn = n1 - int ( real ( n1 - n0, kind = 8 ) / ( 1.0D+00 - f0 / f1 ) )               
    f = envj ( nn, a0 ) - mp
    if ( abs ( nn - n1 ) < 1 ) then
      exit
    end if
    n0 = n1
    f0 = f1
    n1 = nn
    f1 = f
  end do

  msta1 = nn

  return
end
function msta2 ( x, n, mp )

!*****************************************************************************80
!
!! MSTA2 determines a backward recurrence starting point for Jn(x).
!
!  Discussion:
!
!    This procedure determines the starting point for a backward
!    recurrence such that all Jn(x) has MP significant digits.
!
!    Jianming Jin supplied a modification to this code on 12 January 2016.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    14 January 2016
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of Jn(x).
!
!    Input, integer ( kind = 4 ) N, the order of Jn(x).
!
!    Input, integer ( kind = 4 ) MP, the number of significant digits.
!
!    Output, integer ( kind = 4 ) MSTA2, the starting point.
!
  implicit none

  real ( kind = 8 ) a0
  real ( kind = 8 ) ejn
  real ( kind = 8 ) envj
  real ( kind = 8 ) f
  real ( kind = 8 ) f0
  real ( kind = 8 ) f1
  real ( kind = 8 ) hmp
  integer ( kind = 4 ) it
  integer ( kind = 4 ) mp
  integer ( kind = 4 ) msta2
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n0
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) nn
  real ( kind = 8 ) obj
  real ( kind = 8 ) x

  a0 = abs ( x )
  hmp = 0.5D+00 * mp
  ejn = envj ( n, a0 )

  if ( ejn <= hmp ) then
    obj = mp
!
!  Original code:
!
!   n0 = int ( 1.1D+00 * a0 )
!
!  Updated code:
!
    n0 = int ( 1.1D+00 * a0 ) + 1
  else
    obj = hmp + ejn
    n0 = n
  end if

  f0 = envj ( n0, a0 ) - obj
  n1 = n0 + 5
  f1 = envj ( n1, a0 ) - obj

  do it = 1, 20
    nn = n1 - int ( real ( n1 - n0, kind = 8 ) / ( 1.0D+00 - f0 / f1 ) )
    f = envj ( nn, a0 ) - obj
    if ( abs ( nn - n1 ) < 1 ) then
      exit
    end if
    n0 = n1
    f0 = f1
    n1 = nn
    f1 = f
  end do

  msta2 = nn + 10

  return
end
function envj ( n, x )

!*****************************************************************************80
!
!! ENVJ is a utility function used by MSTA1 and MSTA2.
!
!  Discussion:
!
!    ENVJ estimates -log(Jn(x)) from the estimate
!    Jn(x) approx 1/sqrt(2*pi*n) * ( e*x/(2*n))^n
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    14 January 2016
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!    Modifications suggested by Vincent Lafage, 11 January 2016.
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the Bessel function.
!
!    Input, real ( kind = 8 ) X, the absolute value of the argument.
!
!    Output, real ( kind = 8 ) ENVJ, the value.
!
  implicit none

  real ( kind = 8 ) envj
  real ( kind = 8 ) logten
  integer ( kind = 4 ) n
  real ( kind = 8 ) n_r8
  real ( kind = 8 ) r8_gamma_log
  real ( kind = 8 ) x
!
!  Original code
!
  if ( .true. ) then

    envj = 0.5D+00 * log10 ( 6.28D+00 * n ) &
      - n * log10 ( 1.36D+00 * x / n )
!
!  Modification suggested by Vincent Lafage.
!
  else

    n_r8 = real ( n, kind = 8 )
    logten = log ( 10.0D+00 )
    envj = r8_gamma_log ( n_r8 + 1.0D+00 ) / logten - n_r8 * log10 ( x )

  end if

  return
end