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
