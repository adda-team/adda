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