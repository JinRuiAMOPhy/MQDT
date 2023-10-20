module interpolate
!----------------------------------------------------------------------------
! Interpolation Module
!   Ref) Press, Numerical Recipes, Chap. 3.
!
! Abstract data types
! - Linear : linear interpolation
!   . F Linear_create(x(N), y(N))
!   . F Linear_interpolate(L, x)
! - Spline : cubic spline interpolation
!   . F Spline_create(x(N), y(N), yp1, ypn) : yp1 and ypn are optional.
!   . F Spline_interpolate(S, x)
!
! Auxiliary functions and subroutines
! - F locate(table(N), x) : find the right place of x in table
!
! $Id: interpolate.f90 1541 2015-03-24 16:37:08Z ogeffert $
!----------------------------------------------------------------------------
  use constants
  implicit none
  private :: locate
  !
  ! linear interpolation
  !
  type :: Linear
    private
    integer :: n
    real(long), dimension(:), allocatable :: x, y
  end type
  !
  ! Spline interpolation
  !
  type :: Spline
    private
    integer :: n
    real(long), dimension(:), allocatable :: x, y, y2
  end type
  interface Spline_interpolate
    module procedure Spline_interpolate, Spline_interpolate_vector
  end interface
contains

function Linear_create( x, y ) result (L)
  real(long), dimension(:), intent(in) :: x, y
  type(Linear) :: L
  integer :: n
  n = size(x)
  if ( size(x) /= size(y) ) &
    stop "ERROR in Linear_create() : x, y must have the same size"
  L%n = n
  allocate( L%x(n), L%y(n) )
  L%x = x
  L%y = y
end function Linear_create

function Linear_interpolate( L, x ) result (y)
  type(Linear), intent(in) :: L
  real(long), intent(in) :: x
  integer :: j
  real(long) :: y, A, B, h
  j = locate( L%x, x )
  h = L%x(j+1) - L%x(j)
  A = ( L%x(j+1) - x ) / h
  B = ( x - L%x(j) ) / h
  y = A * L%y(j) + B * L%y(j+1)
end function Linear_interpolate

subroutine purge_Linear( L )
  type(Linear), intent(inout) :: L
  deallocate( L%x, L%y )
end subroutine purge_Linear

function Spline_create( x, y, yp1, ypn ) result (S)
  real(long), dimension(:), intent(in) :: x, y
  real(long), intent(in), optional :: yp1, ypn
  type(Spline) :: S
  real(long), dimension(size(x)) :: a, b, c, r, temp
  real(long) :: bet
  integer :: n, j
  n = size(x)
  if ( size(x) /= size(y) ) &
    stop "ERROR in Spline_create() : x, y must have the same size"
  S%n = n
  if ( .not. allocated(S%x) ) allocate( S%x(n) )
  if ( .not. allocated(S%y) ) allocate( S%y(n) )
  if ( .not. allocated(S%y2) ) allocate( S%y2(n) )
  S%x = x
  S%y = y
  c(1:n-1) = x(2:n) - x(1:n-1)          ! c(j) = x(j+1) - x(j)
  r(1:n-1) = 6.0_long * ( y(2:n) - y(1:n-1) ) / c(1:n-1)
  r(2:n-1) = r(2:n-1) - r(1:n-2)
  a(2:n-1) = c(1:n-2)
  b(2:n-1) = 2.0_long * ( c(2:n-1) + a(2:n-1) )
  b(1) = ONE
  b(n) = ONE
  !
  ! lower boundary
  !
  r(1) = ZERO
  c(1) = ZERO
  if ( present(yp1) ) then
    if ( yp1 < .99e30_long ) then
      r(1) = 3.0_long / (x(2)-x(1)) * ( (y(2)-y(1)) / (x(2)-x(1)) - yp1 )
      c(1) = 0.5_long
    end if
  end if
  !
  ! upper boundary
  !
  r(n) = ZERO
  a(n) = ZERO
  if ( present(ypn) ) then
    if ( ypn < .99e30_long ) then
      r(n) = -3.0_long / (x(n)-x(n-1)) * ( (y(n)-y(n-1)) / (x(n)-x(n-1)) - ypn )
      a(n) = 0.5_long
    end if
  end if
  !
  ! solve the tridiagonal matrix of a, b, c for r
  !
  bet = b(1)
  if ( b(1) == ZERO ) stop "TRD_solve: rewrite equations"
  S%y2(1) = r(1) / bet
  do j = 2, n           ! decomposition and forward substitution
    temp(j) = c(j-1) / bet
    bet = b(j) - a(j) * temp(j)
    if ( bet == ZERO ) stop "TRD_solve: failed"
    S%y2(j) = ( r(j) - a(j) * S%y2(j-1) ) / bet
  end do
  do j = n-1, 1, -1     ! backsubstitution
    S%y2(j) = S%y2(j) - temp(j+1) * S%y2(j+1)
  end do
end function Spline_create

subroutine Spline_reuse( S, y, yp1, ypn )
  type(Spline), intent(inout) :: S
  real(long), dimension(:), intent(in) :: y
  real(long), intent(in), optional :: yp1, ypn
  real(long), dimension(S%n) :: a, b, c, r, temp
  real(long) :: bet
  integer :: n, j
  n = S%n
  if ( n /= size(y(:)) ) &
    stop "ERROR in Spline_reuse() : x, y must have the same size"
  S%y(:) = y(:)
  c(1:n-1) = S%x(2:n) - S%x(1:n-1)          ! c(j) = x(j+1) - x(j)
  r(1:n-1) = 6.0_long * ( S%y(2:n) - S%y(1:n-1) ) / c(1:n-1)
  r(2:n-1) = r(2:n-1) - r(1:n-2)
  a(2:n-1) = c(1:n-2)
  b(2:n-1) = 2.0_long * ( c(2:n-1) + a(2:n-1) )
  b(1) = ONE
  b(n) = ONE
  !
  ! lower boundary
  !
  r(1) = ZERO
  c(1) = ZERO
  if ( present(yp1) ) then
    if ( yp1 < .99e30_long ) then
      r(1) = 3.0_long / (S%x(2)-S%x(1)) * &
             ( (S%y(2)-S%y(1)) / (S%x(2)-S%x(1)) - yp1 )
      c(1) = 0.5_long
    end if
  end if
  !
  ! upper boundary
  !
  r(n) = ZERO
  a(n) = ZERO
  if ( present(ypn) ) then
    if ( ypn < .99e30_long ) then
      r(n) = -3.0_long / (S%x(n)-S%x(n-1)) * &
             ( (S%y(n)-S%y(n-1)) / (S%x(n)-S%x(n-1)) - ypn )
      a(n) = 0.5_long
    end if
  end if
  !
  ! solve the tridiagonal matrix of a, b, c for r
  !
  bet = b(1)
  if ( b(1) == ZERO ) stop "TRD_solve: rewrite equations"
  S%y2(1) = r(1) / bet
  do j = 2, n           ! decomposition and forward substitution
    temp(j) = c(j-1) / bet
    bet = b(j) - a(j) * temp(j)
    if ( bet == ZERO ) stop "TRD_solve: failed"
    S%y2(j) = ( r(j) - a(j) * S%y2(j-1) ) / bet
  end do
  do j = n-1, 1, -1     ! backsubstitution
    S%y2(j) = S%y2(j) - temp(j+1) * S%y2(j+1)
  end do
end subroutine Spline_reuse

function Spline_interpolate( S, x ) result (y)
  type(Spline), intent(in) :: S
  real(long), intent(in) :: x
  integer :: j
  real(long) :: y, A, B, C, D, h
  if ( S%n == 1 ) then
    if (  S%x(1) == x ) then
      y = S%y(1)
      return 
    else
      print *,"ERROR in Spline_interpolate(): S%n == 1 at", x
      stop 
    endif
  endif
  j = locate( S%x, x )
  h = S%x(j+1) - S%x(j)
  A = ( S%x(j+1) - x ) / h
  B = ( x - S%x(j) ) / h
  C = ( A**3 - A ) * h**2 / 6.0_long
  D = ( B**3 - B ) * h**2 / 6.0_long
  y = A*S%y(j) + B*S%y(j+1) + C*S%y2(j) + D*S%y2(j+1)
  !
  ! to avoid extrapolation
  !
  if ( x < S%x(1) ) then
    y = S%y(1)
  else if ( x > S%x(S%n) ) then
    y = S%y(S%n)
  end if
end function Spline_interpolate

function Spline_interpolate_vector( S, xx ) result (yy)
  type(Spline), intent(in) :: S
  real(long), intent(in) :: xx(:)
  real(long) :: yy(size(xx))
  integer :: i
  !
  ! It is assumed that x(:) is sorted in ascending order
  !
  do i = 1, size(xx)
    yy(i) = Spline_interpolate( S, xx(i) )
  end do
end function Spline_interpolate_vector

function Spline_2nd_deriv( S, x ) result (y2)
  type(Spline), intent(in) :: S
  real(long), intent(in) :: x
  integer :: j
  real(long) :: y2, A, B, h
  j = locate( S%x, x )
  h = S%x(j+1) - S%x(j)
  A = ( S%x(j+1) - x ) / h
  B = ( x - S%x(j) ) / h
  y2 = A * S%y2(j) + B * S%y2(j+1)
end function Spline_2nd_deriv

subroutine purge_Spline( S )
  type(Spline), intent(inout) :: S
  deallocate( S%x, S%y, S%y2 )
end subroutine purge_Spline

function locate( table, x ) result (j)
!
! find the right place of x in the ordered table by means of bisection
! return j where table(j) < x < table(j+1) (j = 1, ..., n-1)
!   if table(1) == x, j = 1
!   if table(i) == x, j = i
!   but, if table(n) == x, j = n-1
!
  real(long), dimension(:), intent(in) :: table
  real(long), intent(in) :: x
  integer :: j, j_lo, j_hi
  j_lo = 1
  j_hi = size(table)
  do
    if (j_hi - j_lo <= 1) exit
    j = (j_hi + j_lo) / 2
    if ( table(j) > x ) then
      j_hi = j
    else
      j_lo = j
    end if
  end do
  j = j_lo
end function locate

end module interpolate
