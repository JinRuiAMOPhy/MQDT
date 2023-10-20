module envset
   use constants
   implicit none
contains
   function get_NaN() result(nan)
     real(long) :: nan,xx
!     nan = IEEE_VALUE(nan, IEEE_QUIET_NAN)
     xx = -ONE
     nan = sqrt(xx)
   end function get_NaN

   function get_PI() result(PI)
     real(long) :: PI
     PI = asin(ONE) * TWO
   end function get_PI

   subroutine environmentset()
   ! it will be used later via constants module.
      nan = get_NaN() 
      PI = get_PI()
      AU_ACTION  = PLANCK_H_SI / TWO / PI
   end subroutine environmentset
end module envset
