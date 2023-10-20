program main
   implicit none
   real(8) :: a
   real(8), dimension(3) :: b
   integer :: i 
   do i = 0, 11
      a = -2.0 + dble(i) * 0.4
      print '("a ", f8.4, " a( mod 1) ", f8.4 )', a, mod(a, 1.0)
   end do 
   b(1) = 0.1
   b(2) = -0.1+1.0
   b(3) = -0.2+1.0
   print*,(b(i), i = 1, 3)
   b = smooth_Func(b)
   print*,(b(i), i = 1, 3)
contains 
   function smooth_Func(fun) result(f2)
      implicit none
      integer :: fun_end, i
      integer, parameter :: long = 8
      real(8), parameter ::  ONE = 1.d0
      real(8), dimension(:), intent(in) :: fun
      real(8), dimension(size(fun)) :: f2
      real(8) :: dfun
      f2(:) = fun(:)
      fun_end = size(fun)
      do i = 1, fun_end - 1
         dfun = f2(i) - f2(i+1)
         if(dfun > 0.3) then
            f2(i+1) = fun(i+1) + ONE
         else if(dfun < -0.3) then
            f2(i+1) = fun(i+1) - ONE
         end if
      end do
      return
   end function smooth_Func
end program
