module darray
use constants
contains
   subroutine AddToList(list, element)
      IMPLICIT NONE
      integer :: i, isize
      real(long), intent(in) :: element
      real(long), dimension(:), allocatable, intent(inout) :: list
      real(long), dimension(:), allocatable :: clist
      if(allocated(list)) then
         isize = size(list)
         allocate(clist(isize+1))
         do i=1,isize          
            clist(i) = list(i)
         end do
         clist(isize+1) = element
         deallocate(list)
         call move_alloc(clist, list)
      else
         allocate(list(1))
         list(1) = element
      end if
   end subroutine AddToList

   subroutine Add_int_ToList(list, element)
      IMPLICIT NONE
      integer :: i, isize
      integer, intent(in) :: element
      integer, dimension(:), allocatable, intent(inout) :: list
      integer, dimension(:), allocatable :: clist
      if(allocated(list)) then
         isize = size(list)
         allocate(clist(isize+1))
         do i=1,isize          
            clist(i) = list(i)
         end do
         clist(isize+1) = element
         deallocate(list)
         call move_alloc(clist, list)
      else
         allocate(list(1))
         list(1) = element
      end if
   end subroutine Add_int_ToList

   subroutine Add_char_ToList(list, element)
      IMPLICIT NONE
      integer :: i, isize
      character(len = *), intent(in) :: element
      !character(len = short_str), intent(in) :: element
      character(len = short_str), dimension(:), allocatable, intent(inout) :: list
      character(len = short_str), dimension(:), allocatable :: clist
      if(allocated(list)) then
         isize = size(list)
         allocate(clist(isize+1))
         do i=1,isize          
            clist(i) = list(i)
         end do
         clist(isize+1) = element
         deallocate(list)
         call move_alloc(clist, list)
      else
         allocate(list(1))
         list(1) = element
      end if
   end subroutine Add_char_ToList

   subroutine Add_llchar_ToList(list, element)
      IMPLICIT NONE
      integer :: i, isize
      character(len = longlong_str), intent(in) :: element
      character(len = longlong_str), dimension(:), allocatable, intent(inout) :: list
      character(len = longlong_str), dimension(:), allocatable :: clist
      if(allocated(list)) then
         isize = size(list)
         allocate(clist(isize+1))
         do i=1,isize          
            clist(i) = list(i)
         end do
         clist(isize+1) = element
         deallocate(list)
         call move_alloc(clist, list)
      else
         allocate(list(1))
         list(1) = element
      end if
   end subroutine Add_llchar_ToList

   subroutine AddTo2D(mat, element)
      IMPLICIT NONE
      integer :: i, irow, icol
      real(long), dimension(:), intent(in) :: element
      real(long), dimension(:, :), allocatable, intent(inout) :: mat
      real(long), dimension(:, :), allocatable :: cmat
      if(allocated(mat)) then
         irow = size(mat, 1)
         icol = size(mat, 2)
         allocate(cmat(irow+1, icol))
         do i=1,irow          
            cmat(i, :) = mat(i, :)
         end do
         cmat(irow+1, :) = element(:)
         deallocate(mat)
         call move_alloc(cmat, mat)
      else
         allocate(mat(1, size(element)))
         mat(1, :) = element(:)
      end if
   end subroutine AddTo2D

   subroutine Add_int_list_To2D(mat, element)
      IMPLICIT NONE
      integer :: i, irow, icol
      integer, dimension(:), intent(in) :: element
      integer, dimension(:, :), allocatable, intent(inout) :: mat
      integer, dimension(:, :), allocatable :: cmat
      if(allocated(mat)) then
         irow = size(mat, 1)
         icol = size(mat, 2)
         allocate(cmat(irow+1, icol))
         do i=1,irow          
            cmat(i, :) = mat(i, :)
         end do
         cmat(irow+1, :) = element(:)
         deallocate(mat)
         call move_alloc(cmat, mat)
      else
         allocate(mat(1, size(element)))
         mat(1, :) = element(:)
      end if
   end subroutine Add_int_list_To2D
end module darray
