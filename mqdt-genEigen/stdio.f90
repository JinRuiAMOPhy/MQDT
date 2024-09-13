module stdio
   use constants
   use numerical
   use darray
contains
! default if to_advance not present, do advance
! print with advane = 'no' in the end
   subroutine print_vector(V, leng, nw, fmtstr, funit, to_advance)
      implicit none 
      real(long), dimension(:) :: V
      integer, intent(in) :: funit, leng, nw
      integer :: i
      character(len = *), intent(in) :: fmtstr
      character(len = *), intent(in), optional :: to_advance
      logical :: do_advance 
      do_advance = (present(to_advance) &
           .and. (to_advance == 'yes' .or. to_advance == 'YES')) &
            .or.(.not.present(to_advance))
      do i = 1, leng
         write(funit, fmtstr, advance ='no')V(i)
         if((mod(i, nw) == 0 .or. i == leng) .and. do_advance) &
            write(funit, fmtstr, advance ='yes')
      end do 
   end subroutine print_vector

   subroutine print_vector_comment(V, leng, nw, fmtstr, funit, to_advance)
      implicit none 
      real(long), dimension(:) :: V
      integer, intent(in) :: funit, leng, nw
      integer :: i
      character(len = *), intent(in) :: fmtstr
      character(len = *), intent(in), optional :: to_advance
      logical :: do_advance 
      do_advance = (present(to_advance) &
           .and. (to_advance == 'yes' .or. to_advance == 'YES')) &
            .or.(.not.present(to_advance))

      do i = 1, leng
         if(mod(i, nw) == 1) write(funit, '("#")', advance ='no')
         write(funit, fmtstr, advance ='no')V(i)
         if((mod(i, nw) == 0 .or. i == leng) .and. do_advance) &
            write(funit, fmtstr, advance ='yes')
      end do 
   end subroutine print_vector_comment

   subroutine print_int_vector_comment(V, leng, nw, fmtstr, funit, to_advance)
      implicit none 
      integer, dimension(:) :: V
      integer, intent(in) :: funit, leng, nw
      integer :: i
      character(len = *), intent(in) :: fmtstr
      character(len = *), intent(in), optional :: to_advance
      logical :: do_advance 
      do_advance = (present(to_advance) &
           .and. (to_advance == 'yes' .or. to_advance == 'YES')) &
            .or.(.not.present(to_advance))

      do i = 1, leng
         if(mod(i, nw) == 1) write(funit, '("#")', advance ='no')
         write(funit, fmtstr, advance ='no')V(i)
         if((mod(i, nw) == 0 .or. i == leng) .and. do_advance) &
            write(funit, fmtstr, advance ='yes')
      end do 
   end subroutine print_int_vector_comment

   subroutine print_int_vector(V, leng, nw, fmtstr, funit, to_advance)
      implicit none 
      integer, dimension(:) :: V
      integer, intent(in) :: funit, leng, nw
      integer :: i
      character(len = *), intent(in) :: fmtstr
      character(len = *), intent(in), optional :: to_advance
      logical :: do_advance 
      do_advance = (present(to_advance) &
           .and. (to_advance == 'yes' .or. to_advance == 'YES')) &
            .or.(.not.present(to_advance))

      do i = 1, leng
         write(funit, fmtstr, advance ='no')V(i)
         if((mod(i, nw) == 0 .or. i == leng) .and. do_advance) &
            write(funit, fmtstr, advance ='yes')
      end do 
   end subroutine print_int_vector

   subroutine print_matrix_as_vect_by_row(V, nrow, ncol, nw, &
      fmtstr, funit)
      implicit none 
      real(long), dimension(:,:) :: V
      integer, intent(in) :: funit, nrow, ncol, nw
      integer :: i
      character(len = *), intent(in) :: fmtstr
      i = nrow
      if(nrow > 1) then
         do i = 1, nrow - 1
            call print_vector(V(i,:), ncol, nw, fmtstr, funit, 'no')
         end do 
         call print_vector(V(i,:), ncol, nw, fmtstr, funit, 'yes')
      else
         call print_vector(V(i,:), ncol, nw, fmtstr, funit, 'yes')
      end if
   end subroutine print_matrix_as_vect_by_row

   subroutine print_matrix_as_vect_by_col(V, nrow, ncol, nw, &
      fmtstr, funit)
      implicit none 
      real(long), dimension(:,:) :: V
      integer, intent(in) :: funit, nrow, ncol, nw
      integer :: i
      character(len = *), intent(in) :: fmtstr
      i = ncol 
      if(ncol > 1) then
         do i = 1, ncol - 1
            call print_vector(V(:,i), nrow, nw, fmtstr, funit, 'no')
         end do 
         call print_vector(V(:, i), nrow, nw, fmtstr, funit, 'yes')
      else
         call print_vector(V(:, i), nrow, nw, fmtstr, funit, 'yes')
      end if
   end subroutine print_matrix_as_vect_by_col

   subroutine print_matrix(V, nrow, ncol, nw, fmtstr, funit, to_advance)
      implicit none 
      real(long), dimension(:,:) :: V
      integer, intent(in) :: funit, nrow, ncol, nw
      integer :: i
      character(len = *), intent(in) :: fmtstr
      character(len = *), intent(in), optional :: to_advance
      do i = 1, nrow 
         if(present(to_advance)) then
             call print_vector(V(i,:), ncol, nw, fmtstr, funit, to_advance)
         else
             call print_vector(V(i,:), ncol, nw, fmtstr, funit)
         end if
      end do 
   end subroutine print_matrix

   subroutine print_matrix_comment(V, nrow, ncol, nw, fmtstr, funit, to_advance)
      implicit none 
      real(long), dimension(:,:) :: V
      integer, intent(in) :: funit, nrow, ncol, nw
      integer :: i
      character(len = *), intent(in) :: fmtstr
      character(len = *), intent(in), optional :: to_advance
      do i = 1, nrow 
         if(present(to_advance)) then
             call print_vector_comment(V(i,:), ncol, nw, fmtstr, funit, to_advance)
         else
             call print_vector_comment(V(i,:), ncol, nw, fmtstr, funit)
         end if
      end do 
   end subroutine print_matrix_comment

   subroutine print_int_matrix(V, nrow, ncol, nw, fmtstr, funit, to_advance)
      implicit none 
      integer, dimension(:,:) :: V
      integer, intent(in) :: funit, nrow, ncol, nw
      integer :: i
      character(len = *), intent(in) :: fmtstr
      character(len = *), intent(in), optional :: to_advance
      do i = 1, nrow 
         if(present(to_advance)) then
             call print_int_vector(V(i,:), ncol, nw, fmtstr, funit, to_advance)
         else
             call print_int_vector(V(i,:), ncol, nw, fmtstr, funit)
         end if
      end do 
   end subroutine print_int_matrix

   subroutine print_int_matrix_comment(V, nrow, ncol, nw, fmtstr, funit, to_advance)
      implicit none 
      integer, dimension(:,:) :: V
      integer, intent(in) :: funit, nrow, ncol, nw
      integer :: i
      character(len = *), intent(in) :: fmtstr
      character(len = *), intent(in), optional :: to_advance
      do i = 1, nrow 
         if(present(to_advance)) then
             call print_int_vector_comment(V(i,:), ncol, nw, fmtstr, funit, to_advance)
         else
             call print_int_vector_comment(V(i,:), ncol, nw, fmtstr, funit)
         end if
      end do 
   end subroutine print_int_matrix_comment

   subroutine read_key_val_array(funit, keys, vals)
      implicit none
      integer, intent(in) :: funit
      integer :: iost, isep, iend
      character(len = short_str), dimension(:), allocatable :: keys
      character(len = longlong_str), dimension(:), allocatable :: vals
      character(len = longlong_str) :: line, val
      character(len = short_str) :: key

      do 
         read(funit, '(A)', iostat = iost) line
         line = all_lower_case_str(line)
         if( iost > 0 ) then
            exit
         else if(iost < 0 ) then
            exit
         else 
            line = adjustl(line)
            if(line(1:1) == '#') cycle
            if(trim(line) == "") cycle
            isep = index(line, '=')
            key = line(1: isep - 1)
            call Add_char_ToList(keys, key)
            iend = index(line, '#') 
            if(iend == 0) iend = longlong_str
            val = line(isep + 1: iend - 1)
            if(index(val, '{') /= 0 ) then
               call Add_llchar_ToList(vals, read_wrapped_str_array(val, funit))
            else 
               call Add_llchar_ToList(vals, val)
            end if
         end if
      end do 
      return
   end subroutine read_key_val_array

   function read_wrapped_str_array(str_in, funit) result(res)
      implicit none
      character(len = longlong_str), intent(in) :: str_in
      character(len = longlong_str) :: res, line
      integer, intent(in) :: funit
      integer :: iend, roamL, roamR, iost
      line = adjustl(str_in)
      iend = index(line, '}')
      roamL = 1
      res = ''
      if(iend /= 0) then
         roamR = iend - 1 - 2 + 1
         res(roamL:roamR) = line(2: iend - 1)
      else
         iend = index(line, '#') - 1
         if(iend == -1) iend = len_trim(line)
         roamR = iend - 2 + 1
         res(roamL:roamR) = line(2: iend)
         do 
            roamL = roamR + 2
            read(funit, '(A)', iostat = iost) line
            line = adjustl(line)
            if(index(line, '#') /= 0) then
               iend = index(line, '#') - 1
               roamR = roamL + iend 
               res(roamL:roamR) = line(1: iend)
               exit
            else if(index(line, '}') /= 0) then
               iend = index(line, '}') - 1
               roamR = roamL + iend 
               res(roamL:roamR) = line(1: iend)
               exit
            else
               iend = len_trim(line)
               roamR = roamL + iend 
               res(roamL:roamR) = line(1: iend)
            end if
         end do 
      end if
   end function read_wrapped_str_array

   function from_line_get_ncol(line) result(ncol)
      implicit none
      integer :: ncol
      integer :: iost, i
      real(long) :: t
      character(len = longlong_str), intent(in) :: line
      ncol = 0
      t = 0
      do while(ncol <= max_col) 
         read(line, *, iostat = iost) (t, i = 1, ncol)
         if(iost /= 0) then
            exit
         else      
            ncol = ncol + 1
         end if
      end do 
      ncol = ncol - 1 ! !!!
      return
   end function from_line_get_ncol

   function from_line_get_str_ncol(line) result(ncol)
      implicit none
      integer :: ncol
      integer :: iost, i
      character(short_str) :: t
      character(len = longlong_str), intent(in) :: line
      ncol = 0
      t = ''
      do while(ncol <= max_col) 
         read(line, *, iostat = iost) (t, i = 1, ncol)
         if(iost /= 0) then
            exit
         else      
            ncol = ncol + 1
         end if
      end do 
      ncol = ncol - 1 ! !!!
      return
   end function from_line_get_str_ncol

   subroutine read_block_data(U, nrow, ncol, fn)
      implicit none
      real(long), dimension(:,:), allocatable, intent(out) :: U
      integer, intent(out) :: nrow, ncol
      integer, parameter :: funit = 1
      integer :: iost, irow, icol
      character(len = long_str), intent(in) :: fn
      character(len = longlong_str) :: line
      nrow = 0
      ncol = 0
      open(funit, file = fn, action = 'read')
      read(funit, '(A)',iostat = iost) line
      if(iost /= 0) return
      ncol = from_line_get_ncol(trim(line))
      rewind(funit)
      nrow = 0
      do
         read(funit, *,iostat = iost)
         if(iost > 0) then
            exit
         else if(iost < 0) then
            exit
         else      
            nrow = nrow + 1
         end if
      end do
      rewind(funit)
      allocate(U(nrow, ncol))
      do irow = 1, nrow
         read(funit, *,iostat = iost)(U(irow, icol), icol = 1, ncol)
         if(iost > 0) then
            exit
         else if(iost < 0) then
            exit
         end if
      end do
      close(funit)
      return 
   end subroutine read_block_data

   subroutine helper
      implicit none
      call comment_str_among_char("helper", '-', stdout, 80)
      write(*, "('solving mqdt equation in discrete/autoionization region.')")
      write(*, "(a)")
      write(*, "('Author:     Rui Jin')")
      write(*, "(a)")'Command line options:'
      call explain_something_colon(stdout, '--help --h', 'help', 20)
      call explain_something_colon(stdout, '-p -path -mqdt_path', &
           'path for mqdt results to solve, relative directory is also', 20)
      call explain_something_colon(stdout,'','accepted, default &
           &present working directory, if not specified',20)
      call explain_something_colon(stdout,'-debug -DEBUG',&
           'debug print controller, options:',20)
      call explain_something_colon(stdout,'',&
           'Fmat xgrid ygrid solu interp Tir Acoeff &
           &equationsolv para_use', 20)
      call explain_something_colon(stdout,'',&
           'one can combine options by e.g. -debug Fmat+xgrid',20)
      call explain_something_colon(stdout,'',&
           'or simply -debug all to print them all', 20)
      write(*, "(a)")'input control file options:'
      call comment_str_among_char("helper", '-', stdout, 80)
   end subroutine helper
   
   subroutine retrieve_xrange(xrange, ei, ef, n) ! nux_i, nux_f
      implicit none
      character(len = short_str), intent(in) :: xrange
      real(long), intent(out) :: ei, ef
      integer, intent(out), optional :: n
      integer :: isep, idiv, iend
      isep = index(xrange, ':')
      if(isep == 0) stop "range should be given with format begin:end [%number]"
      idiv = index(xrange, '%')
      n = 0
      if(idiv /= 0) then
         iend = idiv - 1
         read(xrange(idiv+1: short_str), *) n
      else 
         iend = short_str
      end if
      read(xrange(1:isep - 1), *) ei
      read(xrange(isep + 1: iend), *) ef
   end subroutine retrieve_xrange

   subroutine retrieve_Erange(arg, ei, ef) ! ei, ef in cm-1.
      implicit none
      character(len = short_str), intent(in) :: arg      
      character(len = short_str) :: Erange, Eunit
      real(long), intent(out) :: ei, ef
      integer :: isep
      call from_str_get_Erange_unit(arg, Erange, Eunit)
      isep = index(Erange, ':')
      read(Erange(1:isep - 1), *) ei
      ei = convert_energy_to_cm(ei, Eunit)
      read(Erange(isep + 1: short_str), *) ef
      ef = convert_energy_to_cm(ef, Eunit)
   end subroutine retrieve_Erange

   subroutine from_str_get_Erange_unit(arg, Erange, Eunit) 
      implicit none
      character(len = short_str), intent(in) :: arg
      character(len = short_str), intent(out) :: Erange, Eunit
      integer :: isep
      isep = index(arg, "_")
      if(isep /= 0) then
         Erange = arg(1:isep - 1)
         Eunit = arg(isep + 1:short_str)
         Eunit = all_lower_case_str(Eunit)
      else
         Erange = arg
         Eunit = 'ryd'
      end if
   end subroutine from_str_get_Erange_unit

   function ordinal(num) result(ord)
      implicit none
      integer, intent(in) :: num
      integer :: last_digit
      character(len = short_str) :: ord
      last_digit = mod(num, 10)
      if(last_digit == 1) then
         write(ord, "(i28, 'st')") num
      else if(last_digit == 2) then
         write(ord, "(i28, 'nd')") num
      else if(last_digit == 3) then
         write(ord, "(i28, 'rd')") num
      else
         write(ord, "(i28, 'th')") num
      end if
      ord = adjustl(ord)
   end function ordinal
   function check_item_in_list(keys, item) result(isthere)
      character(len = short_str), dimension(:), intent(in):: keys
      character(len = *), intent(in):: item
      logical :: isthere
      isthere = .false.
      if(sum(index(keys, item)) /= 0) isthere = .true.
   end function check_item_in_list

   subroutine comment_str_among_char(string, ch, funit, cmtleng) 
      character(len = *), intent(in):: string
      character :: ch
      integer, intent(in) :: funit
      integer, intent(in), optional :: cmtleng
      integer :: leng, before, after, i, cmtl
      integer, parameter :: text_width = 10
      character(len = short_str) :: tmp
      tmp = adjustl(string)
      leng = len_trim(tmp)
      if(.not.present(cmtleng)) then
         cmtl = short_str
      else
         cmtl = cmtleng
      end if
      before = int((cmtl - leng - 1) / 2)
      after = cmtl - before - 1 - leng
      write(funit, "('#')", advance = 'no') 
      write(funit, "(*(A))", advance = 'no') (ch, i = 1, before)
      write(funit, "(A)", advance = 'no') string
      write(funit, "(*(A))", advance = 'no') (ch, i = 1, after)
      write(funit, *)
   end subroutine comment_str_among_char

   subroutine explain_something_colon(funit, stringL, stringR, widthL)
      implicit none
      character(len = *), intent(in):: stringL, stringR
      integer, intent(in) :: funit
      integer, intent(in), optional :: widthL
      integer :: leng, cmtl
      integer, parameter :: text_width = 10
      character(len = long_str) :: tmp
      tmp = adjustl(stringL)
      leng = len_trim(tmp)
      if(.not.present(widthL)) then
         cmtl = short_str
      else
         cmtl = widthL
      end if
      write(funit, "(A, ':')", advance = 'no')tmp(1:cmtl)
      write(funit, "(A)") stringR
   end subroutine explain_something_colon

   function count_file_line(funit) result(nrow)
      implicit none
      integer, intent(in) :: funit
      integer :: nrow, iost
      nrow = 0
      do
         read(funit, *,iostat = iost)
         if(iost > 0) then
            exit
         else if(iost < 0) then
            exit
         else
            nrow = nrow + 1
         end if
      end do
      rewind(funit)
   end function count_file_line
   function all_lower_case_str(str) result(res)
      implicit none
      character(len = *), intent(in):: str
      character(len = len(str)):: res
      integer :: i
      res = str
      do i = 1, len(str)
         if(65<=ichar(str(i:i)) .and. &
            ichar(str(i:i)) <= 90) then
            res(i:i) = char(ichar(str(i:i))+97-65)
         end if
      end do
   end function all_lower_case_str
   function single_plural(word, n) result(word_out)
      implicit none
      character(len = *), intent(in) :: word
      integer, intent(in) :: n
      character(len = len(word) + 1) :: word_out
      if(n > 1) then
         write(word_out, "(A,'s')")word
      else
         write(word_out, "(A)")word
      end if
   end function single_plural
   function float2str(num, fmtstr) result(str)
      implicit none
      real(long), intent(in) :: num
      character(len = *), intent(in), optional :: fmtstr
      character(len = short_str) :: str
      if(present(fmtstr)) then
         write(str, fmtstr) num
      else
         write(str, *) num
      end if
      str = adjustl(str)
   end function float2str
   function int2str(num, fmtstr) result(str)
      implicit none
      integer, intent(in) :: num
      character(len = *), intent(in), optional :: fmtstr
      character(len = short_str) :: str
      if(present(fmtstr)) then
         write(str, fmtstr) num
      else
         write(str, *) num
      end if
      str = adjustl(str)
   end function int2str
end module stdio
