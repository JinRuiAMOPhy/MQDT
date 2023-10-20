module file_cmdl_io
   use type_def
   use constants
   use darray
   use eigensort
   use dfde
   use numerical
contains
   subroutine from_str_get_dbg(CTR)
      implicit none
      type(control), intent(inout) :: CTR
      if(index(CTR%dbgstr, 'os') /= 0) CTR%dbg_os = .true.
      if(index(CTR%dbgstr, 'sort') /= 0) CTR%dbg_sort = .true.
      if(index(CTR%dbgstr, 'interp') /= 0) &
         CTR%dbg_interp = .true.
      if(index(CTR%dbgstr, 'levfind') /= 0) &
         CTR%dbg_levfind = .true.
      if(index(CTR%dbgstr, 'deriv') /= 0) &
         CTR%dbg_deriv = .true.
      if(index(CTR%dbgstr, 'dtaudE') /= 0) &
         CTR%dbg_deriv = .true.
      if(index(CTR%dbgstr, 'all') /= 0) then
         CTR%dbg_os = .true.
         CTR%dbg_sort = .true.
         CTR%dbg_interp = .true.
         CTR%dbg_levfind = .true.
         CTR%dbg_deriv = .true.
         CTR%dbg_deriv = .true.
      end if
      return
   end subroutine from_str_get_dbg

   subroutine from_entry_procedures(CTR)
      implicit none
      type(control), intent(inout) :: CTR
      CTR%yes_eigenchansort = .true. !.false.
      CTR%yes_dfde = .false.
      if(index(CTR%proc, 'dfde') /= 0) then
            CTR%yes_dfde = .true.
      else if(index(CTR%proc, 'eigenchansort') /= 0) then
            CTR%yes_eigenchansort = .true. 
      else if(index(CTR%proc, 'all') /=0 ) then
            CTR%yes_dfde = .true.
            CTR%yes_eigenchansort = .true.
      end if
   end subroutine from_entry_procedures

   function from_str_get_energy_unit(arg) result(E_ryd) 
      implicit none
      character(len = short_str), intent(in) :: arg
      real(long) :: E_ryd
      if(index(arg, "ryd") /= 0 .or. index(arg, "ry") /= 0) then
         read(arg(1:index(arg, "r") - 1), *) E_ryd
      else if(index(arg, "Ry") /= 0 .or. index(arg, "Ryd") /= 0) then
         read(arg(1:index(arg, "R") - 1), *) E_ryd
      else if(index(arg, "au") /= 0 .or. index(arg, "a.u.") /= 0) then
         read(arg(1:index(arg, "a") - 1), *) E_ryd
         E_ryd = TWO * E_ryd
      else if(index(arg, "Hartree") /= 0) then
         read(arg(1:index(arg, "H") - 1), *) E_ryd
         E_ryd = TWO * E_ryd
      else if(index(arg, "hartree") /= 0) then
         read(arg(1:index(arg, "h") - 1), *) E_ryd
         E_ryd = TWO * E_ryd
      else if(index(arg, "eV") /= 0 .or. index(arg, "ev") /= 0) then
         read(arg(1:index(arg, "e") - 1), *) E_ryd
         E_ryd = E_ryd * AU_E_POTENTIAL
      else if(index(arg, "cm-1") /= 0 ) then
         read(arg(1:index(arg, "e") - 1), *) E_ryd
         E_ryd = Energy_convert_cm_to_Ryd(E_ryd)
      else 
         print *, "****"
         print *, "****WARNING: -eg specified only energy value without unit, use default unit ryd "
         print *, "****"
         read(arg, *) E_ryd
      end if
   end function from_str_get_energy_unit
   subroutine read_args(CTR, args)
      implicit none
      type(control), intent(out) :: CTR
      integer :: i
      character(len = short_str) :: arg
      character(len = short_str), allocatable, dimension(:), intent(inout) :: args
      CTR%yes_ignore = .false.
      CTR%path_mqdt = '.'
      CTR%proc = ''
      CTR%figname = 'jhangz.png'
      CTR%couple = 'ls'
      !CTR%egnd_istate = nan
      CTR%egnd_istate = ZERO
      CTR%filter = ZERO ! filter out the fluctuation.
      i = 1

      if(COMMAND_ARGUMENT_COUNT() == 0) then 
         print *, "no arguments specified, use default parameters"
         call Add_char_ToList(args, "")
         return
      else if(COMMAND_ARGUMENT_COUNT() == 1) then
         call get_command_argument( 1, arg)
         arg = all_lower_case_str(arg)
         if( index(arg, '--help') /=0 .or. index(arg, '-h') /=0) &
            call helper
         stop
      end if

      do while( i <= COMMAND_ARGUMENT_COUNT() ) 
         call get_command_argument( i, arg)
         arg = trim(adjustl(arg))
         arg = all_lower_case_str(arg)
         call Add_char_ToList(args, arg(2:short_str))         
         if((trim(arg) .eq. '-path') .or. (trim(arg) .eq. '-mqdt_path') &
            .or. (trim(arg) .eq. '-p'))then
            i = i + 1
            call get_command_argument( i, arg)
            read(arg, '(A)' ) CTR%path_mqdt
         else if((trim(arg) .eq. '-filt') .or. trim(arg) .eq. '-filter')then
            i = i + 1
            call get_command_argument( i, arg)
            read(arg, * ) CTR%filter
            CTR%yes_user_filter = .true.
         else if((trim(arg) .eq. '-ground') .or. (trim(arg) .eq. '-eg'))then
            i = i + 1
            call get_command_argument( i, arg)
            CTR%egnd_istate = from_str_get_energy_unit(arg)
         else if(trim(arg) .eq. '-debug')then
            i = i + 1
            call get_command_argument( i, arg)
            read(arg, * ) CTR%dbgstr
            call from_str_get_dbg(CTR)
         else if((trim(arg) .eq. '-couple') .or. trim(arg) .eq. '-coupling')then
            i = i + 1
            call get_command_argument( i, arg)
            read(arg, * ) CTR%couple
         else if((trim(arg) .eq. '-ignore') .or. trim(arg) .eq. '-ig')then
            CTR%yes_ignore = .true.
         end if
         i = i + 1
      end do 
      call from_entry_procedures(CTR)
   end subroutine read_args

   subroutine run_dfde(FN, CTR, T, S, G)
      implicit none
      type(fname), intent(in) :: FN
      type(control), intent(inout) :: CTR
      type(mqdt), intent(inout) :: T
      type(grid), intent(in) :: G
      type(Smat), intent(inout) :: S
      logical :: file_ex
      inquire(file = FN%f_in_dmat, exist = file_ex) 
      if(CTR%yes_dfde .and. (.not.file_ex)) then
         stop 'intent to calculate dfde but without dalfa.in'
      else
         if(.not.CTR%yes_dfde) then 
            write(*, '(A)')"Don't plan to calculate dfde, but there is dalfa.in&
            & , So I do it any way."
            CTR%yes_dfde = .true.
         end if
         if(T%nop < T%nchan) then
         call dfn(T, S, FN, CTR)
         else 
         call dfde_open(T, S, FN, CTR, G)
         end if
      end if
      return 
   end subroutine run_dfde

   subroutine read_Smat_input(FN, S, CTR)
      implicit none
      type(fname), intent(in) :: FN
      type(Smat), intent(inout) :: S
      type(control), intent(inout) :: CTR
      integer, parameter:: funit = 15
      integer :: i, j, iost, nchan_t, kount
      real(long) :: E_t
      real(long), allocatable, dimension(:) :: u_tmp, mu_tmp, dm_tmp
      logical :: file_ex
      inquire(file = FN%f_in_Smat, exist = file_ex)
      if(.not. file_ex) then
         print *, 'Failed to open '//trim(FN%f_in_Smat)//', Please check using -path '
         stop
      end if
      open(funit, file = FN%f_in_Smat, action = 'read', form = 'formatted')
      read(funit, *,iostat = iost)S%nchan
      rewind(funit)
      S%n_ang = (S%nchan - 1) * S%nchan / 2 
      print *, S%nchan, S%n_ang
      S%n_es = 0
      allocate(mu_tmp(S%nchan), u_tmp(S%n_ang))
      do
         read(funit, *,iostat = iost)nchan_t, E_t, (mu_tmp(i), i = 1, S%nchan),&
             (u_tmp(j), j = 1, S%n_ang)
         if(nchan_t /= S%nchan) stop 'dimension wrong reading miuang.out'
         if(iost > 0) then
            write(*,*) trim(FN%f_in_Smat)//' wrong!'
            exit
         else if(iost < 0) then
            write(*,*)trim(FN%f_in_Smat), S%n_es
            exit
         else
            S%n_es = S%n_es + 1
            call AddToList(S%es, E_t)
            call AddTo2D(S%mu_s, mu_tmp)
            call AddTo2D(S%ang_s, u_tmp)
         end if
      end do 
      close(funit)
      deallocate(mu_tmp, u_tmp)
      inquire(file = FN%f_in_dmat, exist = file_ex)
      if(file_ex) then
         CTR%yes_dfde = .true.
         allocate(dm_tmp(S%nchan))
         open(funit, file = FN%f_in_dmat, action = 'read', form = 'formatted')
         kount = count_file_line(funit)
         if(kount == 1) then
            write(STDERR, "('Dalfa at one energy point is provided.')")
            write(funit_dbg, "('# **** WARNING ****')")
            write(funit_dbg, "('# Dalfa at one energy point is provided.')")
            write(funit_dbg, "('# use constants Dalfa')")
            do i = 1, S%n_es
               call AddTo2D(S%dmat_s, dm_tmp)
            end do 
         else if(kount /=1 .and. kount /= S%n_es) then
            stop 'None constant dalfa mismatch with miuang.out'
         else 
            do j = 1, S%n_es
               read(funit, *,iostat = iost)E_t, (dm_tmp(i), i = 1, S%nchan)
               if(iost > 0) then
                  write(*,*) trim(FN%f_in_dmat)//' wrong!'
                  exit
               else if(iost < 0) then
                  exit
               else
                  call AddTo2D(S%dmat_s, dm_tmp)
               end if
            end do
            print*, size(S%dmat_s, 1), size(S%dmat_s, 2)
         close(funit)
         deallocate(dm_tmp)      
         end if
      end if
   end subroutine read_Smat_input
   subroutine read_tau(FN, T)
      implicit none
      type(fname), intent(in) :: FN
      type(mqdt), intent(inout) :: T
      real(long) :: nu_tmp, E_tmp
      real :: start, finish
      integer, parameter :: funit = 100
      integer :: iost, i, ie
      character(len = *), parameter :: fmtstr = '"(3e20.10)"'
      logical :: file_exist
      T%n_e = 0
      inquire(file = FN%f_in_tau, exist = file_exist)
      if(.not. file_exist) then
         print *, 'Failed to open '//trim(FN%f_in_tau)//', Please check using -path '
         stop
      else if(file_exist) then
         write(*,"('Reading tau in: ')", advance = 'no')
         open(funit, file = FN%f_in_tau, action = 'read', form = 'formatted')
         call cpu_time(start)
         do 
            read( funit, *,iostat = iost)nu_tmp, E_tmp
            if(iost > 0) then
               write(*,*) trim(FN%f_in_tau)//' wrong!'
               exit
            else if(iost < 0) then
               write(*,"(A, i8)", advance = 'no') trim(FN%f_in_tau)//' Total:', T%n_e
               exit 
            else
               T%n_e = T%n_e + 1
            end if
         end do 
         rewind(funit)
         allocate(T%tau(T%n_e, T%nop), T%nux(T%n_e), T%E(T%n_e))
         do ie = 1, T%n_e
            read( funit, *,iostat = iost)T%nux(ie), T%E(ie), (T%tau(ie,i), i = 1, T%nop)
            if(iost > 0) then
               write(*,*) trim(FN%f_in_tau)//' wrong!'
               exit
            else if(iost < 0) then
               exit 
            end if
         end do
         close(funit)
         call cpu_time(finish)
         write(*,"(' takes ', E20.8, ' sec')")finish - start
      end if
      return 
   end subroutine read_tau

   subroutine read_tau_pushback_slow(FN, T)
      implicit none
      type(fname), intent(in) :: FN
      type(mqdt), intent(inout) :: T
      real(long), allocatable, dimension(:) :: tau_tmp
      real(long) :: nu_tmp, E_tmp
      real :: start, finish
      integer, parameter :: funit = 100
      integer :: iost, i
      character(len = *), parameter :: fmtstr = '"(3e20.10)"'
      logical :: file_exist
      inquire(file = FN%f_in_tau, exist = file_exist)
      if(.not. file_exist) then
         print *, 'Failed to open '//trim(FN%f_in_tau)//', Please check using -path '
         stop
      end if
      T%n_e = 0
      if(file_exist) then
         write(*,"('Reading tau in: ')", advance = 'no')
         open(funit, file = FN%f_in_tau, action = 'read', form = 'formatted')
         call cpu_time(start)
         allocate(tau_tmp(T%nop))
         do 
            read( funit, *,iostat = iost)nu_tmp, E_tmp, (tau_tmp(i), i = 1, T%nop)
            if(iost > 0) then
               write(*,*) trim(FN%f_in_tau)//' wrong!'
               exit
            else if(iost < 0) then
               write(*,"(A, i5)", advance = 'no') trim(FN%f_in_tau)//' Total:', T%n_e
               exit 
            else
               T%n_e = T%n_e + 1
               call AddToList(T%nux, nu_tmp)
               call AddToList(T%E, E_tmp)
               call AddTo2D(T%tau, tau_tmp)
            end if
         end do
         close(funit)
         deallocate(tau_tmp)
         call cpu_time(finish)
         write(*,"(' takes ', E20.8, ' sec')")finish - start
      end if
      return 
   end subroutine read_tau_pushback_slow

   subroutine read_A(FN, T, CTR)
      implicit none
      type(fname), intent(in) :: FN
      type(mqdt), intent(inout) :: T
      type(control), intent(in) :: CTR
      real(long) :: E_tmp
      real :: start, finish
      integer, parameter :: funit = 1, funit1 = 2, funit2 = 3, funit3 =  4
      integer :: iost, ie, ichan, io
      character(len = *), parameter :: fmtstr = '"(3e20.10)"'
      logical :: file_ex
      inquire(file = FN%f_in_Ara, exist = file_ex)
      if(.not. file_ex) then
         print*,'Failed to open '//trim(FN%f_in_Ara)//', Please check using -path '
         stop
      end if
      open(funit, file = FN%f_in_Ara, action = 'read', form = 'formatted')
      open(funit1, file = FN%f_in_Anorm, action = 'read', form = 'formatted')
      open(funit2, file = FN%f_in_Dn, action = 'read', form = 'formatted')
      if(CTR%yes_dis)then
         open(funit3, file = FN%f_in_normdos, action = 'read', form = 'formatted')
         allocate(T%norm_dos(T%n_e,T%nop))
      end if
      allocate(T%Ara(T%n_e,T%nchan, T%nop),  T%Anorm(T%n_e,T%nchan, T%nop), T%Dn(T%n_e,T%nop))
      call cpu_time(start)
      write(*,"('Reading A coefficients in: ',A,i8)", advance = 'no')trim(FN%f_in_Ara)//' Total:', T%n_e
      do ie = 1, T%n_e
         read( funit, *,iostat = iost)E_tmp,E_tmp, ((T%Ara(ie, ichan, io), ichan = 1, T%nchan), io = 1, T%nop)
         read( funit1, *,iostat = iost)E_tmp,E_tmp,((T%Anorm(ie, ichan, io), ichan = 1, T%nchan), io = 1, T%nop)
         read( funit2, *,iostat = iost)E_tmp,E_tmp,(T%Dn(ie, io), io = 1, T%nop)
         if(iost > 0) then
            write(*,*) trim(FN%f_in_Ara)//' wrong!'
            exit
         else if(iost < 0) then
            write(*,"(A,i5)", advance = 'no') trim(FN%f_in_Ara)//' Total:', ie
            exit 
         end if
      end do
      close(funit)
      close(funit1)
      close(funit2)
      if(CTR%yes_dis)then
         do ie = 1, T%n_e
            read( funit3, *,iostat = iost) &
               E_tmp,E_tmp,(T%norm_dos(ie, io), io = 1, T%nop)
            if(iost > 0) then
               write(*,*) trim(FN%f_in_Ara)//' wrong!'
               exit
            else if(iost < 0) then
               write(*,"(A,i5)", advance = 'no') trim(FN%f_in_normdos)//' Total:', ie
               exit 
            end if
         end do
         close(funit3)
      end if
      call cpu_time(finish)
      write(*,"(' takes ', E20.8, ' sec')")finish - start
      return 
   end subroutine read_A

   subroutine read_dfde(FN, CTR, T)
      implicit none
      type(fname), intent(in) :: FN
      type(mqdt), intent(inout) :: T
      type(control), intent(inout) :: CTR
      real(long) :: nu_tmp
      logical :: file_exist
      integer, parameter :: funit = 100
      integer :: iost, ie, io
      character(len = *), parameter :: fmtstr = '"(3e20.10)"'
      inquire(file = FN%f_in_dfde, exist = file_exist)
      print*, trim(FN%f_in_dfde), file_exist
      if(file_exist) then
         CTR%yes_dfde =.True.
         open(funit, file = FN%f_in_dfde, action = 'read', form = 'formatted')
         allocate(T%os(T%n_e,T%nop))
         do ie = 1, T%n_e
            read( funit, *,iostat = iost)nu_tmp,(T%os(ie, io), io = 1, T%nop)
            if(iost > 0) then
               write(*,*) trim(FN%f_in_dfde)//' wrong!'
               exit
            else if(iost < 0) then
               write(*,*) trim(FN%f_in_dfde)//' end:', ie
               exit 
            end if
         end do
         close(funit)
      else 
         if(CTR%yes_dfde) stop 'no dfde file'
         CTR%yes_dfde =.false.
      end if
      return 
   end subroutine read_dfde

   subroutine determin_ipx_ipy(T, G)
      implicit none
      type(mqdt), intent(inout) :: T
      type(grid), intent(inout) :: G
      integer :: i, j
      real(long) :: tmax
      if(G%E_cont_i> G%E_cont_f) stop 'exit program: Ei > Ef'
      tmax = T%IPs(1)
      do i = 2, T%n_IP 
         if(T%IPs(i) > tmax) then
            tmax = T%IPs(i)
         else if(T%IPs(i) == tmax) then
            stop 'exit program: repeated IP detected.'
         else 
            stop 'exit program: decreasing IPs.'
         end if
      end do 
      do i = 1, T%n_IP - 1
         if(G%E_cont_i + T%IPs(1) >= T%IPs(i) .and. &
            G%E_cont_i + T%IPs(1) < T%IPs(i+1)) exit
      end do 
      do j = 1, T%n_IP - 1
         if(G%E_cont_f + T%IPs(1) >= T%IPs(j) .and. &
            G%E_cont_f + T%IPs(1) < T%IPs(j+1)) exit
      end do 
      if(i /= j) stop 'exit program: Ei Ef cross ips'
      T%ipx = i + 1
      T%ipy = i
      return 
   end subroutine determin_ipx_ipy

   subroutine determin_nop(T, G)
      implicit none
      type(mqdt), intent(inout) :: T
      type(grid), intent(inout) :: G
      integer :: j
      if(G%E_cont_i> G%E_cont_f) stop 'exit program: Ei > Ef'
      T%nop = 0
      do j = 1, T%nchan 
         if(T%IPs(T%IPy) >= T%IPs(T%IP_seq(j))) then
            T%nop = T%nop + 1
            call Add_int_ToList(T%iopen, j)
         else
            call Add_int_ToList(T%iclose, j)
         end if
      end do 
      T%nclose = T%nchan - T%nop
      return 
   end subroutine determin_nop

   subroutine read_mqdt_input_key(FN, CTR, T, S, G, keys, args)
      implicit none
      type(fname), intent(in) :: FN
      type(control), intent(inout) :: CTR
      type(mqdt), intent(inout) :: T
      type(Smat), intent(inout) :: S
      type(grid), intent(inout) :: G
      integer, parameter :: funit = 15
      integer :: i, ncol, j
      logical :: file_ex
      character(len = short_str), dimension(:), allocatable, intent(inout) :: keys
      character(len = longlong_str), dimension(:), allocatable :: vals
      character(len = short_str), allocatable :: xrange(:)
      character(len = short_str), allocatable, dimension(:), intent(in) :: args
      inquire(file = FN%f_in_mqdt, exist = file_ex)
      if(.not. file_ex) then
         print *, 'Failed to open '//trim(FN%f_in_mqdt)// &
                  ', Please check using -path '
         print *, 'or  -p YOUR_MQDT_DIR'
         stop
      end if

      CTR%omega = 0.999
      G%ygrid_seg_jump = JUMP_THRESH
      G%xgrid_type = 'uniform'
      G%ygrid_type = 'uniform'
      G%ei = nan; G%ef = nan
      G%E_cont_i = nan; G%E_cont_f = nan
      G%nx_flat= 0; G%ny_adapt = 0
      G%nx_spike= 0;
      G%n_x_fine = 0; G%n_y_fine = 0
      G%ny_adapt = 20
      T%Z = nan
      !T%nchan = 0; S%nchan = 0
      !T%n_ang = 0; S%n_ang = 0
      T%nop = -1; T%n_IP = 0
      T%IPx = 0; T%IPy = 0
      T%k_mat_cofact = 1
      T%IPunit = ''
      CTR%discrete_method = 'jhangz' ! jhangz or level
      CTR%au_peak_search_method = 'tot' ! or 'chan'

      open(funit, file = FN%f_in_mqdt, action = 'read', &
           form = 'formatted')
      call read_key_val_array(funit, keys, vals)
      do i = 1, size(keys)
         if(check_item_in_list(args, trim(adjustl(keys(i))))) then
            print *, trim(adjustl(keys(i)))
            stop 'redefining parameters from cmd and mqdt.in'
         end if
         if(keys(i) == 'z') then
            read(vals(i), *) T%Z
         else if(keys(i) == 'nchan') then
            read(vals(i), *) ncol
            if(ncol /= T%nchan .and. T%nchan /= 0)then
               stop 'wrong nchan'
            else
               T%nchan = ncol
            end if
         else if(keys(i) == 'nop') then
            read(vals(i), *) T%nop
         else if(keys(i) == 'nip') then
            read(vals(i), *)  ncol
            if(ncol /= T%n_IP .and. T%n_IP /= 0)then
               stop 'wrong nip'
            else
               T%n_IP = ncol
            end if
         else if(keys(i) == 'k_mat_cofact') then
            read(vals(i), *) T%k_mat_cofact
         else if(keys(i) == 'x_grid_type') then
            read(vals(i), *) G%xgrid_type
         else if(keys(i) == 'y_grid_type') then
            read(vals(i), *) G%ygrid_type
         else if(keys(i) == 'ipy') then
            read(vals(i), *) T%IPy
         else if(keys(i) == 'ipx') then
            read(vals(i), *) T%IPx
         else if(keys(i) == 'twoj') then
            read(vals(i), *) T%twoJ
         else if(keys(i) == 'xrange') then
            call retrieve_xrange(vals(i), G%ei, G%ef, ncol)
         else if(keys(i) == 'ip') then
            ncol = from_line_get_ncol(vals(i))
            if(ncol /= T%n_IP .and. T%n_IP /= 0) then
               stop 'wrong nip'
            else
               T%n_IP = ncol
            end if
            allocate(T%IPs(ncol))
            read(vals(i), *)(T%IPs(j), j = 1, ncol)
            T%IPs(:) = T%IPs(:) - T%IPs(1)
            !print *, T%IPs(:)
         else if(keys(i) == 'ip_seq') then
            ncol = from_line_get_ncol(vals(i))
            if(ncol < T%nchan .and. T%nchan /= 0) then
               stop 'too few IP index or nchan'
            !else
            !   T%nchan = ncol
            end if
            allocate(T%IP_seq(T%nchan))
            read(vals(i), *)(T%IP_seq(j), j = 1, T%nchan)
         else if(keys(i) == 'nx_flat') then
            read(vals(i), *) G%nx_flat
         else if(keys(i) == 'nx_spike') then
            read(vals(i), *) G%nx_spike
         else if(keys(i) == 'ny_adapt') then
            read(vals(i), *) G%ny_adapt
         else if(keys(i) == 'ny_init_guess') then
            read(vals(i), *) G%ny_init_guess
         else if(keys(i) == 'ygrid_seg_jump') then
            read(vals(i), *) G%ygrid_seg_jump
         else if(keys(i) == 'debug') then
            read(vals(i), *) CTR%dbgstr
            call from_str_get_dbg(CTR)
         else if(keys(i) == 'dmu') then
            ncol = from_line_get_ncol(vals(i))
            if(ncol < T%nchan .and. T%nchan /= 0) then
               stop 'too few dmu'
            !else
            !   T%nchan = ncol
            end if
            allocate(S%d_mu(T%nchan))
            read(vals(i), *)(S%d_mu(j), j = 1, T%nchan)
         else if(keys(i) == 'x_fine') then
            ncol = from_line_get_str_ncol(vals(i))
            G%n_x_fine = ncol
            allocate(G%nxdiv(ncol), G%xdiv1(ncol), G%xdiv2(ncol), xrange(ncol))
            read(vals(i), *) (xrange(j), j = 1, ncol)
            do j = 1, ncol
               call retrieve_xrange(xrange(j), G%xdiv1(j), &
                                    G%xdiv2(j), G%nxdiv(j))
            end do 
            deallocate(xrange)
         else if(keys(i) == 'y_fine') then
            ncol = from_line_get_str_ncol(vals(i))
            G%n_y_fine = ncol
            allocate(G%nydiv(ncol), G%ydiv1(ncol), G%ydiv2(ncol), xrange(ncol))
            read(vals(i), *) (xrange(j), j = 1, ncol)
            do j = 1, ncol
               call retrieve_xrange(xrange(j), G%ydiv1(j), &
                                    G%ydiv2(j), G%nydiv(j))
            end do 
            deallocate(xrange)
         else if(keys(i) == 'relax') then
            read(vals(i), *) CTR%omega
         else if(keys(i) == 'eqnsolv_method') then
            read(vals(i), *) CTR%eqnsolv_method
         else if(keys(i) == 'e_continuum') then
            call retrieve_Erange(vals(i), G%E_cont_i, G%E_cont_f)
         else if(keys(i) == 'ip_unit') then
            read(vals(i), *) T%IPunit
         else if(keys(i) == 'discrete_method') then
            read(vals(i), *) CTR%discrete_method ! jhangz / level
         else if(keys(i) == 'au_peak_search_method') then
            read(vals(i), *) CTR%au_peak_search_method
         end if         
      end do 
      deallocate(vals)

      if(index(CTR%discrete_method, 'jhangz') /= 0 )then
         call check_para_jhangz(FN, CTR, T, S, G, keys)
      else if(index(CTR%discrete_method, 'level') /= 0 )then
         call check_para_level(FN, CTR, T, S, G, keys)
      else 
         call comment_str_among_char('INPUT ERROR', '*', stderr)
         write(stderr,*) 'An unacceptable discrete_method was &
         &provided, choose between jhangz/level (no specifying: default jhangz)'
         call comment_str_among_char('INPUT ERROR', '*', stderr)
         stop
      end if
   end subroutine read_mqdt_input_key

   subroutine check_para_jhangz(FN, CTR, T, S, G, keys)
      implicit none
      type(fname), intent(in) :: FN
      type(control), intent(inout) :: CTR
      type(mqdt), intent(inout) :: T
      type(Smat), intent(inout) :: S
      type(grid), intent(inout) :: G
      integer :: i, j
      character(len = short_str), dimension(:), allocatable, intent(inout) :: keys
      real(long) :: tmp1, tmp2
      
      if(check_NaN(T%Z)) stop 'Z not specified.'

      if(.not.check_item_in_list(keys,'eqnsolv_method')) then
         CTR%eqnsolv_method = 'bisect'
      else if(CTR%eqnsolv_method == 'newton' .and. &
         .not.check_item_in_list(keys,'relax')) then
         CTR%omega = 0.618
      else if(CTR%eqnsolv_method == 'hyb' .and. &
         .not.check_item_in_list(keys,'relax')) then
         CTR%omega = 0.618
      end if

      if(.not.allocated(T%IP_seq)) &
         stop 'IP index for each chann not specified.'
      if(.not.allocated(T%IPs)) &
         stop 'IPs not specified.'
      do j = 1, T%n_IP
         T%IPs(j) = convert_energy_to_cm(T%IPs(j), T%IPunit)
      end do 

      if(.not.check_item_in_list(keys, 'xrange') .and. &
         .not.check_item_in_list(keys, 'e_continuum')) then
         stop 'Neither nu range nor E_cont range specified.'
      else if(.not.check_item_in_list(keys, 'xrange')) then ! E_continuum specified but not xrange
         if(G%E_cont_f < ZERO .and. (.not.check_item_in_list(keys, 'ipx') &
            .or. .not.check_item_in_list(keys, 'ipy'))) then
           stop 'DIS detected from E_continuum, but no IPx / IPy specified'
         else if(G%E_cont_f < ZERO) then ! Valid DIS 
            CTR%yes_dis = .true.
            CTR%yes_au = .false.
         else if(G%E_cont_i >= ZERO) then ! AU
            CTR%yes_au = .true.
            CTR%yes_dis = .false.
            if(check_item_in_list(keys, 'ipy') .or. check_item_in_list(keys, 'ipx')) then 
               i = T%IPx; j = T%IPy
               call determin_IPx_IPy(T, G)
               if(i /= T%IPx .or. j /= T%IPy) then
                  write(STDERR,*)&
                  'IPx / IPy specified but mismatch from auto-generated &
                  &please check input'//trim(FN%f_in_mqdt)
                  stop
               end if
             else
               call determin_IPx_IPy(T, G)
               !write(stderr, "('IPx(auto)', i5, 'IPy(auto)', i5)")T%IPx, T%Ipy
             end if
         end if
         G%ei = from_E_cm_get_nu(T, G%E_cont_i + T%IPs(1), T%IPx)
         G%ef = from_E_cm_get_nu(T, G%E_cont_f + T%IPs(1), T%IPx)
         print *,'ei(cm-1)',G%E_cont_i, 'ef', G%E_cont_f
         print *,'Nu(ei)',G%ei, 'Nu(ef)', G%ef
      else if(.not.check_item_in_list(keys, 'e_continuum')) then !xrange specified but not E_continuum
         if(.not.check_item_in_list(keys, 'ipx') &
            .or. .not.check_item_in_list(keys, 'ipy')) then
            stop 'xrange specified but no IPx / IPy.'
         end if
         G%E_cont_i = get_E_cont_ryd(T, G%ei) * Rydberg
         G%E_cont_f = get_E_cont_ryd(T, G%ef) * Rydberg
         if(G%E_cont_f < ZERO) then ! Valid DIS
            CTR%yes_dis = .true.
            CTR%yes_au = .false.
         end if
         if(G%E_cont_i >= ZERO) then ! Valid AU
            CTR%yes_au = .true.
            CTR%yes_dis = .false.
         end if
      else
         if(.not.check_item_in_list(keys, 'ipx') &
            .or. .not.check_item_in_list(keys, 'ipy')) then
            if(G%E_cont_f < ZERO) then !DIS without  IPx
               stop 'xrange specified but no IPx / IPy.'
            else
               call determin_IPx_IPy(T, G)
            end if
         end if
         tmp1 = G%E_cont_i; tmp2 = G%E_cont_f
         G%E_cont_i = get_E_cont_ryd(T, G%ei) * Rydberg
         G%E_cont_f = get_E_cont_ryd(T, G%ef) * Rydberg
         if(tmp1 /= G%E_cont_i .or. tmp2 /= G%E_cont_f) then
            write(STDERR,*) 'conflicting xrange and E_continuum &
               &please check input'//trim(FN%f_in_mqdt)
               stop
         end if
      end if

      S%nchan = T%nchan
      T%n_ang = (T%nchan - 1) * T%nchan / 2 
      S%n_ang = T%n_ang
      if(.not.check_item_in_list(keys, 'nop')) then
         call determin_nop(T, G)
         if(T%nop <= 0) then
            write(STDERR,*)&
              'No nop specified and failed to generate due to &
              &wrong Erange/xrange input '//trim(FN%f_in_mqdt)
            stop 
         end if
      else 
         j = T%nop
         call determin_nop(T, G)
         if(T%nop /= j) then
            write(STDERR,*)&
            'nop specified but different from auto-generated one, please &
            &check input'//trim(FN%f_in_mqdt)
            stop
         end if
      end if
      if(maxval(T%IP_seq) > T%n_IP) stop 'More IP index then n_IP'
      if(size(T%IP_seq) /= T%nchan) stop 'Less IP index than nchan'
      if(.not.check_item_in_list(keys, 'dmu')) then
         allocate(S%d_mu(T%nchan))
         S%d_mu(:) = ZERO
      end if
   end subroutine check_para_jhangz

   subroutine check_para_level(FN, CTR, T, S, G, keys)
      implicit none
      type(fname), intent(in) :: FN
      type(control), intent(inout) :: CTR
      type(mqdt), intent(inout) :: T
      type(Smat), intent(inout) :: S
      type(grid), intent(inout) :: G
      integer :: i, j
      character(len = short_str), dimension(:), allocatable, intent(inout) :: keys
      real(long) :: tmp1, tmp2
      
      if(check_NaN(T%Z)) stop 'Z not specified.'

      if(.not.check_item_in_list(keys,'eqnsolv_method')) then
         CTR%eqnsolv_method = 'bisect'
      else if(CTR%eqnsolv_method == 'newton' .and. &
         .not.check_item_in_list(keys,'relax')) then
         CTR%omega = 0.618
      else if(CTR%eqnsolv_method == 'hyb' .and. &
         .not.check_item_in_list(keys,'relax')) then
         CTR%omega = 0.618
      end if

      if(.not.allocated(T%IP_seq)) &
         stop 'IP index for each chann not specified.'
      if(.not.allocated(T%IPs)) &
         stop 'IPs not specified.'
      do j = 1, T%n_IP
         T%IPs(j) = convert_energy_to_cm(T%IPs(j), T%IPunit)
      end do 

      if(.not.check_item_in_list(keys, 'xrange') .and. &
         .not.check_item_in_list(keys, 'e_continuum')) then
         stop 'Neither nu range nor E_cont range specified.'
      else if(.not.check_item_in_list(keys, 'xrange')) then ! E_continuum specified but not xrange
         if(G%E_cont_f < ZERO .and. (.not.check_item_in_list(keys, 'ipx') &
            .or. .not.check_item_in_list(keys, 'ipy'))) then
           stop 'DIS detected from E_continuum, but no IPx / IPy specified'
         else if(G%E_cont_f < ZERO) then ! Valid DIS 
            CTR%yes_dis = .true.
            CTR%yes_au = .false.
         else if(G%E_cont_i >= ZERO) then ! AU
            CTR%yes_au = .true.
            CTR%yes_dis = .false.
            if(check_item_in_list(keys, 'ipy') .or. check_item_in_list(keys, 'ipx')) then 
               i = T%IPx; j = T%IPy
               call determin_IPx_IPy(T, G)
               if(i /= T%IPx .or. j /= T%IPy) then
                  write(STDERR,*)&
                  'IPx / IPy specified but mismatch from auto-generated &
                  &please check input'//trim(FN%f_in_mqdt)
                  stop
               end if
             else
               call determin_IPx_IPy(T, G)
               !write(stderr, "('IPx(auto)', i5, 'IPy(auto)', i5)")T%IPx, T%Ipy
             end if
         end if
         G%ei = from_E_cm_get_nu(T, G%E_cont_i, T%IPx)
         G%ef = from_E_cm_get_nu(T, G%E_cont_f, T%IPx)
      else if(.not.check_item_in_list(keys, 'e_continuum')) then !xrange specified but not E_continuum
         if(.not.check_item_in_list(keys, 'ipx') &
            .or. .not.check_item_in_list(keys, 'ipy')) then
            stop 'xrange specified but no IPx / IPy.'
         end if
         G%E_cont_i = get_E_cont_ryd(T, G%ei) * Rydberg
         G%E_cont_f = get_E_cont_ryd(T, G%ef) * Rydberg
         if(G%E_cont_f < ZERO) then ! Valid DIS
            CTR%yes_dis = .true.
            CTR%yes_au = .false.
         end if
         if(G%E_cont_i >= ZERO) then ! Valid AU
            CTR%yes_au = .true.
            CTR%yes_dis = .false.
         end if
      else
         if(.not.check_item_in_list(keys, 'ipx') &
            .or. .not.check_item_in_list(keys, 'ipy')) then
            if(G%E_cont_f < ZERO) then !DIS without  IPx
               stop 'xrange specified but no IPx / IPy.'
            else
               call determin_IPx_IPy(T, G)
            end if
         end if
         tmp1 = G%E_cont_i; tmp2 = G%E_cont_f
         G%E_cont_i = get_E_cont_ryd(T, G%ei) * Rydberg
         G%E_cont_f = get_E_cont_ryd(T, G%ef) * Rydberg
         if(tmp1 /= G%E_cont_i .or. tmp2 /= G%E_cont_f) then
            write(STDERR,*) 'conflicting xrange and E_continuum &
               &please check input'//trim(FN%f_in_mqdt)
               stop
         end if
      end if

      S%nchan = T%nchan
      T%n_ang = (T%nchan - 1) * T%nchan / 2 
      S%n_ang = T%n_ang
      if(.not.check_item_in_list(keys, 'nop')) then
         call determin_nop(T, G)
         if(T%nop <= 0) then
            write(STDERR,*)&
              'No nop specified and failed to generate due to &
              &wrong Erange/xrange input'//trim(FN%f_in_mqdt)
            stop 
         end if
      else 
         j = T%nop
         call determin_nop(T, G)
         if(T%nop /= j) then
            write(STDERR,*)&
            'nop specified but different from auto-generated one, please &
            &check input'//trim(FN%f_in_mqdt)
            stop
         end if
      end if
      if(maxval(T%IP_seq) > T%n_IP) stop 'More IP index then n_IP'
      if(size(T%IP_seq) /= T%nchan) stop 'Less IP index than nchan'
      if(.not.check_item_in_list(keys, 'dmu')) then
         allocate(S%d_mu(T%nchan))
         S%d_mu(:) = ZERO
      end if
   end subroutine check_para_level

   subroutine initiate(FN, T, S, CTR, G)
      implicit none
      type(mqdt), intent(out) :: T
      type(fname), intent(out) :: FN
      type(control), intent(inout) :: CTR
      type(Smat), intent(inout) :: S
       type(grid), intent(out) :: G
      character(len = short_str), dimension(:), allocatable :: keys
      real(long) :: start, finish
      character(len = short_str) :: sfmtstr 
      integer :: i, ie, npoint
      character(len = short_str), allocatable, dimension(:) :: args
      integer, parameter :: funit = 1
      real(long) :: nux
      CTR%path_mqdt = '.'
      CTR%couple = 'ls'
      CTR%yes_dis = .false.
      CTR%yes_au = .false.; CTR%yes_dfde = .false.
      CTR%yes_user_filter = .false.
      CTR%omega = nan

      CTR%dbg_os = .false.; CTR%dbg_sort = .false.
      CTR%dbg_interp = .false.
      CTR%dbg_levfind = .false.
      CTR%dbgstr = ''; CTR%dbg_deriv = .false.


      call read_args(CTR, args)

      FN%f_in_tau = trim(CTR%path_mqdt)//"/tau.out"
      FN%f_in_dfde = trim(CTR%path_mqdt)//"/dfde.out"
      FN%f_exp_tau = trim(CTR%path_mqdt)//"/elev.exp"
      FN%f_exp_dfde = trim(CTR%path_mqdt)//"/oscistr.exp"
      FN%f_in_Ara = trim(CTR%path_mqdt)//"/Ara.out"
      FN%f_in_Anorm = trim(CTR%path_mqdt)//"/Anorm.out"
      FN%f_in_Dn = trim(CTR%path_mqdt)//"/Dn.out"
      FN%f_o_os = trim(CTR%path_mqdt)//"/os.out"
      FN%f_in_normdos = trim(CTR%path_mqdt)//"/NormDos.out"
      FN%f_o_os_dos = trim(CTR%path_mqdt)//"/OS_normd_dos.out"
      FN%f_o_sig = trim(CTR%path_mqdt)//"/sig.out"
      FN%f_o_tau_sort = trim(CTR%path_mqdt)//"/tau2-plot.out"
      FN%f_o_os_sort = trim(CTR%path_mqdt)//"/dfde2-plot.out"
      FN%f_o_energyrel = trim(CTR%path_mqdt)//"/nux_nuy.out"
      FN%f_o_phys = trim(CTR%path_mqdt)//"/phystable.out"
      FN%f_o_phys_compare = trim(CTR%path_mqdt)//"/theo-exp-compare.out"
      FN%f_o_exp_plot = trim(CTR%path_mqdt)//"/exp-plot.out"
      FN%f_o_theo_plot = trim(CTR%path_mqdt)//"/theo-plot.out"
      FN%f_in_dmat = trim(CTR%path_mqdt)//"/dalfa.in"
      FN%f_dbg = trim(CTR%path_mqdt)//"/debug-post.out"
      open(funit_dbg, file = FN%f_dbg, action = "write")

      FN%f_in_mqdt = trim(CTR%path_mqdt)//"/mqdt.in"

      FN%f_in_Smat = trim(CTR%path_mqdt)//"/miuang.out"
      call read_Smat_input(FN, S, CTR)
      T%nchan = S%nchan
      T%n_ang = (T%nchan - 1) * T%nchan / 2 
      S%n_ang = T%n_ang
      call read_mqdt_input_key(FN, CTR, T, S, G, keys, args)
      do ie =1, S%n_es
        S%mu_s(ie,:) = S%mu_s(ie,:) + S%d_mu(:)
      end do 


      write(funit_dbg, "('#----- Parameters -------- ')") 
      write(funit_dbg, "('# Z', i3)") int(T%Z)
      write(funit_dbg, "('# nchan', i3)") T%nchan
      write(funit_dbg, "('# IPx', i3)", advance = 'no') T%IPx
      if(.not.check_item_in_list(keys,'ipx')) then
         write(funit_dbg, "(' (internal)')")
      else
         write(funit_dbg, *)
      end if
      write(funit_dbg, "('# IPy', i3)", advance = 'no') T%IPy
      if(.not.check_item_in_list(keys,'ipy')) then
         write(funit_dbg, "(' (internal)')")
      else
         write(funit_dbg, *)
      end if
      write(funit_dbg, "('# E range in terms of nu (IPx =',i3,')',&
            & f8.4, ':', f8.4)", advance = 'no') T%IPx, G%ei, G%ef
      if(.not.check_item_in_list(keys,'xrange')) then
         write(funit_dbg, "(' (internal)')")
      else
         write(funit_dbg, *)
      end if
      write(funit_dbg, "('# Continuum E (Ry)', f8.4, ':', f8.4)",&
         advance = 'no')G%E_cont_i / Rydberg, G%E_cont_f / Rydberg
      if(.not.check_item_in_list(keys,'e_continuum')) then
         write(funit_dbg, "(' (internal)')")
      else
         write(funit_dbg, *)
      end if
      write(funit_dbg, "('# nop', i3)", advance = 'no') T%nop
      if(.not.check_item_in_list(keys,'nop')) then
         write(funit_dbg, "(' (internal)')")
      else
         write(funit_dbg, *)
      end if
      write(funit_dbg, "('# nip: ', i3)", advance = 'no') T%n_IP
      if(.not.check_item_in_list(keys,'nip')) then
         write(funit_dbg, "(' (internal)')")
      else
         write(funit_dbg, *)
      end if
      write(funit_dbg, "('# IP(in cm-1):')", advance = 'yes')
      sfmtstr = '(E20.10)'
      call print_vector_comment(T%IPs, T%n_IP, print_col_l_fmt, &
           sfmtstr, funit_dbg)
      write(funit_dbg, "('# IP index for each chann:')", advance = 'yes')
      sfmtstr = '(i4)'
      call print_int_vector_comment(T%IP_seq, T%nchan, print_col_s_fmt, &
           sfmtstr, funit_dbg)
      write(funit_dbg, "('# x refine', i3)") G%n_x_fine
      do i = 1, G%n_x_fine
         write(funit_dbg, "('# ')", advance = 'no')
         write(funit_dbg, "(2f8.4, i8)")G%xdiv1(i), G%xdiv2(i), G%nxdiv(i)
      end do
      write(funit_dbg, "('# x refine', i3)") G%n_x_fine
      do i = 1, G%n_y_fine
         write(funit_dbg, "('# ')", advance = 'no')
         write(funit_dbg, "(2f8.4, i8)")G%ydiv1(i), G%ydiv2(i), G%nydiv(i)
      end do 
      write(funit_dbg, "('# Smatrix file:', A, ' at ', i3, ' energies')")&
         trim(FN%f_in_Smat), S%n_es

      if(.not.check_item_in_list(keys, 'dmu')) then
         write(funit_dbg, "('# Smatrix correction: ')", advance = 'no')
         write(funit_dbg, "(' [no/off]')")
      else 
         write(funit_dbg, "('# Smatrix correction: ')", advance = 'no')
         sfmtstr = '(f8.4)'
         call print_vector_comment(S%d_mu, T%nchan, print_col_s_fmt, &
              sfmtstr, funit_dbg)
         write(funit_dbg, *)      
      end if
      write(funit_dbg, "('# energy of initial symmetry &
         &state:', E20.10, ' Ry')", advance = 'no')CTR%egnd_istate
      if(.not.check_item_in_list(args, 'eg')) &
         write(funit_dbg, "(' (internal/cmdl)')")

      write(funit_dbg, "('#----- Parameters -------- ')") 
      deallocate(keys)
      deallocate(args)
      if(T%nop < T%nchan) then
        call read_tau(FN, T)
        call read_A(FN, T, CTR)
      else 
        T%n_e = S%n_es
        T%n_e_sort = S%n_es
        allocate(T%tau(T%n_e, T%nop), T%nux(T%n_e), T%E(T%n_e))
        allocate(T%nux_sort(T%n_e_sort), T%tau2(T%n_e_sort,T%nchan),T%os2(T%n_e_sort,T%nchan))
        allocate(T%Ara(T%n_e,T%nchan, T%nop),  T%Anorm(T%n_e,T%nchan, T%nop), T%Dn(T%n_e,T%nop))
        T%Ara = ZERO; T%Anorm = ZERO; T%Dn = ONE
        do i = 1, T%nchan
          T%Ara(:, i, i) = ONE; T%Anorm(:, i,i) = ONE
        end do 
      end if
      print *, "Tau A read<<<<"
      if(CTR%yes_dfde) then
         write(*,"('creating dfde')")
         call cpu_time(start)
         call run_dfde(FN, CTR, T, S, G)
         call cpu_time(finish)
         write(*,"('dfde done, taking ', E20.10, ' sec')") finish - start
         write(*,"('Sort tau and dfde eigenchannel')")
      else if (T%nop < T%nchan) then
         write(*,"('Just sort tau eigenchannel')")
      end if
      call cpu_time(finish)
      if(T%nop < T%nchan) then
        call sort_mqdt(T, CTR, FN)
        write(*,"('sorted: inserted ', i3, ' NaNs to seperate &
           &unwanted connections, taking',E20.10, ' sec')")&
           T%n_e_sort - T%n_e, finish - start
        if(CTR%dbg_deriv .and. index(CTR%au_peak_search_method, 'tot') /=0 ) then
           FN%f_dbg_dtdE = trim(CTR%path_mqdt)//"/dbg_sum_dtaudE.out"
        else if(CTR%dbg_deriv .and. index(CTR%au_peak_search_method, 'chan') /=0 ) then
           FN%f_dbg_dtdE = trim(CTR%path_mqdt)//"/dbg_dtaudE.out"
        end if
      else ! all open channel, nothing to project.
        open(funit, file=FN%f_o_tau_sort, action = 'write')
        sfmtstr = "(E20.10)"
        print *, G%E_cont_i /T%Z**2/ Rydberg,'ei',S%es(1),G%E_cont_f/T%Z**2/ &
         Rydberg, S%es(S%n_es)
        npoint = 0
        do ie = 1, S%n_es
           if((G%E_cont_i /T%Z**2/ Rydberg > S%es(ie)) .or. &
              (S%es(ie) > G%E_cont_f/T%Z**2/ Rydberg)) cycle 
           !print *, 'in range', ie, S%es(ie)
           npoint = npoint + 1
           T%tau(npoint,:) = S%mu_s(ie, :)
           T%tau2(npoint,:) = S%mu_s(ie, :)
           !T%U(:,:) = S%U_s(:,:)
           nux = from_E_cm_get_nu(T, S%es(ie) * Rydberg * T%Z**2 + T%IPs(1), T%IPx)
           T%nux(npoint) = nux
           T%nux_sort(npoint) = nux
           write(funit,"(E20.10)", advance = 'no') nux
           call print_vector(S%mu_s(ie,:), T%nchan, T%nchan, sfmtstr, funit)
        end do 
        close(funit)
        T%n_e = npoint
        T%n_e_sort = npoint
      end if
      return 
   end subroutine initiate

   subroutine finalize(T, S, phys, spec)
      type(mqdt), intent(inout) :: T
      type(Smat), intent(inout) :: S
      type(physobserv), intent(inout) :: phys
      type(spectrodata), intent(inout) :: spec
      if(allocated(S%d_mu)) deallocate(S%d_mu)
      if(allocated(T%IP_seq)) deallocate(T%IP_seq)
      if(allocated(T%IPs)) deallocate(T%IPs)
      if(allocated(T%tau)) deallocate(T%tau)
      if(allocated(T%nux)) deallocate(T%nux)
      if(allocated(T%E)) deallocate(T%E)
      if(allocated(T%Ara)) deallocate(T%Ara)
      if(allocated(T%Anorm)) deallocate(T%Anorm)
      if(allocated(T%Dn)) deallocate(T%Dn)
      if(allocated(T%os)) deallocate(T%os)
      if(allocated(T%os2)) deallocate(T%os2)
      if(allocated(T%tau2)) deallocate(T%tau2)
      if(allocated(T%nux_sort)) deallocate(T%nux_sort)
      if(allocated(T%os_dos)) deallocate(T%os_dos)
      if(allocated(S%mu_s)) deallocate(S%mu_s)
      if(allocated(S%dmat_s)) deallocate(S%dmat_s)
      if(allocated(S%ang_s)) deallocate(S%ang_s)
      if(allocated(S%es)) deallocate(S%es)
      if(allocated(S%spec_jj)) deallocate(S%spec_jj)
      if(allocated(S%spec_ls)) deallocate(S%spec_ls)
      if(allocated(phys%lev_nux)) deallocate(phys%lev_nux)
      if(allocated(phys%lev_nuy)) deallocate(phys%lev_nuy)
      if(allocated(phys%lev_ryd)) deallocate(phys%lev_ryd)
      if(allocated(phys%width)) deallocate(phys%width)
      if(allocated(phys%nux_plot)) deallocate(phys%nux_plot)
      if(allocated(phys%nuy_plot)) deallocate(phys%nuy_plot)
      if(allocated(phys%eigenchan_id)) deallocate(phys%eigenchan_id)
      if(allocated(phys%lev_indx)) deallocate(phys%lev_indx)
      if(allocated(phys%nlev_per_chan)) deallocate(phys%nlev_per_chan)
      if(allocated(phys%elev_spec)) deallocate(phys%elev_spec)
      if(allocated(spec%lev_nux)) deallocate(spec%lev_nux)
      if(allocated(spec%lev_nuy)) deallocate(spec%lev_nuy)
      if(allocated(spec%energy)) deallocate(spec%energy)
      if(allocated(spec%lev_indx)) deallocate(spec%lev_indx)
      if(allocated(spec%width)) deallocate(spec%width)
      if(allocated(spec%eigenchan_id)) deallocate(spec%eigenchan_id)
      if(allocated(spec%elev_spec)) deallocate(spec%elev_spec)
      close(funit_dbg)
   end subroutine finalize
end module file_cmdl_io
