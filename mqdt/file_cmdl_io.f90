module file_cmdl_io
   use type_def
   use constants
   use darray
   use stdio
   use numerical
contains
   subroutine from_str_get_dbg(CTR)
      implicit none
      type(control), intent(inout) :: CTR
      if(index(CTR%dbgstr, 'Fmat') /= 0) then
         CTR%dbg_F = .true.
      end if
      if(index(CTR%dbgstr, 'xgrid') /= 0) then
         CTR%dbg_xgrid = .true.
      end if
      if(index(CTR%dbgstr, 'ygrid') /= 0) then
         CTR%dbg_ygrid = .true.
      end if
      if(index(CTR%dbgstr, 'solu') /= 0) then
         CTR%dbg_solu = .true.
      end if
      if(index(CTR%dbgstr, 'interp') /= 0) then
         CTR%dbg_interp = .true.
      end if
      if(index(CTR%dbgstr, 'Tir') /= 0) then
         CTR%dbg_Tir = .true.
      end if
      if(index(CTR%dbgstr, 'Acoeff') /= 0) then
         CTR%dbg_A = .true.
      end if
      if(index(CTR%dbgstr, 'equationsolv') /= 0) then
         CTR%dbg_eqnsolv = .true.
      end if
      if(index(CTR%dbgstr, 'para_use') /= 0) then
         CTR%dbg_para_use = .true.
      end if

      if(index(CTR%dbgstr, 'all') /= 0) then
         CTR%dbg_F = .true.
         CTR%dbg_xgrid = .true.
         CTR%dbg_ygrid = .true.
         CTR%dbg_solu = .true.
         CTR%dbg_interp = .true.
         CTR%dbg_Tir = .true.
         CTR%dbg_A = .true.
         CTR%dbg_eqnsolv = .true.
         CTR%dbg_para_use = .true.
      end if
      return
   end subroutine from_str_get_dbg

   subroutine read_args(CTR, args)
      implicit none
      type(control), intent(out) :: CTR
      integer :: i
      character(len = short_str) :: arg
      character(len = short_str), allocatable, dimension(:), intent(inout) :: args
      i = 1
      if(COMMAND_ARGUMENT_COUNT() == 0) then
         print *, "no arguments specified, use default parameters"
         call Add_char_ToList(args, "")
         return
      else if(COMMAND_ARGUMENT_COUNT() == 1) then
         call get_command_argument( 1, arg)
         arg = all_lower_case_str(arg)
         arg = trim(adjustl(arg))
         call Add_char_ToList(args, arg(2:short_str))
         if( index(arg, '--help') /=0 .or. index(arg, '--h') /=0) &
            call helper
         stop
      end if

      do while( i <= COMMAND_ARGUMENT_COUNT() ) 
         call get_command_argument( i, arg)
         arg = all_lower_case_str(arg)
         call Add_char_ToList(args, arg)
         if((trim(arg) .eq. '-path') .or. (trim(arg) .eq. '-mqdt_path') &
            .or. (trim(arg) .eq. '-p'))then
            i = i + 1
            call get_command_argument( i, arg)
            read(arg, '(A)' ) CTR%path_mqdt
         else if(trim(arg) .eq. '-debug')then
            i = i + 1
            call get_command_argument( i, arg)
            read(arg, * ) CTR%dbgstr
            call from_str_get_dbg(CTR)
         end if
         i = i + 1
      end do 
   end subroutine read_args

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
      if(i /= j) then
        write(stderr, *) 'exit program: Ei Ef cross ips', T%IPs(i), &
          T%IPs(j)/Rydberg
        stop 
      end if
      T%ipx = i + 1
      T%ipy = i
      print *, 'ipx,ipy', T%ipx, T%ipy, T%nchan
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
      inquire(file = FN%f_mqdt, exist = file_ex)
      if(.not. file_ex) then
         print *, 'Failed to open '//trim(FN%f_mqdt)// &
                  ', Please check using -path '
         print *, 'or  -p YOUR_MQDT_DIR'
         stop
      end if

      CTR%omega = 0.999_long
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
      !T%nchan = 0; S%nchan = 0 ! read from miuang.out
      !T%n_ang = 0; S%n_ang = 0
      T%nop = -1; T%n_IP = 0
      T%IPx = 0; T%IPy = 0
      T%k_mat_cofact = 1
      T%IPunit = ''
      CTR%discrete_method = 'jhangz' ! jhangz or level

      open(funit, file = FN%f_mqdt, action = 'read', &
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
            if(ncol < T%nchan .and. T%nchan /= 0)then
               stop 'wrong nchan'
            !else
               !T%nchan = ncol
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
         else if(keys(i) == 'ip_seq') then
            ncol = from_line_get_ncol(vals(i))
            if(ncol < T%nchan .and. T%nchan /= 0) then
               stop 'wrong IP index or nchan'
            else if(ncol > T%nchan .and. T%nchan /= 0) then
               write(*,*)'only the first', T%nchan ,'will be used'
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
               stop 'wrong dmu'
            !else
            !   T%nchan = ncol
            end if
            allocate(S%d_mu(T%nchan))
            read(vals(i), *)(S%d_mu(j), j = 1, T%nchan)
         else if(keys(i) == 'x_fine') then
            ncol = from_line_get_str_ncol(vals(i))
            G%n_x_fine = ncol
            !print *, 'x_fine :', ncol
            allocate(G%nxdiv(ncol), G%xdiv1(ncol), G%xdiv2(ncol), xrange(ncol))
            read(vals(i), *) (xrange(j), j = 1, ncol)
            !print *, xrange(:)
            do j = 1, ncol
               call retrieve_xrange(xrange(j), G%xdiv1(j), &
                                    G%xdiv2(j), G%nxdiv(j))
            end do 
            !print *, G%xdiv2(:)
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
            !print *, G%E_cont_i, G%E_cont_f
         else if(keys(i) == 'ip_unit') then
            read(vals(i), *) T%IPunit
         else if(keys(i) == 'discrete_method') then
            read(vals(i), *) CTR%discrete_method ! jhangz / level
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
      print *, "read_mqdt_input_key done"
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
         CTR%omega = 0.618_long
      else if(CTR%eqnsolv_method == 'hyb' .and. &
         .not.check_item_in_list(keys,'relax')) then
         CTR%omega = 0.618_long
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
         else if(G%E_cont_f < ZERO) then ! Valid DIS  ! IPs are shifted so that IP1 =0
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
                  &please check input'//trim(FN%f_mqdt)
                  stop
               end if
            else
               call determin_IPx_IPy(T, G)
               !write(stderr, "('IPx(auto)', i5, 'IPy(auto)', i5)")T%IPx, T%Ipy
            end if
         end if
         G%ei = from_E_cm_get_nu(T, G%E_cont_i + T%IPs(1), T%IPx)
         G%ef = from_E_cm_get_nu(T, G%E_cont_f + T%IPs(1), T%IPx)
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
               &please check input'//trim(FN%f_mqdt)
               stop
         end if
      end if

      if(.not.check_item_in_list(keys, 'nop')) then
         call determin_nop(T, G)
         if(T%nop <= 0) then
            write(STDERR,*)&
              'No nop specified and failed to generate due to &
              &wrong Erange/xrange input'//trim(FN%f_mqdt)
            stop 
         end if
      else 
         j = T%nop
         call determin_nop(T, G)
         if(T%nop /= j) then
            write(STDERR,*)&
            'nop specified but different from auto-generated one, please &
            &check input'//trim(FN%f_mqdt)
            stop
         end if
      end if
      if(maxval(T%IP_seq) > T%n_IP) stop 'More IP index then n_IP'
      if(size(T%IP_seq) /= T%nchan) stop 'Less IP index than nchan'
      if(.not.check_item_in_list(keys, 'dmu')) then
         allocate(S%d_mu(T%nchan))
         S%d_mu(:) = ZERO
      end if
      !print *, 'G%E_cont_i',G%E_cont_i/Rydberg
      !print *, 'G%E_cont_f',G%E_cont_f/Rydberg
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
         CTR%omega = 0.618_long
      else if(CTR%eqnsolv_method == 'hyb' .and. &
         .not.check_item_in_list(keys,'relax')) then
         CTR%omega = 0.618_long
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
                  &please check input'//trim(FN%f_mqdt)
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
               &please check input'//trim(FN%f_mqdt)
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
              &wrong Erange/xrange input'//trim(FN%f_mqdt)
            stop 
         end if
      else 
         j = T%nop
         call determin_nop(T, G)
         if(T%nop /= j) then
            write(STDERR,*)&
            'nop specified but different from auto-generated one, please &
            &check input'//trim(FN%f_mqdt)
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


   subroutine read_Smat_input(FN, S)
      implicit none
      type(fname), intent(in) :: FN
      type(Smat), intent(inout) :: S
      integer, parameter:: funit = 15
      integer :: i, j, iost, nchan_t
      real(long) :: E_t
      real(long), allocatable, dimension(:) :: u_tmp, mu_tmp
      logical :: file_ex
      inquire(file = FN%f_Smat, exist = file_ex)
      if(.not. file_ex) then
         print *, 'Failed to open '//trim(FN%f_Smat)//', Please check using -path '
         stop
      end if
      open(funit, file = FN%f_Smat, action = 'read', form = 'formatted')
      read(funit, *,iostat = iost)S%nchan
      rewind(funit)
      S%n_ang = (S%nchan - 1) * S%nchan / 2 
      print *, S%nchan, S%n_ang

      S%n_es = 0
      allocate(mu_tmp(S%nchan), u_tmp(S%n_ang))
      print *, 'allocate(mu_tmp(S%nchan), u_tmp(S%n_ang))'
      do
         read(funit, *,iostat = iost)nchan_t, E_t, (mu_tmp(i), i = 1, S%nchan),&
             (u_tmp(j), j = 1, S%n_ang)
         if(nchan_t /= S%nchan) stop 'dimension wrong reading miuang.out'
         !print *, E_t
         if(iost > 0) then
            write(*,*) trim(FN%f_Smat)//' wrong!'
            exit
         else if(iost < 0) then
            !write(*,*) trim(FN%f_Smat)//' end:',S%n_es
            exit
         else
            S%n_es = S%n_es + 1
            !mu_tmp(:) = mu_tmp(:) + S%d_mu(:)
            call AddToList(S%es, E_t)
            !print * , 'AddToList(S%es, E_t)'
            call AddTo2D(S%mu_s, mu_tmp)
            !print *, 'AddTo2D(S%mu_s, mu_tmp)'
            call AddTo2D(S%ang_s, u_tmp)
         end if
      end do 
      close(funit)
      deallocate(mu_tmp, u_tmp)
   end subroutine read_Smat_input

   subroutine read_Smat_input_block(FN, S)
      implicit none
      type(fname), intent(in) :: FN
      type(Smat), intent(inout) :: S
      integer, parameter:: funit = 15
      integer :: nchan_t, ncol
      real(long), allocatable, dimension(:,:) :: data_block
      logical :: file_ex
      inquire(file = FN%f_Smat, exist = file_ex)
      if(.not. file_ex) then
         print *, 'Failed to open '//trim(FN%f_Smat)//', Please check using -path '
         stop
      end if
      call read_block_data(data_block, S%n_es, ncol, FN%f_Smat)
      nchan_t = int(data_block(1,1))
      if(nchan_t * (nchan_t + 1) / 2 /= ncol - 2) then
         write(STDERR, "('wrong dimension', i5, i5)")nchan_t, ncol -2
         stop 
      end if
      allocate(S%es(S%n_es), S%mu_s(S%n_es, S%nchan), &
               S%ang_s(S%n_es, S%n_ang))
      S%es(:) = data_block(:,2)
      S%mu_s(:, 1:S%nchan) = data_block(:,3:S%nchan+2)
      S%ang_s(:, 1:S%n_ang) = data_block(:,S%nchan+3:ncol)
      if(allocated(data_block)) deallocate(data_block)
   end subroutine read_Smat_input_block

   subroutine initiate(FN, T, S, CTR, G)
      implicit none
      type(mqdt), intent(out) :: T
      type(grid), intent(out) :: G
      type(fname), intent(out) :: FN
      type(control), intent(inout) :: CTR
      type(Smat), intent(out) :: S
      character(len = short_str), dimension(:), allocatable :: keys
      integer :: i, ie
      character(len = short_str), allocatable, dimension(:) :: args

      CTR%path_mqdt = '.'
      CTR%couple = 'ls'
      CTR%yes_dis = .false.
      CTR%yes_au = .false.
      CTR%omega = 0.999_long
      CTR%dbg_F = .false.; CTR%dbg_xgrid = .false.; CTR%dbg_ygrid = .false.
      CTR%dbg_solu = .false.; CTR%dbg_interp = .false.; CTR%dbg_A = .false.
      CTR%dbg_Tir = .false.; CTR%dbg_eqnsolv = .false.; CTR%dbgstr = ''
      CTR%dbg_para_use = .false.

      call read_args(CTR, args)

      FN%f_dbg = trim(CTR%path_mqdt)//"/debug-mqdt.out"
      open(funit_dbg, file = FN%f_dbg, action = "write")
      
      FN%f_mqdt = trim(CTR%path_mqdt)//"/mqdt.in"
      print *, trim(FN%f_mqdt)
      FN%f_Smat = trim(CTR%path_mqdt)//"/miuang.out"
      print *, trim(FN%f_Smat)
      call read_Smat_input(FN, S)
      T%nchan = S%nchan
      T%n_ang = (T%nchan - 1) * T%nchan / 2 
      S%n_ang = T%n_ang
      call read_mqdt_input_key(FN, CTR, T, S, G, keys, args)
      do ie =1, S%n_es
        S%mu_s(ie,:) = S%mu_s(ie,:) + S%d_mu(:)
      end do 


      FN%f_tau = trim(CTR%path_mqdt)//"/tau.out"; FN%funit_tau = 15
      FN%f_Ara = trim(CTR%path_mqdt)//"/Ara.out"; FN%funit_Ara = 16
      FN%f_Anorm = trim(CTR%path_mqdt)//"/Anorm.out"; FN%funit_Anorm =17
      FN%f_Dn = trim(CTR%path_mqdt)//"/Dn.out"; FN%funit_Dn = 18
      FN%f_normdos = trim(CTR%path_mqdt)//"/NormDos.out"; FN%funit_normdos = 19
      FN%f_para_use = trim(CTR%path_mqdt)//"/para_use.out"; FN%funit_para_use = 20

      allocate(T%U(T%nchan, T%nchan), T%Angle(T%n_ang), &
          T%Euler(T%nchan,T%nchan), T%mu(T%nchan), T%F(T%nchan, T%nchan))
      allocate(T%Anorm(T%nchan, T%nop), T%Ara(T%nchan, T%nop), &
          T%Dn(T%nchan), T%Tir(T%nchan, T%nop))
      if(CTR%yes_dis) allocate(T%norm_dos(T%nchan))
      allocate(T%tau(T%nop), T%ty(2,T%nop), T%deriv(2,T%nop))
      call comment_str_among_char("Parameters", '-', funit_dbg)
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
      call print_vector_comment(T%IPs, T%n_IP, print_col_l_fmt, &
           '(E20.10)', funit_dbg)
      write(funit_dbg, "('# IP index for each chann:')", advance = 'yes')
      call print_int_vector_comment(T%IP_seq, T%nchan, print_col_s_fmt, &
           '(i4)', funit_dbg)

      write(funit_dbg, "('# Open chann index :')", advance = 'no')
      call print_int_vector(T%iopen, T%nop, T%nop, '(i5)', funit_dbg)
      write(funit_dbg, "('# x refine', i3)") G%n_x_fine
      do i = 1, G%n_x_fine
         write(funit_dbg, "('# ')", advance = 'no')
         write(funit_dbg, "(2f8.4, i8)")G%xdiv1(i), G%xdiv2(i), G%nxdiv(i)
      end do
      write(funit_dbg, "('# y refine', i3)") G%n_y_fine
      do i = 1, G%n_y_fine
         write(funit_dbg, "('# ')", advance = 'no')
         write(funit_dbg, "(2f8.4, i8)")G%ydiv1(i), G%ydiv2(i), G%nydiv(i)
      end do 
      write(funit_dbg, "('# Smatrix file:', A, ' at ', i3, ' energies')")&
         trim(FN%f_Smat), S%n_es
      if(.not.check_item_in_list(keys, 'dmu')) then
         write(funit_dbg, "('# Smatrix correction: ')", advance = 'no')
         write(funit_dbg, "(' [no/off]')")
      else 
         write(funit_dbg, "('# Smatrix correction: ')", advance = 'yes')
         call print_vector_comment(S%d_mu, T%nchan, print_col_s_fmt, &
              '(f8.4)', funit_dbg)
      end if
      write(funit_dbg, "('# equation solver: ', A)", advance = 'no') &
         CTR%eqnsolv_method
      if(.not.check_item_in_list(keys, 'eqnsolv_method')) then
         write(funit_dbg, "(' (internal)')")
      else
         write(funit_dbg, *)
      end if
      if(CTR%eqnsolv_method == 'newton') then
         write(funit_dbg, "('# relaxing coefficient: ', f8.4)",&
            advance = 'no')CTR%omega
         if(.not.check_item_in_list(keys, 'relax'))then
            write(funit_dbg, "(' (internal)')")
         else
            write(funit_dbg, *)
         end if 
      end if
      call comment_str_among_char("Parameters", '-', funit_dbg)
      write(funit_dbg,*)
      deallocate(keys)
      deallocate(args)
      return 
   end subroutine initiate

   subroutine finalize(T, S, G)
      type(mqdt), intent(inout) :: T
      type(Smat), intent(inout) :: S
      type(grid), intent(inout) :: G
      if(allocated(T%IP_seq)) deallocate(T%IP_seq)
      if(allocated(T%iopen)) deallocate(T%iopen)
      if(allocated(T%IPs)) deallocate(T%IPs)
      if(allocated(T%tau)) deallocate(T%tau)
      if(allocated(T%Ara)) deallocate(T%Ara)
      if(allocated(T%Anorm)) deallocate(T%Anorm)
      if(allocated(T%norm_dos)) deallocate(T%norm_dos)
      if(allocated(T%Dn)) deallocate(T%Dn)
      if(allocated(T%U)) deallocate(T%U)
      if(allocated(T%F)) deallocate(T%F)
      if(allocated(T%Tir)) deallocate(T%Tir)
      if(allocated(T%mu)) deallocate(T%mu)
      if(allocated(T%Angle)) deallocate(T%Angle)
      if(allocated(T%ty)) deallocate(T%ty)
      if(allocated(T%deriv)) deallocate(T%deriv)

      if(allocated(S%mu_s)) deallocate(S%mu_s)
      if(allocated(S%ang_s)) deallocate(S%ang_s)
      if(allocated(S%es)) deallocate(S%es)
      if(allocated(S%d_mu)) deallocate(S%d_mu)
      if(allocated(S%spec_jj)) deallocate(S%spec_jj)
      if(allocated(S%spec_ls)) deallocate(S%spec_ls)
      if(allocated(G%dx)) deallocate(G%dx)
      if(allocated(G%dy)) deallocate(G%dy)
      if(allocated(G%nxdiv)) deallocate(G%nxdiv)
      if(allocated(G%nydiv)) deallocate(G%nydiv)
      if(allocated(G%xdiv1)) deallocate(G%xdiv1)
      if(allocated(G%xdiv2)) deallocate(G%xdiv2)
      if(allocated(G%ydiv1)) deallocate(G%ydiv1)
      if(allocated(G%ydiv2)) deallocate(G%ydiv2)
      if(allocated(G%n_step_seg_y)) deallocate(G%n_step_seg_y)
      if(allocated(G%n_step_seg_x)) deallocate(G%n_step_seg_x)
      if(allocated(G%xseg)) deallocate(G%xseg)
      if(allocated(G%yseg)) deallocate(G%yseg)
      close(funit_dbg)
   end subroutine finalize
end module file_cmdl_io
