module search_root
use constants
use type_def
use numerical
use darray
use interpolate
use stdio
real(long) :: RESONANCE_THRESH
contains
   subroutine findlev(FN, T, CTR, S, phys, spec)
      implicit none
      type(mqdt), intent(inout) :: T
      type(fname), intent(in) :: FN
      type(control), intent(in) :: CTR
      type(Smat), intent(in) :: S
      type(physobserv),intent(inout) :: phys   
      type(spectrodata),intent(inout) :: spec
      if(CTR%yes_dis .and.(.not.CTR%yes_au)) then
          call search_cross(T, phys, CTR%dbg_levfind)  
      else if(CTR%yes_au .and.(.not.CTR%yes_dis)) then
         call resonance_width(T, S, RESONANCE_THRESH)
         if(index(CTR%au_peak_search_method, 'chan') /= 0) then
            print *, 'search_eigen_tau_steep'
            call search_eigen_tau_steep(T, CTR, phys, FN)  
         else if(index(CTR%au_peak_search_method, 'tot') /= 0) then
            print *, 'search_tot_tau_steep'
            call search_tot_tau_steep(T, CTR, phys, FN)  
         end if
      end if
      call get_exp_data(FN,T,spec, phys, CTR)
      call print_phys_tab(FN, CTR, T, phys, spec)
   end subroutine findlev

   subroutine resonance_width(T, S, wi)
      implicit none
      type(mqdt), intent(in) :: T
      type(Smat), intent(in) :: S
      real(long) :: wi
      real(long), dimension(T%nclose) :: width
      real(long), dimension(T%nchan,T%nchan) :: Euler, U
      integer :: ie
      ie = int((1+S%n_es)/2.0)
      call pack_angle(S%ang_s(ie,:), Euler, T%nchan)
      call UANGCON(U, Euler, T%nchan, 1) ! angle ==> U matrix

      call Kmat(T, S%mu_s(ie,:), U, width)
      wi = ONE / minval(width)
   end subroutine resonance_width

   function convert_exp_elev_unit(E_unit0, elev) result(elev2)
      implicit none
      character(short_str), intent(in) :: E_unit0
      character(short_str) :: E_unit
      real(long), intent(in) :: elev
      real(long) :: elev2
      elev2 = elev
      E_unit = all_lower_case_str(E_unit0)
      if((index(E_unit, 'eV') /= 0) .or. (index(E_unit, 'ev') /= 0)) then
         elev2 = Energy_convert_eV_to_cm(elev)
      else if((index(E_unit, 'ryd') /= 0) .or. (index(E_unit, 'Ryd') /= 0) .or. &
         (index(E_unit, 'ry') /= 0) .or. (index(E_unit, 'Ry') /= 0)) then
         elev2 = Energy_convert_Ryd_to_cm(elev)
      else if((index(E_unit, 'a.u.') /= 0) .or. (index(E_unit, 'au') /= 0) .or.&
         (index(E_unit, 'hartree') /= 0) .or. (index(E_unit, 'Hartree') /= 0) ) then
         elev2 = Energy_convert_au_to_cm(elev)
      end if
   end function
   subroutine get_exp_data(FN,T,spec, phys, CTR)
      implicit none
      type(fname), intent(in) :: FN
      type(mqdt), intent(in) :: T
      type(physobserv),intent(in) :: phys   
      type(control), intent(in) :: CTR
      integer, parameter :: funit = 1
      integer :: iost, nlev
      real(long) :: elev
      real(long), dimension(:), allocatable :: tmp
      integer :: chan_id, lev_indx, i
      type(spectrodata),intent(inout) :: spec
      logical :: file_exist
      character(short_str) :: E_unit
      inquire(file = FN%f_exp_tau, exist = file_exist)
      if(.not.file_exist) return
      open(funit, file = FN%f_exp_tau, action = 'read')
      spec%nlev = 0
      read(funit, '(A)') E_unit
      write(funit_dbg, *) "Experimental value in units of "//E_unit
      nlev = 0
      do 
          read(funit, *, iostat = iost) lev_indx, elev, chan_id
          if(iost > 0) then
             write(STDERR,*) trim(FN%f_exp_tau)//' wrong!'
             exit
          else if(iost < 0 ) then
             exit
          else
             elev = convert_exp_elev_unit(E_unit, elev) ! whatever input unit, convert to cm-1
             ! Be careful, exp and theory might be far apart, so we might filter out levels we                        ! need when dealing with a narrow energy range. 
             !if(elev > T%E(T%n_e)*Rydberg .or. elev < T%E(1)*Rydberg) cycle
             nlev = nlev + 1
             call AddToList(spec%energy, elev)
             call Add_int_ToList(spec%lev_indx, lev_indx)
             call Add_int_ToList(spec%eigenchan_id, chan_id)
          end if
      end do 
      close(funit)
      spec%nlev = nlev
      print*,trim(int2str(nlev)),' exp levels loaded and classified:'
      ! to make the exoerimental data index array the same size as theory.
      ! so we can print them together and do arithmic operations.
      if(nlev < phys%nlev) then
         do i = nlev + 1, phys%nlev
            call Add_int_ToList(spec%lev_indx, 0)!, phys%lev_indx(i))
            call Add_int_ToList(spec%eigenchan_id, 0)! phys%eigenchan_id(i))
            call AddToList(spec%energy, nan)
         end do 
         spec%nlev = phys%nlev
         nlev = phys%nlev
      end if
      allocate(spec%lev_nux(nlev), spec%lev_nuy(nlev), tmp(nlev))
      spec%lev_nux(:) =(/(from_E_cm_get_nu(T,spec%energy(i), T%ipx), i = 1, nlev)/)
      if(CTR%yes_dis) then
         tmp(:) =(/(from_E_cm_get_nu(T,spec%energy(i), T%ipy), i = 1, nlev)/)
         spec%lev_nuy(:) =(/(from_nuy_get_tau(tmp(i)), i = 1, nlev)/)
      else
         spec%lev_nuy(:) =(/(ZERO, i = 1, nlev)/)
      end if
      deallocate(tmp)
      print *, "Exp levels : nux chan index"
      do i = 1, nlev
         print "(A5,2x,'root',1x,A,2x,A,2x,A)",ordinal(i),trim(float2str(spec%lev_nux(i),'(f10.5)')),&
         trim(int2str(spec%eigenchan_id(i))),trim(int2str(spec%lev_indx(i)))
      end do
      return 
   end subroutine get_exp_data

   subroutine print_phys_tab(FN, CTR, T, phys, spec)
      implicit none
      type(fname), intent(in) :: FN
      type(control), intent(in) :: CTR
      type(mqdt), intent(in) :: T
      type(physobserv),intent(inout) :: phys   
      type(spectrodata),intent(in) :: spec
      if(CTR%yes_dis) then
         call print_phys_tab_dis(FN, CTR, T, phys, spec)
      else
         call print_phys_tab_au(FN, CTR, T, phys, spec)
      end if
   end subroutine print_phys_tab

   subroutine print_phys_tab_au(FN, CTR, T, phys, spec)
      implicit none
      type(fname), intent(in) :: FN
      type(control), intent(in) :: CTR
      type(mqdt), intent(in) :: T
      type(physobserv),intent(inout) :: phys   
      type(spectrodata),intent(in) :: spec
      integer, parameter :: funit = 1
      real(long), dimension(phys%nlev+spec%nlev, T%nchan) :: spec_x, &
         spec_y, theo_x, theo_y, spec_E
      integer :: i, ny, ichan, lev_ind, ich
      real(long), parameter :: dy = 1.e-3_long
      real(long) :: nuy_min, nuy_max, nuy, xt, yt, Et, xt2
      real(long), dimension(:), allocatable :: nux, tau
      theo_x(:,:) = nan
      theo_y(:,:) = nan
      do i = 1, phys%nlev
         ichan = phys%eigenchan_id(i)
         lev_ind = phys%lev_indx(i)
         theo_x(lev_ind, ichan) = phys%lev_nux(i)
         theo_y(lev_ind, ichan) = phys%lev_nuy(i)
         !print *, 'ich',ichan, 'lev_ind',lev_ind,theo_y(lev_ind, ichan)
      end do
      open(funit, file = FN%f_o_theo_plot, action = 'write')
      write(funit, '("# (nu_x, nu_y) for: ", i5, " channels IPx = ",i3)') T%nchan, T%ipx
      do i = 1, phys%nlev
         write(funit, '(E20.10)', advance = 'no')phys%lev_nux(i)
         lev_ind = phys%lev_indx(i)
         ichan = phys%eigenchan_id(i)
         do ich = 1, T%nchan
            if(ichan == ich) then
               !print *, 'ich',ich, 'lev_ind',lev_ind,theo_y(lev_ind, ich)
               write(funit, '(2E20.10)', advance = 'no') theo_y(lev_ind, ich)
            else
               write(funit, '(2E20.10)', advance = 'no') nan
            end if
         end do 
         write(funit, '(e20.10)', advance = 'yes')
      end do 
      close(funit)

      open(funit, file = FN%f_o_phys, action = 'write')
      write(funit, '("# indx   nu_x, nu_y,  E(ryd)   eigenchan ")')
      write(funit, '("# [* marks the untrustworthy levels, due to ")')
      write(funit, '("# the sparse grid for high energy]")')
      write(funit, '(3i5, "chan", 3E20.10, 2x, A)', advance = 'yes') &
          (i, phys%lev_indx(i), phys%eigenchan_id(i),phys%lev_nux(i), &
          phys%lev_nuy(i), phys%lev_ryd(i), &
          phys%mark(i), i = 1, phys%nlev)
      close(funit)
      if(spec%nlev > 0) then
         spec_x(:,:) = nan
         spec_y(:,:) = nan
         spec_E(:,:) = nan
         ! always spec%nlev >= phys%nlev
         do i = 1, phys%nlev
            ! Here guarantee that spec and theory lev of same order
            ichan = spec%eigenchan_id(i) 
            lev_ind = spec%lev_indx(i) 
            if(ichan == 0 .or. lev_ind == 0) cycle
            spec_x(lev_ind, ichan) = spec%lev_nux(i)   
            spec_y(lev_ind, ichan) = theo_y(lev_ind, ichan)
            spec_E(lev_ind, ichan) = Energy_convert_cm_to_Ryd(spec%energy(i))
         end do
         open(funit, file = FN%f_o_exp_plot, action = 'write')
         write(funit, '("# (nu_x, nu_y) for: ", i5, " channels IPx = ",i3)') T%nchan, T%ipx
         ! always spec%nlev >= phys%nlev
         do i = 1, phys%nlev
            ! Here guarantee that spec and theory lev of same order
            ichan = phys%eigenchan_id(i)
            lev_ind = phys%lev_indx(i)
            write(funit, '(E20.10)', advance = 'no')spec_x(lev_ind, ichan)
            do ich = 1, T%nchan
               if(ichan == ich) then
                  write(funit, '(2E20.10)', advance = 'no') spec_y(lev_ind, ich)
               else 
                  write(funit, '(2E20.10)', advance = 'no') nan
               end if
            end do 
            write(funit, '(e20.10)', advance = 'yes')
         end do 
         close(funit)

         open(funit, file = FN%f_o_phys_compare, action = 'write')
         write(funit, '("# Only print out relevant levels in comparison experiment data, &
            &all other levels can be found in phystable.out")')
         write(funit, '("# indx Theo_x, Theo_y, Exp_x, &
            &[Theo_x-Exp_x], eigenchan")')
         do i = 1, phys%nlev
            lev_ind = phys%lev_indx(i)
            !xt = spec_x(i,phys%eigenchan_id(i))
            !xt = spec_x(lev_ind,phys%eigenchan_id(i))
            ! in order to compare, we print resonance x-position w.r.t. to
            ! corresponding IP 
            xt = from_nux_get_nu(T, xt, T%IP_seq(phys%eigenchan_id(i))) 
            xt2 = from_nux_get_nu(T, phys%lev_nux(i), T%IP_seq(phys%eigenchan_id(i)))
            if(check_NaN(xt)) cycle
            write(funit, '(i5, 4E20.10, i5)', advance = 'yes') i, &
                 xt2, phys%lev_nuy(i), xt, &
                 xt2 - xt, &
                 phys%eigenchan_id(i)
         end do 
         close(funit)
      end if
   end subroutine print_phys_tab_au

   subroutine print_phys_tab_dis(FN, CTR, T, phys, spec)
      implicit none
      type(fname), intent(in) :: FN
      type(control), intent(in) :: CTR
      type(mqdt), intent(in) :: T
      type(physobserv),intent(inout) :: phys   
      type(spectrodata),intent(in) :: spec
      integer, parameter :: funit = 1
      real(long), dimension(phys%nlev+spec%nlev, T%nchan) :: spec_x, &
         spec_y, theo_x, theo_y, spec_E
      integer :: i, ny, ichan, lev_ind, ich
      real(long), parameter :: dy = 1.e-3_long
      real(long) :: nuy_min, nuy_max, nuy, xt, yt, Et, xt2
      real(long), dimension(:), allocatable :: nux, tau
      open(funit, file = FN%f_o_energyrel, action = 'write')
      write(funit, '("#   nu_x nu_y")')
      nuy_min = get_nuy(T, T%nux_sort(1))
      i = T%n_e_sort
      do while(i > 1)
         if(.not.check_NaN(T%nux_sort(i))) exit
         i = i - 1
      end do 
      nuy_max = get_nuy(T, T%nux_sort(i))
      ny = int((nuy_max-nuy_min) / dy)
      allocate(nux(ny), tau(ny))
      do i = 1, ny
         nuy = nuy_min + dy * dble(i)
         nux(i) = from_nuy_get_nux(T, nuy)
         tau(i) = from_nuy_get_tau(nuy)
      end do 
      do i = 1, ny - 1
         write(funit, '(2E20.10)') nux(i), tau(i)
         if(tau(i+1)- tau(i)  > 0.3) then
            write(funit, '(2E20.10)') nan, nan
         end if
      end do 
      i = ny
      write(funit, '(2E20.10)') nux(i), tau(i)
      close(funit)

      theo_x(:,:) = nan
      theo_y(:,:) = nan
      do i = 1, phys%nlev
         ichan = phys%eigenchan_id(i)
         lev_ind = phys%lev_indx(i)
         theo_x(lev_ind, ichan) = phys%lev_nux(i)
         theo_y(lev_ind, ichan) = phys%lev_nuy(i)
         !print *, 'ich',ichan, 'lev_ind',lev_ind,theo_y(lev_ind, ichan)
      end do
      open(funit, file = FN%f_o_theo_plot, action = 'write')
      write(funit, '("# (nu_x, nu_y) for: ", i5, " channels IPx = ",i3)') T%nchan, T%ipx
      do i = 1, phys%nlev
         write(funit, '(E20.10)', advance = 'no')phys%lev_nux(i)
         lev_ind = phys%lev_indx(i)
         ichan = phys%eigenchan_id(i)
         do ich = 1, T%nchan
            if(ichan == ich) then
               !print *, 'ich',ich, 'lev_ind',lev_ind,theo_y(lev_ind, ich)
               write(funit, '(2E20.10)', advance = 'no') theo_y(lev_ind, ich)
            else
               write(funit, '(2E20.10)', advance = 'no') nan
            end if
         end do 
         write(funit, '(e20.10)', advance = 'yes')
      end do 
      close(funit)

      open(funit, file = FN%f_o_phys, action = 'write')
      write(funit, '("# indx   nu_x, nu_y,  E(ryd)   eigenchan ")')
      write(funit, '("# [* marks the untrustworthy levels, due to ")')
      write(funit, '("# the sparse grid for high energy]")')
      write(funit, '(3i5, "chan", 3E20.10, 2x, A)', advance = 'yes') &
          (i, phys%lev_indx(i), phys%eigenchan_id(i),phys%lev_nux(i), &
          phys%lev_nuy(i), phys%lev_ryd(i), &
          phys%mark(i), i = 1, phys%nlev)
      close(funit)
      if(spec%nlev > 0) then
         spec_x(:,:) = nan
         spec_y(:,:) = nan
         spec_E(:,:) = nan
         ! always spec%nlev >= phys%nlev
         do i = 1, phys%nlev
            ! Here guarantee that spec and theory lev of same order
            ichan = spec%eigenchan_id(i) 
            lev_ind = spec%lev_indx(i) 
            if(ichan == 0 .or. lev_ind == 0) cycle
            spec_x(lev_ind, ichan) = spec%lev_nux(i)   
            spec_y(lev_ind, ichan) = spec%lev_nuy(i)
            spec_E(lev_ind, ichan) = Energy_convert_cm_to_Ryd(spec%energy(i))
         end do
         open(funit, file = FN%f_o_exp_plot, action = 'write')
         write(funit, '("# (nu_x, nu_y) for: ", i5, " channels IPx = ",i3)') T%nchan, T%ipx
         ! always spec%nlev >= phys%nlev
         do i = 1, phys%nlev
            ! Here guarantee that spec and theory lev of same order
            ichan = phys%eigenchan_id(i)
            lev_ind = phys%lev_indx(i)
            write(funit, '(E20.10)', advance = 'no')spec_x(lev_ind, ichan)
            do ich = 1, T%nchan
               if(ichan == ich) then
                  write(funit, '(2E20.10)', advance = 'no') spec_y(lev_ind, ich)
               else 
                  write(funit, '(2E20.10)', advance = 'no') nan
               end if
            end do 
            write(funit, '(e20.10)', advance = 'yes')
         end do 
         close(funit)

         open(funit, file = FN%f_o_phys_compare, action = 'write')
         write(funit, '("# Only print out relevant levels in comparison experiment data, &
            &all other levels can be found in phystable.out")')
         write(funit, '("# indx Theo_x Theo_y Exp_x Exp_y &
            &[Theo_x-Exp_x] [Theo_y-Exp_y] eigenchan")')
         do i = 1, phys%nlev
            lev_ind = phys%lev_indx(i)
            !xt = spec_x(i,phys%eigenchan_id(i))
            ! for bound(resonances) we need to check the x-position, w.r.t. its 
            ! corresponding IP, in order to calibrate
            xt = spec_x(lev_ind,phys%eigenchan_id(i))
            xt = from_nux_get_nu(T, xt, T%IP_seq(phys%eigenchan_id(i))) 
            xt2 = from_nux_get_nu(T, phys%lev_nux(i), T%IP_seq(phys%eigenchan_id(i)))
            if(check_NaN(xt)) cycle
            yt = spec_y(lev_ind,phys%eigenchan_id(i))
            write(funit, '(i5, 6E20.10, i5)', advance = 'yes') i, &
                 xt2, phys%lev_nuy(i), &
                 xt, yt, &
                 xt2 - xt, phys%lev_nuy(i) - yt, &
                 phys%eigenchan_id(i)
         end do 
         close(funit)
      end if

   end subroutine print_phys_tab_dis

   function which_eigenchan(T, ie, iop) result (ich)
      implicit none
      type(mqdt), intent(in) :: T
      integer, intent(in) :: iop, ie
      integer :: ich, i
      real(long) :: maxa, A_i_iop
      maxa = abs(T%Anorm(ie, 1, iop))
      ich = 1
      do i = 2, T%nchan
         A_i_iop = abs(T%Anorm(ie, i, iop))
         if(A_i_iop > maxa) then
            maxa = A_i_iop
            ich = i
         end if
      end do 
      return 
   end function which_eigenchan

   function which_eigenchan2(T, ie) result (ich)
      implicit none
      type(mqdt), intent(in) :: T
      integer, intent(in) :: ie
      integer :: ich, i, iop
      real(long) :: maxa, A_i_iop
      maxa = abs(T%Anorm(ie, 1, 1))
      ich = 1
      do iop = 1, T%nop
         do i = 1, T%nchan
            A_i_iop = abs(T%Anorm(ie, i, iop))
            if(A_i_iop > maxa) then
               maxa = A_i_iop
               ich = i
            end if
         end do 
      end do 
      return 
   end function which_eigenchan2

   subroutine get_lev_indx(phys, T)
      type(mqdt), intent(in) :: T
      type(physobserv),intent(inout) :: phys
      integer :: ilev, ich_eig
      allocate(phys%nlev_per_chan(T%nchan))
      phys%nlev_per_chan(:) = 0
      do ilev = 1, phys%nlev
         ich_eig = phys%eigenchan_id(ilev)
         phys%nlev_per_chan(ich_eig) = phys%nlev_per_chan(ich_eig) + 1
         call Add_int_ToList(phys%lev_indx, phys%nlev_per_chan(ich_eig))
      end do 
      return 
   end subroutine get_lev_indx

   function recover_integer(tau) result(tau2)
      implicit none
      real(long), dimension(:), intent(in) :: tau
      real(long), dimension(size(tau)) :: tau2
      integer :: ie, Ne
      real(long) :: principal_n
      Ne = size(tau)
      tau2(:) = tau(:)
      principal_n = ZERO
      do ie = 3, Ne 
         if((tau2(ie) + principal_n -tau2(ie-1)) * (tau2(ie-1) - tau2(ie-2)) < ZERO) &
           principal_n = principal_n + ONE
         tau2(ie) = tau2(ie) + principal_n
      end do 
   end function recover_integer

   function get_resonance_tau2(T, x_l, nux, ich_eig, ibug) result(tau_res) 
      implicit none
      type(mqdt), intent(in) :: T
      integer :: x_l, x_r, ie, ich_eig, x0
      real(long) :: nux, tau_res
      type(Spline) :: SF
      logical, intent(in) :: ibug
      if(ibug)write(funit_dbg,"('Relocate nux ', E20.10, ' in &
          &T%nux_sort around x(', i8,')',E20.10 )", advance = 'yes') &
           nux, x_l,T%nux_sort(x_l)
      ie = 0
      do 
         x0 = x_l + ie
         !if(check_nan(T%nux_sort(x0-1)) .or. check_nan(T%nux_sort(x0))) exit
         if((nux - T%nux_sort(x0-1)) * (T%nux_sort(x0) - nux)>= ZERO) exit
         x0 = x_l - ie
         if(x0 <= 2) then
          print *, 'xl', x_l, 'ie', ie, nux, T%nux_sort(x0)
          exit
         end if
         if((nux - T%nux_sort(x0-1)) * (T%nux_sort(x0) - nux)>= ZERO) exit
         ie = ie + 1
         !if(ie > 1000)stop
      end do 
      x_l = x0 - 1
      !print *, 'ie',ie
      if(ibug)write(funit_dbg, "('relocated: x0', i8, '(shifted )', i8)", advance = 'no')x0, ie
      if(ibug)write(funit_dbg, "(3E20.10)")T%nux_sort(x_l), nux,T%nux_sort(x0)
      do ie = x0-1, 1, -1
         if(ie - x0 <=  -RES_WINDOW) then
            x_l = ie
            exit 
         else if(ie == 1) then
            x_l = ie
         end if
         if(check_NaN(T%tau2(ie, ich_eig))) then
            x_l = ie + 1
            !print *, 'found nan', x_l
            exit
         end if
      end do 
      do ie = x0, T%n_e_sort
         if(ie - x0 >= RES_WINDOW) then
            x_r = ie
            exit 
         else if(ie == T%n_e_sort) then
            x_r = T%n_e_sort
         end if
         if(check_NaN(T%tau2(ie, ich_eig))) then
            x_r = ie - 1
            !print *, 'found nan', x_r
            exit
         end if
      end do 
      print "(2i5,3F8.4)", x_l, x_r, T%nux_sort(x_l), nux, T%nux_sort(x_r)
      print '("E",F8.3, "eV")', get_E_ryd(T, nux) * 13.6058
      if(x_l == x_r) then
         tau_res = T%tau2(x_l,ich_eig)
      else
         SF = Spline_create(T%nux_sort(x_l:x_r), T%tau2(x_l:x_r, ich_eig))
         if(ibug)write(funit_dbg, "(i8,' well behaved sorted openchan tau sample&
            & points from', i8, ' to', i8)") x_r -x_l, x_l, x_r
         if(ibug)write(funit_dbg, "(2E20.10)")&
             (T%nux_sort(ie), T%tau2(ie, ich_eig), ie = x_l, x_r)
         tau_res = Spline_interpolate(SF, nux)
      end if
   end function get_resonance_tau2

!  bound region
   subroutine finer_cross(T, x, y1, y2, i_m, x_m, principal_n_not_same)
      implicit none
      type(mqdt), intent(in) :: T
      real(long), dimension(:), intent(in) :: x, y1, y2
      real(long), dimension(size(y1)) :: f1_smooth, f2_smooth, Func
      integer :: nw, i_m , i
      real(long) :: x_l, x_r, y_l, y_r, y_m
      real(long), intent(out) :: x_m
      real(long) :: principal_n
      type(Spline) :: SPOT
      logical, intent(out) :: principal_n_not_same
      nw = size(x)
      principal_n_not_same = .false.
      ! can only be smoothed locally.
      f1_smooth(:) = smooth_Func(y1(:))
      !f2_smooth(:) = smooth_Func(y2(:))
      f2_smooth(i_m) = from_nuy_get_tau(y2(i_m))
      principal_n = f2_smooth(i_m) + y2(i_m)
      f2_smooth(1:nw) = principal_n - (/(y2(i), i = 1, nw)/) 
      do i = 1, nw - 1
         if(abs(f2_smooth(i) - f2_smooth(i+1))> ONE) principal_n_not_same = .true.
      end do 
      write(funit_dbg, '(5x,"y1               ", 5f8.4)')(y1(i), i = 1, nw)
      write(funit_dbg, '(5x,"mod(y1,1) smooth ", 5f8.4)')(f1_smooth(i), i = 1, nw)
      write(funit_dbg, '(5x,"principal n      ", f8.4)') principal_n
      write(funit_dbg, '(5x,"y2               ", 5f8.4)')(y2(i), i = 1, nw)
      write(funit_dbg, '(5x,"mod(y2,1) smooth ", 5f8.4)')(f2_smooth(i), i = 1, nw)
      if(principal_n_not_same) then
         write(funit_dbg, '("WARNING: the points on the energy relation (and mqdt result tau) are so sparse")')
         write(funit_dbg, '(8x,"that two neighboring Y values on it seperate more than one PI, one should be")')
         write(funit_dbg, '(8x,"very suspecious about this level. I will mark a star for this level in the")')
         write(funit_dbg, '(8x,"summary file: phystable.out. [TODO: increase grid locally and interpolate mqdt")')
         write(funit_dbg, '(8x,"result.]")') 
      end if
      Func(:) = f1_smooth(:) - f2_smooth(:)
      x_l = x(i_m)
      x_r = x(i_m+1)
      y_l = Func(i_m)
      y_r = Func(i_m+1)
      !SPOT = Spline_create(x, Func)
      SPOT = Spline_create(x, f1_smooth)
      do 
         x_m = (x_l + x_r) / TWO
         y_m = Spline_interpolate(SPOT, x_m) &
               - (principal_n - get_nuy(T, x_m))
         if(abs(y_m) <= 1.e-10_long)  return
         if(y_l * y_m > ZERO) then
            x_l = x_m 
            y_l = Spline_interpolate(SPOT, x_l) &
               - (principal_n - get_nuy(T, x_l))
         else if(y_r * y_m > ZERO) then
            x_r = x_m 
            y_r = Spline_interpolate(SPOT, x_r)  &
               - (principal_n - get_nuy(T, x_r))
         else
            return
         end if
      end do 
   end subroutine finer_cross

   function smooth_Func(fun) result(f2)
      implicit none
      integer :: fun_end, i
      real(long), dimension(:), intent(in) :: fun
      real(long), dimension(size(fun)) :: f2
      real(long) :: dfun
      f2(:) = fun(:)
      fun_end = size(fun)
      do i = 1, fun_end - 1
         dfun = f2(i) - f2(i+1)
         if(dfun > 0.8)then
            f2(i+1) = fun(i+1) + ONE
         else if(dfun < -0.8) then
            f2(i+1) = fun(i+1) - ONE
         end if
      end do
   end function smooth_Func

   subroutine rough_cross(Func1, Func2, n_rt, rt_mrks)
      implicit none
      integer, intent(inout) :: n_rt
      integer :: iroam
      real(long) :: principal_n
      real(long), dimension(:), intent(in) :: Func1, Func2
      real(long), dimension(2) :: f1, f2
      integer, dimension(:), allocatable, intent(inout) :: rt_mrks
      n_rt = 0
      do iroam = 1, size(Func1) - 1
         if(check_Nan(Func1(iroam)) .or. check_NaN(Func1(iroam+1))) cycle ! avoid nan cases 
         f1(:) = smooth_Func(Func1(iroam:iroam+1))
         !f2(:) = smooth_Func(Func2(iroam:iroam+1))
         f2(1) = from_nuy_get_tau(Func2(iroam))
         principal_n = f2(1) + Func2(iroam) 
         f2(2) = principal_n - Func2(iroam+1) 
         if( (f1(1)-f2(1)) * (f1(2) - f2(2)) <= ZERO) then
            n_rt = n_rt + 1
            call Add_int_ToList(rt_mrks, iroam)
         end if
      end do 
      return 
   end subroutine rough_cross

   subroutine search_cross(T, phys, ibug) 
      implicit none
      type(mqdt), intent(inout) :: T
      type(physobserv),intent(inout) :: phys
      integer :: ich, n_rt, irt, ie, ilev
      real(long) :: start, finish
      real(long), dimension(T%n_e_sort) :: Func1, Func2, x
      integer, dimension(:), allocatable :: rt_mrks
      integer :: w_l, w_r
      logical :: principal_n_not_same
      character(len = short_str) :: lev_mark
      real(long) :: nux, tau
      logical, intent(in) :: ibug
      call cpu_time(start)

      write(*,"('Searching bound levels ...')")
      write(funit_dbg,"('Searching bound levels ...')")
      phys%nlev = 0
      Func2(:) = (/(get_nuy(T, T%nux_sort(ie)), ie = 1, T%n_e_sort)/) 
      x(:) = T%nux_sort(:)
      loop_chan: do ich = 1, T%nchan
         Func1(:) = T%tau2(:,ich)
         n_rt = 0
         if(allocated(rt_mrks))deallocate(rt_mrks)
         call rough_cross(Func1, Func2, n_rt, rt_mrks)
         write(*,"('Channel', i4)") ich
         write(funit_dbg,"('Channel', i4)") ich
         ilev = 0 ! this should be same as irt, is root search correctly done.
         do irt = 1, n_rt ! very unlikely to how more than one root
               write(funit_dbg,'(A," Cross xi", f8.4, " f1(i) ", f8.4, " f2(i)", f8.4, &
                    &" f1(i+1) ", f8.4, " f2(i+1)", f8.4)') trim(ordinal(irt)),&
                     x(rt_mrks(irt)),Func1(rt_mrks(irt)), Func2(rt_mrks(irt)),&
                     Func1(rt_mrks(irt)+1), Func2(rt_mrks(irt)+1)
               w_l = rt_mrks(irt) - CROSS_WINDOW
               w_r = rt_mrks(irt) + CROSS_WINDOW
               if(w_r >= T%n_e_sort) then
                  w_r = T%n_e_sort
                  w_l = w_r - 2 * CROSS_WINDOW
               else if( w_l < 1) then
                  w_l = 1
                  w_r = w_l + 2 * CROSS_WINDOW
               end if
               if(w_l == w_r) exit
               call finer_cross(T, x(w_l:w_r), Func1(w_l:w_r), &
                      Func2(w_l:w_r), rt_mrks(irt) - w_l + 1, nux, &
                      principal_n_not_same)
               tau = from_nuy_get_tau(get_nuy(T, nux))
               if(ibug)then
                  print '("==> finer search by spline: x ",f8.4, " y ", f8.4)', nux, tau
                  write(funit_dbg, '("==> finer search by spline: x ",f8.4, " y ", f8.4)') nux, tau
               end if
               call AddToList(phys%lev_nux, nux)
               call AddToList(phys%lev_nuy, tau)
               call AddToList(phys%lev_ryd, get_E_ryd(T, nux))
               call Add_int_ToList(phys%eigenchan_id, ich)
               call Add_int_ToList(phys%lev_indx, irt)
               if(principal_n_not_same) then
                  lev_mark = "*"
               else
                  lev_mark = ""
               end if
               call Add_char_ToList(phys%mark, lev_mark)
               phys%nlev = phys%nlev + 1
               ilev = ilev + 1
         end do 
         call Add_int_ToList(phys%nlev_per_chan, ilev) 
      end do loop_chan
      call cpu_time(finish)
      write(*,"(i5, ' levels found, taking ', E20.10, ' sec')") phys%nlev, finish-start
      write(funit_dbg,"(i5, ' levels found, taking ', E20.10, ' sec')") phys%nlev, finish-start
      return 
   end subroutine search_cross

! Autoionization region
   subroutine finer_bisect_root_dtdE(x, y, nw, i_cross, x_m)
      real(long), dimension(:), intent(in) :: x, y
      real(long), dimension(2*nw+1) :: f_smooth
      integer, intent(in) :: i_cross, nw 
      real(long) :: x_l, x_r, y_l, y_r, y_m
      real(long), intent(out) :: x_m
      type(Spline) :: SPOT
      x_l = x(i_cross)
      x_r = x(i_cross+1)
      y_l = y(i_cross)
      y_r = y(i_cross+1)
      f_smooth(:) = smooth_Func(y(i_cross-nw:i_cross+nw))
      SPOT = Spline_create(x(i_cross-nw:i_cross+nw), f_smooth)!y(i_cross-nw:i_cross+nw))
      do 
         x_m = (x_l + x_r) / TWO
         y_m = Spline_interpolate(SPOT, x_m)
         if(abs(y_m) >= EPS)  exit
         if(y_l * y_m > ZERO) then
            x_l = x_m 
            y_l = Spline_interpolate(SPOT, x_l)
         else if(y_r * y_m > ZERO) then
            x_r = x_m 
            y_r = Spline_interpolate(SPOT, x_r)
         else
            exit
         end if
      end do 
   end subroutine finer_bisect_root_dtdE

   subroutine rough_search_eigen_peak_old(CTR, ne, tau, slope, npeak, peaks)
      implicit none
      type(control), intent(in) :: CTR
      real(long), dimension(:) :: tau, slope
      integer, intent(in) :: ne
      integer :: ie
      integer, dimension(:), allocatable :: peaks
      integer, intent(out):: npeak
      real(long) :: filter
      print *, "RESONANCE_THRESH ",RESONANCE_THRESH," (au-1)"
      if(CTR%yes_user_filter) then
         filter = CTR%filter
      else
         !filter = RESONANCE_THRESH / 1.d5
         filter = RESONANCE_THRESH 
      end if
      npeak = 0
      do ie = 2, ne - 1
         if(abs(slope(ie)) > filter .and. & ! filter fluctuation.
            ((slope(ie+1) - slope(ie)) * (slope(ie) - slope(ie-1)) < ZERO) .and. &
            .not. check_NaN(slope(ie)) .and. &
            (tau(ie) - tau(ie -1) ) * (tau(ie+1) - tau(ie)) > ZERO ) then !filter sudden jump.
            npeak = npeak + 1
            call Add_int_ToList(peaks,ie)
         end if
      end do 
      return
   end subroutine rough_search_eigen_peak_old

   subroutine rough_search_eigen_peak(CTR, ne, tau, slope, npeak, peaks)
      implicit none
      type(control), intent(in) :: CTR
      real(long), dimension(:) :: tau, slope
      real(long), dimension(RES_WINDOW * 2 - 1) :: slope_avrg
      integer, intent(in) :: ne
      integer :: ie, ir
      integer, dimension(:), allocatable :: peaks
      integer, intent(out):: npeak
      real(long) :: filter
      print *, "RESONANCE_THRESH ",RESONANCE_THRESH," (au-1)"
      if(CTR%yes_user_filter) then
         filter = CTR%filter
      else
         !filter = RESONANCE_THRESH / 1.d5
         filter = RESONANCE_THRESH 
      end if

      print "(A,e20.10)", "filter out fluctuation below:", filter
      print "(A,I8)", "check dtau/dE in a shifting window ", RES_WINDOW

      npeak = 0
      if(ne <= RES_WINDOW * 2 ) stop "Too few energy points."

      do ie = RES_WINDOW+1, ne - RES_WINDOW - mod(ne, RES_WINDOW), 1
         if(abs(slope(ie)) < filter .or. & ! filter fluctuation.
            ((slope(ie+1) - slope(ie)) * (slope(ie) - slope(ie-1)) > ZERO) .or. & ! 
            check_NaN(slope(ie)) .or. &
            (tau(ie) - tau(ie -1) ) * (tau(ie+1) - tau(ie)) < ZERO ) cycle !filter sudden jump.
         write(stdout, "('loc_max ',i8, e20.10)") ie, slope(ie)
         write(stdout, "('neighbor: ',8e20.10)") &
            slope(ie - RES_WINDOW + 1: ie + RES_WINDOW - 1)

         do ir = -RES_WINDOW + 1, RES_WINDOW -1
            slope_avrg(ir + RES_WINDOW) = sum(slope(ie + ir - RES_WINDOW + 1: ie + ir + RES_WINDOW -1)) / dble(RES_WINDOW * 2-1)
         end do
         write(stdout, "('averaged: ',8e20.10)") &
           slope_avrg
         !if(is_window_center_a_peak(slope(ie - RES_WINDOW + 1: ie + RES_WINDOW -1))) then ! Window always has RES_WINDOW * 2-1 elements.
         if(is_window_center_a_peak(slope_avrg)) then ! Window always has RES_WINDOW * 2-1 elements.
         !if(is_window_center_a_peak(slope_avrg,'force')) then ! Window always has RES_WINDOW * 2-1 elements.
            npeak = npeak + 1
            call Add_int_ToList(peaks,ie)
         end if
      end do
      return
   end subroutine rough_search_eigen_peak
   subroutine rough_search_tot_tau_peak_old(CTR, ne, slope, npeak, peaks)
      implicit none
      type(control), intent(in) :: CTR
      real(long), dimension(:) :: slope
      integer, intent(in) :: ne
      integer :: ie
      integer, dimension(:), allocatable :: peaks
      integer, intent(out):: npeak
      real(long) :: filter
      print *, "RESONANCE_THRESH ",RESONANCE_THRESH," (au-1)"
      if(CTR%yes_user_filter) then 
         filter = CTR%filter
      else
         !filter = RESONANCE_THRESH / 1.d5
         filter = RESONANCE_THRESH 
      end if
      print "(A,e20.10)", "filter out fluctuation below:", filter
      npeak = 0
      do ie = 2, ne - 1
         if(abs(slope(ie)) > filter .and. & ! filter fluctuation.
            ((slope(ie+1) - slope(ie)) * (slope(ie) - slope(ie-1)) < ZERO) .and. &
            .not. check_NaN(slope(ie)) ) then
            npeak = npeak + 1
            call Add_int_ToList(peaks,ie)
         end if
      end do 
      print *, 'rough search:',npeak,"peaks found"
      return
   end subroutine rough_search_tot_tau_peak_old

   subroutine rough_search_tot_tau_peak(CTR, ne, slope, npeak, peaks)
      implicit none
      type(control), intent(in) :: CTR
      real(long), dimension(:) :: slope
      real(long), dimension(RES_WINDOW * 2-1) :: slope_avrg
      integer, intent(in) :: ne
      integer :: ie, ir
      integer, dimension(:), allocatable :: peaks
      integer, intent(out):: npeak
      real(long) :: filter
      print *, "RESONANCE_THRESH ",&
         trim(float2str(RESONANCE_THRESH,'(e16.5)'))," (au-1)"
      if(CTR%yes_user_filter) then 
         filter = CTR%filter
      else
         !filter = RESONANCE_THRESH / 1.d5
         filter = RESONANCE_THRESH 
      end if
      print "(A,A)", "filter out fluctuation below: ", trim(float2str(filter, '(e16.5)'))
      print "(A,A)", "check dtau/dE in a shifting window ", trim(int2str(RES_WINDOW,'(i8)'))

      npeak = 0
      if(ne <= RES_WINDOW * 2 ) stop "Too few energy points."

      do ie = 2 * RES_WINDOW+1, ne - 2 * RES_WINDOW - mod(ne, RES_WINDOW), 1
         ! local maximum at ie. 
         if(abs(slope(ie)) < filter .or. & ! filter fluctuation.
            ((slope(ie+1) - slope(ie)) * (slope(ie) - slope(ie-1)) > ZERO) .or. & 
            check_NaN(slope(ie)) ) cycle
         write(stdout, "('loc_max ',A,1x, A)") trim(int2str(ie)), &
              trim(float2str(slope(ie), '(e14.5)'))
         write(stdout, "('neighbor: ',8e14.5)") &
            slope(ie - RES_WINDOW + 1: ie + RES_WINDOW - 1)
!  Average the function in a window, and shift the window to 
!  look for a window where it's center is a peak.
         do ir = -RES_WINDOW + 1, RES_WINDOW -1
            slope_avrg(ir + RES_WINDOW) = sum(slope(ie + ir - RES_WINDOW + 1: &
            ie + ir + RES_WINDOW -1)) / dble(RES_WINDOW * 2-1)
         end do 
         write(stdout, "('averaged: ',8e14.5)") slope_avrg
         !if(is_window_center_a_peak(slope(ie - RES_WINDOW + 1: ie + RES_WINDOW -1))) then 
         !if(is_window_center_a_peak(slope_avrg, 'force')) then 
         if(is_window_center_a_peak(slope_avrg)) then 
            ! Window always has RES_WINDOW * 2-1 elements.
            npeak = npeak + 1
            call Add_int_ToList(peaks,ie)
         end if
      end do 
      print *, 'rough search:',npeak,"peaks found"
      return
   end subroutine rough_search_tot_tau_peak

!
!  Different from discrete energy region, we search the total original tau (with 
!  integer principal number recovered), so that the sorted function is smooth
!  enough, we then assigen color to the solutions.
!
   subroutine search_tot_tau_steep(T, CTR, phys, FN) 
   ! look for f'' = 0 --> f' max : steepest root 
      implicit none
      type(mqdt), intent(in) :: T
      type(control), intent(in) :: CTR
      type(physobserv),intent(inout) :: phys
      type(fname), intent(in) :: FN
      !real(long), dimension(T%n_e_sort) :: nu, tau, slope, slope2
      real(long), dimension(T%n_e) :: nu, tau, slope, slope2
      integer :: ich, n_rt, irt, x_l, x_r, ie, ich_eig
      integer, dimension(:), allocatable :: rt_mrks
      real(long) :: nux, nuy, fp
      type(Spline) :: SFp
      character(len = short_str) :: lev_mark
      integer :: funit_driv = 100
      phys%nlev = 0
      nu(1:T%n_e) = T%nux(1:T%n_e)
      tau(:) = ZERO
      loop_chan: do ich = 1, T%nop
         tau(1:T%n_e) = tau(1:T%n_e) + T%tau(1:T%n_e, ich)
      end do loop_chan
      tau = recover_integer(tau)
      call difference(T%n_e, nu, tau, slope)
      slope(:) = TWO * slope(:) * nu(:)**3
      call difference(T%n_e, nu, slope, slope2)
      if(CTR%dbg_deriv) then
         open(funit_driv, file = FN%f_dbg_dtdE)
         write(funit_driv, "('# total tau and slope')")
         write(funit_driv,'(3E20.10)')(nu(ie), tau(ie), slope(ie), ie=1, T%n_e)
         close(funit_driv)
      endif 
      write(funit_dbg, "('#RESONANCE_THRESH(internal) ', e20.10)") RESONANCE_THRESH
      call rough_search_tot_tau_peak(CTR, T%n_e, slope, n_rt, rt_mrks)
      write(funit_dbg, "(i5, ' ', A, ' found. (with estimated largest &
            &width: ',E20.10, ' au)' )")n_rt, single_plural('peak',n_rt), &
            ONE / RESONANCE_THRESH
      if(n_rt > 0) then
         do irt = 1, n_rt 
            x_l = rt_mrks(irt) - 2
            x_r = x_l + 4
            write(stdout, "(A,' root ', e20.10)") trim(ordinal(irt)), nu(rt_mrks(irt))
            if(x_r >= T%n_e) exit
            if(CTR%dbg_levfind) then
               write(funit_dbg,"('# Root around ', A, ' point, its 4 &
               &nerighbor points are:')") trim(ordinal(rt_mrks(irt)))
               write(funit_dbg,'(2E20.10)') nu(rt_mrks(irt)-2),slope(rt_mrks(irt)-2)
               write(funit_dbg,'(2E20.10)') nu(rt_mrks(irt)-1),slope(rt_mrks(irt)-1)
               write(funit_dbg,'(2E20.10," <---")') nu(rt_mrks(irt)),slope(rt_mrks(irt))
               write(funit_dbg,'(2E20.10)') nu(rt_mrks(irt)+1),slope(rt_mrks(irt)+1)
               write(funit_dbg,'(2E20.10)') nu(rt_mrks(irt)+2),slope(rt_mrks(irt)+2)
            end if
            SFp = Spline_create(nu(x_l:x_r), slope(x_l:x_r))
            call finer_bisect_root_dtdE(nu, slope2, 5, rt_mrks(irt), nux)
            call AddToList(phys%lev_nux, nux)
            ich_eig = which_eigenchan2(T, rt_mrks(irt))
            call Add_int_ToList(phys%eigenchan_id, ich_eig)
            nuy = get_resonance_tau2(T, rt_mrks(irt), nux, ich_eig, CTR%dbg_levfind)
            write(funit_dbg,'(A," resonance x", f8.4, "tau(",i3,")",f8.4, &
                  &" dtau/dE(i) ", E20.10 )') trim(ordinal(irt)),&
                  nux, ich_eig, nuy, slope(x_l)
                  !nux, ich_eig, nuy, slope(rt_mrks(irt))
            call AddToList(phys%lev_nuy, nuy)
            call AddToList(phys%lev_ryd, get_E_ryd(T, nux))
            lev_mark = ""
            call Add_char_ToList(phys%mark, lev_mark)
            fp = Spline_interpolate(SFp, nux)
            call AddToList(phys%width, ONE/fp)
            phys%nlev = phys%nlev + 1
         end do
         deallocate(rt_mrks)
      end if
      call get_lev_indx(phys, T)
      return
   end subroutine search_tot_tau_steep

!
!  Different from discrete energy region, we don't search the sorted tau,
!  the original tau is smooth enough, we then assigen color to the solutions.
!
   subroutine search_eigen_tau_steep(T, CTR, phys, FN) 
   ! look for f'' = 0 --> f' max : steepest root 
      implicit none
      type(mqdt), intent(in) :: T
      type(control), intent(in) :: CTR
      type(physobserv),intent(inout) :: phys
      type(fname), intent(in) :: FN
      real(long), dimension(T%n_e) :: nu, tau, slope, slope2
      real(long), dimension(T%n_e, T%nop) :: g
      integer :: ich, n_rt, irt, x_l, x_r, ilev, ich_eig, ie
      integer, parameter :: funit_driv = 100
      integer, dimension(:), allocatable :: rt_mrks
      real(long) :: nux, nuy, fp
      type(Spline) :: SF, SFp
      character(len = short_str) :: lev_mark
      phys%nlev = 0
      nu(:) = T%nux(:)
      if(CTR%dbg_deriv) then
         !open(funit_driv, file = FN%f_dbg_dtdE)
      end if
      loop_chan: do ich = 1, T%nop
         write(funit_dbg,"('effective collision channel : rho = ', i4)", advance = 'no') ich
         tau(:) = T%tau(:, ich)
         call difference(T%n_e, nu, tau, slope)
         slope(:) = TWO * slope(:) / nu(:)**3
         g(:,ich) = slope(:)
         call difference(T%n_e, nu, slope, slope2)
         if(CTR%dbg_deriv) then
            open(funit_driv, file = trim(adjustl(FN%f_dbg_dtdE))//'_'//adjustl(ordinal(ich)))
            write(funit_driv, "('# total tau and slope ', A)")trim(ordinal(ich))
            write(funit_driv,'(3E20.10)')(nu(ie), tau(ie), slope(ie), ie=1, T%n_e)
            close(funit_driv)
         endif
         call rough_search_eigen_peak(CTR, T%n_e, tau, slope, n_rt, rt_mrks)
         write(funit_dbg, "(i5, 'peakcs found.')")n_rt
         ilev = 0
         if(n_rt > 0) then
            do irt = 1, n_rt 
               x_l = rt_mrks(irt) - 2
               x_r = x_l + 4
               if(x_r >= T%n_e) exit
               if(CTR%dbg_levfind) then
                  write(funit_dbg,"('Root around ', i8, 'th point, its 4 &
                  &nerighbor points are:')") rt_mrks(irt)
                  write(funit_dbg,'(i8, 2E20.10)') nu(rt_mrks(irt)-2),slope(rt_mrks(irt)-2)
                  write(funit_dbg,'(i8, 2E20.10)') nu(rt_mrks(irt)-1),slope(rt_mrks(irt)-1)
                  write(funit_dbg,'(i8, 2E20.10," <--")') nu(rt_mrks(irt)),slope(rt_mrks(irt))
                  write(funit_dbg,'(i8, 2E20.10)') nu(rt_mrks(irt)+1),slope(rt_mrks(irt)+1)
                  write(funit_dbg,'(i8, 2E20.10)') nu(rt_mrks(irt)+2),slope(rt_mrks(irt)+2)
               end if
               SF = Spline_create(nu(x_l:x_r), tau(x_l:x_r))
               SFp = Spline_create(nu(x_l:x_r), slope(x_l:x_r))
               call finer_bisect_root_dtdE(nu, slope2, 5, rt_mrks(irt), nux)
               write(funit_dbg,'(i3,"th resonance x", f8.4, " dtau/dE(i) ", E20.10 )') irt,&
                     nux, slope(rt_mrks(irt))
               nuy = Spline_interpolate(SF, nux)
               fp = Spline_interpolate(SFp, nux)
               call AddToList(phys%lev_nux, nux)
               call AddToList(phys%lev_nuy, nuy)
               call AddToList(phys%lev_ryd, get_E_ryd(T, nux))
               call AddToList(phys%width, fp)
               ich_eig = which_eigenchan(T, rt_mrks(irt), ich)
               call Add_int_ToList(phys%eigenchan_id, ich_eig)
               lev_mark = ""
               call Add_char_ToList(phys%mark, lev_mark)
               phys%nlev = phys%nlev + 1
            end do
            deallocate(rt_mrks)
         end if
      end do loop_chan
      call get_lev_indx(phys, T)
      if(CTR%dbg_deriv) then
         !close(funit_driv)
      endif
      return
   end subroutine search_eigen_tau_steep
end module search_root
