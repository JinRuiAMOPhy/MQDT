module gengrid
   use constants
   use type_def
   use numerical 
   use stdio
contains
   subroutine uniform_xgrid(G) 
      implicit none
      type(grid), intent(inout) :: G
      real(long) :: dx
      integer :: nstep 
      ! always use T%IPx as the abscissa.
      dx = ONE / dble(G%nx_flat)
      nstep = int((G%ef - G%ei) / dx)
      G%nstep_x = nstep
      if(.not.allocated(G%dx)) allocate(G%dx(nstep))
      G%dx(:) = dx
   end subroutine uniform_xgrid

   subroutine assign_xgrid_in_seg(G)
      implicit none
      type(grid), intent(inout) :: G
      real(long) :: dx
      integer :: i, i2, i3
      i2 = 2; i3 = 2
      G%dx(:) = ZERO
      do i = 2, G%n_seg_x
         dx = (G%xseg(i) - G%xseg(i-1)) / dble(G%n_step_seg_x(i) )
         i3 = i2 + G%n_step_seg_x(i) - 1
         G%dx(i2:i3) = dx
         i2 = i3 + 1
      end do 
   end subroutine assign_xgrid_in_seg

   subroutine uniform_ygrid(G) 
      implicit none
      type(grid), intent(inout) :: G
      real(long) :: dx
      ! always use T%IPx as the abscissa.
      dx = ONE / dble(G%ny_adapt)
      if(.not.allocated(G%dy)) allocate(G%dy(G%ny_adapt))
      G%dy(:) = dx
   end subroutine uniform_ygrid

   subroutine update_ysegment(T, G, CTR, targ)
      implicit none
      type(mqdt), intent(inout) :: T
      type(grid), intent(inout) :: G
      type(control), intent(in) :: CTR
      integer :: i
      real(long), dimension(T%nop), intent(in) :: targ

      G%yseg(1) = ZERO
! The tau should be linear locally, so we extrapolate for the next step, 
! for a better guessing segment than taking the previous tau.
!    
      do i = 2, T%nop + 1
         G%yseg(i) = from_nuy_get_tau(targ(i-1)) 
      end do 
      i = T%nop + 2
      G%yseg(i) = ONE
      call sort_real_vec(G%yseg, 'ascend', G%yseg)
      G%n_step_seg_y(:) = G%ny_adapt
      if(CTR%dbg_ygrid) then
         write(funit_dbg,'("update Y grid at ", f10.5)') T%nux
         call print_vector_comment(G%yseg(:), G%n_seg_y, 15, '(f10.5)', funit_dbg)
         call print_int_vector_comment(G%n_step_seg_y(:), &
            G%n_seg_y, 15, '(i10)', funit_dbg)
      end if
   end subroutine update_ysegment

   subroutine assign_ygrid_in_seg(G)
      implicit none
      type(grid), intent(inout) :: G
      real(long) :: dx
      integer :: i, i2, i3
      i2 = 2; i3 = 2
      do i = 2, G%n_seg_y
         dx = (G%yseg(i) - G%yseg(i-1)) / dble(G%n_step_seg_y(i) )
         i3 = i2 + G%n_step_seg_y(i) - 1
         G%dy(i2:i3) = dx
         i2 = i3 + 1
      end do 
      !print *, G%dy(:)
   end subroutine assign_ygrid_in_seg

   subroutine y_init(T, S, G, CTR, todo)
      implicit none
      type(mqdt), intent(inout) :: T
      type(Smat), intent(in) :: S
      type(grid), intent(inout) :: G
      type(control), intent(in) :: CTR
      character (len = *) , intent(in) :: todo
      integer :: i, iop, i2, n_grid_seg
      logical :: touched
      real(long) :: a, b, c

      if(index(G%ygrid_type, "uniform") /= 0) then
         call uniform_ygrid(G)
         return
      end if

      n_grid_seg = 0
      if(index(todo, "guess") /= 0) then
         call interp_Smat(T%nux, T, S, CTR)
         G%n_seg_y = T%nop + G%n_y_fine * 2 + 2 
         allocate(G%yseg(G%n_seg_y))
         allocate(G%n_step_seg_y(G%n_seg_y))
         n_grid_seg = G%ny_init_guess
      else if(index(todo, "reinitiate") /= 0) then 
         G%n_seg_y = T%nop + 2
         deallocate(G%yseg, G%n_step_seg_y)
         allocate(G%yseg(G%n_seg_y), G%n_step_seg_y(G%n_seg_y))
         n_grid_seg = G%ny_adapt
      end if
      if(index(todo, "renitiate-wider") /= 0) then
        G%yseg(1) = -0.1_long 
      else 
        G%yseg(1) = ZERO
      end if
      do i = 2, T%nop + 1
         iop = T%iopen(i-1)
         G%yseg(i) = from_nuy_get_tau(T%mu(iop)) 
      end do 
      i = T%nop + 2
      if(index(todo, "guess") /= 0) then
         do i2 = 1, G%n_y_fine 
            G%yseg(i) = from_nuy_get_tau(G%ydiv1(i2))
            i = i + 1
            G%yseg(i) = from_nuy_get_tau(G%ydiv2(i2))
            i = i + 1
         end do 
         i = 1 + T%nop + G%n_y_fine * 2 + 1
      else if(index(todo, "reinitiate") /= 0) then
         i = 1 + T%nop + 1
      end if
      G%yseg(i) = ONE

      call sort_real_vec(G%yseg, 'ascend', G%yseg)
      if(index(todo, "guess") /= 0 ) then
         do i = 1, G%n_seg_y
            touched = .false.
            b = G%yseg(i)
            do i2 = 1, G%n_y_fine
               a = from_nuy_get_tau(G%ydiv2(i2))
               c = from_nuy_get_tau(G%ydiv1(i2))
               if(is_ascend(a, b, c)) then
                  G%n_step_seg_y(i) = n_grid_seg * G%nydiv(i2)
                  touched = .true.
                  exit 
               end if
            end do 
            if(.not.touched) G%n_step_seg_y(i) = n_grid_seg 
         end do 
      else if(index(todo, "reinitiate") /= 0 ) then
         G%n_step_seg_y(:) = n_grid_seg
      end if

      G%nstep_y = sum(G%n_step_seg_y(:)) - G%n_step_seg_y(1) + 1
      call comment_str_among_char("Y Grid", '-', funit_dbg)
      if(index(todo, "guess") /= 0 ) then
        write(funit_dbg, "('# Initial guess for y grid &
           &segments: (',i5, 'segments ) ')") G%n_seg_y
      else if(index(todo, "reinitiate") /= 0 ) then
        write(funit_dbg, "('# reconstructing y grid segments &
           &based on solutions : (',i5, 'segments ) ')") G%n_seg_y
      end if
      call print_vector_comment(G%yseg(:), G%n_seg_y, 15, '(f10.5)', funit_dbg)
      call print_int_vector_comment(G%n_step_seg_y(:), &
         G%n_seg_y, 15, '(i10)', funit_dbg)
      write(funit_dbg, "('# The number of y steps', i6,'(&
         &for bisect guess.)')") G%nstep_y
      if(index(todo, "guess") /= 0 )write(funit_dbg, &
         "('# Later the y grids will be updated &
         &according to the previous solution.')")
      if(index(todo, "guess") /= 0) then
         allocate(G%dy(G%nstep_y))
      else if(index(todo, "reinitiate") /= 0) then
         deallocate(G%dy)
         allocate(G%dy(G%nstep_y))
      end if
      if(CTR%dbg_ygrid) then
         open(666,file = trim(CTR%path_mqdt)//"/ygrid.out", action = "write")
         write(666, "('# y grid segments.')")
         do i = 1, G%n_seg_y
            write(666, "(2E20.10)") G%ei, ONE - G%yseg(i), G%ef, ONE - G%yseg(i)
            write(666, "('#')")
         end do
         close(666)
      end if
      if(index(todo, "guess") /= 0 ) then
         call assign_ygrid_in_seg(G)
      else if(index(todo, "reinitiate") /= 0) then
         call update_ysegment(T, G, CTR, T%tau(:))
         call assign_ygrid_in_seg(G)
      end if
      return 
   end subroutine y_init

   subroutine update_ygrid(x_next, T, G, CTR, opt)
      type(mqdt), intent(inout) :: T
      type(grid), intent(inout) :: G
      type(control), intent(in) :: CTR
      real(long), intent(in) :: x_next
      real(long), dimension(T%nop) :: targ
      integer :: i
      character(len = *), intent(in), optional :: opt
      if(G%ygrid_type /= 'adaptive') return
      ! extrapolate to the next step T%ty ==> targ.
      ! and update the y grid segments and readjust the 
      ! y grids.
      ! extrapolating, works good if there is not sudden jump. 
      targ = ZERO
      if(.not. present(opt)) then
         call use_lineint_multi(.false., x_next, T%tx, T%ty, T%nop, targ, '2p-interp')
      else if(index(opt, 'extrap') /= 0) then
         call use_lineint_multi(.false., x_next, T%tx, T%ty, T%nop, targ, '2p-interp')
      else if(index(opt, 'previous') /= 0) then
         targ = T%tau(:)
      else if(index(opt, 'jump') /= 0) then
         call use_lineint_multi(.false., x_next, T%tx, T%ty, T%nop, targ, '2p-interp')
         do i = 1, T%nop
            if(is_at_edge(targ(i)) == 1) then
               !targ(i) = 1.d-5
               targ(i) = 3.0_long * abs(targ(i) - ONE)
            else if(is_at_edge(targ(i)) == 0) then
               !targ(i) = ONE - 1.d-5
               targ(i) = 3.0_long * abs(targ(i))
            end if
         end do 
      end if 
      call update_ysegment(T, G, CTR, targ)
      call assign_ygrid_in_seg(G)
      return
   end subroutine update_ygrid   

   subroutine resonance(T, iclose, wi, mu_eff) 
      implicit none
      type(mqdt), intent(inout) :: T
      integer, intent(in) :: iclose
      real(long) :: mu_eff, wi
      real(long), dimension(T%nclose) :: width, mu
      call Kmat(T, width, mu)
      wi = width(iclose) 
      mu_eff = mu(iclose)
   end subroutine resonance

   subroutine push_nu_eff_in_range(T, G, S, CTR, XL, XR, iclose, mu_R)
      implicit none
      type(mqdt), intent(inout) :: T
      type(grid), intent(inout) :: G
      type(Smat), intent(in) :: S
      type(control), intent(in) :: CTR
      integer, intent(in) :: iclose
      real(long), intent(in) :: XL, XR
      real(long) :: xt, xt2, xx, width, mu_eff
      real(long), dimension(:), allocatable, intent(inout) :: mu_R
      integer :: i_actual_close
      i_actual_close = T%iclose(iclose)
      xx = ONE
      xt = ZERO
      width = zero
      do while(xt < XR)
         call interp_Smat(xx + 0.5_long, T, S, CTR)
         if(index(G%xgrid_type, 'mu-proj') /= 0) then ! seems not so promising.
            call resonance(T, iclose, width, mu_eff) 
            xt = xx + from_nuy_get_tau(mu_eff)
            !print *, from_nuy_get_tau(mu_eff), from_nuy_get_tau(T%mu(i_actual_close))
         else if(index(G%xgrid_type, 'mu-eigen') /= 0) then
            xt = xx + from_nuy_get_tau(T%mu(i_actual_close))
         end if
         width = width * TWO !Ryd. 
         width = width / T%Z**2 * xt**3
         !width = width / T%Z**2 
         !print *, 'width', width
         xt = xt - width
         xt2 = from_nu1_get_nu2(T, xt, T%IP_seq(i_actual_close), T%IPx)
         !write(*, "(2E20.10)", advance = 'no') xt, xt2
         if(is_ascend(XL, xt, XR)) then
            xt2 = from_nu1_get_nu2(T, xt, T%IP_seq(i_actual_close), T%IPx)
            call AddToList(G%xseg, xt2)
         !   write(*,"(' *')")
         !else 
         !   write(*,*)
         end if
         xt = xt + 2 * width
         xt2 = from_nu1_get_nu2(T, xt, T%IP_seq(i_actual_close), T%IPx)
         !write(*, "(2E20.10)", advance = 'no') xt, xt2
         if(is_ascend(XL, xt, XR)) then
            xt2 = from_nu1_get_nu2(T, xt, T%IP_seq(i_actual_close), T%IPx)
            call AddToList(G%xseg, xt2)
            call AddToList(mu_R, xt2)
         !   write(*,"(' *')")
         !else 
         !   write(*,*)
         end if
         xx = xx + ONE
      end do 
   end subroutine push_nu_eff_in_range

   subroutine x_init(T, S, G, CTR)
      implicit none
      type(mqdt), intent(inout) :: T
      type(Smat), intent(in) :: S
      type(grid), intent(inout) :: G
      type(control), intent(in) :: CTR
      integer :: i, iclose, i2
      real(long) :: XL, XR
      real(long) :: a, b, c, xt, xt2
      real(long), dimension(:), allocatable :: mu_R
      logical :: touched

      if(index(G%xgrid_type, "uniform") /= 0) then
         call uniform_xgrid(G)
         return 
      end if

      call AddToList(G%xseg, G%ei)
      call AddToList(G%xseg, G%ef)

      do i = 1, T%nclose 
         iclose = T%iclose(i)
         XL = from_nux_get_nu(T, G%ei, T%IP_seq(iclose))
         XR = from_nux_get_nu(T, G%ef, T%IP_seq(iclose))
         call push_nu_eff_in_range(T, G, S, CTR, XL, XR, findloc(T%iclose,iclose, dim=1), mu_R)
      end do 
      do i = 1, G%n_x_fine 
         xt = G%xdiv1(i)
         xt2 = G%xdiv2(i)
         if(is_ascend(G%ei, xt, G%ef)) call AddToList(G%xseg, xt)
         if(is_ascend(G%ei, xt2, G%ef)) call AddToList(G%xseg, xt2)
      end do 
      call sort_real_vec(G%xseg, 'ascend', G%xseg)
      G%n_seg_x = size(G%xseg)
      allocate(G%n_step_seg_x(G%n_seg_x))
      do i = 1, G%n_seg_x
         touched = .false.
         do i2 = 1, size(mu_R)
            if(G%xseg(i) == mu_R(i2)) then
               G%n_step_seg_x(i) = G%nx_spike
               touched = .true.
               exit 
            end if
         end do 
         b = G%xseg(i) 
         do i2 = 1, G%n_x_fine
            a = G%xdiv1(i2)
            c = G%xdiv2(i2)
            if(is_ascend(a, b,c ) .and. a /= b) then
               G%n_step_seg_x(i) = G%nx_flat * G%nxdiv(i2)
               touched = .true.
               exit
            end if
         end do 
         if(.not.touched) G%n_step_seg_x(i) = G%nx_flat
      end do 
      if(allocated(mu_R))deallocate(mu_R)
      G%nstep_x = sum(G%n_step_seg_x(:) ) - G%n_step_seg_x(1) + 1
      allocate(G%dx(G%nstep_x))
      call assign_xgrid_in_seg(G)

      call comment_str_among_char("X Grid", '-', funit_dbg)
      write(funit_dbg, "('# Initial guess for x grid segments:&
         & (',i5, 'segments ) ')") G%n_seg_x
      call print_vector_comment(G%xseg(:), G%n_seg_x, 15, '(f10.5)', funit_dbg)
      call print_int_vector_comment(G%n_step_seg_x(:), &
         G%n_seg_x, 15, '(i10)', funit_dbg)
      write(funit_dbg, "('# The number of x steps', i10)") G%nstep_x
      write(funit_dbg, "('# Later the x grids will be updated &
         &according to the slope of present solution.')")
      call comment_str_among_char("Grid", '-', funit_dbg)
      if(CTR%dbg_xgrid) then
         open(666,file = trim(CTR%path_mqdt)//"/xgrid.out", action = "write")
         write(666, "('# x grid segments.',i8)")G%n_seg_x
         do i = 1, G%n_seg_x
            !write(666, "(2E20.10)") ONE - G%xseg(i), ONE, ONE - G%xseg(i), ZERO
            write(666, "(2E20.10)") G%xseg(i), ONE, G%xseg(i), ZERO
            write(666, "('#')")
         end do 
         close(666)
      end if
      return 
   end subroutine x_init
end module gengrid
