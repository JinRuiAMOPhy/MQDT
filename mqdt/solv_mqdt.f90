! DISELGF.FOR
! SIMILAR TO DISELUA.FOR, BUT WITH VARY STEP OF ENERGY [DEL(NIU)=CONST.]
! PROGRAM No.: MQDT-07HW

!  ALL RIGHTS RESERVED BY GROUP 309, PHYSICS INSTITUTE, ACADEMIA SINICA
!                  (VERSION 2.6, SEPTEMBER 17, 1994)

!    CALCULATION THE DISCRETE ENERGY LELELS and g-FACTORS ACCORDING U MATRIX
!!   (OR EULER'S ANGLES), QUANTUM DEFECTS MIU AND jj-LS TRANSFORMATION MATRIX
!  |-----------------------------------------------------------------|
!  | NOTICE: THE 'NDIMEN' VALUES(=40) IN MAIN AND SUBROUTINE PROGRAMS|
!  |         SHOULD GREATER then THE REAL CHANNEL No.;               |
!  |         THE FIRST CHANNEL SHOULD BE THE LOWEST CHANNEL WHICH    |
!  |         SERVES AS THE REFERENCE CHANNEL FOR CAL. 'A' AND 'MIU', |
!  |         AS WELL AS THE VARY ENERGY STEP.                        |
!  |-----------------------------------------------------------------|
!
!  Rui Jin modifies into f90
!  replace common blocks with module data structures.
!  

program SOLVMQDT
   use constants
   use type_def
   use numerical
   use file_cmdl_io
   use stdio
   use envset
   use gengrid
   implicit none
   type(fname) :: FN
   type(mqdt) :: T
   type(grid) :: G
   type(Smat) :: S
   type(control) :: CTR
   type(eqnsolv) :: ES
   integer :: counter
   real(long) :: start, finish
   call cpu_time(start)
   call environmentset() ! setting nan et. al.
   call initiate(FN, T, S, CTR, G)
   if(T%nop == T%nchan) then
    write(funit_dbg, "('#------------------')")
    write(funit_dbg, "('# nop == nchan, no need to project')")
    call allchan_open(FN, T, S, CTR, G)
   else
    if(index(G%ygrid_type, "uniform") /= 0) then
       call graphic_solv_Yuniform(FN, T, S, CTR, G)
    else 
       call graphic_solv_Yadaptive(FN, T, S, CTR, G)
    end if
    print *, 'solved'
   end if
   call cpu_time(finish)
   write(stdout, "('took',e20.10,' sec.')") finish - start
    write(funit_dbg, "('#------------------')")
   write(funit_dbg, "('#took ',A,' sec.')") trim(float2str(finish - start, '(f9.3)'))
   call finalize(T, S, G)
contains
   subroutine allchan_open(FN, T, S, CTR, G)
      implicit none
      type(fname), intent(in) :: FN
      type(mqdt), intent(inout) :: T
      type(grid), intent(inout) :: G
      type(Smat), intent(in) :: S
      type(control), intent(in) :: CTR
      integer :: istep, nrt, i
      character(len = short_str) :: fmtstr
      fmtstr = '(2E25.15)'

      open(FN%funit_tau, file = FN%f_tau, action = 'write') 
      open(FN%funit_Ara, file = FN%f_Ara, action = 'write') 
      open(FN%funit_Anorm, file = FN%f_Anorm, action = 'write') 
      open(FN%funit_Dn, file = FN%f_Dn, action = 'write') 
      if(CTR%yes_dis)&
         open(FN%funit_normdos, file = FN%f_normdos, action = 'write') 

      ! calculate on a rather dense y grid based on mu and some manual segments
      ! if provided.
      ! As an initiative guess for the adaptive y grid segments.
      loop_x: do istep = 1, S%n_es 
        write(FN%funit_tau, fmtstr, advance = 'no') S%es(istep)
        call print_vector(S%mu_s(istep,:), T%nop, T%nop, fmtstr, FN%funit_tau)
        write(FN%funit_Ara, fmtstr, advance = 'no') S%es(istep)
        T%Ara(:,:) = ZERO; T%Anorm(:,:) = ZERO; T%Dn = ONE; T%Norm_dos = ONE
        do i = 1, T%nchan
          T%Ara(i,i) = ONE; T%Anorm(i,i) = ONE
        end do 
        call print_matrix_as_vect_by_col(T%Ara, T%nchan, T%nop, T%nchan,&
             fmtstr, FN%funit_Ara)
        write(FN%funit_Anorm, fmtstr, advance = 'no') S%es(istep)
        call print_matrix_as_vect_by_col(T%Anorm, T%nchan, T%nop, T%nchan,&
             fmtstr, FN%funit_Anorm)
        write(FN%funit_Dn, fmtstr, advance = 'no') S%es(istep)
        call print_vector(T%Dn, T%nop, T%nop, fmtstr, FN%funit_Dn)
        if(CTR%yes_dis) then
           write(FN%funit_normdos, fmtstr, advance = 'no') S%es(istep)
           call print_vector(T%Norm_dos, T%nop, T%nop, fmtstr, FN%funit_normdos)
        end if
    !       stop 'missing root(s). resonance(maybe not because of the degeneracy'
      end do loop_x

      close(FN%funit_tau); close(FN%funit_Ara); close(FN%funit_Dn)
      close(FN%funit_Anorm)
      if(CTR%yes_dis) close(FN%funit_normdos)
   end subroutine allchan_open

   subroutine graphic_solv_Yadaptive(FN, T, S, CTR, G)
      implicit none
      type(fname), intent(in) :: FN
      type(mqdt), intent(inout) :: T
      type(grid), intent(inout) :: G
      type(Smat), intent(in) :: S
      type(control), intent(in) :: CTR
      integer :: istep, nrt, i
      real(long) :: x_next
      character(len = short_str) :: fmtstr
      fmtstr = '(2E25.15)'

      open(FN%funit_tau, file = FN%f_tau, action = 'write') 
      open(FN%funit_Ara, file = FN%f_Ara, action = 'write') 
      open(FN%funit_Anorm, file = FN%f_Anorm, action = 'write') 
      open(FN%funit_Dn, file = FN%f_Dn, action = 'write') 
      if(CTR%yes_dis)&
         open(FN%funit_normdos, file = FN%f_normdos, action = 'write') 
      if(CTR%dbg_para_use) &
         open(FN%funit_para_use, file = FN%f_para_use, action = 'write')
      T%nux = G%ei
      T%E_abs = get_E_ryd(T, T%nux)
      ! calculate on a rather dense y grid based on mu and some manual segments
      ! if provided.
      ! As an initiative guess for the adaptive y grid segments.
      call x_init(T, S, G, CTR) 
      call y_init(T, S, G, CTR, 'guess') 
      call search(T%nux, T, S, CTR, G, nrt)
      if(nrt == T%nop) then
         call print_result(T, FN)
      else 
         stop 'Initiative guess grid failed.'
      end if
      ! then adopt previous solution as better guess.
      ! reallocate the y grid segment and redefine grid 
      ! number per segment
      call y_init(T, S, G, CTR, 'reinitiate')

      T%tx(2) = T%nux
      T%ty(2,:) = T%tau(:) ! for calculating derivative istep = 1
      counter = 0
      if(CTR%dbg_para_use) then
         write(FN%funit_para_use, fmtstr, advance='no') T%nux
         call print_vector(T%mu, T%nchan, T%nchan, fmtstr, FN%funit_para_use, 'no')
         call print_vector(T%Angle, T%n_ang, T%n_ang, fmtstr, FN%funit_para_use, 'yes')
      end if
      loop_x: do istep = 2, G%nstep_x - 1
         ! interpolated S matrix are in T data structure.
         if(G%dx(istep) <= ZERO) cycle ! to avoid the repeated grid/segments adaptive grid scheme
         T%nux = T%nux + G%dx(istep)
         T%tx(1) = T%tx(2)
         T%tx(2) = T%nux
         T%ty(1,:) = T%ty(2,:)
         T%ty(2,:) = T%tau(:)
         ! avoid sudden jump by one.
         if(.false.) then
         do i = 1, T%nop
            if(T%ty(1, i) - T%ty(2, i) < -0.3) then
               T%ty(1, i) = T%ty(1, i) + ONE
            else if(T%ty(1, i) - T%ty(2, i) > 0.3) then
               T%ty(1, i) = T%ty(1, i) - ONE
            end if
         end do 
         end if
         !
         T%E_abs = get_E_ryd(T, T%nux)
         call search(T%nux, T, S, CTR, G, nrt)
         if(nrt == T%nop) then
            call print_result(T, FN)
         else if(nrt == 0) then
            print *, "nrt == 0 "
            deallocate(G%yseg, G%n_step_seg_y, G%dy)
            call y_init(T, S, G, CTR, 'guess')
            call search(T%nux, T, S, CTR, G, nrt)
            if(nrt == T%nop) then
               call print_result(T, FN)
               call y_init(T, S, G, CTR, 'reinitiate')
               T%ty(2,:) = T%tau(:) ! for calculating derivative istep = 1
               print *, " re-initiate y grid succeed"
            else
               print *, nrt, "roots found"
               write(funit_dbg, "('try re-initiate y grid failed.',i5)")nrt
               !stop 'try re-initiative y grid failed.'
               print *, 'try re-initiate y grid failed.'
               call y_init(T, S, G, CTR, 'reinitiate-wider')
               call search(T%nux, T, S, CTR, G, nrt)
               if(nrt == T%nop) then
                 print *, " reinitiate-wider y grid succeed"
               else
                 print *, 'try reinitiate-wider y grid failed.'
                 cycle
               end if
               !continue
            end if
         else if(nrt /= T%nop) then ! recalculate using previous results as guessing y segment.
            write(stderr, "('Recalculate based on previous tau:')", advance = 'yes')
            write(funit_dbg, "('#Recalculate based on previous tau:')", advance = 'yes')
            call update_ygrid(T%nux, T, G, CTR, 'previous')
            write(stderr,*)G%yseg(:)
            call print_vector_comment(G%yseg, G%n_seg_y, G%n_seg_y, '(f20.10)', funit_dbg, 'yes')
            call search(T%nux, T, S, CTR, G, nrt)
            if(nrt == T%nop) then
               call print_result(T, FN)
               write(stderr, "('succeed')")
               write(funit_dbg, "('#succeed')")
            else
               write(stderr, "('still fail >>> jumping because of (mod 1)')")
               write(stderr, "('goto the other side with guessed 3*tau')")
               write(funit_dbg, "('#still fail >>> jumping because of (mod 1)')")
               write(funit_dbg, "('#goto the other side with guessed 3*tau')")
               call update_ygrid(T%nux, T, G, CTR, 'jump')
               write(stderr,*)G%yseg(:)
               call print_vector_comment(G%yseg, G%n_seg_y, G%n_seg_y, '(f20.10)', funit_dbg, 'yes')
               call search(T%nux, T, S, CTR, G, nrt)
               if(nrt == T%nop) then
                  call print_result(T, FN)
                  write(stderr, "('succeed')")
                  write(funit_dbg, "('#succeed')")
               else
                  write(stderr, "('still still fail, move on to next x point')")
               !   write(stderr, "('if this point is important try specify &
               !      &suitable ygrid_seg_jump through mqdt.in')")
                  write(funit_dbg, "('#still still fail, move on to next x point')")
               !   write(funit_dbg, "('#if this point is important try specify &
               !      &suitable ygrid_seg_jump through mqdt.in')")
                  !stop
               end if
            end if
         end if
         ! updating the y grid by extrapolation 
         ! works good if there is not sudden jump. 
         x_next = T%nux + G%dx(istep + 1)
         call update_ygrid(x_next, T, G, CTR, 'extrap')

         if(index(G%xgrid_type, 'adaptive') /= 0) then
            ! check slope or extraction. then decide go back and refine
            !
         end if
         if(CTR%dbg_para_use) then
            write(FN%funit_para_use, fmtstr, advance='no') T%nux
            call print_vector(T%mu, T%nchan, T%nchan, fmtstr, FN%funit_para_use, 'no')
            call print_vector(T%Angle, T%n_ang, T%n_ang, fmtstr, FN%funit_para_use, 'yes')
         end if
      end do loop_x

      close(FN%funit_tau); close(FN%funit_Ara); close(FN%funit_Dn)
      close(FN%funit_Anorm)
      if(CTR%yes_dis) close(FN%funit_normdos)
      if(CTR%dbg_para_use) close(FN%funit_para_use)
   end subroutine graphic_solv_Yadaptive

   subroutine graphic_solv_Yuniform(FN, T, S, CTR, G)
      implicit none
      type(fname), intent(in) :: FN
      type(mqdt), intent(inout) :: T
      type(grid), intent(inout) :: G
      type(Smat), intent(in) :: S
      type(control), intent(in) :: CTR
      integer :: istep, nrt
      real(long) :: x_next

      open(FN%funit_tau, file = FN%f_tau, action = 'write') 
      open(FN%funit_Ara, file = FN%f_Ara, action = 'write') 
      open(FN%funit_Anorm, file = FN%f_Anorm, action = 'write') 
      open(FN%funit_Dn, file = FN%f_Dn, action = 'write') 
      if(CTR%yes_dis)&
         open(FN%funit_normdos, file = FN%f_normdos, action = 'write') 

      T%nux = G%ei
      T%E_abs = get_E_ryd(T, T%nux)
      ! calculate on a rather dense y grid based on mu and some manual segments
      ! if provided.
      ! As an initiative guess for the adaptive y grid segments.
      call x_init(T, S, G, CTR) 
      call y_init(T, S, G, CTR, 'uniform') 

      loop_x: do istep = 1, G%nstep_x - 1
         ! interpolated S matrix are in T data structure.
         if(G%dx(istep) <= ZERO) cycle ! to avoid the repeated grid/segments adaptive grid scheme
         T%nux = T%nux + G%dx(istep)
         T%E_abs = get_E_ryd(T, T%nux)
         call search(T%nux, T, S, CTR, G, nrt)
         if(nrt == T%nop) then
            call print_result(T, FN)
         else if(nrt < T%nop) then ! recalculate with previous results as guessing y segment.
            write(stderr, "('miss root')", advance = 'yes')
            write(funit_dbg, "('miss root')", advance = 'yes')
         end if
         x_next = T%nux + G%dx(istep + 1)

         if(index(G%xgrid_type, 'adaptive') /= 0) then
            ! check slope or extraction. then decide go back and refine
            !
         end if
      end do loop_x

      close(FN%funit_tau); close(FN%funit_Ara); close(FN%funit_Dn)
      close(FN%funit_Anorm)
      if(CTR%yes_dis) close(FN%funit_normdos)
   end subroutine graphic_solv_Yuniform

   subroutine print_result(T, FN)
      type(fname), intent(in) :: FN
      type(mqdt), intent(in) :: T
      character(len = short_str) :: fmtstr
      fmtstr = '(2E25.15)'
      write(FN%funit_tau, fmtstr, advance = 'no') T%nux, T%E_abs
      call print_vector(T%tau, T%nop, T%nop, fmtstr, FN%funit_tau)
      write(FN%funit_Ara, fmtstr, advance = 'no') T%nux, T%E_abs
      call print_matrix_as_vect_by_col(T%Ara, T%nchan, T%nop, T%nchan,&
           fmtstr, FN%funit_Ara)
      write(FN%funit_Anorm, fmtstr, advance = 'no') T%nux, T%E_abs
      call print_matrix_as_vect_by_col(T%Anorm, T%nchan, T%nop, T%nchan,&
           fmtstr, FN%funit_Anorm)
      write(FN%funit_Dn, fmtstr, advance = 'no') T%nux, T%E_abs
      call print_vector(T%Dn, T%nop, T%nop, fmtstr, FN%funit_Dn)
      if(CTR%yes_dis) then
         write(FN%funit_normdos, fmtstr, advance = 'no') T%nux, T%E_abs
         call print_vector(T%Norm_dos, T%nop, T%nop, fmtstr, FN%funit_normdos)
      end if
    !     stop 'missing root(s). resonance(maybe not because of the degeneracy'
   end subroutine print_result

   subroutine search(nux, T, S, CTR, G, nrt)
      implicit none
      type(mqdt), intent(inout) :: T
      type(Smat), intent(in) :: S
      type(control), intent(in) :: CTR
      type(grid), intent(in) :: G
      real(long),intent(in) :: nux
      integer, intent(out) :: nrt
      real(long) :: nuy, detF
      integer :: iter, istep, irt
      real(long), parameter :: nuy_thresh = THRESH, &
                               det_thresh = THRESH
      real(long), dimension(2) :: res
      real(long), dimension(T%nop +5 ,2) :: cross, y_cross
      real(long) :: ytmp
      nuy = ZERO 
      ytmp  = nuy
      istep = 1
      nrt = 0
      call interp_Smat(nux, T, S, CTR)
      call prepare_Fmat(nux, nuy, T, CTR)
      detF = T%detF
      do istep = 2, size(G%dy) 
         if(G%dy(istep) <= ZERO) cycle ! to avoid the repeated grid/segments adaptive grid scheme
         call prepare_Fmat(nux, nuy, T, CTR)
         !if(CTR%dbg_eqnsolv)print *, nuy, T%detF
         if(T%detF * detF <= ZERO) then
            nrt = nrt + 1
            cross(nrt, 1) = ytmp
            cross(nrt, 2) = nuy
            y_cross(nrt, 1) = detF
            y_cross(nrt, 2) = T%detF
            if(nrt == T%nop) exit
         end if 
         detF = T%detF
         ytmp = nuy
         nuy = nuy + G%dy(istep)
      end do 
      if(nrt == 0) then
         write(*,*)'nrt == 0 (no roots located).'
         stop
      end if
      if( nrt /= T%nop ) then
         write(funit_dbg, "('# reinitializing Y grids')")
         write(funit_dbg, "('# nux = ', A, 1x,A, ' roots found != nop (',i3,')')") &
            trim(float2str(nux,'(E16.10)')), trim(int2str(nrt,'(i8)')), T%nop
         write(funit_dbg, "('#', 2E20.10)") (ONE - cross(irt, 1), y_cross(irt, 1), &
            ONE - cross(irt, 2), y_cross(irt, 2), irt = 1, nrt)
         write(STDERR, "('# nux = ', E20.10, i3, ' roots found != nop')") nux, nrt
         nrt = 0  ! don't plot, because it's very hard to analyze.
         ! later we will refine y.
         ! The reason for this failure is most probably due to the fact that as tau reaching 
         ! 1, the modulo 1 will lead to discontinuity. 
         ! TODO ; can we get rid of the modulo 1, and solve by a moving window y
         ! range
         ! Or, we can solve with modulo 1, but with slightly wider y-range : 
         ! [0,1] => [0-e,1+e], difficulty is how to provide a universal e.
         write(funit_dbg, "('# ygrid segments')")
         call print_vector_comment(G%yseg, G%n_seg_y, G%n_seg_y, '(f20.10)', funit_dbg, 'yes')
         call print_vector_comment(ONE - G%yseg, G%n_seg_y, G%n_seg_y, '(f20.10)', funit_dbg, 'yes')
         !stop
         return 
      end if
      T%tau(:) = nan
      do irt = 1, nrt
         if(CTR%dbg_eqnsolv) write(funit_dbg, "(A, ' root located between [', A,',',A, '] finer search:')") &
              trim(adjustl(ordinal(irt))), trim(float2str(cross(irt,1),'(E16.10)')), &
              trim(float2str(cross(irt,2),'(E16.10)'))
         ES = create_nonlinear_solver(cross(irt,1), cross(irt,2), &
              y_cross(irt, 1), y_cross(irt, 2), nuy_thresh, det_thresh)
         counter = counter + 1
         if(CTR%eqnsolv_method == 'newton') then
            call newton_relax(T, ES, nux, res, iter, CTR)
         else if(CTR%eqnsolv_method == 'bisect') then
            call bisect(T, ES, nux, res, iter, CTR)
         else if(CTR%eqnsolv_method == 'hyb') then
            call hyb(T, ES, nux, res, iter, CTR)
         end if
! iter < 0 : exceeding maximum iterations
         if(iter < 0) then
            T%tau(irt) = nan
            write(STDERR,"('|F| = 0 reach MAXITER for nux = ', A,1x, A, ' root')")&
                 trim(float2str(nux,'(E16.10)')), trim(adjustl(ordinal(irt)))
            write(funit_dbg,"('# |F| = 0 reach MAXITER for nux =', A,1x, A, ' root')")&
                 trim(float2str(nux,'(E16.10)')), trim(adjustl(ordinal(irt)))
            write(funit_dbg, "('#', 2E20.10)") cross(irt, 1), y_cross(irt, 1), &
               cross(irt, 2), y_cross(irt, 2)
            write(funit_dbg, "('# stuck at', 2E20.10)") res(1), res(2)
            call check_detF(nux,cross(irt, 1),cross(irt, 2) , T, CTR, irt)
            call SOLU(T, CTR, irt)
            nrt = 0 
            return 
         else 
            T%tau(irt) = mod(ONE - res(1), ONE )
            call SOLU(T, CTR, irt)
         end if
      end do 
      if(CTR%dbg_eqnsolv) then
         write(funit_dbg, "('# eigenchannel mu:')")
         call print_vector(T%mu, T%nchan, print_col_s_fmt, '(f9.4)', funit_dbg)
         write(funit_dbg, "('# roots:')")
         call print_vector(T%tau(:), T%nop, print_col_s_fmt, '(f9.4)', funit_dbg)
      end if
   end subroutine search

   subroutine check_detF(nux, yL, yR, T, CTR, irt)
      implicit none      
      type(mqdt), intent(inout) :: T
      type(control), intent(in) :: CTR
      real(long) :: nux, yL, yR, y
      character(len = short_str) :: fn
      integer :: irt
      write(fn, "(i3)") irt
      open(777, file = trim(CTR%path_mqdt)//'/f-test.'//trim(adjustl(fn)))
      y = TWO*yL - yR
      do 
         call prepare_Fmat(nux, y, T, CTR)
         write(777, '(2E20.10)')y, T%detF
         if(y > TWO * yR - yL) exit
         y = y + (3.0_long * (yR - yL) / 50.0_long)
      end do 
      close(777)
   end subroutine check_detF   

! SOLVE THE F EQUATION: F*A = 0
! for any row i: sum_j(Fij * Aj) = 0
! alse |F| = sum_j( Fij * (-1)^(i,j) * cofactor_ij ) = 0 
! so Aj is proportional to (-1)^(i,j) * cofactor_ij
!-------------------------------------------------------------------------
!--SEARCH ith MINOR MATRIXS OF F WITH THE LARGEST SUM OF DETERMINANT**2--|
! now choose the largest determinant are removed, because different row of 
! coffactor may lead to different signs, leading to flipping Ara phases.
! Most of the time, k_mat_cofact can be fixed as any row.
   subroutine SOLU(T, CTR, irho)
      implicit none
      type(mqdt), intent(inout) :: T
      type(control), intent(in) :: CTR
      real(long) :: fdet2, detF
      real(long), dimension(T%nchan-1,T%nchan-1) ::FMIN
      integer, intent(in) :: irho
      integer :: j
      fdet2 = ZERO
      do J = 1, T%nchan
         CALL MINOR(T%F, FMIN, T%k_mat_cofact, j, T%nchan)
         CALL MINV(FMIN, T%nchan - 1, detF)
         fdet2 = fdet2 + detF * detF
         T%Ara(j, irho) = (-ONE)**(T%k_mat_cofact + j) * detF
      end do
      if(CTR%dbg_solu) &
         write(funit_dbg,'("ρ = ", A, " |M(k,ρ)| :")')trim(int2str(irho,'(i8)'))
      if(CTR%dbg_solu) &
         call print_vector(T%Ara(:, irho), T%nchan, print_col_s_fmt, '(E13.5)', funit_dbg) 
      if(CTR%dbg_solu) &
         write(funit_dbg,'("# √ ∑|M(k,ρ)|^2 =", A)')trim(float2str(sqrt(fdet2),'(E16.10)'))
      if(CTR%dbg_solu) &
         write(funit_dbg,'("#   k")')
      T%Ara(:, irho) = T%Ara(:, irho) / sqrt(fdet2)
!
!  T matrix 
! |open      | ==> effective U matrix
! |close     |
      call TIJ(T, CTR, irho)
   end subroutine SOLU
end program SOLVMQDT
