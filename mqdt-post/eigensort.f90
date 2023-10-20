module eigensort
!*************************************
!*  2021 Jin, Rui, CFEL
!*  modified from color.v0.f
!*************************************      
use constants
use type_def
use darray
use stdio
contains
   subroutine sort_mqdt(T, CTR, FN)
      implicit none
      type(mqdt), intent(inout) :: T
      type(control), intent(in) :: CTR
      type(fname), intent(in):: FN
      integer :: ich, ie, nbreaks
      integer, parameter :: funit = 1, funit2 = 2
      integer, dimension(T%nchan) :: MAXNCH, MAXNCH_old
      integer, dimension(:), allocatable :: breakpoint
      integer, dimension(:,:), allocatable :: which_break
      real(long), dimension(T%n_e,T%nchan) :: tau2, os2
      real(long), dimension(T%nchan) :: slope_old
      character(len = short_str) :: fmtstr = '(E20.10)'
      logical :: is_discontinue, is_now_NOT_nan, is_previous_NOT_nan
      write(stderr, *) 'entering sort_mqdt'
      T%n_e_sort = 0
      slope_old(:) = zero !keep the slope
      tau2(:,:) = nan
      os2(:,:) = nan
      nbreaks = 0 
      MAXNCH_old(:) = 0
      loop_E: do ie = 1, T%n_e
         call sort_tau_os(ie, CTR, T, tau2(ie,:), os2(ie,:), MAXNCH)
         if(ie > 2) then
            call break_line(T, tau2, ie, MAXNCH, MAXNCH_old, slope_old, &
                 nbreaks, breakpoint, which_break, CTR%dbg_sort)
         else
            slope_old(:) = (/( (tau2(2,ich) - tau2(1,ich))/(T%nux(2) - T%nux(1)), ich = 1, T%nchan )/)
         end if
         MAXNCH_old(:) = MAXNCH(:)
      end do loop_E
      call push_to_end(T, tau2, os2, nbreaks, breakpoint, which_break)
      if(nbreaks>0 .and. CTR%dbg_sort)then
         write(funit_dbg,"('break points at ')", advance = 'no')
         call print_int_vector(breakpoint, nbreaks, nbreaks, '(i8)', funit_dbg)
      end if
!    COLOR_end

      open(funit, file=FN%f_o_tau_sort, action = 'write')
      if(CTR%yes_dfde) &
         open(funit2, file=FN%f_o_os_sort, action = 'write')
      do ie = 1, T%n_e_sort
         !is_discontinue = .False.
         !if (ie > 2) then
         !  do ich = 1, T%nchan
         !  is_now_NOT_nan = .not. check_NaN(tau2(ie, ich))
         !  is_previous_NOT_Nan = .not.check_NaN(tau2(ie - 1, ich))
         !  end do 
         !  is_discontinue = (is_now_NOT_nan .and. (.not.is_previous_NOT_Nan) ) .or. &
         !                   (is_previous_NOT_Nan .and. (.not.is_now_NOT_nan))
         !end if
         !if(is_discontinue) then
         !  write(funit,"(E20.10)", advance = 'no') nan
         !  call print_vector(T%tau2(ie,:), T%nchan, T%nchan, fmtstr, funit)
         !  if(CTR%yes_dfde) then
         !     write(funit,"(E20.10)", advance = 'no') nan
         !     call print_vector(T%os2(ie,:), T%nchan, T%nchan, fmtstr, funit)
         !  end if
         !end if
         write(funit,"(E20.10)", advance = 'no') T%nux_sort(ie)
         call print_vector(T%tau2(ie,:), T%nchan, T%nchan, fmtstr, funit)
         if(CTR%yes_dfde) then
            write(funit2,"(E20.10)", advance = 'no') T%nux_sort(ie)
            call print_vector(T%os2(ie,:), T%nchan, T%nchan, fmtstr, funit2)
         end if
      end do 
      if(CTR%yes_dfde) &
         close(funit2)
      close(funit)
      if(allocated(breakpoint)) deallocate(breakpoint)
      if(allocated(which_break)) deallocate(which_break)
   end subroutine sort_mqdt

   subroutine break_line(T, tau2, ne_now, MAXNCH, chanord, slope_old, nbreaks, breakpoint, which_break, ibug)
      implicit none
      type(mqdt), intent(inout) :: T
      integer :: ich, ne_now
      integer, intent(inout) :: nbreaks
      integer, dimension(T%nchan), intent(in) :: MAXNCH, chanord
      integer, dimension(:), allocatable :: breakpoint
      integer, dimension(:,:), allocatable :: which_break
      integer, dimension(T%nchan) :: discontinue_chan
      real(long), dimension(T%n_e,T%nchan) :: tau2
      real(long), dimension(T%nchan), intent(inout) :: slope_old
      real(long), dimension(T%nchan) :: slope
      real(long) :: dfun, dx
      logical :: yes_swap, do_add_nan, is_now_NOT_nan, is_previous_NOT_nan
      logical, intent(in) :: ibug
      yes_swap = .false.
      is_now_NOT_nan = .false.
      is_previous_NOT_nan = .false.
      ! when  two channels swapped, and both of them are not NaN, accidental 
      ! connection happens, we definitely want to avoid this
      discontinue_chan = 0
      do ich = 1, T%nchan
         is_now_NOT_nan = .not. check_NaN(tau2(ne_now, ich))
         is_previous_NOT_Nan = .not.check_NaN(tau2(ne_now - 1, ich))
         !if(chanord(ich) /= MAXNCH(ich) .and.is_now_NOT_nan.and.is_previous_NOT_Nan ) then
         if(chanord(ich) /= MAXNCH(ich))then
             yes_swap = .true.
             discontinue_chan(ich) = 1
             if(ibug)write(funit_dbg, "('swap detected at', i8,' for ', i3, 'th &
             &chan, we then check its continuity')", advance = 'no') ne_now, ich 
             exit
         end if
      end do 
      dx = T%nux(ne_now) - T%nux(ne_now - 1)  
      do_add_nan = .false.
      do ich = 1, T%nchan
         dfun = tau2(ne_now,ich) - tau2(ne_now - 1, ich)
         slope(ich) = dfun / dx
         if(abs(dfun) > 0.1) then
            do_add_nan = .true.
            discontinue_chan(ich) = 1
            if(ibug)write(funit_dbg,"(i3, E20.10,i8, ' chan dicontinuity &
            &f(i)-f(i-1)', E20.10)") ich,T%nux(ne_now),ne_now,dfun
         !else if(abs(slope(ich)) > 10.0 * abs(slope_old(ich)) .and. slope(ich) * slope_old(ich) < 0.d0) then
         else if(.true. .and. slope(ich) * slope_old(ich) < 0.d0) then  ! sudden jump. one pi.
            do_add_nan = .true.
            discontinue_chan(ich) = 1
            if(ibug)&
              write(funit_dbg,"(i3, E20.10, i8, ' chan slope discontinuity df/dx(i) ', E20.10,' &
              &df/dx(i-1) ', E20.10)") ich, T%nux(ne_now), ne_now,slope(ich), slope_old(ich)
         end if
         slope_old(ich) = slope(ich)
      end do 
      !if(do_add_nan .and. yes_swap) then
      if(do_add_nan .or. yes_swap) then
         nbreaks = nbreaks + 1
         call Add_int_ToList(breakpoint, ne_now)
         call Add_int_list_To2D(which_break, discontinue_chan)
      end if
   end subroutine break_line

   subroutine sort_tau_os(ie, CTR, T, tau2, os2, MAXNCH) 
      implicit none
      type(mqdt), intent(inout) :: T
      type(control), intent(in) :: CTR
      integer :: ich, iop, ie
      integer, dimension(T%nchan), intent(out) :: MAXNCH
      real(long), dimension(T%nchan) :: MAXA
      real(long), dimension(T%nchan), intent(out) :: tau2, os2
      logical, dimension(T%nchan) :: visited
   
      visited(:) = .false.
      loop_rho: do iop = 1, T%nop
         MAXNCH(iop)=0
         MAXA(iop)=zero
         loop_ch1: do ich = 1, T%nchan
            if( (abs(T%Anorm(ie,ich,iop)) .GT. abs(MAXA(iop))) .and. &
                (.not. visited(ich)) )then
                MAXNCH(iop)=ich
                MAXA(iop)=abs(T%Anorm(ie, ich, iop))
             endif
         end do loop_ch1
         visited(MAXNCH(iop)) = .true.
!    tau_COLOR           
         loop_ch2: do ich = 1, T%nchan
            if(ich .EQ. MAXNCH(iop)) tau2(ich)=T%tau(ie,iop)
         end do loop_ch2
!   DFDE_COLOR
         if(CTR%yes_dfde)then
            loop_chan3: do ich = 1, T%nchan
               if (ich .EQ. MAXNCH(iop)) os2(ich)=T%os(ie,iop)
            end do loop_chan3
         end if
      end do loop_rho
      if(CTR%dbg_sort) then
         write(funit_dbg,'(9i8)')(iop, iop =1, T%nop)
         do ich = 1, T%nchan
            write(funit_dbg,'(i4,9f8.4)')ich, (T%Anorm(ie, ich, iop), iop = 1, T%nop)
         end do 
         write(funit_dbg,'(4x, 9i8)')(MAXNCH(iop), iop = 1, T%nop)
         write(funit_dbg,'(9L8)')visited(:)
      end if
   end subroutine sort_tau_os

   subroutine push_to_end(T, tau2, os2, nbreaks, breakpoint, which_break)
      implicit none
      type(mqdt), intent(inout) :: T
      integer :: ie, ie2, ibreak, iroam_break, ich
      integer, intent(in) :: nbreaks
      integer, dimension(:), allocatable, intent(in) :: breakpoint
      integer, dimension(:,:), allocatable, intent(in) :: which_break
      real(long), dimension(:,:), intent(in) :: tau2, os2
      real(long), dimension(T%nchan)::  array_nan
      array_nan(:) = nan
      if(nbreaks == 0) then
         T%n_e_sort = T%n_e
         allocate(T%nux_sort(T%n_e_sort), T%tau2(T%n_e_sort,T%nchan),T%os2(T%n_e_sort,T%nchan))
         T%nux_sort(:) = T%nux(:)
         T%tau2(:,:) = tau2(:,:)
         T%os2(:,:) = os2(:,:)
      else
         T%n_e_sort = T%n_e + nbreaks
         allocate(T%nux_sort(T%n_e_sort), T%tau2(T%n_e_sort,T%nchan),T%os2(T%n_e_sort,T%nchan))
         ie2 = 0
         iroam_break = 1
         do ie = 1, T%n_e
            ie2 = ie2 + 1
            loop_break: do ibreak = 1, nbreaks
               if(ie == breakpoint(ibreak)) then
                  !T%nux_sort(ie2) = nan
                  T%nux_sort(ie2) = (T%nux(ie-1) + T%nux(ie) ) / TWO
                  !T%tau2(ie2,:) = array_nan(:)
                  !T%os2(ie2,:) = array_nan(:)
                  do ich = 1, T%nchan
                     if(which_break(ibreak,ich) == 1) then
                        T%tau2(ie2,ich) = nan
                        T%os2(ie2,ich) = nan
                        !print "('intert >>> ', i8,' ich', i3, ' tau(i-1) ',f8.4,' tau_m ', &
                        !&f8.4,' tau(i) ', f8.4)",ie, ich,&
                        !tau2(ie-1,ich), T%tau2(ie2,ich), tau2(ie,ich)
                     else 
                        T%tau2(ie2,ich) = (tau2(ie-1,ich) + tau2(ie,ich)) / TWO
                        T%os2(ie2,ich) = (os2(ie-1,ich) + os2(ie,ich)) / TWO
                        !print "('intert >>> ', i8,' ich', i3, ' tau(i-1) ',f8.4,' tau_m ', &
                        !&f8.4,' tau(i) ', f8.4)",ie, ich,&
                        !tau2(ie-1,ich), T%tau2(ie2,ich), tau2(ie,ich)
                     end if
                  end do 
                  ie2 = ie2 + 1
                  exit loop_break
               end if
            end do loop_break
            iroam_break = ibreak
            T%nux_sort(ie2) = T%nux(ie)
            T%tau2(ie2,:) = tau2(ie,:)
            T%os2(ie2,:) = os2(ie,:)
         end do
      end if
      return 
   end subroutine push_to_end
end module eigensort
