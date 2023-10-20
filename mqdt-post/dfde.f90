module dfde
   use constants
   use darray
   use type_def
   use numerical
   use stdio
contains
   subroutine dfn(T, S, FN, CTR)
      implicit none
      type(mqdt), intent(inout):: T
      type(Smat), intent(inout) :: S
      type(fname), intent(in) :: FN
      type(control), intent(in) :: CTR
      integer :: iop, ie
      integer, parameter :: funit1 = 3, funit2 = 4, funit3 = 45, funit4 = 300
      real(long) :: E_ph, Omega, E_chan_ele, sig
      real(long), dimension(T%nop) :: os, DDD
      real(long), dimension(T%nchan) :: dalfa
      character(len = short_str) :: fmtstr = '(E20.10)'

      open(funit1, file = FN%f_o_os)
      open(funit2, file = FN%f_o_sig)
      write(funit1,'("# E_photon  sig*Z^4")')
      write(funit2,'("# nux os")')
      if(CTR%yes_dis) then
         open(funit4, file = FN%f_o_os_dos)
         allocate(T%os_dos(T%n_e, T%nop))
      end if
      allocate(T%os(T%n_e, T%nop))
      do ie = 1, T%n_e
         E_ph = T%E(ie)-CTR%egnd_istate ! T%E (Ryd) wrt IP1, not deduced.
         OMEGA = E_ph
         E_chan_ele = (T%E(ie) -T%IPs(1)/Rydberg )/ T%Z**2 ! reduced
         if(CTR%dbg_interp) write(funit_dbg, "('fit dmat: ', i4, ' E_chan_ele(reduced) =',&
               &f12.4, ' Ryd., E_ph =', f12.4 )")ie, E_chan_ele, E_ph
         call use_lineint_multi(CTR%dbg_interp, E_chan_ele, S%es, S%dmat_s, T%nchan, dalfa)
         if(CTR%dbg_interp) write(funit_dbg, '(200E20.10)')E_chan_ele, (dalfa(iop), iop = 1, T%nchan)
         !stop
         sig = 0.0
      ! reduced energy Ryd. shift with IP(1)
         loop_rho: do iop = 1, T%nop 
            DDD(iop) = sum(dalfa(:) * T%Anorm(ie,:, iop))
            OS(iop) = DDD(iop)**2 * T%Z**4
            sig = sig + DDD(iop) ** 2
         end do loop_rho
         sig = Omega * sig
         os(:) = os(:) * Omega
         T%os(ie, :) = os(:)
         write(funit1,'(2E20.10)')E_ph, sig * T%Z**4
         write(funit2,"(E20.10)", advance = 'no') T%nux(ie)
         call print_vector(os, T%nop, T%nop, fmtstr, funit2)
         if(CTR%yes_dis) then
            T%os_dos(ie, :) = os(:)/T%norm_dos(ie, :)
            write(funit4,"(E20.10)",advance = 'no') T%nux(ie)
            call print_vector(T%os_dos(ie, :), &
              T%nop, T%nop, fmtstr, funit4)
         end if
      end do
      close(funit1); close(funit2); 
      if(CTR%yes_dis) close(funit4)
   end subroutine dfn
   subroutine dfde_open(T, S, FN, CTR, G)
      implicit none
      type(mqdt), intent(inout):: T
      type(Smat), intent(inout) :: S
      type(fname), intent(in) :: FN
      type(grid), intent(in) :: G
      type(control), intent(in) :: CTR
      integer :: iop, ie, i
      integer, parameter :: funit1 = 3, funit2 = 4, funit3 = 45, funit4 = 300
      real(long) :: E_ph, Omega, E_chan_ele, sig
      real(long), dimension(T%nop) :: os, DDD
      real(long), dimension(T%nchan) :: dalfa
      character(len = short_str) :: fmtstr = '(E20.10)'

      open(funit1, file = FN%f_o_os)
      open(funit2, file = FN%f_o_sig)
      write(funit1,'("# E_photon  sig*Z^4")')
      write(funit2,'("# nux os")')
      do ie = 1, S%n_es
         E_ph = S%es(ie)-CTR%egnd_istate ! T%E (Ryd) wrt IP1, not deduced.
         OMEGA = E_ph
         E_chan_ele = (S%es(ie) - T%IPs(1)/Rydberg )/ T%Z**2 ! reduced
         if((G%E_cont_i /T%Z**2/ Rydberg > S%es(ie)) .or. &
              (S%es(ie) > G%E_cont_f/T%Z**2/ Rydberg)) cycle
         if(CTR%dbg_interp) write(funit_dbg, "('fit dmat: ', i4, ' E_chan_ele(reduced) =',&
               &f12.4, ' Ryd. >> E_ph =', f12.4 )")ie, E_chan_ele, E_ph
         dalfa(:) = S%dmat_s(ie, :)
         if(CTR%dbg_interp) write(funit_dbg, '(200E20.10)')E_chan_ele, (dalfa(iop), iop = 1, T%nchan)
         !stop
         sig = ZERO
      ! reduced energy Ryd. shift with IP(1)
         loop_a: do i = 1, T%nchan 
           sig = sig + dalfa(i) ** 2
         end do loop_a
         sig = sig * Omega
         write(funit1,'(2E20.10)')E_ph, sig * T%Z**4
         write(funit2,"(E20.10)", advance = 'no') T%nux(ie)
         call print_vector(os, T%nop, T%nop, fmtstr, funit2)
      end do
      close(funit1); close(funit2); 
   end subroutine dfde_open
end module dfde
