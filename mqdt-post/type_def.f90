module type_def
   use constants
   type mqdt
      integer :: n_e
      integer :: n_e_sort
      integer :: nop, nclose
      integer :: nchan
      integer :: n_ang
      integer :: n_IP
      integer :: twoJ
      integer :: ipy, ipx ! previously ip1, ip2 remember we use higher ip to define x
      integer :: k_mat_cofact
      integer, allocatable :: ip_seq(:), iopen(:), iclose(:)
      character(len = short_str) :: IPunit
      real(long) :: Z
      ! nux previously nu2  
      real(long), allocatable :: E(:), nux(:), &
                                 IPs(:), tau(:,:), Ara(:,:,:), &
                                 Anorm(:,:,:), Dn(:,:), os(:,:), &
                                 os_dos(:,:), tau2(:,:), os2(:,:),&
                                 nux_sort(:), norm_dos(:,:)
   end type mqdt

   type fname
      character(len = long_str + short_str) :: f_in_tau, f_in_Ara,&
                       f_in_Smat, f_in_mqdt, f_in_dfde, f_exp_tau,&
                       f_exp_dfde, f_in_dmat, f_in_Anorm, f_in_Dn,&
                       f_in_XDn, f_o_os, f_o_factor, f_o_os_dos,  &
                       f_o_sig, f_o_tau_sort, f_o_os_sort, f_o_phys, &
                       f_o_energyrel, f_o_exp_plot, f_o_theo_plot,&
                       f_o_phys_compare, f_dbg, f_in_normdos, f_dbg_dtdE
   end type fname
   type control
      !integer :: n_e
      real(long) :: egnd_istate, filter
      character(len = long_str) :: path_mqdt
      character(len = short_str) :: proc 
      character(len = short_str) :: figname
      character(len = short_str) :: couple
      integer :: ibug
      logical :: yes_ignore
      logical :: yes_eigenchansort
      logical :: yes_dfde
      logical :: yes_dis
      logical :: yes_au
      logical :: yes_user_filter
      real(long) :: omega
      logical :: dbg_os, dbg_sort, &
                 dbg_interp, dbg_levfind, &
                 dbg_deriv
      character(len = short_str) :: eqnsolv_method, dbgstr,&
                 discrete_method, au_peak_search_method
   end type control
   type Smat
      integer :: n_es, nchan, n_ang
      real(long), dimension(:,:), allocatable :: ang_s, mu_s, dmat_s
      real(long), dimension(:), allocatable :: es
      character(len = short_str), dimension(:), allocatable :: spec_jj, spec_ls
      real(long), allocatable :: d_mu(:)
   end type Smat
   type physobserv
      integer :: nlev
      real(long), dimension(:), allocatable :: lev_nux, lev_nuy, lev_ryd, width,&
                  nux_plot(:), nuy_plot(:)
      integer, dimension(:), allocatable :: eigenchan_id, lev_indx, nlev_per_chan
      character(len = short_str), dimension(:), allocatable :: elev_spec, mark
   end type physobserv

   type spectrodata
      integer :: nlev
      real(long), dimension(:), allocatable :: lev_nux, lev_nuy, energy, width
      integer, dimension(:), allocatable :: eigenchan_id, lev_indx
      character(len = short_str), dimension(:), allocatable :: elev_spec
   end type spectrodata
   type grid
      real(long), dimension(:), allocatable :: dx, dy
      real(long) :: ei, ef, E_cont_i, E_cont_f, &
                    ygrid_seg_jump
      character(len = short_str) :: xgrid_type, ygrid_type
      integer :: nstep_x, ny_adapt, ny_init_guess,&
                 nx_flat, nx_spike
      integer :: nstep_y, n_x_fine, n_y_fine
      integer, dimension(:), allocatable :: nxdiv, nydiv
      real(long), dimension(:), allocatable :: xdiv1, xdiv2, &
                                ydiv1, ydiv2
      real(long), allocatable :: xseg(:), yseg(:)
      integer, allocatable :: n_step_seg_x(:), n_step_seg_y(:)
      integer :: n_seg_y, n_seg_x
   end type grid   
end module type_def
