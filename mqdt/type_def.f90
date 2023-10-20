module type_def
   use constants
   type mqdt
      integer :: n_e
      integer :: nop, nclose
      integer :: nchan
      integer :: n_ang
      integer :: n_IP
      real(long), allocatable :: IPs(:)
      character(short_str) :: IPunit
      integer, allocatable :: ip_seq(:), iopen(:), iclose(:)
      integer :: twoJ
      integer :: ipy, ipx ! previously ip1, ip2 remember we use higher ip to define x
      integer :: ip3
      integer :: iperiod
      integer :: k_mat_cofact
      real(long) :: Z
      real(long) :: nux
      real(long) :: E_abs ! IPx - Z**2/ nux**2
      real(long), allocatable :: tau(:), Ara(:,:), Anorm(:,:), Dn(:), norm_dos(:)
      real(long), allocatable :: mu(:), Angle(:)
      real(long), allocatable :: U(:,:), F(:,:), Tir(:,:), Euler(:,:)
      real(long), allocatable :: ty(:,:), deriv(:,:)
      real(long), dimension(2) :: tx
      real(long) :: detF
    ! for resonance analysis
      complex(long), allocatable :: K(:,:), Kcc(:,:), Koc(:,:), Kco(:,:), &
                 Koo(:,:), K_eff(:,:)
   end type mqdt

   type fname
      character(len = long_str + short_str) :: f_tau, f_Ara,&
                       f_Smat, f_mqdt, f_Anorm, f_Dn, f_dbg, f_normdos, &
                       f_para_use
      integer :: funit_tau, funit_Ara, funit_Anorm, funit_Dn, funit_normdos, &
                 funit_para_use
   end type fname

   type control
      real(long) :: filter
      character(len = long_str) :: path_mqdt
      character(len = short_str) :: couple
      integer :: ibug
      logical :: yes_dis, yes_au
      logical :: dbg_F, dbg_xgrid, dbg_ygrid, dbg_solu, &
                 dbg_interp, dbg_Tir, dbg_A, dbg_eqnsolv, &
                 dbg_para_use
      real(long) :: omega
      character(len = short_str) :: eqnsolv_method, dbgstr,&
                 discrete_method
   end type control

   type Smat
      integer :: n_es, nchan, n_ang
      real(long), dimension(:,:), allocatable :: ang_s, mu_s
      real(long), dimension(:), allocatable :: es
      real(long), dimension(:), allocatable :: mu
      real(long), allocatable :: d_mu(:)
      character(len = short_str), dimension(:), allocatable :: spec_jj, spec_ls
   end type Smat

   type eqnsolv
      real(long) :: x_l, x_r, y_l, y_r, x_thresh, y_thresh
      integer :: counter
   end type eqnsolv

   type grid
      real(long), dimension(:), allocatable :: dx, dy
      real(long) :: ei, ef, E_cont_i, E_cont_f
      real(long) :: ygrid_seg_jump
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
