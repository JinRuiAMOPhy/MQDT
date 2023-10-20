module numerical 
use constants
use type_def
contains
   function inner_prod(A, B) result(val)
      implicit none
      real(long), dimension(:), intent(in) :: A, B
      real(long) :: val
      integer :: i
      if(size(A) /= size(B)) stop
      val = ZERO
      do i = 1, size(A)
         val = val + A(i) * B(i)
      end do 
   end function inner_prod

   function matXvec(A, V) result(B)
      implicit none
      real(long), dimension(:,:), intent(in) :: A
      real(long), dimension(:), intent(in) :: V
      integer ::  nrow, ncol
      integer :: i
      real(long), dimension(size(V)) :: B
      nrow = size(A, 1)
      ncol = size(A, 2)
      if(ncol /= size(B)) stop 
      do i = 1, nrow
         B(i) = inner_prod(A(i,:), V)
      end do 
   end function matXvec

   function matXmat(A, B) result(C)
      implicit none
      real(long), dimension(:,:), intent(in) :: A, B
      real(long), dimension(size(A, 1),size(B, 2)) :: C
      integer ::  nr1, nc1, nr2, nc2
      integer :: i, j
      nr1 = size(A, 1)
      nc1 = size(A, 2)
      nr2 = size(B, 1)
      nc2 = size(B, 2)
      if(nc1 /= nr2) stop 
      do i = 1, nr1
         do j = 1, nc2
            C(i, j) = inner_prod(A(i,:), B(:,j))
         end do
      end do 
   end function matXmat
   function cinner_prod(A, B) result(val)
      implicit none
      complex(long), dimension(:), intent(in) :: A, B
      complex(long) :: val
      integer :: i
      if(size(A) /= size(B)) then
         stop 'inner_prod dimension ERROR'
      end if
      val = CMPLX_ZERO
      do i = 1, size(A)
         val = val + A(i) * B(i)
      end do
   end function cinner_prod

   function cmatXcvec(A, V) result(B)
      implicit none
      complex(long), dimension(:,:), intent(in) :: A
      complex(long), dimension(:), intent(in) :: V
      integer ::  nrow, ncol
      integer :: i
      complex(long), dimension(size(V)) :: B
      nrow = size(A, 1)
      ncol = size(A, 2)
      if(ncol /= size(B)) then
         stop 'matXvec dimension ERROR'
      end if
      do i = 1, nrow
         B(i) = cinner_prod(A(i,:), V)
      end do
   end function cmatXcvec

   function cmatXcmat(A, B) result(C)
      implicit none
      complex(long), dimension(:,:), intent(in) :: A, B
      complex(long), dimension(size(A, 1),size(B, 2)) :: C
      integer ::  nr1, nc1, nr2, nc2
      integer :: i, j
      nr1 = size(A, 1)
      nc1 = size(A, 2)
      nr2 = size(B, 1)
      nc2 = size(B, 2)
      if(nc1 /= nr2) then
         stop 'matXmat dimension ERROR'
      end if
      do i = 1, nr1
         do j = 1, nc2
            C(i, j) = cinner_prod(A(i,:), B(:,j))
         end do
      end do
   end function cmatXcmat

   function convert_energy_to_cm(Ein, Eunit) result(E_cm) 
      implicit none
      character(len = short_str), intent(in) :: Eunit
      real(long), intent(in) :: Ein
      real(long) :: E_cm
      E_cm = Ein
      if(index(Eunit, "Ry") /= 0 .or. &
         index(Eunit, "Ryd") /= 0 .or. &
         index(Eunit, "ry") /= 0 .or. &
         index(Eunit, "ryd") /= 0 ) then
         E_cm = Ein * Rydberg
      else if(index(Eunit, "au") /= 0 .or. &
         index(Eunit, "a.u.") /= 0 .or. &
         index(Eunit, "Hartree") /= 0 .or. &
         index(Eunit, "hartree") /= 0) then
         E_cm = TWO * Ein * Rydberg
      else if(index(Eunit, "eV") /= 0 .or. &
         index(Eunit, "ev") /= 0) then
         E_cm = Ein / AU_E_POTENTIAL * TWO * Rydberg
      else if(index(Eunit, "cm-1") /= 0 .or. trim(Eunit) == "") then
         E_cm = Ein 
      end if
   end function convert_energy_to_cm

   function from_E_cm_get_nu(T, E_cm, IP) result(nu)
      implicit none
      real(long), intent(in) :: E_cm
      integer, intent(in) :: IP
      type(mqdt), intent(in) :: T
      real(long) :: nu, dIP
      ! by default the E_ryd is wrt (ground state of N+1 system as 0.0, so first IP is positive)
      if(E_cm >= T%IPs(T%ipx)) then 
         write(STDERR, "('Caution: experiment level:', E20.10, &
            &' is greater than your abscissa IP ',E20.10, ' !!!')")E_cm, T%IPs(T%ipx)
         stop
      end if
      dIP = (T%IPs(IP) - E_cm) / Rydberg/ T%Z**2
      nu = ONE / sqrt(dIP)
   end function from_E_cm_get_nu

   function from_nux_get_nu(T, nux, IP) result(nu)
      implicit none
      real(long), intent(in) :: nux
      integer, intent(in) :: IP
      type(mqdt), intent(in) :: T
      real(long) :: nu, dIP
      dIP = (T%IPs(IP) - T%IPs(T%IPx)) / Rydberg/ T%Z**2 &
            + ONE / nux**2
      nu = ONE / sqrt(dIP)
   end function from_nux_get_nu

   function Energy_convert_Ryd_to_cm(E_ryd) result(E_cm)
       implicit none
       real(long), intent(in) :: E_ryd
       real(long) :: E_cm
       E_cm = E_ryd * Rydberg
   end function Energy_convert_Ryd_to_cm

   function Energy_convert_cm_to_Ryd(E_cm) result(E_Ryd)
       implicit none
       real(long), intent(in) :: E_cm
       real(long) :: E_Ryd
       E_Ryd = E_cm / Rydberg
   end function Energy_convert_cm_to_Ryd

   function Energy_convert_eV_to_cm(E_eV) result(E_cm)
       implicit none
       real(long), intent(in) :: E_eV
       real(long) :: E_cm
       E_cm = E_eV * Rydberg / AU_E_POTENTIAL
   end function Energy_convert_eV_to_cm

   function Energy_convert_au_to_cm(E_au) result(E_cm)
       implicit none
       real(long), intent(in) :: E_au
       real(long) :: E_cm
       E_cm = E_au * Rydberg * TWO
   end function Energy_convert_au_to_cm

   function get_E_ryd(T, nux) result(E_ryd)
      implicit none
      real(long) :: E_ryd
      real(long), intent(in) :: nux
      type(mqdt), intent(in) :: T
      E_ryd = T%IPs(T%ipx) / Rydberg - T%Z**2 / nux**2
   end function get_E_ryd

   function get_E_cont_ryd(T, nux) result(E_ryd)
      implicit none
      real(long) :: E_ryd
      real(long), intent(in) :: nux
      type(mqdt), intent(in) :: T
      E_ryd = (T%IPs(T%ipx) - T%IPs(1))/ Rydberg - T%Z**2 / nux**2
   end function get_E_cont_ryd

   function get_nuy(T, nux) result(nuy)
      implicit none
      type(mqdt), intent(in) :: T
      real(long), intent(in) :: nux
      real(long) :: nuy, dIP
      dIP = (T%IPs(T%ipy) - T%IPs(T%ipx)) / Rydberg/ T%Z**2
      nuy = ONE / sqrt(dIP + ONE / nux**2 )
   end function get_nuy

   function from_nuy_get_nux(T, nuy) result(nux)
      implicit none
      type(mqdt), intent(in) :: T
      real(long), intent(in) :: nuy
      real(long) :: nux, dIP
      dIP = (T%IPs(T%ipx) - T%IPs(T%ipy)) / Rydberg/ T%Z**2
      nux = ONE / sqrt(dIP + ONE / nuy**2 )
   end function from_nuy_get_nux

   function from_nuy_get_tau(nuy) result(tau)
      implicit none
      real(long), intent(in) :: nuy
      real(long) :: tau
      tau = mod(-nuy, ONE)
      if(tau < ZERO) then
         tau = tau + ONE
      end if
   end function from_nuy_get_tau

   function is_there_NaN(vec) result(inan)
      implicit none
      integer :: i, inan
      real(long), dimension(:) :: vec
      inan = 0
      do i = 1, size(vec)
         if(check_NaN(vec(i))) then
            inan = i
            return 
         endif
      end do 
      return
   end function is_there_NaN

   subroutine difference(nw, x, f, dfdx)
      implicit none
      integer :: nw, i
      real(long), dimension(:) :: x, f, dfdx
      do i = 1, nw - 1
         if(check_NaN(f(i+1)) .or. check_NaN(f(i))) then
           dfdx(i) = nan
         else
           dfdx(i) = (f(i+1) - f(i)) / (x(i+1) - x(i))
         endif
      end do 
      dfdx(nw) = dfdx(nw-1)
      return 
   end subroutine difference

   subroutine difference2(nw, x, f, dfdx, df2dx2)
      implicit none
      integer :: nw, i
      real(long), dimension(:) :: x, f, dfdx, df2dx2
      do i = 1, nw - 1
         dfdx(i) = (f(i+1) - f(i)) / (x(i+1) - x(i))
      end do 
      do i = 1, nw - 2
         df2dx2(i) = (dfdx(i+1) - dfdx(i)) / (x(i+1) - x(i))
      end do 
      dfdx(nw) = dfdx(nw-1)
      df2dx2(nw-1) = df2dx2(nw-2)
      df2dx2(nw) = df2dx2(nw-1)
      return 
   end subroutine difference2
! TODO : we should make use of the fact that the xs is always increasing,
!        so we can keep track of the current sample data used, instead of
!        always searching from the beginning.
   subroutine use_lineint_multi(ibug, x, Xs, Ys, ndim, Y)
      implicit none
      integer, intent(in) :: ndim
      real(long), dimension(:), intent(in) :: Xs
      real(long), dimension(:,:), intent(in) :: Ys
      real(long), dimension(ndim), intent(out) :: Y
      real(long), intent(in) :: x
      real(long), allocatable, dimension(:) :: y2
      real(long) :: yt
      logical, intent(in) :: ibug
      integer :: i
      allocate(y2(size(Xs)))
      do i = 1, ndim
         y2(:) = Ys(:,i)
         if(ibug)write(funit_dbg, '(5x,"ich", i3)', advance = 'no') i
         call use_lineint_single(ibug, x, Xs, y2, yt)
         if(ibug)write(funit_dbg, '(A)', advance = 'yes') 
         Y(i) = yt
      end do 
      deallocate(y2) 
   end subroutine use_lineint_multi

! TODO : we should make use of the fact that the xs is always increasing,
!        so we can keep track of the current sample data used, instead of
!        always searching from the beginning.
   subroutine use_lineint_single(ibug, x, Xs, Ys, Y)
      implicit none
      integer :: isampl, i, ns, nsampl
      real(long), dimension(:), intent(in) :: Xs
      real(long), dimension(:), intent(in) :: Ys
      real(long), intent(out) :: Y
      real(long), intent(in) :: x
      logical, intent(in) :: ibug
      real(long) :: SA, SB, A, B
      character(len = long_str), parameter :: fmtstr = '(1x,"Xs(",i3,") =", f8.4, &
            &" x = ", f8.4, " Xs(", i3,") =", f8.4, " Ys ", f8.4, " Ys ", f8.4)', &
      fmtstrWarn1 = '(1x, "WARNING Extraplating... x =", f8.4 , " Xs(1) =", f8.4)',&
      fmtstrWarn2 = '(1x, "WARNING Extraplating... x =", f8.4 , " Xs(", i3,") =", f8.4)',&
               fmtstrfit = '(" x = ", f8.4, " A ", f8.4, " B ", f8.4, " Y= ", f8.4)'
      nsampl = 2
      isampl = 1
      ns = size(Xs)
      if(x == Xs(1)) then
         Y = Ys(1)
         return  
      else if(x == Xs(ns)) then
         Y = Ys(ns)
         return
      else if(x < Xs(1)) then
         !if(ibug) write(STDERR, fmtstrWarn1) x, Xs(1)
         !if(ibug) write(funit_dbg, fmtstrWarn1) x, Xs(1)
         write(STDERR, fmtstrWarn1) x, Xs(1)
         write(funit_dbg, fmtstrWarn1) x, Xs(1)
         isampl = 2
      else if(x > Xs(ns)) then
         !if(ibug) write(STDERR, fmtstrWarn2) x, ns, Xs(ns)
         !if(ibug) write(funit_dbg, fmtstrWarn2) x, ns, Xs(ns)
         write(STDERR, fmtstrWarn2) x, ns, Xs(ns)
         write(funit_dbg, fmtstrWarn2) x, ns, Xs(ns)
         isampl = ns
      else 
         do i = 2, ns
           if(x >Xs(i-1) .and. x<= Xs(i)) then 
             isampl = i
             if(ibug)write(funit_dbg,fmtstr, advance = 'no') &
                 i-1, Xs(i-1), x, i, Xs(i), Ys(i-1), Ys(i)
             exit
           end if
         end do 
      end if
      ! linear fit between two points.
      call LINEFIT(Xs(isampl-1:isampl),Ys(isampl-1:isampl), A, B, nsampl, SA, SB) 
      Y = A * x + B         
      if(ibug) write(funit_dbg, fmtstrfit, advance = 'no') x, A, B, Y
      return 
   end subroutine use_lineint_single

!**LINEAR FIT*************************************
!  |---------------------------------------------|
!  |  LINEEAR FIT PROGRAM FOR EQUATION           |
!  |             Y=A*X+B                         |
!  |---------------------------------------------|
!  |  NOTICE: VARIABLE DIMENSION IS USED: N      |
!  |          THE CORRESPONDING ELEMENTS SHOULD  |
!  |          HAVE THE EXACTLY DIMENSION BETWEEN |
!  |          MAIN AND subroutine PROGRAM        |
!  |---------------------------------------------|
!   N: TOTAL POINT NUMBER
   subroutine LINEFIT(x,y,A,B,n,SA,SB)
      implicit none
      integer, intent(in) :: n
      real(long), intent(out) :: A, B, SA, SB
      real(long) :: sigX, sigY, sigXX, sigXY, sigVV, SY, CN, dN
      real(long), dimension(n), intent(in) :: X, Y
      sigX=sum(X(:))
      sigY=sum(Y(:))
      sigXX=sum(X(:)*X(:))
      sigXY=sum(X(:)*Y(:))
      dN = dble(N)

!   Cn: ={N*Sigma(Xi**2)-[Sigma(Xi)]**2}
      CN=dN*sigXX-sigX*sigX
!   A: =[N*Sigma(Xi*Yi)-Sigma(Xi)*Sigma(Yi)]/Cn
!   B: =[N*Sigma(Xi**2)-Sigma(Xi)*Sigma(Xi*Yi)]/Cn
      A=(dN*sigXY-sigX*sigY)/CN
      B=(sigXX*sigY-sigX*sigXY)/CN
!  sigVV = |Y - (Ax + B)|^2
      sigVV = sum((Y(:) - A * X(:) - B)**2)

!  SA: =SQRT[N/Cn]*SY
!  SB: =SQRT[Sigma(Xi**2)/Cn]*SY
!  SY=SQRT"Sigma(Vi**2)/(N-2)"
      SY=0
      if(N.GT.2) SY=sqrt(sigVV/(dN-2.))
      SA=sqrt(dN/CN)*SY
      SB=sqrt(sigXX/CN)*SY
      return
   end subroutine LINEFIT
!
!
!
   subroutine pack_angle(A, Euler, n) 
      implicit none
      integer, intent(in) :: N
      real(long), dimension(:), intent(in) :: A
      real(long), dimension(:, :), intent(out) :: Euler
      integer :: i, j, k
      k = 0
      do i = 1, N - 1
         do j = i + 1, N
            k = k + 1
            Euler(i, j) = A(k)
         end do
      end do
      return 
   end subroutine pack_angle
!

!**TRANSFORM BETWEEN ORTHOGONAL MATRIX AND EULER'S ANGLE**********
!    A: EULER'S ANGLE WITH SUBSCRIPTS (I,J);
!    U: ORTHOGONAL MATRIX U(n,n)
!    N: DIMENSION OF MATRIX
   subroutine ang2U(A, U, N)
      implicit none
      integer, intent(in) :: N
      real(long), dimension(n,n), intent(inout) :: U, A
      real(long), dimension(n,n) :: DS, DC
      integer :: i, j, m
      real(long) :: umi, umj
      call unitmat(U,N)
      do J=2,N
         do I=1,J-1
            DC(I,J)=cos(A(I,J))
            DS(I,J)=sin(A(I,J))
            do M=1,N
               UMI=U(M,I)
               UMJ=U(M,J)
               U(M,I)=UMI*DC(I,J)+UMJ*DS(I,J)
               U(M,J)=-UMI*DS(I,J)+UMJ*DC(I,J)
            end do 
         end do
      end do
   end subroutine ang2U


!    A: EULER'S ANGLE WITH SUBSCRIPTS (I,J);
!  IND: MATRIX-->ANGLE: IND=0; ANGLE-->MATRIX: IND=1
!    U: ORTHOGONAL MATRIX U(n,n)
!    K: NUMBER OF ANGLE;
!    N: DIMENSION OF MATRIX;
!    R: a copy of U
!    R: SHOULD BE UNIT MATRIX AFTER TERMINATING THE PROCESS OF MATRIX-->ANGLE;

   subroutine U2ang(A, U, N)
      implicit none
      real(long) :: T
      integer :: I, J, ia, m
      integer, intent(in) :: N
      real(long), dimension(n,n) :: U, A, DS, DC, R
      real(long) :: rmi, rmj
! initialize A, and copy U --> R
      A(:,:)  = 0.
      R(:,:) = U(:,:)
      loop_J: do J=N,2,-1 !30 
         T = ONE
!  CALCULATE A(1,J)...A(J-1,J) OF ORTHOGONAL MATRIX "do 10"
         loop_I: do I=1,J-1 !10
            DS(I,J)=R(J,I)/T
!cc Add the following one sent. xxp: 99/5/21 to ensure correct
            if(abs(ds(i,j))>1) ds(i,j) = ds(i,j)/abs(ds(i,j))
            A(I,J)=asin(DS(I,J))
            DC(I,J)=cos(A(I,J))
            T=T*DC(I,J)
!--IN CASE SOME ANGLE A(I,J)=90 DEGREE--|
            if(abs(DC(I,J)).LT.1e-10_long.and.I.LT.J-1) then
               WRITE(funit_dbg,*) 'A(',I,J,')-90 DEGREE<1e-10'
               do IA=I+1,J-1
                  A(IA,J)=0
                  DS(IA,J)=0
                  DC(IA,J)=1
               end do
               exit loop_I
            end if
        end do loop_I !10
!--IN CASE THE DIAGONAL ELEMENT OF THE ORTHOGONAL MATRIX IS NEGATIVE--|
        if(R(J,J).LT.0) then
           A(J-1,J)=PI-A(J-1,J)
           if(A(J-1,J).GT.PI) A(J-1,J)=A(J-1,J)-PI*TWO
           DC(J-1,J)=-DC(J-1,J)
        end if
!  ELIMINATE A(1,J)...A(J-1,J) FROM ORTHOGONAL MATRIX "do 20"
         do I=J-1,1,-1
            do M=1,N
              RMI=R(M,I)
              RMJ=R(M,J)
              R(M,I)=RMI*DC(I,J)-RMJ*DS(I,J)
              R(M,J)=RMI*DS(I,J)+RMJ*DC(I,J)
            end do 
         end do
      end do  loop_J
      return
   end subroutine U2ang

   subroutine Uangcon(U,A,N,IND)
!  |---------------------------------------------|
!  |  NOTICE: VARIABLE DIMENSION IS USED: N      |
!  |          THE CORRESPONDING ELEMENTS SHOULD  |
!  |          HAVE THE EXACTLY DIMENSION BETWEEN |
!  |          MAIN AND subroutine PROGRAM        |
!  |---------------------------------------------|
!    A: EULER'S ANGLE WITH SUBSCRIPTS (I,J);
!  IND: MATRIX-->ANGLE: IND=0; ANGLE-->MATRIX: IND=1
!    U: ORTHOGONAL MATRIX U(n,n)
!    N: DIMENSION OF MATRIX;
      implicit none
      integer, intent(in) ::ind, n
      real(long), dimension(n,n) :: U, A

      if(n.EQ.1) then
         U(1,1)=1
      else if(IND.EQ.1) then
         call ang2U(A, U, N)
      else
         call U2ang(A, U, N)
      end if
      return 
   end subroutine Uangcon
!**UNIT MATRIX WITH Nth DIMENSION************************************
   subroutine unitmat(um, n)
      implicit none
      integer, intent(in) :: n
      integer :: i, j
      real(long), dimension(n, n), intent(inout) :: UM
      do i = 1, n
         do j = 1, n
            um(i,j) = ZERO
         end do 
         um(i,i) = ONE
      end do 
   end subroutine unitmat
!-------------------------------------------------------------
! CAL. DETERMINANT OF RECTANGULAR MATRIX u BY TRANSFORMING IT
! INTO A TRIAGNULAR MATRIX u' WITH GAUSSIAN ELIMINATION (G)
! u' = G*u
! det{u} = det{u'} = \Prod_i{u'_{ii}}
! copy u before operation, so no need to copy before passing it
! to this routine
   subroutine MINV(u_in, N, w)
      implicit none
      real(long), dimension(:,:), intent(inout) :: u_in
      real(long), intent(out) :: w
      real(long) :: t, x
      integer, intent(in) :: N
      integer :: i, j, l, k, m
      real(long), allocatable :: u(:,:)
      allocate(u(N,N))
      u(:,:) = u_in(1:N, 1:N)
      w = ONE
      l = 1
      i = 1
      loop_row: do while(i < n .and. w /= 0) 
         if(u(i,i) < 1.e-32_long) then
         j = i + 1
         loop_pivot: do while( j <= n .and. l > 0)  
            if(abs(u(j,i)) >= 1.e-32_long) then
               do k = i, n
                  t = u(i,k)
                  u(i,k) = u(j,k)
                  u(j,k) = t
               end do 
               w = -w
               l = 0
            end if
            j = j + 1
         end do loop_pivot
         if(l /= 0) w = 0 ! all zero
         end if
         if(w /= ZERO) then
            do m = i + 1, n
               x = u(m,i) / u(i,i)
               do l = i + 1, n
                  u(m,l)=u(m,l)-x*u(i,l)
               end do 
            end do
         end if
         i=i+1
      end do loop_row

      if( w /= ZERO) then
         do i = 1, n
            w = w * u(i,i)
         end do
      end if
      deallocate(u)
      return
   end subroutine MINV

!**PRODUCE MINOR MATRIX OF RECTANGULAR MATRIX************************
   subroutine MINOR(F,FMIN,I,J,N)
      implicit none
      real(long), dimension(:,:), intent(in) :: F
      real(long), dimension(:,:), intent(out) :: FMIN
      integer, intent(in) :: I, J, N
      integer :: L, K, L1, K1
      do L = 1, N - 1
        do K = 1, N - 1
           if(L.GE.I) then
              L1 = L + 1
           else
              L1 = L
           endif
           if(K.GE.J) then
              K1 = K + 1
           else
              K1 = K
           end if
           FMIN(L,K) = F(L1,K1)
        end do 
     end do
     return 
   end subroutine MINOR

   subroutine sort_real_vec(f, decend_or_ascend, g2, order2)
      implicit none
      real(long), dimension(:), intent(in) :: f
      real(long), dimension(size(f)), intent(inout) :: g2
      integer, dimension(:), intent(inout), optional :: order2
      real(long), dimension(size(f)) :: g
      integer, dimension(size(f)) :: order
      real(long) :: tmp
      character(len = *) :: decend_or_ascend
      integer :: n, i, j, it
      n = size(f)
      g(:) = f(:)
      order(:) = (/(i , i = 1, n)/)
      do i = 1, n
         do j = 1, n
            if(g(i) > g(j)) then
               tmp = g(j)
               it = order(j)
               g(j) = g(i)
               order(j) = order(i)
               g(i) = tmp
               order(i) = it
            endif
         end do
      end do
      if(index(decend_or_ascend, 'decend') /= 0) then
         g2(:) = g(:)
         if(present(order2))order2(:) = order(:)
      else if(index(decend_or_ascend, 'ascend') /= 0) then
         do i=1,n
            g2(i) = g(n-i+1)
            if(present(order2))order2(i) = order(n-i+1)
         enddo
      end if
   end subroutine sort_real_vec
   function cumsum(V) result(G)
      implicit none
      real(long), dimension(:), intent(in) :: V
      real(long), dimension(size(V)) :: G
      integer :: i
      do i = 1, size(V)
         G(i) = sum(V(1:i))
      end do
   end function cumsum
   function is_ascend(a, b, c) result(res)
      implicit none
      real(long), intent(in) :: a, b, c
      logical :: res
      if(a<=b .and. b <= c) then
         res = .true.
      else
         res = .false.
      end if
   end function is_ascend
   function is_array_ascend(a, n) result(res)
      implicit none
      real(long), dimension(:), intent(in) :: a
      integer, intent(in) :: n
      logical :: res
      integer :: i
      res = .true.
      do i = 2, n
         if(a(i) < a(i-1)) then
            res = .false.
            exit
         end if
      end do

   end function is_array_ascend

   function is_window_center_a_peak(a, force) result(res)
   ! the window width is 2 * RES_WINDOW -1
      implicit none
      real(long), dimension(:), intent(in) :: a
      integer :: imax
      character(len=*), intent(in), optional :: force
      logical :: res
      if(present(force) .and. index(force, 'force') /= 0) then
         res = .true.
         return 
      else 
         imax =  maxloc(a, dim =1)
         if(is_array_ascend(a(1:RES_WINDOW), RES_WINDOW) .and. & 
            is_array_ascend(a(RES_WINDOW * 2 - 1 : RES_WINDOW:-1), RES_WINDOW)) then
            res = .true.
         else if(RES_WINDOW == imax) then
            res = .true.
         else 
            res = .false.
         end if
      end if
   end function is_window_center_a_peak

!***********************************************
!* SOLVING A COMPLEX LINEAR MATRIX SYSTEM AX=B *
!* with Gauss-Jordan method using full pivoting*
!* at each step. During the process, original  *
!* A and B matrices are destroyed to spare     *
!* storage location.                           *
!* ------------------------------------------- *
!* INPUTS:    A    COMPLEX MATRIX N*N          *
!*            B    COMPLEX MATRIX N*M          *
!* ------------------------------------------- *
!* OUTPUTS:   A    INVERSE OF A N*N            *
!*            DET  COMPLEX DETERMINANT OF A    *
!*            B    SOLUTION MATRIX N*M         *
!* ------------------------------------------- *
!* NOTA - If M=0 inversion of A matrix only.   *
!***********************************************
   function CUMAT(N) result(U)
      implicit none
      integer, intent(in) :: N
      complex, dimension(N, N) :: U
      integer :: i
      U(:,:) = CMPLX_ZERO
      do i = 1, n
         U(i, i) = REAL_ONE
      end do 
   end function
   subroutine CMATINV(N,M,AA,BB,det)
      implicit none
      integer, intent(in) :: N, M
      real(long), PARAMETER :: EPSMACH=2.E-12_long
      complex(long), dimension(N, N), intent(inout) :: AA
      complex(long), dimension(N, M), intent(inout) :: BB
      complex(long), intent(out) :: det
      INTEGER, pointer :: PC(:), PL(:)
      complex(long), pointer :: CS(:)
      complex(long) :: PV, TT
      real(long) :: PAV
      integer :: i, j, k, ik, jk, ialloc
                                                    
!Initializations :                       
      allocate(PC(1:N),stat=ialloc)
      allocate(PL(1:N),stat=ialloc)
      allocate(CS(1:N),stat=ialloc)

      det=REAL_ONE
      do I=1,N
         PC(I)=0
         PL(I)=0
         CS(I)=CMPLX_ZERO
      end do
!main loop                                                                   
      do K=1,N                                                                  
!Searching greatest pivot:                                               
        PV=AA(K,K)                                                              
        IK=K                                                                    
        JK=K                                                                    
        PAV=ABS(PV)                                                            
        do I=K,N                                                                
          do J=K,N                                                              
            if (ABS(AA(I,J)).GT.PAV) then                                      
              PV=AA(I,J)                                                        
              PAV=ABS(PV)                                                      
              IK=I                                                              
              JK=J                                                              
            endif                                                               
          enddo                                                                 
        enddo                                                                   
                                                                               
!Search terminated, the pivot is in location I=IK, J=JK.
!Memorizing pivot location: :                                        
        PC(K)=JK                                                                
        PL(K)=IK                                                                
                                                                               
!Determinant det is actualised
!If det=0, ERROR MESSAGE and STOP
!Machine dependant EPSMACH equals here 2.E-12                                        
                                                                               
        if (IK.NE.K) det=-det                                                   
        if (JK.NE.K) det=-det                                                   
        det=det*PV                                                              
        if (ABS(det).LT.EPSMACH) &
           stop 'DETERMINANT EQUALS ZERO, NO SOLUTION!'
!POSITIONNING PIVOT IN K,K:                                              
        if(IK.NE.K) then                                                        
          do I=1,N                                                              
!EXCHANGE LINES IK and K:                                            
            TT=AA(IK,I)                                                         
            AA(IK,I)=AA(K,I)                                                    
            AA(K,I)=TT                                                          
          enddo                                                                 
        endif                                                                   
      if (M.NE.0) then                                                          
        do I=1,M                                                               
          TT=BB(IK,I)                                                           
          BB(IK,I)=BB(K,I)                                                      
          BB(K,I)=TT                                                            
        enddo                                                                   
      endif                                                                     
!Pivot is at correct line                                                
        if(JK.NE.K) then                                                        
          do I=1,N                                                              
!Exchange columns JK and K of matrix AA                                         
            TT=AA(I,JK)                                                         
            AA(I,JK)=AA(I,K)                                                    
            AA(I,K)=TT                                                          
          enddo                                                                 
        endif                                                                   
!Pivot is at correct column and located in K,K                                              
!Store column K in vector CS then set column K to zero                                             
        do I=1,N                                                                
          CS(I)=AA(I,K)                                                         
          AA(I,K)=CMPLX(0.0,0.0)                                                          
        enddo                                                                   
!                                                                               
        CS(K)=CMPLX(0.0,0.0)                                                                
        AA(K,K)=CMPLX(1.0,0.0)                                                              
!Modify line K :                                            
        if(ABS(PV).LT.EPSMACH) &
           stop '  PIVOT TOO SMALL - STOP'
        do I=1,N                                                                
          AA(K,I)=AA(K,I)/PV                                                    
        enddo                                                                   
        if (M.NE.0) then                                                        
          do I=1,M                                                             
            BB(K,I)=BB(K,I)/PV                                                  
          enddo                                                                 
        endif                                                                   
!Modify other lines of matrix AA:                                        
        do J=1,N                                                                
          if (J.EQ.K) continue                                                  
          do I=1,N                                                              
!Modify line J of matrix AA :                                            
            AA(J,I)=AA(J,I)-CS(J)*AA(K,I)                                       
          enddo                                                                 
          if (M.NE.0) then                                                      
            do I=1,M                                                          
              BB(J,I)=BB(J,I)-CS(J)*BB(K,I)                                     
            enddo                                                               
          endif                                                                 
        enddo                                                                   
!Line K is ready.                                                
      enddo                                                                     
!End of K loop                                                              
!The matrix AA is inverted - Rearrange AA                         
!Exchange lines                                                            
      do I=N,1,-1                                                               
        IK=PC(I)                                                                
        if (IK.EQ.I) continue                                                   
!EXCHANGE LINES I AND PC(I) OF AA:                                         
        do J=1,N                                                                
          TT=AA(I,J)                                                            
          AA(I,J)=AA(IK,J)                                                      
          AA(IK,J)=TT                                                           
        enddo                                                                   
        if (M.NE.0) then                                                        
          do J=1,M                                                             
            TT=BB(I,J)                                                          
            BB(I,J)=BB(IK,J)                                                    
            BB(IK,J)=TT                                                         
          enddo                                                                 
        endif                                                                   
!NO MORE EXCHANGE NEEDED                                                      
!GO TO NEXT LINE                                                  
      enddo                                                                     
                                                                               
!EXCHANGE COLUMNS                                                          
      do J=N,1,-1                                                               
        JK=PL(J)                                                                
        if (JK.EQ.J) continue                                                   
!EXCHANGE COLUMNS J AND PL(J) OF AA :                                       
        do I=1,N                                                                
          TT=AA(I,J)                                                            
          AA(I,J)=AA(I,JK)                                                      
          AA(I,JK)=TT                                                           
        enddo                                                                   
!NO MORE EXCHANGE NEEDED                                                      
!GO TO NEXT COLUMN   
      enddo                                                                     
!REARRANGEMENT TERMINATED.                                                        
      return 
   end subroutine CMATINV                                                                      

!*******************************************
   subroutine Kmat(T, mu, U, width)
      implicit none
      type(mqdt), intent(in) :: T
      real(long), dimension(T%nchan,T%nchan), intent(in) :: U
      real(long), dimension(T%nchan), intent(in) :: mu
      complex(long), dimension(T%nchan, T%nchan) :: K
      complex(long), dimension(T%nop, T%nop) :: Koo, invmat, CUnit
      complex(long), dimension(T%nop, T%nclose) :: Koc
      complex(long), dimension(T%nclose, T%nop) :: Kco
      complex(long), dimension(T%nclose, T%nclose) :: Kcc, Keff, Kr
      real(long), dimension(T%nclose),intent(out) :: width 
      complex(long) :: det
      integer :: i, j, ia
      K(:,:) = CMPLX_ZERO
      do i = 1, T%nchan
         do j = 1, T%nchan
            do ia = 1, T%nchan
              K(i, j) = K(i, j) + U(i,ia) * tan(PI * mu(ia)) * U(j, ia)
            end do 
         end do 
      end do 
      Koo(:,:) = K(1:T%nop, 1:T%nop)
      CUnit = CUMAT(T%nop)
      invmat(:,:) = Koo(:,:) + CUnit(:,:) * (ZERO, ONE)
      call CMATINV(T%nop, T%nop, invmat, CUnit, det)
      Kcc(:,:) = K(T%nop+1:T%nchan, T%nop+1:T%nchan)
      Koc(:,:) = K(1:T%nop, T%nop+1:T%nchan)
      Kco(:,:) = K(T%nop+1:T%nchan, 1:T%nop)
      Keff(:,:) = Kcc(:,:) - cmatXcmat(Kco, cmatXcmat(invmat, Koc))
      Kr(:,:) = real(Keff(:,:))
      do i = 1, T%nclose
         width(i) = AIMAG(Keff(i, i)) / PI * TWO  
      end do 
      return 
   end subroutine Kmat
end module numerical
