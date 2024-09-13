! lande g-factor
   subroutine lande(Ara, g, twoJ, N)
      implicit none
      real(long), dimension(:), intent(in) :: Ara
      integer, intent(in) :: N, twoJ
      real(long), intent(out)  :: g
      real(long), dimension(N) :: gaa, gls, Ukk
      real(long), dimension(N, N) :: Uia, Ula, Uka, Ujl, Ujk, &
          Ukl, Uls, Ujj
      integer, dimension(N) :: Nl, Ns
      integer :: i, j

      G = ZERO
      loop_row: do i = 1, N
         UlS(i) = ZERO
         Ujj(i) = ZERO
         Ukk(i)=0.
         loop_col: do j = 1, N
            UlA(i,j) = ZERO
            UkA(i,j) = ZERO
            loop_inter: do k = 1, N
               Uka(i,j) = Uka(i,j) + Uia(k,j) * Ujk(k,i)
               Ula(i,j)=Ula(i,j) + Uia(k,j) * Ujl(k,i)
            end do inter
            Ujj(i)=Ujj(i)+UiA(i,j)*Ara(j)
            Ukk(i)=Ukk(i)+UkA(i,j)*Ara(j)
            UlS(i)=UlS(i)+UlA(i,j)*Ara(j)
         end do loop_col
        if(twoJ == 0) cycle
        G=G+UlS(i)*UlS(i)*GlS(i)
      end do loop_row


        if(twoJ == 0) GOTO 40
        do 38 i=1,N
        GAA(i)=0
        do 38 j=1,N
38      GAA(i)=GAA(i)+UlA(j,i)*GlS(j)*UlA(j,i)

        UjjM=0.
        UkkM=0.
        UlSM=0.
        AraM=0.

        do i=1,N
        if(UjjM.lT.Ujj(i)*Ujj(i)) then
          UjjM=Ujj(i)*Ujj(i)
          iUjjM=i
        endif
        if(UkkM.lT.Ukk(i)*Ukk(i)) then
          UkkM=Ukk(i)*Ukk(i)
          iUkkM=i
        endif

        if(UlSM.lT.UlS(i)*UlS(i)) then
          UlSM=UlS(i)*UlS(i)
          iUlSM=i
        endif

        if(AraM.lT.Ara(i)*Ara(i)) then
          AraM=Ara(i)*Ara(i)
          iUAlM=i
        endif
        enddo
   end subroutine lande
