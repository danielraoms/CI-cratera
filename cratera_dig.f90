!*************************************************************
module cratera_dig
	use parametros
	implicit none
	integer :: Ixcell_hole, Iycell_hole
	integer :: left_hole, right_hole, bottom_hole, upper_hole !define as fronteiras da região a ser cavada
	double precision :: Rcell

contains
	
subroutine cavar_cratera(N_resto, flag_dummy)
	use parametros
	implicit none
	integer :: contdig
	integer, intent(out) :: N_resto
	integer, dimension(N), intent(out) :: flag_dummy

	contdig = 0 !contador de partículas que ainda permanecem no sistema

	Rcell = 2.0d0*raiomax !o raio da célula é o diâmetro das partículas das paredes (fronteiras) do domínio 

	!definindo que todas as partículas móveis serão evoluídas e aparecerão no .eps
	do j = 1, N
		flag_dummy(j) = 1
	end do


	!identificador para a coluna de células que define a fronteira esquerda do buraco
	left_hole = 30 !um terço do comprimento da parede inferior das células à esquerda não serão cavadas

	!identificador para a coluna de células que define a fronteira direita do buraco
	right_hole = 99 !um terço do comprimento da parede inferior das células à direita não serão cavadas

	!identificador para a linha de células que define a fronteira inferior do buraco
	bottom_hole = 1 !a altura de aproximadamente 20 partículas da fronteira, partindo da parede inferior, não serão cavadas

	!identificador para a linha de células que define a fronteira superior do buraco
	upper_hole = 50 !acima da altura de 45 partículas, partindo da parede inferior, não há partículas

!		++             +                                        +
!		||            ++                                        |
!		||            ||                                        |
!		||            |-----------------------------------------------+ Upper hole
!		||            ||                                        |
!		||            ||                                        |
!		||            ||                                        |
!		||            ||                                        |
!		||            ||                                        |
!		||            ||                                        |
!		||            ||                                        |
!		||            ||                                        |
!		||            ||                                        |
!		||            ||                                        |
!		|--------------------------------------------------------------+ Bottom Hole
!		+-------------------------------------------------------+
!			       |                                        |
!			       +                                        +
!			    Left hole                                 Right hole
!É UM BATALHA NAVAL! MUDE U, B, R, L COMO QUISER


	!identificando quais partículas serão cavadas
	do j = 1, N
		!identificando as células das partículas
		Ixcell_hole = int((xold(j) - (xinicial-raiomax))/(Rcell) + 1) 
		Iycell_hole = int((yold(j) - (yinicial-raiomax))/(Rcell) + 1) 

		!se a partícula está dentro da região a ser cavada, mude seu flag_dig
		!se a partícula está dentro das colunas de células da região a ser cavada e
		if ((Ixcell_hole .ge. left_hole) .AND. (Ixcell_hole .le. right_hole)) then
			!se a partícula está dentro das linhas de células da região a ser cavada
			if ((Iycell_hole .gt. bottom_hole) .AND. (Iycell_hole .lt. upper_hole)) then
				flag_dummy(j) = -1
				contdig = contdig + 1
			end if
		end if
	end do

	N_resto = N - contdig

	write(*,*) "terminou dig"

	return
end subroutine cavar_cratera

end module
