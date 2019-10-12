!cria o setup da Lista de Verlet
subroutine lista_verlet(tdl_dummy)
	use parametros
	integer :: k, row, col 		!r for row, c for column
	integer :: Ixcell, Iycell 	!coluna da célula em questão, linha da célula em questão
	integer, dimension(N,N), intent(inout) :: tdl_dummy	   !última partícula colocada na célula em questão (Tête de Liste)
	double precision :: Rcell

	!iniciar a Tête de Liste (tdl)
	do col = 1, maxIxcell + 1
		do row = 1, maxIycell + 1
			tdl_dummy(row,col) = 0
		end do
	enddo	

	Rcell = (2.1d0*raiomed) !o raio da célula é o diâmetro das partículas das paredes (fronteiras) do domínio 
	
	!Localizar partículas e por nas respectivas células
	do k = 1, N
	
		!definindo a LINHA de células em que se encontra a partícula
		Iycell = int((yold(k) - (yinicial-raiomax))/(Rcell) + 1) 

		!definindo a COLUNA de células em que se encontra a partícula
		Ixcell = int((xold(k) - (xinicial-raiomax))/(Rcell) + 1)

		!a partícula k encontra-se em tdl(linha de células da partícula k, coluna de célula da partícula k) 
		links(k) = tdl_dummy(Iycell, Ixcell) 

		!a partícula k é a nova tête de liste da célula (Iycell, Ixcell)
		tdl_dummy(Iycell, Ixcell) = k 

		
		if (Ixcell .GT. maxIxcell) then
			maxIxcell = Ixcell
		end if
		
		if (Iycell .GT. maxIycell) then
			maxIycell = Iycell
		end if

		!redefine coluna e linha de células para cálculo subsequente
		Ixcell = 0
		Iycell = 0
	enddo

	return
end subroutine lista_verlet
