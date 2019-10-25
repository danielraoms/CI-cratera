!define o número total de partículas do sistema
subroutine parts_number(qtdec,qtdel)
	use parametros
	implicit none
	double precision :: distx !distância horizontal entre dois centros de partículas adjacentes, considerando a folga entre partículas
	double precision :: disty !distância vertical entre dois centros de partículas adjacentes, considerando a folga entre partículas
	integer :: qtdel !quantidade máxima de partículas em cada linha de partículas alinhadas
	integer :: qtdec !quantidade máxima de partículas em cada coluna de partículas alinhadas

	distx = 2.0d0*raiomax + folgax
	disty = 2.0d0*raiomax + folgay

	qtdec = floor((xold(bw) - xold(1) - distx)/distx) 
	qtdel = floor((yold(bw+lw) - yold(1))/disty) 

	!escreve o número de linhas e colunas de partículas na geração da CI
	write(*,*) "qtde", qtdel, qtdec

	N = paredes + qtdel*(2*qtdec)
	
	!escreve o número de partículas na tela
	write(*,*) "N", N

	return
end subroutine parts_number

!*********************************************************************************************************

!define parâmetros da PAREDE (raio das partículas na configuração inicial e forças, posições e velocidades iniciais) 
subroutine condicoes_iniciais_wall()
	use parametros
	implicit none

	!define parâmetros gerais das partículas da parede
	do i = 1, paredes
		r(i) = raiomax !raio da partícula i
		forcax(i) = 0.0d0
		forcay(i) = 0.0d0
	enddo
	
      !todas as paredes fixadas
	!todas as partículas móveis começam do repouso
    	do i = 1, N
		vxold(i) = 0.0d0
		vyold(i) = 0.0d0
	end do

	!parede inferior
	do i = 1, bw
		xold(i) = xinicial + 2.0*(i-1)*raiomax
		xnew(i) = xinicial + 2.0*(i-1)*raiomax
		yold(i) = yinicial
		ynew(i) = yinicial
	enddo
	
	!parede esquerda
	do i = (bw + 1), (bw + lw)
		xold(i) = xinicial
		xnew(i) = xinicial
		yold(i) = yinicial + 2.0d0*(i-bw)*raiomax
		ynew(i) = yinicial + 2.0d0*(i-bw)*raiomax
	enddo
	
	!parede direita
	do i = (bw + lw + 1),(bw + lw + rw)
		xold(i) = xold(bw)
		xnew(i) = xnew(bw)
		yold(i) = yinicial + 2.0d0*(i-(bw+lw))*raiomax
		ynew(i) = yinicial + 2.0d0*(i-(bw+lw))*raiomax
	enddo
	
	return
end subroutine condicoes_iniciais_wall

!*********************************************************************************************************

!subroutine que define parâmetros das partículas MÓVEIS (raio, massa e inércia das partículas na configuração inicial e forças, posições e velocidades iniciais) 
subroutine condicoes_iniciais_position_moving()
	use parametros
	implicit none
	double precision :: rho !densidade do material de que é feito uma partícula
	double precision :: p1, p2 !intermediários
	double precision :: distx !distância horizontal entre dois centros de partículas adjacentes, considerando a folga entre partículas
	double precision :: disty !distância vertical entre dois centros de partículas adjacentes, considerando a folga entre partículas
	integer :: qtdel !quantidade máxima de partículas em cada linha de partículas alinhadas
	integer :: qtdec !quantidade máxima de partículas em cada coluna de partículas alinhadas
	integer :: k, O, parts
	integer :: clock
	integer, dimension(:), allocatable :: iseed
	real :: randD1(N), randD2(N), rand1(N), rand2(N), u_lg1(N), u_lg2(N)
	real :: rand3, rand4, rand5, s_aux
	real :: rand_aux(N) !variável auxiliar para gerar tamanhos de grãos com distribuição log-normal
	real :: sqrt_aux
	real, parameter :: v_small = TINY(1.0)
	double precision :: mu_aux, sigma_aux !mean and deviation for log-normal distribution

	!inicia o gerador de random numbers
	call random_seed(size=O)
	allocate(iseed(O))
	
	call system_clock(COUNT=clock)

	iseed = clock + 37*[(i, i=1,O)]
	call random_seed(PUT=iseed)
	deallocate(iseed)

	allocate(x_dummy(qtde_l,4*qtde_c))
	allocate(y_dummy(qtde_l,4*qtde_c))
	allocate(r_dummy(qtde_l,4*qtde_c))	

	distx = 2.0d0*raiomax + folgax
	disty = 2.0d0*raiomax + folgay

	!parâmetros para distribuição de tamanho das partículas
	mu_aux = -0.5
	sigma_aux = 0.9

	!posicionando as partículas MÓVEIS dentro do recipiente com um leve desvio
	k = paredes + 1
	do j = 1, 2*qtde_c
		do i = 1, qtde_l
		if (mod(j,2) .EQ. 0) then
			x_dummy(i,j) = xinicial + i*distx
		else 	
			call random_number(rand1(k))
			call random_number(randD1(k))
			x_dummy(i,j) = 0.005*rand1(k)*(-1)**(int(10*randD1(k))) + xinicial + i*distx
		end if 

		if (mod(i,2) .EQ. 0) then
			y_dummy(i,j) = yinicial + j*disty
		else 	
			call random_number(rand2(k))
			call random_number(randD2(k))
			y_dummy(i,j) = 0.005*rand2(k)*(-1)**(int(10*randD2(k))) + yinicial + j*disty
		end if 
		k = k+1
		end do
	end do

	!gerando distribuição normal de números randômicos
	!tome N-(paredes) números randômicos
	do k = (paredes + 1), N
		!chamar números aleatórios até que a soma deles dê próxima de 1
		do			
			call random_number(rand3)
			call random_number(rand4)

			rand3 = scale(rand3,1) - 1.0
			rand4 = scale(rand4,1) - 1.0

			s_aux = rand3**2.0d0 + rand4**2.0d0 + v_small*2.0d0 !somar v_small impede calcularmos log de zero/0

			if (s_aux .lt. 1.0) then
					exit
			end if
		end do
			sqrt_aux = sqrt(-2.0d0*log(s_aux)/s_aux)
			u_lg1(k) = rand3*sqrt_aux
			u_lg2(k) = rand4*sqrt_aux	
	end do

	!gerando a distribuição log-normal de números randômicos a partir de uma distribuição normal
	do k = (paredes + 1), N
		u_lg1(k) = exp(mu_aux + sigma_aux*u_lg1(k))
		u_lg2(k) = exp(mu_aux + sigma_aux*u_lg2(k))

		u_lg1(k) = u_lg1(k) - floor(u_lg1(k))
		u_lg2(k) = u_lg2(k) - floor(u_lg2(k))
	end do

	!renomenando as posições das partículas MÓVEIS (array bidimensional --> array unidimensional)
	k = (paredes + 1)
	do j = 1, 2*qtde_c
		do i = 1, qtde_l 
		xold(k) = x_dummy(i,j)
		xnew(k) = x_dummy(i,j)
		yold(k) = y_dummy(i,j) 
		ynew(k) = y_dummy(i,j)
		theta_old(k) = 0.0d0
		theta_new(k) = 0.0d0
		omega_old(k) = 0.0d0
		omega_new(k) = 0.0d0
		r(k) = 0.75d0*raiomax + (0.1*raiomed)*u_lg1(k)*(-1)**(int(10*u_lg2(k)))
		k = k+1
		end do
	end do

	!definindo massas, raios e momentos de inércia de cada partícula do SISTEMA
	rho = 8000.0d0 !densidade do acrílico
	p1 = (4.0d0/3.0d0)*pi
	do i = 1, N
		p2 = (r(i))**(3.0d0)
		m(i) = rho*(p1*p2)
		inertia(i) = 0.5d0*m(i)*((r(i))**2.0d0)
	end do

	deallocate(x_dummy)
	deallocate(y_dummy)
	deallocate(r_dummy)

	return
end subroutine condicoes_iniciais_position_moving
