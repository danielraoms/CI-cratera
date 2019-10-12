include 'parametros.f90'		!parametros necessários para todo o main
include 'coefficients.f90'		!define os coeficientes necessários para cálculos no module forces
include 'forces.f90'			!define as subroutines de forças a serem usadas
include 'initial_conditions.f90'	!define as condições iniciais das partículas
include 'verlet_list.f90'		!contém subroutine para início da lista de verlet
include 'integracao_eps.f90'		!contém subroutines para evolução das partículas e salvar eps
include 'cratera_dig.f90'		!contém subroutine para cavar uma região a partir da condição inicial

!*************************************************************
!*************************************************************
PROGRAM Colisaoverlet2D
	use parametros 		
	use integracao_eps	
	use cratera_dig	

	implicit none
	integer :: l
	integer :: a, b !contadores para célula teste (a,b)
	integer :: p, q !contadores para partículas teste (p) e comparada (q)
	integer :: cont !contador do tempo 
	integer :: cont_chains !contador para force chains
	double precision :: start, finish !auxiliares para contar tempo de simulação

	!arquivo para dados de entrada
	write(entrada,"(a,i0,a)") "dados_entrada.dat"

	!arquivo para dados dos coeficientes
	write(coef,"(a,i0,a)") "coefficients.dat"
	
	!arquivo para dados da condição inicial
	write(dados_condini, "(a,i0,a)") "dados_condini.dat"
	
	!arquivo para quantidades da condição inicial
	write(condini, "(a,i0,a)") "condini.dat"

	!arquivo para quantidades da cratera
	write(cratera, "(a,i0,a)") "cratera.dat"

	!arquivo para história do atrito
	write(historia_atrito, "(a,i0,a)") "historia_atrito.dat"

	!arquivo para energias
	write(kinenergy, "(a,i0,a)") "energiacinetica.dat"
	write(potgenergy, "(a,i0,a)") "energiapotg.dat"
	write(rotenergy, "(a,i0,a)") "energiarot.dat"
	write(totenergy, "(a,i0,a)") "energiatotal.dat"

	!arquivo para validação da força normal
	write(normalforces, "(a,i0,a)") "normalforces.dat"

	!arquivo para validação das posições para colisões normais
	write(posicoes_val, "(a,i0,a)") "posicoes.dat"

	!abre arquivo para dados de entrada
	open(unit = 1, file = trim(entrada), status = "unknown")

	!abre arquivo para energias
	open (unit = 50, file = trim(kinenergy), status = "unknown")
	open (unit = 60, file = trim(potgenergy), status = "unknown")
	open (unit = 70, file = trim(rotenergy), status = "unknown")
	open (unit = 80, file = trim(totenergy), status = "unknown")

	!abre arquivo para validação da força normal
	open (unit = 90, file = trim(normalforces), status = "unknown")	

	!abre arquivo para validação das posições para colisões normais
	open (unit = 91, file = trim(posicoes_val), status = "unknown")

	!escreve e abre arquivo para validaçao da força de atrito
	write(validation_friction, "(a,i0,a)") "atrito.dat"
	open (unit = 1500, file = trim(validation_friction), status="unknown")

	


call cpu_time(start)

!*******************************************************************************************************************************************************************************************
!*****PRÉ-EVOLUÇÃO DO SISTEMA (GERANDO A CONFIGURAÇÃO INICIAL)******************************************************************************************************************************

	!puxa dados do arquivo de entrada para definir parâmetros de entrada da SIMULAÇÃO
	read(1,*) bw			!número de partículas na BOTTOM WALL
	read(1,*) lw			!número de partículas na LEFT WALL, sem contar a bottom wall
	read(1,*) rw			!número de partículas na RIGHT WALL, sem contar a bottom wall e a left wall
	read(1,*) contmod		!salva resultados a cada contmod passos
	read(1,*) contmodeps		!salva arquivo .eps a cada contmodeps passos
	read(1,*) contmodsafe		!salva energia a cada contmodsafe passos
	read(1,*) dt			!passo de tempo
	read(1,*) g 			!constante gravitacional
	read(1,*) xinicial		!posição x do centro da partícula mais perto da origem (0,0)
	read(1,*) yinicial		!posição y do centro da partícula mais perto da origem (0,0)
	read(1,*) raiomed		!raio médio das partículas do SISTEMA

	write(*,*) "entrada", bw, lw, rw, contmod, contmodeps, contmodsafe, dt, g, xinicial, yinicial, raiomed
	
	!define raio máximo das partículas do SISTEMA
	raiomax = 1.05d0*raiomed 

	!define folgas horizontais e verticais entre duas partículas adjacentes da CONFIGURAÇÃO INICIAL
	folgax = 0.5d0*raiomax !folga horizontal
	folgay = 0.5d0*raiomax  !folga vertical	
 
	!define o número de partículas na PAREDE
	paredes = bw + lw + rw

	!alocando posições, importante para definir número total de partículas no SISTEMA
	allocate(r(paredes))
	allocate(xold(paredes), yold(paredes))
	allocate(xnew(paredes), ynew(paredes))
	allocate(vxold(paredes), vyold(paredes))
	allocate(forcax(paredes), forcay(paredes))

	!definindo as posições iniciais das partículas da PAREDE para cálculo do número total de partículas do SISTEMA
	call condicoes_iniciais_wall()

	!calcula o número total de partículas do SISTEMA (importante para a alocação dos arrays)
	call parts_number(qtde_l, qtde_c)

	!desalocando posições, agora que sabemos o número total de partículas no SISTEMA
	deallocate(r)	
	deallocate(xold, yold)
	deallocate(xnew, ynew)
	deallocate(vxold, vyold)
	deallocate(forcax,forcay)

	!alocando posições e outros parâmetros físicos importantes do SISTEMA corretamente:
	!arrays que devem estar alocados a todo o momento
	allocate(r(N), m(N), inertia(N))
	allocate(xold(N), xnew(N))
	allocate(yold(N), ynew(N))
	allocate(xnewer(N), ynewer(N))
	allocate(vxold(N), vyold(N))
	allocate(vxnew(N), vynew(N))
	allocate(theta_old(N), theta_new(N))
	allocate(omega_old(N), omega_new(N))
	allocate(forcax(N), forcay(N), torque(N))
	allocate(tdl(N,N), links(N)) 

	!arrays necessários para gerar arquivos .eps
	allocate(F_elastica(N))
	allocate(velocidade_total(2,N))
	
	!define a configuração inicial das partículas da PAREDE
	call condicoes_iniciais_wall() 		  

	!define a configuração inicial das partículas MÓVEIS
	call condicoes_iniciais_position_moving()

	!definindo as dimensões dos arrays para a história do atrito
	allocate(detector_old(N,N), detector_new(N,N))
	allocate(dx_history_x(N,N), dx_history_y(N,N))
	allocate(E_young(N), v_poisson(N), G_shear(N)) 

	write(*,*) "AQUI", paredes+1, vxold(paredes+1), vyold(paredes+1) 

	!inicializando número máximo de células de Verlet
	maxIxcell = bw + 1
	maxIycell = bw + 1

	!escreve o número total de partículas do SISTEMA, número máximo de colunas de células de Verlet, número máximo de linhas de células de Verlet
	write(*,*) "N", N, maxIycell, maxIxcell

	!define a lista de Verlet para calcular corretamente número máximo de células de Verlet
	call lista_verlet(tdl)

	!calcula os coeficientes físicos para a simulação
	call coeficientes()

	!escreve os coeficientes da SIMULAÇÃO
	write(*,*) "coeficientes", gama_n, mi_t, gama_s, mi_roll_1, mi_roll_2

	!define o modelo de rolling friction
	mi_roll = mi_roll_1

	!a tolerância para o equilíbrio de energia cinética 
	!é igual à energia cinética de uma partícula da parede quase no repouso multiplicado pelo número de partículas móveis
	toleranciakin = (N - paredes)*m(1)*0.4d0*(1.0e-1)/2.0d0

	!a tolerância para o equilíbrio de energia rotacional
	!é igual à energia rotacional de uma partícula da parede quase no repouso multiplicado pelo número de partículas móveis
	toleranciarot = (N - paredes)*m(1)*0.4d0*(1.0e-1)/2.0d0

	!definindo módulos de Young, razão de Poisson e Shear modulus
	E_young(:) = 7.0e8
	v_poisson(:) = 0.2
	G_shear(:) = E_young(:)/(2.0d0*(1.0d0 + v_poisson(:)))

	!escreve a tolerância para as energias cinética e rotacional
	write(*,*) "tolerancia energia cinetica", toleranciakin
	write(*,*) "tolerancia energia rotacional", toleranciarot
	!read(*,*)
	

	!Definindo placa plana inferior - para validaçao do atrito
	xwall = xold(paredes+1)
	ywall = 0.5d0 + 2.0d0*raiomax
	vxwall = 0.0d0
	vywall = 0.0d0
	
	!iniciar a Friction History 
	detector_old(:,:) = 0
	detector_new(:,:) = 0


!**********************************************************************************************************************************************************************************************
!*****EVOLUÇÃO DO SISTEMA (GERANDO A CONDIÇÃO INICIAL)*****************************************************************************************************************************************
!loop do tempo
DO cont = 1, int(50/dt)

	!calcula as energias do sistema no tempo atual 
	call energias()

	if (mod(cont,100) .eq. 0) then
		write(*,*) "t", cont*dt, kinetic, rotational_en 
	end if


	!escreve o tempo atual na tela
	!escreve as energias em seus respectivos arquivos .dat
	if (mod(cont,10) .eq. 0) then
		write(50,*) cont*dt, kinetic
		write(60,*) cont*dt, potentialg
		write(70,*) cont*dt, rotational_en
		write(80,*) cont*dt, totalenergy
	end if

	!escreve os parâmetros e dados da configuração do sistema a cada contmodsafe passos, para manter a simulação segura de shutdowns
	if (mod(cont,contmodsafe) .eq. 0) then
			!escreve num arquivo .dat os parâmetros da condição inicial
			open (unit = 100, file=trim(dados_condini), status = "unknown")

			write(100, *) cont*dt,N,bw,lw,rw,xinicial,yinicial,raiomed,raiomax,dt,&
				      gama_n,mi_t,gama_s,mi_roll_1

			close(unit=100)

			!escreve num arquivo .dat as quantidades das partículas na condição inicial
			open (unit = 101, file=trim(condini), status = "unknown")

			do i = 1, N
				write(101,*) r(i),m(i),inertia(i),xold(i),xnew(i),vxold(i),vxnew(i),forcax(i),yold(i),ynew(i),&
			    		    vyold(i),vynew(i),forcay(i),theta_old(i),theta_new(i),omega_old(i),omega_new(i),torque(i),&
					    E_young(i),v_poisson(i),G_shear(i)
			end do

			close(unit=101)

			!escreve num arquivo .dat a história do atrito 
			open (unit = 102, file=trim(historia_atrito), status = "unknown" )

			do i = 1, N
			do j = 1, N 
				write(102,*) dx_history_x(i,j), dx_history_y(i,j), detector_old(i,j), detector_new(i,j)
			end do
			end do

			close(unit=102)		

			!itere enquanto (kinetic .LT. toleranciastop)
			if ((kinetic .LT. toleranciakin) .AND. (rotational_en .LT. toleranciarot) .AND. (cont .GT. 300000)) then
			
			write(*,*) "Equilíbrio:", t, kinetic, rotational_en
			!call execute_command_line ("gfortran cratera.f90")

			
			!termine a simulação
			EXIT
		end if
	end if

	!redefine as forças e torques do SISTEMA sofrendo apenas forças de campo para novo cálculo da iteração atual
	forcax(:) = 0.0d0
	forcay(:) = -g*m(:)
	torque(:) = 0.0d0

	!reinicia a lista de Verlet e reposiciona as partículas nas células
	!if (mod(cont,1) .EQ. 0) then
		call lista_verlet(tdl) 
	!end if

	!used for validation of friction of a single moving particle moving over a fixed plate 
	!xwall = xnew(paredes+1)

	
!	!itere todas as células interiores do sistema
	do b = 1, maxIxcell
	do a = 1, maxIycell
	
		!tome a partícula head da célula teste (a,b) (Head Of Cell, Hoc; Tête de Liste, TDL) como a partícula teste
		p = tdl(a,b)

		!itere todas as partículas da célula teste (a,b) 
		do while (p .GT. 0)

			!compare a célula atual (a,b) com ela mesma e com as células adjacentes
			do j = b-1, b+1
				!ignore as ghost cells (colunas) que seriam iteradas além das fronteiras do sistema
				if ((j .le. 0) .OR. (j .gt. maxIxcell)) then 
					cycle
				end if	
			do i = a-1, a+1
				!ignore as ghost cells (linhas) que seriam iteradas além das fronteiras do sistema
				if ((i .le. 0) .OR. (i .gt. maxIycell)) then
					cycle
				end if

			!tome a partícula head da célula (i,j) como a partícula comparada
			q = tdl(i,j)		
				
				!compare todas as partículas da célula comparada (i,j) 
				do while (q .NE. 0)

					!não calcule as quantidades das partículas da fronteira
					if (p .gt. paredes) then

					!não compare a partícula teste com ela mesma
					!critério de colisão entre a partícula teste (p) e a partícula comparada (q)
					if ((p.ne.q).AND.((((xnew(p)-xnew(q))**2.0d0)+(ynew(p)-ynew(q))**2.0d0).LE.((r(p)+r(q))**2.0))) then

						detector_new(p,q) = 1
		
						!Friction History Detection 
						!caso 1: início de uma colisão
						if ((detector_old(p,q) .eq. 0) .and. (detector_new(p,q) .eq. 1)) then						
							dx_history_x(p,q) = 0.0d0
							dx_history_y(p,q) = 0.0d0
						!caso 2: no meio de uma colisão
						else if ((detector_old(p,q) .eq. 1) .and. (detector_new(p,q) .eq. 1)) then
							dx_history_x(p,q) = dx_history_x(p,q)
							dx_history_y(p,q) = dx_history_y(p,q)
							!mantenha dx_history_x(p,q) e dx_history_y(p,q) unchanged
						end if
							 
						if ((p .le. paredes) .and. (q .le. paredes)) then
							!não calcule a força resultante para qualquer partícula da parede
							continue
						else
						  !calcule as forças 
						call all_forces(gama_n,mi_t,gama_s,mi_roll,E_young(p),E_young(q),v_poisson(p),v_poisson(q),&
									G_shear(p),G_shear(q),m(p),m(q),r(p),r(q),xold(p),xold(q),yold(p),yold(q),&
									vxold(p),vxold(p),vyold(p),vyold(q),omega_old(p),omega_old(q),Fx_elastica,Fy_elastica,&
									Fx_viscosa,Fy_viscosa,Fat_x,Fat_y,Fs_tangencial,T_rolling,ndx,ndy,&
									dx_history_x(p,q),dx_history_y(p,q),sinal_vrel)

						end if

						!soma as contribuições de força que a partícula comparada (q) realiza na partícula teste (p) à força total 
						!que a partícula teste (p) sofre de todas suas partículas vizinhas
						forcax(p) = forcax(p) + Fx_elastica + Fx_viscosa + Fat_x
						forcay(p) = forcay(p) + Fy_elastica + Fy_viscosa + Fat_y
						torque(p) = torque(p) - sinal_vrel*(Fs_tangencial)*r(p) + T_rolling

						!força elástica resultante na partícula p
						F_elastica(p) = F_elastica(p) + dsqrt(Fx_elastica**2.0d0 + Fy_elastica**2.0d0)

						!reseta as contribuições que a partícula comparada (q) realiza na partícula teste (p) - já foram contabilizadas!
						Fx_elastica = 0.0d0
						Fy_elastica = 0.0d0
						Fx_viscosa = 0.0d0
						Fy_viscosa = 0.0d0
						Fs_tangencial = 0.0d0
						ndx = 0.0d0
						ndy = 0.0d0

					else
						detector_new(p,q) = 0
					end if !fecha condicional de critério de colisão

					end if

					!atualiza a partícula comparada (q) como sendo a próxima da lista de vizinhos da célula comparada (i,j)
					q = links(q)

				end do !fecha loop de partículas da célula comparada (i,j)

			end do !fecha double loop de células comparadas (i,j), vizinhas da célula teste (a,b)
			end do

			!atualiza a partícula atual (p) como sendo a próxima da lista de vizinhos da célula teste (a,b)
			p = links(p)

		end do !fecha loop de partículas da célula teste (a,b)

	end do !fecha double loop de células teste (a,b)
	end do

	!evolui no tempo as posições, velocidades e ângulos das partículas
	call integracao_verlet()

	if (mod(cont,10) .eq. 0) then
		write(91,*) cont*dt, yold(12)
	end if

	do i = 1, N
		forcax(i) = 0.0d0
		forcay(i) = 0.0d0
		torque(i) = 0.0d0
	end do

	velocidade_total(1,:) = vxold(:)
	velocidade_total(2,:) = vyold(:)  

	!gera arquivo .eps com configuração do sistema na iteração atual
	!calcula as energias da configuração do sistema na iteração atual
	if (mod(cont,contmodeps) .eq. 0) then

		!chama subroutine para gerar arquivo .eps com configuração do sistema 
		call salva_eps(int(cont/contmodeps),ywall,N,paredes,r,xnew,ynew,0,theta_new,F_elastica,velocidade_total)

		call energias()
	end if

	!atualizando friction history detection arrays
	detector_old(:,:) = detector_new(:,:)

END DO

	write(*,*) "CI gerada."
	write(*,*) "Cavando cratera..."

!**********************************************************************************************************************************************************************************************
!*****CAVANDO A CRATERA (GERANDO A CONDIÇÃO INICIAL)*****************************************************************************************************************************************

	allocate (flag_dig(N))

	!cava uma região a partir da condição inicial
	call cavar_cratera(N_restante, flag_dig) 

	write(*,*) "N_restante", N, N_restante

	write(*,*) "Cratera cavada. Evoluindo..."

!loop do tempo
DO cont = floor(t/dt), 10000000	

	!calcula as energias do sistema no tempo atual 
	call energias()

	if (mod(cont,100) .eq. 0) then
		write(*,*) "t", cont*dt, kinetic, rotational_en 
	end if


	!escreve o tempo atual na tela
	!escreve as energias em seus respectivos arquivos .dat
	if (mod(cont,10) .eq. 0) then
		write(50,*) cont*dt, kinetic
		write(60,*) cont*dt, potentialg
		write(70,*) cont*dt, rotational_en
		write(80,*) cont*dt, totalenergy
	end if

	!escreve os parâmetros e dados da configuração do sistema a cada contmodsafe passos, para manter a simulação segura de shutdowns
	if (mod(cont,contmodsafe) .eq. 0) then

			!escreve num arquivo .dat as quantidades das partículas na condição inicial
			open (unit = 105, file=trim(cratera), status = "unknown")

			do i = 1, N
				write(105,*) r(i),m(i),inertia(i),xold(i),xnew(i),vxold(i),vxnew(i),forcax(i),yold(i),ynew(i),&
			    		    vyold(i),vynew(i),forcay(i),theta_old(i),theta_new(i),omega_old(i),omega_new(i),torque(i),&
					    E_young(i),v_poisson(i),G_shear(i)
			end do

			close(unit=105)

			!itere enquanto (kinetic .LT. toleranciastop)
			if ((kinetic .LT. toleranciakin/10.0d0) .AND. (rotational_en .LT. toleranciarot/10.0d0) .AND. (cont .GT. 3000000)) then
			
			write(*,*) "Equilíbrio:", t, kinetic, rotational_en
			!call execute_command_line ("gfortran cratera.f90")

			
			!termine a simulação
			EXIT
		end if
	end if

	!redefine as forças e torques do SISTEMA sofrendo apenas forças de campo para novo cálculo da iteração atual
	forcax(:) = 0.0d0
	forcay(:) = -g*m(:)
	torque(:) = 0.0d0

	do j = 1, N
		if (flag_dig(j) .eq. -1) then
			forcax(j) = 0.0d0
			forcay(j) = 0.0d0		
			torque(j) = 0.0d0
		end if
	end do

	!reinicia a lista de Verlet e reposiciona as partículas nas células
	!if (mod(cont,1) .EQ. 0) then
		call lista_verlet(tdl) 
	!end if

	!used for validation of friction of a single moving particle moving over a fixed plate 
	!xwall = xnew(paredes+1)

	
	!itere todas as células interiores do sistema
	do b = 1, maxIxcell
	do a = 1, maxIycell
	
		!tome a partícula head da célula teste (a,b) (Head Of Cell, Hoc; Tête de Liste, TDL) como a partícula teste
		p = tdl(a,b)

		!itere todas as partículas da célula teste (a,b) 
		do while (p .GT. 0)

			!compare a célula atual (a,b) com ela mesma e com as células adjacentes
			do j = b-1, b+1
				!ignore as ghost cells (colunas) que seriam iteradas além das fronteiras do sistema
				if ((j .le. 0) .OR. (j .gt. maxIxcell)) then 
					cycle
				end if	
			do i = a-1, a+1
				!ignore as ghost cells (linhas) que seriam iteradas além das fronteiras do sistema
				if ((i .le. 0) .OR. (i .gt. maxIycell)) then
					cycle
				end if

			!tome a partícula head da célula (i,j) como a partícula comparada
			q = tdl(i,j)		
				
				!compare todas as partículas da célula comparada (i,j) 
				do while (q .NE. 0)

					!não calcule as quantidades das partículas da fronteira
					if (p .gt. paredes) then

					if (((flag_dig(p) .eq. 1) .AND. (flag_dig(q) .eq. 1)) .OR.& !duas partículas móveis não-cavadas
					   ((flag_dig(p) .eq. 1) .AND. (q .le. paredes))      .OR.& !uma partícula móvel e uma da parede
					   ((p .le. paredes) .AND. (flag_dig(q) .eq. 1))) then	    !uma partícula da parede e uma partícula móvel


					!não compare a partícula teste com ela mesma
					!critério de colisão entre a partícula teste (p) e a partícula comparada (q)
					if ((p.ne.q).AND.((((xnew(p)-xnew(q))**2.0d0)+(ynew(p)-ynew(q))**2.0d0).LE.((r(p)+r(q))**2.0))) then

						detector_new(p,q) = 1
		
						!Friction History Detection 
						!caso 1: início de uma colisão
						if ((detector_old(p,q) .eq. 0) .and. (detector_new(p,q) .eq. 1)) then						
							dx_history_x(p,q) = 0.0d0
							dx_history_y(p,q) = 0.0d0
						!caso 2: no meio de uma colisão
						else if ((detector_old(p,q) .eq. 1) .and. (detector_new(p,q) .eq. 1)) then
							dx_history_x(p,q) = dx_history_x(p,q)
							dx_history_y(p,q) = dx_history_y(p,q)
							!mantenha dx_history_x(p,q) e dx_history_y(p,q) unchanged
						end if
							 
						if ((p .le. paredes) .and. (q .le. paredes)) then
							!não calcule a força resultante para qualquer partícula da parede
							continue
						else
						  !calcule as forças 
						call all_forces(gama_n,mi_t,gama_s,mi_roll,E_young(p),E_young(q),v_poisson(p),v_poisson(q),&
									G_shear(p),G_shear(q),m(p),m(q),r(p),r(q),xold(p),xold(q),yold(p),yold(q),&
									vxold(p),vxold(p),vyold(p),vyold(q),omega_old(p),omega_old(q),Fx_elastica,Fy_elastica,&
									Fx_viscosa,Fy_viscosa,Fat_x,Fat_y,Fs_tangencial,T_rolling,ndx,ndy,&
									dx_history_x(p,q),dx_history_y(p,q),sinal_vrel)

						end if

						!soma as contribuições de força que a partícula comparada (q) realiza na partícula teste (p) à força total 
						!que a partícula teste (p) sofre de todas suas partículas vizinhas
						forcax(p) = forcax(p) + Fx_elastica + Fx_viscosa + Fat_x
						forcay(p) = forcay(p) + Fy_elastica + Fy_viscosa + Fat_y
						torque(p) = torque(p) - sinal_vrel*(Fs_tangencial)*r(p) + T_rolling

						!força elástica resultante na partícula p
						F_elastica(p) = F_elastica(p) + dsqrt(Fx_elastica**2.0d0 + Fy_elastica**2.0d0)

						!reseta as contribuições que a partícula comparada (q) realiza na partícula teste (p) - já foram contabilizadas!
						Fx_elastica = 0.0d0
						Fy_elastica = 0.0d0
						Fx_viscosa = 0.0d0
						Fy_viscosa = 0.0d0
						Fs_tangencial = 0.0d0
						ndx = 0.0d0
						ndy = 0.0d0

					else
						detector_new(p,q) = 0
					end if !fecha condicional de critério de colisão

					end if

					end if

					!atualiza a partícula comparada (q) como sendo a próxima da lista de vizinhos da célula comparada (i,j)
					q = links(q)

				end do !fecha loop de partículas da célula comparada (i,j)

			end do !fecha double loop de células comparadas (i,j), vizinhas da célula teste (a,b)
			end do

			!atualiza a partícula atual (p) como sendo a próxima da lista de vizinhos da célula teste (a,b)
			p = links(p)

		end do !fecha loop de partículas da célula teste (a,b)

	end do !fecha double loop de células teste (a,b)
	end do

	!evolui no tempo as posições, velocidades e ângulos das partículas
	call integracao_verlet_cratera()

	do i = 1, N
		forcax(i) = 0.0d0
		forcay(i) = 0.0d0
		torque(i) = 0.0d0
	end do

	velocidade_total(1,:) = vxold(:)
	velocidade_total(2,:) = vyold(:)  

	!gera arquivo .eps com configuração do sistema na iteração atual
	!calcula as energias da configuração do sistema na iteração atual
	if (mod(cont,contmodeps) .eq. 0) then

		sum_vcell(:,:,:) = 0.0d0
		cont_cell(:,:) = 0
		highest_of_cell(:,:) = 0.0d0
		sum_cell = 0

		!calculando a soma das velocidades das partículas para cada célula do sistema
		!itere todas as células do sistema
		do b = 1, maxIxcell
		do a = 1, maxIycell

			!tome a partícula head da célula teste (a,b) como a partícula teste
			p = tdl(a,b)

			do while (p .gt. 0)

				if (flag_dig(p) .eq. 1) then

					sum_vcell(1,a,b) = sum_vcell(1,a,b) + vxold(p)
					sum_vcell(2,a,b) = sum_vcell(2,a,b) + vyold(p)

					!calculando a altura máxima para cada célula
					if (xold(p) .gt. highest_of_cell(a,b)) then
						highest_of_cell(a,b) = xold(p)
					else
						highest_of_cell(a,b) = highest_of_cell(a,b)
					end if

					cont_cell(a,b) = cont_cell(a,b) + 1

					p = links(p)
				else 
					sum_vcell(1,a,b) = sum_vcell(1,a,b)
					sum_vcell(2,a,b) = sum_vcell(2,a,b)
					p = links(p)
					continue
				end if
			end do	
		end do
		end do

		mean_velocity_cell(:,:,:) = 0.0d0

		!calculando a velocidade média de cada célula
		do b = 1, maxIxcell
		do a = 1, maxIycell
			sum_cell = 0
			sum_cell = sum_cell +&
			           cont_cell(a+1,b+1)  + cont_cell(a+1,b)   + cont_cell(a+1,b)&
				  +cont_cell(a,b-1)    + cont_cell(a,b)     + cont_cell(a,b+1)&
				  +cont_cell(a-1,b-1)  + cont_cell(a-1,b)   + cont_cell(a-1,b+1)

			if (sum_cell .ne. 0) then

				mean_velocity_cell(1,a,b) = (sum_vcell(1,a+1,b-1) + sum_vcell(1,a+1,b) + sum_vcell(1,a+1,b+1)&
							   +sum_vcell(1,a,b-1)   + sum_vcell(1,a,b)   + sum_vcell(1,a,b+1)&
							   +sum_vcell(1,a-1,b-1) + sum_vcell(1,a-1,b) + sum_vcell(1,a-1,b-1))/sum_cell

				mean_velocity_cell(2,a,b) = (sum_vcell(2,a+1,b-1) + sum_vcell(2,a+1,b) + sum_vcell(2,a+1,b+1)&
							   +sum_vcell(2,a,b-1)   + sum_vcell(2,a,b)   + sum_vcell(2,a,b+1)&
							   +sum_vcell(2,a-1,b-1) + sum_vcell(2,a-1,b) + sum_vcell(2,a-1,b-1))/sum_cell

			else
				mean_velocity_cell(1,a,b) = 0.0d0
				mean_velocity_cell(2,a,b) = 0.0d0
			end if

			position_cell(1,a,b) = 2.0d0*raiomed*b
			position_cell(2,a,b) = 2.0d0*raiomed*a
		end do
		end do


		!calculando velocidades médias da posição x para gráfico espaço-temporal
		do b = 1, maxIxcell
			space_time_vel(1,b) = sum(mean_velocity_cell(1,:,b))/(2.0d0*raiomax*39)
			space_time_vel(2,b) = sum(mean_velocity_cell(2,:,b))/(2.0d0*raiomax*39)

			position_cell_st(b) = b*2.0d0*raiomax
		end do

		!calculando altura máxima para cada coluna (para cada posição x)
		do b = 1, maxIxcell
			highest_height(b) = maxval(highest_of_cell(:,b))
		end do

		!arquivo para campo de velocidades
		write(campo_continuo, '("campo_continuo-", I4, ".dat")') cont/contmodeps
	
		!abre arquivo para campo de velocidades
		open(unit=337,file=campo_continuo,status='unknown')

		!escreve dados do campo de velocidades no instante atual
		do b = 1, maxIxcell
		do a = 1, maxIycell
			write(337,*) cont*dt, a, b, position_cell(1,a,b), position_cell(2,a,b), mean_velocity_cell(1,a,b), mean_velocity_cell(2,a,b)
			write(*,*) "fake writing velocity field"
		end do
		end do

		!abre arquivo para diagramas espaço-temporais
		open(unit=437, file=espaco_temporal,status='unknown')

		!escreve dados para diagramas espaço-temporais 
		do b = 1, maxIxcell
			write(437,*) cont*dt, position_cell_st(b), space_time_vel(1,b), space_time_vel(2,b), highest_height(b)
		end do
			write(437,*) ''


		!chama subroutine para gerar arquivo .eps com configuração do sistema 
		call salva_eps_cratera(int(cont/contmodeps),ywall,N,paredes,r,xnew,ynew,0,theta_new,F_elastica,velocidade_total,flag_dig)

		call energias()
	end if

	!atualizando friction history detection arrays
	detector_old(:,:) = detector_new(:,:)

END DO

	write(*,*) "Cratera ou degrau obtido/a."

	call cpu_time(finish)

	!escreve o tempo de simulação na tela
	write(101,*) '("Time = ", f6.3," minutes.")', (finish-start)/60, (finish-start)/3600
	write(*,*) '("Time = ", f6.3," minutes.")', (finish-start)/60, (finish-start)/3600

	write(*,*) "forca máxima", maximum_force

	write(*,*) "The simulation reached its end time."
end program Colisaoverlet2D
