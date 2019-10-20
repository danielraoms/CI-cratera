module parametros
	implicit none
	
	!arquivos
	character(100) :: entrada, posvel, overlaps, coef, kinenergy, potgenergy, rotenergy, totenergy, gainenergy, lossenergy !arquivos .dat	
	character(100) :: kinenergyc, potgenergyc, rotenergyc, totenergyc, cratera
	character(100) :: teste, forcasx, forcasy, forcas, dados_condini, condini, verlet_debug, simulation_time      !arquivos .dat
	character(100) :: validation_friction
	character(100) :: normalforces, posicoes_val, historia_atrito, cratera_evolucao


	!parâmetros de entrada 
	integer :: bw 							!número de partículas na BOTTOM WALL
	integer :: lw  							!número de partículas na LEFT WALL, sem contar a bottom wall
	integer :: rw  							!número de partículas na RIGHT WALL, sem contar a bottom wall nem left wall
	integer :: contmod  						!salva resultados a cada contmod passos
	integer :: contmodeps  						!salva arquivo .eps a cada contmodeps passos
	integer :: contmodsafe 						!salva energia a cada contmodsafe passos
	double precision :: dt 						!passo de tempo 
	double precision :: g  						!constante gravitacional terrestre
	double precision :: xinicial  					!posição x do centro da partícula mais perto da origem (0,0)
	double precision :: yinicial  					!posição y do centro da partícula mais perto da origem (0,0)
	double precision :: raiomed  					!raio médio das partículas do SISTEMA
	double precision :: raiomax 	 				!raio máximo, das partículas da parede, é 105% do raio médio
	double precision :: folgax 					!folga horizontal entre duas partículas adjacentes da CONFIGURAÇÃO INICIAL
	double precision :: folgay 		 			!folga vertical entre duas partículas adjacentes da CONFIGURAÇÃO INICIAL 
	double precision :: xwall, ywall				!placa horizontal plana - chão
	double precision :: vxwall, vywall				!velocidade da placa horizontal plana - chao
	double precision :: forca_x, forca_y, torque_s			!arquivos auxiliares para validaçao

	!parâmetros principais
	integer :: qtde_l 						!quantidade máxima de partículas em cada linha de partículas alinhadas
	integer :: qtde_c 						!quantidade máxima de partículas em cada linha de partículas alinhadas
	integer :: N 							!número de partículas na SIMULAÇÃO
	integer :: paredes 						!número de partículas nas PAREDES
	double precision :: t 						!tempo
	double precision :: gama_n, mi_t, gama_s, mi_roll, mi_roll_1, mi_roll_2 !coeficientes viscoso e de atrito dinâmico, coeficiente de regularização da força de atrito, coeficiente de rolling friction (modelo 1), coeficiente de rolling friction (modelo 2)
	double precision :: sinal_tan
	integer :: N_restante !número de partículas restante após cavar

	!lista de verlet
	integer ::  maxIxcell, maxIycell	    											 !número máximo de células na lista de Verlet
	integer, allocatable, dimension (:,:) :: tdl 										 !Tête de Liste
	integer, allocatable, dimension (:) :: links 										 !next tdl
	double precision, allocatable, dimension(:) :: pos_x_1, pos_y_1, chains_x, chains_y
	double precision, allocatable, dimension(:) :: normald_x, normald_y, raio_chains

	!para condições iniciais
	integer :: i, j																			!contadores de loop
	integer :: nummax  																		!número máximo de partículas na SIMULAÇÃO
	double precision, allocatable, dimension (:) :: r, m, inertia 							!raio, massa e momento de inércia de cada partícula
	double precision, allocatable, dimension (:) :: xold, yold, xnew, ynew, xnewer, ynewer  !posições x e y de cada partícula
	double precision, allocatable, dimension (:) :: vxold, vyold, vxnew, vynew				!velocidades x e y de cada partícula
	double precision, allocatable, dimension (:) :: theta_old, theta_new	       			!deslocamento angular de cada partícula
 	double precision, allocatable, dimension (:) :: omega_old, omega_new					!velocidade angular de cada partícula
	double precision, allocatable, dimension (:) :: forcax, forcay 							!forças x e y atuantes em cada partícula
	double precision, allocatable, dimension (:) :: torque 									!torque atuante em cada partícula
	double precision, allocatable, dimension (:,:) :: x_dummy, y_dummy, r_dummy				!arrays intermediários para posicionar partículas móveis na caixa
	
	!para a cratera
	integer, allocatable, dimension(:) :: flag_dig !flag para sinalizar as partículas a serem cavadas
	
	!post-processing
	double precision :: massa_media  										!massa média das partículas móveis
	double precision :: vel_maxima   										!velocidade máxima da simulação
	double precision :: forca_maxima 										!força máxima atingida considerando toda a evolução
	double precision :: kinetic 											!energia cinética total do sistema
	double precision :: potentialg 											!energia potencial gravitacional total do sistema
	double precision :: rotational_en 										!energia rotacional total do sistema
	double precision :: totalenergy 										!energia total do sistema
	double precision :: aux_pack											!auxiliar para calcular packing fraction
	integer, allocatable, dimension (:,:) :: cont_cell						!número de partículas em cada célula de Verlet
	integer :: sum_cell														!
	double precision, allocatable, dimension (:,:,:) :: sum_vcell 			!soma das velocidades das partículas para cada célula do sistema
	double precision, allocatable, dimension (:,:) :: highest_of_cell   	!altura máxima entre as partículas para cada célula do sistema
	double precision, allocatable, dimension (:,:,:) :: mean_velocity_cell  !velocidade média de cada célula
	double precision, allocatable, dimension (:,:,:) :: position_cell 		!posição do vetor velocidade em cada célula 
	double precision, allocatable, dimension (:,:) :: sum_packcell 			!soma das áreas das partículas de uma célula, para cada célula
	double precision, allocatable, dimension (:,:) :: packing_fraction 		!packing fraction de cada célula
	
	!para cálculo das forças
	double precision :: Fx_elastica, Fy_elastica  									!força normal elástica nas posições x e y
	double precision :: Fx_viscosa, Fy_viscosa    									!força normal viscosa nas posições x e y
	double precision :: Fs_tangencial	  											!força tangencial 
	double precision :: Fat_x, Fat_y												!força de atrito nas posições x e y
	double precision, allocatable, dimension (:) :: F_elastica 						!força elástica total de cada partícula (para .eps)
	double precision :: ndx, ndy													!componentes x e y da direção normal do contato
	double precision :: T_rolling     												!torque de rolling friction
	double precision, allocatable, dimension (:,:) :: dx_history_x, dx_history_y 	!valores da história do atrito para um dado contato 
	integer, allocatable, dimension(:,:) :: detector_old, detector_new				!detectores de colisão para a história do atrito 
	double precision ::  sinal_vrel 												!sinal da velocidade relativa tangencial

	!dados das partículas
	double precision, allocatable, dimension (:) :: E_young						!Young's modulus de cada partícula
	double precision, allocatable, dimension (:) :: v_poisson					!razão de Poisson de cada partícula
	double precision, allocatable, dimension (:) :: G_shear 					!shear modulus de cada partícula (para cálculo do atrito)

	!para gerar eps
	double precision, allocatable, dimension (:,:) :: velocidade_total


	!outros parâmetros 
	double precision, parameter :: pi = 4*atan(1.0d0)		!pi
	double precision :: toleranciakin 				!tolerância da energia cinética para parada - equilíbrio do SISTEMA
	double precision :: toleranciarot 				!tolerância da energia rotacional para parada - equilíbrio do SISTEMA

end module parametros
